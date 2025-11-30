##  This script generates Figure 4a and related panels.
##  This script performs:
##    1) PCA of gut content samples (Colon/Cecum/Rectum × WT/DB).
##    2) PLS-DA for:
##         - 6 groups (Colon/Cecum/Rectum × WT/DB),
##         - WT vs DB (all gut segments),
##         - Colon vs Cecum vs Rectum (WT only).
##    3) Differential analysis (DB vs WT) per tissue using limma.
##    4) UpSet plot of significant metabolites across tissues.
##    5) Faceted volcano plots (per tissue) and highlighting
##       core metabolites.
##
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(limma)
  library(factoextra)
  library(pheatmap)
  library(ggrepel)
  library(glue)
  library(viridis)
  library(ropls)
  library(UpSetR)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggplot2)
})

# ------------------ 1. Paths & basic setup ------------------
# Root working dir for this figure (adapt to your local path)
setwd("path/to/project_root/Figure4")

# Input Excel files (unchanged names / formats)
serum_file   <- "Metabolomics_Serum.xlsx"
content_file <- "Metabolomics_CecalContent.xlsx"

# Output folder for all plots / results
dir.create("plots", showWarnings = FALSE)

# ------------------ 2. Read & unify data --------------------
serum   <- read_xlsx(serum_file)
content <- read_xlsx(content_file)

# Fix a common header typo (keep original logic)
names(serum) <- str_replace_all(names(serum), "^Seurm", "Serum")

# Feature columns (kept exactly as in original script)
feature_cols <- c(
  "Type","ID","MZ","RT","MS2.name","MS2_score",
  "Formula","SuperClass","Class","Subclass",
  "HMDB","kegg","kegg_pathway","MS1.name"
)

# Long-format intensity tables
serum_long <- serum %>%
  pivot_longer(
    -all_of(feature_cols),
    names_to  = "Sample",
    values_to = "Intensity"
  )

content_long <- content %>%
  pivot_longer(
    -all_of(feature_cols),
    names_to  = "Sample",
    values_to = "Intensity"
  )

# Combine all (serum + content)
all_long <- bind_rows(serum_long, content_long)

# ------------------ 3. Sample metadata ----------------------
# New column name convention:
#   Colon-DB1, Colon-WT2, Cecum-DB3, Rectum-WT4, Serum-WT1, Serum-DB6, ...
# We parse:
#   Tissue_raw: part before first "-"
#   Genotype  : "WT" or "DB"
#   Rep       : trailing digits
#   Tissue    : copy of Tissue_raw (Colon/Cecum/Rectum/Serum)
meta <- all_long %>%
  distinct(Sample) %>%
  mutate(
    # Part before first dash, e.g. "Colon", "Cecum", "Rectum", "Serum"
    Tissue_raw = sub("-.*$", "", Sample),

    # Genotype encoded as "-WT" or "-DB" in the sample name
    Genotype = case_when(
      str_detect(Sample, "-WT") ~ "WT",
      str_detect(Sample, "-DB") ~ "DB",
      TRUE                      ~ NA_character_
    ),

    # Replicate ID = trailing digits (works for 1 or 2 digits)
    Rep = as.integer(str_extract(Sample, "\\d+$")),

    # Unified tissue label (already in "Colon/Cecum/Rectum/Serum" format)
    Tissue = Tissue_raw
  ) %>%
  dplyr::filter(!str_detect(Sample, "^QC"))  # drop QC samples if present

# Restrict main table to well-defined samples
all_long <- all_long %>%
  dplyr::filter(Sample %in% meta$Sample)

# ------------------ 4. Expression matrix --------------------
expr_mat <- all_long %>%
  dplyr::filter(!is.na(MS2.name)) %>%
  group_by(MS2.name, Sample) %>%
  summarise(Intensity = max(Intensity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Intensity) %>%
  distinct(MS2.name, .keep_all = TRUE) %>%
  column_to_rownames("MS2.name") %>%
  as.matrix()

write.csv(expr_mat, "expr_mat.csv")

# Log2 transform
expr_log <- log2(expr_mat + 1)

# ------------------ 5. Median imputation helper -------------
impute_median <- function(mat) {
  t(apply(mat, 1, function(x) {
    if (all(is.na(x))) return(x)
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  }))
}

# ============================================================
# 6. PCA – gut contents only (6 groups)
#    Colon/Cecum/Rectum × WT/DB
# ============================================================

content_samples <- meta %>%
  dplyr::filter(Tissue %in% c("Colon", "Cecum", "Rectum")) %>%
  mutate(Group = paste(Tissue, Genotype, sep = "_"))

mat_raw <- expr_log[, content_samples$Sample, drop = FALSE]
mat_raw <- mat_raw[rowSums(is.na(mat_raw)) < ncol(mat_raw), , drop = FALSE]

mat_imp <- t(apply(mat_raw, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}))
mat_imp <- mat_imp[apply(mat_imp, 1, sd) > 0, , drop = FALSE]

pca_res <- prcomp(t(mat_imp), scale. = TRUE)

group_colors <- RColorBrewer::brewer.pal(6, "Set1")
names(group_colors) <- unique(content_samples$Group)

p <- fviz_pca_ind(
  pca_res,
  geom        = "point",
  col.ind     = content_samples$Group,
  addEllipses = TRUE,
  ellipse.type = "confidence"
) +
  scale_color_manual(values = group_colors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title    = element_text(face = "bold"),
    legend.title  = element_blank()
  ) +
  ggtitle("PCA - Gut Content (6 Groups)")
p

ggsave("plots/PCA_GutContent_6groups.png", p, width = 7, height = 6, dpi = 300)
pdf("plots/PCA_GutContent_6groups.pdf", 7, 6); print(p); dev.off()

message("PCA (6 groups) saved.")

# ============================================================
# 7. PLS-DA – 6 groups (Colon/Cecum/Rectum × WT/DB)
# ============================================================

content_meta <- meta %>%
  dplyr::filter(Tissue %in% c("Colon", "Cecum", "Rectum")) %>%
  mutate(Group = paste(Tissue, Genotype, sep = "_"))

samples <- content_meta$Sample

expr_sub <- expr_log[, samples, drop = FALSE]
expr_sub <- expr_sub[rowSums(is.na(expr_sub)) < ncol(expr_sub), ]
expr_sub_imp <- t(apply(expr_sub, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}))
expr_sub_imp <- expr_sub_imp[apply(expr_sub_imp, 1, sd) > 0, ]

X <- t(expr_sub_imp)
Y <- content_meta$Group

set.seed(123)
if (any(table(Y) < 2)) {
  warning("PLS-DA (6 groups): At least one class has <2 samples; model may not run.")
}

plsda_mod <- opls(X, Y, predI = 2, orthoI = 0, permI = 100)

scores <- data.frame(
  PLS1  = plsda_mod@scoreMN[, 1],
  PLS2  = plsda_mod@scoreMN[, 2],
  Group = Y
)

colors <- brewer.pal(length(unique(Y)), "Set1")
names(colors) <- unique(Y)

p_pls <- ggplot(scores, aes(PLS1, PLS2, color = Group, shape = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(linetype = 2, color = "grey40") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title   = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "PLS-DA - Gut Content (Colon/Cecum/Rectum x WT/DB)",
    x     = "PLS-DA Component 1",
    y     = "PLS-DA Component 2"
  )
p_pls

ggsave("plots/PLSDA_Content_6groups.png", p_pls, width = 6, height = 5, dpi = 300)
pdf("plots/PLSDA_Content_6groups.pdf", 6, 5); print(p_pls); dev.off()

# ============================================================
# 8. Additional PLS-DA:
#    (1) WT vs DB – all contents
#    (2) Colon vs Cecum vs Rectum – WT only
# ============================================================

plot_plsda <- function(score_df, out_prefix, title) {
  group_colors <- brewer.pal(length(unique(score_df$Group)), "Set1")
  names(group_colors) <- unique(score_df$Group)
  
  p <- ggplot(score_df, aes(x = PLS1, y = PLS2, color = Group, shape = Group)) +
    geom_point(size = 3) +
    stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, color = NA) +
    stat_ellipse(linetype = 2, color = "grey40") +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x     = "PLS Component 1",
      y     = "PLS Component 2"
    ) +
    theme(
      legend.title = element_blank(),
      plot.title   = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(glue::glue("plots/{out_prefix}_score.png"), plot = p, width = 6, height = 5, dpi = 500)
  pdf(glue::glue("plots/{out_prefix}_score.pdf"), 6, 5); print(p); dev.off()
}

# --- 8.1 WT vs DB (all contents) ---
meta_sub1 <- meta %>%
  filter(Tissue %in% c("Colon", "Cecum", "Rectum")) %>%
  mutate(Group = Genotype)

samples1 <- meta_sub1$Sample
mat1 <- expr_log[, samples1, drop = FALSE]
mat1 <- mat1[rowSums(is.na(mat1)) < ncol(mat1), ]
mat1 <- impute_median(mat1)
mat1 <- mat1[apply(mat1, 1, sd) > 0, ]

X1 <- t(mat1)
Y1 <- meta_sub1$Group

set.seed(123)
if (any(table(Y1) < 2)) {
  warning("PLS-DA (WT vs DB): some classes have <2 samples.")
}
pls1 <- opls(X1, Y1, predI = 2, orthoI = 0, permI = 100)

score1 <- data.frame(
  PLS1  = pls1@scoreMN[, 1],
  PLS2  = pls1@scoreMN[, 2],
  Group = Y1
)
plot_plsda(score1, out_prefix = "PLSDA_WTvsDB", title = "PLS-DA: WT vs DB")

# --- 8.2 Colon vs Cecum vs Rectum (WT only) ---
meta_sub2 <- meta %>%
  filter(Tissue %in% c("Colon", "Cecum", "Rectum") & Genotype == "WT") %>%
  mutate(Group = Tissue)

samples2 <- meta_sub2$Sample
mat2 <- expr_log[, samples2, drop = FALSE]
mat2 <- mat2[rowSums(is.na(mat2)) < ncol(mat2), ]
mat2 <- impute_median(mat2)
mat2 <- mat2[apply(mat2, 1, sd) > 0, ]

X2 <- t(mat2)
Y2 <- meta_sub2$Group

set.seed(123)
if (any(table(Y2) < 2)) {
  warning("PLS-DA (WT tissues): some classes have <2 samples.")
}
pls2 <- opls(X2, Y2, predI = 2, orthoI = 0, permI = 100)

score2 <- data.frame(
  PLS1  = pls2@scoreMN[, 1],
  PLS2  = pls2@scoreMN[, 2],
  Group = Y2
)
plot_plsda(score2, out_prefix = "PLSDA_WT_tissues",
           title = "PLS-DA: Colon vs Cecum vs Rectum (WT only)")

# ============================================================
# 9. Differential metabolites (DB vs WT) per tissue
#    + UpSet plot + faceted volcano
# ============================================================

tissues   <- c("Serum", "Colon", "Cecum", "Rectum")
diff_list <- list()

for (tis in tissues) {
  cat("Processing:", tis, "\n")
  
  samp <- meta %>%
    dplyr::filter(Tissue == tis) %>%
    pull(Sample)
  mat  <- expr_log[, samp, drop = FALSE]
  
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), ]
  mat <- t(apply(mat, 1, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  }))
  mat <- mat[apply(mat, 1, sd) > 0, ]
  
  geno <- factor(
    meta$Genotype[match(colnames(mat), meta$Sample)],
    levels = c("WT", "DB")
  )
  design <- model.matrix(~ 0 + geno)
  colnames(design) <- c("WT", "DB")
  
  fit  <- lmFit(mat, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(DB - WT, levels = design)))
  
  tab <- topTable(fit2, number = Inf) %>%
    rownames_to_column("Metabolite")
  
  write_csv(tab, glue("plots/{tis}_DBvsWT_limma.csv"))
  
  sig <- tab %>%
    dplyr::filter(P.Value < 0.05, abs(logFC) > 0.5) %>%
    pull(Metabolite)
  diff_list[[tis]] <- sig
  cat(" Significant:", length(sig), "metabolites\n")
}

up_in <- fromList(diff_list)

png("plots/DiffMetabolites_UpSet.png", width = 2600, height = 2200, res = 350)
UpSetR::upset(
  up_in,
  sets           = tissues,
  order.by       = "freq",
  keep.order     = TRUE,
  sets.bar.color = "#4DBBD5",
  main.bar.color = "#E64B35",
  text.scale     = 2
)
dev.off()

pdf("plots/DiffMetabolites_UpSet.pdf", width = 6, height = 5)
UpSetR::upset(
  up_in,
  sets           = tissues,
  order.by       = "freq",
  keep.order     = TRUE,
  sets.bar.color = "#4DBBD5",
  main.bar.color = "#E64B35",
  text.scale     = 2
)
dev.off()

volcano_df <- bind_rows(lapply(tissues, function(tis) {
  read_csv(glue("plots/{tis}_DBvsWT_limma.csv"), show_col_types = FALSE) %>%
    mutate(Tissue = tis)
}))

volcano_df <- volcano_df %>%
  mutate(
    sig = case_when(
      P.Value < 0.05 & logFC >  0.5 ~ "Up",
      P.Value < 0.05 & logFC < -0.5 ~ "Down",
      TRUE                          ~ "NS"
    )
  )

label_df <- volcano_df %>%
  dplyr::filter(sig != "NS") %>%
  group_by(Tissue) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 5) %>%
  ungroup()

p_vol <- ggplot(volcano_df, aes(logFC, -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(Up = "#D73027", Down = "#4575B4", NS = "grey70")) +
  facet_grid(. ~ Tissue, scales = "free_x") +
  geom_text_repel(
    data            = label_df,
    aes(label       = Metabolite),
    size            = 3,
    max.overlaps    = 50,
    box.padding     = 0.5,
    point.padding   = 0.4,
    segment.color   = "grey30",
    force           = 2,
    force_pull      = 1,
    min.segment.length = 0
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot of Differential Metabolites (DB vs WT)",
    x     = "log2 Fold Change (DB / WT)",
    y     = "-log10(P value)"
  ) +
  theme(
    strip.text        = element_text(face = "bold", size = 14),
    legend.title      = element_blank(),
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )

ggsave("plots/Volcano_Facet_raw.png", p_vol, width = 12, height = 4, dpi = 400)
pdf("plots/Volcano_Facet_raw.pdf", 12, 4); print(p_vol); dev.off()

message("limma/UpSet/Faceted Volcano exported to plots/")

# ------------------ 10. Core-metabolite highlight ----------

core_metabolites <- c(
  "Cholic acid",
  "Taurodeoxycholic acid",
  "Succinic acid",
  "Pyruvic acid"
)

label_df_core <- volcano_df %>%
  dplyr::filter(Metabolite %in% core_metabolites, sig != "NS")

p_vol_core <- ggplot(volcano_df, aes(logFC, -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(Up = "#D73027", Down = "#4575B4", NS = "grey70")) +
  facet_grid(. ~ Tissue, scales = "free_x") +
  geom_text_repel(
    data            = label_df_core,
    aes(label       = Metabolite),
    size            = 5,
    max.overlaps    = 50,
    box.padding     = 0.5,
    point.padding   = 0.4,
    segment.color   = "grey30",
    force           = 2,
    force_pull      = 1,
    min.segment.length = 0
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot of Differential Metabolites (DB vs WT)",
    x     = "log2 Fold Change (DB / WT)",
    y     = "-log10(P value)"
  ) +
  theme(
    strip.text        = element_text(face = "bold", size = 14),
    legend.title      = element_blank(),
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )

ggsave("plots/Volcano_Facet.png",  p_vol_core, width = 12, height = 5, dpi = 400)
pdf("plots/Volcano_Facet.pdf", 12, 5); print(p_vol_core); dev.off()
