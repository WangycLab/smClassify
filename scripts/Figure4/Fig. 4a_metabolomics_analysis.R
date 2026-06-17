############################################################
## Figure 4 metabolomics analysis and visualisation

# ============================================================
# 0. Load packages
# ============================================================

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
})

# ============================================================
# 1. Input files and output directory
# ============================================================
base_dir   <- "Figure4"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

serum_file   <- file.path(input_dir, "serum_norm_ms2_metabolites_intensity.xlsx")
content_file <- file.path(input_dir, "content_norm_ms2_metabolites_intensity.xlsx")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. Shared palettes, symbols and plotting theme
# ============================================================
pal_pheno <- c(
  WT = "#4DBBD5",
  DB = "#E64B35"
)

pal_region <- c(
  Cecum  = "#1f77b4",
  Colon  = "#2ca02c",
  Rectum = "#ff7f0e",
  Serum  = "#9467bd"
)

pal_volcano <- c(
  Up = "#E64B35",
  Down = "#4DBBD5",
  NS = "grey70"
)

shape_region <- c(
  Cecum  = 16,
  Colon  = 17,
  Rectum = 15,
  Serum  = 18
)

shape_pheno <- c(
  WT = 16,
  DB = 17
)

theme_nature <- function(base_size = 14){
  theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# 3. Read metabolomics tables and combine serum/content data
# ============================================================
serum   <- read_xlsx(serum_file)
content <- read_xlsx(content_file)

names(serum) <- str_replace_all(names(serum), "^Seurm", "Serum")

feature_cols <- c(
  "Type","ID","MZ","RT","MS2.name","MS2_score",
  "Formula","SuperClass","Class","Subclass",
  "HMDB","kegg","kegg_pathway","MS1.name"
)

serum_long <- serum %>%
  pivot_longer(-all_of(feature_cols), names_to = "Sample", values_to = "Intensity")

content_long <- content %>%
  pivot_longer(-all_of(feature_cols), names_to = "Sample", values_to = "Intensity")

all_long <- bind_rows(serum_long, content_long)

# ============================================================
# 4. Parse sample metadata from sample names
# ============================================================

meta <- all_long %>%
  distinct(Sample) %>%
  mutate(
    Tissue_raw = str_extract(Sample, "^[A-Za-z]+"),
    Genotype   = str_extract(Sample, "(WT|DB)"),
    Rep        = as.integer(str_extract(Sample, "\\d{1,2}$")),
    Tissue = Tissue_raw
  ) %>%
  filter(!str_detect(Sample, "^QC")) %>%
  filter(!is.na(Tissue), !is.na(Genotype)) %>%
  mutate(
    Genotype = factor(Genotype, levels = c("WT", "DB")),
    Tissue   = factor(Tissue, levels = c("Cecum", "Colon", "Rectum", "Serum"))
  )

all_long <- all_long %>%
  filter(Sample %in% meta$Sample)

# Export parsed metadata for checking.
write_csv(meta, file.path(output_dir, "sample_metadata_parsed.csv"))

# ============================================================
# 5. Build metabolite intensity matrix
# ============================================================
expr_mat <- all_long %>%
  filter(!is.na(MS2.name)) %>%
  group_by(MS2.name, Sample) %>%
  summarise(Intensity = max(Intensity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Intensity) %>%
  distinct(MS2.name, .keep_all = TRUE) %>%
  column_to_rownames("MS2.name") %>%
  as.matrix()

write.csv(expr_mat, file.path(output_dir, "expr_mat.csv"))

expr_log <- log2(expr_mat + 1)

# ============================================================
# 6. Helper functions
# ============================================================
impute_median <- function(mat){
  t(apply(mat, 1, function(x){
    if (all(is.na(x))) return(x)
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  }))
}

clean_matrix <- function(mat){
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  mat <- impute_median(mat)
  mat <- mat[apply(mat, 1, sd) > 0, , drop = FALSE]
  mat
}

save_plot_pdf_png <- function(plot_obj, prefix, width = 6, height = 5, dpi = 300){
  ggsave(file.path(output_dir, glue("{prefix}.png")), plot_obj, width = width, height = height, dpi = dpi)
  pdf(file.path(output_dir, glue("{prefix}.pdf")), width = width, height = height)
  print(plot_obj)
  dev.off()
}


# ============================================================
# 7. PCA of gut-content metabolite profiles
# ============================================================
content_samples <- meta %>%
  filter(Tissue %in% c("Colon", "Cecum", "Rectum")) %>%
  mutate(Group = paste(Tissue, Genotype, sep = "_"))

mat_raw <- expr_log[, content_samples$Sample, drop = FALSE]
mat_raw <- clean_matrix(mat_raw)

pca_res <- prcomp(t(mat_raw), scale. = TRUE)

pca_scores <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("Sample") %>%
  left_join(content_samples, by = "Sample")

pc1_var <- round(100 * summary(pca_res)$importance[2, 1], 1)
pc2_var <- round(100 * summary(pca_res)$importance[2, 2], 1)

# 7.1 PCA coloured by phenotype.
p_pca_pheno <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = Genotype, shape = Tissue)
) +
  geom_point(size = 3.2) +
  stat_ellipse(
    aes(group = Genotype, fill = Genotype),
    geom = "polygon", alpha = 0.16, color = NA
  ) +
  stat_ellipse(
    aes(group = Genotype),
    linetype = 2, color = "grey40", linewidth = 0.5
  ) +
  scale_color_manual(values = pal_pheno) +
  scale_fill_manual(values = pal_pheno) +
  scale_shape_manual(values = shape_region[c("Cecum", "Colon", "Rectum")]) +
  labs(
    title = "PCA â€? Gut Content (colored by phenotype)",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_nature(14)

save_plot_pdf_png(p_pca_pheno, "PCA_GutContent_byPhenotype", width = 6, height = 5, dpi = 300)

# 7.2 PCA coloured by gut region.
p_pca_region <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = Tissue, shape = Genotype)
) +
  geom_point(size = 3.2) +
  stat_ellipse(
    aes(group = Tissue, fill = Tissue),
    geom = "polygon", alpha = 0.16, color = NA
  ) +
  stat_ellipse(
    aes(group = Tissue),
    linetype = 2, color = "grey40", linewidth = 0.5
  ) +
  scale_color_manual(values = pal_region[c("Cecum", "Colon", "Rectum")]) +
  scale_fill_manual(values = pal_region[c("Cecum", "Colon", "Rectum")]) +
  scale_shape_manual(values = shape_pheno) +
  labs(
    title = "PCA â€? Gut Content (colored by region)",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_nature(14)

save_plot_pdf_png(p_pca_region, "PCA_GutContent_byRegion", width = 6, height = 5, dpi = 300)

# 7.3 PCA combined six-group view.
p_pca_6g <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = Genotype, shape = Tissue)
) +
  geom_point(size = 3.2) +
  stat_ellipse(
    aes(group = interaction(Tissue, Genotype), fill = Genotype),
    geom = "polygon", alpha = 0.14, color = NA
  ) +
  stat_ellipse(
    aes(group = interaction(Tissue, Genotype)),
    linetype = 2, color = "grey40", linewidth = 0.5
  ) +
  scale_color_manual(values = pal_pheno) +
  scale_fill_manual(values = pal_pheno) +
  scale_shape_manual(values = shape_region[c("Cecum", "Colon", "Rectum")]) +
  labs(
    title = "PCA â€? Gut Content (Cecum / Colon / Rectum Ã— WT / DB)",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_nature(14)

save_plot_pdf_png(p_pca_6g, "PCA_GutContent_6groups", width = 7, height = 6, dpi = 300)
message("PCA plots saved.")

############################################################
## 8. PLS-DA for six gut-content groups
##    Colon/Cecum/Rectum x WT/DB
############################################################
content_meta <- meta %>%
  filter(Tissue %in% c("Colon", "Cecum", "Rectum")) %>%
  mutate(Group = paste(Tissue, Genotype, sep = "_"))

samples <- content_meta$Sample
expr_sub <- expr_log[, samples, drop = FALSE]
expr_sub <- clean_matrix(expr_sub)

X <- t(expr_sub)
Y <- content_meta$Group

set.seed(123)
if (any(table(Y) < 2)) {
  warning("PLS-DA (6 groups): At least one class has <2 samples; model may not run.")
}

plsda_mod <- opls(X, Y, predI = 2, orthoI = 0, permI = 100)

scores_6g <- data.frame(
  Sample = rownames(plsda_mod@scoreMN),
  PLS1 = plsda_mod@scoreMN[, 1],
  PLS2 = plsda_mod@scoreMN[, 2]
) %>%
  left_join(content_meta, by = "Sample")

p_pls_6g <- ggplot(
  scores_6g,
  aes(x = PLS1, y = PLS2, color = Genotype, shape = Tissue)
) +
  geom_point(size = 3.2) +
  stat_ellipse(
    aes(group = interaction(Tissue, Genotype), fill = Genotype),
    geom = "polygon", alpha = 0.18, color = NA
  ) +
  stat_ellipse(
    aes(group = interaction(Tissue, Genotype)),
    linetype = 2, color = "grey40", linewidth = 0.5
  ) +
  scale_color_manual(values = pal_pheno) +
  scale_fill_manual(values = pal_pheno) +
  scale_shape_manual(values = shape_region[c("Cecum", "Colon", "Rectum")]) +
  labs(
    title = "PLS-DA:Cecum/Colon/Rectum Ã— WT/DB",
    x = "PLS-DA Component 1",
    y = "PLS-DA Component 2"
  ) +
  theme_nature(14)

save_plot_pdf_png(p_pls_6g, "PLSDA_Content_6groups", width = 6, height = 5, dpi = 300)

############################################################
## 9. PLS-DA plotting helper
############################################################
plot_plsda <- function(score_df, out_prefix, title, mode = c("pheno", "region")){
  mode <- match.arg(mode)
  
  if (mode == "pheno") {
    p <- ggplot(score_df, aes(x = PLS1, y = PLS2, color = Group, shape = Group)) +
      geom_point(size = 3.2) +
      stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.18, color = NA) +
      stat_ellipse(linetype = 2, color = "grey40", linewidth = 0.5) +
      scale_color_manual(values = pal_pheno) +
      scale_fill_manual(values = pal_pheno)
  } else {
    use_cols <- pal_region[names(pal_region) %in% unique(as.character(score_df$Group))]
    p <- ggplot(score_df, aes(x = PLS1, y = PLS2, color = Group, shape = Group)) +
      geom_point(size = 3.2) +
      stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.18, color = NA) +
      stat_ellipse(linetype = 2, color = "grey40", linewidth = 0.5) +
      scale_color_manual(values = use_cols) +
      scale_fill_manual(values = use_cols)
  }
  
  p <- p +
    labs(title = title, x = "PLS Component 1", y = "PLS Component 2") +
    theme_nature(14)
  
  save_plot_pdf_png(p, paste0(out_prefix, "_score"), width = 6, height = 5, dpi = 500)
}

############################################################
## 10. PLS-DA: WT vs DB across all gut segments
############################################################
meta_sub1 <- meta %>%
  filter(Tissue %in% c("Colon", "Cecum", "Rectum")) %>%
  mutate(Group = Genotype)

samples1 <- intersect(meta_sub1$Sample, colnames(expr_log))
mat1 <- expr_log[, samples1, drop = FALSE]
mat1 <- clean_matrix(mat1)

meta_sub1_aligned <- meta_sub1 %>%
  filter(Sample %in% colnames(mat1)) %>%
  mutate(Sample = factor(Sample, levels = colnames(mat1))) %>%
  arrange(Sample)

Y1 <- droplevels(factor(meta_sub1_aligned$Group, levels = c("WT", "DB")))
X1 <- t(mat1[, as.character(meta_sub1_aligned$Sample), drop = FALSE])

set.seed(123)
if (any(table(Y1) < 2)) warning("PLS-DA (WT vs DB): some classes have <2 samples.")

pls1 <- opls(X1, Y1, predI = 2, orthoI = 0, permI = 100)

score1 <- data.frame(
  PLS1 = pls1@scoreMN[, 1],
  PLS2 = pls1@scoreMN[, 2],
  Group = factor(Y1, levels = c("WT", "DB"))
)

plot_plsda(score1, out_prefix = "PLSDA_WTvsDB", title = "PLS-DA: WT vs DB", mode = "pheno")

############################################################
## 11. PLS-DA: Colon vs Cecum vs Rectum in WT samples only
############################################################
meta_sub2 <- meta %>%
  filter(Tissue %in% c("Colon", "Cecum", "Rectum") & Genotype == "WT") %>%
  mutate(Group = Tissue)

samples2 <- intersect(meta_sub2$Sample, colnames(expr_log))
mat2 <- expr_log[, samples2, drop = FALSE]
mat2 <- clean_matrix(mat2)

meta_sub2_aligned <- meta_sub2 %>%
  filter(Sample %in% colnames(mat2)) %>%
  mutate(Sample = factor(Sample, levels = colnames(mat2))) %>%
  arrange(Sample)

Y2 <- droplevels(factor(meta_sub2_aligned$Group, levels = c("Cecum", "Colon", "Rectum")))
print(table(Y2))

keep_groups2 <- names(table(Y2))[table(Y2) >= 2]
meta_sub2_aligned <- meta_sub2_aligned %>%
  filter(Group %in% keep_groups2)

if (nrow(meta_sub2_aligned) == 0 || length(unique(meta_sub2_aligned$Group)) < 2) {
  warning("PLS-DA (WT tissues) skipped: fewer than 2 groups with >=2 samples.")
} else {
  mat2 <- mat2[, as.character(meta_sub2_aligned$Sample), drop = FALSE]
  Y2 <- droplevels(factor(meta_sub2_aligned$Group, levels = keep_groups2))
  X2 <- t(mat2)
  
  print(table(Y2))
  
  set.seed(123)
  pls2 <- opls(X2, Y2, predI = 2, orthoI = 0, permI = 100)
  
  score2 <- data.frame(
    PLS1 = pls2@scoreMN[, 1],
    PLS2 = pls2@scoreMN[, 2],
    Group = factor(Y2, levels = keep_groups2)
  )
  
  plot_plsda(
    score2,
    out_prefix = "PLSDA_WT_tissues",
    title = "PLS-DA: Colon vs Cecum vs Rectum (WT only)",
    mode = "region"
  )
}

############################################################
## 12. Differential metabolite analysis: DB vs WT per tissue
############################################################
tissues <- c("Serum", "Colon", "Cecum", "Rectum")
diff_list <- list()

for (tis in tissues) {
  cat("Processing:", tis, "\n")
  
  samp <- meta %>%
    filter(Tissue == tis) %>%
    pull(Sample)
  
  mat <- expr_log[, samp, drop = FALSE]
  mat <- clean_matrix(mat)
  
  geno <- factor(meta$Genotype[match(colnames(mat), meta$Sample)], levels = c("WT", "DB"))
  design <- model.matrix(~ 0 + geno)
  colnames(design) <- c("WT", "DB")
  
  fit <- lmFit(mat, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(DB - WT, levels = design)))
  tab <- topTable(fit2, number = Inf) %>%
    rownames_to_column("Metabolite")
  
  write_csv(tab, file.path(output_dir, glue("{tis}_DBvsWT_limma.csv")))
  
  sig <- tab %>%
    filter(P.Value < 0.05, abs(logFC) > 0.5) %>%
    pull(Metabolite)
  
  diff_list[[tis]] <- sig
  cat(" -> Significant:", length(sig), "metabolites\n")
}

# ============================================================
# 13. UpSet plot of significant differential metabolites
# ============================================================
up_in <- fromList(diff_list)

png(file.path(output_dir, "DiffMetabolites_UpSet.png"), width = 2600, height = 2200, res = 350)
UpSetR::upset(
  up_in,
  sets = tissues,
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = pal_pheno["WT"],
  main.bar.color = pal_pheno["DB"],
  text.scale = 2
)
dev.off()

pdf(file.path(output_dir, "DiffMetabolites_UpSet.pdf"), width = 6, height = 5)
UpSetR::upset(
  up_in,
  sets = tissues,
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = pal_pheno["WT"],
  main.bar.color = pal_pheno["DB"],
  text.scale = 2
)
dev.off()

# ============================================================
# 14. Faceted volcano plot across tissues
# ============================================================
volcano_df <- bind_rows(lapply(tissues, function(tis) {
  read_csv(file.path(output_dir, glue("{tis}_DBvsWT_limma.csv")), show_col_types = FALSE) %>%
    mutate(Tissue = tis)
})) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Serum", "Colon", "Cecum", "Rectum")),
    sig = case_when(
      P.Value < 0.05 & logFC >  0.5 ~ "Up",
      P.Value < 0.05 & logFC < -0.5 ~ "Down",
      TRUE ~ "NS"
    )
  )


label_df <- volcano_df %>%
  filter(sig != "NS") %>%
  group_by(Tissue) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 5) %>%
  ungroup()

p_vol <- ggplot(volcano_df, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = pal_volcano) +
  facet_grid(. ~ Tissue, scales = "free_x") +
  geom_text_repel(
    data = label_df,
    aes(label = Metabolite),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.color = "grey30",
    force = 2,
    force_pull = 1,
    min.segment.length = 0
  ) +
  labs(
    title = "Volcano Plot of Differential Metabolites (DB vs WT)",
    x = "log2 Fold Change (DB / WT)",
    y = "-log10(P value)"
  ) +
  theme_nature(14) +
  theme(strip.text = element_text(face = "bold", size = 14))

save_plot_pdf_png(p_vol, "Volcano_Facet_all_labels", width = 12, height = 4, dpi = 400)

# ============================================================
# 15. Faceted volcano plot with selected core metabolite labels
# ============================================================
core_metabolites <- c(
  "Cholic acid",
  "Taurodeoxycholic acid",
  "Succinic acid",
  "Pyruvic acid",
  "Leucine",
  "Isoleucine"
)

label_df_core <- volcano_df %>%
  filter(Metabolite %in% core_metabolites, sig != "NS")

p_vol_core <- ggplot(volcano_df, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = pal_volcano) +
  facet_grid(. ~ Tissue, scales = "free_x") +
  geom_text_repel(
    data = label_df_core,
    aes(label = Metabolite),
    size = 5,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.color = "grey30",
    force = 2,
    force_pull = 1,
    min.segment.length = 0
  ) +
  labs(
    title = "Volcano Plot of Differential Metabolites (DB vs WT)",
    x = "log2 Fold Change (DB / WT)",
    y = "-log10(P value)"
  ) +
  theme_nature(14) +
  theme(strip.text = element_text(face = "bold", size = 14))

save_plot_pdf_png(p_vol_core, "Volcano_Facet_core_labels", width = 12, height = 5, dpi = 400)

message("All PCA / PLS-DA / limma / UpSet / volcano plots exported to: ", output_dir)


# ============================================================
# 16. Source data for Fig. 4a volcano plots
# ============================================================
# Export full source data underlying Fig. 4a, including adjusted P values.

fig4a_source_data <- volcano_df %>%
  mutate(
    Direction = case_when(
      sig == "Up"   ~ "Up in DB",
      sig == "Down" ~ "Down in DB",
      TRUE          ~ "Not significant"
    ),
    Is_labelled_core_metabolite = Metabolite %in% core_metabolites,
    log2FC_cutoff = 0.5,
    P_value_cutoff = 0.05,
    Significant_for_volcano = P.Value < 0.05 & abs(logFC) > 0.5
  ) %>%
  select(
    Tissue,
    Metabolite,
    log2_fold_change_DB_vs_WT = logFC,
    average_log2_intensity = AveExpr,
    moderated_t_statistic = t,
    P_value = P.Value,
    adjusted_P_value_BH = adj.P.Val,
    B_statistic = B,
    Direction,
    Significant_for_volcano,
    Is_labelled_core_metabolite,
    log2FC_cutoff,
    P_value_cutoff
  ) %>%
  arrange(Tissue, P_value)

write_csv(
  fig4a_source_data,
  file.path(output_dir, "Fig4a_source_data_volcano_metabolites_DB_vs_WT.csv")
)
