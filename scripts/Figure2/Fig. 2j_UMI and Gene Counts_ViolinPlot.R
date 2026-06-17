############################################################
# Fig. 2j 
############################################################

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(scales)

# ============================================================
# Input and output paths
# ============================================================

work_dir <- "Figure2"
plot_dir <- "Fig2j_smClassify_output_summary"

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

sce_rds <- file.path(work_dir, "sce_bac_annotated.rds")

# ============================================================
# Load species-assigned Seurat object
# ============================================================

sce_bac_annot <- readRDS(sce_rds)


# ============================================================
# Plotting theme and colors
# ============================================================

theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      axis.line   = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.4, colour = "black"),
      axis.title  = element_text(size = base_size),
      axis.text   = element_text(size = base_size - 1, colour = "black"),
      plot.title  = element_text(face = "bold", hjust = 0, size = base_size + 1),
      legend.title= element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      plot.margin = margin(6, 6, 6, 6)
    )
}

pal_pheno <- c(
  WT = "#4DBBD5",
  DB = "#E64B35"
)


meta_df <- sce_bac_annot@meta.data %>%
  as.data.frame()

nonr_umi_col <- dplyr::case_when(
  "nonr_total" %in% colnames(meta_df) ~ "nonr_total",
  "nonr_umi"   %in% colnames(meta_df) ~ "nonr_umi",
  TRUE ~ NA_character_
)


sample_order <- c(
  "Cecum-WT-1",
  "Cecum-WT-2",
  "Cecum-DB-1",
  "Cecum-DB-2",
  "Colon-WT-1",
  "Colon-WT-2",
  "Colon-DB-1",
  "Colon-DB-2",
  "Rectum-WT-1",
  "Rectum-WT-2",
  "Rectum-DB-1",
  "Rectum-DB-2"
)

fig2j_qc <- meta_df %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(
    orig.ident = as.character(orig.ident),
    loc = if ("loc" %in% colnames(.)) as.character(loc) else "(all)",
    nonr_umi = as.numeric(.data[[nonr_umi_col]]),
    nonr_features = as.numeric(nonr_features),
    pheno = ifelse(grepl("WT", orig.ident), "WT", "DB")
  ) %>%
  mutate(
    orig.ident = factor(orig.ident, levels = sample_order),
    pheno = factor(pheno, levels = c("WT", "DB"))
  )

fig2j_qc <- fig2j_qc %>%
  dplyr::filter(!is.na(orig.ident))

fwrite(
  as.data.table(fig2j_qc),
  file = file.path(plot_dir, "Fig2j_nonrRNA_recovery_source_data.csv"),
  bom = TRUE
)

# ============================================================
# Fig. 2j: non-rRNA UMI counts per cell
# ============================================================

p_fig2j_nonr_umi <- ggplot(
  fig2j_qc,
  aes(x = orig.ident, y = nonr_umi, fill = pheno)
) +
  geom_violin(colour = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  scale_y_log10() +
  scale_fill_manual(values = pal_pheno) +
  labs(
    title = "smClassify: non-rRNA UMI counts per cell",
    y = "Non-rRNA UMI (log10)",
    x = NULL
  ) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_fig2j_nonr_umi

# ============================================================
# Fig. 2j: detected non-rRNA genes per cell
# ============================================================

p_fig2j_nonr_gene <- ggplot(
  fig2j_qc,
  aes(x = orig.ident, y = nonr_features, fill = pheno)
) +
  geom_violin(colour = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  scale_y_log10() +
  scale_fill_manual(values = pal_pheno) +
  labs(
    title = "smClassify: detected non-rRNA genes per cell",
    y = "Detected non-rRNA genes (log10)",
    x = NULL
  ) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_fig2j_nonr_gene

# ============================================================
# Export Fig. 2j panels
# ============================================================

ggsave(
  filename = file.path(plot_dir, "Fig2j_nonrRNA_UMI_counts_per_cell.png"),
  plot = p_fig2j_nonr_umi,
  width = 5,
  height = 3,
  dpi = 600
)

ggsave(
  filename = file.path(plot_dir, "Fig2j_nonrRNA_UMI_counts_per_cell.pdf"),
  plot = p_fig2j_nonr_umi,
  width = 5,
  height = 3
)

ggsave(
  filename = file.path(plot_dir, "Fig2j_detected_nonrRNA_genes_per_cell.png"),
  plot = p_fig2j_nonr_gene,
  width = 5,
  height = 3,
  dpi = 600
)

ggsave(
  filename = file.path(plot_dir, "Fig2j_detected_nonrRNA_genes_per_cell.pdf"),
  plot = p_fig2j_nonr_gene,
  width = 5,
  height = 3
)

pdf(
  file.path(plot_dir, "Fig2j_nonrRNA_UMI_and_detected_genes_per_cell.pdf"),
  width = 5,
  height = 3
)
print(p_fig2j_nonr_umi)
print(p_fig2j_nonr_gene)
dev.off()
