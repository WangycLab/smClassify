# ============================================================
# Figure 3e. Weighted KEGG pathway score heatmap
library(Seurat)
library(dplyr)
library(tidyr)
library(emmeans)
library(purrr)
library(data.table)
library(ggplot2)
library(pheatmap)

# ============================================================
# 1. Input and output paths
# ============================================================

base_dir   <- "Figure3"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

fig3b_dir <- file.path(output_dir, "Figure3b_KEGG_enrichment")
fig3e_dir <- file.path(output_dir, "Figure3e_WeightedPathwayHeatmap")

dir.create(fig3e_dir, showWarnings = FALSE, recursive = TRUE)

rds_file <- file.path(output_dir, "Figure3d_functional_cluster_composition_input_annotated.rds")
kegg_sig_file <- file.path(fig3b_dir, "Figure3b_kegg_sig.csv")

# ============================================================
# 2. Load data
# ============================================================

seu <- readRDS(rds_file)

samples <- unique(seu$group)

sample_meta <- seu@meta.data %>%
  dplyr::select(group, phenotype, loc) %>%
  dplyr::distinct() %>%
  as.data.frame()
row.names(sample_meta) <- sample_meta$group
sample_meta$group <- factor(sample_meta$group)

seurat_clusters <- unique(seu$seurat_clusters)

kegg_sig <- read.csv(kegg_sig_file, check.names = FALSE)

enrichment_results_filtered <- kegg_sig

enrichment_results_filtered$cluster <- as.character(enrichment_results_filtered$cluster)
enrichment_results_filtered$RichFactor <- as.numeric(enrichment_results_filtered$RichFactor)

seu$seurat_clusters <- as.character(seu$seurat_clusters)
seu$group <- as.character(seu$group)

# ============================================================
# 3. Initial matrix
# ============================================================

sample_pathway_matrix <- matrix(
  NA,
  nrow = length(samples),
  ncol = length(unique(enrichment_results_filtered$Description))
)

rownames(sample_pathway_matrix) <- samples
colnames(sample_pathway_matrix) <- unique(enrichment_results_filtered$Description)

write.csv(
  sample_pathway_matrix,
  file.path(fig3e_dir, "Figure3e_initial_sample_pathway_matrix.csv")
)

# ============================================================
# 4. Cell proportions
# ============================================================

cell_proportions_dt <- seu@meta.data %>%
  dplyr::mutate(
    seurat_clusters = as.character(seurat_clusters),
    group = as.character(group)
  ) %>%
  group_by(seurat_clusters, group) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup()

write.csv(
  cell_proportions_dt,
  file.path(fig3e_dir, "Figure3e_cluster_group_cell_proportions.csv"),
  row.names = FALSE
)

# ============================================================
# 5. Weighted score calculation
# ============================================================

all_pw <- unique(enrichment_results_filtered$Description)

all_weighted_scores <- list()

for (pw in all_pw) {

  pathway_scores <- enrichment_results_filtered[
    enrichment_results_filtered$Description == pw, ]

  weighted_scores <- rep(0, length(samples))
  names(weighted_scores) <- samples

  for (cluster in unique(pathway_scores$cluster)) {

    cluster_score <- pathway_scores$RichFactor[
      pathway_scores$cluster == cluster][1]

    cell_proportion <- cell_proportions_dt[
      cell_proportions_dt$seurat_clusters == cluster, ]

    for (i in seq_along(samples)) {

      group <- samples[i]

      group_proportion <- cell_proportion$proportion[
        cell_proportion$group == group]

      if (length(group_proportion) == 1 &&
          !is.na(group_proportion) &&
          group_proportion > 0) {

        weighted_scores[i] <- weighted_scores[i] +
          (cluster_score * group_proportion)
      }
    }
  }

  all_weighted_scores[[pw]] <- weighted_scores
}

# ============================================================
# 6. Matrix + normalization
# ============================================================

all_weighted_scores_df <- as.data.frame(t(do.call(rbind, all_weighted_scores)))
colnames(all_weighted_scores_df) <- names(all_weighted_scores)

all_weighted_scores_df <- all_weighted_scores_df[as.character(samples), , drop = FALSE]

write.csv(
  all_weighted_scores_df,
  file.path(fig3e_dir, "Figure3e_weighted_scores_original_structure_group_by_pathway.csv")
)

normalized_data <- scale(all_weighted_scores_df)

write.csv(
  normalized_data,
  file.path(fig3e_dir, "Figure3e_weighted_scores_normalized_data.csv")
)

write.csv(
  normalized_data,
  file.path(fig3e_dir, "Heatmap of Weighted Scores-normalized_data.csv")
)

# ============================================================
# 7. Heatmap
# ============================================================

pdf(file.path(fig3e_dir, "heatmap_margin_demo.pdf"), width = 6, height = 7)
pheatmap(
  t(normalized_data),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  main = "Heatmap of Weighted Scores",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

png(file.path(fig3e_dir, "heatmap_margin_demo.png"), width = 1200, height = 1400, res = 200)
pheatmap(
  t(normalized_data),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  main = "Heatmap of Weighted Scores",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

pdf(file.path(fig3e_dir, "Figure3e_weighted_pathway_score_heatmap.pdf"), width = 6, height = 7)
pheatmap(
  t(normalized_data),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  main = "Heatmap of Weighted Scores",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

png(file.path(fig3e_dir, "Figure3e_weighted_pathway_score_heatmap.png"), width = 1200, height = 1400, res = 200)
pheatmap(
  t(normalized_data),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  main = "Heatmap of Weighted Scores",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

