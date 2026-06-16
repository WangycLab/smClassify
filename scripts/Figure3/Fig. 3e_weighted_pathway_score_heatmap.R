# ============================================================
# Figure 3e. Weighted KEGG pathway score heatmap

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(emmeans)
  library(purrr)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# ============================================================
# 1. Input and output paths
# ============================================================

base_dir   <- "Figure3"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

fig3b_dir <- file.path(output_dir, "Figure3b_KEGG_enrichment")
fig3e_dir <- file.path(output_dir, "Figure3e_WeightedPathwayHeatmap")

dir.create(fig3e_dir, showWarnings = FALSE, recursive = TRUE)

rds_file <- file.path(output_dir, "Figure3d_Sankey_Composition/Figure3d_functional_cluster_composition_input_annotated.rds")
kegg_sig_file <- file.path(fig3b_dir, "Figure3b_kegg_sig.csv")



# ============================================================
# 2. Load data
# ============================================================

seu <- readRDS(rds_file)


# Keep the original group order from the Seurat object.
# This is important for reproducing the original heatmap column order.
samples <- unique(seu$group)

# Extract sample metadata as in the original script.
sample_meta <- seu@meta.data %>%
  dplyr::select(group, phenotype, loc) %>%
  dplyr::distinct() %>%
  as.data.frame()
row.names(sample_meta) <- sample_meta$group
sample_meta$group <- factor(sample_meta$group)

# Get all Seurat clusters as in the original script.
seurat_clusters <- unique(seu$seurat_clusters)

# Read KEGG enrichment results.
kegg_sig <- read.csv(kegg_sig_file, check.names = FALSE)


enrichment_results_filtered <- kegg_sig


# Match original type behavior as closely as possible.
enrichment_results_filtered$cluster <- as.character(enrichment_results_filtered$cluster)
enrichment_results_filtered$RichFactor <- as.numeric(enrichment_results_filtered$RichFactor)
seu$seurat_clusters <- as.character(seu$seurat_clusters)
seu$group <- as.character(seu$group)

# ============================================================
# 3. Initialize pathway matrix, as in the original script
# ============================================================

sample_pathway_matrix <- matrix(
  NA,
  nrow = length(samples),
  ncol = length(unique(enrichment_results_filtered$Description))
)
rownames(sample_pathway_matrix) <- samples
colnames(sample_pathway_matrix) <- unique(enrichment_results_filtered$Description)

# Export for source-data checking.
write.csv(
  sample_pathway_matrix,
  file.path(fig3e_dir, "Figure3e_initial_sample_pathway_matrix.csv")
)

# ============================================================
# 4. Calculate cell proportions by Seurat cluster and group
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

print(cell_proportions_dt)

write.csv(
  cell_proportions_dt,
  file.path(fig3e_dir, "Figure3e_cluster_group_cell_proportions.csv"),
  row.names = FALSE
)

# ============================================================
# 5. Original pathway-wise weighted-score calculation
# ============================================================

all_pw <- unique(enrichment_results_filtered$Description)
print(paste("All pathways:", paste(all_pw, collapse = ", ")))

for (pw in all_pw) {
  print(paste("Processing pathway:", pw))
  pathway_scores <- enrichment_results_filtered[enrichment_results_filtered$Description == pw, ]
  print(paste("Clusters with RichFactor for this pathway:", paste(pathway_scores$cluster, collapse = ", ")))
}

print(paste("Group list:", paste(samples, collapse = ", ")))

# Keep the original diagnostic block for the last pathway_scores object.
weighted_scores <- rep(0, length(samples))

for (cluster in unique(pathway_scores$cluster)) {
  print(paste("Processing cluster", cluster))
  cluster_score <- pathway_scores$RichFactor[pathway_scores$cluster == cluster]
  if (length(cluster_score) == 0 || is.na(cluster_score[1])) {
    next
  }
  cell_proportion <- cell_proportions_dt[cell_proportions_dt$seurat_clusters == as.character(cluster), ]
  print(paste("Cell proportions in cluster", cluster, ":"))
  print(cell_proportion)
}

for (cluster in unique(pathway_scores$cluster)) {
  cluster_score <- pathway_scores$RichFactor[pathway_scores$cluster == cluster]
  if (length(cluster_score) == 0 || is.na(cluster_score[1])) {
    next
  }
  # If multiple rows exist for the same pathway-cluster pair, keep the first value.
  # In the original script this was assumed to be length 1.
  cluster_score <- cluster_score[1]
  cell_proportion <- cell_proportions_dt[cell_proportions_dt$seurat_clusters == as.character(cluster), ]

  for (i in seq_along(samples)) {
    group <- samples[i]
    group_proportion <- cell_proportion$proportion[cell_proportion$group == group]
    if (length(group_proportion) == 1 && !is.na(group_proportion) && group_proportion > 0) {
      weighted_scores[i] <- weighted_scores[i] + (cluster_score * group_proportion)
    }
  }
}

print(paste("Weighted scores for this pathway:", paste(weighted_scores, collapse = ", ")))

# Main original calculation.
all_weighted_scores <- list()

for (pw in all_pw) {
  print(paste("Calculating pathway:", pw))

  pathway_scores <- enrichment_results_filtered[enrichment_results_filtered$Description == pw, ]
  print(paste("Clusters with RichFactor for this pathway:", paste(pathway_scores$cluster, collapse = ", ")))

  weighted_scores <- rep(0, length(samples))
  names(weighted_scores) <- samples

  for (cluster in unique(pathway_scores$cluster)) {
    print(paste("Processing cluster", cluster))

    cluster_score <- pathway_scores$RichFactor[pathway_scores$cluster == cluster]
    print(paste("RichFactor for cluster", cluster, ":", paste(cluster_score, collapse = ", ")))

    if (length(cluster_score) == 0 || is.na(cluster_score[1])) {
      print(paste("Cluster", cluster, "has no valid RichFactor; skipped."))
      next
    }

    # Keep the original assumption of one score per pathway-cluster pair.
    # This line prevents vector-length issues while retaining the original result
    # when each pair is unique.
    cluster_score <- cluster_score[1]

    cell_proportion <- cell_proportions_dt[cell_proportions_dt$seurat_clusters == cluster, ]
    print(paste("Cell proportions in cluster", cluster, ":", paste(cell_proportion$proportion, collapse = ", ")))

    for (i in seq_along(samples)) {
      group <- samples[i]
      group_proportion <- cell_proportion$proportion[cell_proportion$group == group]
      print(paste("Group", group, "cell proportion:", paste(group_proportion, collapse = ", ")))

      if (length(group_proportion) == 1 && !is.na(group_proportion) && group_proportion > 0) {
        weighted_scores[i] <- weighted_scores[i] + (cluster_score * group_proportion)
      }
    }

    print(paste("Weighted scores after cluster", cluster, ":", paste(weighted_scores, collapse = ", ")))
  }

  print(paste("Final weighted scores for pathway:", paste(weighted_scores, collapse = ", ")))
  all_weighted_scores[[pw]] <- weighted_scores
}

# ============================================================
# 6. Convert weighted scores to data frame using original structure
# ============================================================
# Original:
#   all_weighted_scores_df <- as.data.frame(t(do.call(rbind, all_weighted_scores)))
#   colnames(all_weighted_scores_df) <- names(all_weighted_scores)
#   normalized_data <- scale(all_weighted_scores_df)
#   pheatmap(t(normalized_data), ...)
# ============================================================

all_weighted_scores_df <- as.data.frame(t(do.call(rbind, all_weighted_scores)))
colnames(all_weighted_scores_df) <- names(all_weighted_scores)

# Keep original row order from samples.
all_weighted_scores_df <- all_weighted_scores_df[as.character(samples), , drop = FALSE]

write.csv(
  all_weighted_scores_df,
  file.path(fig3e_dir, "Figure3e_weighted_scores_original_structure_group_by_pathway.csv")
)

head(all_weighted_scores_df)

# Z-score normalization exactly as in the original script.
normalized_data <- scale(all_weighted_scores_df)

write.csv(
  normalized_data,
  file.path(fig3e_dir, "Figure3e_weighted_scores_normalized_data.csv")
)

# Keep original output name as a compatibility file, but save inside Figure 3e output dir.
write.csv(
  normalized_data,
  file.path(fig3e_dir, "Heatmap of Weighted Scores-normalized_data.csv")
)

# ============================================================
# 7. Plot heatmap exactly as original
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

# Updated Figure 3e output names with identical plotting object and parameters.
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


