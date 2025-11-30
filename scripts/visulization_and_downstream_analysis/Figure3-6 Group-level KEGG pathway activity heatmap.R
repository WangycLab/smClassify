# This script was used to generate Figure 3e
#
# Description
#   This script aggregates cluster-level KEGG enrichment (RichFactor)
#   into group-level “weighted pathway scores”, and visualizes them
#   as a heatmap.
#
#   For each group and each KEGG pathway:
#     - Take cluster-wise KEGG RichFactor from kegg_sig.csv.
#     - Weight each cluster's RichFactor by the proportion of cells
#       that group contributes to that cluster.
#     - Sum over clusters → group × pathway weighted score.
#
#   Finally, z-score normalization and heatmap visualization are used
#   to generate Figure 3e.
#
# Inputs (assumed in project root):
#   - seu_symbol_SCT_clustered_annotated.RData
#       Seurat object named `seu` with:
#         * meta.data columns: group, phenotype, loc, seurat_clusters
#   - kegg_sig.csv
#       Cluster-level significant KEGG terms with at least:
#         * Description (pathway name)
#         * cluster (cluster ID)
#         * RichFactor (enrichment strength per cluster)
#
# Outputs:
#   - heatmap_margin_demo.png
#   - heatmap_margin_demo.pdf
#
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(emmeans)
  library(purrr)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# ------------------------------------------------------------
# 0. Project directory and data loading
# ------------------------------------------------------------

project_dir <- "path/to/project_root"  # <-- modify to your own path
setwd(project_dir)

load("seu_symbol_SCT_clustered_annotated.RData")  # loads `seu`

# ------------------------------------------------------------
# 1. Sample metadata and basic objects
# ------------------------------------------------------------

# Extract unique sample-level metadata from Seurat object
sample_meta <- seu@meta.data %>%
  dplyr::select(group, phenotype, loc) %>%
  dplyr::distinct() %>%
  as.data.frame()

row.names(sample_meta) <- sample_meta$group
sample_meta$group <- factor(sample_meta$group)

# For reference (not directly used in the downstream calculations below)
seurat_clusters <- unique(seu$seurat_clusters)

# ------------------------------------------------------------
# 2. Load KEGG enrichment results (cluster-level)
# ------------------------------------------------------------

# kegg_sig.csv should come from your prior KEGG analysis script
kegg_sig <- read.csv("kegg_sig.csv")
enrichment_results_filtered <- kegg_sig   # keep the original object name

# ------------------------------------------------------------
# 3. Initialize group × pathway matrix
# ------------------------------------------------------------

samples <- unique(seu$group)
all_pw  <- unique(enrichment_results_filtered$Description)

sample_pathway_matrix <- matrix(
  NA_real_,
  nrow = length(samples),
  ncol = length(all_pw)
)
rownames(sample_pathway_matrix) <- samples
colnames(sample_pathway_matrix) <- all_pw

# ------------------------------------------------------------
# 4. Cluster-wise cell proportions per sample
# ------------------------------------------------------------

# For each (cluster, group), compute the fraction of cells contributed
# by that group to the given cluster.
cell_proportions_dt <- seu@meta.data %>%
  dplyr::group_by(seurat_clusters, group) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%
  dplyr::ungroup()

print(cell_proportions_dt)

# ------------------------------------------------------------
# 5. Compute weighted KEGG pathway scores per sample
# ------------------------------------------------------------

# all_weighted_scores:
#   a list of named numeric vectors (one vector per pathway),
#   with names corresponding to `samples`.
all_weighted_scores <- list()

for (pw in all_pw) {
  message("Calculating pathway: ", pw)
  
  # Cluster-level entries for this pathway
  pathway_scores <- enrichment_results_filtered[
    enrichment_results_filtered$Description == pw, 
  ]
  message("Clusters with this pathway: ",
          paste(pathway_scores$cluster, collapse = ", "))
  
  # Initialize pathway scores for each sample (group)
  weighted_scores <- rep(0, length(samples))
  names(weighted_scores) <- samples
  
  # Loop over clusters contributing to this pathway
  for (cluster in unique(pathway_scores$cluster)) {
    message("  Processing cluster: ", cluster)
    
    cluster_score <- pathway_scores$RichFactor[pathway_scores$cluster == cluster]
    message("    RichFactor for cluster ", cluster, ": ", paste(cluster_score, collapse = ", "))
    
    if (is.na(cluster_score)) {
      message("    No RichFactor available, skipping this cluster.")
      next
    }
    
    # Cell fraction of each sample within this cluster
    cell_proportion <- cell_proportions_dt[
      cell_proportions_dt$seurat_clusters == cluster, 
    ]
    
    message("    Cluster ", cluster, " cell proportions: ",
            paste(cell_proportion$proportion, collapse = ", "))
    
    # For each group, add weighted contribution
    for (i in seq_along(samples)) {
      group <- samples[i]
      group_proportion <- cell_proportion$proportion[cell_proportion$group == group]
      message("      Group ", group, " proportion: ", group_proportion)
      
      if (!is.na(group_proportion) && group_proportion > 0) {
        weighted_scores[i] <- weighted_scores[i] + (cluster_score * group_proportion)
      }
    }
    
    message("    Weighted scores after cluster ", cluster, ": ",
            paste(weighted_scores, collapse = ", "))
  }
  
  message("Final pathway (", pw, ") weighted scores: ",
          paste(weighted_scores, collapse = ", "))
  
  all_weighted_scores[[pw]] <- weighted_scores
}

# Convert list → sample × pathway matrix (samples as rows)
all_weighted_scores_df <- as.data.frame(
  t(do.call(rbind, all_weighted_scores))
)
colnames(all_weighted_scores_df) <- names(all_weighted_scores)
head(all_weighted_scores_df)

# ------------------------------------------------------------
# 6. Z-score normalization across pathways
# ------------------------------------------------------------

# Normalize scores per pathway (column-wise z-score)
normalized_data <- scale(all_weighted_scores_df)

# ------------------------------------------------------------
# 7. Heatmap visualization (Figure 3e)
# ------------------------------------------------------------

# PNG output
png("heatmap_margin_demo.png", width = 1200, height = 1400, res = 200)
pheatmap(
  t(normalized_data),
  cluster_rows = TRUE,   # cluster pathways
  cluster_cols = FALSE,  # keep sample order as-is
  scale        = "column",
  main         = "Heatmap of Weighted Scores",
  color        = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

# PDF output
pdf("heatmap_margin_demo.pdf", width = 6, height = 7)
pheatmap(
  t(normalized_data),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale        = "column",
  main         = "Heatmap of Weighted Scores",
  color        = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()
