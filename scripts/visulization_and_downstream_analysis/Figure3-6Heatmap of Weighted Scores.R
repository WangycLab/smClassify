library(dplyr)
library(tidyr)
library(emmeans)
library(purrr)
library(data.table)
library(ggplot2)
library(pheatmap)

setwd("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009")
load("seu_symbol_SCT_clustered_annotated.RData")

# Get metadata from Seurat object
sample_meta <- seu@meta.data %>%
  dplyr::select(group, phenotype, loc) %>%
  dplyr::distinct() %>%
  as.data.frame()
row.names(sample_meta) <- sample_meta$group
sample_meta$group <- factor(sample_meta$group)
seurat_clusters <- unique(seu$seurat_clusters)

# Get KEGG Result
kegg_sig <- read.csv("kegg_sig.csv")
enrichment_results_filtered <- kegg_sig

# Create a matrix to store the weighted enrichment scores of each group in each cluster
samples <- unique(seu$group)
samples
sample_pathway_matrix <- matrix(
  NA,
  nrow = length(samples),
  ncol = length(unique(enrichment_results_filtered$Description))
)
rownames(sample_pathway_matrix) <- samples
colnames(sample_pathway_matrix) <- unique(enrichment_results_filtered$Description)
sample_pathway_matrix

# Calculate ratio of each cell type
cell_proportions_dt <- seu@meta.data %>%
  group_by(seurat_clusters, group) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup()

print(cell_proportions_dt)

# Get all pathways
all_pw <- unique(enrichment_results_filtered$Description)

# Initialize weighted score
all_weighted_scores <- list()

# Loop
for (pw in all_pw) {
  print(paste("Calculating pathway:", pw))
  
  # Get scores of all clusters
  pathway_scores <- enrichment_results_filtered[enrichment_results_filtered$Description == pw, ]
  print(paste("The richfactor of this pathway is:", paste(pathway_scores$cluster, collapse = ", ")))
  
  # Initialize weighted score in each cluster
  weighted_scores <- rep(0, length(samples)) 
  names(weighted_scores) <- samples  
  
  # Loop through each cluster
  for (cluster in unique(pathway_scores$cluster)) {
    print(paste("Processing cluster:", cluster))  
    cluster_score <- pathway_scores$RichFactor[pathway_scores$cluster == cluster]
    
    print(paste("Cluster", cluster, "Richfactor:", cluster_score))  
    
    if (is.na(cluster_score)) {
      print(paste("Cluster", cluster, "No richfactor available, skipping"))
      next  
    }
    
    # Get the cell proportion of each group within this cluster
    cell_proportion <- cell_proportions_dt[cell_proportions_dt$seurat_clusters == cluster, ]
    
    print(paste("cluster", cluster, "cell proportion", paste(cell_proportion$proportion, collapse = ", ")))  
    # Loop through each group and calculate the weighted enrichment score.
    for (i in 1:length(samples)) {
      group <- samples[i]
      group_proportion <- cell_proportion$proportion[cell_proportion$group == group]
      print(paste("Group", group, "cell proportion:", group_proportion))  
      
      if (!is.na(group_proportion) && group_proportion > 0) {
        weighted_scores[i] <- weighted_scores[i] + (cluster_score * group_proportion)
      }
    }
    
    print(paste("cluster", cluster, "weighted score after clustering", paste(weighted_scores, collapse = ", ")))  
  }
  
  print(paste("Pathway weighted scores:", paste(weighted_scores, collapse = ", ")))
  
  all_weighted_scores[[pw]] <- weighted_scores
}

all_weighted_scores_df <- as.data.frame(t(do.call(rbind, all_weighted_scores)))
colnames(all_weighted_scores_df) <- names(all_weighted_scores)
head(all_weighted_scores_df)

# Z-score Normalization
normalized_data <- scale(all_weighted_scores_df)


# Output png and pdf files
png("heatmap_margin_demo.png", width = 1200, height = 1400, res = 200)

pheatmap(t(normalized_data), 
         cluster_rows = T,   
         cluster_cols = F,    
         scale = "column",          
         main = "Heatmap of Weighted Scores",
         color = colorRampPalette(c("blue", "white", "red"))(50)
)

dev.off()


pdf("heatmap_margin_demo.pdf", width = 6, height = 7)
pheatmap(t(normalized_data), 
         cluster_rows = T,   
         cluster_cols = F,    
         scale = "column",         
         main = "Heatmap of Weighted Scores",
         color = colorRampPalette(c("blue", "white", "red"))(50)
)

dev.off()

