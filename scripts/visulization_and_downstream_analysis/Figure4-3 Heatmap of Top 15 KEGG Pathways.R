# =============================================================================
# This script generates Figure 4b
# Input:
#   - cpd_DE_enrich.csv  # KEGG enrichment results from metabolomics (per group)
#
# Output:
#   - cluster_top15_heat_matrix.csv  # matrix of -log10(qvalue) for selected pathways
#   - Top15_KEGG_Pathway_Heatmap.pdf # heatmap of top 15 KEGG pathways per group
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(igraph)
  library(ggraph)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
})

# -----------------------------------------------------------------------------
# 1. Working directory and input
# -----------------------------------------------------------------------------
# Set your working directory to the folder containing:
#   - cpd_DE_enrich.csv
# Example:
#   setwd("/path/to/project/Figure4")
# Here we assume the working directory has already been set correctly.

# Read metabolomics KEGG enrichment table
kegg_meta <- read.csv("cpd_DE_enrich.csv", stringsAsFactors = FALSE)
head(kegg_meta)

# -----------------------------------------------------------------------------
# 2. Build matrix of -log10(qvalue) for top KEGG pathways
# -----------------------------------------------------------------------------
# Step 2.1: For each metabolomics group (Cluster column), keep top 15 pathways
#           by RichFactor (as in original code).
top_paths <- kegg_meta %>%
  mutate(log_q = -log10(qvalue)) %>%
  group_by(Cluster) %>%
  arrange(RichFactor) %>%         # same ordering rule as original
  slice_head(n = 15) %>%
  ungroup()

# Step 2.2: Collect the set of all unique pathway descriptions (row candidates)
top_descriptions <- unique(top_paths$Description)

# Step 2.3: For these pathways, get -log10(qvalue) across all clusters/groups
heat_data <- kegg_meta %>%
  mutate(log_q = -log10(qvalue)) %>%
  filter(Description %in% top_descriptions) %>%
  select(Cluster, Description, log_q) %>%
  pivot_wider(
    names_from  = Cluster,
    values_from = log_q,
    values_fill = 0
  )

# Step 2.4: Use Description as row names to obtain a pathway × group matrix
heat_matrix <- as.data.frame(heat_data)
rownames(heat_matrix) <- heat_matrix$Description
heat_matrix <- heat_matrix[, -1, drop = FALSE]

# Save the raw matrix for inspection / downstream reuse
write.csv(heat_matrix, "cluster_top15_heat_matrix.csv", row.names = TRUE)

# -----------------------------------------------------------------------------
# 3. Define pathway categories (for row annotation)
# -----------------------------------------------------------------------------
# Build a table that assigns each pathway Description to a broad functional
# category (e.g. lipid metabolism, amino acid metabolism, etc.)
pathway_names <- rownames(heat_matrix)

pathway_category <- tibble(
  Description = pathway_names,
  Category = case_when(
    # Lipid metabolism
    Description %in% c(
      "Biosynthesis of unsaturated fatty acids",
      "Linoleic acid metabolism",
      "Arachidonic acid metabolism",
      "Steroid hormone biosynthesis",
      "Retinol metabolism",
      "Ovarian steroidogenesis"
    ) ~ "Lipid metabolism",
    
    # Amino acid metabolism
    Description %in% c(
      "Biosynthesis of amino acids",
      "Valine, leucine and isoleucine biosynthesis",
      "Valine, leucine and isoleucine degradation",
      "Alanine, aspartate and glutamate metabolism",
      "Cysteine and methionine metabolism",
      "Arginine and proline metabolism",
      "Tryptophan metabolism",
      "Tyrosine metabolism",
      "Arginine biosynthesis"
    ) ~ "Amino acid metabolism",
    
    # Energy & carbohydrate metabolism
    Description %in% c(
      "2-Oxocarboxylic acid metabolism",
      "Butanoate metabolism",
      "Carbohydrate digestion and absorption"
    ) ~ "Energy & carbohydrate metabolism",
    
    # Nucleotide & vitamin metabolism
    Description %in% c(
      "Nucleotide metabolism",
      "Purine metabolism",
      "Vitamin digestion and absorption",
      "Sulfur metabolism"
    ) ~ "Nucleotide & vitamin metabolism",
    
    # Hormonal & signaling pathways
    Description %in% c(
      "Serotonergic synapse",
      "Aldosterone synthesis and secretion",
      "Parathyroid hormone synthesis, secretion and action",
      "Prolactin signaling pathway",
      "Sphingolipid signaling pathway"
    ) ~ "Hormonal & signaling pathways",
    
    # Secondary metabolite biosynthesis
    Description %in% c(
      "Biosynthesis of phenylpropanoids",
      "Glucosinolate biosynthesis",
      "Biosynthesis of alkaloids derived from shikimate pathway",
      "Tropane, piperidine and pyridine alkaloid biosynthesis",
      "Degradation of flavonoids",
      "Furfural degradation",
      "Aminobenzoate degradation"
    ) ~ "Secondary metabolite biosynthesis",
    
    # Barrier & immune function
    Description %in% c(
      "Bile secretion",
      "Cornified envelope formation"
    ) ~ "Barrier & immune function",
    
    # Transporters
    Description %in% c(
      "ABC transporters"
    ) ~ "Transporters",
    
    # Default
    TRUE ~ "Unclassified"
  )
)

pathway_category <- as.data.frame(pathway_category)
rownames(pathway_category) <- pathway_category$Description

# -----------------------------------------------------------------------------
# 4. Reorder heat_matrix by category, then by pathway name
# -----------------------------------------------------------------------------
heat_matrix_sorted <- heat_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Description") %>%
  left_join(pathway_category, by = "Description") %>%
  arrange(Category, Description)

heat_matrix_sorted_mat <- heat_matrix_sorted %>%
  column_to_rownames("Description") %>%
  select(-Category)

# -----------------------------------------------------------------------------
# 5. Define annotation colors for pathway categories
# -----------------------------------------------------------------------------
categories_used <- unique(pathway_category$Category)

category_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(categories_used)),
  categories_used
)

ann_colors <- list(Category = category_colors)

# -----------------------------------------------------------------------------
# 6. Plot heatmap (Figure 4b)
# -----------------------------------------------------------------------------
out_pdf <- "Top15_KEGG_Pathway_Heatmap.pdf"

pheatmap::pheatmap(
  mat = as.matrix(heat_matrix_sorted_mat),
  color = colorRampPalette(c("white", "red4"))(100),
  cluster_rows = FALSE,   # keep category-based ordering of pathways
  cluster_cols = TRUE,    # allow clustering of groups/columns
  annotation_row = data.frame(
    Category = pathway_category[rownames(heat_matrix_sorted_mat), "Category", drop = FALSE]
  ),
  annotation_colors = ann_colors,
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 11,
  cellwidth = 30,
  cellheight = 12,
  main = "Top 15 KEGG Pathways per Group (Grouped by Category)",
  filename = out_pdf
)
