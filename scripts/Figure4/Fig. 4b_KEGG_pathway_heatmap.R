############################################################
## Figure 4b KEGG pathway enrichment heatmap
# ============================================================
# 1. Load packages
# ============================================================
library(tidyverse)
library(Seurat)
library(igraph)
library(ggraph)
library(scales)
library(pheatmap)
library(RColorBrewer)


# ============================================================
# 2. Input and output paths
# ============================================================
base_dir   <- "Figure4"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

kegg_enrichment_csv <- file.path(input_dir, "Metabolomics_KEGG_pathway_enrichment_by_cluster.csv")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

kegg_meta <- read.csv(kegg_enrichment_csv, stringsAsFactors = FALSE)
head(kegg_meta)

# ============================================================
# 3. Select top 15 pathways per cluster
# ============================================================
top_paths <- kegg_meta %>%
  mutate(log_q = -log10(qvalue)) %>%
  group_by(Cluster) %>%
  arrange(RichFactor) %>%
  slice_head(n = 15) %>%
  ungroup()

top_descriptions <- unique(top_paths$Description)
top_descriptions

heat_data <- kegg_meta %>%
  mutate(log_q = -log10(qvalue)) %>%
  filter(Description %in% top_descriptions) %>%
  select(Cluster, Description, log_q) %>%
  pivot_wider(names_from = Cluster, values_from = log_q, values_fill = 0)

heat_matrix <- as.data.frame(heat_data)
rownames(heat_matrix) <- heat_matrix$Description
heat_matrix <- heat_matrix[, -1, drop = FALSE]

write.csv(
  heat_matrix,
  file.path(output_dir, "cluster_top15_heat_matrix.csv"),
  row.names = TRUE
)

# ============================================================
# 4. Create pathway category annotation
# ============================================================
pathway_names <- rownames(heat_matrix)
pathway_names

pathway_category <- tibble(
  Description = pathway_names,
  Category = case_when(
    Description %in% c("Biosynthesis of unsaturated fatty acids",
                       "Linoleic acid metabolism",
                       "Arachidonic acid metabolism",
                       "Steroid hormone biosynthesis",
                       "Retinol metabolism",
                       "Ovarian steroidogenesis") ~ "Lipid metabolism",
    
    Description %in% c("Biosynthesis of amino acids",
                       "Valine, leucine and isoleucine biosynthesis",
                       "Valine, leucine and isoleucine degradation",
                       "Alanine, aspartate and glutamate metabolism",
                       "Cysteine and methionine metabolism",
                       "Arginine and proline metabolism",
                       "Tryptophan metabolism",
                       "Tyrosine metabolism",
                       "Arginine biosynthesis") ~ "Amino acid metabolism",
    
    Description %in% c("2-Oxocarboxylic acid metabolism",
                       "Butanoate metabolism",
                       "Carbohydrate digestion and absorption") ~ "Energy & carbohydrate metabolism",
    
    Description %in% c("Nucleotide metabolism",
                       "Purine metabolism",
                       "Vitamin digestion and absorption",
                       "Sulfur metabolism") ~ "Nucleotide & vitamin metabolism",
    
    Description %in% c("Serotonergic synapse",
                       "Aldosterone synthesis and secretion",
                       "Parathyroid hormone synthesis, secretion and action",
                       "Prolactin signaling pathway",
                       "Sphingolipid signaling pathway") ~ "Hormonal & signaling pathways",
    
    Description %in% c("Biosynthesis of phenylpropanoids",
                       "Glucosinolate biosynthesis",
                       "Biosynthesis of alkaloids derived from shikimate pathway",
                       "Tropane, piperidine and pyridine alkaloid biosynthesis",
                       "Degradation of flavonoids",
                       "Furfural degradation",
                       "Aminobenzoate degradation") ~ "Secondary metabolite biosynthesis",
    
    Description %in% c("Bile secretion",
                       "Cornified envelope formation") ~ "Barrier & immune function",
    
    Description %in% c("ABC transporters") ~ "Transporters",
    
    TRUE ~ "Unclassified"
  )
)

pathway_category <- as.data.frame(pathway_category)
row.names(pathway_category) <- pathway_category$Description

# ============================================================
# 5. Sort heatmap matrix by pathway category
# ============================================================
heat_matrix_sorted <- heat_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Description") %>%
  left_join(pathway_category, by = "Description") %>%
  arrange(Category, Description)

heat_matrix_sorted_mat <- heat_matrix_sorted %>%
  column_to_rownames("Description") %>%
  select(-Category)

# ============================================================
# 6. Create annotation colours
# ============================================================
categories_used <- unique(pathway_category$Category)

category_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(categories_used)),
  categories_used
)

ann_colors <- list(Category = category_colors)

# ============================================================
# 7. Draw heatmap
# ============================================================
out_pdf <- file.path(output_dir, "Top15_KEGG_Pathway_Heatmap.pdf")

pheatmap::pheatmap(
  mat = as.matrix(heat_matrix_sorted_mat),
  color = colorRampPalette(c("white", "red4"))(100),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
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
