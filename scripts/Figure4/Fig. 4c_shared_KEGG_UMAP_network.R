############################################################
## Figure 4c shared KEGG pathway network over UMAP

############################################################

# ============================================================
# 1. Load packages
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(igraph)
  library(ggraph)
  library(scales)
  library(ggrepel)
  library(viridis)
})

library(dplyr)

# ============================================================
# 2. Input and output paths
# ============================================================
base_dir   <- "Figure4"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")


annotated_seurat_rds <- file.path(
  "Figure3/Output/Figure3_mouse_bacterial_functional_cluster_annotated.rds"
)

cluster_summary_csv <- file.path(input_dir, "Functional_cluster_summary.csv")
sc_kegg_csv         <- file.path(input_dir, "Single_cell_KEGG_enrichment.csv")
meta_kegg_csv       <- file.path(input_dir, "Metabolomics_KEGG_pathway_enrichment_by_cluster.csv")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 3. Load data
# ============================================================
seu <- readRDS(annotated_seurat_rds)

seu@active.ident <- seu$celltype

lookup    <- read_csv(cluster_summary_csv, show_col_types = FALSE)
kegg_sc   <- read_csv(sc_kegg_csv, show_col_types = FALSE)
kegg_meta <- read_csv(meta_kegg_csv, show_col_types = FALSE)

# ============================================================
# 4. Format functional cluster labels
# ============================================================
cluster_levels <- lookup$Cluster
cluster_labels <- lookup$`Figure label (short)`

kegg_sc <- kegg_sc %>%
  mutate(Cluster = paste0("C", cluster)) %>%
  mutate(
    Cluster = factor(
      Cluster,
      levels = cluster_levels,
      labels = cluster_labels
    )
  )

# ============================================================
# 5. Extract UMAP coordinates
# ============================================================
umap_df <- Embeddings(seu, "umap") %>%
  as_tibble(rownames = "CellBarcode") %>%
  mutate(Cluster = factor(Idents(seu)))

cluster_centroid <- umap_df %>%
  group_by(Cluster) %>%
  summarise(across(c(umap_1, umap_2), mean), .groups = "drop")

# ============================================================
# 6. Harmonise enrichment tables
# ============================================================
kegg_meta <- kegg_meta %>%
  dplyr::rename(
    seuGroup = Cluster,
    Pathway  = Description
  )

kegg_sc <- kegg_sc %>%
  dplyr::rename(Pathway = Description)

# ============================================================
# 7. Define metabolomics nodes
# ============================================================
meta_nodes <- tibble(
  name = unique(kegg_meta$seuGroup),
  x    = max(cluster_centroid$umap_1, na.rm = TRUE) + 8,
  y    = seq(
    min(cluster_centroid$umap_2, na.rm = TRUE),
    max(cluster_centroid$umap_2, na.rm = TRUE),
    length.out = n_distinct(kegg_meta$seuGroup)
  ),
  type = "Metabolomics"
)

# ============================================================
# 8. Build shared pathway edges
# ============================================================
edge_df <- inner_join(kegg_sc, kegg_meta, by = "Pathway")

edge_df2 <- edge_df %>%
  mutate(
    rf  = RichFactor.x,
    sig = if ("p.adjust.x" %in% names(.)) -log10(p.adjust.x)
    else if ("pvalue.x" %in% names(.)) -log10(pvalue.x)
    else NA_real_
  ) %>%
  mutate(
    rf_s  = rescale(rf,  to = c(0, 1), from = range(rf,  na.rm = TRUE)),
    sig_s = if (!all(is.na(sig))) rescale(sig, to = c(0, 1), from = range(sig, na.rm = TRUE)) else NA_real_,
    score = if (!all(is.na(sig_s))) 2 * (rf_s * sig_s) / (rf_s + sig_s) else rf_s
  )

edge_tbl <- edge_df2 %>%
  transmute(
    from    = as.character(Cluster),
    to      = as.character(seuGroup),
    score,
    rf      = RichFactor.x,
    pathway = Pathway
  ) %>%
  dplyr::filter(!is.na(from), !is.na(to), !is.na(score)) %>%
  group_by(from, to) %>%
  summarise(
    weight  = mean(score, na.rm = TRUE),
    n_path  = n_distinct(pathway),
    mean_rf = mean(rf, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 9. Build network nodes
# ============================================================
cluster_nodes <- cluster_centroid %>%
  transmute(
    name = as.character(Cluster),
    x = umap_1,
    y = umap_2,
    type = "Cluster"
  )

node_df <- bind_rows(cluster_nodes, meta_nodes) %>%
  distinct(name, .keep_all = TRUE)

# ============================================================
# 10. Generate graph layout
# ============================================================
g <- graph_from_data_frame(d = edge_tbl, vertices = node_df, directed = FALSE)

V(g)$x <- node_df$x[match(V(g)$name, node_df$name)]
V(g)$y <- node_df$y[match(V(g)$name, node_df$name)]

layout_manual <- create_layout(g, layout = "manual", x = V(g)$x, y = V(g)$y)

# ============================================================
# 11. Plot network with labels
# ============================================================
p1 <- ggraph(layout_manual) +
  geom_point(
    data = umap_df, aes(umap_1, umap_2),
    colour = "grey85", size = 0.2, alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_edge_link(aes(width = weight),
                 colour = "#B2182B", alpha  = 0.5, lineend = "round") +
  scale_edge_width(range = c(0.15, 2.8), guide = "none") +
  geom_node_point(aes(colour = type), size = 4.6, stroke = 0) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.6, fontface = "bold") +
  scale_colour_manual(values = c(Cluster = "#2166AC", Metabolomics = "#4DAF4A")) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin = margin(14, 18, 14, 18, "pt")
  ) +
  ggtitle("Shared KEGG pathways")

p1

ggsave(
  file.path(output_dir, "UMAP_Linked_KEGG.pdf"),
  p1,
  width = 5.5,
  height = 5.0,
  dpi = 450
)

# ============================================================
# 12. Plot network without labels
# ============================================================
p1 <- ggraph(layout_manual) +
  geom_point(
    data = umap_df, aes(umap_1, umap_2),
    colour = "grey85", size = 0.2, alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_edge_link(aes(width = weight),
                 colour = "#B2182B", alpha  = 0.5, lineend = "round") +
  scale_edge_width(range = c(0.15, 2.8), guide = "none") +
  geom_node_point(aes(colour = type), size = 4.6, stroke = 0) +
  geom_node_text(aes(label = name), repel = TRUE, size = 0, fontface = "bold") +
  scale_colour_manual(values = c(Cluster = "#2166AC", Metabolomics = "#4DAF4A")) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin = margin(14, 18, 14, 18, "pt")
  ) +
  ggtitle("Shared KEGG pathways")

p1

ggsave(
  file.path(output_dir, "UMAP_Linked_KEGG.png"),
  p1,
  width = 5.5,
  height = 5.0,
  dpi = 450,
  bg = "white"
)
