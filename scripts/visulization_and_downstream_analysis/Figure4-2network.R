# =============================================================================
# Shared KEGG Pathway Network: scRNA-seq Clusters x Metabolomics Groups
# Visualization + Quantification 
# =============================================================================

# 1) Dependencies -------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)   # Data wrangling & plotting
  library(Seurat)      # scRNA-seq object access
  library(igraph)      # Graph data structures
  library(ggraph)      # Graph visualization (ggplot2 grammar)
  library(scales)      # Scaling utilities (rescale, palettes)
  library(ggrepel)     # Non-overlapping text labels
  library(viridis)     # Color palettes (optional)
})
library(dplyr)

setwd("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009")
# 2) Output directory ---------------------------------------------------------
plot_dir <- file.path("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009")

# 3) Load Seurat object -------------------------------------------------------
# Expect an object named `seu` with UMAP embeddings and `celltype` meta.
load("seu_symbol_SCT_clustered_annotated.RData")
seu@active.ident <- seu$celltype

# 4) Load lookup & enrichment results ----------------------------------------
# - lookup: cluster label mapping table
# - kegg_sc: significant KEGG terms from scRNA-seq (per cluster)
# - kegg_meta: significant KEGG terms from metabolomics groups
lookup    <- read_csv("Cluster_functional_summary.csv", show_col_types = FALSE)
kegg_sc   <- read_csv("kegg_sig.csv", show_col_types = FALSE)
kegg_meta <- read_csv("cpd_DE_enrich.csv", show_col_types = FALSE)

# 5) Map cluster short labels for plotting -----------------------------------
cluster_levels <- paste0("C", lookup$Cluster)
cluster_labels <- lookup$`Figure label (short)`

kegg_sc <- kegg_sc %>%
  mutate(Cluster = paste0("C", cluster)) %>%            # unify "C1/C2/..." style
  mutate(Cluster = factor(Cluster,
                          levels = cluster_levels,
                          labels = cluster_labels))     # map to short labels

# 6) UMAP coordinates & cluster centroids ------------------------------------
umap_df <- Embeddings(seu, "umap") %>%
  as_tibble(rownames = "CellBarcode") %>%
  mutate(Cluster = factor(Idents(seu)))

# cluster centroids used as node locations
cluster_centroid <- umap_df %>%
  group_by(Cluster) %>%
  summarise(across(c(umap_1, umap_2), mean), .groups = "drop")

# 7) Harmonize column names ---------------------------------------------------
# - kegg_meta 'Cluster' column holds metabolomics groups; rename to 'seuGroup'
# - both sc & meta use 'Pathway' for the pathway descriptor
kegg_meta <- kegg_meta %>%
  dplyr::rename(
    seuGroup = Cluster,
    Pathway  = Description
  )
kegg_sc <- kegg_sc %>%
  dplyr::rename(Pathway = Description)

# 8) Create metabolomics side nodes ------------------------------------------
# Place metabolomics group nodes to the right of UMAP cloud (x-shifted),
# and vertically distribute them across UMAP y-range for readability.
meta_nodes <- tibble(
  name = unique(kegg_meta$seuGroup),
  x    = max(cluster_centroid$umap_1, na.rm = TRUE) + 8,
  y    = seq(min(cluster_centroid$umap_2, na.rm = TRUE),
             max(cluster_centroid$umap_2, na.rm = TRUE),
             length.out = n_distinct(kegg_meta$seuGroup)),
  type = "Metabolomics"
)

# 9) Join scRNA & metabolomics by pathway ------------------------------------
# edge_df contains matched rows (shared pathways)
edge_df <- inner_join(kegg_sc, kegg_meta, by = "Pathway")

# Compute composite score per shared pathway pair:
# - rf: enrichment strength proxy (RichFactor.x)
# - sig: significance proxy (-log10 p.adjust.x / pvalue.x)
# - rf_s/sig_s: [0,1] rescaled
# - score: harmonic mean of rf_s & sig_s when possible; fallback to rf_s
edge_df2 <- edge_df %>%
  mutate(
    rf  = RichFactor.x,
    sig = if ("p.adjust.x" %in% names(.)) -log10(p.adjust.x)
    else if ("pvalue.x" %in% names(.)) -log10(pvalue.x)
    else NA_real_
  ) %>%
  mutate(
    rf_s  = rescale(rf,  to = c(0,1), from = range(rf,  na.rm = TRUE)),
    sig_s = if (!all(is.na(sig))) rescale(sig, to = c(0,1), from = range(sig, na.rm = TRUE)) else NA_real_,
    score = if (!all(is.na(sig_s))) 2*(rf_s*sig_s)/(rf_s+sig_s) else rf_s
  )

# Aggregate to one edge per (Cluster x MetabolomicsGroup) ---------------------
# - weight: mean score over shared pathways
# - n_path: number of distinct shared pathways
# - mean_rf: average RichFactor for reference
edge_tbl <- edge_df2 %>%
  transmute(
    from    = as.character(Cluster),   # sc cluster
    to      = as.character(seuGroup),  # metabolomics group
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

# 10) Build node table for graph layout --------------------------------------
cluster_nodes <- cluster_centroid %>%
  transmute(name = as.character(Cluster),
            x = umap_1, y = umap_2, type = "Cluster")

node_df <- bind_rows(cluster_nodes, meta_nodes) %>%
  distinct(name, .keep_all = TRUE)

# 11) Build graph + manual layout --------------------------------------------
g <- graph_from_data_frame(d = edge_tbl, vertices = node_df, directed = FALSE)
V(g)$x <- node_df$x[match(V(g)$name, node_df$name)]
V(g)$y <- node_df$y[match(V(g)$name, node_df$name)]
layout_manual <- create_layout(g, layout = "manual", x = V(g)$x, y = V(g)$y)

# 12) Plot A: Network over UMAP (Nature-style) --------------------------------
# - Grey UMAP points as background
# - Edges weighted by mean shared-pathway score
# - Node colors by type (Cluster / Metabolomics)
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

# Save network figure (PNG + PDF). 

ggsave(file.path(plot_dir, "UMAP_Linked_KEGG.pdf"), p1,
       width = 5.5, height = 5.0, dpi = 450)


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
ggsave(file.path(plot_dir, "UMAP_Linked_KEGG.png"), p1,
       width = 5.5, height = 5.0, dpi = 450, bg = "white")
# 13) Extract KO-level stats (if KO column exists) ----------------------------
# - edge_df is sc meta joined by pathway
# - Expect KO column from sc enrichment (assumed "ID.x"). If different, change name only.
colnames(edge_df)
if ("ID.x" %in% names(edge_df)) {
  ko_df <- edge_df %>%
    transmute(
      Cluster  = as.character(Cluster),
      seuGroup = as.character(seuGroup),
      Pathway  = Pathway,
      KO       = ID.x,                  # <-- adjust if your KO column differs
      rf       = RichFactor.x,
      sig      = if ("p.adjust.x" %in% names(.)) -log10(p.adjust.x)
      else if ("pvalue.x" %in% names(.)) -log10(pvalue.x)
      else NA_real_
    ) %>%
    mutate(
      rf_s  = rescale(rf,  to = c(0,1), from = range(rf,  na.rm = TRUE)),
      sig_s = if (!all(is.na(sig))) rescale(sig, to = c(0,1), from = range(sig, na.rm = TRUE)) else NA_real_,
      score = if (!all(is.na(sig_s))) 2*(rf_s*sig_s)/(rf_s+sig_s) else rf_s
    ) %>%
    dplyr::filter(!is.na(KO))
  
  # For each Cluster x Group, keep top-5 KOs by composite score (unchanged)
  top_ko <- ko_df %>%
    group_by(Cluster, seuGroup) %>%
    arrange(desc(score)) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  write_csv(top_ko, file.path(plot_dir, "TopKO_perClusterGroup.csv"))
  cat("\n>>> Top-KO statistics saved to:",
      file.path(plot_dir, "TopKO_perClusterGroup.csv"), "\n")
} else {
  cat("!!! KO column not found in edge_df. Please check your sc enrichment file (e.g., GeneID/ID column).\n")
}

# 14) Chord diagram: Cluster x KO x Group (tri-partite) -----------------------
# - Build two edge sets: Cluster x KO and KO x Group, then bind
library(circlize)

edges_cluster_ko <- top_ko %>%
  dplyr::select(from = Cluster, to = KO, value = score)

edges_ko_group <- top_ko %>%
  dplyr::select(from = KO, to = seuGroup, value = score)

edges_all <- bind_rows(edges_cluster_ko, edges_ko_group)

# Colors:
# - Clusters: Set2
# - Groups:   Set3
# - KOs:      neutral grey
clusters <- unique(top_ko$Cluster)
groups   <- unique(top_ko$seuGroup)
kos      <- unique(top_ko$KO)

grid.col <- c(
  setNames(RColorBrewer::brewer.pal(n = length(clusters), "Set2"), clusters),
  setNames(RColorBrewer::brewer.pal(n = length(groups), "Set3"), groups),
  setNames(rep("grey60", length(kos)), kos)
)

# Initialize circos canvas & draw chord diagram
circos.clear()
circos.par(start.degree = 90, gap.degree = 3,
           track.margin = c(0.01, 0.01), canvas.xlim = c(-1, 1), canvas.ylim = c(-1, 1))

chordDiagram(
  x = edges_all,
  grid.col = grid.col,
  transparency = 0.25,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = list(track.height = 0.08)
)

# Label styling (bold sector names on the first track)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(2),
              sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),
              cex = 0.6, font = 2)
}, bg.border = NA)

# Optional inspection
table(edges_all$to)

# Save PNG (pixel units; keep as is)
png(file.path("Cluster_KO_Group_Chord.png"),
    width = 2200, height = 2200, res = 300, bg = "white")
chordDiagram(
  x = edges_all,
  grid.col = grid.col,
  transparency = 0.25,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = list(track.height = 0.08)
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(2),
              sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),
              cex = 0.6, font = 2)
}, bg.border = NA)
dev.off()

# Save PDF (inches; NOTE: width/height are in inches, not pixels)
# If you find the PDF huge, set width/height to ~8-12.
pdf(file.path("Cluster_KO_Group_Chord.pdf"),
    width = 12, height = 12, bg = "white")  
chordDiagram(
  x = edges_all,
  grid.col = grid.col,
  transparency = 0.25,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = list(track.height = 0.08)
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(2),
              sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),
              cex = 0.6, font = 2)
}, bg.border = NA)
dev.off()

