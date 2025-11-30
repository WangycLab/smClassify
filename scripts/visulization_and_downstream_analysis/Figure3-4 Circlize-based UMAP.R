# ============================================================
# Circlize-based UMAP visualization & subtype summaries
# Generates Figure 3c, 3d and related panels
#
# Major steps:
#   1) Load pre-processed Seurat object (stepF_seu_SCT_clustered.rds).
#   2) Add group / position / phenotype metadata and functional cluster labels.
#   3) Prepare circular UMAP layout using plot1cell::prepare_circlize_data().
#   4) Draw circlize UMAP plots for cell types, plus outer tracks for region
#      and phenotype.
#   5) Build a Sankey diagram (Function ↔ Cluster).
#   6) Draw stacked barplots:
#        - Cell-type composition by phenotype
#        - Cell-type composition by loc
#        - Total cell counts per celltype
#
# Notes:
#   - All analysis parameters (pt.size, kde2d.n, colors, sizes, etc.)
#     are identical to the original script.
#   - Only the root directory is generalized to `project_dir`.
# ============================================================

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

req_pkgs <- c(
  "Seurat", "tidyverse", "Matrix", "data.table", "ggrepel",
  "RColorBrewer", "gplots", "rtracklayer", "scales",
  "pheatmap", "ComplexHeatmap", "circlize", "paletteer", "plot1cell"
)

for (pkg in req_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

library(plot1cell)
library(circlize)
library(networkD3)
library(htmlwidgets)
library(webshot2)   # Requires Chrome/Chromium for html → image export
library(patchwork)

# ------------------------------------------------------------
# 2. Global parameters and colors
# ------------------------------------------------------------

# Root project directory (modify this for your own environment)
project_dir <- "path/to/project_root"  # <-- set to your project root
ROOT_DIR    <- project_dir
out_dir     <- ROOT_DIR

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

nature_palette <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85"
)

nature_palette2 <- c(
  "#3C5488", "#F39B7F", "#8491B4", "#91D1C2",
  "#F5A0A1", "#C2B5D8", "#FFB6C1", "#7E6148",
  "#E64B35", "#4DBBD5", "#00A087", "#B09C85"
)

# Expand a base palette to length n (kept exactly as original behavior)
expand_palette <- function(base_colors, n) {
  if (n <= length(base_colors)) base_colors[1:n] else colorRampPalette(base_colors)(n)
}

# Simple helper to save ggplots as PNG + PDF
save_plot <- function(p, name, w = 7, h = 5, dpi = 600) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
}

# ------------------------------------------------------------
# 3. Load Seurat object and add metadata
# ------------------------------------------------------------

set.seed(1234)
setwd(ROOT_DIR)

# Pre-processed object from earlier steps
seu <- readRDS("stepF_seu_SCT_clustered.rds")

# group: collapse suffix "-digit" so that e.g. "Cecum-WT-1" & "Cecum-WT-2" share a group
seu$group <- sub("-[0-9]$", "", seu$orig.ident)

# position: infer region from orig.ident prefix
seu$position <- dplyr::case_when(
  grepl("^Cecum",  seu$orig.ident) ~ "Cecum",
  grepl("^Colon",  seu$orig.ident) ~ "Colon",
  grepl("^Rectum", seu$orig.ident) ~ "Rectum",
  TRUE ~ NA_character_
)

# phenotype: WT vs T2DM (kept identical to original)
seu$phenotype <- ifelse(grepl("WT", seu$orig.ident), "WT", "T2DM")

# region: identical to position (for compatibility with earlier code)
seu$region <- seu$position

# Cluster-level functional annotation:
#   - Cluster_functional_summary.csv should map clusters to functional labels.
anno_df <- read.csv("Cluster_functional_summary.csv", header = TRUE)

# Mapping: old cluster levels -> new functional IDs (2nd column of anno_df)
new.ids <- setNames(as.character(anno_df[[2]]), levels(seu))

seu@active.ident <- seu$seurat_clusters
seu <- RenameIdents(seu, new.ids)
seu$celltype <- seu@active.ident

# Save / reload annotated object (same logic as original)
save(seu, file = "seu_symbol_SCT_clustered_annotated.RData")
load("seu_symbol_SCT_clustered_annotated.RData")

# ------------------------------------------------------------
# 4. Circlize UMAP data preparation and plotting
# ------------------------------------------------------------

# Use functional celltype as active identity
seu@active.ident <- as.factor(seu$celltype)

# Prepare circlize layout data from Seurat object
circ_data <- prepare_circlize_data(seu, scale = 0.7)

# Colors for celltypes and outer tracks (region, phenotype)
celltype_colors <- c(
  "#4DBBD5", "#E64B35", "#00A087", "#3C5488",
  "#F39B7F", "#91D1C2", "#C2B5D8", "#7E6148"
)

region_cols    <- expand_palette(nature_palette,  length(unique(seu$region)))
phenotype_cols <- expand_palette(brewer.pal(6, "Set2"), length(unique(seu$phenotype)))

# ---- 4.1 PNG version: circlize UMAP by celltype ----
png(file.path(out_dir, "Circlize_Celltype.png"), width = 8000, height = 8000, res = 1200)
plot_circlize(
  circ_data,
  do.label  = FALSE,
  pt.size   = 0.4,
  col.use   = celltype_colors,
  bg.color  = "white",
  kde2d.n   = 200,
  repel     = TRUE,
  label.cex = 0.8
)
dev.off()

# ---- 4.2 PDF version (with tracks for region & phenotype) ----
cairo_pdf(
  file   = file.path(out_dir, "Circlize_Celltype.pdf"),
  width  = 5,
  height = 5,
  onefile = FALSE,
  family  = "sans"
)
op <- par(mar = rep(0, 4))
plot_circlize(
  circ_data,
  do.label  = FALSE,
  pt.size   = 0.4,
  col.use   = celltype_colors,
  bg.color  = "white",
  kde2d.n   = 200,
  repel     = TRUE,
  label.cex = 1
)
add_track(circ_data, group = "region",    colors = region_cols,    track_num = 2)
add_track(circ_data, group = "phenotype", colors = phenotype_cols, track_num = 3)
par(op)
dev.off()

# ---- 4.3 Alternative PDF (without extra tracks, larger canvas) ----
cairo_pdf(
  file   = file.path(out_dir, "Circlize_Celltype.pdf"),
  width  = 8,
  height = 8,
  onefile = FALSE,
  family  = "sans"
)
op <- par(mar = rep(0, 4))
plot_circlize(
  circ_data,
  do.label  = FALSE,
  pt.size   = 0.5,
  col.use   = celltype_colors,
  bg.color  = "white",
  kde2d.n   = 200,
  repel     = TRUE,
  label.cex = 1
)
par(op)
dev.off()

# ============================================================
# 5. Sankey network: Function ↔ Group
# ============================================================

# The following assumes:
#   - Idents(seu) = celltype (functional label)
#   - seu$group   = sample / group label

seu@active.ident <- seu$celltype

# Proportions of (celltype, group)
data <- as.data.frame(prop.table(table(Idents(seu), seu$group), margin = 1))
head(data)

df <- data.frame(
  Function = data$Var2,
  Cluster  = data$Var1,
  Freq     = data$Freq
)

# Nodes = unique Functions + unique Clusters
nodes <- data.frame(name = c(unique(df$Function), unique(df$Cluster)))

# Links with numeric IDs (0-based indices for networkD3)
links <- df |>
  mutate(
    IDsource = match(Function, nodes$name) - 1,
    IDtarget = match(Cluster,  nodes$name) - 1
  ) |>
  dplyr::select(IDsource, IDtarget, Freq)

# Draw Sankey diagram (interactive HTML widget)
p_sankey <- sankeyNetwork(
  Links     = links,
  Nodes     = nodes,
  Source    = "IDsource",
  Target    = "IDtarget",
  Value     = "Freq",
  NodeID    = "name",
  nodeWidth = 20,
  fontSize  = 30,
  sinksRight = FALSE
)
p_sankey  # Will display in an interactive R session

# ============================================================
# 6. Stacked barplots and cell count summary
# ============================================================

meta_df <- seu@meta.data

# -------------------------------
# 6.1 Stacked barplot: celltype × phenotype
# -------------------------------

frac_df <- meta_df %>%
  dplyr::group_by(celltype, phenotype) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(celltype) %>%
  dplyr::mutate(freq = n / sum(n))

p1 <- ggplot(frac_df, aes(
  x    = freq,
  y    = forcats::fct_rev(factor(celltype)),
  fill = phenotype
)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Fraction", y = NULL) +
  scale_fill_manual(values = c("T2DM" = "#E74C3C", "WT" = "#3498DB")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y   = element_text(size = 12),
    axis.text.x   = element_text(size = 12),
    axis.title.x  = element_text(size = 14),
    legend.title  = element_blank(),
    legend.position = "right",
    panel.grid    = element_blank()
  )

# -------------------------------
# 6.2 Stacked barplot: celltype × loc
# -------------------------------

frac_df <- meta_df %>%
  dplyr::group_by(celltype, loc) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(celltype) %>%
  dplyr::mutate(freq = n / sum(n))

p1loc <- ggplot(frac_df, aes(
  x    = freq,
  y    = forcats::fct_rev(factor(celltype)),
  fill = loc
)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Fraction", y = NULL) +
  scale_fill_manual(values = c(
    "Cecum"  = "#1f77b4",
    "Colon"  = "#2ca02c",
    "Rectum" = "#ff7f0e"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y   = element_blank(),
    axis.text.x   = element_text(size = 12),
    axis.title.x  = element_text(size = 14),
    legend.title  = element_blank(),
    legend.position = "right",
    panel.grid    = element_blank()
  )

# -------------------------------
# 6.3 Total cell numbers per celltype
# -------------------------------

meta_df2 <- as.data.frame(seu@meta.data)

total_df <- meta_df2 %>%
  dplyr::count(celltype)

p2 <- ggplot(total_df, aes(
  x = n / 1000,
  y = forcats::fct_rev(factor(celltype))
)) +
  geom_col(fill = "#BCBD22", width = 0.7) +
  labs(x = "Cell Nums (x10^3)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y   = element_blank(),
    axis.text.x   = element_text(size = 12),
    axis.title.x  = element_text(size = 14),
    panel.grid    = element_blank()
  )

# -------------------------------
# 6.4 Combine the three panels side by side
# -------------------------------

p_combined <- p1 + p1loc + p2 + plot_layout(widths = c(1, 1, 0.8))
p_combined

ggsave(
  "seu_Cellsubtype_Phenotype_Stacked.pdf",
  p_combined,
  width = 8,
  height = 2.5,
  units = "in",
  dpi = 600
)
ggsave(
  "seu_Cellsubtype_Phenotype_Stacked.png",
  p_combined,
  width = 8,
  height = 2.5,
  units = "in",
  dpi = 600
)
