# This script:
#   1) Loads a clustered Seurat object (stepF_seu_SCT_clustered.rds).
#   2) Runs FindAllMarkers() to identify cluster marker genes.
#   3) Selects top 50 marker genes per cluster, then top 10 for heatmap.
#   4) Joins marker genes to a wide GTF annotation table.
#   5) Computes average expression (SCT assay) per cluster.
#   6) Builds and saves a ComplexHeatmap heatmap of Top10 marker genes.
#
# Inputs (expected in project_dir):
#   - stepF_seu_SCT_clustered.rds
#       Seurat object with SCT assay and clustering (seurat_clusters).
#   - gtf_gene_annotation.csv
#       GTF-wide gene annotation, containing at least a "symbol" column.
#
# Outputs:
#   - Total cells_markers_min.pct0.05.csv
#       All marker genes from FindAllMarkers() (min.pct = 0.05).
#   - Total cells_markers_top50.csv
#       Top 50 marker genes per cluster (by avg_log2FC).
#   - Total cells_markers_top50_gtf.csv
#       Top 50 markers joined to GTF (symbol + product).
#   - All_genes_symbol_gtf.csv
#       All genes mapped to GTF annotation (symbol / product / gene_id / gene.x).
#   - Top10_marker_genes_heatmap.png
#   - Top10_marker_genes_heatmap.pdf
#       ComplexHeatmap heatmap of Top10 marker genes across clusters.
#
# NOTE:
#   - All numeric parameters and thresholds are kept exactly as in the
#     original script (no changes).
#   - Only comments, structure, and path handling have been made
#     more GitHub-friendly.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(future)
  library(tibble)
  library(paletteer)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)          # grid.text, gpar
  library(RColorBrewer)  # brewer.pal
  library(ggplot2)
  library(scales)        # hue_pal
})

# ------------------------------------------------------------
# 0) Global options and project directory
# ------------------------------------------------------------
set.seed(1234)

# User should set this to the project root directory when running the script.
project_dir <- "path/to/project_root"   # <-- modify this to your environment
dir.create(project_dir, showWarnings = FALSE, recursive = TRUE)
setwd(project_dir)

# ------------------------------------------------------------
# 0.1) Check required input files
# ------------------------------------------------------------
if (!file.exists("stepF_seu_SCT_clustered.rds")) {
  stop("Missing file: stepF_seu_SCT_clustered.rds (clustered Seurat object).")
}
if (!file.exists("gtf_gene_annotation.csv")) {
  stop("Missing file: gtf_gene_annotation.csv (GTF-wide annotation).")
}

# Load clustered Seurat object and GTF-wide annotation
seu <- readRDS("stepF_seu_SCT_clustered.rds")
sce <- seu   # for compatibility with original code using 'sce'
gtf_wide <- fread("gtf_gene_annotation.csv", check.names = FALSE)

# ============================================================
# 1) Find cluster markers using FindAllMarkers()
# ============================================================

## =========================
## Step H: Add markers
## =========================
markers <- FindAllMarkers(
  seu,
  only.pos        = TRUE,
  logfc.threshold = 1,
  min.pct         = 0.05
)
write.csv(markers, "Total cells_markers_min.pct0.05.csv")
table(markers$cluster)

# Select Top 50 marker genes per cluster
topN <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()
table(topN$cluster)
write.csv(topN, "Total cells_markers_top50.csv")

# ============================================================
# 2) Join marker genes with GTF annotation
# ============================================================

# Original convention: replace "-" with "_" to match gtf_wide$symbol
topN$symbol <- gsub("-", "_", topN$gene)

# Left join with GTF-wide annotation using 'symbol'
topN_gtf <- topN %>% left_join(gtf_wide, by = "symbol")

# Keep core columns and remove duplicates
topN_gtf_subset <- unique(topN_gtf[, c(
  "p_val", "avg_log2FC", "pct.1", "pct.2",
  "p_val_adj", "cluster", "symbol", "product"
)])
write.csv(topN_gtf_subset, "Total cells_markers_top50_gtf.csv", row.names = FALSE)

# ============================================================
# 3) Map all genes to GTF annotation (for reference)
# ============================================================

## =========================
## Add all genes to GTF
## =========================
all_genes <- rownames(sce)
all_genes_df <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

# Match naming convention used above: "-" -> "_"
all_genes_df$symbol <- gsub("-", "_", all_genes_df$gene)

# Left-join to GTF to obtain product / gene_id if available
all_genes_df <- all_genes_df %>% left_join(gtf_wide, by = "symbol")

# Deduplicate and keep selected columns
all_genes_df <- unique(all_genes_df[, c("symbol", "product", "gene_id", "gene.x")])
write.csv(all_genes_df, "All_genes_symbol_gtf.csv", row.names = FALSE)
colnames(all_genes_df)  # quick sanity check

# ============================================================
# 4) Select Top10 marker genes per cluster for heatmap
# ============================================================

## =========================
## Top10 Marker genes
## =========================
top10 <- topN_gtf_subset %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()

# Remove duplicated symbols so each gene appears at most once
top10 <- top10[!duplicated(top10$symbol), ]

# For heatmap labels: convert "_" back to "-" in symbol
top10$symbol <- gsub("_", "-", top10$symbol)
top10

genes <- unique(top10$symbol)
genes  # the set of genes to be used in the heatmap

# ============================================================
# 5) Average expression matrix across clusters (SCT assay)
# ============================================================

## =========================
## Average expression across clusters
## =========================
aver_dt <- AverageExpression(
  seu,
  assays   = "SCT",
  features = genes,
  group.by = "seurat_clusters",
  slot     = "data"
)$SCT %>% as.data.frame()

# Save gene symbol as a column for merging
aver_dt$symbol <- rownames(aver_dt)

# Merge basic annotation (product, cluster) from Top10 table
top102 <- top10[, c("symbol", "product", "cluster")]
aver_dt <- left_join(aver_dt, top102, by = "symbol")

# Build combined row labels: "symbol__product"
aver_dt$combined <- paste(aver_dt$symbol, aver_dt$product, sep = "__")
rownames(aver_dt) <- make.unique(aver_dt$combined)

# Extract pure expression matrix (all columns except the last 4)
expr_mat <- aver_dt[, 1:(ncol(aver_dt) - 4)]

# ============================================================
# 6) Prepare annotations for ComplexHeatmap
# ============================================================

## =========================
## Heatmap annotations
## =========================

# Row annotation: original cluster of each marker gene
gene_anno <- data.frame(gene_anno = top10$cluster, row.names = top10$symbol)

# Column annotation: cluster IDs (seurat_clusters)
celltypes <- colnames(expr_mat)
cell_anno <- data.frame(cell_anno = celltypes, row.names = celltypes)

# Color mapping for clusters (columns)
celltype_col <- setNames(
  brewer.pal(length(unique(cell_anno$cell_anno)), "Set3"),
  unique(cell_anno$cell_anno)
)

# Color mapping for gene cluster annotation (rows)
gene_anno_col <- setNames(
  brewer.pal(length(unique(gene_anno$gene_anno)), "Set1"),
  unique(gene_anno$gene_anno)
)

# ============================================================
# 7) Build ComplexHeatmap object
# ============================================================

## =========================
## ComplexHeatmap
## =========================

# Row-wise Z-score scaling (standardization per gene)
expr_scaled <- t(scale(t(expr_mat)))

# Column annotation (cluster labels on top of heatmap)
col_anno <- HeatmapAnnotation(
  cell = cell_anno$cell_anno,
  col  = list(cell = celltype_col),
  show_annotation_name = FALSE,
  gp   = gpar(col = 'white', lwd = 1.5)
)

# Row-side text annotation: full gene label (symbol__product)
row_anno <- rowAnnotation(
  foo = anno_text(
    rownames(expr_scaled),
    location = 0,
    just = "left",
    gp = gpar(col = "black", fontface = "italic"),
    width = max_text_width(rownames(expr_scaled)) * 1.25
  )
)

# Diverging color scale for (scaled) expression
mycol2 <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))

# Main heatmap (no row/column clustering to keep cluster order fixed)
ht <- Heatmap(
  expr_scaled,
  name = "expression",
  col  = mycol2,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 60,
  row_names_gp     = gpar(fontsize = 10, fontface = "italic"),
  rect_gp          = gpar(col = "white", lwd = 1.5),
  top_annotation   = col_anno
) + row_anno

# Draw heatmap to current device (useful in interactive sessions)
ht

# ============================================================
# 8) Save heatmap to PNG and PDF
# ============================================================

# High-resolution PNG output
png("Top10_marker_genes_heatmap.png", width = 6500, height = 8000, res = 500)
draw(ht)
grid.text(
  "Top10 marker genes across clusters",
  x = 0.5, y = 1.5,
  gp = gpar(fontsize = 18, fontface = "bold")
)
dev.off()

# PDF output (print-friendly / for publication)
pdf("Top10_marker_genes_heatmap.pdf", width = 9, height = 12)
draw(ht)
grid.text(
  "Top10 marker genes across clusters",
  x = 0.5, y = 0.98,
  gp = gpar(fontsize = 18, fontface = "bold")
)
dev.off()
