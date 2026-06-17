# ============================================================
# Supplementary Figure 9c marker-gene heatmap

# =================


library(data.table)
library(Seurat)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)

set.seed(1234)

# ============================================================
# Input and output paths
# ======================================
base_dir   <- "Figure3"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

seurat_rds <- file.path(output_dir, "Figure3_mouse_bacterial_functional_cluster_annotated.rds")
gtf_csv    <- file.path(input_dir, "Mouse_gtf_gene_annotation.csv")

out_marker_all_csv     <- file.path(output_dir, "FigureS9c_all_cluster_markers_min_pct0.05.csv")
out_marker_top50_csv   <- file.path(output_dir, "FigureS9c_top50_cluster_markers.csv")
out_marker_top50_gtf   <- file.path(output_dir, "FigureS9c_top50_cluster_markers_with_annotation.csv")
out_all_genes_gtf_csv  <- file.path(output_dir, "FigureS9c_all_genes_with_gtf_annotation.csv")
out_top10_gtf_csv      <- file.path(output_dir, "FigureS9c_top10_marker_genes_with_annotation.csv")
out_heatmap_source_csv <- file.path(output_dir, "FigureS9c_marker_heatmap_source_data.csv")
out_heatmap_png        <- file.path(output_dir, "FigureS9c_top10_marker_gene_heatmap.png")
out_heatmap_pdf        <- file.path(output_dir, "FigureS9c_top10_marker_gene_heatmap.pdf")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# Load input data
# ============================================================
seu <- readRDS(seurat_rds)
gtf_wide <- data.table::fread(gtf_csv, check.names = FALSE) %>%
  as.data.frame()


seu$seurat_clusters <- factor(seu$seurat_clusters)
Idents(seu) <- seu$seurat_clusters

# ============================================================
# Marker-gene detection
# ============================================================
markers <- FindAllMarkers(
  object = seu,
  only.pos = TRUE,
  logfc.threshold = 1,
  min.pct = 0.05
)

write.csv(markers, out_marker_all_csv, row.names = FALSE)

# Select the top 50 marker genes per cluster by average log2 fold-change.
top50 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50, with_ties = FALSE) %>%
  ungroup()

write.csv(top50, out_marker_top50_csv, row.names = FALSE)

# ============================================================
# Marker-gene annotation
# ============================================================
top50$symbol <- gsub("-", "_", top50$gene)
top50_gtf <- top50 %>%
  left_join(gtf_wide, by = "symbol")

top50_gtf_subset <- top50_gtf %>%
  dplyr::select(any_of(c(
    "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj",
    "cluster", "gene", "symbol", "product", "gene_id", "ko"
  ))) %>%
  distinct()

write.csv(top50_gtf_subset, out_marker_top50_gtf, row.names = FALSE)

# Export GTF annotation for all genes present in the Seurat object.
all_genes_df <- data.frame(gene = rownames(seu), stringsAsFactors = FALSE)
all_genes_df$symbol <- gsub("-", "_", all_genes_df$gene)
all_genes_df <- all_genes_df %>%
  left_join(gtf_wide, by = "symbol") %>%
  dplyr::select(any_of(c("gene", "symbol", "product", "gene_id", "ko", "gene.x", "gene.y"))) %>%
  distinct()

write.csv(all_genes_df, out_all_genes_gtf_csv, row.names = FALSE)

# ============================================================
# Select Top10 marker genes per cluster
# ============================================================
top10 <- top50_gtf_subset %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  distinct(symbol, .keep_all = TRUE)

# Convert symbols back to the Seurat feature style.
top10$symbol <- gsub("_", "-", top10$symbol)

genes <- unique(top10$symbol)
genes <- genes[genes %in% rownames(seu)]


write.csv(top10, out_top10_gtf_csv, row.names = FALSE)

# ============================================================
# Average expression by Seurat cluster
# ============================================================
aver_dt <- AverageExpression(
  object = seu,
  assays = DefaultAssay(seu),
  features = genes,
  group.by = "seurat_clusters",
  slot = "data"
)[[DefaultAssay(seu)]] %>%
  as.data.frame()

aver_dt$symbol <- rownames(aver_dt)

top10_annot <- top10 %>%
  dplyr::select(any_of(c("symbol", "product", "cluster"))) %>%
  distinct(symbol, .keep_all = TRUE)

aver_dt <- left_join(aver_dt, top10_annot, by = "symbol")


aver_dt$combined <- paste(aver_dt$symbol, aver_dt$product, sep = "__")
rownames(aver_dt) <- make.unique(aver_dt$combined)

expr_mat <- aver_dt[, setdiff(colnames(aver_dt), c("symbol", "product", "cluster", "combined")), drop = FALSE]
expr_mat <- as.matrix(expr_mat)

# Keep row order according to the selected Top10 marker list.
expr_scaled <- t(scale(t(expr_mat)))
expr_scaled[is.na(expr_scaled)] <- 0

heatmap_source <- as.data.frame(expr_scaled) %>%
  rownames_to_column("marker_gene_annotation")
write.csv(heatmap_source, out_heatmap_source_csv, row.names = FALSE)

# ============================================================
# Heatmap annotation and plotting
# ============================================================
celltypes <- colnames(expr_scaled)
cell_anno <- data.frame(cell_anno = celltypes, row.names = celltypes)

celltype_col <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(unique(cell_anno$cell_anno))),
  unique(cell_anno$cell_anno)
)

col_anno <- HeatmapAnnotation(
  cell = cell_anno$cell_anno,
  col = list(cell = celltype_col),
  show_annotation_name = FALSE,
  gp = gpar(col = "white", lwd = 1.5)
)

row_anno <- rowAnnotation(
  marker = anno_text(
    rownames(expr_scaled),
    location = 0,
    just = "left",
    gp = gpar(col = "black", fontface = "italic", fontsize = 8),
    width = max_text_width(rownames(expr_scaled)) * 1.25
  )
)

heatmap_col <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))

ht <- Heatmap(
  expr_scaled,
  name = "expression",
  col = heatmap_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 60,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  rect_gp = gpar(col = "white", lwd = 1.5),
  top_annotation = col_anno
) + row_anno

png(out_heatmap_png, width = 6500, height = 8000, res = 500)
draw(ht)
grid.text(
  "Top10 marker genes across clusters",
  x = 0.5,
  y = 0.98,
  gp = gpar(fontsize = 18, fontface = "bold")
)
dev.off()
message("[output] Heatmap PNG saved: ", normalizePath(out_heatmap_png, winslash = "/", mustWork = FALSE))

pdf(out_heatmap_pdf, width = 9, height = 12)
draw(ht)
grid.text(
  "Top10 marker genes across clusters",
  x = 0.5,
  y = 0.98,
  gp = gpar(fontsize = 18, fontface = "bold")
)
dev.off()
