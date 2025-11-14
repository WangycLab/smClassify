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

set.seed(1234)
work_dir <- "D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009"
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

# ========= 0) Check files =========
if (!file.exists("stepF_seu_SCT_clustered.rds")) stop("Missing file: stepF_seu_SCT_clustered.rds")
if (!file.exists("gtf_gene_annotation.csv")) stop("Missing file: gtf_gene_annotation.csv")

seu <- readRDS("stepF_seu_SCT_clustered.rds")
gtf_wide <- fread("gtf_gene_annotation.csv", check.names = FALSE)

# ========= 1) Get markers =========
## =========================
## Step H: Add markers
## =========================
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  logfc.threshold = 1,
  min.pct = 0.05
)
write.csv(markers, "Total cells_markers_min.pct0.05.csv")
table(markers$cluster)

topN <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()
table(topN$cluster)
write.csv(topN, "Total cells_markers_top50.csv")

# Add markers to GTF
topN$symbol <- gsub("-", "_", topN$gene)
topN_gtf <- topN %>% left_join(gtf_wide, by = "symbol")
topN_gtf_subset <- unique(topN_gtf[, c(
  "p_val", "avg_log2FC", "pct.1", "pct.2",
  "p_val_adj", "cluster", "symbol", "product"
)])
write.csv(topN_gtf_subset, "Total cells_markers_top50_gtf.csv", row.names = FALSE)

## =========================
## Add all genes to GTF
## =========================
all_genes <- rownames(sce)
all_genes_df <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
all_genes_df$symbol <- gsub("-", "_", all_genes_df$gene)
all_genes_df <- all_genes_df %>% left_join(gtf_wide, by = "symbol")
all_genes_df <- unique(all_genes_df[, c("symbol", "product", "gene_id", "gene.x")])
write.csv(all_genes_df, "All_genes_symbol_gtf.csv", row.names = FALSE)
colnames(all_genes_df)

## =========================
## Top10 Marker genes
## =========================
top10 <- topN_gtf_subset %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()
top10 <- top10[!duplicated(top10$symbol), ]
top10$symbol <- gsub("_", "-", top10$symbol)
top10

genes <- unique(top10$symbol)
genes

## =========================
## Average expression across clusters
## =========================
aver_dt <- AverageExpression(
  seu,
  assays = "SCT",
  features = genes,
  group.by = "seurat_clusters",
  slot = "data"
)$SCT %>% as.data.frame()

aver_dt$symbol <- rownames(aver_dt)

# Merge annotations
top102 <- top10[, c("symbol", "product", "cluster")]
aver_dt <- left_join(aver_dt, top102, by = "symbol")
aver_dt$combined <- paste(aver_dt$symbol, aver_dt$product, sep = "__")
rownames(aver_dt) <- make.unique(aver_dt$combined)

# Extract expression matrix
expr_mat <- aver_dt[, 1:(ncol(aver_dt) - 4)]

## =========================
## Heatmap annotations
## =========================
gene_anno <- data.frame(gene_anno = top10$cluster, row.names = top10$symbol)

celltypes <- colnames(expr_mat)
cell_anno <- data.frame(cell_anno = celltypes, row.names = celltypes)

celltype_col <- setNames(
  brewer.pal(length(unique(cell_anno$cell_anno)), "Set3"),
  unique(cell_anno$cell_anno)
)
gene_anno_col <- setNames(
  brewer.pal(length(unique(gene_anno$gene_anno)), "Set1"),
  unique(gene_anno$gene_anno)
)

## =========================
## ComplexHeatmap
## =========================
expr_scaled <- t(scale(t(expr_mat)))

col_anno <- HeatmapAnnotation(
  cell = cell_anno$cell_anno,
  col  = list(cell = celltype_col),
  show_annotation_name = FALSE,
  gp   = gpar(col = 'white', lwd = 1.5)
)

row_anno <- rowAnnotation(
  foo = anno_text(
    rownames(expr_scaled),
    location = 0,
    just = "left",
    gp = gpar(col = "black", fontface = "italic"),
    width = max_text_width(rownames(expr_scaled)) * 1.25
  )
)

mycol2 <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))

ht <- Heatmap(
  expr_scaled,
  name = "expression",
  col  = mycol2,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 60,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  rect_gp = gpar(col = "white", lwd = 1.5),
  top_annotation = col_anno
) + row_anno
ht

# Plot and save heatmap
png("Top10_marker_genes_heatmap.png", width = 6500, height = 8000, res = 500)
draw(ht)
grid.text(
  "Top10 marker genes across clusters",
  x = 0.5, y = 1.5,
  gp = gpar(fontsize = 18, fontface = "bold")
)
dev.off()

pdf("Top10_marker_genes_heatmap.pdf", width = 9, height = 12)
draw(ht)
grid.text(
  "Top10 marker genes across clusters",
  x = 0.5, y = 0.98,
  gp = gpar(fontsize = 18, fontface = "bold")
)
dev.off()
