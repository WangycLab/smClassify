suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(future)
  library(tibble)
  library(paletteer)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

set.seed(1234)   # Fix random seed for reproducibility
work_dir <- "D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009"
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

# Global palette used in multiple plots
color21 <- c(paletteer::paletteer_d("ggthemes::Classic_20", n = 20), "black")

## =========================
## Helper functions
## =========================

.blankish <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; !nzchar(trimws(x)) }
.norm_id  <- function(x) toupper(gsub("_","-", as.character(x), fixed = TRUE))

# Build rRNA/tmRNA gene set
build_rrna_set <- function(attr_df) {
  if ("is_rRNA" %in% names(attr_df) || "is_tmRNA" %in% names(attr_df)) {
    isR <- attr_df$is_rRNA %in% TRUE
    isT <- attr_df$is_tmRNA %in% TRUE
    return(.norm_id(unique(na.omit(attr_df$locus_tag[isR | isT]))))
  }
  prod <- as.character(attr_df$product)
  mask <- (!is.na(prod)) &
    (grepl("(?i)(ribosomal\\s*RNA|[0-9]+S\\s*rRNA|RNASE\\s*P\\s*RNA)", prod, perl = TRUE) |
       grepl("(?i)tmRNA|transfer-messenger", prod, perl = TRUE)) &
    (!grepl("(?i)tRNA|riboswitch", prod, perl = TRUE))
  .norm_id(unique(na.omit(attr_df$locus_tag[mask])))
}

# Align a wide GTF-like annotation table to the expression matrix row names
align_gtf_to_genes <- function(gtf_wide, genes){
  setDT(gtf_wide)
  genes_norm <- .norm_id(genes)
  key_order <- c("locus_tag","ID","gene_id","Name")
  gtf_wide[, gene_id := as.character(gene_id)]
  for (k in key_order) {
    to_fill <- .blankish(gtf_wide$gene_id)
    cand    <- .norm_id(gtf_wide[[k]])
    hit     <- to_fill & !.blankish(cand) & (cand %in% genes_norm)
    gtf_wide[hit, gene_id := cand[hit]]
  }
  keep_cols <- intersect(c("gene_id","symbol","ko","product"), names(gtf_wide))
  unique(gtf_wide[! .blankish(gene_id), ..keep_cols])
}

# Fast sparse row aggregation by group
fast_rowsum_by_group <- function(M, gvec, min_total = 1L, verbose = TRUE){
  if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  stopifnot(length(gvec) == nrow(M))
  rs <- Matrix::rowSums(M)
  keep <- !is.na(gvec) & rs >= min_total
  M2 <- M[keep,,drop = FALSE]; grp <- gvec[keep]
  f <- factor(grp, levels = unique(grp))
  G <- Matrix::sparseMatrix(i = as.integer(f), j = seq_along(f), x = 1L,
                            dims = c(nlevels(f), length(f)))
  out <- G %*% M2
  rownames(out) <- levels(f)
  if (verbose) {
    message(sprintf("Aggregation complete: %d groups ?? %d cells; nnz=%d",
                    nrow(out), ncol(out), Matrix::nnzero(out)))
  }
  out
}

## =========================
## Main workflow
## =========================

## Step A: Remove rRNA/tmRNA features
sce <- readRDS("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009/sce_bac_annotated.rds")
attr_df <- fread("mouse_annot_feature_attributes_rrna_flags.csv")
rrna_ids <- build_rrna_set(attr_df)
sce <- sce[!(rownames(sce) %in% rrna_ids), ]
saveRDS(sce, "mouse_20250927_sce_no_rrna.rds")

## Step B: Align GTF-wide annotation
gtf_wide <- fread("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20250910/gtf_wide.csv")
genes_nr <- rownames(GetAssayData(sce, "RNA", "counts"))
sym_map  <- align_gtf_to_genes(gtf_wide, genes_nr)
sym_map  <- sym_map[!(tolower(product) == "hypothetical protein"), ]
fwrite(sym_map, "mouse_20250927_gtf_symbol_map_filtered.tsv", sep = "\t")

## Step C: Aggregate counts by gene symbol
sym_map <- sym_map %>%
  mutate(symbol = ifelse(grepl("^MGYG", symbol), symbol, gsub("[-_][0-9]+$", "", symbol)))
gvec <- setNames(sym_map$symbol, sym_map$gene_id)[.norm_id(genes_nr)]
M_symbol <- fast_rowsum_by_group(GetAssayData(sce, "RNA", "counts"), gvec, min_total = 5)
saveRDS(M_symbol, "mouse_20250927_counts_aggregated_by_symbol.rds")

## Step D: Build Seurat object and carry over metadata
seu <- CreateSeuratObject(counts = M_symbol, meta.data = sce@meta.data[colnames(M_symbol), ])
Idents(seu) <- seu$species
saveRDS(seu, "mouse_20250927_seurat_symbol_raw.rds")

## Step E: Downsampling & light QC
# Keep features expressed (>0) in >= 5 cells
keep_features <- names(which(Matrix::rowSums(GetAssayData(seu, "RNA", "counts") > 0) >= 5))
seu <- subset(seu, features = keep_features)

# 7000 cells per sample (by orig.ident)
cells_keep <- seu@meta.data %>% rownames_to_column("cell") %>%
  group_by(orig.ident) %>%
  group_modify(~ slice_sample(.x, n = min(nrow(.x), 7000))) %>%
  pull(cell)
seu <- subset(seu, cells = cells_keep)

# Basic nCount_RNA range filter
seu <- subset(seu, subset = nCount_RNA > 30 & nCount_RNA < 1000)
saveRDS(seu, "mouse_20250927_seurat_symbol_qc.rds")

## Step F: SCTransform, neighbors, clustering, UMAP
plan(sequential)
options(future.globals.maxSize = 32 * 1024^3)
seu <- SCTransform(
  seu, vst.flavor = "v2", vars.to.regress = "nCount_RNA",
  variable.features.n = 5000, conserve.memory = TRUE,
  method = if (requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else "poisson"
)
seu <- RunPCA(seu, npcs = 50, assay = "SCT")
seu <- FindNeighbors(seu, dims = 1:6)
seu <- FindClusters(seu, resolution = 0.1)
seu <- RunUMAP(seu, dims = 1:6)

DimPlot(
  seu, reduction = "umap", group.by = "seurat_clusters",
  pt.size = .1, raster = FALSE, cols = color21
)
saveRDS(seu, "mouse_20250927_seurat_sct_umap_clusters.rds")

## Step G: UMAP colored by top-20 species
top_names <- seu@meta.data %>% as_tibble() %>%
  count(species, orig.ident, name = "n") %>%
  group_by(species) %>%
  summarise(total_count = sum(n), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 20) %>%
  pull(species)

seu@meta.data <- seu@meta.data %>%
  mutate(top20_species = ifelse(species %in% top_names, species, "Others"))

top20_species <- DimPlot(
  seu, reduction = "umap", group.by = "top20_species",
  pt.size = .1, raster = FALSE, cols = color21
)
top20_species
ggsave("mouse_20250927_umap_top20_species.png", top20_species, width = 12, height = 5, dpi = 600)

## Step H: Marker detection and functional annotation join
markers <- FindAllMarkers(
  seu, only.pos = TRUE,
  logfc.threshold = 1,
  min.pct = 0.05
)
write.csv(markers, "mouse_20250927_markers_minpct0.05_logfc1.csv", row.names = FALSE)
table(markers$cluster)

topN <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()
table(topN$cluster)
write.csv(topN, "mouse_20250927_markers_top50_per_cluster.csv", row.names = FALSE)

# Join marker genes to GTF-wide table
topN$symbol <- gsub("-", "_", topN$gene)
topN_gtf <- topN %>% left_join(gtf_wide, by = "symbol")
topN_gtf_subset <- unique(
  topN_gtf[, c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","symbol","product")]
)
write.csv(topN_gtf_subset, "mouse_20250927_markers_top50_with_gtf.csv", row.names = FALSE)

# Export all genes with symbol/product mapping
all_genes <- rownames(sce)
all_genes_df <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
all_genes_df$symbol <- gsub("-", "_", all_genes_df$gene)
all_genes_df <- all_genes_df %>% left_join(gtf_wide, by = "symbol")
all_genes_df <- unique(all_genes_df[, c("symbol","product","gene_id","gene.x")])
write.csv(all_genes_df, "mouse_20250927_all_genes_symbol_gtf_map.csv", row.names = FALSE)
colnames(all_genes_df)

# Select top-10 per cluster for heatmap labels
top10 <- topN_gtf_subset %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()
top10 <- top10[!duplicated(top10$symbol), ]
top10$symbol <- gsub("_", "-", top10$symbol)

genes <- unique(top10$symbol)
genes

# Also collect top non-MGYG genes (optional showcase list used for averaging)
top5_nonMGYG <- markers %>%
  filter(!grepl("^MGYG", gene)) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20) %>%
  ungroup()
genes <- unique(top5_nonMGYG$gene)
genes

# Average expression per cluster (SCT assay, log-normalized "data" slot)
aver_dt <- AverageExpression(
  seu,
  assays = "SCT",
  features = genes,
  group.by = "seurat_clusters",
  slot = "data"
)$SCT %>% as.data.frame()
aver_dt$symbol <- rownames(aver_dt)

# Merge product info for composite row labels
top102 <- top10[, c("symbol", "product", "cluster")]
aver_dt <- left_join(aver_dt, top102, by = "symbol")
aver_dt$combined <- paste(aver_dt$symbol, aver_dt$product, sep = "__")
rownames(aver_dt) <- make.unique(aver_dt$combined)

# Keep only the expression matrix part
expr_mat <- aver_dt[, 1:(ncol(aver_dt) - 4)]

## ============ Annotations for ComplexHeatmap ============
# Row annotation: marker cluster
gene_anno <- data.frame(gene_anno = top10$cluster, row.names = top10$symbol)

# Column annotation: cluster IDs
celltypes <- colnames(expr_mat)
cell_anno <- data.frame(cell_anno = celltypes, row.names = celltypes)

# Colors for annotations
celltype_col <- setNames(brewer.pal(length(unique(cell_anno$cell_anno)), "Set3"),
                         unique(cell_anno$cell_anno))
gene_anno_col <- setNames(brewer.pal(length(unique(gene_anno$gene_anno)), "Set1"),
                          unique(gene_anno$gene_anno))

## ============ ComplexHeatmap ============
# Z-score by gene (row-wise)
expr_scaled <- t(scale(t(expr_mat)))

# Column annotation (clusters)
col_anno <- HeatmapAnnotation(
  cell = cell_anno$cell_anno,
  col  = list(cell = celltype_col),
  show_annotation_name = FALSE,
  gp   = gpar(col = 'white', lwd = 1.5)
)

# Row-side labels
row_anno <- rowAnnotation(
  foo = anno_text(
    rownames(expr_scaled),
    location = 0, just = "left",
    gp = gpar(col = "black", fontface = "italic"),
    width = max_text_width(rownames(expr_scaled)) * 1.25
  )
)

# Heatmap color scale
mycol2 <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))

# Compose heatmap
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

# Export heatmap (PNG/PDF)
png("mouse_20250927_heatmap_expr_selected_genes.png", width = 6500, height = 5000, res = 500)
draw(ht)
grid.text("COE1butyricum hallii", x = 0.5, y = 0.98, gp = gpar(fontsize = 18, fontface = "bold"))
dev.off()

pdf("mouse_20250927_heatmap_expr_selected_genes.pdf", width = 9, height = 8)
draw(ht)
grid.text("COE1butyricum hallii", x = 0.5, y = 0.98, gp = gpar(fontsize = 18, fontface = "bold"))
dev.off()
