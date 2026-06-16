## ======================================================================
## Figure 6b-g. Functional subtypes and pseudotime analysis of
## Muribaculum gordoncarteri cells

library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(tibble)
library(Matrix)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(patchwork)

## ----------------------------------------------------------------------
## 1. Input and output paths
## ----------------------------------------------------------------------

base_dir   <- "Figure6"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


murigor_rds <- file.path(output_dir, "Fig6b_g_Murigor_processed.rds")

annotated_seurat_rds <- "Figure3/Output/Figure3_mouse_bacterial_functional_cluster_annotated.rds"

# ============================================================
# 3. Load annotated Seurat object
# ============================================================
seu <- readRDS(annotated_seurat_rds)

## ----------------------------------------------------------------------
## 3. Subset M. gordoncarteri cells and perform clustering
## ----------------------------------------------------------------------

Idents(seu) <- seu$species

Murigor <- subset(
  seu,
  idents = "Muribaculum gordoncarteri"
)

DefaultAssay(Murigor) <- "RNA"

Murigor <- SCTransform(
  Murigor,
  vars.to.regress = "nCount_RNA",
  verbose = FALSE
)

Murigor <- RunPCA(Murigor, verbose = FALSE)

Murigor <- FindNeighbors(
  Murigor,
  dims = 1:6
)

Murigor <- FindClusters(
  Murigor,
  resolution = 0.1
)

Murigor <- RunUMAP(
  Murigor,
  dims = 1:6
)

## ----------------------------------------------------------------------
## 4. Assign functional subtype labels
## ----------------------------------------------------------------------

Murigor$subtype_id <- paste0("Subtype ", Murigor$seurat_clusters)

subtype_labels <- c(
  "0" = "Subtype 0: nitrogen and glutamate metabolism",
  "1" = "Subtype 1: carbohydrate and fatty acid utilization",
  "2" = "Subtype 2: stress adaptation"
)

Murigor$cellsubtype <- unname(subtype_labels[as.character(Murigor$seurat_clusters)])

Murigor$cellsubtype <- factor(
  Murigor$cellsubtype,
  levels = subtype_labels
)

subtype_cols <- c(
  "Subtype 0: nitrogen and glutamate metabolism" = "#4DBBD5",
  "Subtype 1: carbohydrate and fatty acid utilization" = "#00A087",
  "Subtype 2: stress adaptation" = "#E64B35"
)

theme_umap <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

## ----------------------------------------------------------------------
## 5. Figure 6b: UMAP of functional subtypes
## ----------------------------------------------------------------------

p_fig6b <- DimPlot(
  Murigor,
  reduction = "umap",
  group.by = "cellsubtype",
  label = FALSE,
  pt.size = 0.1
) +
  scale_color_manual(values = subtype_cols) +
  labs(color = "Subtype") +
  ggtitle("M. gordoncarteri subtypes") +
  theme_umap
p_fig6b
ggsave(
  file.path(output_dir, "Fig6b_Murigor_UMAP_subtypes.pdf"),
  p_fig6b,
  width = 7,
  height = 4
)

ggsave(
  file.path(output_dir, "Fig6b_Murigor_UMAP_subtypes.png"),
  p_fig6b,
  width = 7,
  height = 4,
  dpi = 300
)

## ----------------------------------------------------------------------
## 6. Figure 6c: Relative abundance across WT and db/db tissue groups
## ----------------------------------------------------------------------

Murigor$group <- factor(
  Murigor$group,
  levels = c("Cecum-WT", "Cecum-DB", "Colon-WT", "Colon-DB", "Rectum-WT", "Rectum-DB")
)

pal_group <- c(
  "Cecum-WT"  = "#6BAED6",
  "Cecum-DB"  = "#08519C",
  "Colon-WT"  = "#74C476",
  "Colon-DB"  = "#006D2C",
  "Rectum-WT" = "#FDAE6B",
  "Rectum-DB" = "#E6550D"
)

fig6c_source <- Murigor@meta.data %>%
  dplyr::filter(!is.na(group), !is.na(cellsubtype)) %>%
  dplyr::count(group, cellsubtype, name = "cell_number") %>%
  group_by(group) %>%
  dplyr::mutate(relative_abundance = cell_number / sum(cell_number)) %>%
  ungroup()

p_fig6c <- ggplot(
  fig6c_source,
  aes(x = group, y = relative_abundance, fill = cellsubtype)
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = subtype_cols) +
  labs(x = NULL, y = "Relative abundance", fill = "Subtype") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
p_fig6c
ggsave(
  file.path(output_dir, "Fig6c_Murigor_subtype_relative_abundance.pdf"),
  p_fig6c,
  width = 6,
  height = 4
)

ggsave(
  file.path(output_dir, "Fig6c_Murigor_subtype_relative_abundance.png"),
  p_fig6c,
  width = 6,
  height = 4,
  dpi = 300
)

write.csv(
  fig6c_source,
  file.path(output_dir, "Fig6c_Murigor_subtype_relative_abundance_source_data.csv"),
  row.names = FALSE
)

## ----------------------------------------------------------------------
## 7. Figure 6d: Monocle3 pseudotime trajectory
## ----------------------------------------------------------------------

seurat_obj <- Murigor

expr_matrix <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  slot = "counts"
)

cell_metadata <- seurat_obj@meta.data

gene_annotation <- data.frame(
  gene_short_name = rownames(expr_matrix),
  row.names = rownames(expr_matrix)
)

cds <- new_cell_data_set(
  expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

reducedDims(cds)$UMAP <- Embeddings(seurat_obj, reduction = "umap")

cds <- cluster_cells(
  cds,
  reduction_method = "UMAP"
)

cds <- learn_graph(
  cds,
  use_partition = TRUE
)

## Select the principal graph node most enriched for cecum cells as root,
## consistent with a cecum-to-colon-to-rectum trajectory.
get_root_pr_nodes <- function(cds, root_region = "Cecum") {
  cell_ids <- colnames(cds)[colData(cds)$region == root_region]
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])
  
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
  
  root_pr_nodes
}

root_pr_nodes <- get_root_pr_nodes(cds, root_region = "Cecum")

cds <- order_cells(
  cds,
  root_pr_nodes = root_pr_nodes
)

p_fig6d <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  cell_size = 0.5,
  alpha = 0.5
) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
p_fig6d
ggsave(
  file.path(output_dir, "Fig6d_Murigor_monocle3_pseudotime.pdf"),
  p_fig6d,
  width = 6,
  height = 4
)

ggsave(
  file.path(output_dir, "Fig6d_Murigor_monocle3_pseudotime.png"),
  p_fig6d,
  width = 6,
  height = 4,
  dpi = 300
)

## Add pseudotime back to Seurat metadata.
Murigor$pseudotime <- pseudotime(cds)[colnames(Murigor)]

## ----------------------------------------------------------------------
## 8. Figure 6e: Density distributions along pseudotime
## ----------------------------------------------------------------------

meta_df <- Murigor@meta.data %>%
  rownames_to_column("cell_id") %>%
  filter(!is.na(pseudotime))

phenotype_cols <- c(
  "WT" = "#3498DB",
  "DB" = "#E74C3C",
  "T2DM" = "#E74C3C",
  "db/db" = "#E74C3C"
)

region_cols <- c(
  "Cecum" = "#1F77B4",
  "Colon" = "#2CA02C",
  "Rectum" = "#FF7F0E"
)

p_fig6e_phenotype <- ggplot(
  meta_df,
  aes(x = pseudotime, fill = phenotype)
) +
  geom_density(alpha = 0.65, color = NA) +
  facet_wrap(~phenotype, ncol = 1) +
  scale_fill_manual(values = phenotype_cols) +
  labs(x = "Pseudotime", y = "Cell density") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )

p_fig6e_region <- ggplot(
  meta_df,
  aes(x = pseudotime, fill = region)
) +
  geom_density(alpha = 0.65, color = NA) +
  facet_wrap(~region, ncol = 1) +
  scale_fill_manual(values = region_cols) +
  labs(x = "Pseudotime", y = "Cell density") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )

ggsave(
  file.path(output_dir, "Fig6e_Murigor_pseudotime_density_phenotype.pdf"),
  p_fig6e_phenotype,
  width = 4,
  height = 3.5
)

ggsave(
  file.path(output_dir, "Fig6e_Murigor_pseudotime_density_phenotype.png"),
  p_fig6e_phenotype,
  width = 4,
  height = 3.5,
  dpi = 600
)

ggsave(
  file.path(output_dir, "Fig6e_Murigor_pseudotime_density_region.pdf"),
  p_fig6e_region,
  width = 4,
  height = 3.5
)

ggsave(
  file.path(output_dir, "Fig6e_Murigor_pseudotime_density_region.png"),
  p_fig6e_region,
  width = 4,
  height = 3.5,
  dpi = 600
)

write.csv(
  meta_df,
  file.path(output_dir, "Fig6d_e_Murigor_pseudotime_source_data.csv"),
  row.names = FALSE
)

## ----------------------------------------------------------------------
## 9. Figure 6f: Heatmap of pseudotime-associated genes
## ----------------------------------------------------------------------

deg_pseudotime <- graph_test(
  cds,
  neighbor_graph = "principal_graph",
  cores = 4
)

deg_pseudotime <- deg_pseudotime %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  arrange(q_value)

write.csv(
  deg_pseudotime,
  file.path(output_dir, "Fig6f_Murigor_pseudotime_graph_test_all_genes.csv"),
  row.names = FALSE
)

sig_genes <- deg_pseudotime %>%
  filter(!is.na(gene_id), !grepl("^MGYG", gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  slice_head(n = 100) %>%
  pull(gene_id)

sig_genes <- intersect(sig_genes, rownames(Murigor[["SCT"]]))

cell_order <- meta_df %>%
  arrange(pseudotime) %>%
  pull(cell_id)

cell_order <- intersect(cell_order, colnames(Murigor))

log_expr <- GetAssayData(
  Murigor,
  assay = "SCT",
  slot = "data"
)[sig_genes, cell_order, drop = FALSE]

log_expr <- as.matrix(log_expr)
log_expr_scaled <- t(scale(t(log_expr)))
log_expr_scaled[is.na(log_expr_scaled)] <- 0

row_clust <- hclust(
  dist(log_expr_scaled),
  method = "complete"
)

gene_modules <- cutree(row_clust, k = 4)

gene_info <- data.frame(
  gene_id = rownames(log_expr_scaled),
  Module = paste0("Module", gene_modules),
  stringsAsFactors = FALSE
)

gene_info$sd <- apply(log_expr_scaled, 1, sd)

top_genes_by_module <- gene_info %>%
  group_by(Module) %>%
  slice_max(order_by = sd, n = 10, with_ties = FALSE) %>%
  ungroup()

row_anno <- data.frame(
  Module = gene_info$Module,
  row.names = gene_info$gene_id
)

module_levels <- sort(unique(gene_info$Module))
module_colors <- RColorBrewer::brewer.pal(length(module_levels), "Set2")
names(module_colors) <- module_levels

anno_colors <- list(Module = module_colors)

labels_vec <- ifelse(
  gene_info$gene_id %in% top_genes_by_module$gene_id,
  gene_info$gene_id,
  ""
)

heatmap_cols <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
)(100)

heatmap_breaks <- seq(-2, 2, length.out = 101)

pdf(
  file.path(output_dir, "Fig6f_Murigor_pseudotime_heatmap.pdf"),
  width = 4.5,
  height = 8
)

pheatmap(
  log_expr_scaled,
  cluster_rows = row_clust,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = labels_vec,
  show_colnames = FALSE,
  color = heatmap_cols,
  breaks = heatmap_breaks,
  annotation_row = row_anno,
  annotation_colors = anno_colors,
  main = "Pseudotime-associated genes",
  use_raster = FALSE
)

dev.off()

png(
  file.path(output_dir, "Fig6f_Murigor_pseudotime_heatmap.png"),
  width = 4.5,
  height = 8,
  units = "in",
  res = 600
)

pheatmap(
  log_expr_scaled,
  cluster_rows = row_clust,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = labels_vec,
  show_colnames = FALSE,
  color = heatmap_cols,
  breaks = heatmap_breaks,
  annotation_row = row_anno,
  annotation_colors = anno_colors,
  main = "Pseudotime-associated genes",
  use_raster = FALSE
)

dev.off()

write.csv(
  gene_info,
  file.path(output_dir, "Fig6f_Murigor_pseudotime_gene_modules.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(log_expr_scaled) %>%
    rownames_to_column("gene_id"),
  file.path(output_dir, "Fig6f_Murigor_pseudotime_heatmap_source_data.csv"),
  row.names = FALSE
)

## ----------------------------------------------------------------------
## 10. Figure 6g: Selected gene trends along pseudotime
## ----------------------------------------------------------------------

target_genes <- c(
  "fabF",
  "gdhA", "glnA", "gadB", "gadC",
  "sodB", "rbr",
  "ureG"
)

target_genes <- intersect(target_genes, rownames(cds))

if (length(target_genes) == 0) {
  stop("None of the selected Figure 6g genes were found in cds.")
}

rowData(cds)$gene_short_name <- rownames(cds)

p_fig6g <- plot_genes_in_pseudotime(
  cds[target_genes, ],
  ncol = 4,
  min_expr = 0.5,
  color_cells_by = "pseudotime",
  cell_size = 1
) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave(
  file.path(output_dir, "Fig6g_Murigor_selected_gene_pseudotime_trends.pdf"),
  p_fig6g,
  width = 14,
  height = 5
)

ggsave(
  file.path(output_dir, "Fig6g_Murigor_selected_gene_pseudotime_trends.png"),
  p_fig6g,
  width = 14,
  height = 5,
  dpi = 300
)

fig6g_expr <- FetchData(
  Murigor,
  vars = c("pseudotime", "region", "phenotype", "group", "cellsubtype", target_genes)
) %>%
  rownames_to_column("cell_id")

write.csv(
  fig6g_expr,
  file.path(output_dir, "Fig6g_Murigor_selected_gene_trends_source_data.csv"),
  row.names = FALSE
)

saveRDS(
  Murigor,
  murigor_rds
)