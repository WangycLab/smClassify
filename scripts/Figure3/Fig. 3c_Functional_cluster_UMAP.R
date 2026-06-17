# ============================================================
# Figure 3c functional-cluster UMAP

# ============================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(scales)


set.seed(1234)


base_dir   <- "Figure3"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

seurat_rds <- file.path(input_dir, "Mouse_bacterial_SCT_clustered_seurat.rds")
anno_csv   <- file.path(input_dir, "Cluster_functional_annotation.csv")

out_png    <- file.path(output_dir, "Figure3c_UMAP_functional_cluster.png")
out_pdf    <- file.path(output_dir, "Figure3c_UMAP_functional_cluster.pdf")
out_source <- file.path(output_dir, "Figure3c_functional_cluster_source_data.csv")
out_rds    <- file.path(output_dir, "Figure3_mouse_bacterial_functional_cluster_annotated.rds")


seu <- readRDS(seurat_rds)
anno_df <- read.csv(anno_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

colnames(anno_df)[1:2] <- c("cluster_id", "functional_cluster")

anno_df <- anno_df %>%
  mutate(
    cluster_id_clean = gsub("^C", "", as.character(cluster_id)),
    functional_cluster = as.character(functional_cluster)
  )

cluster_map <- setNames(
  anno_df$functional_cluster,
  anno_df$cluster_id_clean
)

seu$celltype <- unname(
  cluster_map[as.character(seu$seurat_clusters)]
)


functional_levels <- anno_df$functional_cluster[match(
  sort(as.numeric(anno_df$cluster_id_clean)),
  as.numeric(anno_df$cluster_id_clean)
)]
functional_levels <- unique(functional_levels)

seu$celltype <- factor(seu$celltype, levels = functional_levels)
Seurat::Idents(seu) <- seu$celltype

print(table(seu$celltype, useNA = "ifany"))

saveRDS(seu, out_rds)

# ============================================================
# Extract UMAP source data
# ============================================================
umap_df <- Seurat::Embeddings(seu, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell")

colnames(umap_df)[1:3] <- c("cell", "UMAP_1", "UMAP_2")

meta_df <- seu@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, seurat_clusters, celltype, dplyr::everything())

plot_df <- left_join(umap_df, meta_df, by = "cell")

# ============================================================
# Plot Figure 3c
# ============================================================
celltype_colors<- c("#0071B1","#AF79A0", "#E69F00", "#4E79A7", "#ECC847",
  "#58A04F", "#999999", "#D55E00")
celltype_colors <- celltype_colors[seq_along(functional_levels)]
names(celltype_colors) <- functional_levels
celltype_colors
theme_umap <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.title = element_text(size = base_size * 0.85, face = "bold"),
      legend.text = element_text(size = base_size * 0.75),
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5)
    )
}

p_fig3c <- ggplot(
  plot_df,
  aes(x = UMAP_1, y = UMAP_2, color = celltype)
) +
  geom_point(size = 0.3, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = celltype_cols ) +
  labs(
    title = "UMAP by functional cluster identity",
    x = "UMAP_1",
    y = "UMAP_2",
    color = "Functional cluster"
  ) +
  theme_umap(base_size = 14) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p_fig3c)

ggsave(out_png, p_fig3c, width = 6.2, height = 4.5, dpi = 600)
ggsave(out_pdf, p_fig3c, width = 6.2, height = 4.5)

