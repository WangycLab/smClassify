############################################################
## Figure 5f,g and supplementary Duncaniella analyses
##

# ============================================================
# 1. Load packages
# ============================================================
library(Seurat)
library(dplyr)
library(ggplot2)


# ============================================================
# 2. Input and output paths
# ============================================================
base_dir   <- "Figure5"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

annotated_seurat_rds <- "Figure3/Output/Figure3_mouse_bacterial_functional_cluster_annotated.rds"
gtf_wide_csv         <- file.path(input_dir, "gtf_wide.csv")
dun_type_csv         <- file.path(input_dir, "Dun_type.csv")

# ============================================================
# 3. Load annotated Seurat object
# ============================================================
seu <- readRDS(annotated_seurat_rds)

# ============================================================
# 4. Define colours
# ============================================================
cell_type_cols <- c(
  "#1F77B4", "#2CA02C", "#FF7F0E", "#6A5ACD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22",
  "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1",
  "#3CB371", "#9ACD32", "#8B0000", "#FFD700",
  "#DC143C", "#228B22", "#FF6347", "#483D8B",
  "#BDB76B", "#20B2AA", "#FF1493", "#FF4500",
  "#32CD32", "#3E8E41", "#20B2AA"
)

cell_type_cols2 <- c(
  "#1F77B4", "#8C564B", "#E377C2", "#7F7F7F",
  "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8",
  "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000",
  "#FFD700", "#DC143C", "#228B22", "#FF6347",
  "#483D8B", "#BDB76B", "#FF1493", "#FF4500",
  "#32CD32", "#3E8E41", "#20B2AA"
)

pal_group <- c(
  "Cecum-WT"  = "#6baed6",
  "Cecum-DB"  = "#08519c",
  "Colon-WT"  = "#74c476",
  "Colon-DB"  = "#006d2c",
  "Rectum-WT" = "#fdae6b",
  "Rectum-DB" = "#e6550d"
)


# ============================================================
# 8. Re-cluster Duncaniella cells
# ============================================================
Dun_seu <- SCTransform(
  Dun_seu,
  verbose = FALSE,
  vars.to.regress = "nCount_RNA"
)

Dun_seu <- RunPCA(Dun_seu, verbose = FALSE)
Dun_seu <- FindNeighbors(Dun_seu, dims = 1:6)
Dun_seu <- FindClusters(Dun_seu, resolution = 0.1)
Dun_seu <- RunUMAP(Dun_seu, dims = 1:6)

saveRDS(
  Dun_seu,
  file = file.path(output_dir, "Figure5fg_Duncaniella_reclustered.rds")
)

# ============================================================
# 9. Supplementary UMAPs before subtype annotation
# ============================================================
p_cluster <- DimPlot(
  Dun_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "seurat_clusters"
) +
  scale_color_manual(values = cell_type_cols) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )
p_cluster 
ggsave(
  file.path(output_dir, "Supplementary_Duncaniella_UMAP_by_cluster.png"),
  p_cluster,
  width = 4,
  height = 3,
  dpi = 300
)

p_region <- DimPlot(
  Dun_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "region"
) +
  scale_color_manual(values = c("Cecum" = "#1f77b4", "Colon" = "#2ca02c", "Rectum" = "#ff7f0e")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )

ggsave(
  file.path(output_dir, "Supplementary_Duncaniella_UMAP_by_region.png"),
  p_region,
  width = 4,
  height = 3,
  dpi = 300
)

p_phenotype <- DimPlot(
  Dun_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "phenotype"
) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )

ggsave(
  file.path(output_dir, "Supplementary_Duncaniella_UMAP_by_phenotype.png"),
  p_phenotype,
  width = 4,
  height = 3,
  dpi = 300
)

p_species <- DimPlot(
  Dun_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "species"
) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )

ggsave(
  file.path(output_dir, "Supplementary_Duncaniella_UMAP_by_species.png"),
  p_species,
  width = 6,
  height = 3,
  dpi = 300
)

# ============================================================
# 10. Supplementary marker analysis
# ============================================================
Idents(Dun_seu) <- Dun_seu$seurat_clusters

markers <- FindAllMarkers(
  object = Dun_seu,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.15
)

write.csv(
  markers,
  file = file.path(output_dir, "Supplementary_Duncaniella_markers_minpct0.05_logfc0.15.csv"),
  row.names = FALSE
)

topN <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()

write.csv(
  topN,
  file = file.path(output_dir, "Supplementary_Duncaniella_markers_top50.csv"),
  row.names = FALSE
)

gtf_wide <- fread(gtf_wide_csv)

topN$symbol <- gsub("-", "_", topN$gene)
topN_gtf <- topN %>% left_join(gtf_wide, by = "symbol")

topN_gtf_subset <- unique(
  topN_gtf[, c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "symbol", "product")]
)

write.csv(
  topN_gtf_subset,
  file = file.path(output_dir, "Supplementary_Duncaniella_markers_top50_gtf.csv"),
  row.names = FALSE
)

top5_nonMGYG <- topN %>%
  filter(!grepl("^MGYG", gene)) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  ungroup()

write.csv(
  top5_nonMGYG,
  file = file.path(output_dir, "Supplementary_Duncaniella_top5_nonMGYG_markers.csv"),
  row.names = FALSE
)

genes <- unique(top5_nonMGYG$gene)

p_marker <- DotPlot(Dun_seu, features = genes) +
  scale_color_gradient(low = "#2ECC71", high = "#8E44AD") +
  theme_minimal(base_size = 13) +
  scale_size_continuous(range = c(1, 8)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "Functional Subtypes", title = "DotPlot of Top Marker Genes")

ggsave(
  file.path(output_dir, "Supplementary_Duncaniella_top_marker_gene_dotplot.pdf"),
  p_marker,
  width = 12,
  height = 3.5
)

ggsave(
  file.path(output_dir, "Supplementary_Duncaniella_top_marker_gene_dotplot.png"),
  p_marker,
  width = 12,
  height = 3.5,
  dpi = 300
)

# ============================================================
# 11. Add Duncaniella subpopulation annotation
# ============================================================
dun_type <- read.csv(dun_type_csv)

new.cluster.ids <- as.character(dun_type[, 2])

Dun_seu@active.ident <- as.factor(Dun_seu$seurat_clusters)
names(new.cluster.ids) <- levels(Dun_seu)

Dun_seu <- RenameIdents(Dun_seu, new.cluster.ids)
Dun_seu$cellsubtype <- Dun_seu@active.ident

write.csv(
  as.data.frame(table(Dun_seu$cellsubtype)),
  file = file.path(output_dir, "Figure5f_Duncaniella_subpopulation_cell_counts.csv"),
  row.names = FALSE
)

saveRDS(
  Dun_seu,
  file = file.path(output_dir, "Figure5fg_Duncaniella_subpopulation_annotated.rds")
)

# ============================================================
# 12. Figure 5f: Duncaniella UMAP colored by subpopulation
# ============================================================
p_5f <- DimPlot(
  Dun_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "cellsubtype"
) +
  scale_color_manual(values = cell_type_cols2) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )

ggsave(
  file.path(output_dir, "Figure5f_Duncaniella_subpopulation_UMAP.png"),
  p_5f,
  width = 8,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(output_dir, "Figure5f_Duncaniella_subpopulation_UMAP.pdf"),
  p_5f,
  width = 8,
  height = 4
)

# ============================================================
# 13. Figure 5g: Duncaniella subpopulation proportions
# ============================================================
meta_all <- seu@meta.data
meta_dun <- Dun_seu@meta.data

group_levels <- c(
  "Cecum-WT", "Cecum-DB",
  "Colon-WT", "Colon-DB",
  "Rectum-WT", "Rectum-DB"
)

meta_all$group <- factor(meta_all$group, levels = group_levels)
meta_dun$group <- factor(meta_dun$group, levels = group_levels)

group_total <- meta_all %>%
  group_by(group) %>%
  summarise(total_cells = n(), .groups = "drop")

dun_counts <- meta_dun %>%
  group_by(cellsubtype, group) %>%
  summarise(dun_cells = n(), .groups = "drop")

df_ratio <- dun_counts %>%
  left_join(group_total, by = "group") %>%
  mutate(freq = dun_cells / total_cells)

write.csv(
  df_ratio,
  file = file.path(output_dir, "Figure5g_Duncaniella_subpopulation_fraction_source_data.csv"),
  row.names = FALSE
)

p_5g <- ggplot(df_ratio, aes(x = group, y = freq, fill = group)) +
  geom_bar(stat = "identity", width = 0.6) +
  facet_wrap(~ cellsubtype, scales = "free_y", ncol = 5) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    y = "Proportion of total cells in group",
    x = "Group",
    title = "Duncaniella subproportions across groups"
  ) +
  scale_fill_manual(values = pal_group)

ggsave(
  file.path(output_dir, "Figure5g_Duncaniella_subpopulation_fraction_by_group.png"),
  p_5g,
  width = 18,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(output_dir, "Figure5g_Duncaniella_subpopulation_fraction_by_group.pdf"),
  p_5g,
  width = 18,
  height = 6
)

# ============================================================
# 14. Supplementary Duncaniella species proportions
# ============================================================
dun_species_counts <- meta_dun %>%
  group_by(species, group) %>%
  summarise(dun_cells = n(), .groups = "drop")

df_species_ratio <- dun_species_counts %>%
  left_join(group_total, by = "group") %>%
  mutate(freq = dun_cells / total_cells)

write.csv(
  df_species_ratio,
  file = file.path(output_dir, "Supplementary_Duncaniella_species_fraction_source_data.csv"),
  row.names = FALSE
)

p_species_fraction <- ggplot(df_species_ratio, aes(x = group, y = freq, fill = group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ species, scales = "free_y", ncol = 5) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    y = "Proportion of total cells in group",
    x = "Group",
    title = "Duncaniella species proportions across groups"
  ) +
  scale_fill_manual(values = pal_group) +
  theme(
    strip.background = element_rect(fill = "#91D1C2", color = "black"),
    strip.text = element_text(color = "black", face = "bold")
  )

ggsave(
  file.path(output_dir, "Supplementary_Duncaniella_species_fraction_by_group.png"),
  p_species_fraction,
  width = 16,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(output_dir, "Supplementary_Duncaniella_species_fraction_by_group.pdf"),
  p_species_fraction,
  width = 16,
  height = 4
)
