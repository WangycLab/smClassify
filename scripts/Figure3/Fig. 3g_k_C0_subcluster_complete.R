# ============================================================
# Figure 3g-k. Subclustering and functional module analysis of cluster 0
# ============================================================

library(Seurat)
library(tidyverse)
library(Matrix)
library(data.table)
library(ggrepel)
library(RColorBrewer)
library(gplots)
library(rtracklayer)
library(scales)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(paletteer)
library(plot1cell)
library(grid)

# ============================================================
# 1. Input and output paths
# ============================================================

base_dir   <- "Figure3"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")
out_dir    <- file.path(output_dir, "Figure3g_k_C0_subcluster")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_file <- file.path(input_dir, "Figure3_mouse_bacterial_functional_cluster_annotated.rds")
seu <- readRDS(rds_file)

# ============================================================
# 2. Color palette
# ============================================================

nature_palette2 <- c(
  "#3C5488", "#F39B7F", "#8491B4", "#91D1C2",
  "#F5A0A1", "#C2B5D8", "#FFB6C1", "#7E6148",
  "#E64B35", "#4DBBD5", "#00A087", "#B09C85"
)

expand_palette <- grDevices::colorRampPalette(nature_palette2)
cell_type_cols <- expand_palette(80)

mycol2 <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#0da9ce", "white", "#e74a32")
)

# ============================================================
# 3. Subclustering of cluster 0
# ============================================================

seu@active.ident <- as.factor(seu$seurat_clusters)

glycan_seu <- subset(seu, idents = "0")

glycan_seu <- SCTransform(
  glycan_seu,
  verbose = FALSE,
  vars.to.regress = "nCount_RNA"
)

glycan_seu <- RunPCA(glycan_seu, verbose = FALSE)

p_elbow <- ElbowPlot(glycan_seu)

ggsave(file.path(out_dir, "C0_glycan_ElbowPlot.png"), p_elbow, width = 6, height = 5, dpi = 600)
ggsave(file.path(out_dir, "C0_glycan_ElbowPlot.pdf"), p_elbow, width = 6, height = 5)

glycan_seu <- FindNeighbors(glycan_seu, dims = 1:4)
glycan_seu <- FindClusters(glycan_seu, resolution = 0.15)
glycan_seu <- RunUMAP(glycan_seu, dims = 1:4)

saveRDS(glycan_seu, file.path(out_dir, "Figure3g_k_C0_subcluster_seurat.rds"))

# ============================================================
# 4. UMAP plots
# ============================================================

p1 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 0.5
) +
  scale_color_manual(values = nature_palette2) +
  theme_minimal()

ggsave(file.path(out_dir, "C0_glycan_UMAP_Clusters.png"), p1, width = 8, height = 7, dpi = 1200)

p2 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  group.by = "group",
  pt.size = 1
) +
  scale_color_manual(values = cell_type_cols) +
  theme_minimal()

ggsave(file.path(out_dir, "C0_glycan_UMAP_Group.png"), p2, width = 8, height = 8, dpi = 1200)

p3 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  group.by = "phenotype",
  pt.size = 0.5
) +
  scale_color_manual(values = c("#E64B35", "#4DBBD5")) +
  theme_minimal()

ggsave(file.path(out_dir, "C0_glycan_UMAP_Phenotype.png"), p3, width = 8, height = 6, dpi = 300)

p4 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  group.by = "loc",
  pt.size = 0.5
) +
  scale_color_manual(values = c("#1F77B4", "#2CA02C", "#FF7F0E")) +
  theme_minimal()

ggsave(file.path(out_dir, "C0_glycan_UMAP_Location.png"), p4, width = 8, height = 6, dpi = 300)

# ============================================================
# 5. Gene sets + module scoring
# ============================================================

genesets <- list(
  Plant_PUL = c("susC","susD","susE","susF","GH13","GH32","GH28","GH10","GH11"),
  Butyrate_core = c("thl","hbd","crt","bcd","etfA","etfB"),
  Vitamin_B12_MK = c("cobA","cobB","cobC","cobD","cobQ","cobS"),
  Antioxidant = c("katA","katE","sodA","sodB","ahpC","ahpF","trxA","trxB"),
  Mucin_GH = c("nanA","nanE","nanK","nagA","nagB","GH33","GH29"),
  Bile_total = c("cbh","bsh","tolC","acrA","acrB","emrA","emrB"),
  SCFA_total = c("thl","hbd","crt","bcd","ptb","buk","pta","ackA"),
  Flagella = c("fliC","fliD","fliE","fliF","fliG","fliH","fliI"),
  T6SS = c("hcp","vgrG","tssA","tssB","tssC","tssD","clpV"),
  LPS_Core = c("lpxA","lpxB","lpxC","lpxD","kdsA","kdsB")
)

all_genes <- rownames(glycan_seu)
genesets_clean <- lapply(genesets, function(x) intersect(x, all_genes))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]

glycan_seu <- AddModuleScore(
  glycan_seu,
  features = genesets_clean,
  name = names(genesets_clean),
  nbin = 24,
  seed = 123
)

write.csv(glycan_seu@meta.data, file.path(out_dir, "module_scores.csv"))

# ============================================================
# 6. DotPlots
# ============================================================

p_protective <- DotPlot(
  glycan_seu,
  features = c("Plant_PUL1","Vitamin_B12_MK3","Antioxidant4"),
  group.by = "seurat_clusters",
  scale = TRUE
) + theme_bw()

ggsave(file.path(out_dir, "Dotplot_Protective.png"), p_protective)

p_risk <- DotPlot(
  glycan_seu,
  features = c("LPS_Core10","Ethanol_Acetaldehyde18","Sugar_Overload19"),
  group.by = "seurat_clusters",
  scale = TRUE
) + theme_bw()

ggsave(file.path(out_dir, "Dotplot_Risk.png"), p_risk)

p_scfa <- DotPlot(
  glycan_seu,
  features = c("Acetate22","Propionate21","Butyrate_core2"),
  group.by = "seurat_clusters",
  scale = TRUE
) + theme_bw()

ggsave(file.path(out_dir, "Dotplot_SCFA.png"), p_scfa)

# ============================================================
# 7. Composition
# ============================================================

glycan_seu@active.ident <- as.factor(glycan_seu$group)

prop_df_group <- as.data.frame(prop.table(table(Idents(glycan_seu), glycan_seu$seurat_clusters)))
colnames(prop_df_group) <- c("Group","Cluster","Freq")

p_bar_group <- ggplot(prop_df_group, aes(x = Cluster, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  theme_classic()

ggsave(file.path(out_dir, "composition.png"), p_bar_group)

write.csv(prop_df_group, file.path(out_dir, "composition_source.csv"))

# ============================================================
# 8. Marker genes + heatmap
# ============================================================

Idents(glycan_seu) <- glycan_seu$seurat_clusters

markers <- FindAllMarkers(glycan_seu, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write.csv(markers, file.path(out_dir, "markers.csv"))

top_genes <- markers %>%
  filter(!grepl("^MGYG", gene)) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  pull(gene)

expr <- AverageExpression(glycan_seu, features = top_genes, group.by = "seurat_clusters")$SCT

expr_scaled <- scale(t(expr))

ht <- Heatmap(expr_scaled, name = "expression", col = mycol2)

png(file.path(out_dir, "heatmap.png"), width = 8000, height = 6000, res = 500)
draw(ht)
grid.text("Top marker genes", x = 0.5, y = 1.1)
dev.off()

pdf(file.path(out_dir, "heatmap.pdf"))
draw(ht)
grid.text("Top marker genes", x = 0.5, y = 1.1)
dev.off()
