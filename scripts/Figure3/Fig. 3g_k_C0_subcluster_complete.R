# ============================================================
# Figure 3g-k. Subclustering and functional module analysis of cluster 0
# ===========================================================
# ============================================================

suppressPackageStartupMessages({
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
})

# ============================================================
# 1. Input and output paths
# ============================================================

base_dir   <- "Figure3"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")
out_dir    <- file.path(output_dir, "Figure3g_k_C0_subcluster")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

rds_file <- file.path(input_dir, "Figure3_mouse_bacterial_functional_cluster_annotated.rds")
seu <- readRDS(rds_file)

# ============================================================
# 2. Shared plotting objects from previous Figure 3 scripts
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
# 3. Subset global cluster 0 and rerun subclustering
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
ggsave("C0_glycan_ElbowPlot.png", p_elbow, width = 6, height = 5, dpi = 600)
ggsave("C0_glycan_ElbowPlot.pdf", p_elbow, width = 6, height = 5, device = "pdf")

glycan_seu <- FindNeighbors(glycan_seu, dims = 1:4)
glycan_seu <- FindClusters(glycan_seu, resolution = 0.15)
glycan_seu <- RunUMAP(glycan_seu, dims = 1:4)

saveRDS(glycan_seu, "Figure3g_k_C0_subcluster_seurat.rds")

# ============================================================
# 4. UMAP visualization
# ============================================================

p1 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  group.by = "seurat_clusters"
) +
  scale_color_manual(values = nature_palette2) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave("C0_glycan_UMAP_Clusters.png", p1, width = 8, height = 7, dpi = 1200)
ggsave("C0_glycan_UMAP_Clusters.pdf", p1, width = 8, height = 7, device = "pdf")

p2 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 1,
  group.by = "group"
) +
  scale_color_manual(values = cell_type_cols) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave("C0_glycan_UMAP_Group.png", p2, width = 8, height = 8, dpi = 1200)
ggsave("C0_glycan_UMAP_Group.pdf", p2, width = 8, height = 8, device = "pdf")

p3 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  group.by = "phenotype"
) +
  scale_color_manual(values = c("#E64B35", "#4DBBD5")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave("C0_glycan_UMAP_Phenotype.png", p3, width = 8, height = 6, dpi = 300)
ggsave("C0_glycan_UMAP_Phenotype.pdf", p3, width = 8, height = 6, device = "pdf")

p4 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  group.by = "loc"
) +
  scale_color_manual(values = c("#1F77B4", "#2CA02C", "#FF7F0E")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave("C0_glycan_UMAP_Location.png", p4, width = 8, height = 6, dpi = 300)
ggsave("C0_glycan_UMAP_Location.pdf", p4, width = 8, height = 6, device = "pdf")

# ============================================================
# 5. Functional gene sets and module scores
# ============================================================

genesets <- list(
  Plant_PUL = c("susC", "susD", "susE", "susF",
                "GH13", "GH32", "GH28", "GH10", "GH11", "CE1", "CE4",
                "amyA", "xynA", "xynB", "araF", "araG", "lacZ", "malE", "malF"),
  
  Butyrate_core = c("thl", "hbd", "crt", "bcd", "etfA", "etfB", "croR", "bdhA", "paaH", "paaF", "paaZ"),
  Butyrate_alt  = c("ptb", "buk"),
  
  Vitamin_B12_MK = c("cobA", "cobB", "cobC", "cobD", "cobQ", "cobS", "cobT", "cobU",
                     "cbiA", "cbiC", "cbiD", "cbiK", "menA", "menB", "menC", "ubiE"),
  
  Antioxidant = c("katA", "katE", "sodA", "sodB", "sodC", "ahpC", "ahpF", "trxA", "trxB",
                  "dps", "perR", "oxyR", "msrA", "msrB", "tsaA", "gor", "grxA", "tpx", "bcp"),
  
  Mucin_GH = c("nanA", "nanE", "nanK", "nanT", "nagA", "nagB", "nagC", "nagZ",
               "fucA", "fucI", "fucK", "fucU", "fucP", "GH33", "GH29", "GH95", "GH20", "sulfatase"),
  Bile_acid_secondary = c("baiA", "baiCD", "baiF"),
  
  Bile_Hydrolase = c("cbh", "bsh", "bshA", "bshB", "bshC"),
  
  Efflux_Resistance = c("tolC", "acrA", "acrB", "acrD", "ompR", "ompF", "ompC",
                        "emrA", "emrB", "norM", "mdtK", "mdtF", "cusA", "cusB", "mexB"),
  
  Bile_total = c("cbh", "bsh", "bshA", "bshB", "bshC",
                 "tolC", "acrA", "acrB", "acrD", "ompR", "ompF", "ompC",
                 "emrA", "emrB", "norM", "mdtK", "mdtF", "cusA", "cusB", "mexB"),
  
  LPS_Core = c("lpxA", "lpxB", "lpxC", "lpxD", "lpxK", "lpxL", "lpxM", "lpxP", "kdsA", "kdsB",
               "waaC", "waaF", "waaL", "msbA", "lptA", "lptB", "lptC"),
  
  LPS_Mod_PEtN = c("eptA", "pagP", "lpxT", "arnT", "pmrA", "pmrB", "phoP", "phoQ"),
  
  OM_Lipid_Transport = c("mlaA", "mlaB", "mlaC", "mlaD", "mlaE", "mlaF", "lptD", "lptE", "vacJ"),
  
  Ethanolamine = c("eutA", "eutB", "eutC", "eutD", "eutE", "eutG", "eutH", "eutJ", "eutK", "eutL",
                   "eutM", "eutN", "eutP", "eutQ", "eutR", "eutS", "eutT", "eutV", "eutW", "eutX"),
  
  Succinate_Block = c("frdA", "frdB", "frdC", "frdD", "sdhA", "sdhB", "sdhC", "sdhD",
                      "mdh", "fumA", "fumB", "fumC", "sucA", "sucB", "sucC", "sucD"),
  
  T6SS = c("hcp", "vgrG", "tssA", "tssB", "tssC", "tssD", "tssE", "tssF", "tssG", "tssH",
           "tssI", "tssJ", "tssK", "tssL", "tssM", "tssN", "clpV", "dotU", "icmF"),
  
  Flagella = c("fliC", "fliD", "fliE", "fliF", "fliG", "fliH", "fliI", "fliJ", "fliK", "fliL",
               "fliM", "fliN", "fliO", "fliP", "fliQ", "fliR", "fliS", "fliT", "flgA", "flgB"),
  
  TMA_TMAO = c("cutC", "cutD", "cntA", "cntB", "yeaW", "yeaX",
               "caiA", "caiB", "caiC", "caiD", "caiE", "caiF", "caiT",
               "betT", "torA", "torZ", "torC", "dmsA", "dmsB", "dmsC"),
  
  Ethanol_Acetaldehyde = c("pdc", "adhE", "adhA", "adhB", "adhP", "adhC", "aldA", "aldB", "aldH",
                           "fucO", "yqhD", "yiaY", "frmA", "frmB", "exaA", "exaB", "exaC", "mdh"),
  
  Sugar_Overload = c("ptsI", "ptsH", "lacC", "lacE", "lacF", "lacG", "galT"),
  
  Lipid_Biosynthesis = c("fabD", "fabF"),
  
  Propionate = c("pccA", "pccB", "mutA", "mutB", "epi", "prpE", "prpC", "prpD", "prpB", "acsA",
                 "mmdA", "lcdA", "mmcA", "prpF", "methylcitrate synthase", "methylcitrate dehydratase"),
  
  Acetate = c("pta", "ackA", "acs", "poxB", "yfiD", "ldhA", "pflB", "pflA",
              "pykF", "pykA", "eno", "pgk", "gapA", "tpiA", "pgm"),
  
  SCFA_total = c("thl", "hbd", "crt", "bcd", "etfA", "etfB", "croR", "bdhA", "paaH", "paaF", "paaZ",
                 "ptb", "buk", "atoA", "atoD", "fadE", "etfQ", "etfR",
                 "pccA", "pccB", "mutA", "mutB", "epi", "prpE", "prpC", "prpD", "prpB", "acsA",
                 "mmdA", "lcdA", "mmcA", "prpF",
                 "pta", "ackA", "acs", "poxB", "yfiD", "ldhA", "pflB", "pflA",
                 "pykF", "pykA", "eno", "pgk", "gapA", "tpiA", "pgm"),
  
  Bile_acid_secondary = c("baiA", "baiCD", "baiF")
)

all_genes <- rownames(glycan_seu)

genesets_clean <- lapply(genesets, function(g) intersect(g, all_genes))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]

write.csv(
  data.frame(
    module = names(genesets_clean),
    n_genes = sapply(genesets_clean, length),
    genes = sapply(genesets_clean, paste, collapse = ";")
  ),
  "Figure3g_k_module_gene_sets_used.csv",
  row.names = FALSE
)

modules_to_remove <- grep(
  "Plant_PUL|Butyrate|Vitamin_B12|Antioxidant|Mucin_GH|Bile|LPS|Ethanolamine|Succinate|T6SS|Flagella|TMA_TMAO|Ethanol_Acetaldehyde|Propionate|Acetate|SCFA|Sugar_Overload|Lipid_Biosynthesis",
  colnames(glycan_seu@meta.data),
  value = TRUE
)

glycan_seu@meta.data <- glycan_seu@meta.data[, !(colnames(glycan_seu@meta.data) %in% modules_to_remove)]

cols_to_remove <- c(
  "OM_Lipid_Transport9",
  "Efflux_Resistance8",
  "OM_Lipid_Transport12",
  "Plant_PUL1",
  "Butyrate_core2",
  "Butyrate_total3",
  "Vitamin_B12_MK4",
  "Antioxidant5",
  "Mucin_GH6",
  "Bile_Hydrolase7",
  "Bile_total9",
  "LPS_Core10",
  "LPS_Mod_PEtN11",
  "LPS_total13",
  "Ethanolamine14",
  "Succinate_Block15",
  "T6SS16",
  "Flagella17",
  "TMA_TMAO18",
  "Ethanol_Acetaldehyde19",
  "Sugar_Overload20",
  "Lipid_Biosynthesis21",
  "Propionate22",
  "Acetate23",
  "SCFA_total24",
  "Efflux_Resistance7"
)

glycan_seu@meta.data <- glycan_seu@meta.data[, !(colnames(glycan_seu@meta.data) %in% cols_to_remove)]

glycan_seu <- AddModuleScore(
  glycan_seu,
  features = genesets_clean,
  name = names(genesets_clean),
  nbin = 24,
  seed = 123
)

write.csv(glycan_seu@meta.data, "Figure3g_k_C0_subcluster_metadata_with_module_scores.csv")

# ============================================================
# 6. Module-score DotPlots
# ============================================================

p_protective <- DotPlot(
  glycan_seu,
  features = c("Plant_PUL1", "Vitamin_B12_MK3", "Antioxidant4"),
  group.by = "seurat_clusters",
  scale = TRUE
) +
  scale_color_gradient(low = "white", high = "red3") +
  coord_flip() +
  scale_size_continuous(range = c(4, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Protective")

ggsave("Dotplot_Protective.pdf", p_protective, width = 7, height = 5.6, device = "pdf")
ggsave("Dotplot_Protective.png", p_protective, width = 7, height = 5.6, dpi = 600)

p_risk <- DotPlot(
  glycan_seu,
  features = c(
    "LPS_Mod_PEtN11",
    "LPS_Core10",
    "Ethanol_Acetaldehyde18",
    "Sugar_Overload19",
    "Lipid_Biosynthesis20",
    "Mucin_GH5"
  ),
  group.by = "seurat_clusters",
  scale = TRUE
) +
  scale_color_gradient(low = "white", high = "#1F77B4") +
  coord_flip() +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Risk")

ggsave("Dotplot_Risk.pdf", p_risk, width = 7.5, height = 5.5, device = "pdf")
ggsave("Dotplot_Risk.png", p_risk, width = 7.5, height = 5.5, dpi = 600)

p_scfa <- DotPlot(
  glycan_seu,
  features = c("Acetate22", "Propionate21", "Butyrate_core2", "SCFA_total23"),
  group.by = "seurat_clusters",
  scale = TRUE
) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_flip() +
  scale_size_continuous(range = c(4, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("SCFA_ContextFlux")

ggsave("Dotplot_SCFA_ContextFlux.pdf", p_scfa, width = 6.5, height = 5.5, device = "pdf")
ggsave("Dotplot_SCFA_ContextFlux.png", p_scfa, width = 6.5, height = 5.5, dpi = 600)

# ============================================================
# 7. Composition bar plots
# ============================================================

glycan_seu@active.ident <- as.factor(glycan_seu$group)
prop_df_group <- as.data.frame(prop.table(table(Idents(glycan_seu), glycan_seu$seurat_clusters)))
colnames(prop_df_group) <- c("Group", "Cluster", "Freq")
prop_df_group$Group <- factor(prop_df_group$Group, levels = unique(prop_df_group$Group))

p_bar_group <- ggplot(prop_df_group, aes(x = Cluster, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c("#1F77B4", "#17BECF", "#2CA02C", "#9ACD32", "#FF7F0E", "#F5A0A1")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(x = "Cluster", y = "Proportion", title = "Group Composition by clsuter")

ggsave("C0_group_composition_by_cluster.pdf", p_bar_group, width = 7, height = 5, device = "pdf")
ggsave("C0_group_composition_by_cluster.png", p_bar_group, width = 7, height = 5, dpi = 600)
write.csv(prop_df_group, "C0_group_composition_by_cluster_source_data.csv", row.names = FALSE)

# ============================================================
# 8. Marker genes and ComplexHeatmap
# ============================================================

Idents(glycan_seu) <- glycan_seu$seurat_clusters

markers <- FindAllMarkers(
  object = glycan_seu,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.1
)
write.csv(markers, "glycan_seu_markers.csv")

top5_nonMGYG <- markers %>%
  dplyr::filter(!grepl("^MGYG", gene)) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()

write.csv(top5_nonMGYG, "glycan_seu_top10_nonMGYG_markers.csv", row.names = FALSE)

genes <- unique(top5_nonMGYG$gene)

aver_dt <- AverageExpression(
  glycan_seu,
  assays = "SCT",
  features = genes,
  group.by = "seurat_clusters",
  slot = "data"
)$SCT %>%
  as.data.frame()

write.csv(aver_dt, "Subcluster0_Top10_nonMGYG_marker_average_expression.csv")

expr_mat <- aver_dt

celltypes <- colnames(expr_mat)
cell_anno <- data.frame(cell_anno = celltypes, row.names = celltypes)

celltype_col <- setNames(
  brewer.pal(length(unique(cell_anno$cell_anno)), "Set3"),
  unique(cell_anno$cell_anno)
)

expr_scaled_t <- scale(t(expr_mat))

clusters <- rownames(expr_scaled_t)
row_anno_df <- data.frame(cluster = clusters, row.names = clusters)

cluster_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(length(unique(row_anno_df$cluster))),
  unique(row_anno_df$cluster)
)

row_anno_t <- rowAnnotation(
  cluster = row_anno_df$cluster,
  col = list(cluster = cluster_colors),
  annotation_label = "Cluster",
  gp = gpar(col = "white", lwd = 1.5)
)

ht <- Heatmap(
  expr_scaled_t,
  name = "expression",
  col = mycol2,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_gp = gpar(fontsize = 10, fontface = "italic"),
  rect_gp = gpar(col = "white", lwd = 1.5),
  left_annotation = row_anno_t
)

png("Subcluster0_Top5_marker_genes_heatmap.png", width = 10000, height = 6500, res = 500)
draw(ht)
grid.text(
  "Top10 marker genes across clusters",
  x = 0.5,
  y = 1.5,
  gp = gpar(fontsize = 18, fontface = "bold")
)
dev.off()

pdf("Subcluster0_Top5_marker_genes_heatmap.pdf", width = 8, height = 3.5)
draw(ht)
grid.text(
  "Top10 marker genes across clusters",
  x = 0.5,
  y = 1.2,
  gp = gpar(fontsize = 18, fontface = "bold")
)
dev.off()

