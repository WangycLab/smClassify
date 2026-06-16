############################################################
## Figure 5c,d,e correlated-species subpopulation analysis
##
############################################################

# ============================================================
# 1. Load packages
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ============================================================
# 2. Input and output paths
# ============================================================
base_dir   <- "Figure5"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

annotated_seurat_rds <- "Figure3/Output/Figure3_mouse_bacterial_functional_cluster_annotated.rds"
cor_species_type_csv <- file.path(input_dir, "cor_species_type.csv")

dpi_out <- 300

# ============================================================
# 3. Load data
# ============================================================
seu <- readRDS(annotated_seurat_rds)

# ============================================================
# 4. Define colours
# ============================================================
species_cols <- c(
  "Duncaniella muris" = "#91D1C2",
  "Muribaculum gordoncarteri" = "#3C5488",
  "Bacteroides muris" = "#FF6347"
)

region_cols <- c(
  "Cecum" = "#1f77b4",
  "Colon" = "#2ca02c",
  "Rectum" = "#ff7f0e"
)

cluster_cols <- c(
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85", "#E64B35", "#4DBBD5",
  "#00A087", "#3C5488"
)

subtype_cols <- c(
  "#6A5ACD", "#2CA02C", "#F5A0A1", "#BDB76B", "#DC143C",
  "#6A5ACD", "#C2B5D8", "#1F77B4", "#BCBD22", "#17BECF",
  "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000",
  "#FFD700", "#DC143C", "#228B22", "#483D8B", "#BDB76B",
  "#20B2AA", "#FF1493", "#FF4500", "#32CD32", "#3E8E41"
)

# ============================================================
# 5. Subset correlated species and recluster
# ============================================================
cor_species <- c(
  "Duncaniella muris",
  "Muribaculum gordoncarteri",
  "Bacteroides muris"
)

seu@active.ident <- as.factor(seu$species)
cor_seu <- subset(seu, idents = cor_species)

write.csv(
  as.data.frame(table(cor_seu$group)),
  file.path(output_dir, "Figure5cde_correlated_species_group_cell_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(cor_seu$species)),
  file.path(output_dir, "Figure5c_correlated_species_cell_counts.csv"),
  row.names = FALSE
)

cor_seu <- SCTransform(cor_seu, verbose = FALSE, vars.to.regress = "nCount_RNA")
cor_seu <- RunPCA(cor_seu, verbose = FALSE)
cor_seu <- FindNeighbors(cor_seu, dims = 1:10)
cor_seu <- FindClusters(cor_seu, resolution = 0.1)
cor_seu <- RunUMAP(cor_seu, dims = 1:10)

saveRDS(
  cor_seu,
  file.path(output_dir, "Figure5cde_correlated_species_reclustered.rds")
)

# ============================================================
# 6. Figure 5c: UMAP colored by species
# ============================================================
p_5c <- DimPlot(
  cor_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "species"
) +
  scale_color_manual(values = species_cols) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  ) +
  labs(title = "Highly correlated species")

ggsave(file.path(output_dir, "Figure5c_correlated_species_UMAP_by_species.png"),
       p_5c, width = 8.5, height = 4.5, dpi = dpi_out)
ggsave(file.path(output_dir, "Figure5c_correlated_species_UMAP_by_species.pdf"),
       p_5c, width = 8.5, height = 4.5)

# ============================================================
# 7. Supplementary UMAPs before subtype annotation
# ============================================================
p_region <- DimPlot(
  cor_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "region"
) +
  scale_color_manual(values = region_cols) +
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

ggsave(file.path(output_dir, "Supplementary_correlated_species_UMAP_by_region.png"),
       p_region, width = 6, height = 4, dpi = dpi_out)

p_phenotype <- DimPlot(
  cor_seu,
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

ggsave(file.path(output_dir, "Supplementary_correlated_species_UMAP_by_phenotype.png"),
       p_phenotype, width = 6, height = 4, dpi = dpi_out)

p_cluster <- DimPlot(
  cor_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "seurat_clusters"
) +
  scale_color_manual(values = cluster_cols) +
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

ggsave(file.path(output_dir, "Supplementary_correlated_species_UMAP_by_cluster.png"),
       p_cluster, width = 6, height = 4, dpi = dpi_out)

# ============================================================
# 8. Supplementary cluster species-composition bar plot
# ============================================================
cluster_species_df <- data.frame(
  cluster = cor_seu$seurat_clusters,
  species = cor_seu$species
)

comp_tab <- cluster_species_df %>%
  group_by(cluster, species) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(freq = count / sum(count))

write.csv(
  comp_tab,
  file.path(output_dir, "Supplementary_correlated_species_cluster_species_composition.csv"),
  row.names = FALSE
)

p_comp <- ggplot(comp_tab, aes(x = cluster, y = freq, fill = species)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = species_cols) +
  labs(x = "Seurat Cluster", y = "Proportion", fill = "Species") +
  theme_minimal(base_size = 14) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

ggsave(file.path(output_dir, "Supplementary_correlated_species_cluster_species_composition.png"),
       p_comp, width = 5, height = 4, dpi = dpi_out)
ggsave(file.path(output_dir, "Supplementary_correlated_species_cluster_species_composition.pdf"),
       p_comp, width = 5, height = 4)

# ============================================================
# 9. Marker genes and subtype annotation
# ============================================================
Idents(cor_seu) <- cor_seu$seurat_clusters

message("Finding markers ...")
markers <- FindAllMarkers(
  object = cor_seu,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.15
)

write.csv(
  markers,
  file.path(output_dir, "Supplementary_correlated_species_markers.csv"),
  row.names = FALSE
)

cor_type <- read.csv(cor_species_type_csv)

new.cluster.ids <- as.character(cor_type[, 2])
cor_seu@active.ident <- as.factor(cor_seu$seurat_clusters)
names(new.cluster.ids) <- levels(cor_seu)

cor_seu <- RenameIdents(cor_seu, new.cluster.ids)
cor_seu$cellsubtype <- cor_seu@active.ident

write.csv(
  as.data.frame(table(cor_seu$cellsubtype)),
  file.path(output_dir, "Figure5d_correlated_species_subtype_cell_counts.csv"),
  row.names = FALSE
)

# ============================================================
# 10. Figure 5d: UMAP colored by subtype
# ============================================================
p_5d <- DimPlot(
  cor_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "cellsubtype"
) +
  scale_color_manual(values = subtype_cols) +
  ggtitle("Cell Subtype") +
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

ggsave(file.path(output_dir, "Figure5d_correlated_species_UMAP_by_subtype.png"),
       p_5d, width = 9, height = 4, dpi = dpi_out)
ggsave(file.path(output_dir, "Figure5d_correlated_species_UMAP_by_subtype.pdf"),
       p_5d, width = 9, height = 4, device = cairo_pdf)

saveRDS(
  cor_seu,
  file.path(output_dir, "Figure5cde_correlated_species_subtype_annotated.rds")
)

# ============================================================
# 11. Supplementary top marker-gene dot plot
# ============================================================
topN <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()

write.csv(
  topN,
  file.path(output_dir, "Supplementary_correlated_species_markers_top50.csv"),
  row.names = FALSE
)

top5_nonMGYG <- markers %>%
  filter(!grepl("^MGYG", gene)) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  ungroup()

write.csv(
  top5_nonMGYG,
  file.path(output_dir, "Supplementary_correlated_species_top5_nonMGYG_markers.csv"),
  row.names = FALSE
)

marker_genes <- unique(top5_nonMGYG$gene)

cor_seu@active.ident <- as.factor(cor_seu$cellsubtype)

p_marker <- DotPlot(cor_seu, features = marker_genes) +
  scale_color_gradient(low = "#2ECC71", high = "#8E44AD") +
  theme_minimal(base_size = 13) +
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

ggsave(file.path(output_dir, "Supplementary_correlated_species_top_marker_gene_dotplot.pdf"),
       p_marker, width = 12, height = 4)
ggsave(file.path(output_dir, "Supplementary_correlated_species_top_marker_gene_dotplot.png"),
       p_marker, width = 12, height = 4, dpi = dpi_out)

# ============================================================
# 12. Figure 5e: functional module activity dot plot
# ============================================================
genesets <- list(
  CellEnvelope_LPSBiogenesis = c(
    "lpxA", "lpxC", "lpxD", "kdsA", "waaA", "rfbB", "gmd", "fcl", "lptD", "ompA"
  ),
  Mucosal_GlycanUtilization = c(
    "susC", "susD", "pulA", "amyE", "xynD", "cbh", "axe1-6A", "algA", "tagO"
  ),
  Carbohydrate_Uptake = c(
    "xylE", "mglB", "sglT", "glpT", "msmX", "glpK", "glpA", "glpB", "glpC", "lacZ"
  ),
  Fermentation_Pathways = c(
    "pgi", "pfkA", "gapA", "tpiA", "pgk", "pyk", "pflB", "pta", "ackA", "adhE"
  ),
  Motility_Chemotaxis = c(
    "fliC", "fljB", "flgG", "hag", "fliA", "cheA", "swrC"
  ),
  Redox_StressDefense = c(
    "sodB", "katA", "ahpC", "tpx", "trxA", "trxB", "fprA", "rbr", "norV", "hcp"
  ),
  AcidResistance_pH = c(
    "gadB", "gadC", "ureA", "ureB", "ureC", "ureD", "ureE", "ureF", "ureG"
  ),
  Efflux_IntrinsicResistance = c(
    "mexB", "acrF", "mdtA", "macB", "oprM", "ttgD", "bepE", "bepF", "bepG"
  )
)

all_genes <- rownames(cor_seu)
genesets_clean <- lapply(genesets, function(g) intersect(g, all_genes))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]

write.csv(
  data.frame(
    Module = rep(names(genesets_clean), lengths(genesets_clean)),
    Gene = unlist(genesets_clean, use.names = FALSE)
  ),
  file.path(output_dir, "Figure5e_functional_modules_used.csv"),
  row.names = FALSE
)

cor_seu <- AddModuleScore(
  cor_seu,
  features = genesets_clean,
  name = names(genesets_clean),
  nbin = 6
)

cor_seu@active.ident <- as.factor(cor_seu$species)

modules <- c(
  "Fermentation_Pathways4",
  "Carbohydrate_Uptake3",
  "AcidResistance_pH7",
  "Efflux_IntrinsicResistance8",
  "Mucosal_GlycanUtilization2",
  "Redox_StressDefense6",
  "Motility_Chemotaxis5",
  "CellEnvelope_LPSBiogenesis1"
)

write.csv(
  cor_seu@meta.data[
    ,
    c("orig.ident", "species", "seurat_clusters", "cellsubtype", modules),
    drop = FALSE
  ],
  file.path(output_dir, "Figure5e_functional_module_scores_source_data.csv"),
  row.names = FALSE
)

p_5e <- DotPlot(
  cor_seu,
  features = modules,
  group.by = "cellsubtype",
  scale = TRUE
) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red3",
    midpoint = 0
  )

ggsave(file.path(output_dir, "Figure5e_functional_module_activity_dotplot.png"),
       p_5e, width = 9, height = 4, dpi = dpi_out)
ggsave(file.path(output_dir, "Figure5e_functional_module_activity_dotplot.pdf"),
       p_5e, width = 9, height = 4)

# ============================================================
# 13. Supplementary individual functional-gene dot plot
# ============================================================
genesets <- list(
  Fermentation_Pathways = c(
    "pgi", "pfkA", "gapA", "tpiA", "pgk", "pyk", "pflB", "pta", "ackA", "adhE"
  ),
  Carbohydrate_Uptake = c(
    "xylE", "mglB", "sglT", "glpT", "msmX", "glpK", "glpA", "glpB", "glpC", "lacZ"
  ),
  AcidResistance_pH = c(
    "gadB", "gadC", "ureA", "ureB", "ureC", "ureD", "ureE", "ureF", "ureG"
  ),
  Efflux_IntrinsicResistance = c(
    "mexB", "acrF", "mdtA", "macB", "oprM", "ttgD", "bepE", "bepF", "bepG"
  ),
  Mucosal_GlycanUtilization = c(
    "susC", "susD", "pulA", "amyE", "xynD", "cbh", "axe1-6A", "algA", "tagO"
  ),
  Redox_StressDefense = c(
    "sodB", "katA", "ahpC", "tpx", "trxA", "trxB", "fprA", "rbr", "norV", "hcp"
  ),
  Motility_Chemotaxis = c(
    "fliC", "fljB", "flgG", "hag", "fliA", "cheA", "swrC"
  ),
  CellEnvelope_LPSBiogenesis = c(
    "lpxA", "lpxC", "lpxD", "kdsA", "waaA", "rfbB", "gmd", "fcl", "lptD", "ompA"
  )
)

all_functional_genes <- unique(unlist(genesets))

p_gene <- DotPlot(
  cor_seu,
  features = all_functional_genes,
  group.by = "cellsubtype",
  scale = TRUE
) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red3",
    midpoint = 0
  )

ggsave(file.path(output_dir, "Supplementary_functional_gene_activity_dotplot.png"),
       p_gene, width = 20, height = 3.5, dpi = dpi_out)

saveRDS(
  cor_seu,
  file.path(output_dir, "Figure5cde_correlated_species_final_object.rds")
)