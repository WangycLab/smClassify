# This script performs a genus-level deep dive on Duncaniella:
#   1. Load pre-processed Seurat object and metadata
#   2. Subset Duncaniella cells and compute pathway module scores
#   3. Run dimensionality reduction, clustering, and UMAP visualization
#   4. Identify marker genes and annotate them with GTF-derived metadata
#   5. Re-label clusters into functional subtypes
#   6. Plot subtype/species proportions across experimental groups
# ================================================================

library(networkD3)
library(dplyr)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(gplots)
library(rtracklayer)
library(scales)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# ---------------------------------------------------------------
# 0. Project directory & data loading
# ---------------------------------------------------------------
# Set the project directory to the root of this analysis.
# Adapt this line when using the script on a new machine or path.
project_dir <- "path/to/project_root"

# Load Seurat object (must contain `seu`)
load(file.path(project_dir, "seu_symbol_SCT_clustered_annotated.RData"))

# Set working directory for all outputs (figures, tables, RDS, etc.)
setwd(file.path(project_dir, "Figure5"))

# ---------------------------------------------------------------
# 1. Color palettes for UMAPs and barplots
# ---------------------------------------------------------------
nature_palette <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85"
)

cell_type_cols <- c(
  "#1F77B4", "#2CA02C", "#FF7F0E", "#6A5ACD", "#8C564B", "#E377C2", "#7F7F7F",
  "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32",
  "#8B0000", "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B",
  "#20B2AA", "#FF1493", "#FF4500", "#32CD32", "#3E8E41", "#20B2AA"
)

# ---------------------------------------------------------------
# 2. Subset Seurat object to genus Duncaniella
# ---------------------------------------------------------------
seu@active.ident <- as.factor(seu$Tax_Genus)
table(seu$Tax_Genus)

Dun_seu  <- subset(seu, idents = "Duncaniella")
table(Dun_seu$species)

# ---------------------------------------------------------------
# 3. Pathway module scoring (AddModuleScore)
#    - Requires `genesets` object defined in the environment
# ---------------------------------------------------------------

# Score modules using pre-defined gene sets (genesets_clean)
Dun_seu <- AddModuleScore(
  Dun_seu,
  features = genesets_clean,
  name = names(genesets_clean),
  nbin = 6
)

# Inspect metadata and set identity to species for plotting
colnames(Dun_seu@meta.data)
Dun_seu@active.ident <- as.factor(Dun_seu$species)

# Filter gene sets to genes present in this object
all_genes <- rownames(Dun_seu)

# Keep only genes that overlap with Dun_seu and remove empty sets
genesets_clean <- lapply(genesets, function(g) intersect(g, all_genes))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]
genesets_clean

# ---------------------------------------------------------------
# 4. DotPlot of key metabolic and stress pathways by species
# ---------------------------------------------------------------
p <- DotPlot(Dun_seu, features = c(
  "Acetate2",
  "Propionate-MCM19",
  "Propionate-MCC20",
  "Propionate-Acrylate21",
  "Butyrate1",
  "LCFA-Beta-Oxidation13",
  "LCFA-Other-Aux16",
  "LCFA-Regulation15",
  "LCFA-Auxiliary12",
  "LCFA-Transport-Efflux14",
  "LCFA-Activation11",
  "PUL-Core4",
  "Nitrogen-Metabolism10",
  "Mucin5",
  "Bile-Total8",
  "Bile-Hydrolase6",
  "Bile-Acid-Secondary7",
  "Stress17",
  "Efflux-Resistance18",
  "Stickland9"
),
group.by = "species",
scale = TRUE) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3",
    midpoint = 0
  )

p
ggsave("Duncaniella_Key Metabolic and Stress Pathways_top20 species.png",
       p, width = 10, height = 5, dpi = 600)

# ---------------------------------------------------------------
# 5. SCTransform, PCA, neighbors, clustering & UMAP
# ---------------------------------------------------------------
Dun_seu <- SCTransform(Dun_seu, verbose = FALSE, vars.to.regress = "nCount_RNA")
Dun_seu <- RunPCA(Dun_seu, verbose = FALSE)
Dun_seu <- FindNeighbors(Dun_seu, dims = 1:6)
Dun_seu <- FindClusters(Dun_seu, resolution = 0.1)
Dun_seu <- RunUMAP(Dun_seu, dims = 1:6)

# ---------------------------------------------------------------
# 6. UMAP plots colored by cluster, region, phenotype, species
# ---------------------------------------------------------------

# 6a. By Seurat clusters
p <- DimPlot(Dun_seu, reduction = "umap", label = FALSE, pt.size = 0.1,
             group.by = "seurat_clusters") +
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
p
ggsave("Dun_DimPlot_seurat_clusters.png", p, width = 4, height = 3, dpi = 300)

# Save & reload object (for reproducible checkpoint)
saveRDS(Dun_seu, "Dun_seu.Rds")
Dun_seu <- readRDS("Dun_seu.Rds")

# 6b. By region (Cecum/Colon/Rectum)
p <- DimPlot(Dun_seu, reduction = "umap", label = FALSE, pt.size = 0.1,
             group.by = "region") +
  theme_minimal() +
  scale_color_manual(values = c(
    "Cecum" = "#1f77b4",
    "Colon" = "#2ca02c",
    "Rectum" = "#ff7f0e"
  )) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )
p
ggsave("Dun_DimPlot_seurat_region.png", p, width = 4, height = 3, dpi = 300)

# 6c. By phenotype (e.g., WT vs T2DM)
p <- DimPlot(Dun_seu, reduction = "umap", label = FALSE, pt.size = 0.1,
             group.by = "phenotype") +
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
p
ggsave("Dun_DimPlot_seurat_phenotype.png", p, width = 4, height = 3, dpi = 300)

# 6d. By species
p <- DimPlot(Dun_seu, reduction = "umap", label = FALSE, pt.size = 0.1,
             group.by = "species") +
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
p
ggsave("Dun_DimPlot_species.png", p, width = 6, height = 3, dpi = 300)

# ---------------------------------------------------------------
# 7. Differential gene analysis (FindAllMarkers)
# ---------------------------------------------------------------
Idents(Dun_seu) <- Dun_seu$seurat_clusters

## 7a. Differential gene screening
message("Finding markers")
markers <- FindAllMarkers(
  object          = Dun_seu,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.15
)
markers

# Save full marker table
write.csv(markers, "Dun_seu_markers_min.pct0.05.csv")
table(markers$cluster)

## 7b. Top 50 marker genes per cluster
topN <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()
table(topN$cluster)

write.csv(topN, "Dun_seu_markers_top50.csv")

# ---------------------------------------------------------------
# 8. Annotation of marker genes using wide-format GTF
# ---------------------------------------------------------------

# Annotation table (symbol/product columns) under project_dir
gtf_wide <- fread(file.path(project_dir, "gtf_wide.csv"))

# Harmonize gene symbol naming and join with annotation
topN$symbol <- gsub("-", "_", topN$gene)
topN_gtf <- topN %>% left_join(gtf_wide, by = "symbol")

topN_gtf_subset <- unique(topN_gtf[, c(
  "p_val", "avg_log2FC", "pct.1", "pct.2",
  "p_val_adj", "cluster", "symbol", "product"
)])
write.csv(topN_gtf_subset, "Dun_seu_markers_top50_gtf.csv", row.names = FALSE)

# ---------------------------------------------------------------
# 9. Select top marker genes for visualization (DotPlot)
# ---------------------------------------------------------------

# 9a. Top 5 per cluster (annotated) for inspection
top10 <- topN_gtf_subset |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 5) |>
  ungroup()

top10 <- top10[!duplicated(top10$symbol), ]
top10$symbol <- gsub("_", "-", top10$symbol)
top10

genes <- unique(top10$symbol)
genes

# 9b. Top 5 non-MGYG genes per cluster
top5_nonMGYG <- topN %>%
  dplyr::filter(!grepl("^MGYG", gene)) %>%     # Remove genes starting with MGYG
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%            # Select top 5 per cluster
  ungroup()
top5_nonMGYG

genes <- unique(top5_nonMGYG$gene)
genes

# 9c. DotPlot of top marker genes
p <- DotPlot(Dun_seu, features = genes) +
  scale_color_gradient(low = "#2ECC71", high = "#8E44AD") +
  theme_minimal(base_size = 13) +
  scale_size_continuous(range = c(1, 8)) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title   = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.grid   = element_blank()
  ) +
  labs(x = "", y = "Functional Subtypes", title = "DotPlot of Top Marker Genes")

print(p)

ggsave("Dun_seu_dotplot_topgenes.pdf", p, width = 12, height = 3.5)
ggsave("Dun_seu_dotplot_topgenes.png", p, width = 12, height = 3.5)

# ---------------------------------------------------------------
# 10. Manual annotation of clusters into functional subtypes
# ---------------------------------------------------------------

# Mapping table: columns = cluster, cellsubtype (from 'Dun_type.csv')
a <- read.csv("Dun_type.csv")
a
levels(Dun_seu)

new.cluster.ids <- as.character(a[, 2])
new.cluster.ids
levels(Dun_seu)

Dun_seu@active.ident <- as.factor(Dun_seu$seurat_clusters)
names(new.cluster.ids) <- levels(Dun_seu)

Dun_seu <- RenameIdents(Dun_seu, new.cluster.ids)
Dun_seu$cellsubtype <- Dun_seu@active.ident

cell_type_cols2 <- c(
  "#1F77B4", "#8C564B", "#E377C2", "#7F7F7F",
  "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32",
  "#8B0000", "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B",
  "#FF1493", "#FF4500", "#32CD32", "#3E8E41", "#20B2AA"
)

p <- DimPlot(Dun_seu, reduction = "umap", label = FALSE, pt.size = 0.1,
             group.by = "cellsubtype") +
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

p
ggsave("Dun_DimPlot_celltype.pdf", p, width = 8, height = 4)

# ---------------------------------------------------------------
# 11. Stacked barplots of Duncaniella proportions across groups
# ---------------------------------------------------------------

# Color palettes
pal_pheno <- c(
  "WT"   = "#4DBBD5",
  "T2DM" = "#E64B35"
)

pal_region <- c(
  "Cecum"  = "#1f77b4",
  "Colon"  = "#2ca02c",
  "Rectum" = "#ff7f0e"
)

# Color palette for 6 group combinations
pal_group <- c(
  "Cecum-WT"   = "#6baed6",
  "Cecum-DB"   = "#08519c",
  "Colon-WT"   = "#74c476",
  "Colon-DB"   = "#006d2c",
  "Rectum-WT"  = "#fdae6b",
  "Rectum-DB"  = "#e6550d"
)

# Meta data for all cells and Duncaniella cells
meta_all <- seu@meta.data
meta_dun <- Dun_seu@meta.data

# ---------------------------------------------------------------
# 11a. Proportions per group, faceted by Duncaniella subtypes
# ---------------------------------------------------------------

# 1. Total cell count per group (denominator)
group_total <- meta_all %>%
  group_by(group) %>%
  summarise(total_cells = n(), .groups = "drop")

# 2. Cell counts of each Duncaniella subtype per group (numerator)
dun_counts <- meta_dun %>%
  group_by(cellsubtype, group) %>%
  summarise(dun_cells = n(), .groups = "drop")

# 3. Merge and compute proportions
df_ratio <- dun_counts %>%
  left_join(group_total, by = "group") %>%
  mutate(freq = dun_cells / total_cells)

# 4. Visualization (faceted by subtype)
p <- ggplot(df_ratio, aes(x = group, y = freq, fill = group)) +
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

p
ggsave("Duncaniella subproportions across groups.pdf", p,
       width = 18, height = 6)

# ---------------------------------------------------------------
# 11b. Proportions per group, faceted by Duncaniella species
# ---------------------------------------------------------------

# 1. Total cell count per group (denominator)
group_total <- meta_all %>%
  group_by(group) %>%
  summarise(total_cells = n(), .groups = "drop")

# 2. Cell counts of each Duncaniella species per group (numerator)
dun_counts <- meta_dun %>%
  group_by(species, group) %>%
  summarise(dun_cells = n(), .groups = "drop")

# 3. Merge and compute proportions
df_ratio <- dun_counts %>%
  left_join(group_total, by = "group") %>%
  mutate(freq = dun_cells / total_cells)

# 4. Visualization (faceted by species)
p <- ggplot(df_ratio, aes(x = group, y = freq, fill = group)) +
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

p
ggsave("Duncaniella subproportions across groups.png", p,
       width = 16, height = 4)
