# =======================================================================
# Correlated-Species Single-cell Analysis (B. muris / M. gordoncarteri /
# D. muris)
# -----------------------------------------------------------------------
# This script performs a focused Seurat analysis on three correlated
# species:
#   1. Load pre-processed Seurat object and define top-20 species
#   2. Subset three species of interest and run SCTransform + UMAP + clustering
#   3. Generate UMAP plots by species, region, phenotype, cluster
#   4. Assess species composition across Seurat clusters (stacked barplot)
#   5. Identify marker genes, re-annotate clusters into cell subtypes
#   6. Visualize top marker genes via DotPlot
#   7. Define species-specific functional gene modules and score them
#   8. Plot module-level and gene-level DotPlots across subtypes
# =======================================================================
# 0. Libraries
# -----------------------------------------------------------------------
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

# -----------------------------------------------------------------------
# 1. Project directory, data loading, and working directory
# -----------------------------------------------------------------------
# Set this to the root directory of the project when running the script.
project_dir <- "path/to/project_root"

# Load Seurat object containing `seu`
load(file.path(project_dir, "seu_symbol_SCT_clustered_annotated.RData"))

# Set working directory for all outputs (figures, tables, RDS, etc.)
setwd(file.path(project_dir, "Figure5"))

# -----------------------------------------------------------------------
# 2. Define top 20 species & subset correlated species
# -----------------------------------------------------------------------

# Top 20 species:
# 1) Count and sort
sp_counts <- sort(table(seu$species), decreasing = TRUE)
n_top <- min(20, length(sp_counts))

top_species <- names(sp_counts)[seq_len(n_top)]
seu$top20_species <- ifelse(
  seu$species %in% top_species,
  as.character(seu$species),
  "Others"
)

seu$top20_species <- factor(
  seu$top20_species,
  levels = c(top_species, "Others")
)

n_other <- sum(!(seu$species %in% top_species))
levels(seu$top20_species)[levels(seu$top20_species) == "Others"]

# Quick check
table(seu$top20_species)

# Color palettes
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

# Correlated species panel
cor_species <- c(
  "Duncaniella muris",
  "Muribaculum gordoncarteri",
  "Bacteroides muris"
)

seu@active.ident <- as.factor(seu$species)

cor_seu <- subset(seu, idents = cor_species)

# Basic counts
table(cor_seu$group)
table(seu$group)
table(seu$top20_species)
table(cor_seu$species)

# -----------------------------------------------------------------------
# 3. SCTransform, PCA, neighbors, clustering & UMAP on correlated species
# -----------------------------------------------------------------------
cor_seu <- SCTransform(cor_seu, verbose = FALSE, vars.to.regress = "nCount_RNA")
cor_seu <- RunPCA(cor_seu, verbose = FALSE)
cor_seu <- FindNeighbors(cor_seu, dims = 1:10)
cor_seu <- FindClusters(cor_seu, resolution = 0.1)
cor_seu <- RunUMAP(cor_seu, dims = 1:10)

# -----------------------------------------------------------------------
# 4. UMAP plots colored by species, region, phenotype and clusters
# -----------------------------------------------------------------------

# 4a. UMAP by species
nature_palette <- c(
  "#91D1C2", "#3C5488", "#FF6347"
)

p <- DimPlot(
  cor_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "species"
) +
  scale_color_manual(values = nature_palette) +
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
ggsave("Cor_DimPlot_species.png", p, width = 8.5, height = 4.5, dpi = 300)
ggsave("Cor_DimPlot_species.pdf", p, width = 8.5, height = 4.5, dpi = 300)

# Save and reload for consistency
saveRDS(cor_seu, "cor_seu.Rds")
cor_seu <- readRDS("cor_seu.Rds")

# 4b. UMAP by region (Cecum / Colon / Rectum)
p <- DimPlot(
  cor_seu, reduction = "umap",
  label = FALSE, pt.size = 0.1,
  group.by = "region"
) +
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
ggsave("cor_seu _DimPlot_seurat_region.png", p, width = 6, height = 4, dpi = 300)

# 4c. UMAP by phenotype (e.g., WT vs T2DM)
p <- DimPlot(
  cor_seu, reduction = "umap",
  label = FALSE, pt.size = 0.1,
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
p
ggsave("cor_seu _DimPlot_seurat_phenotype.png", p, width = 6, height = 4, dpi = 300)

# 4d. UMAP by Seurat clusters
nature_palette <- c(
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85", "#E64B35", "#4DBBD5", "#00A087", "#3C5488"
)

p <- DimPlot(
  cor_seu, reduction = "umap",
  label = FALSE, pt.size = 0.1,
  group.by = "seurat_clusters"
) +
  theme_minimal() +
  scale_color_manual(values = nature_palette) +
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
ggsave("cor_seu _DimPlot_seurat_clusters.png", p, width = 6, height = 4, dpi = 300)

# -----------------------------------------------------------------------
# 5. Species composition per cluster (stacked proportional barplot)
# -----------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

# Extract cluster and species information
df <- data.frame(
  cluster = cor_seu$seurat_clusters,
  species = cor_seu$species
)

# Count number of each species in each cluster and compute proportions
comp_tab <- df %>%
  group_by(cluster, species) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(freq = count / sum(count))   # proportion per cluster

# Plot (proportional stacked bar chart)
ggplot(comp_tab, aes(x = cluster, y = freq, fill = species)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = nature_palette) +
  labs(
    x = "Seurat Cluster",
    y = "Proportion",
    fill = "Species"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

# Ensure unique color mapping (extra palette if needed)
cell_type_cols <- c(
  "#1F77B4", "#2CA02C", "#FF7F0E", "#F0E68C", "#6A5ACD", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#F5A0A1",
  "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000", "#FFD700",
  "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B", "#20B2AA",
  "#FF1493", "#FF4500", "#32CD32", "#3E8E41"
)

# -----------------------------------------------------------------------
# 6. Differential gene screening (FindAllMarkers) on correlated species
# -----------------------------------------------------------------------
Idents(cor_seu) <- cor_seu$seurat_clusters

## 3. Differential gene screening and annotation
message("Finding markers")
markers <- FindAllMarkers(
  object          = cor_seu,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.15
)
markers

# -----------------------------------------------------------------------
# 7. Manual annotation of clusters into cell subtypes
# -----------------------------------------------------------------------
a <- read.csv("cor_species_type.csv")
a
levels(cor_seu)

new.cluster.ids <- as.character(a[, 2])
new.cluster.ids
levels(cor_seu)

cor_seu@active.ident <- as.factor(cor_seu$seurat_clusters)
names(new.cluster.ids) <- levels(cor_seu)

cor_seu <- RenameIdents(cor_seu, new.cluster.ids)
cor_seu$cellsubtype <- cor_seu@active.ident

cell_type_cols <- c(
  "#6A5ACD", "#2CA02C", "#F5A0A1", "#BDB76B", "#DC143C", "#6A5ACD", "#C2B5D8",
  "#1F77B4", "#BCBD22", "#17BECF",
  "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000", "#FFD700",
  "#DC143C", "#228B22", "#483D8B", "#BDB76B", "#20B2AA",
  "#FF1493", "#FF4500", "#32CD32", "#3E8E41"
)

p_Subtype <- DimPlot(
  cor_seu,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "cellsubtype"
) +
  scale_color_manual(values = cell_type_cols) +
  ggtitle("Cell Subtype") +
  theme_minimal() +
  scale_color_manual(values = cell_type_cols) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )

p_Subtype
ggsave("cor_seu_UMAP_Subtype.png", p_Subtype, width = 9, height = 4, dpi = dpi_out)
ggsave("cor_seu_UMAP_Subtype.pdf", p_Subtype, width = 9, height = 4, device = cairo_pdf)

# -----------------------------------------------------------------------
# 8. Top-50 markers per cluster and subset of non-MGYG markers
# -----------------------------------------------------------------------
topN <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()
table(topN$cluster)
write.csv(topN, "cor_seu_markers_top50.csv")

top5_nonMGYG <- markers %>%
  filter(!grepl("^MGYG", gene)) %>%    # Remove genes starting with "MGYG"
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%     # Take the top 5 genes for each cluster
  ungroup()

genes <- unique(top5_nonMGYG$gene)
genes

cor_seu@active.ident <- as.factor(cor_seu$cellsubtype)

# Plotting top marker genes across functional subtypes
p <- DotPlot(cor_seu, features = genes) +
  scale_color_gradient(low = "#2ECC71", high = "#8E44AD") +
  theme_minimal(base_size = 13) +
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
  labs(
    x = "",
    y = "Functional Subtypes",
    title = "DotPlot of Top Marker Genes"
  )

print(p)

ggsave("cor_seu_dotplot_topgenes.pdf", p, width = 12, height = 4)
ggsave("cor_seu_dotplot_topgenes.png", p, width = 12, height = 4, dpi = 300)

# -----------------------------------------------------------------------
# 9. Define species-specific functional gene modules and score them
# -----------------------------------------------------------------------
genesets <- list(
  # B. muris
  Bmuris_Sus = c("susC", "susD", "susE", "susF"),
  Bmuris_CAZy = c("GH13", "GH32", "GH28", "GH10", "GH11", "CE1", "CE4"),
  Bmuris_SugarTransport = c("araF", "araG", "malE", "malF", "lacZ"),

  # M. gordoncarteri
  Mgordon_LCFA = c("fadD", "fadN", "atoD", "crt", "hbd", "mutB"),
  Mgordon_Efflux = c("bepF", "bepG", "tolC", "acrB", "mexB"),
  Mgordon_Bile = c("bsh", "baiA", "baiB", "baiCD"),

  # D. muris
  Dmuris_Efflux = c("acrB", "mexB", "emrB", "mdtB", "tolC"),
  Dmuris_Detox = c("katA", "ahpC", "trxA", "rbr", "norV", "hcp"),
  Dmuris_Motility = c("fliC", "flgG", "fljB", "cheA", "swrC"),

  # Common
  Common_OMP = c("ompA", "omp41", "omp40", "inlJ"),
  Common_Detox = c("norV", "hcp", "katA", "ahpC", "trxA", "rbr", "fprA")
)

# Clean gene sets, keep only existing genes
all_genes <- rownames(cor_seu)
genesets_clean <- lapply(genesets, function(g) intersect(g, all_genes))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]
genesets_clean

# AddModuleScore for each module
cor_seu <- AddModuleScore(
  cor_seu,
  features = genesets_clean,
  name = names(genesets_clean),
  nbin = 6
)

# -----------------------------------------------------------------------
# 10. DotPlot for module scores across subtypes and genes across subtypes
# -----------------------------------------------------------------------

# Metadata overview (sanity check)
colnames(cor_seu@meta.data)
cor_seu@active.ident <- as.factor(cor_seu$species)

# 10a. DotPlot of module scores per subtype
p <- DotPlot(
  cor_seu,
  features = c(
    modules <- c(
      modules <- c(
        "Bmuris_Sus1", "Bmuris_SugarTransport2",
        "Mgordon_LCFA3", "Mgordon_Efflux4", "Mgordon_Bile5",
        "Dmuris_Efflux6", "Dmuris_Detox7", "Dmuris_Motility8",
        "Common_OMP9", "Common_Detox10"
      )
    )
  ),
  group.by = "cellsubtype",
  scale = TRUE
) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3",
    midpoint = 0
  )

p
ggsave(
  "Cor_Dot plot of species-specific functional gene module activity.png",
  p, width = 9, height = 3.5, dpi = 300
)

ggsave(
  "Cor_Dot plot of species-specific functional gene module activity.pdf",
  p, width = 9, height = 3.5
)

# 10b. DotPlot of all genes in modules per subtype
all_genes <- unique(unlist(genesets))

# Check merged gene sets
all_genes

p <- DotPlot(
  cor_seu,
  features = all_genes,
  group.by = "cellsubtype",
  scale = TRUE
) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3",
    midpoint = 0
  )

p
ggsave(
  "Cor_Dot plot of species-specific functional gene module activit_genes.png",
  p, width = 13, height = 3.5, dpi = 300
)
