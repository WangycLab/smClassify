# This script performs a genus-level deep dive on Parabacteroides:
#   1. Load pre-processed Seurat object (seu) and set project directory
#   2. Subset Parabacteroides cells and run SCTransform + PCA + UMAP + clustering
#   3. Generate UMAP plots colored by cluster, phenotype, region, species
#   4. Identify cluster markers and visualize top markers via DotPlot
#   5. Manually annotate clusters into Parabacteroides cell subtypes
#   6. Define functional modules related to culture requirements
#   7. Score modules and gene-level signatures per species
#   8. Summarize culture “lever” scores and visualize as heatmap (+/-/●)
#   9. Plot all genes used for culture modules across species

# ============================================================================

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

# ---------------------------------------------------------------------------
# 0. Project configuration & data loading
# ---------------------------------------------------------------------------
# Set this to the root directory of your project when running the script.
project_dir <- "path/to/project_root"

# Load Seurat object (must contain `seu`)
load(file.path(project_dir, "seu_symbol_SCT_clustered_annotated.RData"))

# Set working directory for this figure panel and all downstream outputs
setwd(file.path(project_dir, "Figure5"))

# ---------------------------------------------------------------------------
# 1. Color palettes
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 2. Subset Parabacteroides genus and basic clustering/UMAP
# ---------------------------------------------------------------------------
sample <- seu
sample@active.ident <- as.factor(sample$Tax_Genus)

Parabac <- subset(sample, idents = "Parabacteroides")
table(Parabac$species)

Parabac <- SCTransform(Parabac, verbose = FALSE, vars.to.regress = "nCount_RNA")
Parabac <- RunPCA(Parabac, verbose = FALSE)
Parabac <- FindNeighbors(Parabac, dims = 1:5)
Parabac <- FindClusters(Parabac, resolution = 0.15)
Parabac <- RunUMAP(Parabac, dims = 1:5)

# ---------------------------------------------------------------------------
# 3. UMAP visualizations (clusters / phenotype / region / species)
# ---------------------------------------------------------------------------

# 3a. UMAP colored by Seurat clusters
cell_type_cols <- c(
  "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32",
  "#8B0000", "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B",
  "#FF1493", "#FF4500", "#32CD32", "#3E8E41", "#20B2AA"
)

p <- DimPlot(
  Parabac,
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
p
ggsave("Para_DimPlot_seurat_clusters.png", p, width = 4, height = 3, dpi = 300)

# 3b. UMAP colored by phenotype (first version)
p <- DimPlot(
  Parabac,
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
p
ggsave("Parabac_DimPlot_seurat_phenotype.png", p, width = 4, height = 3, dpi = 300)

# Save & reload object as checkpoint
saveRDS(Parabac, "Parabac.Rds")
Parabac <- readRDS("Parabac.Rds")

# 3c. UMAP colored by region
p <- DimPlot(
  Parabac,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "region"
) +
  theme_minimal() +
  scale_color_manual(values = c(
    "Cecum"  = "#1f77b4",
    "Colon"  = "#2ca02c",
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
ggsave("Parabac_DimPlot_seurat_region.png", p, width = 4.5, height = 3, dpi = 300)

# 3d. UMAP colored by phenotype (second version, same filename)
p <- DimPlot(
  Parabac,
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
p
ggsave("Parabac_DimPlot_seurat_phenotype.png", p, width = 4, height = 3, dpi = 300)

# 3e. UMAP colored by Parabacteroides species
cell_type_cols <- c("#BCBD22", "#228B22", "#7F7F7F")

p <- DimPlot(
  Parabac,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  group.by = "species"
) +
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
p
ggsave("Parabac_DimPlot_species.png", p, width = 6, height = 3, dpi = 300)
ggsave("Parabac_DimPlot_species.pdf", p, width = 6, height = 3)

# ---------------------------------------------------------------------------
# 4. Simple DotPlot for selected genes (mutB, epi) across groups
# ---------------------------------------------------------------------------
DotPlot(
  Parabac,
  features = c("mutB", "epi"),
  group.by = "group"
) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text  = element_text(face = "bold")
  ) +
  labs(title = "BCAA-related genes - Parabacteroides distasonis")

# ---------------------------------------------------------------------------
# 5. Figure layout parameters & common themes (for paper-ready plots)
# ---------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)

## 1. Basic parameters for paper layout
fig_width  <- 10   # inches
fig_height <- 7    # inches
dpi_out    <- 300  # PNG resolution
base_size  <- 10   # base font size

## 2. Color scheme for cell types (generic palette)
cell_type_cols <- c(
  "#1F77B4", "#2CA02C", "#FF7F0E", "#F0E68C", "#6A5ACD", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#F5A0A1",
  "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000", "#FFD700",
  "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B", "#20B2AA",
  "#FF1493", "#FF4500", "#32CD32", "#3E8E41"
)

## 3. General themes for Nature-style UMAPs
nature_theme <- theme_minimal(base_family = "Helvetica", base_size = base_size) +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.title   = element_blank(),
    legend.title = element_text(size = base_size, face = "bold"),
    legend.text  = element_text(size = base_size),
    plot.title   = element_text(size = base_size + 2, face = "bold", hjust = 0.5)
  )

common_theme <- theme_minimal(base_family = "Helvetica", base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10)
  )

# ---------------------------------------------------------------------------
# 6. Differential gene screening (FindAllMarkers) for Parabac clusters
# ---------------------------------------------------------------------------
Idents(Parabac) <- Parabac$seurat_clusters

message("Finding markers")
markers <- FindAllMarkers(
  object          = Parabac,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.15
)
markers

# ---------------------------------------------------------------------------
# 7. Top 50 markers per cluster
# ---------------------------------------------------------------------------
markers$gene

topN <- markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 50) |>
  ungroup()
topN
markers

# ---------------------------------------------------------------------------
# 8. Manual annotation: map clusters to Parabacteroides cell subtypes
# ---------------------------------------------------------------------------
a <- read.csv("Parabac_type.csv")
a
levels(Parabac)

new.cluster.ids <- as.character(a[, 2])
new.cluster.ids
levels(Parabac)

Parabac@active.ident <- as.factor(Parabac$seurat_clusters)
names(new.cluster.ids) <- levels(Parabac)

Parabac <- RenameIdents(Parabac, new.cluster.ids)
Parabac$cellsubtype <- Parabac@active.ident

cell_type_cols <- c(
  "#8C564B", "#2CA02C", "#F5A0A1", "#FF7F0E", "#6A5ACD", "#8B0000", "#C2B5D8",
  "#1F77B4", "#BCBD22", "#17BECF",
  "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000", "#FFD700",
  "#DC143C", "#228B22", "#483D8B", "#BDB76B", "#20B2AA",
  "#FF1493", "#FF4500", "#32CD32", "#3E8E41"
)

cell_type_cols2 <- c("#6A5ACD", "#F5A0A1", "#3CB371", "#FFD700")

p <- DimPlot(
  Parabac,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
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

p
ggsave("Parabac_DimPlot_celltype.pdf", p, w = 7, h = 3)

# ---------------------------------------------------------------------------
# 9. DotPlot of top non-MGYG marker genes per cluster
# ---------------------------------------------------------------------------
Idents(Parabac) <- "seurat_clusters"

top5_nonMGYG <- topN %>%
  dplyr::filter(!grepl("^MGYG", gene)) %>%   # Remove genes starting with MGYG
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%          # Select top 5 per cluster
  ungroup()
top5_nonMGYG
genes <- unique(top5_nonMGYG$gene)
genes

p <- DotPlot(Parabac, features = genes) +
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
  ) +
  scale_size_continuous(range = c(1, 8))

print(p)

ggsave("Parabac_dotplot_topgenes.pdf", p, width = 12, height = 3)
ggsave("Parabac_dotplot_topgenes.png", p, width = 12, height = 3, dpi = 300)

# ---------------------------------------------------------------------------
# 10. Culture requirement analysis across Parabacteroides species
#     (module-level and gene-level scoring)
# ---------------------------------------------------------------------------

# Compare culture requirements among Parabacteroides species
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)

# =========================================================
# Culture levers across Parabacteroides species (final)
# =========================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(reshape2)
  library(stringr)
})

# 10.1 Configuration: functional module gene sets & weights
culture_modules <- list(
  # 1. Inulin/Amylose (plant polysaccharide utilization)
  Inulin_Amylose = c("susC", "susD", "malQ", "glgA"),

  # 2. Mucin/NAG
  Mucin_NAG = c("nanA", "nagA", "nagB", "nagC"),

  # 3. Fucose
  Fucose = c("fucI", "fucO"),

  # 4. Xylose / Mannose / Galactose
  Xylose_Man_Gal = c("xylA", "xylB", "xylE", "manA", "galE", "galK"),

  # 5. Ethanolamine (empty in this object)
  Ethanolamine = character(0),

  # 6. Vitamin B12
  Vitamin_B12 = c("cobB", "cobD", "cbiA", "cbiD"),

  # 7. Fumarate (electron acceptor)
  Fumarate_acceptor = c("frdA", "frdB"),

  # 8. Bile salts tolerance
  Bile_salts = c("cbh", "acrF", "tolC", "ompR"),

  # 9. Branched SCFAs
  Branched_SCFAs = c("ilvE", "leuA", "leuB", "leuC"),

  # 10. Broad vitamin/cofactor metabolism
  Vitamins_broad = c(
    # Biotin
    "bioA", "bioB", "bioC", "bioD", "bioF",
    # Thiamine
    "thiC", "thiD", "thiE",
    # Folate
    "folB", "folC", "folD",
    # NAD
    "nadA", "nadB", "nadC", "nadD", "nadE",
    # Menaquinone
    "menA", "menB", "menD", "menE", "menF"
  )
)

# 10.2 Gene keywords (gene-level evidence)
gene_keywords <- list(
  Mucin_NAG = c("nagA", "nagB", "nanA"),
  Fucose    = c("fucI", "fucP"),
  XMG       = c("xylA", "xylB", "manA", "galE", "galK"),
  Prop_core = c("mutA", "mutB", "pccA", "pccB", "mcmA", "mcmB", "prpB", "prpC", "prpD"),
  B12_core  = c("cobB", "cobD", "cbiA", "cbiD")
)

lever_module_weights <- list(
  "Mucin/NAG"            = c(PUL_sugars = 1),
  "Fucose"               = c(PUL_sugars = 1),
  "Xylose/Man/Gal"       = c(PUL_sugars = 1),
  "Fumarate (e- acceptor)" = c(Propionate = 1),
  "Ethanolamine"         = c(Ethanolamine_B12 = 1),
  "Vitamin B12"          = c(Ethanolamine_B12 = 1),
  "Bile salts (tolerate)"= c(Bile_tolerance = 1)
)

# 10.3 Collapse into six module-score feature sets
module_sets <- list(
  PUL_sugars         = unique(unlist(c(
    culture_modules$Inulin_Amylose,
    culture_modules$Mucin_NAG,
    culture_modules$Fucose,
    culture_modules$Xylose_Man_Gal
  ))),
  Ethanolamine_B12   = unique(unlist(c(
    culture_modules$Ethanolamine,
    culture_modules$Vitamin_B12
  ))),
  Propionate         = culture_modules$Fumarate_acceptor,
  Bile_tolerance     = culture_modules$Bile_salts,
  AA_bios_deg        = culture_modules$Branched_SCFAs,
  Vitamins_cofactors = culture_modules$Vitamins_broad
)

# 10.4 Module-level AddModuleScore per cell
DefaultAssay(Parabac) <- DefaultAssay(Parabac)
for (nm in names(module_sets)) {
  Parabac <- AddModuleScore(
    Parabac,
    features = list(module_sets[[nm]]),
    name = paste0(nm, "_score"),
    vst = FALSE,
    assay = DefaultAssay(Parabac),
    slot = "data",
    nbin = 6
  )
}

# 10.5 Aggregate module scores per species
scores_species <- Parabac@meta.data %>%
  group_by(species) %>%
  summarise(across(ends_with("_score1"), mean, na.rm = TRUE)) %>%
  arrange(species)

scores_species

mat_mod <- as.matrix(scores_species[, -1])
rownames(mat_mod) <- scores_species$species
colnames(mat_mod) <- gsub("_score1$", "", colnames(mat_mod))
z_mod <- scale(mat_mod)

# 10.6 Gene-level scoring helpers
DefaultAssay(Parabac) <- DefaultAssay(Parabac)
all_genes <- rownames(Parabac)

strict_match <- function(vec, all_genes) {
  vec <- unique(vec)
  ag <- all_genes
  ag_low <- tolower(ag)
  res <- unlist(lapply(vec, function(g) {
    g_low <- tolower(g)
    g_low2 <- tolower(gsub("-", "_", g))
    idx <- which(ag_low == g_low | ag_low == g_low2)
    if (length(idx) > 0) ag[idx] else character(0)
  }))
  unique(res)
}

avg_expr_by_species <- function(gene_set) {
  feats <- unique(unlist(lapply(gene_set, strict_match, all_genes = all_genes)))
  feats <- feats[feats %in% rownames(Parabac)]
  if (length(feats) == 0) return(NULL)
  avg <- AverageExpression(
    Parabac,
    features = feats,
    group.by = "species",
    assays = DefaultAssay(Parabac)
  )[[1]]
  rowMeans(avg, na.rm = TRUE)
}

gene_scores <- list()
for (n in names(gene_keywords)) {
  sc <- avg_expr_by_species(gene_keywords[[n]])
  if (!is.null(sc)) {
    names(sc) <- gsub("-", "_", names(sc))
    gene_scores[[n]] <- sc
  }
}

if (length(gene_scores) > 0) {
  gene_df <- do.call(cbind, lapply(gene_scores, function(x) {
    x <- x[rownames(z_mod)]
    return(x)
  }))
  colnames(gene_df) <- names(gene_scores)
  z_gene <- scale(gene_df)
} else {
  z_gene <- NULL
}

# 10.7 Merge module & gene evidence into culture “lever” scores
lever_score <- matrix(0, nrow = nrow(z_mod), ncol = length(lever_module_weights))
rownames(lever_score) <- rownames(z_mod)
colnames(lever_score) <- names(lever_module_weights)

for (lv in names(lever_module_weights)) {
  w <- lever_module_weights[[lv]]
  mod_cols <- intersect(names(w), colnames(z_mod))

  mod_part <- 0
  if (length(mod_cols) > 0) {
    mod_part <- as.matrix(z_mod[, mod_cols, drop = FALSE]) %*%
      matrix(w[mod_cols], ncol = 1)
    mod_part <- as.numeric(mod_part)
  }

  gene_part <- rep(0, nrow(z_mod))
  if (!is.null(z_gene)) {
    if (lv == "Mucin/NAG" && "Mucin_NAG" %in% colnames(z_gene)) gene_part <- z_gene[, "Mucin_NAG"]
    if (lv == "Fucose"    && "Fucose"    %in% colnames(z_gene)) gene_part <- z_gene[, "Fucose"]
    if (lv == "Xylose/Man/Gal" && "XMG"  %in% colnames(z_gene)) gene_part <- z_gene[, "XMG"]
    if (lv == "Fumarate (e- acceptor)" && "Prop_core" %in% colnames(z_gene)) gene_part <- z_gene[, "Prop_core"]
    if (lv %in% c("Ethanolamine", "Vitamin B12") && "B12_core" %in% colnames(z_gene)) gene_part <- z_gene[, "B12_core"]
  }

  gene_part[is.na(gene_part)] <- 0

  if (all(gene_part == 0)) {
    lever_score[, lv] <- mod_part
  } else {
    lever_score[, lv] <- 0.5 * mod_part + 0.5 * gene_part
  }
}

# Additional derived lever: "Avoid bile" = negative of "Bile salts (tolerate)"
lever_score <- cbind(
  lever_score,
  "Avoid bile" = -lever_score[, "Bile salts (tolerate)"]
)

# ---------------------------------------------------------------------------
# 11. Convert lever scores to symbolic matrix (+ / - / ●) and plot
# ---------------------------------------------------------------------------
to_symbol <- function(x, pos_thr = 0.5, neg_thr = -0.5) {
  ifelse(x > pos_thr, "+", ifelse(x < neg_thr, "-", "●"))
}
sym_mat <- apply(lever_score, 2, to_symbol)
sym_mat

sym_long <- melt(sym_mat, varnames = c("species", "lever"), value.name = "sym")
num_long <- melt(lever_score, varnames = c("species", "lever"), value.name = "score") %>%
  left_join(sym_long, by = c("species", "lever"))

p <- ggplot(num_long, aes(x = lever, y = species, fill = score)) +
  geom_tile() +
  geom_text(aes(label = sym), size = 6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(
    title = "Culture lever evidence (+ / - / ●)",
    x = "Culture levers",
    y = "Parabacteroides species"
  ) +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p
ggsave("Para_Culture lever evidence.pdf", p, width = 8, height = 3.5)
ggsave("Para_Culture lever evidence.png", p, width = 8, height = 3.5)

# ---------------------------------------------------------------------------
# 12. Save culture lever tables
# ---------------------------------------------------------------------------
write.csv(as.data.frame(lever_score), "Parabacteroides_species_culture_scores.csv")
write.csv(as.data.frame(sym_mat), "Parabacteroides_species_culture_symbols.csv")

# ---------------------------------------------------------------------------
# 13. DotPlot of all genes used in culture modules across species
# ---------------------------------------------------------------------------
all_genes <- unique(unlist(module_sets))

# Check merged gene sets
all_genes

p <- DotPlot(
  Parabac,
  features = all_genes,
  group.by = "species",
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
  "Parabacteroides_species_culture_genes.png",
  p, width = 15, height = 3.5, dpi = 300
)
