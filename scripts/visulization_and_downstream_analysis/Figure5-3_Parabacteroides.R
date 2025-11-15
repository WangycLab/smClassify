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
# Set working directory
load("seu_symbol_SCT_clustered_annotated.RData")

setwd("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009/Figure5")

nature_palette <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85"
)
cell_type_cols <- c("#1F77B4",  "#2CA02C","#FF7F0E","#6A5ACD", "#8C564B", "#E377C2", "#7F7F7F", 
                    "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32",
                    "#8B0000",  "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B", 
                    "#20B2AA", "#FF1493", "#FF4500", "#32CD32", "#3E8E41", "#20B2AA")


sample<-seu
sample@active.ident <- as.factor(sample$Tax_Genus)
Parabac <- subset(sample, idents = "Parabacteroides")
table(Parabac$species)


Parabac <- SCTransform(Parabac, verbose = FALSE,vars.to.regress = "nCount_RNA")
Parabac <- RunPCA(Parabac, verbose = FALSE)
Parabac <- FindNeighbors(Parabac, dims = 1:5)
Parabac <- FindClusters(Parabac, resolution = 0.15)
Parabac <- RunUMAP(Parabac, dims = 1:5)

cell_type_cols<-c( "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32",
                  "#8B0000",  "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B", 
                  "#FF1493", "#FF4500", "#32CD32", "#3E8E41", "#20B2AA")
p<-DimPlot(Parabac, reduction = "umap", label = F, pt.size = 0.1,
        group.by = "seurat_clusters") +
  scale_color_manual(values = cell_type_cols)+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("Para_DimPlot_seurat_clusters.png", p, width = 4, height = 3, dpi = 300)

p<-DimPlot(Parabac, reduction = "umap", label = F, pt.size = 0.1,
           group.by = "phenotype") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("Parabac_DimPlot_seurat_phenotype.png", p, width = 4, height = 3, dpi = 300)

saveRDS(Parabac,"Parabac.Rds")
Parabac<-readRDS("Parabac.Rds")

p<-DimPlot(Parabac, reduction = "umap", label = F, pt.size = 0.1,
           group.by = "region") +
  theme_minimal() +
  scale_color_manual(values=c( "Cecum" = "#1f77b4",
                               "Colon" = "#2ca02c",
                               "Rectum" = "#ff7f0e"))+
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 14)) 
p
ggsave("Parabac_DimPlot_seurat_region.png", p, width = 4.5, height = 3, dpi = 300)

p<-DimPlot(Parabac, reduction = "umap", label = F, pt.size = 0.1,
           group.by = "phenotype") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("Parabac_DimPlot_seurat_phenotype.png", p, width = 4, height = 3, dpi = 300)


cell_type_cols <- c( "#BCBD22","#228B22","#7F7F7F")
p<-DimPlot(Parabac, reduction = "umap", label = F, pt.size = 0.5,
           group.by = "species") +
  theme_minimal() +
  scale_color_manual(values=cell_type_cols)+
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("Parabac_DimPlot_species.png", p, width = 6, height = 3, dpi = 300)
ggsave("Parabac_DimPlot_species.pdf", p, width = 6, height = 3)


DotPlot(
  Parabac,
  features = c("mutB","epi"),
  group.by = "group"   
) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text  = element_text(face = "bold")
  ) +
  labs(title = "BCAA-related genes - Parabacteroides distasonis")

## ------------------------------------------------------------
## Figure 1. UMAP overview - Acetyl-phosphate metabolism dataset
## ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)

## 1. Basic Parameters
# Paper layout size
fig_width  <- 10             # inches
fig_height <- 7             # square
dpi_out    <- 300              # PNG resolution
base_size  <- 10                # base font size for annotations/ticks (8 pt)

## 2. Color Scheme
cell_type_cols <- c(
  "#1F77B4", "#2CA02C", "#FF7F0E", "#F0E68C", "#6A5ACD", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#F5A0A1",
  "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000", "#FFD700",
  "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B", "#20B2AA",
  "#FF1493", "#FF4500", "#32CD32", "#3E8E41"
)

## 3. General Themes
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
# Set global theme
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

Idents(Parabac) <- Parabac$seurat_clusters      

# 4. Differential Gene Screening + Annotation
message("Finding markers")
markers <- FindAllMarkers(
  object          = Parabac,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.15
) 

markers


# 5. Plot Heatmap
# 1. Gene and expression data preparation
markers$gene
topN <- markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 50) |>
  ungroup()
topN
markers


a <- read.csv('Parabac_type.csv')
a
levels(Parabac)
new.cluster.ids <- as.character(a[,2])
new.cluster.ids
levels(Parabac)
Parabac@active.ident <- as.factor(Parabac$seurat_clusters)
names(new.cluster.ids) <- levels(Parabac)
Parabac <- RenameIdents(Parabac, new.cluster.ids)
Parabac$cellsubtype <- Parabac@active.ident


cell_type_cols <- c(
  "#8C564B", "#2CA02C", "#F5A0A1", "#FF7F0E",  "#6A5ACD", "#8B0000","#C2B5D8",
  "#1F77B4",  "#BCBD22", "#17BECF",
  "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000", "#FFD700",
  "#DC143C", "#228B22", "#483D8B", "#BDB76B", "#20B2AA",
  "#FF1493", "#FF4500", "#32CD32", "#3E8E41")


cell_type_cols2 <- c( "#6A5ACD","#F5A0A1", "#3CB371", "#FFD700")

p<-DimPlot(Parabac, reduction = "umap", label = F, pt.size = 0.5,
           group.by = "cellsubtype") +
  scale_color_manual(values =  cell_type_cols2)+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 

p
ggsave( "Parabac_DimPlot_celltype.pdf", p, w=7, h=3)



# Set it as the Seurat object identity
Idents(Parabac) <- "seurat_clusters"

top5_nonMGYG <- topN %>%
  dplyr::filter(!grepl("^MGYG", gene)) %>%      # Remove genes starting with MGYG
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%       # Select top 5 per cluster
  ungroup()
top5_nonMGYG 
genes <- unique(top5_nonMGYG$gene)
genes

# Plotting
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
  labs(x = "", y = "Functional Subtypes", title = "DotPlot of Top Marker Genes")+
  scale_size_continuous(range = c(1, 8)) 
print(p)

ggsave("Parabac_dotplot_topgenes.pdf", p, width = 12, height = 3)
ggsave("Parabac_dotplot_topgenes.png", p, width = 12, height = 3, dpi = 300)




# Compare culture requirements among Parabacteroides species

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)

# ----------------------------
# 1. Define modules related to culture requirements
# ----------------------------
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

# ----------------------------
# 0) Configuration: module/gene keywords & weights
# ----------------------------
# 1) Original functional module gene sets
# ----------------------------
culture_modules <- list(
  # 1. Inulin/Amylose (plant polysaccharide utilization; removed missing susE/susF/gh13/gh32/amyA)
  Inulin_Amylose = c("susC","susD","malQ","glgA"),
  
  # 2. Mucin/NAG (removed missing nanK/nanE/nanH/nagE)
  Mucin_NAG = c("nanA","nagA","nagB","nagC"),
  
  # 3. Fucose (removed missing fucA/fucK/fucU)
  Fucose = c("fucI","fucO"),
  
  # 4. Xylose / Mannose / Galactose (removed missing xylF/manB/galM)
  Xylose_Man_Gal = c("xylA","xylB","xylE","manA","galE","galK"),
  
  # 5. Ethanolamine (all missing in this object: eutA/B/C/D/E/G/H empty)
  Ethanolamine = character(0),
  
  # 6. Vitamin B12 (removed missing cobA/cobC/cobI/cobU/cbiB/cbiC/cbiT)
  Vitamin_B12 = c("cobB","cobD","cbiA","cbiD"),
  
  # 7. Fumarate (electron acceptor) (removed missing frdC/frdD)
  Fumarate_acceptor = c("frdA","frdB"),
  
  # 8. Bile salts (removed missing bsh/acrA/acrB/ompF/ompC)
  Bile_salts = c("cbh","acrF","tolC","ompR"),
  
  # 9. Branched SCFAs (removed missing bkdA/bkdB/bkdC/leuD)
  Branched_SCFAs = c("ilvE","leuA","leuB","leuC"),
  
  # 10. Vitamins (broad)
  Vitamins_broad = c(
    # Biotin
    "bioA","bioB","bioC","bioD","bioF",
    # Thiamine (kept)
    "thiC","thiD","thiE",
    # Folate (kept except folA)
    "folB","folC","folD",
    # NAD
    "nadA","nadB","nadC","nadD","nadE",
    # Menaquinone (kept except menC/menG)
    "menA","menB","menD","menE","menF"
  )
)
# ----------------------------
# 2) Gene keywords (for ��gene-level�� evidence; filtered)
# ----------------------------
gene_keywords <- list(
  Mucin_NAG = c("nagA","nagB","nanA"),
  Fucose    = c("fucI","fucP"),
  XMG       = c("xylA","xylB","manA","galE","galK"),
  Prop_core = c("mutA","mutB","pccA","pccB","mcmA","mcmB","prpB","prpC","prpD"),
  B12_core  = c("cobB","cobD","cbiA","cbiD")
)
lever_module_weights <- list(
  "Mucin/NAG" = c(PUL_sugars = 1),
  "Fucose" = c(PUL_sugars = 1),
  "Xylose/Man/Gal" = c(PUL_sugars = 1),
  "Fumarate (e- acceptor)" = c(Propionate = 1),
  "Ethanolamine" = c(Ethanolamine_B12 = 1),
  "Vitamin B12" = c(Ethanolamine_B12 = 1),
  "Bile salts (tolerate)" = c(Bile_tolerance = 1)
)
# ----------------------------
# 3) Summarize into 6 module score feature sets (consistent with lever_module_weights; based on filtered sets)
# ----------------------------
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

# 3) Use Seurat AddModuleScore2 to generate *_score1 columns one by one
#    Each name generates a column name1 (i.e., *_score1)
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

# ----------------------------
# 2) Aggregate ModuleScore by species
# ----------------------------
scores_species <- Parabac@meta.data %>%
  group_by(species) %>%
  summarise(across(ends_with("_score1"), mean, na.rm = TRUE)) %>%
  arrange(species)

scores_species

mat_mod <- as.matrix(scores_species[,-1])
rownames(mat_mod) <- scores_species$species
colnames(mat_mod) <- gsub("_score1$", "", colnames(mat_mod))
z_mod <- scale(mat_mod)

# ----------------------------
# 3) Gene-level scoring
# ----------------------------
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
  avg <- AverageExpression(Parabac, features = feats, group.by = "species", assays = DefaultAssay(Parabac))[[1]]
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

# ----------------------------
# 4) Merge into culture lever scores
# ----------------------------
lever_score <- matrix(0, nrow = nrow(z_mod), ncol = length(lever_module_weights))
rownames(lever_score) <- rownames(z_mod)
colnames(lever_score) <- names(lever_module_weights)

for (lv in names(lever_module_weights)) {
  w <- lever_module_weights[[lv]]
  mod_cols <- intersect(names(w), colnames(z_mod))
  mod_part <- 0
  if (length(mod_cols) > 0) {
    mod_part <- as.matrix(z_mod[, mod_cols, drop = FALSE]) %*% matrix(w[mod_cols], ncol = 1)
    mod_part <- as.numeric(mod_part)
  }
  gene_part <- rep(0, nrow(z_mod))
  if (!is.null(z_gene)) {
    if (lv == "Mucin/NAG" && "Mucin_NAG" %in% colnames(z_gene)) gene_part <- z_gene[, "Mucin_NAG"]
    if (lv == "Fucose"    && "Fucose"    %in% colnames(z_gene)) gene_part <- z_gene[, "Fucose"]
    if (lv == "Xylose/Man/Gal" && "XMG"  %in% colnames(z_gene)) gene_part <- z_gene[, "XMG"]
    if (lv == "Fumarate (e- acceptor)" && "Prop_core" %in% colnames(z_gene)) gene_part <- z_gene[, "Prop_core"]
    if (lv %in% c("Ethanolamine","Vitamin B12") && "B12_core" %in% colnames(z_gene)) gene_part <- z_gene[, "B12_core"]
  }
  gene_part[is.na(gene_part)] <- 0
  if (all(gene_part == 0)) {
    lever_score[, lv] <- mod_part
  } else {
    lever_score[, lv] <- 0.5 * mod_part + 0.5 * gene_part
  }
}

# Additional derived column: avoid bile salts
lever_score <- cbind(lever_score,
                     "Avoid bile" = -lever_score[, "Bile salts (tolerate)"])

# ----------------------------
# 5) Symbol matrix
# ----------------------------
to_symbol <- function(x, pos_thr=0.5, neg_thr=-0.5) {
  ifelse(x > pos_thr, "+", ifelse(x < neg_thr, "-", "●"))
}
sym_mat <- apply(lever_score, 2, to_symbol)
sym_mat 
# ----------------------------
# 6) Visualization
# ----------------------------
# Symbol matrix
sym_long <- melt(sym_mat, varnames = c("species","lever"), value.name = "sym")
num_long <- melt(lever_score, varnames = c("species","lever"), value.name = "score") %>%
  left_join(sym_long, by = c("species","lever"))
p<-ggplot(num_long, aes(x = lever, y = species, fill = score)) +
  geom_tile() +
  geom_text(aes(label = sym), size = 6) +
  scale_fill_gradient2(low = "blue",  mid = "white",high = "red") +
  labs(title = "Culture lever evidence (+ / - / ●)",
       x = "Culture levers", y = "Parabacteroides species") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("Para_Culture lever evidence.pdf", p, width = 8, height = 3.5)
ggsave("Para_Culture lever evidence.png", p, width = 8, height = 3.5)

# ----------------------------
# 7) Save results
# ----------------------------
write.csv(as.data.frame(lever_score), "Parabacteroides_species_culture_scores.csv")
write.csv(as.data.frame(sym_mat), "Parabacteroides_species_culture_symbols.csv")


all_genes <- unique(unlist(module_sets))

# Check merged gene sets
all_genes
p<-DotPlot(Parabac, features = all_genes,
           group.by = "species",
           scale = T) + 
  RotatedAxis()+
  scale_size_continuous(range = c(2, 8))+
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3",
    midpoint = 0  
  ) 
p
ggsave("Parabacteroides_species_culture_genes.png", p, 
       width = 15, height = 3.5, dpi = 300)
