# This script performs the following analyses for Figure 5a, 5b:
#
# 1. Construct Sample × Species count matrix from Seurat object
# 2. Convert counts to relative abundance per sample
# 3. Compute Spearman correlation among microbial species
#    - Global heatmap of all species
#    - Focused heatmap of top 20 abundant species
#    - Extract species highly correlated with *Muribaculum gordoncarteri*
#
# 4. Derive species order from hierarchical clustering (top 20)
#    - Used downstream for consistent DotPlot ordering
#
# 5. Define microbial metabolic pathway gene sets
#    (SCFA, LCFA, nitrogen metabolism, PULs, mucin, bile acids,
#     propionate subpathways, stress & resistance modules)
#
# 6. Calculate module scores for the top 20 bacterial species
#    using Seurat::AddModuleScore
#
# 7. Visualize pathway activities across species
#    - DotPlot grouped by species
#    - Color indicates module expression level
#    - Dot size indicates percent of expressing cells
#
# ================================================================
suppressPackageStartupMessages({
  library(networkD3)
  library(dplyr)
  library(Matrix)
  library(data.table)
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
})

## 0) Basic settings: project_dir + plot_dir
project_dir <- "YOUR/PATH/HERE"   
setwd(project_dir)
plot_dir <- file.path(project_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

## 1) Load Seurat object
load("seu_symbol_SCT_clustered_annotated.RData")

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

## 1. Build Sample x Species matrix
species_abund <- as.data.frame(table(seu$orig.ident, seu$species)) %>%
  tidyr::spread(Var2, Freq)

rownames(species_abund) <- species_abund$Var1
species_abund <- species_abund[, -1]

## 2. Convert to relative abundance matrix
species_abund_rel <- sweep(species_abund, 1, rowSums(species_abund), FUN = "/")

## Calculate correlation matrix
corr_matrix <- cor(species_abund_rel, method = "spearman")

## 5. Heatmap visualization
# All species (interactive)
pheatmap(
  corr_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation Heatmap (Species Relative Abundance)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 6,
  color = colorRampPalette(c("blue", "white", "red"))(50)
)

# All species (saved)
pheatmap::pheatmap(
  corr_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation Heatmap (Species Relative Abundance)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 6,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  filename = file.path(plot_dir, "Correlation_Heatmap_Species_Relative_Abundance_pheatmap.png"),
  width = 8, height = 6
)

# High-res PNG using pheatmap (不是 ComplexHeatmap)
png(
  file.path(plot_dir, "Correlation_Heatmap_Species_Relative_Abundance.png"),
  width = 8, height = 8, units = "in", res = 600,
  type = "cairo", bg = "white"
)
pheatmap::pheatmap(
  corr_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation Heatmap (Species Relative Abundance)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 16,
  color = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

# Subset with |cor| > 0.6
sel_species <- names(cor_res[abs(cor_res) > 0.6 & names(cor_res) != "Muribaculum gordoncarteri"])
if (length(sel_species) > 0) {
  mat_sel <- species_abund_rel[, c("Muribaculum gordoncarteri", sel_species)]
  pheatmap(
    cor(mat_sel, method = "spearman"),
    cluster_rows = TRUE, cluster_cols = TRUE,
    main = "|cor| > 0.6 with Muribaculum gordoncarteri",
    fontsize = 8,
    color = colorRampPalette(c("blue", "white", "red"))(50)
  )
} else {
  message("missing |cor| > 0.6 species")
}

## Top 20 species correlation heatmap
top20_species <- c(
  "Alloprevotella sp002933955",
  "Bacteroides acidifaciens",
  "Bacteroides muris",
  "Cryptobacteroides sp910585445",
  "Duncaniella freteri",
  "Duncaniella muricolitica",
  "Duncaniella muris",
  "Eubacterium_J sp910574915",
  "Eubacterium_J sp910580075",
  "Kineothrix sp910578405",
  "Muribaculum gordoncarteri",
  "Parabacteroides distasonis",
  "Parabacteroides sp910577325",
  "Prevotella sp002933775",
  "Uncultured 14-2",
  "Uncultured 1XD42-69",
  "Uncultured CAG-873",
  "Uncultured COE1",
  "Uncultured Roseburia",
  "Uncultured UBA3402"
)

## 2. Extract data of these 20 species from abundance matrix
mat_top20 <- species_abund_rel[, intersect(top20_species, colnames(species_abund_rel)), drop = FALSE]

## 3. Compute correlation matrix
corr_top20 <- cor(mat_top20, method = "spearman")

## 4. Draw heatmap (interactive + saving)
ph <- pheatmap(
  corr_top20,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation Heatmap (Top 20 Species)",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  color = colorRampPalette(c("blue", "white", "red"))(50)
)

pheatmap::pheatmap(
  corr_top20,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation Heatmap (Top 20 Species)",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  filename = file.path(plot_dir, "Correlation_Top20_Species.pdf"),
  width = 8, height = 7
)

pheatmap::pheatmap(
  corr_top20,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation Heatmap (Top 20 Species)",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  filename = file.path(plot_dir, "Correlation_Top20_Species.png"),
  width = 2400 / 300, height = 2400 / 300
)

## Order from heatmap
row_order <- ph$tree_row$order
ordered_species <- rownames(corr_top20)[row_order]
print(ordered_species)

col_order <- ph$tree_col$order
ordered_species_col <- colnames(corr_top20)[col_order]
print(ordered_species_col)

suppressPackageStartupMessages({
  library(purrr)
  library(patchwork)
  library(gridExtra)
})

## 1) Define gene sets (names use hyphens, gene order preserved exactly)
genesets <- list(
  # Short-chain fatty acid (SCFA) metabolism
  Butyrate = c("hbd", "crt", "bcd", "etfB", "buk", "atoD"),
  Acetate  = c("pta", "ackA", "acs"),
  Propionate = c("pccB", "epi", "mutA", "mutB", "prpC", "prpB",
                 "acnA", "icd", "acsA", "ldhA"),

  # Polysaccharide utilization (PUL)
  "PUL-Core" = c(
    "susC", "susD", "susE", "susF",
    "MGYG000321753-02998", "MGYG000321753-02040",
    "MGYG000413837-01743", "MGYG000413837-01568"
  ),
  Mucin = c(
    "nanA", "nanK", "nanE", "nanT",
    "nagA", "nagB", "nagE",
    "fucI", "fucA", "fucO", "fucU", "fucP",
    "atsA"
  ),

  # Bile acid metabolism
  "Bile-Hydrolase"      = c("cbh", "bshA"),
  "Bile-Acid-Secondary" = c("baiA", "baiCD", "baiF"),
  "Bile-Total"          = c("cbh", "bshA", "tolC", "acrA", "acrB",
                            "ompR", "emrA", "norM", "mdtK", "cusB", "mexB"),

  # Protein/amino acid metabolism
  Stickland = c("grdA", "grdB", "grdC", "grdD", "fldA", "fldB"),
  "Nitrogen-Metabolism" = c(
    "gdhA", "gltB", "gltD", "glnA", "glnB", "glnD", "glnK",
    "arcA", "arcB", "argF", "argI", "argG", "argH", "speF"
  ),

  # Long-chain fatty acid (LCFA)
  "LCFA-Activation" = c(
    "fadD", "lcfB", "fadK",
    "MGYG000321702-00384", "MGYG000321710-04598", "MGYG000321739-01041",
    "MGYG000321748-01952", "MGYG000321753-00482", "MGYG000321753-01135",
    "MGYG000321759-00152", "MGYG000321759-00407", "MGYG000321767-00431",
    "MGYG000321767-03315", "MGYG000321767-04474", "MGYG000321780-01043",
    "MGYG000321780-01305", "MGYG000321860-00614", "MGYG000321881-00212",
    "MGYG000321881-00512", "MGYG000322022-00940", "MGYG000322033-00047",
    "MGYG000322051-00257", "MGYG000322051-01863", "MGYG000322086-01485",
    "MGYG000322086-01760", "MGYG000322091-00676", "MGYG000322091-01355",
    "MGYG000322227-00008", "MGYG000322400-01304", "MGYG000322419-00056",
    "MGYG000322419-00736", "MGYG000322532-02317", "MGYG000324532-00302",
    "MGYG000324532-01730", "MGYG000325900-00921", "MGYG000325900-01392",
    "MGYG000326783-00086", "MGYG000326940-00507", "MGYG000326940-01738",
    "MGYG000326940-02329", "MGYG000328850-02706", "MGYG000331430-00871",
    "MGYG000348300-00625", "MGYG000348300-00626", "MGYG000352581-00020",
    "MGYG000352581-00893", "MGYG000352581-01087", "MGYG000410738-01264",
    "MGYG000410792-00098", "MGYG000410901-02309", "MGYG000411321-01353",
    "MGYG000411368-02148", "MGYG000413837-01999", "MGYG000413837-02122",
    "MGYG000424354-01582", "MGYG000424354-02200", "MGYG000321702-01376",
    "MGYG000321777-03348", "MGYG000322033-00839", "MGYG000322033-01839",
    "MGYG000322086-01698", "MGYG000322124-01460", "MGYG000322545-02390",
    "MGYG000325507-00136", "MGYG000326481-01377", "MGYG000328854-01397",
    "MGYG000330350-00451", "MGYG000330350-00941", "MGYG000330676-02392",
    "MGYG000348300-01229", "MGYG000410792-01708", "MGYG000321701-01279",
    "MGYG000321739-00875", "MGYG000322081-01083", "MGYG000324299-00020",
    "MGYG000324299-00442", "MGYG000325351-00772", "MGYG000329649-02417",
    "MGYG000367640-01082", "MGYG000411309-01219", "MGYG000411319-00015",
    "MGYG000321701-00226", "MGYG000321748-00531", "MGYG000322081-00959",
    "MGYG000322137-01884", "MGYG000328540-02350", "MGYG000321701-00821",
    "MGYG000321860-01919", "MGYG000325799-01281", "MGYG000326481-01570",
    "MGYG000328854-01260", "MGYG000328854-01608", "MGYG000321707-00023",
    "MGYG000328540-00790", "MGYG000413837-00394", "MGYG000326609-00586"
  ),
  "LCFA-Auxiliary" = c("fadD3", "acsA", "acsC", "acsD", "acsE",
                       "accA", "accB", "accC", "accD", "fadM"),
  "LCFA-Beta-Oxidation" = c("hbd", "crt", "mutB", "epi", "etfB", "fadN"),
  "LCFA-Transport-Efflux" = c("acrB", "bepC", "bepF", "bepG", "tolC",
                              "mexB", "mdtA", "mdtB", "mdtC"),
  "LCFA-Regulation" = c("fadR", "fnr", "crp", "arcA", "arcB"),
  "LCFA-Other-Aux" = c("caiA", "caiT", "acdA", "pflB", "fixA", "fixB", "fixX"),

  # Stress response / resistance
  Stress = c("sodA", "sodB", "sodC", "katA", "katE", "rbr",
             "dnaK", "groL", "groS", "clpB", "clpP", "grpE"),
  "Efflux-Resistance" = c("tolC", "acrA", "acrB", "ompR", "emrA",
                          "norM", "mdtK", "cusB", "mexB"),

  # Propionate sub-pathways
  "Propionate-MCM" = c("pccB", "epi", "mutA", "mutB"),
  "Propionate-MCC" = c("prpC", "prpB", "acnA", "icd",
                       "emrA", "norM", "mdtK", "cusB", "mexB"),
  "Propionate-Acrylate" = c("acsA", "ldhA")
)

# Filter gene sets by genes present in Seurat object
all_genes <- rownames(seu)
genesets_clean <- lapply(genesets, function(g) intersect(g, all_genes))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]
genesets_clean

## Re-define top20_species (for safety)
top20_species <- c(
  "Alloprevotella sp002933955",
  "Bacteroides acidifaciens",
  "Bacteroides muris",
  "Cryptobacteroides sp910585445",
  "Duncaniella freteri",
  "Duncaniella muricolitica",
  "Duncaniella muris",
  "Eubacterium_J sp910574915",
  "Eubacterium_J sp910580075",
  "Kineothrix sp910578405",
  "Muribaculum gordoncarteri",
  "Parabacteroides distasonis",
  "Parabacteroides sp910577325",
  "Prevotella sp002933775",
  "Uncultured 14-2",
  "Uncultured 1XD42-69",
  "Uncultured CAG-873",
  "Uncultured COE1",
  "Uncultured Roseburia",
  "Uncultured UBA3402"
)

seu@active.ident <- as.factor(seu$species)
sel_seu  <- subset(seu, idents = top20_species)
table(sel_seu$species)

# Scoring
sel_seu <- AddModuleScore(
  sel_seu,
  features = genesets_clean,
  name = names(genesets_clean),
  nbin = 6
)

colnames(sel_seu@meta.data)

# Use species order from correlation heatmap (top20)
ordered_species <- rev(ordered_species)
print(ordered_species)

sel_seu$species <- factor(sel_seu$species, levels = ordered_species)

## Define module meta columns used in DotPlot
modules <- c(
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
)

p <- DotPlot(
  sel_seu,
  features = modules,
  group.by = "species",
  scale = TRUE
) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3",
    midpoint = 0
  ) +
  labs(
    title = "Key Metabolic and Stress Pathways in Top 20 Species",
    x = "Pathways", y = "Species"
  )

p

ggsave(
  file.path(plot_dir, "Key_Metabolic_and_Stress_Pathways_top20_species.png"),
  p, width = 10, height = 7.5, dpi = 600
)

ggsave(
  file.path(plot_dir, "Key_Metabolic_and_Stress_Pathways_top20_species.pdf"),
  p, width = 10, height = 7.5
)
