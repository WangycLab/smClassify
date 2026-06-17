############################################################
## Figure 5a,b species correlation and functional pathway dotplot
##

############################################################

# ============================================================
# 1. Load packages
# ============================================================
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(Matrix)
library(scales)
library(grid)


# ============================================================
# 2. Input and output paths
# ============================================================
base_dir   <- "Figure5"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

annotated_seurat_rds <- "Figure3/Output/Figure3_mouse_bacterial_functional_cluster_annotated.rds"

# ============================================================
# 3. Load annotated Seurat object
# ============================================================
seu <- readRDS(annotated_seurat_rds)


# ============================================================
# 4. Define top 20 species
# ============================================================
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

# ============================================================
# 5. Figure 5a: top 20 species correlation heatmap
# ============================================================
species_abund <- as.data.frame(table(seu$orig.ident, seu$species)) %>%
  tidyr::pivot_wider(
    names_from = Var2,
    values_from = Freq,
    values_fill = 0
  )

rownames(species_abund) <- species_abund$Var1
species_abund <- species_abund[, -1, drop = FALSE]

species_abund_rel <- sweep(
  species_abund,
  1,
  rowSums(species_abund),
  FUN = "/"
)

top20_present <- intersect(top20_species, colnames(species_abund_rel))

mat_top20 <- species_abund_rel[, top20_present, drop = FALSE]
corr_top20 <- cor(mat_top20, method = "spearman")

write.csv(
  mat_top20,
  file = file.path(output_dir, "Figure5a_top20_species_relative_abundance_source_data.csv")
)

write.csv(
  corr_top20,
  file = file.path(output_dir, "Figure5a_top20_species_spearman_correlation_source_data.csv")
)

ph <- pheatmap::pheatmap(
  corr_top20,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation Heatmap (Top 20 Species)",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  filename = file.path(output_dir, "Figure5a_top20_species_correlation_heatmap.pdf"),
  width = 8,
  height = 7
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
  filename = file.path(output_dir, "Figure5a_top20_species_correlation_heatmap.png"),
  width = 8,
  height = 8
)

ordered_species <- rownames(corr_top20)[ph$tree_row$order]

write.csv(
  data.frame(ordered_species = ordered_species),
  file = file.path(output_dir, "Figure5a_top20_species_heatmap_row_order.csv"),
  row.names = FALSE
)

# ============================================================
# 6. Define functional gene modules
# ============================================================
genesets <- list(
  "Butyrate Production" = c("hbd", "crt", "bcd", "etfB"),
  "Acetate Production" = c("pta", "ackA", "acs"),
  "Propionate Metabolism" = c("pccB", "epi", "mutA", "mutB", "prpC", "prpB", "acnA", "icd", "acsA", "ldhA"),
  
  "PUL Core Complex" = c("susC", "susD", "susE"),
  "Mucin Degradation" = c("nagA", "nagB", "nagE", "fucI", "fucA", "fucO", "fucU", "fucP", "atsA"),
  
  "Bile Salt Hydrolase" = c("cbh"),
  "Bile Acid Transport" = c("cbh", "tolC", "acrA", "acrB", "ompR", "mdtK", "cusB", "mexB"),
  
  "Stickland Fermentation" = c("fldA"),
  "Nitrogen Assimilation" = c("gdhA", "gltB", "gltD", "glnA", "glnB", "glnD", "glnK", "arcA", "arcB", "argF", "argG", "argH"),
  
  "LCFA Activation" = c("fadD", "lcfB"),
  "LCFA Elongation" = c("fadD3", "acsA", "accA", "accB", "accC", "accD"),
  "LCFA Beta-Oxidation" = c("hbd", "crt", "mutB", "epi", "etfB", "fadN"),
  "LCFA Transport & Efflux" = c("acrB", "bepF", "bepG", "tolC", "mexB", "mdtA", "mdtB", "mdtC"),
  "LCFA Regulation" = c("crp", "arcA", "arcB"),
  "LCFA Supporting Pathways" = c("caiA", "caiT", "pflB"),
  
  "Stress Response" = c("sodB", "katA", "rbr", "dnaK", "groL", "groS", "clpB", "clpP", "grpE"),
  "Efflux & Antibiotic Resistance" = c("tolC", "acrA", "acrB", "ompR", "mdtK", "cusB", "mexB"),
  
  "Propionate (MCM Pathway)" = c("pccB", "epi", "mutA", "mutB"),
  "Propionate (MCC Pathway)" = c("prpC", "prpB", "acnA", "icd", "mdtK", "cusB", "mexB"),
  "Propionate (Acrylate Pathway)" = c("acsA", "ldhA")
)

genesets_clean <- lapply(genesets, function(g) intersect(g, rownames(seu)))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]

write.csv(
  data.frame(
    Module = rep(names(genesets_clean), lengths(genesets_clean)),
    Gene = unlist(genesets_clean, use.names = FALSE)
  ),
  file = file.path(output_dir, "Figure5b_functional_gene_modules_used.csv"),
  row.names = FALSE
)

# ============================================================
# 7. Figure 5b: pathway module-score dotplot
# ============================================================
seu@active.ident <- as.factor(seu$species)

sel_seu <- subset(seu, idents = top20_present)

sel_seu <- AddModuleScore(
  sel_seu,
  features = genesets_clean,
  name = names(genesets_clean),
  nbin = 6
)

ordered_species <- rev(ordered_species)
sel_seu$species <- factor(sel_seu$species, levels = ordered_species)

module_order <- c(
  "Acetate Production2",
  "Propionate (MCM Pathway)18",
  "Propionate (MCC Pathway)19",
  "Propionate (Acrylate Pathway)20",
  "Butyrate Production1",
  "LCFA Beta-Oxidation12",
  "LCFA Supporting Pathways15",
  "LCFA Regulation14",
  "LCFA Elongation11",
  "LCFA Transport & Efflux13",
  "LCFA Activation10",
  "PUL Core Complex4",
  "Nitrogen Assimilation9",
  "Mucin Degradation5",
  "Bile Acid Transport7",
  "Bile Salt Hydrolase6",
  "Stress Response16",
  "Efflux & Antibiotic Resistance17",
  "Stickland Fermentation8"
)

module_order <- module_order[module_order %in% colnames(sel_seu@meta.data)]

write.csv(
  sel_seu@meta.data[, c("orig.ident", "species", module_order), drop = FALSE],
  file = file.path(output_dir, "Figure5b_top20_species_module_scores_source_data.csv")
)

p_5b <- DotPlot(
  sel_seu,
  features = module_order,
  group.by = "species",
  scale = TRUE
) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red3",
    midpoint = 0
  ) +
  labs(
    title = "Key Metabolic and Stress Pathways in Top 20 Species",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(output_dir, "Figure5b_top20_species_key_metabolic_stress_pathway_dotplot.png"),
  p_5b,
  width = 10,
  height = 8,
  dpi = 600
)

ggsave(
  file.path(output_dir, "Figure5b_top20_species_key_metabolic_stress_pathway_dotplot.pdf"),
  p_5b,
  width = 10,
  height = 8
)
