#   This script loads the processed Seurat object (SCT-normalized, clustered, and annotated at species level), computes module scores for multiple
#   microbial metabolic pathways, and visualizes pathway activity using
#   DotPlots across:
#       - host groups (WT vs DB)
#       - top bacterial species (top20_species)
#
# Pathways included:
#   1) Nitrogen metabolism:
#       - Glutamine–glutamate axis
#       - Ornithine/citrulline/arginine interconversion
#       - Polyamine-related enzymes
#
#   2) Bile acid–related genes:
#       - Deconjugation (cbh)
#       - Efflux pumps (acrA, mdtK)
#       - Osmotic/bile-stress regulators (ompR)
#       - Oxidative/metal-stress (cusB)
#
#   3) Long-chain fatty acid (LCFA) related genes:
#       - Uptake/activation (fadD)
#       - Transport (tolC)
#       - β-oxidation–related enzymes (fabF, crt, hbd)
#
# ======================================================
# ---- 0. Load packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(glue)
  library(ggrepel)
  library(ggpubr)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ---- 1. Paths & data ----
# Set your project directory to the Figure4 folder
project_dir <- "YOUR/PATH/HERE"   # e.g. "D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009/Figure4"
setwd(project_dir)

# Plot output directory
plot_dir <- file.path(project_dir, "plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
}

# Load Seurat object (stored one level above Figure4)
seurat_file <- file.path(project_dir, "..", "seu_symbol_SCT_clustered_annotated.RData")
load(seurat_file)

table(seu$species)

# ---- 2. Define Nitrogen Metabolism Gene Sets ----
NitrogenMetabolismGenesets <- list(
  Glutamate                = c("gdhA", "gltB", "gltD"),
  Glutamine                = c("glnA", "glnB", "glnD", "glnK"),
  Ornithine_to_Citrulline  = c("arcB", "argF", "argI"),
  Citrulline_to_Arginine   = c("argG", "argH"),
  Arginine_to_Citrulline   = c("arcA"),
  Ornithine_to_Polyamine   = c("speF")
)

seu$group

# ---- 3. Clean gene sets to match Seurat object ----
all_genes <- rownames(seu)
genesets_clean <- lapply(NitrogenMetabolismGenesets, function(g) intersect(g, all_genes))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]
genesets_clean

# ---- 4. Calculate Module Scores ----
seuAddModule <- AddModuleScore(
  object   = seu,
  features = genesets_clean,
  name     = names(genesets_clean),
  nbin     = 24,
  seed     = 123
)

# ---- 5. Check meta.data columns ----
colnames(seuAddModule@meta.data)

# ======================================================
# 6. Visualization
# ======================================================

## Helper function: general DotPlot style
theme_dot <- function() {
  theme_bw(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 15),
      axis.text.y  = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      strip.text   = element_text(face = "bold")
    )
}

# ---- 6.4. DotPlot by Group for individual genes ----
## Glutamate / Glutamine / Arginine / Ornithine genes
p_dot <- DotPlot(
  seuAddModule,
  features = c(
    # Arginine / Citrulline
    "arcA", "argG", "argH",
    # Citrulline / Ornithine
    "arcB", "argI", "argF",
    # Ornithine / upstream Glutamate-Glutamine
    "glnA", "glnB", "glnD", "glnK",
    "gdhA", "gltB", "gltD"
  ),
  group.by = "group",
  scale = TRUE
) +
  coord_flip() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3", midpoint = 0
  ) +
  scale_size_continuous(range = c(2, 8)) +
  theme_dot() +
  ggtitle("Amino acid metabolism")

p_dot
ggsave(
  file.path(plot_dir, "Dotplot_Amino acid metabolism by Group.pdf"),
  p_dot, width = 5, height = 5.5
)

p_dot <- DotPlot(
  seuAddModule,
  features = c(
    # Arginine / Citrulline
    "argH",
    # Citrulline / Ornithine
    "argF",
    # Ornithine / upstream Glutamate-Glutamine
    "glnA",
    "gdhA", "gltB", "gltD"
  ),
  group.by = "top20_species",
  scale = TRUE
) +
  coord_flip() +
  scale_color_gradient(
    low = "white", high = "red3"
  ) +
  scale_size_continuous(range = c(1, 8)) +
  theme_dot() +
  ggtitle("Amino acid metabolism")

p_dot
ggsave(
  file.path(plot_dir, "Dotplot_Amino acid metabolism by top20_species.pdf"),
  p_dot, width = 9, height = 5
)

# ------------------------------------------------------
# A) Bile acid related genes
# ------------------------------------------------------
genes_bile <- c(
  "cbh",   # Directly deconjugates bile acids (including taurine/glycine conjugates of TDCA/TCDCA).
  "ompR",  # Global response regulator, involved in osmotic/outer membrane/bile salt resistance networks; species-specific direction.
  "acrA",  # Core complex of RND efflux pump
  "mdtK",  # RND/MFS/MATE efflux pump family
  "cusB"   # Cu efflux system, mainly reflecting metal/oxidative stress
)

p_dot <- DotPlot(
  seuAddModule,
  features = genes_bile,
  group.by = "group",
  scale = TRUE
) +
  coord_flip() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3", midpoint = 0
  ) +
  scale_size_continuous(range = c(4, 8)) +
  theme_dot() +
  ggtitle("Bile acid related genes") +
  scale_x_discrete(limits = rev)

p_dot
ggsave(
  file.path(plot_dir, "Dotplot_Bile acid.pdf"),
  p_dot, width = 5.5, height = 4
)

p_dot <- DotPlot(
  seuAddModule,
  features = genes_bile,
  group.by = "top20_species",
  scale = TRUE
) +
  coord_flip() +
  scale_color_gradient2(
    low = "white", high = "red3"
  ) +
  scale_size_continuous(range = c(1, 8)) +
  theme_dot() +
  ggtitle("Bile acid related genes") +
  scale_x_discrete(limits = rev) +
  theme(plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "mm"))

p_dot
ggsave(
  file.path(plot_dir, "Dotplot_Bile acid related genes_top20_species.pdf"),
  p_dot, width = 12, height = 6
)
ggsave(
  file.path(plot_dir, "Dotplot_Bile acid related genes_top20_species.png"),
  p_dot, width = 12, height = 6, dpi = 300
)

# ------------------------------------------------------
# B) LFCA related genes
# ------------------------------------------------------
genelist <- c(
  "fadD",  # LCFA uptake/activation
  "tolC",  # Transport/regulation
  "fabF",
  "crt", "hbd"
)

p_dot <- DotPlot(
  seuAddModule,
  features = genelist,
  group.by = "group",
  scale = TRUE
) +
  coord_flip() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3", midpoint = 0
  ) +
  scale_size_continuous(range = c(4, 8)) +
  theme_dot() +
  ggtitle("LFCA related genes") +
  scale_x_discrete(limits = rev)

p_dot
ggsave(
  file.path(plot_dir, "Dotplot_LFCA related genes.pdf"),
  p_dot, width = 5.5, height = 4
)

p_dot <- DotPlot(
  seuAddModule,
  features = c(
    "fadD", # LCFA uptake/activation
    "tolC", # Transport/regulation
    "fabF",
    "crt", "hbd"
  ),
  group.by = "top20_species",
  scale = TRUE
) +
  coord_flip() +
  scale_color_gradient2(
    low = "white", high = "red3"
  ) +
  scale_size_continuous(range = c(1, 8)) +
  theme_dot() +
  ggtitle("LFCA related genes") +
  scale_x_discrete(limits = rev) +
  theme(plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "mm"))

p_dot
ggsave(
  file.path(plot_dir, "Dotplot_LFCA related genes_top20_species.pdf"),
  p_dot, width = 12, height = 6
)
ggsave(
  file.path(plot_dir, "Dotplot_LFCA related genes_top20_species.png"),
  p_dot, width = 12, height = 6, dpi = 300
)
