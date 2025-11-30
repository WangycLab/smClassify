# This script:
#   - Loads the Seurat object (seu)
#   - Annotates top 20 species (unchanged behavior)
#   - Generates FeaturePlots for 4 genes:
#         BoGH43B  → MGYG000413837-01743
#         BoGH31A  → MGYG000413837-01569
#         bepF
#         bepG
# ============================================================================

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

# ----------------------------------------------------------------------------
# 0. Project directory & data
# ----------------------------------------------------------------------------
project_dir <- "path/to/project_root"

load(file.path(project_dir, "seu_symbol_SCT_clustered_annotated.RData"))
setwd(file.path(project_dir, "Figure6"))

# ----------------------------------------------------------------------------
# Helper: unified FeaturePlot theme
# ----------------------------------------------------------------------------
fp_theme <- theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 14)
  )

# ============================================================================
#  FeaturePlots for 4 genes
#    BoGH43B → MGYG000413837-01743
#    BoGH31A → MGYG000413837-01569
#    bepF
#    bepG
# ============================================================================

# --- BoGH43B (MGYG000413837-01743) -----------------------------------------
p <- FeaturePlot(
  seu,
  features = c("MGYG000413837-01743"),
  cols     = c("lightgrey", "red"),
  pt.size  = 0.1
) + fp_theme

ggsave("FeaturePlot_seu_MGYG000413837-01743.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_MGYG000413837-01743.png", p, width = 7.5, height = 6, dpi = 300)

# --- BoGH31A (MGYG000413837-01569) -----------------------------------------
p <- FeaturePlot(
  seu,
  features = c("MGYG000413837-01569"),
  cols     = c("lightgrey", "red"),
  pt.size  = 0.1
) + fp_theme

ggsave("FeaturePlot_seu_MGYG000413837_01569.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_MGYG000413837_01569.png", p, width = 7.5, height = 6, dpi = 300)

# --- bepF -------------------------------------------------------------------
p <- FeaturePlot(
  seu,
  features = c("bepF"),
  cols     = c("lightgrey", "red"),
  pt.size  = 0.1
) + fp_theme

ggsave("FeaturePlot_seu_bepF.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_bepF.png", p, width = 7.5, height = 6, dpi = 300)

# --- bepG -------------------------------------------------------------------
p <- FeaturePlot(
  seu,
  features = c("bepG"),
  cols     = c("lightgrey", "red"),
  pt.size  = 0.1
) + fp_theme

ggsave("FeaturePlot_seu_bepG.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_bepG.png", p, width = 7.5, height = 6, dpi = 300)
