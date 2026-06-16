## ======================================================================
## Figure 6a. Feature plots of selected microbial genes


library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)

## ----------------------------------------------------------------------
## 1. Input and output paths
## ----------------------------------------------------------------------

base_dir   <- "Figure6"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


annotated_seurat_rds <- "Figure3/Output/Figure3_mouse_bacterial_functional_cluster_annotated.rds"

# ============================================================
# 3. Load annotated Seurat object
# ============================================================
seu <- readRDS(annotated_seurat_rds)

## ----------------------------------------------------------------------
## 3. Define genes shown in Figure 6a
## ----------------------------------------------------------------------

fig6a_genes <- tibble::tribble(
  ~gene_id,                ~gene_label, ~gene_class,
  "MGYG000413837-01743",   "BoGH43B",   "Glycoside hydrolase",
  "MGYG000413837-01568",   "BoGH31A",   "Glycoside hydrolase",
  "bepF",                 "bepF",      "Efflux pump",
  "bepG",                 "bepG",      "Efflux pump"
)


## ----------------------------------------------------------------------
## 4. Plotting function
## ----------------------------------------------------------------------

theme_feature <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

plot_feature_gene <- function(gene_id, gene_label) {
  FeaturePlot(
    seu,
    features = gene_id,
    cols = c("lightgrey", "red"),
    pt.size = 0.1
  ) +
    ggtitle(gene_label) +
    theme_feature
}

## ----------------------------------------------------------------------
## 5. Generate and save individual feature plots
## ----------------------------------------------------------------------

fig6a_plots <- list()

for (i in seq_len(nrow(fig6a_genes))) {
  gene_id <- fig6a_genes$gene_id[i]
  gene_label <- fig6a_genes$gene_label[i]
  
  p <- plot_feature_gene(gene_id, gene_label)
  
  fig6a_plots[[gene_label]] <- p
  
  ggsave(
    filename = file.path(output_dir, paste0("Fig6a_FeaturePlot_", gene_label, ".pdf")),
    plot = p,
    width = 7.5,
    height = 6
  )
  
  ggsave(
    filename = file.path(output_dir, paste0("Fig6a_FeaturePlot_", gene_label, ".png")),
    plot = p,
    width = 7.5,
    height = 6,
    dpi = 300
  )
}

