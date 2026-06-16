## ======================================================================
## Figure 6h-i. Dot plots of nitrogen metabolism genes

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)

## ----------------------------------------------------------------------
## 1. Input and output paths
## ----------------------------------------------------------------------

base_dir   <- "Figure6"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)



murigor_rds <- file.path(
  output_dir,
  "Fig6b_g_Murigor_processed.rds"
)

annotated_seurat_rds <- "Figure3/Output/Figure3_mouse_bacterial_functional_cluster_annotated.rds"

# ============================================================
# 3. Load annotated Seurat object
# ============================================================
seu <- readRDS(annotated_seurat_rds)

## ----------------------------------------------------------------------
## 3. Prepare species order for Figure 6h
## ----------------------------------------------------------------------

sp_counts <- sort(table(seu$species), decreasing = TRUE)

n_top_species <- min(20, length(sp_counts))
top_species <- names(sp_counts)[seq_len(n_top_species)]

seu$top_species <- ifelse(
  seu$species %in% top_species,
  as.character(seu$species),
  "Others"
)

seu$top_species <- factor(
  seu$top_species,
  levels = c(top_species, "Others")
)

## ----------------------------------------------------------------------
## 4. Define urease genes
## ----------------------------------------------------------------------

fig6h_genes <- c(
  "ureG", "ureF","ureE","ureD", "ureC","ureB","ureA"
)

fig6i_genes <- c(
  "ureB", "ureC", "ureD", "ureE", "ureF", "ureG"
)

fig6h_genes <- intersect(fig6h_genes, rownames(seu))

## ----------------------------------------------------------------------
## 5. Figure 6h: ureA-ureG across species
## ----------------------------------------------------------------------

p_fig6h <-DotPlot (seu,
            features = fig6h_genes  , 
            group.by = "top20_species",
            scale = TRUE) +
  scale_color_gradient(low =  "white", high = "red3") +
  coord_flip() +
  scale_size_continuous(range = c(0, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold")
  ) 
p_fig6h
ggsave(
  file.path(output_dir, "Fig6h_urease_genes_across_species.pdf"),
  p_fig6h,
  width = 8,
  height = 6
)

ggsave(
  file.path(output_dir, "Fig6h_urease_genes_across_species.png"),
  p_fig6h,
  width = 8,
  height = 6,
  dpi = 300
)

fig6h_source <- p_fig6h$data

write.csv(
  fig6h_source,
  file.path(output_dir, "Fig6h_urease_genes_across_species_source_data.csv"),
  row.names = FALSE
)

## ----------------------------------------------------------------------
## 6. Prepare M. gordoncarteri object for Figure 6i
## ----------------------------------------------------------------------


Murigor <- readRDS(murigor_rds)


fig6i_genes <- intersect(fig6i_genes, rownames(Murigor))


Murigor$group <- factor(
  Murigor$group,
  levels = c("Cecum-WT", "Cecum-DB", "Colon-WT", "Colon-DB", "Rectum-WT", "Rectum-DB")
)

## ----------------------------------------------------------------------
## 7. Figure 6i
## ----------------------------------------------------------------------

p_fig6i <- DotPlot(
  Murigor,
  features = c(
    "ureG", "ureF", "ureE", "ureD", "ureC", "ureB",
    "bepF", "bepG",
    "MGYG000413837-01569",
    "MGYG000413837-01743",
    "pulA",
    "susD",
    "susC"
  ),
  group.by = "group",
  scale = TRUE
) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red3",
    midpoint = 0
  ) +
  coord_flip() +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Polysaccharide utilization/\nStress adaptation/\nNitrogen buffering")

p_fig6i

ggsave(
  file.path(output_dir, "Fig6i_Murigor_polysaccharide_stress_nitrogen_group.pdf"),
  p_fig6i,
  width = 6.5,
  height = 5
)

ggsave(
  file.path(output_dir, "Fig6i_Murigor_polysaccharide_stress_nitrogen_group.png"),
  p_fig6i,
  width = 6.5,
  height = 5,
  dpi = 300
)

