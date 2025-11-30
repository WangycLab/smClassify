# This script was used to generate Figure 3f.
#   Visualize, for each (celltype, species) pair, the composition
#   of WT vs DB as a "split pie" (polar arc) where:
#     - Each dot sits on a 2D grid: x = species, y = celltype.
#     - Dot radius encodes total abundance of that (celltype, species).
#     - The dot is split into arcs by phenotype (WT / DB).
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggforce)
  library(tibble)
  library(forcats)
})

# ------------------------------------------------------------
# 0. Project directory and data loading
# ------------------------------------------------------------

project_dir <- "path/to/project_root"   # <-- modify for your local environment
setwd(project_dir)

load("seu_symbol_SCT_clustered_annotated.RData")

# Ensure required metadata exists
stopifnot(
  all(c("celltype", "species", "phenotype") %in% colnames(seu@meta.data))
)

# Identity = celltype
Idents(seu) <- seu$celltype

# ============================================================
# 1. Compute per-celltype/species stats (counts, props, angles)
# ============================================================

data_pie <- seu@meta.data %>%
  dplyr::group_by(celltype, species, phenotype) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(celltype, species) %>%
  dplyr::mutate(
    total = sum(n),
    prop  = n / total
  ) %>%
  dplyr::arrange(celltype, species, phenotype) %>%
  dplyr::mutate(
    start = 2 * pi * cumsum(dplyr::lag(prop, default = 0)),
    end   = 2 * pi * cumsum(prop)
  ) %>%
  dplyr::ungroup()

# Keep top 5 species for each celltype
top5_species <- data_pie %>%
  dplyr::group_by(celltype, species) %>%
  dplyr::summarise(total_species_count = sum(n), .groups = "drop") %>%
  dplyr::group_by(celltype) %>%
  dplyr::slice_max(order_by = total_species_count, n = 5, with_ties = FALSE) %>%
  dplyr::ungroup()

data_pie <- dplyr::inner_join(data_pie, top5_species, by = c("celltype", "species"))

# Grid positions
x_lab <- data_pie %>%
  dplyr::distinct(species) %>%
  dplyr::arrange(species) %>%
  dplyr::mutate(mid_x = dplyr::row_number())

y_lab <- data_pie %>%
  dplyr::distinct(celltype) %>%
  dplyr::arrange(celltype) %>%
  dplyr::mutate(mid_y = dplyr::row_number())

max_total <- max(data_pie$total)

data_pie <- data_pie %>%
  dplyr::left_join(x_lab, by = "species") %>%
  dplyr::left_join(y_lab, by = "celltype") %>%
  dplyr::mutate(
    r = (total / max_total)^(1 / 4) * 0.6
  )

# ============================================================
# 2. SplitDotPie plot (WT vs DB)
# ============================================================

cols <- c(
  WT = "#4DBBD5",
  DB = "#E64B35"  
)

p_splitpie <- ggplot(data_pie) +
  ggforce::geom_arc_bar(
    aes(
      x0 = mid_x,
      y0 = mid_y,
      r0 = 0,
      r = r,
      start = start,
      end = end,
      fill = phenotype,
      group = interaction(celltype, species, phenotype)
    ),
    alpha = 0.8,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(
    breaks = x_lab$mid_x,
    labels = x_lab$species,
    expand = expansion(mult = .1)
  ) +
  scale_y_reverse(
    breaks = y_lab$mid_y,
    labels = y_lab$celltype,
    expand = expansion(mult = .1)
  ) +
  coord_fixed() +
  labs(
    title = "Species Composition per Cell Type (WT vs DB)",
    x = "Species",
    y = "Cell Type",
    fill = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    legend.position = "right"
  )

print(p_splitpie)

ggsave("SplitDotPie_WTvsDB.png", p_splitpie, width = 13, height = 7, dpi = 600)
ggsave("SplitDotPie_WTvsDB.pdf", p_splitpie, width = 13, height = 7, dpi = 600)
