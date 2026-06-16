# ============================================================
# Figure 3f. Split pie plot of top species by functional cluster
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggforce)
  library(tibble)
  library(forcats)
  library(scales)
})

set.seed(1234)

# ------------------------------------------------------------
# 1. Input and output paths
# ------------------------------------------------------------

base_dir   <- "Figure3"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")
out_dir    <- file.path(output_dir, "Figure3f_SplitPie")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_file <- file.path(input_dir, "Figure3_mouse_bacterial_functional_cluster_annotated.rds")


seu <- readRDS(rds_file)

# ------------------------------------------------------------
# 2. Check and harmonize metadata
# ------------------------------------------------------------

functional_cluster_levels <- c(
  "0_Motile chemotactic and stress-adaptive",
  "1_Butyrate metabolism",
  "2_Polysaccharide uptake",
  "3_Organic acid metabolism",
  "4_Anaerobic energy metabolism",
  "5_Carbohydrate degradation",
  "6_Polysaccharide degradation",
  "7_Oxidative stress response"
)
functional_cluster_levels <- functional_cluster_levels[
  functional_cluster_levels %in% unique(seu$celltype)
]

if (length(functional_cluster_levels) > 0) {
  seu$celltype <- factor(seu$celltype, levels = functional_cluster_levels)
} else {
  seu$celltype <- factor(seu$celltype)
}

Idents(seu) <- seu$celltype

# ------------------------------------------------------------
# 3. Summarize cell counts by functional cluster, species and phenotype
# ------------------------------------------------------------

data_pie <- seu@meta.data %>%
  dplyr::filter(
    !is.na(celltype),
    !is.na(species),
    !is.na(phenotype)
  ) %>%
  dplyr::mutate(
    celltype = factor(celltype, levels = levels(seu$celltype)),
    phenotype = factor(phenotype, levels = c("WT", "DB"))
  ) %>%
  dplyr::group_by(celltype, species, phenotype) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
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

# ------------------------------------------------------------
# 4. Retain the top species for each functional cluster
# ------------------------------------------------------------

top_n_species <- 5

top_species <- data_pie %>%
  dplyr::group_by(celltype, species) %>%
  dplyr::summarise(total_species_count = sum(n), .groups = "drop") %>%
  dplyr::group_by(celltype) %>%
  dplyr::slice_max(
    order_by = total_species_count,
    n = top_n_species,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

data_pie <- data_pie %>%
  dplyr::inner_join(top_species, by = c("celltype", "species"))

# ------------------------------------------------------------
# 5. Define plot coordinates and pie radii
# ------------------------------------------------------------

x_lab <- data_pie %>%
  dplyr::distinct(species) %>%
  dplyr::arrange(species) %>%
  dplyr::mutate(mid_x = dplyr::row_number())

y_lab <- data_pie %>%
  dplyr::distinct(celltype) %>%
  dplyr::arrange(celltype) %>%
  dplyr::mutate(mid_y = dplyr::row_number())

max_total <- max(data_pie$total, na.rm = TRUE)

if (!is.finite(max_total) || max_total <= 0) {
  stop("No valid celltype-species counts were found for plotting.")
}

data_pie <- data_pie %>%
  dplyr::left_join(x_lab, by = "species") %>%
  dplyr::left_join(y_lab, by = "celltype") %>%
  dplyr::mutate(
    r = (total / max_total)^(1 / 4) * 0.6
  )

# ------------------------------------------------------------
# 6. Build an in-panel circle-size key
# ------------------------------------------------------------

size_key_values <- c(10, 100, 1000, 10000)
size_key_values <- size_key_values[size_key_values <= max_total]
if (!length(size_key_values)) {
  size_key_values <- signif(max_total, digits = 1)
}

legend_x <- max(x_lab$mid_x) + 1.4
legend_y_start <- min(y_lab$mid_y) + 0.5

size_key <- tibble::tibble(
  total = size_key_values,
  r = (total / max_total)^(1 / 4) * 0.6,
  x0 = legend_x,
  y0 = legend_y_start + seq_along(size_key_values) * 0.9,
  label = scales::comma(total)
)

# ------------------------------------------------------------
# 7. Plot settings
# ------------------------------------------------------------

group_cols <- c(
  WT   = "#4DBBD5",
  DB = "#E64B35"
)

x_limits <- c(0.5, max(x_lab$mid_x) + 2.8)

# ------------------------------------------------------------
# 8. Draw split pie plot
# ------------------------------------------------------------

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
    alpha = 1,
    color = "black",
    linewidth = 0.2
  ) +
  ggforce::geom_circle(
    data = size_key,
    aes(x0 = x0, y0 = y0, r = r),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    data = size_key,
    aes(
      x = x0 + 0.75,
      y = y0,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.2,
    color = "black"
  ) +
  annotate(
    "text",
    x = legend_x,
    y = legend_y_start,
    label = "Cell count",
    hjust = 0.5,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = group_cols,
    name = "Group",
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = x_lab$mid_x,
    labels = x_lab$species,
    limits = x_limits,
    expand = expansion(mult = 0.02)
  ) +
  scale_y_reverse(
    breaks = y_lab$mid_y,
    labels = y_lab$celltype,
    expand = expansion(mult = 0.10)
  ) +
  coord_fixed(clip = "off") +
  labs(
    title = "Species composition per functional cluster",
    x = "Species",
    y = "Functional cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black",
      size = 12
    ),
    axis.text.y = element_text(
      color = "black",
      size = 12
    ),
    axis.title = element_text(
      color = "black",
      size = 13
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey85"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 8, r = 35, b = 8, l = 8, unit = "mm")
  )

print(p_splitpie)

# ------------------------------------------------------------
# 9. Export figure and source data
# ------------------------------------------------------------

ggsave(
  filename = file.path(out_dir, "Figure3f_split_pie_top_species_WT_vs_DB.png"),
  plot = p_splitpie,
  width = 13,
  height = 7,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Figure3f_split_pie_top_species_WT_vs_DB.pdf"),
  plot = p_splitpie,
  width = 13,
  height = 7,
  dpi = 600
)

write.csv(
  data_pie,
  file.path(out_dir, "Figure3f_split_pie_source_data.csv"),
  row.names = FALSE
)

write.csv(
  top_species,
  file.path(out_dir, "Figure3f_top_species_by_functional_cluster.csv"),
  row.names = FALSE
)

