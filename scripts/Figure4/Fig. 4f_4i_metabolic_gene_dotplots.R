############################################################
## Figure 4f and 4i metabolic gene DotPlots

# ============================================================
# 1. Load packages
# ============================================================
library(tidyverse)
library(Seurat)
library(ggplot2)

# ============================================================
# 2. Input and output paths
# ============================================================
base_dir   <- "Figure4"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

annotated_seurat_rds <- file.path(
  "Figure3/Output/Figure3_mouse_bacterial_functional_cluster_annotated.rds"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 3. Load Seurat object
# ============================================================
seu <- readRDS(annotated_seurat_rds)


# ============================================================
# 4. Define top 20 species plus Others
# ============================================================
sp_counts <- sort(table(seu$species), decreasing = TRUE)
n_top <- min(20, length(sp_counts))
top_species <- names(sp_counts)[seq_len(n_top)]

seu$top20_species <- ifelse(
  seu$species %in% top_species,
  as.character(seu$species),
  "Others"
)

species_levels <- c(top_species, "Others")

seu$top20_species <- factor(
  seu$top20_species,
  levels = species_levels
)

species_count_summary <- tibble(
  Species = names(sp_counts),
  Cell_number = as.integer(sp_counts),
  Fig4f_group = ifelse(names(sp_counts) %in% top_species, names(sp_counts), "Others")
)

write_csv(
  species_count_summary,
  file.path(output_dir, "Fig4f_top20_species_cell_counts.csv")
)

# ============================================================
# 5. Define gut region-phenotype group order
# ============================================================
group_order <- c(
  "Cecum-WT", "Cecum-DB",
  "Colon-WT", "Colon-DB",
  "Rectum-WT", "Rectum-DB"
)

seu$group <- paste0(
  seu@meta.data[[region_col]], "-",
  seu@meta.data[[phenotype_col]]
)

seu$group <- factor(as.character(seu$group), levels = group_order)

write_csv(
  as.data.frame(table(seu$group)) %>%
    rename(Group = Var1, Cell_number = Freq),
  file.path(output_dir, "Fig4i_group_cell_counts.csv")
)

# ============================================================
# 6. Define Figure 4f metabolic genes and categories
# ============================================================
fig4f_genes <- c(
  "gdhA",
  "gltB", "gltD",
  "glnA",
  "argF",
  "argH",
  "MGYG000321881-01184",
  "MGYG000321780-01489",
  "glsA1",
  "asD",
  "fabD",
  "fabH",
  "fabG",
  "fabF",
  "fabI",
  "ompR",
  "acrA",
  "tolC",
  "cbh"
)

gene_category <- tibble(
  Gene = fig4f_genes,
  Category = c(
    rep("Ammonia\nassimilation", 6),
    rep("Organic acid\nmetabolism", 4),
    rep("Fatty acid\nsynthesis", 5),
    rep("Bile acid\nmetabolism", 4)
  ),
  Category_fill = c(
    rep("#E8E8F4", 6),
    rep("#F7D6A3", 4),
    rep("#CDEFF2", 5),
    rep("#E7DDF3", 4)
  )
)

gene_check <- tibble(
  Gene = fig4f_genes,
  Present_in_object = fig4f_genes %in% rownames(seu)
)

write_csv(
  gene_check,
  file.path(output_dir, "Fig4f_metabolic_gene_check.csv")
)

features_use <- intersect(fig4f_genes, rownames(seu))


# ============================================================
# 7. Generate Figure 4f DotPlot data
# ============================================================
p_dot_raw <- DotPlot(
  seu,
  features = features_use,
  group.by = "top20_species",
  scale = TRUE
)

dot_data <- p_dot_raw$data %>%
  mutate(
    Gene = as.character(features.plot),
    Species = as.character(id),
    Gene = factor(Gene, levels = fig4f_genes),
    Species = factor(Species, levels = species_levels),
    gene_index = match(as.character(Gene), fig4f_genes),
    gene_y = length(fig4f_genes) - gene_index + 1,
    species_x = match(as.character(Species), species_levels)
  ) %>%
  left_join(gene_category, by = c("Gene" = "Gene")) %>%
  arrange(gene_index, species_x)

write_csv(
  dot_data,
  file.path(output_dir, "Fig4f_metabolic_gene_top20_species_dotplot_source_data.csv")
)

# ============================================================
# 8. Prepare Figure 4f category bands
# ============================================================
category_band <- gene_category %>%
  mutate(
    gene_index = match(Gene, fig4f_genes),
    gene_y = length(fig4f_genes) - gene_index + 1
  ) %>%
  group_by(Category, Category_fill) %>%
  summarise(
    ymin = min(gene_y) - 0.5,
    ymax = max(gene_y) + 0.5,
    ymid = mean(range(gene_y)),
    .groups = "drop"
  ) %>%
  arrange(desc(ymid))

# ============================================================
# 9. Plot Figure 4f
# ============================================================
plot_x_min <- -2.7
plot_x_max <- length(species_levels) + 0.5

p_fig4f <- ggplot() +
  geom_rect(
    data = category_band,
    aes(xmin = plot_x_min, xmax = -0.55, ymin = ymin, ymax = ymax, fill = Category),
    colour = NA,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = dot_data,
    aes(x = species_x, y = gene_y, size = pct.exp, colour = avg.exp.scaled),
    alpha = 1
  ) +
  geom_text(
    data = category_band,
    aes(x = -2.15, y = ymid, label = Category),
    angle = 90,
    size = 5,
    lineheight = 0.9,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = setNames(category_band$Category_fill, category_band$Category),
    guide = "none"
  ) +
  scale_colour_gradient(
    low = "white",
    high = "red3",
    name = "Average Expression"
  ) +
  scale_size_continuous(
    range = c(1, 8),
    breaks = c(0, 10, 20, 30, 40, 50),
    limits = c(0, max(50, max(dot_data$pct.exp, na.rm = TRUE))),
    name = "Percent Expressed"
  ) +
  scale_x_continuous(
    breaks = seq_along(species_levels),
    labels = species_levels,
    limits = c(plot_x_min, plot_x_max),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = rev(seq_along(fig4f_genes)),
    labels = fig4f_genes,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Expression of Metabolic Genes Across the Top 20 Species",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
    axis.text.y = element_text(size = 13, face = "italic", colour = "black"),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.margin = margin(t = 8, r = 12, b = 8, l = 10, unit = "mm")
  )

ggsave(
  file.path(output_dir, "Fig4f_Expression_of_Metabolic_Genes_Across_Top20_Species.pdf"),
  p_fig4f,
  width = 10,
  height = 8
)

ggsave(
  file.path(output_dir, "Fig4f_Expression_of_Metabolic_Genes_Across_Top20_Species.png"),
  p_fig4f,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ============================================================
# 10. Figure 4i bile acid metabolism DotPlot
# ============================================================
fig4i_genes <- c("mdtK", "acrA", "ompR", "cbh")

fig4i_gene_check <- tibble(
  Gene = fig4i_genes,
  Present_in_object = fig4i_genes %in% rownames(seu)
)

write_csv(
  fig4i_gene_check,
  file.path(output_dir, "Fig4i_bile_acid_metabolism_gene_check.csv")
)

fig4i_features_use <- intersect(fig4i_genes, rownames(seu))


p_fig4i_raw <- DotPlot(
  seu,
  features = fig4i_features_use,
  group.by = "group",
  scale = TRUE
)

fig4i_data <- p_fig4i_raw$data %>%
  mutate(
    Gene = as.character(features.plot),
    Group = as.character(id),
    Gene = factor(Gene, levels = fig4i_genes),
    Group = factor(Group, levels = group_order)
  ) %>%
  arrange(match(as.character(Gene), fig4i_genes), match(as.character(Group), group_order))

write_csv(
  fig4i_data,
  file.path(output_dir, "Fig4i_bile_acid_metabolism_dotplot_source_data.csv")
)

p_fig4i <- ggplot(
  fig4i_data,
  aes(x = Group, y = Gene, size = pct.exp, colour = avg.exp.scaled)
) +
  geom_point(alpha = 1) +
  scale_x_discrete(limits = group_order, drop = FALSE) +
  scale_y_discrete(limits = rev(fig4i_genes), drop = FALSE) +
  scale_colour_gradient2(
    low = "#7B61FF",
    mid = "white",
    high = "#D7301F",
    midpoint = 0,
    name = "Average Expression"
  ) +
  scale_size_continuous(
    range = c(4, 8),
    breaks = c(2.5, 5.0, 7.5, 10.0),
    name = "Percent Expressed"
  ) +
  labs(
    title = "Bile acid metabolism",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 14, face = "italic", colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.margin = margin(t = 8, r = 12, b = 8, l = 8, unit = "mm")
  )

ggsave(
  file.path(output_dir, "Fig4i_Bile_acid_metabolism_dotplot.pdf"),
  p_fig4i,
  width = 5.5,
  height = 3.5
)

ggsave(
  file.path(output_dir, "Fig4i_Bile_acid_metabolism_dotplot.png"),
  p_fig4i,
  width = 5.5,
  height = 3.5,
  dpi = 300,
  bg = "white"
)
