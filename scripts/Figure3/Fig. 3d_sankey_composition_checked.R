# ============================================================
# Figure 3d | Sankey plot and functional-cluster composition
#

# ============================================================


library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(scales)
library(stringr)
library(forcats)
library(networkD3)
library(htmlwidgets)
library(webshot2)
library(patchwork)
library(jsonlite)


set.seed(1234)
options(stringsAsFactors = FALSE)

# ============================================================
# 1. Paths and input files
# ============================================================

base_dir   <- "Figure3"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")
out_dir    <- file.path(output_dir, "Figure3d_Sankey_Composition")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Preferred input from Figure 3c.
rds_file <- file.path(input_dir, "Figure3_mouse_bacterial_functional_cluster_annotated.rds")

seu<-readRDS(rds_file )
# ============================================================
# 3. Prepare metadata
# ============================================================

meta <- seu@meta.data

seu$group <- sub("-[0-9]$", "", as.character(meta$orig.ident))

seu$position <- dplyr::case_when(
  grepl("^Cecum",  as.character(meta$orig.ident)) ~ "Cecum",
  grepl("^Colon",  as.character(meta$orig.ident)) ~ "Colon",
  grepl("^Rectum", as.character(meta$orig.ident)) ~ "Rectum",
  TRUE ~ NA_character_
)

seu$region <- seu$position
seu$phenotype <- ifelse(grepl("WT", as.character(meta$orig.ident)), "WT", "DB")
seu$replicate <- sub("^.*-([0-9]+)$", "\\1", as.character(meta$orig.ident))

# Fixed factor levels for consistent plotting.
sample_levels <- c(
  "Cecum-WT-1", "Cecum-WT-2", "Cecum-DB-1", "Cecum-DB-2",
  "Colon-WT-1", "Colon-WT-2", "Colon-DB-1", "Colon-DB-2",
  "Rectum-WT-1", "Rectum-WT-2", "Rectum-DB-1", "Rectum-DB-2"
)
sample_levels <- sample_levels[sample_levels %in% unique(as.character(seu$orig.ident))]

region_levels <- c("Cecum", "Colon", "Rectum")
phenotype_levels <- c("WT", "DB")
group_levels <- c(
  "Cecum-WT", "Cecum-DB",
  "Colon-WT", "Colon-DB",
  "Rectum-WT", "Rectum-DB"
)
group_levels <- group_levels[group_levels %in% unique(as.character(seu$group))]

seu$orig.ident <- factor(as.character(seu$orig.ident), levels = sample_levels)
seu$phenotype  <- factor(seu$phenotype, levels = phenotype_levels)
seu$position   <- factor(seu$position, levels = region_levels)
seu$region     <- factor(seu$region, levels = region_levels)
seu$replicate  <- factor(seu$replicate, levels = c("1", "2"))
seu$group      <- factor(seu$group, levels = group_levels)


celltype_levels_expected <- c(
  "0_Motile chemotactic and stress-adaptive",
  "1_Butyrate metabolism",
  "2_Polysaccharide uptake",
  "3_Organic acid metabolism",
  "4_Anaerobic energy metabolism",
  "5_Carbohydrate degradation",
  "6_Polysaccharide degradation",
  "7_Oxidative stress response"
)
celltype_levels <- celltype_levels_expected[celltype_levels_expected %in% unique(seu$celltype)]


extra_celltypes <- setdiff(sort(unique(seu$celltype)), celltype_levels)
celltype_levels <- c(celltype_levels, extra_celltypes)

seu$celltype <- factor(seu$celltype, levels = celltype_levels)
Idents(seu) <- seu$celltype


meta_df <- as.data.frame(seu@meta.data) %>%
  tibble::rownames_to_column("cell")

# ============================================================
# 4. Colors
# ============================================================

pal_pheno <- c(
  WT = "#3498DB",
  DB = "#D94F45"
)

pal_region <- c(
  Cecum  = "#4C78A8",
  Colon  = "#59A14F",
  Rectum = "#F28E2B"
)


cluster_cols_nature <- c(
  "0" = "#0071B1",
  "1" = "#AF79A0",
  "2" = "#E69F00",
  "3" = "#4E79A7",
  "4" = "#ECC847",
  "5" = "#58A04F",
  "6" = "#999999",
  "7" = "#D55E00"
)

cellnum_col <- "#BCBD22"


celltype_cluster_id <- sub("^([0-9]+).*", "\\1", celltype_levels)
names(celltype_cluster_id) <- celltype_levels


pal_celltype <- setNames(
  cluster_cols_nature[celltype_cluster_id],
  celltype_levels
)

pal_group <- c(
  "Cecum-WT"  = unname(pal_pheno["WT"]),
  "Cecum-DB"  = unname(pal_pheno["DB"]),
  "Colon-WT"  = unname(pal_pheno["WT"]),
  "Colon-DB"  = unname(pal_pheno["DB"]),
  "Rectum-WT" = unname(pal_pheno["WT"]),
  "Rectum-DB" = unname(pal_pheno["DB"])
)
pal_group <- pal_group[names(pal_group) %in% group_levels]

# ============================================================
# 5. Plot theme and save helper
# ============================================================

theme_fig3d_bar <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text.x = element_text(size = base_size * 0.85, colour = "black"),
      axis.text.y = element_text(size = base_size * 0.85, colour = "black"),
      legend.title = element_text(size = base_size * 0.9, face = "bold"),
      legend.text = element_text(size = base_size * 0.85),
      legend.key.size = unit(0.45, "cm"),
      panel.grid = element_blank(),
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5)
    )
}

save_plot <- function(p, filename, width, height, dpi = 600) {
  ggsave(
    filename = file.path(out_dir, paste0(filename, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = dpi
  )
  ggsave(
    filename = file.path(out_dir, paste0(filename, ".pdf")),
    plot = p,
    width = width,
    height = height
  )
}

# ============================================================
# 6. Sankey plot: group -> functional cluster
# ============================================================

sankey_counts <- as.data.frame(
  table(Group = seu$group, Celltype = seu$celltype)
)
colnames(sankey_counts) <- c("Group", "Celltype", "n")

sankey_frac <- as.data.frame(
  prop.table(
    table(Group = seu$group, Celltype = seu$celltype),
    margin = 2
  )
)
colnames(sankey_frac) <- c("Group", "Celltype", "Freq")

sankey_df <- sankey_frac %>%
  dplyr::left_join(sankey_counts, by = c("Group", "Celltype")) %>%
  dplyr::filter(Freq > 0) %>%
  dplyr::mutate(
    Group = as.character(Group),
    Celltype = as.character(Celltype),
    Phenotype = ifelse(grepl("-WT$", Group), "WT", "DB")
  )

nodes <- data.frame(
  name = c(group_levels, celltype_levels),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(name %in% c(sankey_df$Group, sankey_df$Celltype)) %>%
  dplyr::mutate(
    node_type = dplyr::case_when(
      name %in% names(pal_group) ~ "group",
      name %in% names(pal_celltype) ~ "celltype",
      TRUE ~ "other"
    ),
    cluster_id = dplyr::case_when(
      node_type == "celltype" ~ sub("^([0-9]+).*", "\\1", name),
      TRUE ~ NA_character_
    ),
    node_color = dplyr::case_when(
      name %in% names(pal_group) ~ pal_group[name],
      name %in% names(pal_celltype) ~ pal_celltype[name],
      TRUE ~ "#BDBDBD"
    ),
    node_group = name
  )

links <- sankey_df %>%
  dplyr::mutate(
    IDsource = match(Group, nodes$name) - 1,
    IDtarget = match(Celltype, nodes$name) - 1,
    link_group = Phenotype,
    link_color = dplyr::case_when(
      Phenotype == "WT" ~ unname(pal_pheno["WT"]),
      Phenotype == "DB" ~ unname(pal_pheno["DB"]),
      TRUE ~ "#BDBDBD"
    )
  ) %>%
  dplyr::select(IDsource, IDtarget, Freq, n, link_group, link_color)


node_color_map <- c(pal_group, pal_celltype)
node_color_map <- node_color_map[!is.na(names(node_color_map))]
node_color_map <- node_color_map[!is.na(node_color_map)]

link_color_map <- c(
  WT = unname(pal_pheno["WT"]),
  DB = unname(pal_pheno["DB"])
)

node_color_json <- jsonlite::toJSON(as.list(node_color_map), auto_unbox = TRUE)
link_color_json <- jsonlite::toJSON(as.list(link_color_map), auto_unbox = TRUE)

colourScale <- sprintf(
  'd3.scaleOrdinal().domain(["WT", "DB"]).range(["%s", "%s"])',
  unname(pal_pheno["WT"]),
  unname(pal_pheno["DB"])
)

p_sankey <- networkD3::sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "IDsource",
  Target = "IDtarget",
  Value = "Freq",
  NodeID = "name",
  NodeGroup = "node_group",
  LinkGroup = "link_group",
  colourScale = htmlwidgets::JS(colourScale),
  nodeWidth = 20,
  nodePadding = 12,
  fontSize = 24,
  sinksRight = FALSE
)

# Force exact node and link colors by visible node name.
p_sankey <- htmlwidgets::onRender(
  p_sankey,
  sprintf(
    '
    function(el, x) {
      var nodeColor = %s;
      var linkColor = %s;

      d3.select(el).selectAll(".node rect")
        .style("fill", function(d) {
          return nodeColor[d.name] || "#BDBDBD";
        })
        .style("stroke", function(d) {
          return d3.rgb(nodeColor[d.name] || "#BDBDBD").darker(0.6);
        });

      d3.select(el).selectAll(".link")
        .style("stroke", function(d) {
          var pheno = d.source.name.endsWith("-WT") ? "WT" : "DB";
          return linkColor[pheno] || "#BDBDBD";
        })
        .style("stroke-opacity", 0.55);
    }
    ',
    node_color_json,
    link_color_json
  )
)

html_file <- file.path(out_dir, "Figure3d_sankey_group_to_functional_cluster.html")
png_file  <- file.path(out_dir, "Figure3d_sankey_group_to_functional_cluster.png")

htmlwidgets::saveWidget(
  p_sankey,
  file = html_file,
  selfcontained = TRUE
)
p_sankey

webshot2::webshot(
  url = html_file,
  file = png_file,
  vwidth = 1800,
  vheight = 1200,
  zoom = 2
)
# ============================================================
# 7. Composition plot 1: phenotype fraction within each cluster
# ============================================================

phenotype_frac_df <- meta_df %>%
  dplyr::filter(!is.na(celltype), !is.na(phenotype)) %>%
  dplyr::count(celltype, phenotype, name = "n") %>%
  group_by(celltype) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

p_phenotype <- ggplot(
  phenotype_frac_df,
  aes(
    x = freq,
    y = fct_rev(factor(celltype, levels = celltype_levels)),
    fill = phenotype
  )
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = pal_pheno) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Fraction", y = NULL, fill = "Phenotype") +
  theme_fig3d_bar(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )
p_phenotype
# ============================================================
# 8. Composition plot 2: region fraction within each cluster
# ============================================================

region_frac_df <- meta_df %>%
  dplyr::filter(!is.na(celltype), !is.na(region)) %>%
  dplyr::count(celltype, region, name = "n") %>%
  group_by(celltype) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

p_region <- ggplot(
  region_frac_df,
  aes(
    x = freq,
    y = fct_rev(factor(celltype, levels = celltype_levels)),
    fill = region
  )
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = pal_region) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Fraction", y = NULL, fill = "Region") +
  theme_fig3d_bar(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

# ============================================================
# 9. Composition plot 3: total cell number per cluster
# ============================================================

cellnum_df <- meta_df %>%
  dplyr::filter(!is.na(celltype)) %>%
  dplyr::count(celltype, name = "n")

p_cellnum <- ggplot(
  cellnum_df,
  aes(
    x = n / 1000,
    y = fct_rev(factor(celltype, levels = celltype_levels))
  )
) +
  geom_col(fill = cellnum_col, width = 0.75, color = "black", linewidth = 0.15) +
  labs(x = "Cell number (脳10鲁)", y = NULL) +
  theme_fig3d_bar(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# ============================================================
# 10. Combine and save composition plots
# ============================================================

p_composition <- p_phenotype + p_region + p_cellnum +
  patchwork::plot_layout(widths = c(1.15, 1.15, 0.8))

save_plot(
  p_composition,
  filename = "Figure3d_functional_cluster_composition_phenotype_region_cellnumber",
  width = 8,
  height = 2.5,
  dpi = 600
)

save_plot(
  p_phenotype,
  filename = "Figure3d_composition_phenotype",
  width = 3.2,
  height = 2.5,
  dpi = 600
)

save_plot(
  p_region,
  filename = "Figure3d_composition_region",
  width = 3.2,
  height = 2.5,
  dpi = 600
)

save_plot(
  p_cellnum,
  filename = "Figure3d_composition_cell_number",
  width = 2.5,
  height = 2.5,
  dpi = 600
)

# ============================================================
# 11. Export source data and color-check tables
# ============================================================

write.csv(
  sankey_df,
  file.path(out_dir, "Figure3d_sankey_group_to_functional_cluster_source_data.csv"),
  row.names = FALSE
)

write.csv(
  links,
  file.path(out_dir, "Figure3d_sankey_links_source_data.csv"),
  row.names = FALSE
)

write.csv(
  nodes,
  file.path(out_dir, "Figure3d_sankey_nodes_source_data.csv"),
  row.names = FALSE
)

write.csv(
  phenotype_frac_df,
  file.path(out_dir, "Figure3d_functional_cluster_phenotype_fraction_source_data.csv"),
  row.names = FALSE
)

write.csv(
  region_frac_df,
  file.path(out_dir, "Figure3d_functional_cluster_region_fraction_source_data.csv"),
  row.names = FALSE
)

write.csv(
  cellnum_df,
  file.path(out_dir, "Figure3d_functional_cluster_total_cell_number_source_data.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(name = names(node_color_map), color = unname(node_color_map)),
  file.path(out_dir, "Figure3d_sankey_node_color_check.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(name = names(link_color_map), color = unname(link_color_map)),
  file.path(out_dir, "Figure3d_sankey_link_color_check.csv"),
  row.names = FALSE
)

saveRDS(
  seu,
  file.path(out_dir, "Figure3d_functional_cluster_composition_input_annotated.rds")
)

