# Supplementary Figure 1: m-MGnify KEGG annotation


# 1. Load packages ---------------------------------------------------------

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)


# 2. Load common data ------------------------------------------------------

kegg_level <- read.delim2("~/Mouse_gut/KEGG_levels_output.tsv", header = TRUE)
eggnog <- read.delim2("~/Mouse_gut/M1_eggNOG.tsv", header = TRUE, check.names = FALSE)


# 3. Supplementary Fig. 1a: KEGG pathway circular bar plot -----------------

sample_q <- readRDS("Expression_martrix.rds")

genes_in_sample <- rownames(sample_q)
genes_in_sample <- gsub("-", "_", genes_in_sample)

eggnog_in_sample <- eggnog[eggnog$`#query` %in% genes_in_sample, ]

KEGG_table <- eggnog_in_sample %>%
  dplyr::select(KEGG_Pathway, `#query`) %>%
  filter(KEGG_Pathway != "-") %>%
  mutate(KEGG_Pathway = strsplit(KEGG_Pathway, ",")) %>%
  unnest(KEGG_Pathway)

KEGG_table <- KEGG_table[grep("map", KEGG_table$KEGG_Pathway), ]

pathway_count <- table(KEGG_table$KEGG_Pathway) %>% as.data.frame()
pathway_count <- merge(pathway_count, kegg_level, by.x = "Var1", by.y = "kegg")
pathway_count <- pathway_count %>%
  group_by(Level1) %>%
  slice_max(n = 20, order_by = Freq) %>%
  ungroup()

write.csv(pathway_count, file = "pathway_gene_freq.csv", row.names = FALSE)

coord_data <- pathway_count
colnames(coord_data) <- c("KEGG", "count", "Level1", "Level2", "Level3")
coord_data$id <- 1:nrow(coord_data)

number <- nrow(coord_data)
angle <- 90 - 360 * (coord_data$id - 0.5) / number
coord_data$hjust <- ifelse(angle < -90, 1, 0)
coord_data$angle <- ifelse(angle < -90, angle + 180, angle)

coord_theme <- theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.25, "cm")
  )

p_kegg_circle <- ggplot(coord_data, aes(x = as.factor(id), y = count, fill = Level1)) +
  geom_bar(stat = "identity") +
  ylim(-500, 1500) +
  coord_polar(start = 0) +
  coord_theme +
  scale_fill_brewer(palette = "Set3") +
  geom_text(
    data = coord_data,
    aes(
      x = id,
      y = count + 5,
      label = paste0(Level3, "(", count, ")"),
      hjust = hjust
    ),
    color = "black",
    size = 2.6,
    angle = coord_data$angle,
    inherit.aes = FALSE
  )

p_kegg_circle


# 4. Supplementary Fig. 1b: KEGG pathway heatmap by species ----------------

mgnify_info <- read.delim2("~/Mouse_gut/genomes-all_metadata.tsv", header = TRUE)
mgnify_info <- mgnify_info %>%
  separate(Lineage, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(across(everything(), ~ sub("^[a-z]__", "", .)))

species_info <- mgnify_info %>% distinct(Genome, p, s, .keep_all = FALSE)

eggnog_species <- eggnog
eggnog_species$species <- substr(eggnog_species$`#query`, 1, 13)
eggnog_species <- merge(eggnog_species, species_info, by.x = "species", by.y = "Genome")

eggnog_freq <- table(eggnog_species$p) %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  head(8)

top_phylum_species <- species_info[species_info$p %in% eggnog_freq$Var1, ]

eggnog_species <- eggnog_species[eggnog_species$species %in% top_phylum_species$Genome, ]
eggnog_species$`#query` <- gsub("_", "-", eggnog_species$`#query`)

KEGG_species_table <- eggnog_species %>%
  dplyr::select(KEGG_Pathway, `#query`, species) %>%
  filter(KEGG_Pathway != "-") %>%
  mutate(KEGG_Pathway = strsplit(KEGG_Pathway, ",")) %>%
  unnest(KEGG_Pathway)

KEGG_species_table <- KEGG_species_table[grepl("map", KEGG_species_table$KEGG_Pathway), ]
KEGG_species_table <- merge(KEGG_species_table, species_info, by.x = "species", by.y = "Genome")
KEGG_species_table <- merge(KEGG_species_table, kegg_level, by.x = "KEGG_Pathway", by.y = "kegg")
KEGG_species_table$s <- ifelse(KEGG_species_table$s == "", KEGG_species_table$species, KEGG_species_table$s)

interested_pathway <- c(
  "Metabolic pathways",
  "Biosynthesis of amino acids",
  "Carbon metabolism",
  "Pyruvate metabolism",
  "Starch and sucrose metabolism",
  "Galactose metabolism",
  "Butanoate metabolism",
  "Propanoate metabolism",
  "Fatty acid biosynthesis",
  "Valine, leucine and isoleucine biosynthesis",
  "Fatty acid metabolism"
)

species_description_count <- KEGG_species_table %>%
  count(KEGG_Pathway, Level3, s, p)

species_description_count <- species_description_count[
  species_description_count$Level3 %in% interested_pathway,
]
species_description_count$p <- factor(species_description_count$p, levels = eggnog_freq$Var1)
species_description_count <- arrange(species_description_count, p)

data_matrix <- species_description_count %>%
  select(Level3, s, n) %>%
  pivot_wider(names_from = s, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Level3") %>%
  as.matrix()

log_matrix <- log10(data_matrix + 1)

annotation_df <- species_description_count %>% distinct(s, p)
annotation_df <- column_to_rownames(annotation_df, var = "s")

annotation_colors <- structure(
  brewer.pal(length(unique(annotation_df$p)), "Set3"),
  names = as.character(unique(annotation_df$p))
)

col_anno <- columnAnnotation(
  Phylum = annotation_df[colnames(data_matrix), "p"],
  col = list(Phylum = annotation_colors),
  annotation_legend_param = list(title = "Phylum")
)

ht_kegg <- Heatmap(
  log_matrix,
  name = "Log10 count",
  col = colorRamp2(c(1, max(log_matrix)), c("white", "red")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  use_raster = FALSE,
  show_column_names = FALSE,
  top_annotation = col_anno,
  row_names_gp = gpar(fontsize = 10)
)

draw(ht_kegg)


# 5. Combined Supplementary Figure 1 --------------------------------------

p_heatmap <- ggdraw() + draw_grob(grid.grabExpr(draw(ht_kegg)))

figS1 <- plot_grid(
  p_kegg_circle,
  p_heatmap,
  ncol = 1,
  labels = c("a", "b"),
  label_size = 16,
  label_fontface = "bold",
  rel_heights = c(2.3, 1)
)

figS1
