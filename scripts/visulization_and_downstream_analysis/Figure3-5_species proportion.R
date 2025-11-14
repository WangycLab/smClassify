library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggforce)
library(tibble)
library(forcats)

setwd("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009")
load("seu_symbol_SCT_clustered_annotated.RData")

stopifnot(
  all(c("celltype", "species", "phenotype") %in% colnames(seu@meta.data))
)

Idents(seu) <- seu$celltype

## Calculate angle and radius
data_pie <- seu@meta.data %>%
  group_by(celltype, species, phenotype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype, species) %>%
  mutate(total = sum(n),
         prop  = n / total) %>%
  arrange(celltype, species, phenotype) %>%
  mutate(start = 2 * pi * cumsum(lag(prop, default = 0)),
         end   = 2 * pi * cumsum(prop)) %>%
  ungroup()

## Only keep top 3 species of each celltype
top3_species <- data_pie %>%
  group_by(celltype, species) %>%
  summarise(total_species_count = sum(n), .groups = "drop") %>%
  group_by(celltype) %>%
  slice_max(order_by = total_species_count, n = 5, with_ties = FALSE) %>%
  ungroup()

data_pie <- inner_join(data_pie, top3_species, by = c("celltype", "species"))

## Set mid point of coordinates 
x_lab <- data_pie %>% distinct(species)  %>% arrange(species)  %>% mutate(mid_x = row_number())
y_lab <- data_pie %>% distinct(celltype) %>% arrange(celltype) %>% mutate(mid_y = row_number())

max_total <- max(data_pie$total)

data_pie <- data_pie %>%
  left_join(x_lab, by = "species") %>%
  left_join(y_lab, by = "celltype") %>%
  mutate(r = sqrt(total / max_total))

data_pie <- data_pie %>%
  mutate(r = (total / max_total)^(1/4) * 0.6)


cols <- c(WT = "#4DBBD5", T2DM = "#E64B35")

p_splitpie <- ggplot(data_pie) +
  geom_arc_bar(
    aes(x0 = mid_x, y0 = mid_y, r0 = 0, r = r,
        start = start, end = end,
        fill = phenotype,
        group = interaction(celltype, species, phenotype)),
    alpha = 0.8, color = "black", linewidth = 0.2
  ) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = x_lab$mid_x, labels = x_lab$species, expand = expansion(mult = .1)) +
  scale_y_reverse(breaks = y_lab$mid_y, labels = y_lab$celltype, expand = expansion(mult = .1)) +
  coord_fixed() +
  labs(title = "Species Composition per Cell Type (WT vs T2DM)",
       x = "Species", y = "Cell Type", fill = "Group") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "right")


print(p_splitpie)
ggsave("SplitDotPie_WTvsT2DM.png", p_splitpie, width = 13, height = 7, dpi = 600)
ggsave("SplitDotPie_WTvsT2DM.pdf", p_splitpie, width = 13, height = 7, dpi = 600)

