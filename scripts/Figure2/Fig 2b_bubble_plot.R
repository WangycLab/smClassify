# Fig2b annotation summary


# 1. Load packages ---------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(cowplot)


# 2. Load annotation data --------------------------------------------------

eggnog <- read.delim2("~/Mouse_gut/M1_eggNOG.tsv", header = TRUE, check.names = FALSE)
gtfinfo <- read.csv("~/Mouse_gut/M1_genome_all.fix_info.csv", header = TRUE, stringsAsFactors = FALSE)

annotation_cols <- c("COG" = "COG_category", "KEGG" = "KEGG_Pathway", "GO" = "GOs", "CAZy" = "CAZy")
annotation_df <- eggnog[, annotation_cols]
colnames(annotation_df) <- names(annotation_cols)

annotation_df[annotation_df == "-"] <- NA


# 3. Pie charts: annotated gene proportion --------------------------------

mgnify_anno_num <- sum(apply(annotation_df, 1, function(row) any(!is.na(row))))

mgnify_pie_data <- data.frame(
  Database = "m-MGnify",
  category = c("Annotated", "Unannotated"),
  count = c(mgnify_anno_num, nrow(gtfinfo) - mgnify_anno_num)
)
mgnify_pie_data$percentage <- round(mgnify_pie_data$count / sum(mgnify_pie_data$count) * 100, 1)

refseq_pie_data <- data.frame(
  Database = "p-RefSeq",
  category = c("Annotated", "Unannotated"),
  percentage = c(39.1, 60.9)
)

p_refseq_pie <- ggplot(refseq_pie_data, aes(x = "", y = percentage, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("Annotated" = "#1F78B4", "Unannotated" = "lightgray")) +
  geom_text(
    aes(label = paste0(percentage, "%")),
    color = "black",
    size = 2.5,
    position = position_stack(vjust = 0.5)
  ) +
  labs(title = "p-RefSeq") +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10)
  )

p_mgnify_pie <- ggplot(mgnify_pie_data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("Annotated" = "#1F78B4", "Unannotated" = "lightgray")) +
  geom_text(
    aes(label = paste0(percentage, "%")),
    color = "black",
    size = 2.5,
    position = position_stack(vjust = 0.5)
  ) +
  labs(title = "m-MGnify") +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10)
  )

p_pie_legend <- ggplot(mgnify_pie_data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("Annotated" = "#1F78B4", "Unannotated" = "lightgray")) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  )

pie_legend <- get_legend(p_pie_legend)

p_pies <- plot_grid(
  pie_legend,
  plot_grid(p_refseq_pie, p_mgnify_pie, nrow = 1),
  ncol = 1,
  rel_heights = c(0.2, 1)
)

p_pies


# 4. Bubble plot: annotation database combinations -------------------------

annotation_logic <- annotation_df %>%
  mutate(across(everything(), ~ !is.na(.)))

database_counts <- colSums(annotation_logic)

combinations <- annotation_logic %>%
  summarise(
    COG_KEGG = sum(COG & KEGG, na.rm = TRUE),
    COG_GO = sum(COG & GO, na.rm = TRUE),
    COG_CAZy = sum(COG & CAZy, na.rm = TRUE),
    KEGG_GO = sum(KEGG & GO, na.rm = TRUE),
    KEGG_CAZy = sum(KEGG & CAZy, na.rm = TRUE),
    GO_CAZy = sum(GO & CAZy, na.rm = TRUE),
    COG_KEGG_GO = sum(COG & KEGG & GO, na.rm = TRUE),
    COG_KEGG_CAZy = sum(COG & KEGG & CAZy, na.rm = TRUE),
    COG_GO_CAZy = sum(COG & GO & CAZy, na.rm = TRUE),
    KEGG_GO_CAZy = sum(KEGG & GO & CAZy, na.rm = TRUE),
    COG_KEGG_GO_CAZy = sum(COG & KEGG & GO & CAZy, na.rm = TRUE)
  )

bubble_data <- data.frame(
  COG_KEGG = c(TRUE, TRUE, FALSE, FALSE),
  COG_GO = c(TRUE, FALSE, TRUE, FALSE),
  COG_CAZy = c(TRUE, FALSE, FALSE, TRUE),
  KEGG_GO = c(FALSE, TRUE, TRUE, FALSE),
  KEGG_CAZy = c(FALSE, TRUE, FALSE, TRUE),
  GO_CAZy = c(FALSE, FALSE, TRUE, TRUE),
  COG_KEGG_GO = c(TRUE, TRUE, TRUE, FALSE),
  COG_KEGG_CAZy = c(TRUE, TRUE, FALSE, TRUE),
  COG_GO_CAZy = c(TRUE, FALSE, TRUE, TRUE),
  KEGG_GO_CAZy = c(FALSE, TRUE, TRUE, TRUE),
  COG_KEGG_GO_CAZy = c(TRUE, TRUE, TRUE, TRUE)
)
bubble_data$database <- names(annotation_cols)

bubble_long <- melt(bubble_data, id = "database")
bubble_long$database <- factor(bubble_long$database, levels = rev(bubble_data$database))
bubble_long$database_label <- paste0(
  bubble_long$database,
  " ",
  database_counts[as.character(bubble_long$database)]
)

database_label_levels <- paste0(
  rev(names(annotation_cols)),
  " ",
  database_counts[rev(names(annotation_cols))]
)
bubble_long$database_label <- factor(bubble_long$database_label, levels = database_label_levels)

combination_counts <- t(combinations) %>% as.data.frame()
combination_counts$combination <- rownames(combination_counts)
combination_counts$log10_count <- log10(combination_counts$V1)
combination_counts$combination <- factor(combination_counts$combination, levels = combination_counts$combination)

p_combination_bar <- ggplot(combination_counts, aes(x = combination, y = log10_count)) +
  geom_bar(stat = "identity", fill = "gray35", width = 0.8) +
  labs(title = "m-MGnify", x = "", y = "log10 Count") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 6, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 8),
    legend.position = "none"
  )

p_bubble_matrix <- ggplot(bubble_long, aes(x = variable, y = database_label, size = value, color = value)) +
  geom_point() +
  scale_size_manual(values = c(`FALSE` = 0, `TRUE` = 2.2)) +
  scale_color_manual(values = c(`FALSE` = "white", `TRUE` = "red")) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 5, color = "black"),
    legend.position = "none"
  )

p_bubble <- plot_grid(
  p_combination_bar,
  p_bubble_matrix,
  ncol = 1,
  align = "v",
  rel_heights = c(2, 1)
)

p_bubble


# 5. Fig2b combined panel --------------------------------------------------

fig2b <- plot_grid(
  p_pies,
  p_bubble,
  nrow = 1,
  rel_widths = c(1.2, 1)
)

fig2b
