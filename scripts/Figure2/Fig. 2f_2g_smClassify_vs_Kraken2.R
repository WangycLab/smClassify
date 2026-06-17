############################################################
## Figure 2f-g smClassify versus Kraken2 taxonomic comparison

############################################################

# ============================================================
# 0. Load packages
# ============================================================
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(VennDiagram)
library(grid)

set.seed(1234)

# ============================================================
# 1. Input files and output directory
# ============================================================

base_dir   <- "Figure2"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")
plot_dir   <- file.path(output_dir, "smClassify_vs_Kraken2_MGnify")

rds_file <- file.path(
  input_dir,
  "Mouse_microbial_taxonomy_smClassify_Kraken2_MGnify.rds"
)

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# 2. Read metadata and check required columns
# ============================================================

sce_bac_annot <- readRDS(rds_file)
meta <- sce_bac_annot@meta.data

required_cols <- c(
  "orig.ident",
  "smClassify_species", "MGnify_species",
  "smClassify_genus", "MGnify_genus"
)


meta$Group <- sub("-[0-9]+$", "", meta$orig.ident)

# ============================================================
# 3. Theme and thresholds
# ============================================================

theme_cns <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

threshold_df <- data.frame(
  Threshold = c("0.001%", "0.01%", "0.1%", "1%", "5%"),
  value = c(1e-5, 1e-4, 1e-3, 1e-2, 5e-2),
  stringsAsFactors = FALSE
)

overlap_cols <- c(
  "Shared" = "#56B4E9",
  "smClassify_only" = "#E69F00",
  "Kraken2_only" = "#009E73"
)

# ============================================================
# 4. Species-level abundance, overlap and top-10 summaries
# ============================================================

species_df <- meta %>%
  select(
    sm = smClassify_species,
    kr = MGnify_species
  )

species_sm_ab <- species_df %>%
  filter(!is.na(sm), sm != "") %>%
  count(sm, name = "n") %>%
  mutate(abundance = n / sum(n))

species_kr_ab <- species_df %>%
  filter(!is.na(kr), kr != "") %>%
  count(kr, name = "n") %>%
  mutate(abundance = n / sum(n))

write.csv(
  species_sm_ab,
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_smClassify_abundance.csv"),
  row.names = FALSE
)
write.csv(
  species_kr_ab,
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_Kraken2_abundance.csv"),
  row.names = FALSE
)

species_overlap_summary <- lapply(seq_len(nrow(threshold_df)), function(i) {
  th_lab <- threshold_df$Threshold[i]
  th_val <- threshold_df$value[i]

  species_sm_set <- species_sm_ab %>%
    filter(abundance > th_val) %>%
    pull(sm) %>%
    unique()

  species_kr_set <- species_kr_ab %>%
    filter(abundance > th_val) %>%
    pull(kr) %>%
    unique()

  data.frame(
    Threshold = th_lab,
    Shared = length(intersect(species_sm_set, species_kr_set)),
    smClassify_only = length(setdiff(species_sm_set, species_kr_set)),
    Kraken2_only = length(setdiff(species_kr_set, species_sm_set))
  )
}) %>%
  bind_rows()

write.csv(
  species_overlap_summary,
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_overlap_summary_by_abundance_threshold.csv"),
  row.names = FALSE
)

species_overlap_long <- species_overlap_summary %>%
  pivot_longer(
    cols = c(Shared, smClassify_only, Kraken2_only),
    names_to = "Category",
    values_to = "N"
  )

species_overlap_long$Threshold <- factor(
  species_overlap_long$Threshold,
  levels = threshold_df$Threshold
)
species_overlap_long$Category <- factor(
  species_overlap_long$Category,
  levels = c("Shared", "smClassify_only", "Kraken2_only")
)

write.csv(
  species_overlap_long,
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_overlap_summary_by_abundance_threshold_long.csv"),
  row.names = FALSE
)

p_species_overlap_threshold <- ggplot(
  species_overlap_long,
  aes(x = Threshold, y = N, fill = Category)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = N),
    position = position_dodge(width = 0.75),
    vjust = -0.25,
    size = 3.2
  ) +
  scale_fill_manual(
    values = overlap_cols,
    labels = c("Shared", "smClassify only", "Kraken2 only")
  ) +
  labs(
    title = "Species overlap across abundance thresholds",
    x = "Relative abundance threshold",
    y = "Number of Species"
  ) +
  theme_cns() +
  theme(legend.position = "top") +
  coord_cartesian(clip = "off")

ggsave(
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_overlap_by_threshold_grouped_barplot.pdf"),
  p_species_overlap_threshold,
  width = 6,
  height = 4
)
ggsave(
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_overlap_by_threshold_grouped_barplot.png"),
  p_species_overlap_threshold,
  width = 6,
  height = 4,
  dpi = 600
)

species_sm_set_all <- unique(species_sm_ab$sm)
species_kr_set_all <- unique(species_kr_ab$kr)

species_overlap_df <- data.frame(
  Category = c("Shared", "smClassify_only", "Kraken2_only"),
  Count = c(
    length(intersect(species_sm_set_all, species_kr_set_all)),
    length(setdiff(species_sm_set_all, species_kr_set_all)),
    length(setdiff(species_kr_set_all, species_sm_set_all))
  )
)

write.csv(
  species_overlap_df,
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_overlap.csv"),
  row.names = FALSE
)

p_species_overlap <- ggplot(species_overlap_df, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.6, colour = "black") +
  scale_fill_manual(values = overlap_cols) +
  labs(
    title = "Species_smClassify_vs_Kraken2_MGnify overlap",
    x = NULL,
    y = "Number of taxa"
  ) +
  theme_cns() +
  theme(legend.position = "none")

ggsave(
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_overlap.pdf"),
  p_species_overlap,
  width = 4,
  height = 3
)
ggsave(
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_overlap.png"),
  p_species_overlap,
  width = 4,
  height = 3,
  dpi = 600
)

species_top_sm <- species_sm_ab %>%
  arrange(desc(abundance)) %>%
  slice_head(n = 10) %>%
  mutate(Method = "smClassify", taxon = sm) %>%
  select(taxon, Method, abundance)

species_top_kr <- species_kr_ab %>%
  arrange(desc(abundance)) %>%
  slice_head(n = 10) %>%
  mutate(Method = "Kraken2", taxon = kr) %>%
  select(taxon, Method, abundance)

species_top_df <- bind_rows(species_top_sm, species_top_kr)

write.csv(
  species_top_df,
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_top10.csv"),
  row.names = FALSE
)

p_species_top <- ggplot(
  species_top_df,
  aes(x = reorder(taxon, abundance), y = abundance, fill = Method)
) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.3) +
  facet_wrap(~Method, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = c(
    "smClassify" = "#E69F00",
    "Kraken2" = "#009E73"
  )) +
  labs(
    title = "Species_smClassify_vs_Kraken2_MGnify Top10",
    x = NULL,
    y = "Cell fraction"
  ) +
  theme_cns() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_top10.pdf"),
  p_species_top,
  width = 6.5,
  height = 4
)
ggsave(
  file.path(plot_dir, "Species_smClassify_vs_Kraken2_MGnify_top10.png"),
  p_species_top,
  width = 6.5,
  height = 4,
  dpi = 600
)

# ============================================================
# 5. Genus-level abundance, overlap and top-10 summaries
# ============================================================

genus_df <- meta %>%
  select(
    sm = smClassify_genus,
    kr = MGnify_genus
  )

genus_sm_ab <- genus_df %>%
  filter(!is.na(sm), sm != "") %>%
  count(sm, name = "n") %>%
  mutate(abundance = n / sum(n))

genus_kr_ab <- genus_df %>%
  filter(!is.na(kr), kr != "") %>%
  count(kr, name = "n") %>%
  mutate(abundance = n / sum(n))

write.csv(
  genus_sm_ab,
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_smClassify_abundance.csv"),
  row.names = FALSE
)
write.csv(
  genus_kr_ab,
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_Kraken2_abundance.csv"),
  row.names = FALSE
)

genus_overlap_summary <- lapply(seq_len(nrow(threshold_df)), function(i) {
  th_lab <- threshold_df$Threshold[i]
  th_val <- threshold_df$value[i]

  genus_sm_set <- genus_sm_ab %>%
    filter(abundance > th_val) %>%
    pull(sm) %>%
    unique()

  genus_kr_set <- genus_kr_ab %>%
    filter(abundance > th_val) %>%
    pull(kr) %>%
    unique()

  data.frame(
    Threshold = th_lab,
    Shared = length(intersect(genus_sm_set, genus_kr_set)),
    smClassify_only = length(setdiff(genus_sm_set, genus_kr_set)),
    Kraken2_only = length(setdiff(genus_kr_set, genus_sm_set))
  )
}) %>%
  bind_rows()

write.csv(
  genus_overlap_summary,
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_overlap_summary_by_abundance_threshold.csv"),
  row.names = FALSE
)

genus_overlap_long <- genus_overlap_summary %>%
  pivot_longer(
    cols = c(Shared, smClassify_only, Kraken2_only),
    names_to = "Category",
    values_to = "N"
  )

genus_overlap_long$Threshold <- factor(
  genus_overlap_long$Threshold,
  levels = threshold_df$Threshold
)
genus_overlap_long$Category <- factor(
  genus_overlap_long$Category,
  levels = c("Shared", "smClassify_only", "Kraken2_only")
)

write.csv(
  genus_overlap_long,
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_overlap_summary_by_abundance_threshold_long.csv"),
  row.names = FALSE
)

p_genus_overlap_threshold <- ggplot(
  genus_overlap_long,
  aes(x = Threshold, y = N, fill = Category)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = N),
    position = position_dodge(width = 0.75),
    vjust = -0.25,
    size = 3.2
  ) +
  scale_fill_manual(
    values = overlap_cols,
    labels = c("Shared", "smClassify only", "Kraken2 only")
  ) +
  labs(
    title = "Genus overlap across abundance thresholds",
    x = "Relative abundance threshold",
    y = "Number of Genera"
  ) +
  theme_cns() +
  theme(legend.position = "top") +
  coord_cartesian(clip = "off")

ggsave(
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_overlap_by_threshold_grouped_barplot.pdf"),
  p_genus_overlap_threshold,
  width = 6,
  height = 4
)
ggsave(
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_overlap_by_threshold_grouped_barplot.png"),
  p_genus_overlap_threshold,
  width = 6,
  height = 4,
  dpi = 600
)

genus_sm_set_all <- unique(genus_sm_ab$sm)
genus_kr_set_all <- unique(genus_kr_ab$kr)

genus_overlap_df <- data.frame(
  Category = c("Shared", "smClassify_only", "Kraken2_only"),
  Count = c(
    length(intersect(genus_sm_set_all, genus_kr_set_all)),
    length(setdiff(genus_sm_set_all, genus_kr_set_all)),
    length(setdiff(genus_kr_set_all, genus_sm_set_all))
  )
)

write.csv(
  genus_overlap_df,
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_overlap.csv"),
  row.names = FALSE
)

p_genus_overlap <- ggplot(genus_overlap_df, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.6, colour = "black") +
  scale_fill_manual(values = overlap_cols) +
  labs(
    title = "Genus_smClassify_vs_Kraken2_MGnify overlap",
    x = NULL,
    y = "Number of taxa"
  ) +
  theme_cns() +
  theme(legend.position = "none")

ggsave(
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_overlap.pdf"),
  p_genus_overlap,
  width = 4,
  height = 3
)
ggsave(
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_overlap.png"),
  p_genus_overlap,
  width = 4,
  height = 3,
  dpi = 600
)

genus_top_sm <- genus_sm_ab %>%
  arrange(desc(abundance)) %>%
  slice_head(n = 10) %>%
  mutate(Method = "smClassify", taxon = sm) %>%
  select(taxon, Method, abundance)

genus_top_kr <- genus_kr_ab %>%
  arrange(desc(abundance)) %>%
  slice_head(n = 10) %>%
  mutate(Method = "Kraken2", taxon = kr) %>%
  select(taxon, Method, abundance)

genus_top_df <- bind_rows(genus_top_sm, genus_top_kr)

write.csv(
  genus_top_df,
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_top10.csv"),
  row.names = FALSE
)

p_genus_top <- ggplot(
  genus_top_df,
  aes(x = reorder(taxon, abundance), y = abundance, fill = Method)
) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.3) +
  facet_wrap(~Method, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = c(
    "smClassify" = "#E69F00",
    "Kraken2" = "#009E73"
  )) +
  labs(
    title = "Genus_smClassify_vs_Kraken2_MGnify Top10",
    x = NULL,
    y = "Cell fraction"
  ) +
  theme_cns() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_top10.pdf"),
  p_genus_top,
  width = 6.5,
  height = 4
)
ggsave(
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_top10.png"),
  p_genus_top,
  width = 6.5,
  height = 4,
  dpi = 600
)

# ============================================================
# 6. Figure 2f Venn diagrams
# ============================================================

genus_sm_set <- meta %>%
  filter(!is.na(smClassify_genus), smClassify_genus != "") %>%
  pull(smClassify_genus) %>%
  unique()

genus_kr_set <- meta %>%
  filter(!is.na(MGnify_genus), MGnify_genus != "") %>%
  pull(MGnify_genus) %>%
  unique()

species_sm_set <- meta %>%
  filter(!is.na(smClassify_species), smClassify_species != "") %>%
  pull(smClassify_species) %>%
  unique()

species_kr_set <- meta %>%
  filter(!is.na(MGnify_species), MGnify_species != "") %>%
  pull(MGnify_species) %>%
  unique()

genus_shared <- intersect(genus_sm_set, genus_kr_set)
genus_venn_grob <- draw.pairwise.venn(
  area1 = length(genus_sm_set),
  area2 = length(genus_kr_set),
  cross.area = length(genus_shared),
  category = c("smClassify", "Kraken2"),
  fill = c("#E69F00", "#009E73"),
  alpha = c(0.6, 0.6),
  lty = "solid",
  lwd = 1.2,
  col = c("black", "black"),
  cex = 1.4,
  cat.cex = 1.2,
  cat.col = c("black", "black"),
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  scaled = TRUE
)

pdf(
  file.path(plot_dir, "Figure2f_Genus_smClassify_vs_Kraken2_MGnify_overlap_venn.pdf"),
  width = 2,
  height = 2
)
grid.newpage()
grid.draw(genus_venn_grob)
dev.off()

png(
  file.path(plot_dir, "Figure2f_Genus_smClassify_vs_Kraken2_MGnify_overlap_venn.png"),
  width = 1500,
  height = 1500,
  res = 600
)
grid.newpage()
grid.draw(genus_venn_grob)
dev.off()

figure2f_genus_source <- data.frame(
  Category = c("Shared", "smClassify_only", "Kraken2_only"),
  Count = c(
    length(genus_shared),
    length(setdiff(genus_sm_set, genus_kr_set)),
    length(setdiff(genus_kr_set, genus_sm_set))
  )
)

write.csv(
  figure2f_genus_source,
  file.path(plot_dir, "Figure2f_Genus_smClassify_vs_Kraken2_MGnify_overlap_venn_SourceData.csv"),
  row.names = FALSE
)

species_shared <- intersect(species_sm_set, species_kr_set)
species_venn_grob <- draw.pairwise.venn(
  area1 = length(species_sm_set),
  area2 = length(species_kr_set),
  cross.area = length(species_shared),
  category = c("smClassify", "Kraken2"),
  fill = c("#E69F00", "#009E73"),
  alpha = c(0.6, 0.6),
  lty = "solid",
  lwd = 1.2,
  col = c("black", "black"),
  cex = 1.4,
  cat.cex = 1.2,
  cat.col = c("black", "black"),
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  scaled = TRUE
)

pdf(
  file.path(plot_dir, "Figure2f_Species_smClassify_vs_Kraken2_MGnify_overlap_venn.pdf"),
  width = 1.5,
  height = 1.5
)
grid.newpage()
grid.draw(species_venn_grob)
dev.off()

png(
  file.path(plot_dir, "Figure2f_Species_smClassify_vs_Kraken2_MGnify_overlap_venn.png"),
  width = 1500,
  height = 1500,
  res = 600
)
grid.newpage()
grid.draw(species_venn_grob)
dev.off()

figure2f_species_source <- data.frame(
  Category = c("Shared", "smClassify_only", "Kraken2_only"),
  Count = c(
    length(species_shared),
    length(setdiff(species_sm_set, species_kr_set)),
    length(setdiff(species_kr_set, species_sm_set))
  )
)

write.csv(
  figure2f_species_source,
  file.path(plot_dir, "Figure2f_Species_smClassify_vs_Kraken2_MGnify_overlap_venn_SourceData.csv"),
  row.names = FALSE
)

# ============================================================
# 6. Figure 2g genus-level abundance correlations
# ============================================================

get_region_from_group <- function(x) {
  dplyr::case_when(
    grepl("Cecum|MC", x, ignore.case = TRUE) ~ "Cecum",
    grepl("Colon|JC", x, ignore.case = TRUE) ~ "Colon",
    grepl("Rectum|ZC", x, ignore.case = TRUE) ~ "Rectum",
    TRUE ~ x
  )
}

get_phenotype_from_group <- function(x) {
  dplyr::case_when(
    grepl("WT", x, ignore.case = TRUE) ~ "WT",
    grepl("DB|T2DM|db/db", x, ignore.case = TRUE) ~ "DB",
    TRUE ~ NA_character_
  )
}

groups <- unique(meta$Group)
cor_list <- list()
stat_list <- list()

for (g in groups) {
  genus_group_df <- meta %>%
    filter(Group == g) %>%
    select(
      sm = smClassify_genus,
      kr = MGnify_genus
    )

  genus_group_sm_ab <- genus_group_df %>%
    filter(!is.na(sm), sm != "") %>%
    count(sm, name = "n_sm") %>%
    mutate(sm_abundance = n_sm / sum(n_sm)) %>%
    rename(taxon = sm)

  genus_group_kr_ab <- genus_group_df %>%
    filter(!is.na(kr), kr != "") %>%
    count(kr, name = "n_kr") %>%
    mutate(kr_abundance = n_kr / sum(n_kr)) %>%
    rename(taxon = kr)

  genus_group_cor <- inner_join(genus_group_sm_ab, genus_group_kr_ab, by = "taxon")

  genus_group_cor <- genus_group_cor %>%
    mutate(
      Group = g,
      Region = get_region_from_group(g),
      Phenotype = get_phenotype_from_group(g),
      log_sm = log10(sm_abundance + 1e-8),
      log_kr = log10(kr_abundance + 1e-8)
    )

  pearson <- cor(genus_group_cor$log_sm, genus_group_cor$log_kr, method = "pearson")
  spearman <- cor(genus_group_cor$log_sm, genus_group_cor$log_kr, method = "spearman")

  stat_list[[g]] <- data.frame(
    Group = g,
    Region = unique(genus_group_cor$Region),
    Phenotype = unique(genus_group_cor$Phenotype),
    Pearson = pearson,
    Spearman = spearman,
    N_taxa = nrow(genus_group_cor)
  )

  cor_list[[g]] <- genus_group_cor
}

cor_genus_all <- bind_rows(cor_list)
cor_summary_all <- bind_rows(stat_list)

write.csv(
  cor_genus_all,
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_correlation_by_group_SourceData.csv"),
  row.names = FALSE
)
write.csv(
  cor_summary_all,
  file.path(plot_dir, "Genus_smClassify_vs_Kraken2_MGnify_correlation_by_group_summary.csv"),
  row.names = FALSE
)

cor_genus_wt <- cor_genus_all %>%
  filter(Phenotype == "WT") %>%
  mutate(Region = factor(Region, levels = c("Cecum", "Colon", "Rectum")))

cor_summary_wt <- cor_summary_all %>%
  filter(Phenotype == "WT") %>%
  mutate(Region = factor(Region, levels = c("Cecum", "Colon", "Rectum"))) %>%
  arrange(Region)

write.csv(
  cor_genus_wt,
  file.path(plot_dir, "Figure2g_Genus_WT_regions_correlation_SourceData.csv"),
  row.names = FALSE
)
write.csv(
  cor_summary_wt,
  file.path(plot_dir, "Figure2g_Genus_WT_regions_correlation_summary.csv"),
  row.names = FALSE
)

lim_min <- min(c(cor_genus_wt$log_sm, cor_genus_wt$log_kr), na.rm = TRUE)
lim_max <- max(c(cor_genus_wt$log_sm, cor_genus_wt$log_kr), na.rm = TRUE)

label_df <- cor_summary_wt %>%
  mutate(
    x = lim_max,
    y = lim_min,
    label = paste0(
      "Pearson r = ", round(Pearson, 3),
      "\nSpearman rho = ", round(Spearman, 3),
      "\nN = ", N_taxa
    )
  )

p_fig2g <- ggplot(cor_genus_wt, aes(x = log_sm, y = log_kr)) +
  geom_point(
    size = 2,
    alpha = 0.8,
    colour = "black",
    fill = "#56B4E9",
    shape = 21,
    stroke = 0.3
  ) +
  geom_smooth(method = "lm", se = FALSE, colour = "#D55E00", linewidth = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 0,
    size = 3.1
  ) +
  facet_wrap(~Region, nrow = 1) +
  coord_fixed(xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max)) +
  labs(
    x = "log10(smClassify genus abundance)",
    y = "log10(Kraken2 genus abundance)"
  ) +
  theme_cns(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  file.path(plot_dir, "Figure2g_Genus_WT_regions_correlation.pdf"),
  p_fig2g,
  width = 8.5,
  height = 3.2
)
ggsave(
  file.path(plot_dir, "Figure2g_Genus_WT_regions_correlation.png"),
  p_fig2g,
  width = 8.5,
  height = 3.2,
  dpi = 600
)

for (region_i in levels(cor_genus_wt$Region)) {
  cor_region <- cor_genus_wt %>%
    filter(Region == region_i)

  stat_region <- cor_summary_wt %>%
    filter(Region == region_i)

  write.csv(
    cor_region,
    file.path(
      plot_dir,
      paste0("Figure2g_Genus_WT_", region_i, "_correlation_SourceData.csv")
    ),
    row.names = FALSE
  )

  label_region <- stat_region %>%
    mutate(
      x = lim_max,
      y = lim_min,
      label = paste0(
        "Pearson r = ", round(Pearson, 3),
        "\nSpearman rho = ", round(Spearman, 3),
        "\nN = ", N_taxa
      )
    )

  p_region <- ggplot(cor_region, aes(x = log_sm, y = log_kr)) +
    geom_point(
      size = 2,
      alpha = 0.8,
      colour = "black",
      fill = "#56B4E9",
      shape = 21,
      stroke = 0.3
    ) +
    geom_smooth(method = "lm", se = FALSE, colour = "#D55E00", linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_text(
      data = label_region,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 3.1
    ) +
    coord_fixed(xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max)) +
    labs(
      title = region_i,
      x = "log10(smClassify genus abundance)",
      y = "log10(Kraken2 genus abundance)"
    ) +
    theme_cns(base_size = 11) +
    theme(
      legend.position = "none"
    )

  ggsave(
    file.path(plot_dir, paste0("Figure2g_Genus_WT_", region_i, "_correlation.pdf")),
    p_region,
    width = 3.2,
    height = 3.2
  )
  ggsave(
    file.path(plot_dir, paste0("Figure2g_Genus_WT_", region_i, "_correlation.png")),
    p_region,
    width = 3.2,
    height = 3.2,
    dpi = 600
  )
}

