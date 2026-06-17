## =========================================================
## Supplementary Figure 7a
## Metagenome Bracken composition barplot (SCI style cleaned)
## =========================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)

## =========================================================
## 0. paths
## =========================================================

base_dir   <- "Figure2"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

bracken_long_file <- file.path(input_dir, "bracken_long.csv")

plot_dir <- output_dir

sample_order_all <- c(
  "WT1","WT2","WT3","WT4","WT5","WT6",
  "DB1","DB2","DB3","DB4","DB5","DB6"
)

topN_map <- list(
  Family = 10,
  Genus = 15,
  Species = 20,
  Species_harmonized = 20
)

## =========================================================
## 1. helper
## =========================================================

`%||%` <- function(x, y) if (!is.null(x)) x else y

make_named_palette <- function(df, base_cols = NULL) {

  taxa <- unique(df$taxon)
  taxa <- taxa[!is.na(taxa)]

  taxa_main <- if ("Others" %in% taxa) setdiff(taxa, "Others") else taxa
  n <- length(taxa_main)

  if (is.null(base_cols)) {
    cols <- if (n <= 8) {
      brewer.pal(max(3, n), "Set2")[seq_len(n)]
    } else if (n <= 12) {
      brewer.pal(12, "Set3")[seq_len(n)]
    } else {
      colorRampPalette(c(
        "#1F77B4","#2CA02C","#FF7F0E","#6A5ACD","#8C564B",
        "#E377C2","#7F7F7F","#BCBD22","#17BECF"
      ))(n)
    }
  } else {
    cols <- colorRampPalette(base_cols)(n)
  }

  pal <- setNames(cols, taxa_main)

  if ("Others" %in% taxa) {
    pal <- c(pal, Others = "grey80")
  }

  pal
}

## =========================================================
## 2. load data + SCI-style sample harmonization
## =========================================================

bracken_long <- read.csv(bracken_long_file)

mgx_all <- bracken_long %>%
  mutate(
    sample = as.character(sample),

    ## SCI-style inline mapping (NO function)
    sample_plot = sample %>%
      str_remove("^MC_") %>%
      str_replace("^WT8([1-6])$", "WT\\1") %>%
      str_replace("^DB8([1-6])$", "DB\\1"),

    phenotype = case_when(
      str_detect(sample, "^WT") ~ "WT",
      str_detect(sample, "^DB") ~ "DB",
      TRUE ~ NA_character_
    ),

    reads_input = if ("new_est_reads" %in% colnames(.)) new_est_reads else NA_real_
  ) %>%
  filter(!is.na(phenotype))

## =========================================================
## 3. summarize
## =========================================================

summarize_rank_bracken <- function(df, level_name) {

  rank_df <- df %>%
    filter(level == level_name,
           !is.na(phenotype),
           !is.na(taxon),
           taxon != "")

  sample_tax <- rank_df %>%
    group_by(sample, sample_plot, phenotype, taxon) %>%
    summarise(
      abundance = sum(fraction_total_reads, na.rm = TRUE),
      reads = sum(reads_input, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(sample, sample_plot, phenotype) %>%
    mutate(abundance = abundance / sum(abundance, na.rm = TRUE)) %>%
    ungroup()

  group_tax <- sample_tax %>%
    group_by(phenotype, taxon) %>%
    summarise(
      abundance = mean(abundance, na.rm = TRUE),
      reads = sum(reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(group = phenotype) %>%
    arrange(group, desc(abundance))

  list(sample = sample_tax, group = group_tax)
}

family_res  <- summarize_rank_bracken(mgx_all, "family")
genus_res   <- summarize_rank_bracken(mgx_all, "genus")
species_res <- summarize_rank_bracken(mgx_all, "species")
species_hm_res <- summarize_rank_bracken(mgx_all, "species_harmonized")

## =========================================================
## 4. output tables
## =========================================================

make_wide_group <- function(df) {
  df %>%
    pivot_wider(names_from = group,
                values_from = c(abundance, reads),
                values_fill = 0)
}

make_wide_sample <- function(df) {
  df %>%
    select(sample_plot, taxon, abundance, reads) %>%
    pivot_wider(names_from = sample_plot,
                values_from = c(abundance, reads),
                values_fill = 0)
}

write_csv(species_res$sample, file.path(output_dir, "MGX_species_sample_abundance.csv"))
write_csv(genus_res$sample,   file.path(output_dir, "MGX_genus_sample_abundance.csv"))
write_csv(family_res$sample,  file.path(output_dir, "MGX_family_sample_abundance.csv"))
write_csv(species_hm_res$sample, file.path(output_dir, "MGX_species_harmonized_sample_abundance.csv"))

write_csv(species_res$group, file.path(output_dir, "MGX_species_group_abundance_all.csv"))
write_csv(genus_res$group,   file.path(output_dir, "MGX_genus_group_abundance_all.csv"))
write_csv(family_res$group,  file.path(output_dir, "MGX_family_group_abundance_all.csv"))
write_csv(species_hm_res$group, file.path(output_dir, "MGX_species_harmonized_group_abundance_all.csv"))

## =========================================================
## 5. composition
## =========================================================

build_mgx_comp_sample <- function(df, topN = 10) {

  df2 <- df %>%
    mutate(
      sample_plot = sample %>%
        str_remove("^MC_") %>%
        str_replace("^WT8([1-6])$", "WT\\1") %>%
        str_replace("^DB8([1-6])$", "DB\\1")
    )

  top_taxa <- df2 %>%
    group_by(taxon) %>%
    summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = topN) %>%
    pull(taxon)

  df2 %>%
    mutate(
      taxon = ifelse(taxon %in% top_taxa, taxon, "Others"),
      sample_plot = factor(sample_plot, levels = sample_order_all)
    ) %>%
    group_by(sample_plot, phenotype, taxon) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(sample_plot, phenotype) %>%
    mutate(Percent = abundance / sum(abundance) * 100) %>%
    ungroup() %>%
    rename(sample = sample_plot, position = phenotype)
}

mgx_comp_list <- list(
  Family = build_mgx_comp_sample(family_res$sample, topN_map$Family),
  Genus  = build_mgx_comp_sample(genus_res$sample, topN_map$Genus),
  Species = build_mgx_comp_sample(species_res$sample, topN_map$Species),
  Species_harmonized = build_mgx_comp_sample(species_hm_res$sample,
                                              topN_map$Species_harmonized)
)

## =========================================================
## 6. plotting
## =========================================================

plot_rank_nice <- function(df, rank_name,
                           sample_order = NULL,
                           position_order = NULL,
                           show_pct_labels = FALSE) {

  df <- df %>%
    mutate(
      position = factor(position, levels = position_order %||% unique(position)),
      sample = factor(sample, levels = sample_order %||% unique(sample)),
      taxon = fct_relevel(taxon, "Others", after = Inf)
    )

  pal <- make_named_palette(df)

  p <- ggplot(df, aes(sample, Percent, fill = taxon)) +
    geom_col(width = 0.85) +
    facet_wrap(~position, scales = "free_x") +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(
      title = paste0(rank_name, " composition (Top + Others)"),
      x = "Sample",
      y = "Percentage"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  if (show_pct_labels) {
    p <- p + geom_text(
      aes(label = ifelse(Percent >= 5, paste0(round(Percent,1), "%"), "")),
      position = position_stack(vjust = 0.5),
      size = 3,
      color = "white"
    )
  }

  p
}

p_family <- plot_rank_nice(mgx_comp_list$Family, "Family",
                            sample_order_all, c("WT","DB"))

p_genus <- plot_rank_nice(mgx_comp_list$Genus, "Genus",
                          sample_order_all, c("WT","DB"))

p_species <- plot_rank_nice(mgx_comp_list$Species, "Species",
                            sample_order_all, c("WT","DB"))

p_species_hm <- plot_rank_nice(mgx_comp_list$Species_harmonized,
                               "Species_harmonized",
                               sample_order_all, c("WT","DB"))

## =========================================================
## 7. save
## =========================================================

ggsave(file.path(plot_dir, "MGX_Family.png"), p_family, width = 5, height = 4.5, dpi = 600)
ggsave(file.path(plot_dir, "MGX_Genus.png"), p_genus, width = 5, height = 4.5, dpi = 600)
ggsave(file.path(plot_dir, "MGX_Species.png"), p_species, width = 6, height = 4.5, dpi = 600)
ggsave(file.path(plot_dir, "MGX_Species_harmonized.png"), p_species_hm,
       width = 6, height = 4.5, dpi = 600)

pdf(file.path(plot_dir, "MGX_all.pdf"), width = 6, height = 4)
print(p_family)
print(p_genus)
print(p_species)
print(p_species_hm)
dev.off()
