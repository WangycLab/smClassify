#!/usr/bin/env Rscript

## =========================================================
## Supplementary Figure 7a
## Metagenome Bracken composition barplot
##
## Input:
##   bracken_long.csv
##
## Output:
##   1) abundance summary tables
##   2) Family / Genus / Species / Species_harmonized composition barplots
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(RColorBrewer)
})

## =========================================================
## 0. paths
## =========================================================
base_dir   <- "Figure2"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

bracken_long_file <- file.path(input_dir, "bracken_long.csv")

outdir <- output_dir
plot_dir <- output_dir

## Fixed plotting order
sample_order_all <- c(
  "WT1", "WT2", "WT3", "WT4", "WT5", "WT6",
  "DB1", "DB2", "DB3", "DB4", "DB5", "DB6"
)

topN_map <- list(
  Family = 10,
  Genus = 15,
  Species = 20,
  Species_harmonized = 20
)

## =========================================================
## 1. helpers
## =========================================================
`%||%` <- function(x, y) if (!is.null(x)) x else y

rename_sample_simple <- function(x) {
  case_when(
    x == "WT81" ~ "WT1",
    x == "WT82" ~ "WT2",
    x == "WT83" ~ "WT3",
    x == "WT84" ~ "WT4",
    x == "WT85" ~ "WT5",
    x == "WT86" ~ "WT6",
    x == "DB81" ~ "DB1",
    x == "DB82" ~ "DB2",
    x == "DB83" ~ "DB3",
    x == "DB84" ~ "DB4",
    x == "DB85" ~ "DB5",
    x == "DB86" ~ "DB6",
    x == "MC_WT81" ~ "WT1",
    x == "MC_WT82" ~ "WT2",
    x == "MC_WT83" ~ "WT3",
    x == "MC_WT84" ~ "WT4",
    x == "MC_WT85" ~ "WT5",
    x == "MC_WT86" ~ "WT6",
    x == "MC_DB81" ~ "DB1",
    x == "MC_DB82" ~ "DB2",
    x == "MC_DB83" ~ "DB3",
    x == "MC_DB84" ~ "DB4",
    x == "MC_DB85" ~ "DB5",
    x == "MC_DB86" ~ "DB6",
    TRUE ~ x
  )
}

make_named_palette <- function(df, base_cols = NULL) {
  taxa <- unique(df$taxon)
  taxa <- taxa[!is.na(taxa)]
  
  if ("Others" %in% taxa) {
    taxa_main <- setdiff(taxa, "Others")
  } else {
    taxa_main <- taxa
  }
  
  n <- length(taxa_main)
  
  if (is.null(base_cols)) {
    if (n <= 8) {
      cols <- brewer.pal(max(3, n), "Set2")[seq_len(n)]
    } else if (n <= 12) {
      cols <- brewer.pal(12, "Set3")[seq_len(n)]
    } else {
      cols <- colorRampPalette(c(
          "#1F77B4", "#2CA02C", "#FF7F0E", "#6A5ACD", "#8C564B", "#E377C2", "#7F7F7F",
          "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32",
          "#8B0000", "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B",
          "#20B2AA", "#FF1493", "#FF4500", "#32CD32", "#3E8E41", "#20B2AA"
      ))(n)
    }
  } else {
    if (length(base_cols) < n) {
      cols <- colorRampPalette(base_cols)(n)
    } else {
      cols <- base_cols[seq_len(n)]
    }
  }
  
  pal <- setNames(cols, taxa_main)
  
  if ("Others" %in% taxa) {
    pal <- c(pal, Others = "grey80")
  }
  
  pal
}
bracken_long <- read.csv(bracken_long_file)
## =========================================================
## 2. prepare bracken_long


mgx_all <- bracken_long %>%
  mutate(
    sample = as.character(sample),
    sample_plot = rename_sample_simple(sample),
    phenotype = case_when(
      str_detect(sample, "^WT") ~ "WT",
      str_detect(sample, "^DB") ~ "DB",
      str_detect(sample, "WT") ~ "WT",
      str_detect(sample, "DB") ~ "DB",
      TRUE ~ NA_character_
    ),
    reads_input = if ("new_est_reads" %in% colnames(.)) new_est_reads else NA_real_
  ) %>%
  filter(!is.na(phenotype))

## =========================================================
## 3. summarize by rank
## =========================================================
summarize_rank_bracken <- function(df, level_name) {
  stopifnot(all(c("sample", "sample_plot", "phenotype", "level", "taxon", "fraction_total_reads") %in% colnames(df)))
  
  rank_df <- df %>%
    dplyr::filter(.data$level == level_name) %>%
    dplyr::filter(!is.na(.data$phenotype)) %>%
    dplyr::filter(!is.na(.data$taxon), .data$taxon != "")
  
  # sample-level abundance
  sample_tax <- rank_df %>%
    dplyr::group_by(.data$sample, .data$sample_plot, .data$phenotype, .data$taxon) %>%
    dplyr::summarise(
      abundance = sum(.data$fraction_total_reads, na.rm = TRUE),
      reads = if ("reads_input" %in% colnames(rank_df)) {
        if (all(is.na(rank_df$reads_input))) NA_real_ else sum(.data$reads_input, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$sample, .data$sample_plot, .data$phenotype) %>%
    dplyr::mutate(abundance = .data$abundance / sum(.data$abundance, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  # group-level abundance
  group_tax <- sample_tax %>%
    dplyr::group_by(.data$phenotype, .data$taxon) %>%
    dplyr::summarise(
      abundance = mean(.data$abundance, na.rm = TRUE),
      reads = if (all(is.na(.data$reads))) NA_real_ else sum(.data$reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(group = .data$phenotype) %>%
    dplyr::arrange(.data$group, dplyr::desc(.data$abundance))
  
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
    pivot_wider(
      names_from = group,
      values_from = c(abundance, reads),
      values_fill = 0
    )
}

make_wide_sample <- function(df) {
  df %>%
    select(sample_plot, taxon, abundance, reads) %>%
    pivot_wider(
      names_from = sample_plot,
      values_from = c(abundance, reads),
      values_fill = 0
    )
}

## sample-level
write_csv(species_res$sample, file.path(outdir, "MGX_species_sample_abundance.csv"))
write_csv(genus_res$sample,   file.path(outdir, "MGX_genus_sample_abundance.csv"))
write_csv(family_res$sample,  file.path(outdir, "MGX_family_sample_abundance.csv"))
write_csv(species_hm_res$sample, file.path(outdir, "MGX_species_harmonized_sample_abundance.csv"))

## group-level
write_csv(species_res$group, file.path(outdir, "MGX_species_group_abundance_all.csv"))
write_csv(genus_res$group,   file.path(outdir, "MGX_genus_group_abundance_all.csv"))
write_csv(family_res$group,  file.path(outdir, "MGX_family_group_abundance_all.csv"))
write_csv(species_hm_res$group, file.path(outdir, "MGX_species_harmonized_group_abundance_all.csv"))

## WT / DB split group-level
write_csv(species_res$group %>% filter(group == "WT") %>% select(-group),
          file.path(outdir, "MGX_WT_species_abundance.csv"))
write_csv(species_res$group %>% filter(group == "DB") %>% select(-group),
          file.path(outdir, "MGX_DB_species_abundance.csv"))

write_csv(genus_res$group %>% filter(group == "WT") %>% select(-group),
          file.path(outdir, "MGX_WT_genus_abundance.csv"))
write_csv(genus_res$group %>% filter(group == "DB") %>% select(-group),
          file.path(outdir, "MGX_DB_genus_abundance.csv"))

write_csv(family_res$group %>% filter(group == "WT") %>% select(-group),
          file.path(outdir, "MGX_WT_family_abundance.csv"))
write_csv(family_res$group %>% filter(group == "DB") %>% select(-group),
          file.path(outdir, "MGX_DB_family_abundance.csv"))

write_csv(species_hm_res$group %>% filter(group == "WT") %>% select(-group),
          file.path(outdir, "MGX_WT_species_harmonized_abundance.csv"))
write_csv(species_hm_res$group %>% filter(group == "DB") %>% select(-group),
          file.path(outdir, "MGX_DB_species_harmonized_abundance.csv"))

## WT / DB split sample-level
write_csv(species_res$sample %>% filter(phenotype == "WT"),
          file.path(outdir, "MGX_WT_species_sample_abundance.csv"))
write_csv(species_res$sample %>% filter(phenotype == "DB"),
          file.path(outdir, "MGX_DB_species_sample_abundance.csv"))

write_csv(genus_res$sample %>% filter(phenotype == "WT"),
          file.path(outdir, "MGX_WT_genus_sample_abundance.csv"))
write_csv(genus_res$sample %>% filter(phenotype == "DB"),
          file.path(outdir, "MGX_DB_genus_sample_abundance.csv"))

write_csv(family_res$sample %>% filter(phenotype == "WT"),
          file.path(outdir, "MGX_WT_family_sample_abundance.csv"))
write_csv(family_res$sample %>% filter(phenotype == "DB"),
          file.path(outdir, "MGX_DB_family_sample_abundance.csv"))

write_csv(species_hm_res$sample %>% filter(phenotype == "WT"),
          file.path(outdir, "MGX_WT_species_harmonized_sample_abundance.csv"))
write_csv(species_hm_res$sample %>% filter(phenotype == "DB"),
          file.path(outdir, "MGX_DB_species_harmonized_sample_abundance.csv"))

## wide tables
write_csv(make_wide_group(species_res$group), file.path(outdir, "MGX_species_abundance_wide_group.csv"))
write_csv(make_wide_group(genus_res$group),   file.path(outdir, "MGX_genus_abundance_wide_group.csv"))
write_csv(make_wide_group(family_res$group),  file.path(outdir, "MGX_family_abundance_wide_group.csv"))
write_csv(make_wide_group(species_hm_res$group), file.path(outdir, "MGX_species_harmonized_abundance_wide_group.csv"))

write_csv(make_wide_sample(species_res$sample), file.path(outdir, "MGX_species_abundance_wide_sample.csv"))
write_csv(make_wide_sample(genus_res$sample),   file.path(outdir, "MGX_genus_abundance_wide_sample.csv"))
write_csv(make_wide_sample(family_res$sample),  file.path(outdir, "MGX_family_abundance_wide_sample.csv"))
write_csv(make_wide_sample(species_hm_res$sample), file.path(outdir, "MGX_species_harmonized_abundance_wide_sample.csv"))

## =========================================================
## 5. build composition data
## =========================================================
build_mgx_comp_sample <- function(df, topN = 10) {
  stopifnot(all(c("sample", "phenotype", "taxon", "abundance") %in% colnames(df)))
  
  df2 <- df %>%
    mutate(
      sample_plot = if ("sample_plot" %in% colnames(.)) {
        as.character(sample_plot)
      } else {
        rename_sample_simple(as.character(sample))
      }
    )
  
  top_taxa <- df2 %>%
    dplyr::group_by(.data$taxon) %>%
    dplyr::summarise(total = sum(.data$abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$total)) %>%
    dplyr::slice_head(n = topN) %>%
    dplyr::pull(.data$taxon)
  
  df2 %>%
    dplyr::mutate(
      taxon = ifelse(.data$taxon %in% top_taxa, .data$taxon, "Others"),
      sample_plot = factor(.data$sample_plot, levels = sample_order_all)
    ) %>%
    dplyr::group_by(.data$sample_plot, .data$phenotype, .data$taxon) %>%
    dplyr::summarise(abundance = sum(.data$abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data$sample_plot, .data$phenotype) %>%
    dplyr::mutate(Percent = .data$abundance / sum(.data$abundance) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::rename(sample = .data$sample_plot, position = .data$phenotype)
}

mgx_comp_list <- list(
  Family  = build_mgx_comp_sample(family_res$sample,  topN = topN_map$Family),
  Genus   = build_mgx_comp_sample(genus_res$sample,   topN = topN_map$Genus),
  Species = build_mgx_comp_sample(species_res$sample, topN = topN_map$Species),
  Species_harmonized = build_mgx_comp_sample(
    species_hm_res$sample,
    topN = topN_map$Species_harmonized
  )
)
## =========================================================
## 6. plotting function
## =========================================================


plot_rank_nice <- function(df, rank_name, base_cols = NULL,
                           sample_order = NULL, position_order = NULL,
                           show_pct_labels = FALSE) {
  df <- df %>%
    mutate(
      position = factor(position, levels = position_order %||% unique(position)),
      sample   = factor(sample, levels = sample_order %||% unique(sample)),
      taxon    = fct_relevel(taxon, "Others", after = Inf)
    ) %>%
    arrange(position, sample)
  
  pal <- make_named_palette(df, base_cols)
  
  p <- ggplot(df, aes(sample, Percent, fill = taxon)) +
    geom_col(width = 0.85) +
    facet_wrap(~position, scales = "free_x") +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(
      title = paste0(rank_name, " composition (Top", topN_map[[rank_name]], " + Others)"),
      x = "Sample",
      y = "Percentage",
      fill = "Taxon"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  if (show_pct_labels) {
    p <- p + geom_text(
      aes(label = ifelse(Percent >= 5, paste0(round(Percent, 1), "%"), "")),
      position = position_stack(vjust = 0.5),
      size = 3,
      color = "white"
    )
  }
  
  p
}

## =========================================================
## 7. plot
## =========================================================
p_family_all <- plot_rank_nice(
  mgx_comp_list$Family,
  "Family",
  sample_order = sample_order_all,
  position_order = c("WT", "DB")
)

p_genus_all <- plot_rank_nice(
  mgx_comp_list$Genus,
  "Genus",
  sample_order = sample_order_all,
  position_order = c("WT", "DB")
)

p_species_all <- plot_rank_nice(
  mgx_comp_list$Species,
  "Species",
  sample_order = sample_order_all,
  position_order = c("WT", "DB")
)

p_species_hm_all <- plot_rank_nice(
  mgx_comp_list$Species_harmonized,
  "Species_harmonized",
  sample_order = sample_order_all,
  position_order = c("WT", "DB")
)
p_species_hm_all
## =========================================================
## 8. save plots
## =========================================================
pdf(file.path(plot_dir, "MGX_Family_Genus_Species_all_samples.pdf"), width = 5, height = 4)
print(p_family_all)
print(p_genus_all)
print(p_species_all)
dev.off()

pdf(file.path(plot_dir, "MGX_Family_all_samples.pdf"), width = 5, height = 4)
print(p_family_all)
dev.off()

pdf(file.path(plot_dir, "MGX_Genus_all_samples.pdf"), width = 5, height = 4)
print(p_genus_all)
dev.off()

pdf(file.path(plot_dir, "MGX_Species_all_samples.pdf"), width = 6, height = 4)
print(p_species_all)
dev.off()

ggsave(file.path(plot_dir, "MGX_Family_all_samples.png"), p_family_all, width = 5, height = 4.5, dpi = 600)
ggsave(file.path(plot_dir, "MGX_Genus_all_samples.png"), p_genus_all, width = 5, height = 4.5, dpi = 600)
ggsave(file.path(plot_dir, "MGX_Species_all_samples.png"), p_species_all, width = 8, height = 4.5, dpi = 600)

pdf(file.path(plot_dir, "MGX_Species_harmonized_all_samples.pdf"), width = 6, height = 4)
print(p_species_hm_all)
dev.off()

ggsave(file.path(plot_dir, "MGX_Species_harmonized_all_samples.png"),
       p_species_hm_all, width = 8, height = 4.5, dpi = 600)
