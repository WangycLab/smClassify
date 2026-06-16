
# =========================================================
# Figure 2i
# Cecum single-cell vs metagenome (Bracken MGnify) comparison

# =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(scales)
})

# =========================================================
# 0. paths
# =========================================================

base_dir   <- "Figure2"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sce_file          <- file.path(input_dir, "Mouse_microbial_taxonomy_smClassify_Kraken2_MGnify.rds")
mgnify_meta_file  <- file.path(input_dir, "mgnify_metadata.csv")
bracken_long_file <- file.path(input_dir, "Metagenome_taxonomic_abundance.csv")

outdir <- output_dir

# =========================================================
# 1. theme
# =========================================================

theme_cns <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# =========================================================
# 2. read data
# =========================================================

sce_bac_annot <- readRDS(sce_file)
mgnify_meta <- read_csv(mgnify_meta_file, show_col_types = FALSE)
bracken_long <- read.csv(bracken_long_file)

mgx_bracken_long <- bracken_long %>%
  mutate(
    sample = as.character(sample),
    group = case_when(
      str_detect(sample, "^WT") ~ "WT",
      str_detect(sample, "^DB") ~ "DB",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group))

# =========================================================
# 3. prepare single-cell meta (Cecum)
# =========================================================

sc_meta_cecum <- sce_bac_annot@meta.data %>%
  dplyr::filter(loc == "Cecum") %>%
  mutate(
    group = case_when(
      str_detect(orig.ident, "WT") ~ "WT",
      str_detect(orig.ident, "DB") ~ "DB",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(group))

# =========================================================
# 4. MGnify metadata mapping
# =========================================================

mgnify_family_map <- mgnify_meta %>%
  select(
    mgnify_species,
    MGnify_Family = Family,
    MGnify_Genus  = Genus
  ) %>%
  distinct()

sc_meta_cecum <- sc_meta_cecum %>%
  left_join(
    mgnify_family_map,
    by = c("MGnify_species" = "mgnify_species")
  )

# =========================================================
# 5. helper: summarize abundance
# =========================================================

summarize_sc_group <- function(df, tax_col) {
  df %>%
    select(group, taxon = all_of(tax_col)) %>%
    mutate(taxon = as.character(taxon)) %>%
    filter(!is.na(group), !is.na(taxon), taxon != "") %>%
    count(group, taxon, name = "n") %>%
    group_by(group) %>%
    mutate(abundance = n / sum(n)) %>%
    ungroup()
}

summarize_mgx_group <- function(df, level_name) {
  df %>%
    filter(level == level_name) %>%
    mutate(taxon = as.character(taxon)) %>%
    filter(!is.na(taxon), taxon != "") %>%
    group_by(group, taxon) %>%
    summarise(
      abundance = sum(fraction_total_reads),
      .groups = "drop"
    ) %>%
    group_by(group) %>%
    mutate(abundance = abundance / sum(abundance)) %>%
    ungroup()
}

# =========================================================
# 6. build abundance tables
# =========================================================

# smClassify
sc_sm_family  <- summarize_sc_group(sc_meta_cecum, "Tax_Family")
sc_sm_genus   <- summarize_sc_group(sc_meta_cecum, "smClassify_genus")
sc_sm_species <- summarize_sc_group(sc_meta_cecum, "smClassify_species")

# sc-Kraken2(MGnify)
sc_kr_family  <- summarize_sc_group(sc_meta_cecum, "MGnify_Family")
sc_kr_genus   <- summarize_sc_group(sc_meta_cecum, "MGnify_Genus")
sc_kr_species <- summarize_sc_group(sc_meta_cecum, "MGnify_species")

# metagenome Bracken
mgx_family     <- summarize_mgx_group(mgx_bracken_long, "family")
mgx_genus      <- summarize_mgx_group(mgx_bracken_long, "genus")
mgx_species    <- summarize_mgx_group(mgx_bracken_long, "species")
mgx_species_hm <- summarize_mgx_group(mgx_bracken_long, "species_harmonized")

write_csv(sc_sm_family,  file.path(outdir, "Cecum_smClassify_family_group_abundance.csv"))
write_csv(sc_sm_genus,   file.path(outdir, "Cecum_smClassify_genus_group_abundance.csv"))
write_csv(sc_sm_species, file.path(outdir, "Cecum_smClassify_species_group_abundance.csv"))

write_csv(sc_kr_family,  file.path(outdir, "Cecum_sm_Kraken2_MGnify_family_group_abundance.csv"))
write_csv(sc_kr_genus,   file.path(outdir, "Cecum_sm_Kraken2_MGnify_genus_group_abundance.csv"))
write_csv(sc_kr_species, file.path(outdir, "Cecum_sm_Kraken2_MGnify_species_group_abundance.csv"))

# =========================================================
# 7. shared-taxa correlation by rank
# =========================================================

run_shared_correlation <- function(sc_df, mgx_df, prefix, tax_label, sc_label, comparison_label) {

  groups_use <- sort(intersect(unique(sc_df$group), unique(mgx_df$group)))
  stats_list <- list()

  for (g in groups_use) {

    sc_sub <- sc_df %>%
      dplyr::filter(group == g) %>%
      mutate(taxon = as.character(taxon)) %>%
      dplyr::select(taxon, sc_abundance = abundance)

    mgx_sub <- mgx_df %>%
      dplyr::filter(group == g, abundance > 1e-8) %>%
      mutate(taxon = as.character(taxon)) %>%
      dplyr::select(taxon, mgx_abundance = abundance)

    cor_df <- inner_join(sc_sub, mgx_sub, by = "taxon") %>%
      mutate(
        log_sc = log10(sc_abundance + 1e-8),
        log_mgx = log10(mgx_abundance + 1e-8)
      ) %>%
      arrange(desc(sc_abundance + mgx_abundance))

    write_csv(
      cor_df,
      file.path(outdir, paste0(prefix, "_", g, "_shared_taxa_table.csv"))
    )

    if (nrow(cor_df) < 3) next

    pearson_raw <- suppressWarnings(cor.test(
      cor_df$sc_abundance,
      cor_df$mgx_abundance,
      method = "pearson"
    ))

    pearson_log <- suppressWarnings(cor.test(
      cor_df$log_sc,
      cor_df$log_mgx,
      method = "pearson"
    ))

    spearman_res <- suppressWarnings(cor.test(
      cor_df$sc_abundance,
      cor_df$mgx_abundance,
      method = "spearman",
      exact = FALSE
    ))

    lims_all <- c(cor_df$log_sc, cor_df$log_mgx)
    lim_min <- min(lims_all)
    lim_max <- max(lims_all)

    p <- ggplot(cor_df, aes(x = log_sc, y = log_mgx)) +
      geom_point(
        size = 2.4,
        alpha = 0.9,
        colour = "black",
        fill = "#1F77B4",
        shape = 21,
        stroke = 0.3
      ) +
      geom_smooth(method = "lm", se = FALSE, colour = "#D55E00", linewidth = 0.7) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      coord_fixed(
        xlim = c(lim_min, lim_max),
        ylim = c(lim_min, lim_max)
      ) +
      labs(
        title = paste0(tax_label, " correlation - ", g),
        x = paste0("log10(", sc_label, ")"),
        y = "log10(Metagenome Kraken2)"
      ) +
      annotate(
        "text",
        x = lim_max,
        y = lim_min,
        label = paste0(
          "Pearson(log) = ", round(unname(pearson_log$estimate), 3),
          "\nSpearman = ", round(unname(spearman_res$estimate), 3),
          "\nN = ", nrow(cor_df)
        ),
        hjust = 1,
        vjust = 0,
        size = 4
      ) +
      theme_cns()

    ggsave(
      file.path(outdir, paste0(prefix, "_", g, "_shared_correlation.pdf")),
      p,
      width = 4.3,
      height = 4
    )
    ggsave(
      file.path(outdir, paste0(prefix, "_", g, "_shared_correlation.png")),
      p,
      width = 4.3,
      height = 4,
      dpi = 600
    )

    stats_list[[g]] <- data.frame(
      Comparison = comparison_label,
      Condition = g,
      Taxonomic_rank = tax_label,
      Number_of_taxa = nrow(cor_df),
      Pearson_r_relative_abundance = unname(pearson_raw$estimate),
      Pearson_P_value_relative_abundance = pearson_raw$p.value,
      Pearson_r_log_abundance = unname(pearson_log$estimate),
      Pearson_P_value_log_abundance = pearson_log$p.value,
      Spearman_rho = unname(spearman_res$estimate),
      Spearman_P_value = spearman_res$p.value
    )
  }

  stats_df <- bind_rows(stats_list)

  write_csv(
    stats_df,
    file.path(outdir, paste0(prefix, "_shared_correlation_summary.csv"))
  )

  stats_df
}

# =========================================================
# 8. run correlation comparisons
# =========================================================

stats_sm_family <- run_shared_correlation(
  sc_df = sc_sm_family,
  mgx_df = mgx_family,
  prefix = "Cecum_smClassify_vs_MGX_family",
  tax_label = "Family",
  sc_label = "smClassify",
  comparison_label = "smClassify_vs_Metagenome"
)

stats_sm_genus <- run_shared_correlation(
  sc_df = sc_sm_genus,
  mgx_df = mgx_genus,
  prefix = "Cecum_smClassify_vs_MGX_genus",
  tax_label = "Genus",
  sc_label = "smClassify",
  comparison_label = "smClassify_vs_Metagenome"
)

stats_sm_species <- run_shared_correlation(
  sc_df = sc_sm_species,
  mgx_df = mgx_species_hm,
  prefix = "Cecum_smClassify_vs_MGX_species",
  tax_label = "Species",
  sc_label = "smClassify",
  comparison_label = "smClassify_vs_Metagenome"
)

stats_kr_family <- run_shared_correlation(
  sc_df = sc_kr_family,
  mgx_df = mgx_family,
  prefix = "Cecum_sm_Kraken2MGnify_vs_MGX_family",
  tax_label = "Family",
  sc_label = "sc-Kraken2(MGnify)",
  comparison_label = "sm_Kraken2_MGnify_vs_Metagenome"
)

stats_kr_genus <- run_shared_correlation(
  sc_df = sc_kr_genus,
  mgx_df = mgx_genus,
  prefix = "Cecum_sm_Kraken2MGnify_vs_MGX_genus",
  tax_label = "Genus",
  sc_label = "sc-Kraken2(MGnify)",
  comparison_label = "sm_Kraken2_MGnify_vs_Metagenome"
)

stats_kr_species <- run_shared_correlation(
  sc_df = sc_kr_species,
  mgx_df = mgx_species,
  prefix = "Cecum_sm_Kraken2MGnify_vs_MGX_species",
  tax_label = "Species",
  sc_label = "sc-Kraken2(MGnify)",
  comparison_label = "sm_Kraken2_MGnify_vs_Metagenome"
)

all_stats <- bind_rows(
  stats_sm_family,
  stats_sm_genus,
  stats_sm_species,
  stats_kr_family,
  stats_kr_genus,
  stats_kr_species
)

write_csv(
  all_stats,
  file.path(outdir, "Cecum_All_SC_vs_MGX_shared_correlation_summary.csv")
)

all_stats_ordered <- all_stats %>%
  mutate(
    Taxonomic_rank = factor(Taxonomic_rank, levels = c("Family", "Genus", "Species")),
    Condition = factor(Condition, levels = c("WT", "DB"))
  ) %>%
  arrange(Comparison, Condition, Taxonomic_rank)

write_csv(
  all_stats_ordered,
  file.path(outdir, "Cecum_All_SC_vs_MGX_shared_correlation_summary_ordered.csv")
)

# =========================================================
# 9. overlap comparison
# =========================================================

run_overlap_comparison <- function(sc_df, mgx_df, prefix, tax_label, comparison_label, sc_only_label = "smClassify_only") {

  groups_use <- sort(intersect(unique(sc_df$group), unique(mgx_df$group)))
  res_list <- list()

  for (g in groups_use) {
    sc_set <- sc_df %>%
      dplyr::filter(group == g) %>%
      mutate(taxon = as.character(taxon)) %>%
      pull(taxon) %>%
      unique()

    mgx_set <- mgx_df %>%
      dplyr::filter(group == g, abundance > 1e-8) %>%
      mutate(taxon = as.character(taxon)) %>%
      pull(taxon) %>%
      unique()

    n_shared <- length(intersect(sc_set, mgx_set))
    n_sc_only <- length(setdiff(sc_set, mgx_set))
    n_mgx_only <- length(setdiff(mgx_set, sc_set))

    df_plot <- data.frame(
      Category = c("Shared", sc_only_label, "Metagenome_only"),
      Count = c(n_shared, n_sc_only, n_mgx_only)
    )

    write_csv(
      df_plot,
      file.path(outdir, paste0(prefix, "_", g, "_overlap_counts.csv"))
    )

    p <- ggplot(df_plot, aes(x = Category, y = Count, fill = Category)) +
      geom_col(width = 0.6, colour = "black") +
      geom_text(
        aes(label = Count),
        position = position_stack(vjust = 0.5),
        size = 4,
        colour = "black"
      ) +
      scale_fill_manual(values = c(
        "Shared" = "#009E73",
        "smClassify_only" = "#BCBD22",
        "smKraken2_only" = "#BCBD22",
        "Metagenome_only" = "#F5A0A1"
      )) +
      labs(
        title = paste0(tax_label, " overlap - ", g),
        x = NULL,
        y = "Number of taxa"
      ) +
      theme_cns() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14)
      )

    ggsave(
      file.path(outdir, paste0(prefix, "_", g, "_overlap.pdf")),
      p,
      width = 3.2,
      height = 4
    )
    ggsave(
      file.path(outdir, paste0(prefix, "_", g, "_overlap.png")),
      p,
      width = 3.2,
      height = 4,
      dpi = 600
    )

    res_list[[g]] <- data.frame(
      Comparison = comparison_label,
      Condition = g,
      Taxonomic_rank = tax_label,
      Shared = n_shared,
      SingleCell_only = n_sc_only,
      Metagenome_only = n_mgx_only
    )
  }

  res_df <- bind_rows(res_list)

  write_csv(
    res_df,
    file.path(outdir, paste0(prefix, "_overlap_summary.csv"))
  )

  return(res_df)
}

ov_sm_family <- run_overlap_comparison(
  sc_sm_family, mgx_family,
  "Cecum_smClassify_vs_MGX_family",
  "Family",
  "smClassify_vs_Metagenome",
  sc_only_label = "smClassify_only"
)

ov_sm_genus <- run_overlap_comparison(
  sc_sm_genus, mgx_genus,
  "Cecum_smClassify_vs_MGX_genus",
  "Genus",
  "smClassify_vs_Metagenome",
  sc_only_label = "smClassify_only"
)

ov_sm_species <- run_overlap_comparison(
  sc_sm_species, mgx_species_hm,
  "Cecum_smClassify_vs_MGX_species",
  "Species",
  "smClassify_vs_Metagenome",
  sc_only_label = "smClassify_only"
)

ov_kr_family <- run_overlap_comparison(
  sc_kr_family, mgx_family,
  "Cecum_sm_Kraken2MGnify_vs_MGX_family",
  "Family",
  "sm_Kraken2_MGnify_vs_Metagenome",
  sc_only_label = "smKraken2_only"
)

ov_kr_genus <- run_overlap_comparison(
  sc_kr_genus, mgx_genus,
  "Cecum_sm_Kraken2MGnify_vs_MGX_genus",
  "Genus",
  "sm_Kraken2_MGnify_vs_Metagenome",
  sc_only_label = "smKraken2_only"
)

ov_kr_species <- run_overlap_comparison(
  sc_kr_species, mgx_species,
  "Cecum_sm_Kraken2MGnify_vs_MGX_species",
  "Species",
  "sm_Kraken2_MGnify_vs_Metagenome",
  sc_only_label = "smKraken2_only"
)

overlap_all <- bind_rows(
  ov_sm_family, ov_sm_genus, ov_sm_species,
  ov_kr_family, ov_kr_genus, ov_kr_species
)

write_csv(
  overlap_all,
  file.path(outdir, "Cecum_all_overlap_summary.csv")
)

# =========================================================
# 10. top10-union log2FC between single-cell and Metagenome
# =========================================================

run_top10_log2fc_union <- function(sc_df, mgx_df, prefix, tax_label, sc_vs_mgx_label = "smClassify / Metagenome") {

  groups_use <- sort(intersect(unique(sc_df$group), unique(mgx_df$group)))
  pseudo <- 1e-8

  for (g in groups_use) {
    sc_sub <- sc_df %>%
      dplyr::filter(group == g) %>%
      mutate(taxon = as.character(taxon)) %>%
      dplyr::select(taxon, sc_abundance = abundance)

    mgx_sub <- mgx_df %>%
      dplyr::filter(group == g, abundance > 1e-8) %>%
      mutate(taxon = as.character(taxon)) %>%
      dplyr::select(taxon, mgx_abundance = abundance)

    top10_sc <- sc_sub %>%
      arrange(desc(sc_abundance)) %>%
      slice_head(n = 10) %>%
      pull(taxon)

    top10_mgx <- mgx_sub %>%
      arrange(desc(mgx_abundance)) %>%
      slice_head(n = 10) %>%
      pull(taxon)

    top_union <- union(top10_sc, top10_mgx)

    df_fc <- full_join(sc_sub, mgx_sub, by = "taxon") %>%
      mutate(
        sc_abundance = ifelse(is.na(sc_abundance), 0, sc_abundance),
        mgx_abundance = ifelse(is.na(mgx_abundance), 0, mgx_abundance)
      ) %>%
      dplyr::filter(taxon %in% top_union) %>%
      mutate(
        Category = case_when(
          sc_abundance > 0 & mgx_abundance > 0 ~ "Shared",
          sc_abundance > 0 & mgx_abundance == 0 ~ "SingleCell_only",
          sc_abundance == 0 & mgx_abundance > 0 ~ "Metagenome_only"
        ),
        log2FC = log2((sc_abundance + pseudo) / (mgx_abundance + pseudo))
      )

    write_csv(
      df_fc,
      file.path(outdir, paste0(prefix, "_", g, "_top10_union_log2FC_table.csv"))
    )

    df_plot <- df_fc %>%
      arrange(log2FC) %>%
      mutate(
        taxon = factor(taxon, levels = taxon),
        Category = factor(
          Category,
          levels = c("Shared", "SingleCell_only", "Metagenome_only")
        )
      )

    p <- ggplot(df_plot, aes(x = taxon, y = log2FC, fill = Category)) +
      geom_col(
        width = 0.65,
        colour = "black",
        linewidth = 0.3
      ) +
      geom_text(
        aes(label = round(log2FC, 2)),
        hjust = ifelse(df_plot$log2FC > 0, -0.2, 1.2),
        size = 3
      ) +
      scale_fill_manual(values = c(
        "Shared" = "#009E73",
        "SingleCell_only" = "#BCBD22",
        "Metagenome_only" = "#F5A0A1"
      )) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      coord_flip() +
      labs(
        title = paste0(tax_label, " top10-union log2FC - ", g),
        x = NULL,
        y = paste0("log2FC [", sc_vs_mgx_label, "]"),
        fill = NULL
      ) +
      theme_cns() +
      theme(
        legend.position = "top"
      )

    ggsave(
      file.path(outdir, paste0(prefix, "_", g, "_top10_union_log2FC.pdf")),
      p,
      width = 7,
      height = 5.5
    )

    ggsave(
      file.path(outdir, paste0(prefix, "_", g, "_top10_union_log2FC.png")),
      p,
      width = 7,
      height = 5.5,
      dpi = 600
    )
  }
}

run_top10_log2fc_union(
  sc_sm_family,
  mgx_family,
  "Cecum_smClassify_vs_MGX_family",
  "Family",
  "smClassify / Metagenome"
)

run_top10_log2fc_union(
  sc_sm_genus,
  mgx_genus,
  "Cecum_smClassify_vs_MGX_genus",
  "Genus",
  "smClassify / Metagenome"
)

run_top10_log2fc_union(
  sc_sm_species,
  mgx_species_hm,
  "Cecum_smClassify_vs_MGX_species",
  "Species",
  "smClassify / Metagenome"
)

run_top10_log2fc_union(
  sc_kr_family,
  mgx_family,
  "Cecum_sm_Kraken2MGnify_vs_MGX_family",
  "Family",
  "scKraken2 / Metagenome"
)

run_top10_log2fc_union(
  sc_kr_genus,
  mgx_genus,
  "Cecum_sm_Kraken2MGnify_vs_MGX_genus",
  "Genus",
  "scKraken2 / Metagenome"
)

run_top10_log2fc_union(
  sc_kr_species,
  mgx_species,
  "Cecum_sm_Kraken2MGnify_vs_MGX_species",
  "Species",
  "scKraken2 / Metagenome"
)

cat("\nDone. Results written to:\n", outdir, "\n")
