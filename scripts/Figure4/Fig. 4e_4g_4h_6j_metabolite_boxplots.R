############################################################
## Figure 4e, 4g, 4h and Figure 6j metabolite boxplots
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggpubr)
})

# ======================================================
# 0) Input and output paths
# ======================================================

base_dir   <- "Figure4"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

serum_file   <- file.path(input_dir, "serum_norm_ms2_metabolites_intensity.xlsx")
content_file <- file.path(input_dir, "content_norm_ms2_metabolites_intensity.xlsx")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

source_out_file <- file.path(
  output_dir,
  "Fig4e_4g_4h_6j_SourceData_long_ALL8.csv"
)

# ======================================================
# 1) Targets
# ======================================================

# Fig. 4e ratios
fig4e_ratio_targets <- c(
  "Gln/Glu",
  "NAG/Glu"
)

fig4e_component_metabolites <- c(
  "Glutamine",
  "Glutamic acid",
  "N-Acetyl-L-glutamic acid"
)

# Fig. 4g components
# Even LCFA = C16:0 + C16:1 + C18:1
# Odd-chain FA = C15:0 + C17:0
fig4g_component_metabolites <- c(
  "Palmitic acid",        # C16:0
  "Palmitoleic acid",     # C16:1
  "Oleic acid",           # C18:1
  "Pentadecanoic acid",   # C15:0
  "Heptadecanoic acid"    # C17:0
)

# Fig. 4h conjugated bile acids
fig4h_metabolites <- c(
  "Taurodeoxycholic acid",
  "Taurochenodesoxycholic acid"
)

# Fig. 6j polyamine intermediates
fig6j_metabolites <- c(
  "N-Acetylputrescine",
  "N-Acetylcadaverine",
  "Sinapoylspermine",
  "N8-Acetylspermidine"
)

all_target_metabolites <- unique(c(
  fig4e_component_metabolites,
  fig4g_component_metabolites,
  fig4h_metabolites,
  fig6j_metabolites
))

# ======================================================
# 2) Fixed x-axis groups
# ======================================================

group_levels_all <- c(
  "Cecum-WT",  "Cecum-DB",
  "Colon-WT",  "Colon-DB",
  "Rectum-WT", "Rectum-DB",
  "Serum-WT",  "Serum-DB"
)

compare_list_all <- list(
  c("Cecum-WT",  "Cecum-DB"),
  c("Colon-WT",  "Colon-DB"),
  c("Rectum-WT", "Rectum-DB"),
  c("Serum-WT",  "Serum-DB")
)

# ======================================================
# 3) Read and reshape metabolomics data
# ======================================================

feature_cols <- c(
  "Type", "ID", "MZ", "RT", "MS2.name", "MS2_score",
  "Formula", "SuperClass", "Class", "Subclass",
  "HMDB", "kegg", "kegg_pathway", "MS1.name"
)

serum   <- read_xlsx(serum_file)
content <- read_xlsx(content_file)

names(serum) <- str_replace_all(names(serum), "^Seurm", "Serum")

serum_long <- serum %>%
  pivot_longer(
    cols = -all_of(feature_cols),
    names_to = "Sample",
    values_to = "Intensity"
  )

content_long <- content %>%
  pivot_longer(
    cols = -all_of(feature_cols),
    names_to = "Sample",
    values_to = "Intensity"
  )

all_long <- bind_rows(serum_long, content_long)

dat_long <- all_long %>%
  rename(Metabolite = MS2.name) %>%
  mutate(
    Metabolite = ifelse(is.na(Metabolite) | !nzchar(Metabolite), MS1.name, Metabolite),
    Metabolite = str_trim(Metabolite),
    Treatment = case_when(
      str_detect(Sample, "WT") ~ "WT",
      str_detect(Sample, "DB|T2DM") ~ "DB",
      TRUE ~ NA_character_
    ),
    Tissue = case_when(
      str_detect(Sample, "^Serum") ~ "Serum",
      str_detect(Sample, "^MC")    ~ "Cecum",
      str_detect(Sample, "^JC")    ~ "Colon",
      str_detect(Sample, "^ZC")    ~ "Rectum",
      TRUE ~ NA_character_
    ),
    Intensity = as.numeric(Intensity)
  ) %>%
  filter(!is.na(Treatment), !is.na(Tissue)) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("WT", "DB")),
    Tissue    = factor(Tissue, levels = c("Cecum", "Colon", "Rectum", "Serum")),
    Group     = paste0(as.character(Tissue), "-", as.character(Treatment)),
    Group     = factor(Group, levels = group_levels_all)
  )

# ======================================================
# 4) Robust metabolite name matching
# ======================================================

norm_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[\u00A0]", " ") %>%
    str_trim() %>%
    str_squish() %>%
    str_to_lower()
}

target_map <- tibble(
  Target_name = all_target_metabolites,
  Metabolite_norm = norm_name(all_target_metabolites)
)

dat_target <- dat_long %>%
  mutate(Metabolite_norm = norm_name(Metabolite)) %>%
  inner_join(target_map, by = "Metabolite_norm") %>%
  mutate(Metabolite = Target_name) %>%
  select(-Target_name, -Metabolite_norm) %>%
  group_by(Tissue, Sample, Treatment, Group, Metabolite) %>%
  summarise(
    Intensity = max(Intensity, na.rm = TRUE),
    .groups = "drop"
  )

missing_targets <- setdiff(
  all_target_metabolites,
  unique(dat_target$Metabolite)
)

if (length(missing_targets) > 0) {
  warning(
    "The following target metabolites were not matched: ",
    paste(missing_targets, collapse = ", ")
  )
}

target_match_summary <- tibble(
  Target_metabolite = all_target_metabolites,
  Matched = all_target_metabolites %in% unique(dat_target$Metabolite)
)

write_csv(
  target_match_summary,
  file.path(output_dir, "Fig4e_4g_4h_6j_target_match_summary.csv")
)

# ======================================================
# 5) Helper functions
# ======================================================

safe_div <- function(num, den) {
  num <- as.numeric(num)
  den <- as.numeric(den)
  out <- num / den
  out[!is.finite(out) | !is.finite(num) | !is.finite(den) | den <= 0] <- NA_real_
  out
}

row_sum_na <- function(...) {
  mat <- cbind(...)
  out <- rowSums(mat, na.rm = TRUE)
  all_na <- rowSums(is.finite(mat)) == 0
  out[all_na] <- NA_real_
  out
}

ensure_cols <- function(df, cols) {
  for (cc in cols) {
    if (!cc %in% colnames(df)) {
      df[[cc]] <- NA_real_
    }
  }
  df
}

sanitize <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

pad_limits <- function(y, ypad_ratio = 0.20, min_span_abs = 1e-6) {
  y <- y[is.finite(y)]
  if (!length(y)) return(c(0, 1))
  
  ymin <- min(y, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)
  span <- ymax - ymin
  
  if (!is.finite(span) || span <= 0) {
    ymin2 <- ymin - 0.5
    ymax2 <- ymax + 0.5
    return(c(ymin2, ymax2 * (1 + ypad_ratio)))
  } else {
    span_use <- max(span, min_span_abs)
    ymax2 <- ymax + span_use * ypad_ratio
    return(c(ymin, ymax2))
  }
}

safe_welch_stats <- function(df) {
  x_db <- df %>%
    filter(Treatment == "DB") %>%
    pull(Value)
  
  x_wt <- df %>%
    filter(Treatment == "WT") %>%
    pull(Value)
  
  x_db <- x_db[is.finite(x_db)]
  x_wt <- x_wt[is.finite(x_wt)]
  
  n_db <- length(x_db)
  n_wt <- length(x_wt)
  
  base_out <- tibble(
    n_WT = n_wt,
    n_DB = n_db,
    mean_WT = mean(x_wt, na.rm = TRUE),
    mean_DB = mean(x_db, na.rm = TRUE),
    median_WT = median(x_wt, na.rm = TRUE),
    median_DB = median(x_db, na.rm = TRUE),
    sd_WT = sd(x_wt, na.rm = TRUE),
    sd_DB = sd(x_db, na.rm = TRUE)
  )
  
  if (n_wt < 2 || n_db < 2) {
    return(
      base_out %>%
        mutate(
          mean_difference_DB_minus_WT = NA_real_,
          t_statistic = NA_real_,
          df = NA_real_,
          p_value = NA_real_,
          conf_low_95 = NA_real_,
          conf_high_95 = NA_real_,
          cohen_d_DB_vs_WT = NA_real_
        )
    )
  }
  
  tt <- suppressWarnings(
    t.test(
      x_db,
      x_wt,
      var.equal = FALSE,
      alternative = "two.sided"
    )
  )
  
  pooled_sd <- sqrt(
    ((n_db - 1) * var(x_db) + (n_wt - 1) * var(x_wt)) /
      (n_db + n_wt - 2)
  )
  
  cohen_d <- ifelse(
    is.finite(pooled_sd) && pooled_sd > 0,
    (mean(x_db) - mean(x_wt)) / pooled_sd,
    NA_real_
  )
  
  base_out %>%
    mutate(
      mean_difference_DB_minus_WT = mean(x_db) - mean(x_wt),
      t_statistic = unname(tt$statistic),
      df = unname(tt$parameter),
      p_value = tt$p.value,
      conf_low_95 = tt$conf.int[1],
      conf_high_95 = tt$conf.int[2],
      cohen_d_DB_vs_WT = cohen_d
    )
}

make_stats <- function(source_df) {
  source_df %>%
    group_by(Panel, Feature, Tissue) %>%
    group_modify(~ safe_welch_stats(.x)) %>%
    ungroup() %>%
    group_by(Panel) %>%
    mutate(BH_adjusted_P_value_within_panel = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    group_by(Panel, Tissue) %>%
    mutate(BH_adjusted_P_value_within_tissue = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      Test = "Two-sided Welch's t-test",
      Comparison = "DB vs WT",
      Unit_of_study = "biologically independent male mouse"
    )
}

# ======================================================
# 6) Plot theme
# ======================================================

nature_theme <- theme_classic(base_size = 16) +
  theme(
    axis.title.y = element_text(size = 17, face = "bold"),
    axis.text.x  = element_text(size = 17, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 17),
    strip.text   = element_text(size = 20, face = "bold"),
    legend.position = "none",
    plot.title   = element_text(size = 21),
    plot.margin  = margin(t = 8, r = 35, b = 8, l = 35, unit = "mm")
  )

col_fill  <- c("WT" = "#444444", "DB" = "#D62728")
col_point <- c("WT" = "#222222", "DB" = "#D62728")
shape_group <- c("WT" = 16, "DB" = 17)

plot_each_feature_all8 <- function(source_df,
                                   y_label = "Value",
                                   file_prefix,
                                   width = 9.6,
                                   height = 6) {
  
  all_features <- unique(source_df$Feature)
  
  for (ff in all_features) {
    df_plot <- source_df %>%
      filter(Feature == ff) %>%
      mutate(
        Group = factor(Group, levels = group_levels_all),
        Treatment = factor(Treatment, levels = c("WT", "DB"))
      )
    
    if (nrow(df_plot) == 0) next
    
    present_groups <- df_plot %>%
      filter(is.finite(Value)) %>%
      distinct(Group) %>%
      pull(Group) %>%
      as.character()
    
    cmp_list <- compare_list_all[
      sapply(compare_list_all, function(z) all(z %in% present_groups))
    ]
    
    y_limits <- pad_limits(df_plot$Value)
    
    p <- ggplot(df_plot, aes(x = Group, y = Value, fill = Treatment)) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.6,
        color = "black",
        alpha = 0.7,
        size = 0.6
      ) +
      geom_jitter(
        aes(color = Treatment, shape = Treatment),
        width = 0.15,
        size = 1.6,
        alpha = 0.75
      ) +
      scale_x_discrete(limits = group_levels_all, drop = FALSE) +
      scale_fill_manual(values = col_fill, drop = FALSE) +
      scale_color_manual(values = col_point, drop = FALSE) +
      scale_shape_manual(values = shape_group, drop = FALSE) +
      coord_cartesian(ylim = y_limits) +
      nature_theme +
      labs(y = y_label, x = NULL, title = ff) +
      {
        if (length(cmp_list) > 0)
          stat_compare_means(
            comparisons = cmp_list,
            method = "t.test",
            method.args = list(var.equal = FALSE, alternative = "two.sided"),
            label = "p.signif",
            hide.ns = FALSE,
            label.y.npc = 0.96,
            step.increase = 0,
            tip.length = 0.02,
            bracket.size = 0.5,
            size = 6
          )
        else NULL
      }
    
    fn <- paste0(file_prefix, "_", sanitize(ff))
    
    ggsave(
      file.path(output_dir, paste0(fn, ".pdf")),
      p,
      width = width,
      height = height,
      dpi = 300
    )
    
    ggsave(
      file.path(output_dir, paste0(fn, ".png")),
      p,
      width = width,
      height = height,
      dpi = 300
    )
  }
}

# ======================================================
# 7) Fig. 4e source data and plots
#    Gln/Glu and NAG/Glu
# ======================================================

wide_4e <- dat_target %>%
  filter(Metabolite %in% fig4e_component_metabolites) %>%
  select(Tissue, Sample, Treatment, Group, Metabolite, Intensity) %>%
  pivot_wider(
    names_from = Metabolite,
    values_from = Intensity,
    values_fill = NA_real_
  )

wide_4e <- ensure_cols(wide_4e, fig4e_component_metabolites)

fig4e_source <- wide_4e %>%
  mutate(
    `Gln/Glu` = safe_div(`Glutamine`, `Glutamic acid`),
    `NAG/Glu` = safe_div(`N-Acetyl-L-glutamic acid`, `Glutamic acid`)
  ) %>%
  pivot_longer(
    cols = all_of(fig4e_ratio_targets),
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  filter(is.finite(Value)) %>%
  mutate(
    Panel = "Fig4e",
    Feature_class = "Ratio"
  ) %>%
  select(Panel, Feature_class, Feature, Tissue, Sample, Treatment, Group, Value)

plot_each_feature_all8(
  source_df = fig4e_source,
  y_label = "Relative abundance ratio",
  file_prefix = "Fig4e",
  width = 9.6,
  height = 6
)

# ======================================================
# 8) Fig. 4g source data and plots
#    Even LCFA and Odd-chain FA
#    X-axis also includes serum
# ======================================================

wide_4g <- dat_target %>%
  filter(Metabolite %in% fig4g_component_metabolites) %>%
  select(Tissue, Sample, Treatment, Group, Metabolite, Intensity) %>%
  pivot_wider(
    names_from = Metabolite,
    values_from = Intensity,
    values_fill = NA_real_
  )

wide_4g <- ensure_cols(wide_4g, fig4g_component_metabolites)

fig4g_source <- wide_4g %>%
  mutate(
    `Even LCFA (C16:0+C16:1+C18:1)` = row_sum_na(
      `Palmitic acid`,
      `Palmitoleic acid`,
      `Oleic acid`
    ),
    `Odd-chain FA (C15:0+C17:0)` = row_sum_na(
      `Pentadecanoic acid`,
      `Heptadecanoic acid`
    )
  ) %>%
  pivot_longer(
    cols = c(
      `Even LCFA (C16:0+C16:1+C18:1)`,
      `Odd-chain FA (C15:0+C17:0)`
    ),
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  filter(is.finite(Value)) %>%
  mutate(
    Panel = "Fig4g",
    Feature_class = case_when(
      Feature == "Even LCFA (C16:0+C16:1+C18:1)" ~
        "Even-chain long-chain fatty acids",
      Feature == "Odd-chain FA (C15:0+C17:0)" ~
        "Odd-chain fatty acids",
      TRUE ~ "Fatty acids"
    )
  ) %>%
  select(Panel, Feature_class, Feature, Tissue, Sample, Treatment, Group, Value)

plot_each_feature_all8(
  source_df = fig4g_source,
  y_label = "Relative abundance",
  file_prefix = "Fig4g",
  width = 9.6,
  height = 6
)

# ======================================================
# 9) Fig. 4h source data and plots
#    Taurodeoxycholic acid and taurochenodeoxycholic acid
# ======================================================

fig4h_source <- dat_target %>%
  filter(Metabolite %in% fig4h_metabolites) %>%
  transmute(
    Panel = "Fig4h",
    Feature_class = "Conjugated bile acid",
    Feature = Metabolite,
    Tissue,
    Sample,
    Treatment,
    Group,
    Value = Intensity
  ) %>%
  filter(is.finite(Value))

plot_each_feature_all8(
  source_df = fig4h_source,
  y_label = "Relative abundance",
  file_prefix = "Fig4h",
  width = 9.6,
  height = 6
)

# ======================================================
# 10) Fig. 6j source data and plots
#     Polyamine intermediates
# ======================================================

fig6j_source <- dat_target %>%
  filter(Metabolite %in% fig6j_metabolites) %>%
  transmute(
    Panel = "Fig6j",
    Feature_class = "Polyamine intermediate",
    Feature = Metabolite,
    Tissue,
    Sample,
    Treatment,
    Group,
    Value = Intensity
  ) %>%
  filter(is.finite(Value))

plot_each_feature_all8(
  source_df = fig6j_source,
  y_label = "Relative abundance",
  file_prefix = "Fig6j",
  width = 9.6,
  height = 6
)

# ======================================================
# 11) Combined Source Data table with Welch's statistics
#     Wide format:
#     - One row per metabolite / ratio / metric
#     - Sample values are separate columns
#     - Feature is unique and used as row identifier
#     - Other metadata and statistics are saved as columns
# ======================================================

all_source <- bind_rows(
  fig4e_source,
  fig4g_source,
  fig4h_source,
  fig6j_source
)

all_stats <- make_stats(all_source)

write_csv(
  all_source,
  source_out_file
)

write_csv(
  all_stats,
  file.path(output_dir, "Fig4e_4g_4h_6j_Welch_stats_long_ALL8.csv")
)

# ----------------------------
# 11.1 Sample values as separate columns
# ----------------------------

sample_value_wide <- all_source %>%
  mutate(
    Feature = as.character(Feature),
    Sample_col = paste0(
      as.character(Tissue), "_",
      as.character(Treatment), "_",
      make.names(as.character(Sample))
    )
  ) %>%
  select(
    Feature,
    Panel,
    Feature_class,
    Sample_col,
    Value
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = Sample_col,
    values_from = Value
  )

# ----------------------------
# 11.2 Statistics as columns
#      Statistics are tissue-specific, so each tissue gets its own columns
# ----------------------------

stats_wide <- all_stats %>%
  mutate(
    Feature = as.character(Feature),
    Tissue = as.character(Tissue)
  ) %>%
  select(
    Feature,
    Tissue,
    Test,
    Comparison,
    Unit_of_study,
    n_WT,
    n_DB,
    mean_WT,
    mean_DB,
    median_WT,
    median_DB,
    sd_WT,
    sd_DB,
    mean_difference_DB_minus_WT,
    t_statistic,
    df,
    p_value,
    BH_adjusted_P_value_within_panel,
    BH_adjusted_P_value_within_tissue,
    conf_low_95,
    conf_high_95,
    cohen_d_DB_vs_WT
  ) %>%
  pivot_wider(
    id_cols = Feature,
    names_from = Tissue,
    values_from = c(
      n_WT,
      n_DB,
      mean_WT,
      mean_DB,
      median_WT,
      median_DB,
      sd_WT,
      sd_DB,
      mean_difference_DB_minus_WT,
      t_statistic,
      df,
      p_value,
      BH_adjusted_P_value_within_panel,
      BH_adjusted_P_value_within_tissue,
      conf_low_95,
      conf_high_95,
      cohen_d_DB_vs_WT
    ),
    names_glue = "{.value}_{Tissue}"
  ) %>%
  left_join(
    all_stats %>%
      distinct(Feature, Test, Comparison, Unit_of_study),
    by = "Feature"
  )

# ----------------------------
# 11.3 Combine sample values + statistics
# ----------------------------

source_data_wide <- sample_value_wide %>%
  left_join(stats_wide, by = "Feature") %>%
  relocate(
    Feature,
    Panel,
    Feature_class,
    Test,
    Comparison,
    Unit_of_study
  ) %>%
  arrange(Panel, Feature)

# ----------------------------
# 11.4 Save wide Source Data table
# ----------------------------

source_out_file_wide <- file.path(
  output_dir,
  "Fig4e_4g_4h_6j_SourceData_with_Welch_stats_ALL8_wide.csv"
)

write_csv(
  source_data_wide,
  source_out_file_wide
)
