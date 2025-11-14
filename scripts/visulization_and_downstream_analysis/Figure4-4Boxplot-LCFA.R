# ======================================================
# LCFA metabolism - Per-metabolite boxplots (WT vs DB)
# - Welch's t-test (two-sided)
# - Serum: x = "Serum-WT", "Serum-DB"
# - Gut:   x = "Cecum-WT","Cecum-DB","Colon-WT","Colon-DB","Rectum-WT","Rectum-DB"
# - Export stats + source data + per-metabolite ylim
# - y-axis auto padding (+20% top) to avoid clipping p-value labels
# - Gut significance labels at the SAME height (no step increase)
# ======================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggpubr)
})

## 0) Basic settings
setwd("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009/Figure4")
dir.create("plots", showWarnings = FALSE)
plot_dir <- file.path(getwd(), "plots")

#  Targets (standardize & deduplicate)
targets_raw <- c(
  # ---- SFA ----
  "Palmitic acid",
  "Stearic acid",
  "Pentadecanoic acid",
  "Heptadecanoic acid",
  
  # ---- MUFA ----
  "Palmitoleic acid",
  "Oleic acid",
  
  # ---- PUFA n-6 ----
  "Linoleic acid",
  "Gamma-Linolenic acid",
  "Dihomo-gamma-linolenic acid",
  "Eicosadienoic acid",
  "Arachidonic acid",
  
  # ---- PUFA n-3 ----
  "Eicosapentaenoic acid",
  "Docosapentaenoic acid",
  "Docosahexaenoic acid",
  
  "Conjugated linoleic acid",
  "9,10-EpOME",
  "12,13-EpOME",
  "9,10-DiHOME",
  "12,13-DiHOME",
  "5-HETE",
  "12-HETE",
  "15-HETE"
)

# Deduplicate (case-insensitive), keep the first occurrence as the standard name
targets_lc <- tolower(targets_raw)
targets     <- targets_raw[match(unique(targets_lc), targets_lc)]

# File name prefix 
fname_prefix <- "LCFA"

## 1) Read data & convert to long format
serum_file   <- "serum_norm_ms2_metabolites_intensity.xlsx"
content_file <- "content_norm_ms2_metabolites_intensity.xlsx"

serum   <- read_xlsx(serum_file)
content <- read_xlsx(content_file)
names(serum) <- str_replace_all(names(serum), "^Seurm", "Serum")

feature_cols <- c("Type","ID","MZ","RT","MS2.name","MS2_score",
                  "Formula","SuperClass","Class","Subclass",
                  "HMDB","kegg","kegg_pathway","MS1.name")

serum_long   <- serum   |> pivot_longer(-all_of(feature_cols), names_to = "Sample", values_to = "Intensity")
content_long <- content |> pivot_longer(-all_of(feature_cols), names_to = "Sample", values_to = "Intensity")
all_long <- bind_rows(serum_long, content_long)

## 2) Grouping (WT vs DB; Serum/Cecum/Colon/Rectum)
dat_long <- all_long %>%
  dplyr::rename(Metabolite = MS2.name) %>%
  mutate(
    Metabolite = ifelse(is.na(Metabolite) | !nzchar(Metabolite), MS1.name, Metabolite),
    Treatment  = case_when(
      str_detect(Sample, "WT") ~ "WT",
      str_detect(Sample, "DB|T2DM") ~ "DB",
      TRUE ~ NA_character_
    ),
    Tissue = case_when(
      str_detect(Sample, "^Serum") ~ "Serum",
      str_detect(Sample, "^MC")    ~ "Cecum",
      str_detect(Sample, "^JC")    ~ "Colon",
      str_detect(Sample, "^ZC")    ~ "Rectum",
      TRUE ~ "Other"
    )
  ) %>%
  filter(!is.na(Treatment), Tissue != "Other") %>%
  mutate(
    Treatment = factor(Treatment, levels = c("WT","DB")),
    Tissue    = factor(Tissue, levels = c("Serum","Cecum","Colon","Rectum"))
  )

## 3) Select only target metabolites (case-insensitive; standardized names)
targets_regex <- paste0("(", paste(stringr::str_replace_all(targets, "([\\W])","\\\\\\1"), collapse="|"), ")")

dat_target <- dat_long %>%
  filter(str_detect(Metabolite, regex(targets_regex, ignore_case = TRUE))) %>%
  rowwise() %>%
  mutate(Metabolite_std = {
    i <- match(tolower(Metabolite), tolower(targets))
    if (!is.na(i)) targets[i] else Metabolite
  }) %>%
  ungroup() %>%
  mutate(Metabolite = Metabolite_std) %>%
  select(-Metabolite_std)

stopifnot(nrow(dat_target) > 0)

## 4) Welch t-test (two-sided)
safe_t_p <- function(x, y, var_equal = FALSE, alternative = "two.sided") {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(t.test(y, x, var.equal = var_equal, alternative = alternative)$p.value)
}

stats_all <- dat_target %>%
  group_by(Tissue, Metabolite) %>%
  summarise(
    n_WT   = sum(Treatment=="WT" & is.finite(Intensity)),
    n_DB   = sum(Treatment=="DB" & is.finite(Intensity)),
    mean_WT= mean(Intensity[Treatment=="WT"], na.rm = TRUE),
    mean_DB= mean(Intensity[Treatment=="DB"], na.rm = TRUE),
    sd_WT  = sd(Intensity[Treatment=="WT"], na.rm = TRUE),
    sd_DB  = sd(Intensity[Treatment=="DB"], na.rm = TRUE),
    p_val  = safe_t_p(Intensity[Treatment=="WT"], Intensity[Treatment=="DB"],
                      var_equal = FALSE, alternative = "two.sided"),
    .groups = "drop"
  ) %>%
  group_by(Tissue) %>%
  mutate(FDR = p.adjust(p_val, method = "fdr")) %>%
  ungroup()

write.csv(stats_all, file.path(plot_dir, paste0(fname_prefix, "_WT_vs_DB_stats_ttest_all.csv")), row.names = FALSE)

## 5) Source Data + unified group labels
dat_target <- dat_target %>%
  mutate(Group = paste0(as.character(Tissue), "-", as.character(Treatment)))

source_serum <- dat_target %>%
  filter(Tissue == "Serum") %>%
  select(Metabolite, Tissue, Sample, Treatment, Group, Intensity)

source_gut <- dat_target %>%
  filter(Tissue %in% c("Cecum","Colon","Rectum")) %>%
  select(Metabolite, Tissue, Sample, Treatment, Group, Intensity)

write.csv(source_serum, file.path(plot_dir, paste0(fname_prefix, "_SourceData_SERUM_points.csv")), row.names = FALSE)
write.csv(source_gut,   file.path(plot_dir, paste0(fname_prefix, "_SourceData_GUT_points.csv")),   row.names = FALSE)

## 5.1) y-axis limits (including 20% top padding)
ypad_ratio   <- 0.20
min_span_abs <- 1e-6

pad_limits <- function(y) {
  y <- y[is.finite(y)]
  if (!length(y)) return(c(0, 1))
  ymin <- min(y, na.rm = TRUE); ymax <- max(y, na.rm = TRUE)
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

serum_ylim <- source_serum %>%
  group_by(Metabolite) %>%
  summarise(
    ymin = min(Intensity, na.rm = TRUE),
    ymax = max(Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(ymin_pad = pad_limits(c(ymin, ymax))[1],
         ymax_pad = pad_limits(c(ymin, ymax))[2]) %>%
  ungroup()

gut_ylim <- source_gut %>%
  group_by(Metabolite) %>%
  summarise(
    ymin = min(Intensity, na.rm = TRUE),
    ymax = max(Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(ymin_pad = pad_limits(c(ymin, ymax))[1],
         ymax_pad = pad_limits(c(ymin, ymax))[2]) %>%
  ungroup()

write.csv(serum_ylim, file.path(plot_dir, paste0(fname_prefix, "_SourceData_SERUM_ylim_by_metabolite.csv")), row.names = FALSE)
write.csv(gut_ylim,   file.path(plot_dir, paste0(fname_prefix, "_SourceData_GUT_ylim_by_metabolite.csv")),   row.names = FALSE)

## 6) Theme (sufficient left/right margin)
nature_theme <- theme_classic(base_size = 16) +
  theme(
    axis.title.y = element_text(size = 17, face = "bold"),
    axis.text.x  = element_text(size = 19, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 17),
    strip.text   = element_text(size = 20, face = "bold"),
    legend.position = "none",
    plot.title   = element_text(size = 21),
    plot.margin  = margin(t = 8, r = 35, b = 8, l = 35, unit = "mm")
  )
col_fill  <- c("WT"="#444444", "DB"="#D62728")
col_point <- c("WT"="#222222", "DB"="#D62728")

sanitize <- function(x) gsub("[^A-Za-z0-9]+","_", x)
all_mets <- unique(dat_target$Metabolite)

serum_levels <- c("Serum-WT","Serum-DB")
gut_levels   <- c("Cecum-WT","Cecum-DB","Colon-WT","Colon-DB","Rectum-WT","Rectum-DB")

## 7) Plotting
# Serum: single plot per metabolite 
for (m in all_mets) {
  df_serum <- source_serum %>% filter(Metabolite == m)
  if (nrow(df_serum) == 0) next
  
  df_serum <- df_serum %>%
    mutate(Group = factor(Group, levels = serum_levels),
           Treatment = factor(Treatment, levels = c("WT","DB")))
  
  lim_ser <- serum_ylim %>% filter(Metabolite == m)
  y_limits <- c(lim_ser$ymin_pad, lim_ser$ymax_pad)
  
  p_serum <- ggplot(df_serum, aes(x = Group, y = Intensity, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", alpha = 0.7, size = 0.6) +
    geom_jitter(aes(color = Treatment), width = 0.12, size = 1.6, alpha = 0.6) +
    scale_fill_manual(values = col_fill) +
    scale_color_manual(values = col_point) +
    coord_cartesian(ylim = y_limits) +
    nature_theme +
    labs(y = "Intensity", x = NULL, title = m) +
    stat_compare_means(
      comparisons = list(c("Serum-WT","Serum-DB")),
      method = "t.test",
      method.args = list(var.equal = FALSE, alternative = "two.sided"),
      label = "p.signif",
      hide.ns = FALSE,
      label.y.npc = 0.98,  
      tip.length = 0.02,
      bracket.size = 0.5,
      size = 6
    )
  
  fn <- paste0("Nature_", fname_prefix, "_SERUM_", sanitize(m))
  ggsave(file.path(plot_dir, paste0(fn, ".pdf")), p_serum, width = 5.4, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, paste0(fn, ".png")), p_serum, width = 5.4, height = 6, dpi = 300)
}

#  Gut: single plot per metabolite
for (m in all_mets) {
  df_gut <- source_gut %>% filter(Metabolite == m)
  if (nrow(df_gut) == 0) next
  
  present_levels <- gut_levels[gut_levels %in% df_gut$Group]
  df_gut <- df_gut %>%
    mutate(
      Group = factor(Group, levels = present_levels),
      Treatment = factor(Treatment, levels = c("WT","DB"))
    )
  
  lim_gut <- gut_ylim %>% filter(Metabolite == m)
  y_limits <- c(lim_gut$ymin_pad, lim_gut$ymax_pad)
  
  # Compare only WT vs DB per tissue (only if both groups are present)
  tissues_in_df <- c("Cecum","Colon","Rectum")
  cmp_list <- list()
  for (tiss in tissues_in_df) {
    g1 <- paste0(tiss, "-WT")
    g2 <- paste0(tiss, "-DB")
    if (all(c(g1, g2) %in% present_levels)) {
      cmp_list <- append(cmp_list, list(c(g1, g2)))
    }
  }
  
  width_n <- max(6, 1.2 * length(present_levels))
  
  p_gut <- ggplot(df_gut, aes(x = Group, y = Intensity, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", alpha = 0.7, size = 0.6) +
    geom_jitter(aes(color = Treatment), width = 0.15, size = 1.6, alpha = 0.6) +
    scale_fill_manual(values = col_fill) +
    scale_color_manual(values = col_point) +
    coord_cartesian(ylim = y_limits) +
    nature_theme +
    labs(y = "Intensity", x = NULL, title = m) +
    { if (length(cmp_list) > 0)
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
      ) else NULL }
  
  fn <- paste0("Nature_", fname_prefix, "_GUT_", sanitize(m))
  ggsave(file.path(plot_dir, paste0(fn, ".pdf")), p_gut, width = width_n, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, paste0(fn, ".png")), p_gut, width = width_n, height = 6, dpi = 300)
}

message("Done. Plots + source data saved to: ", plot_dir)



## 8)Calculate key ratios


# Vectorized division
v_div <- function(num, den) {
  num <- as.numeric(num); den <- as.numeric(den)
  out <- num / den
  bad <- !is.finite(out) | !is.finite(num) | !is.finite(den) | is.na(num) | is.na(den) | den <= 0
  out[bad] <- NA_real_
  out
}

# 8.1 Wide table (one row per sample; columns are metabolites)
wide_target <- dat_target %>%
  mutate(
    # Standardize DPA naming
    Metabolite = str_replace(
      Metabolite,
      regex("^Docosapentaenoic acid \\(22n-3\\)$", ignore_case = TRUE),
      "Docosapentaenoic acid"
    )
  ) %>%
  select(Tissue, Sample, Treatment, Metabolite, Intensity) %>%
  mutate(Metabolite = factor(Metabolite, levels = unique(Metabolite))) %>%
  tidyr::pivot_wider(names_from = Metabolite, values_from = Intensity, values_fill = NA)

# 8.2 Ensure required columns exist for calculation
ensure_cols <- function(df, cols) {
  for (cc in cols) if (!cc %in% names(df)) df[[cc]] <- NA_real_
  df
}
need_cols <- c(
  "Arachidonic acid",              # AA
  "Dihomo-gamma-linolenic acid",   # DGLA
  "Eicosapentaenoic acid",         # EPA
  "Docosapentaenoic acid"          # DPA
)
wide_target <- ensure_cols(wide_target, need_cols)

# 8.3 Calculate available ratios
wide_ratio <- wide_target %>%
  mutate(
    # n-6 conversion efficiency (5-desaturase direction)
    `AA/DGLA` = v_div(`Arachidonic acid`, `Dihomo-gamma-linolenic acid`),
    `DGLA/AA` = v_div(`Dihomo-gamma-linolenic acid`, `Arachidonic acid`),
    
    # n-3 / n-6 balance
    `EPA/AA`  = v_div(`Eicosapentaenoic acid`, `Arachidonic acid`),
    
    # n-3 metabolic flux (EPA DPA)
    `DPA/EPA` = v_div(`Docosapentaenoic acid`, `Eicosapentaenoic acid`)
  )


# 8.4 Convert to long format
ratio_long <- wide_ratio %>%
  pivot_longer(
    cols = c(`AA/DGLA`,`DGLA/AA`,`EPA/AA`,`DPA/EPA`),
    names_to = "Ratio", values_to = "Value"
  ) %>%
  filter(is.finite(Value)) %>%
  mutate(Group = paste0(as.character(Tissue), "-", as.character(Treatment)))

# 8.5 Statistics (Welch's t-test), FDR by tissue
ratio_stats <- ratio_long %>%
  group_by(Tissue, Ratio) %>%
  summarise(
    n_WT   = sum(Group %in% paste0(unique(Tissue), "-WT") & is.finite(Value)),
    n_DB   = sum(Group %in% paste0(unique(Tissue), "-DB") & is.finite(Value)),
    mean_WT= mean(Value[Group %in% paste0(unique(Tissue), "-WT")], na.rm = TRUE),
    mean_DB= mean(Value[Group %in% paste0(unique(Tissue), "-DB")], na.rm = TRUE),
    sd_WT  = sd(Value[Group %in% paste0(unique(Tissue), "-WT")], na.rm = TRUE),
    sd_DB  = sd(Value[Group %in% paste0(unique(Tissue), "-DB")], na.rm = TRUE),
    p_val  = safe_t_p(Value[Group %in% paste0(unique(Tissue), "-WT")],
                      Value[Group %in% paste0(unique(Tissue), "-DB")],
                      var_equal = FALSE, alternative = "two.sided"),
    .groups = "drop_last"
  ) %>%
  group_by(Tissue) %>%
  mutate(FDR = p.adjust(p_val, method = "fdr")) %>%
  ungroup()

write.csv(ratio_stats, file.path(plot_dir, paste0(fname_prefix, "_RATIOS_WT_vs_DB_stats_ttest_all.csv")), row.names = FALSE)
write.csv(ratio_long,  file.path(plot_dir, paste0(fname_prefix, "_RATIOS_SourceData_points.csv")),      row.names = FALSE)
# RATIO-specific color scheme (Okabe�CIto colorblind-friendly)
col_fill_ratio  <- c("WT" = "#1A1A1A",  
                     "DB" = "#E69F00")  
col_point_ratio <- c("WT" = "#1A1A1A",  
                     "DB" = "#E69F00") 
shape_ratio     <- c("WT" = 16, "DB" = 17)  

# 8.6 y-axis limits for ratios (including 20% top padding)
ratio_ylim <- ratio_long %>%
  group_by(Tissue, Ratio) %>%
  summarise(
    ymin = min(Value, na.rm = TRUE),
    ymax = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(ymin_pad = pad_limits(c(ymin, ymax))[1],
         ymax_pad = pad_limits(c(ymin, ymax))[2]) %>%
  ungroup()

write.csv(ratio_ylim, file.path(plot_dir, paste0(fname_prefix, "_RATIOS_ylim_by_ratio.csv")), row.names = FALSE)

# 8.7 Plotting LCFA ratios - Serum
all_ratios <- unique(ratio_long$Ratio)
for (r in all_ratios) {
  df_s <- ratio_long %>% filter(Tissue == "Serum", Ratio == r)
  if (nrow(df_s) == 0) next
  
  df_s <- df_s %>%
    mutate(
      Group = factor(Group, levels = c("Serum-WT","Serum-DB")),
      Treatment = factor(ifelse(str_detect(Group, "WT$"), "WT", "DB"),
                         levels = c("WT","DB"))
    )
  
  lim_s <- ratio_ylim %>% filter(Tissue == "Serum", Ratio == r)
  y_limits <- c(lim_s$ymin_pad, lim_s$ymax_pad)
  
  p_s <- ggplot(df_s, aes(x = Group, y = Value, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", alpha = 0.7, size = 0.6) +
    geom_jitter(aes(color = Treatment, shape = Treatment), width = 0.12, size = 1.6, alpha = 0.75) +
    scale_fill_manual(values = col_fill_ratio) +
    scale_color_manual(values = col_point_ratio) +
    scale_shape_manual(values = shape_ratio) +
    coord_cartesian(ylim = y_limits) +
    nature_theme +
    labs(y = r, x = NULL, title = r) +
    stat_compare_means(
      comparisons = list(c("Serum-WT","Serum-DB")),
      method = "t.test",
      method.args = list(var.equal = FALSE, alternative = "two.sided"),
      label = "p.signif",
      hide.ns = FALSE,
      label.y.npc = 0.98,
      tip.length = 0.02,
      bracket.size = 0.5,
      size = 6
    )
  
  fn <- paste0("Nature_", fname_prefix, "_RATIO_SERUM_", sanitize(r))
  ggsave(file.path(plot_dir, paste0(fn, ".pdf")), p_s, width = 5.4, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, paste0(fn, ".png")), p_s, width = 5.4, height = 6, dpi = 300)
}

# 8.8 Plotting LCFA ratios - Gut (Cecum/Colon/Rectum; equal height for significance brackets)
for (r in all_ratios) {
  df_g <- ratio_long %>% filter(Tissue %in% c("Cecum","Colon","Rectum"), Ratio == r)
  if (nrow(df_g) == 0) next
  
  present_levels <- gut_levels[gut_levels %in% df_g$Group]
  df_g <- df_g %>%
    mutate(
      Group = factor(Group, levels = present_levels),
      Treatment = factor(ifelse(str_detect(Group, "WT$"), "WT", "DB"),
                         levels = c("WT","DB"))
    )
  
  lim_g <- ratio_ylim %>% filter(Tissue %in% c("Cecum","Colon","Rectum"), Ratio == r)
  y_limits <- c(min(lim_g$ymin_pad, na.rm = TRUE), max(lim_g$ymax_pad, na.rm = TRUE))
  
  cmp_list <- list()
  for (tiss in c("Cecum","Colon","Rectum")) {
    g1 <- paste0(tiss, "-WT"); g2 <- paste0(tiss, "-DB")
    if (all(c(g1, g2) %in% present_levels)) cmp_list <- append(cmp_list, list(c(g1, g2)))
  }
  
  width_n <- max(6, 1.2 * length(present_levels))
  
  p_g <- ggplot(df_g, aes(x = Group, y = Value, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", alpha = 0.7, size = 0.6) +
    geom_jitter(aes(color = Treatment, shape = Treatment), width = 0.15, size = 1.6, alpha = 0.75) +
    scale_fill_manual(values = col_fill_ratio) +
    scale_color_manual(values = col_point_ratio) +
    scale_shape_manual(values = shape_ratio) +
    coord_cartesian(ylim = y_limits) +
    nature_theme +
    labs(y = r, x = NULL, title = r) +
    { if (length(cmp_list) > 0)
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
      ) else NULL }
  
  fn <- paste0("Nature_", fname_prefix, "_RATIO_GUT_", sanitize(r))
  ggsave(file.path(plot_dir, paste0(fn, ".pdf")), p_g, width = width_n, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, paste0(fn, ".png")), p_g, width = width_n, height = 6, dpi = 300)
}

message("Ratio analysis (simplified) done. Files saved to: ", plot_dir)
