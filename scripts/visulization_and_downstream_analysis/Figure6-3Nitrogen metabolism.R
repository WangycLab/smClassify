# This script performs metabolite-level differential analysis across serum
# and intestinal segments (Cecum, Colon, Rectum) in WT and DB mice.

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggpubr)
})

## ---------- 0) Basic settings ----------
# Root project directory (edit if project path changes)
project_dir <- "D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009"

# Working directory for this figure
setwd(file.path(project_dir, "Figure6"))

# Directory for all plots
dir.create("plots", showWarnings = FALSE)
plot_dir <- file.path(getwd(), "plots")

# ====== Targets (canonical names & de-duplication) ======
targets_raw <- c(
  # Urea cycle core
  "Urea","Citrulline","Ornithine","Arginine","Glutamine","Glutamic acid","D-Glutamic acid",
  "N-Acetyl-L-glutamic acid","Homo-L-arginine","DL-Arginine","L-Phosphoarginine",
  
  # Uric acid branch
  "Uric acid","1,7-Dimethyluric acid","7-Methyluric acid",
  "Hypoxanthine","Xanthine","1-Methylxanthine","Allantoin",
  
  # Polyamine branch
  "N-Acetylputrescine","N8-Acetylspermidine","N8-Acetylspermine","N1,N12-Diacetylspermine","Sinapoylspermine",
  "N-Acetylcadaverine",
  
  # Related conjugates / derivatives
  "Formiminoglutamic acid","N-Phenylacetylglutamic acid","N5-Acetylornithine","N2-Acetylornithine",
  "Creatine","Creatinine","Phosphocreatine",
  "Asymmetric dimethylarginine","Symmetric dimethylarginine",
  "Urocanic acid",
  "N-Acetylputrescine","N1,N12-Diacetylspermine","Sinapoylspermine","N-Acetylcadaverine"
)

# Case-insensitive de-duplication, keeping the first canonical spelling
targets_lc <- tolower(targets_raw)
targets     <- targets_raw[match(unique(targets_lc), targets_lc)]
targets

# Filename prefix (change to short code like "AA" if necessary)
fname_prefix <- "Nitrogen metabolism"

## ---------- 1) Read input tables & pivot to long format ----------
serum_file   <- "serum_norm_ms2_metabolites_intensity.xlsx"
content_file <- "content_norm_ms2_metabolites_intensity.xlsx"

serum   <- read_xlsx(serum_file)
content <- read_xlsx(content_file)

# Fix typo in column names if present
names(serum) <- str_replace_all(names(serum), "^Seurm", "Serum")

# Feature columns that are not sample intensities
feature_cols <- c(
  "Type","ID","MZ","RT","MS2.name","MS2_score",
  "Formula","SuperClass","Class","Subclass",
  "HMDB","kegg","kegg_pathway","MS1.name"
)

serum_long   <- serum   |> pivot_longer(-all_of(feature_cols), names_to = "Sample", values_to = "Intensity")
content_long <- content |> pivot_longer(-all_of(feature_cols), names_to = "Sample", values_to = "Intensity")

# Merge serum + gut into one long table
all_long <- bind_rows(serum_long, content_long)

## ---------- 2) Define groups (WT vs DB; Serum / Cecum / Colon / Rectum) ----------
dat_long <- all_long %>%
  dplyr::rename(Metabolite = MS2.name) %>%
  mutate(
    # Use MS1.name as backup if MS2.name is empty
    Metabolite = ifelse(is.na(Metabolite) | !nzchar(Metabolite), MS1.name, Metabolite),
    
    # Treatment group (WT vs DB/T2DM)
    Treatment  = case_when(
      str_detect(Sample, "WT")        ~ "WT",
      str_detect(Sample, "DB|T2DM")   ~ "DB",
      TRUE                            ~ NA_character_
    ),
    
    # Tissue, based on sample naming pattern
    Tissue = case_when(
      str_detect(Sample, "^Serum") ~ "Serum",
      str_detect(Sample, "^MC")    ~ "Cecum",
      str_detect(Sample, "^JC")    ~ "Colon",
      str_detect(Sample, "^ZC")    ~ "Rectum",
      TRUE                         ~ "Other"
    )
  ) %>%
  # Keep only WT/DB and known tissues
  filter(!is.na(Treatment), Tissue != "Other") %>%
  mutate(
    Treatment = factor(Treatment, levels = c("WT","DB")),
    Tissue    = factor(Tissue,    levels = c("Serum","Cecum","Colon","Rectum"))
  )

## ---------- 3) Keep only target metabolites (case-insensitive; normalized names) ----------
targets_regex <- paste0(
  "(",
  paste(stringr::str_replace_all(targets, "([\\W])","\\\\\\1"), collapse = "|"),
  ")"
)

dat_target <- dat_long %>%
  filter(str_detect(Metabolite, regex(targets_regex, ignore_case = TRUE))) %>%
  rowwise() %>%
  mutate(
    # Map to canonical name in `targets` (case-insensitive)
    Metabolite_std = {
      i <- match(tolower(Metabolite), tolower(targets))
      if (!is.na(i)) targets[i] else Metabolite
    }
  ) %>%
  ungroup() %>%
  mutate(Metabolite = Metabolite_std) %>%
  select(-Metabolite_std)

stopifnot(nrow(dat_target) > 0)

## ---------- 4) Welch's t-test (two-sided, safe wrapper) ----------
safe_t_p <- function(x, y, var_equal = FALSE, alternative = "two.sided") {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(
    t.test(y, x, var.equal = var_equal, alternative = alternative)$p.value
  )
}

stats_all <- dat_target %>%
  group_by(Tissue, Metabolite) %>%
  summarise(
    n_WT    = sum(Treatment == "WT" & is.finite(Intensity)),
    n_DB    = sum(Treatment == "DB" & is.finite(Intensity)),
    mean_WT = mean(Intensity[Treatment == "WT"], na.rm = TRUE),
    mean_DB = mean(Intensity[Treatment == "DB"], na.rm = TRUE),
    sd_WT   = sd(Intensity[Treatment == "WT"], na.rm = TRUE),
    sd_DB   = sd(Intensity[Treatment == "DB"], na.rm = TRUE),
    p_val   = safe_t_p(
      Intensity[Treatment == "WT"],
      Intensity[Treatment == "DB"],
      var_equal   = FALSE,
      alternative = "two.sided"
    ),
    .groups = "drop"
  ) %>%
  group_by(Tissue) %>%
  mutate(FDR = p.adjust(p_val, method = "fdr")) %>%
  ungroup()

write.csv(
  stats_all,
  file.path(plot_dir, paste0(fname_prefix, "_WT_vs_DB_stats_ttest_all.csv")),
  row.names = FALSE
)

## ---------- 5) Source data tables + unified group labels ----------
dat_target <- dat_target %>%
  mutate(Group = paste0(as.character(Tissue), "-", as.character(Treatment)))

# Serum source data (point-level)
source_serum <- dat_target %>%
  filter(Tissue == "Serum") %>%
  select(Metabolite, Tissue, Sample, Treatment, Group, Intensity)

# Gut source data (Cecum / Colon / Rectum; point-level)
source_gut <- dat_target %>%
  filter(Tissue %in% c("Cecum","Colon","Rectum")) %>%
  select(Metabolite, Tissue, Sample, Treatment, Group, Intensity)

write.csv(
  source_serum,
  file.path(plot_dir, paste0(fname_prefix, "_SourceData_SERUM_points.csv")),
  row.names = FALSE
)
write.csv(
  source_gut,
  file.path(plot_dir, paste0(fname_prefix, "_SourceData_GUT_points.csv")),
  row.names = FALSE
)

## ---------- 5.1) Per-metabolite y-axis limits (with 20% padding on top) ----------
ypad_ratio   <- 0.20     # fraction of the data span to add at the top
min_span_abs <- 1e-6     # minimal span to avoid zero-range issues

pad_limits <- function(y) {
  y <- y[is.finite(y)]
  if (!length(y)) return(c(0, 1))
  ymin <- min(y, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)
  span <- ymax - ymin
  
  # If no span (flat data), create a small artificial range
  if (!is.finite(span) || span <= 0) {
    ymin2 <- ymin - 0.5
    ymax2 <- ymax + 0.5
    return(c(ymin2, ymax2 * (1 + ypad_ratio)))
  } else {
    span_use <- max(span, min_span_abs)
    ymax2    <- ymax + span_use * ypad_ratio
    return(c(ymin, ymax2))
  }
}

# Serum y-limits per metabolite
serum_ylim <- source_serum %>%
  group_by(Metabolite) %>%
  summarise(
    ymin = min(Intensity, na.rm = TRUE),
    ymax = max(Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    ymin_pad = pad_limits(c(ymin, ymax))[1],
    ymax_pad = pad_limits(c(ymin, ymax))[2]
  ) %>%
  ungroup()

# Gut y-limits per metabolite
gut_ylim <- source_gut %>%
  group_by(Metabolite) %>%
  summarise(
    ymin = min(Intensity, na.rm = TRUE),
    ymax = max(Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    ymin_pad = pad_limits(c(ymin, ymax))[1],
    ymax_pad = pad_limits(c(ymin, ymax))[2]
  ) %>%
  ungroup()

write.csv(
  serum_ylim,
  file.path(plot_dir, paste0(fname_prefix, "_SourceData_SERUM_ylim_by_metabolite.csv")),
  row.names = FALSE
)
write.csv(
  gut_ylim,
  file.path(plot_dir, paste0(fname_prefix, "_SourceData_GUT_ylim_by_metabolite.csv")),
  row.names = FALSE
)

## ---------- 6) Plot theme (with sufficient margins for p-value brackets) ----------
nature_theme <- theme_classic(base_size = 16) +
  theme(
    axis.title.y   = element_text(size = 17, face = "bold"),
    axis.text.x    = element_text(size = 19, angle = 45, hjust = 1, vjust = 1),
    axis.text.y    = element_text(size = 17),
    strip.text     = element_text(size = 20, face = "bold"),
    legend.position = "none",
    plot.title     = element_text(size = 21),
    plot.margin    = margin(t = 8, r = 35, b = 8, l = 35, unit = "mm")
  )

col_fill  <- c("WT" = "#444444", "DB" = "#D62728")
col_point <- c("WT" = "#222222", "DB" = "#D62728")

# Utility to sanitize metabolite names for filenames
sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

all_mets <- unique(dat_target$Metabolite)

# Fixed group order for Serum and Gut
serum_levels <- c("Serum-WT","Serum-DB")
gut_levels   <- c("Cecum-WT","Cecum-DB","Colon-WT","Colon-DB","Rectum-WT","Rectum-DB")

## ---------- 7) Plotting ----------
# ---- Serum: one plot per metabolite ----
for (m in all_mets) {
  df_serum <- source_serum %>% filter(Metabolite == m)
  if (nrow(df_serum) == 0) next
  
  df_serum <- df_serum %>%
    mutate(
      Group     = factor(Group, levels = serum_levels),
      Treatment = factor(Treatment, levels = c("WT","DB"))
    )
  
  # Per-metabolite y-limits
  lim_ser  <- serum_ylim %>% filter(Metabolite == m)
  y_limits <- c(lim_ser$ymin_pad, lim_ser$ymax_pad)
  
  p_serum <- ggplot(df_serum, aes(x = Group, y = Intensity, fill = Treatment)) +
    geom_boxplot(
      outlier.shape = NA,
      width         = 0.6,
      color         = "black",
      alpha         = 0.7,
      size          = 0.6
    ) +
    geom_jitter(
      aes(color = Treatment),
      width = 0.12,
      size  = 1.6,
      alpha = 0.6
    ) +
    scale_fill_manual(values = col_fill) +
    scale_color_manual(values = col_point) +
    coord_cartesian(ylim = y_limits) +
    nature_theme +
    labs(y = "Intensity", x = NULL, title = m) +
    # Single comparison: Serum-WT vs Serum-DB
    stat_compare_means(
      comparisons  = list(c("Serum-WT","Serum-DB")),
      method       = "t.test",
      method.args  = list(var.equal = FALSE, alternative = "two.sided"),
      label        = "p.signif",
      hide.ns      = FALSE,
      label.y.npc  = 0.98,   # Fixed relative height close to top of panel
      tip.length   = 0.02,
      bracket.size = 0.5,
      size         = 6
    )
  
  fn <- paste0("Nature_", fname_prefix, "_SERUM_", sanitize(m))
  ggsave(file.path(plot_dir, paste0(fn, ".pdf")), p_serum, width = 5.4, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, paste0(fn, ".png")), p_serum, width = 5.4, height = 6, dpi = 300)
}

# ---- Gut: one plot per metabolite (same p-label height across three gut sites) ----
for (m in all_mets) {
  df_gut <- source_gut %>% filter(Metabolite == m)
  if (nrow(df_gut) == 0) next
  
  # Keep only gut groups that actually exist in this metabolite
  present_levels <- gut_levels[gut_levels %in% df_gut$Group]
  
  df_gut <- df_gut %>%
    mutate(
      Group     = factor(Group, levels = present_levels),
      Treatment = factor(Treatment, levels = c("WT","DB"))
    )
  
  # Per-metabolite y-limits
  lim_gut  <- gut_ylim %>% filter(Metabolite == m)
  y_limits <- c(lim_gut$ymin_pad, lim_gut$ymax_pad)
  
  # For each gut site, compare WT vs DB only if both exist
  tissues_in_df <- c("Cecum","Colon","Rectum")
  cmp_list <- list()
  for (tiss in tissues_in_df) {
    g1 <- paste0(tiss, "-WT")
    g2 <- paste0(tiss, "-DB")
    if (all(c(g1, g2) %in% present_levels)) {
      cmp_list <- append(cmp_list, list(c(g1, g2)))
    }
  }
  
  # Figure width scales with number of groups
  width_n <- max(6, 1.2 * length(present_levels))
  
  p_gut <- ggplot(df_gut, aes(x = Group, y = Intensity, fill = Treatment)) +
    geom_boxplot(
      outlier.shape = NA,
      width         = 0.6,
      color         = "black",
      alpha         = 0.7,
      size          = 0.6
    ) +
    geom_jitter(
      aes(color = Treatment),
      width = 0.15,
      size  = 1.6,
      alpha = 0.6
    ) +
    scale_fill_manual(values = col_fill) +
    scale_color_manual(values = col_point) +
    coord_cartesian(ylim = y_limits) +
    nature_theme +
    labs(y = "Intensity", x = NULL, title = m) +
    {
      # Add p-value annotations only if we have comparisons
      if (length(cmp_list) > 0)
        stat_compare_means(
          comparisons   = cmp_list,
          method        = "t.test",
          method.args   = list(var.equal = FALSE, alternative = "two.sided"),
          label         = "p.signif",
          hide.ns       = FALSE,
          label.y.npc   = 0.96,  # All brackets at the same relative height
          step.increase = 0,     # Do not stagger brackets upward
          tip.length    = 0.02,
          bracket.size  = 0.5,
          size          = 6
        ) else NULL
    }
  
  fn <- paste0("Nature_", fname_prefix, "_GUT_", sanitize(m))
  ggsave(file.path(plot_dir, paste0(fn, ".pdf")), p_gut, width = width_n, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, paste0(fn, ".png")), p_gut, width = width_n, height = 6, dpi = 300)
}

