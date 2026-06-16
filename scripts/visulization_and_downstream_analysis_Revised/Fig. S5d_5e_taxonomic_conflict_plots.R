############################################################
## Supplementary Figure 5d-e taxonomic conflict plots

############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(tibble)
  library(purrr)
})

set.seed(1234)

# ============================================================
# 1. Input and output paths
# ============================================================

base_dir   <- "Figure2"
output_dir <- file.path(base_dir, "Output")

taxonomy_rds <- file.path(
  output_dir,
  "Mouse_microbial_taxonomy_smClassify_Kraken2_MGnify.rds"
)

plot_dir <- file.path(output_dir, "SupplementaryFigure5d_5e_taxonomic_conflict")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(taxonomy_rds)) {
  stop("Input RDS not found: ", taxonomy_rds, call. = FALSE)
}

# ============================================================
# 2. Read annotated Seurat object
# ============================================================

sce_bac_annot <- readRDS(taxonomy_rds)
meta <- sce_bac_annot@meta.data %>%
  rownames_to_column("cell")

required_cols <- c(
  "smClassify_genus", "MGBC_genus", "Refseq_genus", "MGnify_genus",
  "smClassify_species", "MGBC_species", "Refseq_species", "MGnify_species"
)

missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# ============================================================
# 3. Helper functions
# ============================================================

.norm_taxon <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "unclassified", "Unclassified", "unknown", "Unknown")] <- NA
  x
}

pair_status <- function(a, b) {
  a <- .norm_taxon(a)
  b <- .norm_taxon(b)
  case_when(
    is.na(a) & is.na(b) ~ "both_na",
    !is.na(a) & is.na(b) ~ "only_first",
    is.na(a) & !is.na(b) ~ "only_second",
    a == b ~ "agree",
    TRUE ~ "conflict"
  )
}

make_pair <- function(a, b) {
  a <- .norm_taxon(a)
  b <- .norm_taxon(b)
  if (is.na(a) | is.na(b)) return(NA_character_)
  
  pair <- sort(c(a, b))
  paste(pair, collapse = " | ")
}

save_plot <- function(p, filename, width, height, dpi = 600) {
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = p, width = width, height = height,
    device = cairo_pdf, bg = "white"
  )
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".png")),
    plot = p, width = width, height = height,
    dpi = dpi, bg = "white"
  )
}

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}

# ============================================================
# 4. Build pairwise conflict source data
# ============================================================

conflict_df <- meta %>%
  transmute(
    cell = cell,
    orig.ident = orig.ident,
    loc = if ("loc" %in% colnames(meta)) loc else NA_character_,
    
    smClassify_genus   = .norm_taxon(smClassify_genus),
    MGBC_genus         = .norm_taxon(MGBC_genus),
    Refseq_genus       = .norm_taxon(Refseq_genus),
    MGnify_genus       = .norm_taxon(MGnify_genus),
    
    smClassify_species = .norm_taxon(smClassify_species),
    MGBC_species       = .norm_taxon(MGBC_species),
    Refseq_species     = .norm_taxon(Refseq_species),
    MGnify_species     = .norm_taxon(MGnify_species)
  ) %>%
  mutate(
    sm_vs_MGBC_genus      = pair_status(smClassify_genus, MGBC_genus),
    sm_vs_Refseq_genus    = pair_status(smClassify_genus, Refseq_genus),
    sm_vs_MGnify_genus    = pair_status(smClassify_genus, MGnify_genus),
    MGBC_vs_Refseq_genus  = pair_status(MGBC_genus, Refseq_genus),
    MGBC_vs_MGnify_genus  = pair_status(MGBC_genus, MGnify_genus),
    Refseq_vs_MGnify_genus = pair_status(Refseq_genus, MGnify_genus),
    
    sm_vs_MGBC_species      = pair_status(smClassify_species, MGBC_species),
    sm_vs_Refseq_species    = pair_status(smClassify_species, Refseq_species),
    sm_vs_MGnify_species    = pair_status(smClassify_species, MGnify_species),
    MGBC_vs_Refseq_species  = pair_status(MGBC_species, Refseq_species),
    MGBC_vs_MGnify_species  = pair_status(MGBC_species, MGnify_species),
    Refseq_vs_MGnify_species = pair_status(Refseq_species, MGnify_species)
  )

pair_order_genus <- c(
  "sm_vs_MGBC_genus",
  "sm_vs_Refseq_genus",
  "sm_vs_MGnify_genus",
  "MGBC_vs_Refseq_genus",
  "MGBC_vs_MGnify_genus",
  "Refseq_vs_MGnify_genus"
)

pair_order_species <- c(
  "sm_vs_MGBC_species",
  "sm_vs_Refseq_species",
  "sm_vs_MGnify_species",
  "MGBC_vs_Refseq_species",
  "MGBC_vs_MGnify_species",
  "Refseq_vs_MGnify_species"
)

pair_label_map <- c(
  "sm_vs_MGBC_genus" = "smClassify vs MGBC",
  "sm_vs_Refseq_genus" = "smClassify vs RefSeq",
  "sm_vs_MGnify_genus" = "smClassify vs MGnify",
  "MGBC_vs_Refseq_genus" = "MGBC vs RefSeq",
  "MGBC_vs_MGnify_genus" = "MGBC vs MGnify",
  "Refseq_vs_MGnify_genus" = "RefSeq vs MGnify",
  "sm_vs_MGBC_species" = "smClassify vs MGBC",
  "sm_vs_Refseq_species" = "smClassify vs RefSeq",
  "sm_vs_MGnify_species" = "smClassify vs MGnify",
  "MGBC_vs_Refseq_species" = "MGBC vs RefSeq",
  "MGBC_vs_MGnify_species" = "MGBC vs MGnify",
  "Refseq_vs_MGnify_species" = "RefSeq vs MGnify"
)

summary_pairwise_all <- bind_rows(
  conflict_df %>% count(level = "sm_vs_MGBC_genus", status = sm_vs_MGBC_genus),
  conflict_df %>% count(level = "sm_vs_Refseq_genus", status = sm_vs_Refseq_genus),
  conflict_df %>% count(level = "sm_vs_MGnify_genus", status = sm_vs_MGnify_genus),
  conflict_df %>% count(level = "MGBC_vs_Refseq_genus", status = MGBC_vs_Refseq_genus),
  conflict_df %>% count(level = "MGBC_vs_MGnify_genus", status = MGBC_vs_MGnify_genus),
  conflict_df %>% count(level = "Refseq_vs_MGnify_genus", status = Refseq_vs_MGnify_genus),
  conflict_df %>% count(level = "sm_vs_MGBC_species", status = sm_vs_MGBC_species),
  conflict_df %>% count(level = "sm_vs_Refseq_species", status = sm_vs_Refseq_species),
  conflict_df %>% count(level = "sm_vs_MGnify_species", status = sm_vs_MGnify_species),
  conflict_df %>% count(level = "MGBC_vs_Refseq_species", status = MGBC_vs_Refseq_species),
  conflict_df %>% count(level = "MGBC_vs_MGnify_species", status = MGBC_vs_MGnify_species),
  conflict_df %>% count(level = "Refseq_vs_MGnify_species", status = Refseq_vs_MGnify_species)
) %>%
  group_by(level) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

write_csv_table <- function(x, file) {
  data.table::fwrite(as.data.table(x), file = file, bom = TRUE)
}

write_csv_table(
  summary_pairwise_all,
  file.path(plot_dir, "SupplementaryFigure5d_pairwise_taxonomic_disagreement_source_data.csv")
)

# ============================================================
# 5. Supplementary Figure 5d
# ============================================================

plot_pair_conflict <- summary_pairwise_all %>%
  filter(status == "conflict") %>%
  mutate(
    tax_level = ifelse(grepl("_genus$", level), "Genus", "Species"),
    level = factor(level, levels = c(pair_order_genus, pair_order_species)),
    pair_label = factor(
      pair_label_map[as.character(level)],
      levels = c(
        "smClassify vs MGBC", "smClassify vs RefSeq", "smClassify vs MGnify",
        "MGBC vs RefSeq", "MGBC vs MGnify", "RefSeq vs MGnify"
      )
    )
  )

p_sfig5d <- ggplot(plot_pair_conflict, aes(x = pair_label, y = fraction)) +
  geom_col(width = 0.8, fill = "#BDB76B", color = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(fraction, accuracy = 0.1)), vjust = -0.4, size = 3) +
  facet_wrap(~tax_level, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  labs(
    x = NULL,
    y = "Conflict fraction",
    title = "Pairwise taxonomic disagreements"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(20, 20, 20, 20)
  )

save_plot(p_sfig5d, "SupplementaryFigure5d_pairwise_taxonomic_disagreement", width = 5.5, height = 4.5)

# ============================================================
# 6. Supplementary Figure 5e
# ============================================================

get_top_pairs_conflict <- function(df, col_a, col_b, status_col, pair_name, top_n = 5) {
  df %>%
    filter(.data[[status_col]] == "conflict") %>%
    mutate(species_pair = map2_chr(.data[[col_a]], .data[[col_b]], make_pair)) %>%
    filter(!is.na(species_pair)) %>%
    count(species_pair, sort = TRUE) %>%
    mutate(
      fraction = n / sum(n),
      comparison = pair_name
    ) %>%
    slice_head(n = top_n)
}

top_pairs_conflict <- bind_rows(
  get_top_pairs_conflict(conflict_df, "smClassify_species", "MGBC_species",
                         "sm_vs_MGBC_species", "smClassify vs MGBC"),
  get_top_pairs_conflict(conflict_df, "smClassify_species", "Refseq_species",
                         "sm_vs_Refseq_species", "smClassify vs RefSeq"),
  get_top_pairs_conflict(conflict_df, "smClassify_species", "MGnify_species",
                         "sm_vs_MGnify_species", "smClassify vs MGnify"),
  get_top_pairs_conflict(conflict_df, "MGBC_species", "Refseq_species",
                         "MGBC_vs_Refseq_species", "MGBC vs RefSeq"),
  get_top_pairs_conflict(conflict_df, "MGBC_species", "MGnify_species",
                         "MGBC_vs_MGnify_species", "MGBC vs MGnify"),
  get_top_pairs_conflict(conflict_df, "Refseq_species", "MGnify_species",
                         "Refseq_vs_MGnify_species", "RefSeq vs MGnify")
)

write_csv_table(
  top_pairs_conflict,
  file.path(plot_dir, "SupplementaryFigure5e_top5_conflict_species_pairs_source_data.csv")
)

plot_top_pairs <- top_pairs_conflict %>%
  mutate(
    species_pair = fct_reorder(species_pair, n),
    comparison = factor(
      comparison,
      levels = c(
        "smClassify vs MGBC", "smClassify vs RefSeq", "smClassify vs MGnify",
        "MGBC vs RefSeq", "MGBC vs MGnify", "RefSeq vs MGnify"
      )
    )
  )

p_sfig5e <- ggplot(plot_top_pairs, aes(x = species_pair, y = fraction)) +
  geom_col(fill = "#BDB76B", color = "black", linewidth = 0.2) +
  facet_wrap(~comparison, scales = "free_y") +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Fraction within conflicts",
    title = "Top conflicting species pairs"
  ) +
  theme_pub()

save_plot(p_sfig5e, "SupplementaryFigure5e_top5_conflict_species_pairs", width = 16, height = 4)

message("Supplementary Figure 5d-e outputs saved to: ", plot_dir)
