############################################################
## Supplementary Figure 7d species support comparison
## Kraken2/Bracken metagenomics versus smClassify
############################################################

# ============================================================
# 0. Load packages
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

set.seed(1234)

# ============================================================
# 1. Input files and output directory
# ============================================================

base_dir   <- "Figure2"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")
plot_dir   <- file.path(output_dir, "Fig_S7d_species_support_Kraken2_smClassify")

kraken_report_dir <- "D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20260319/smClassify-metagenome/Metagenome_kraken2/reports"

sce_file <- file.path(
  input_dir,
  "Mouse_microbial_taxonomy_smClassify_Kraken2_MGnify.rds"
)

crosswalk_file <- file.path(
  input_dir,
  "Mouse_taxonomy_crosswalk_MGnify_RefSeq.csv"
)

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

kraken_report_files <- list.files(
  kraken_report_dir,
  pattern = "\\.report$",
  full.names = TRUE
)

stopifnot(length(kraken_report_files) > 0)
stopifnot(file.exists(sce_file))
stopifnot(file.exists(crosswalk_file))

# ============================================================
# 2. Utility functions
# ============================================================

.blankish <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  !nzchar(trimws(x))
}

finalize_species_mgnify_logic <- function(cw) {
  stopifnot(is.data.frame(cw))
  
  cw <- as.data.table(cw)
  
  ln <- tolower(names(cw))
  ln <- sub("^refseq[ ._]?species$", "refseq.species", ln)
  names(cw) <- ln
  
  if (!"mgnify_species" %in% names(cw)) {
    cw[, mgnify_species := NA_character_]
  }
  
  if (!"refseq.species" %in% names(cw)) {
    cw[, refseq.species := NA_character_]
  }
  
  mg <- cw[["mgnify_species"]]
  rf <- cw[["refseq.species"]]
  
  mg_blank <- .blankish(mg)
  rf_blank <- .blankish(rf)
  
  GENUS_SP_DOT <- "^\\s*([A-Za-z0-9_-]+)\\s+sp\\.$"
  GENUS_SP_NUM <- "^\\s*([A-Za-z0-9_-]+)\\s+sp[0-9]+$"
  
  mg_is_spdot <- !mg_blank & grepl(GENUS_SP_DOT, mg, perl = TRUE)
  mg_is_spnum <- !mg_blank & grepl(GENUS_SP_NUM, mg, perl = TRUE)
  
  rf_low <- !rf_blank & (
    grepl(GENUS_SP_DOT, rf, perl = TRUE) |
      grepl(GENUS_SP_NUM, rf, perl = TRUE)
  )
  
  get_genus <- function(x, pat) {
    sub(pat, "\\1", x, perl = TRUE)
  }
  
  genus_allcaps <- function(genus) {
    grepl("^[A-Z0-9-]+$", genus)
  }
  
  chosen <- mg
  reason <- ifelse(mg_blank, "MGnify_blank", "MGnify_used")
  
  use_ref_blank <- mg_blank & !rf_blank & !rf_low
  chosen[use_ref_blank] <- rf[use_ref_blank]
  reason[use_ref_blank] <- "Use_RefSeq_due_to_MGnify_blank"
  
  if (any(mg_is_spdot, na.rm = TRUE)) {
    use_ref <- mg_is_spdot & !rf_blank & !rf_low
    chosen[use_ref] <- rf[use_ref]
    reason[use_ref] <- "Use_RefSeq_due_to_MGnify_spdot"
    
    still_sp <- mg_is_spdot & (rf_blank | rf_low | .blankish(rf))
    
    if (any(still_sp, na.rm = TRUE)) {
      genus_tok <- get_genus(mg[still_sp], GENUS_SP_DOT)
      chosen[still_sp] <- paste0("Uncultured ", genus_tok)
      reason[still_sp] <- "MGnify_spdot_rewritten_to_Uncultured_Genus"
    }
  }
  
  if (any(mg_is_spnum, na.rm = TRUE)) {
    use_ref <- mg_is_spnum & !rf_blank & !rf_low
    chosen[use_ref] <- rf[use_ref]
    reason[use_ref] <- "Use_RefSeq_due_to_MGnify_spnum"
    
    still_sp <- mg_is_spnum & (rf_blank | rf_low | .blankish(rf))
    
    if (any(still_sp, na.rm = TRUE)) {
      genus_tok <- get_genus(mg[still_sp], GENUS_SP_NUM)
      is_allcaps <- genus_allcaps(genus_tok)
      
      chosen[still_sp] <- ifelse(
        is_allcaps,
        paste0("Uncultured ", genus_tok),
        mg[still_sp]
      )
      
      reason[still_sp] <- ifelse(
        is_allcaps,
        "MGnify_spnum_allcaps_rewritten_to_Uncultured_Genus",
        "Keep_MGnify_spnum"
      )
    }
  }
  
  chosen[.blankish(chosen)] <- "Unknown sp."
  
  cw[, Species := chosen]
  cw[, Species_reason := reason]
  
  cw
}

read_kraken_species_report <- function(file) {
  df <- fread(file, header = FALSE, fill = TRUE)
  
  colnames(df) <- c(
    "percent",
    "clade_reads",
    "taxon_reads",
    "rank",
    "taxid",
    "name"
  )
  
  sample_name <- basename(file) %>%
    str_remove("\\.report$")
  
  df[, name := str_trim(name)]
  df[, sample := sample_name]
  
  df <- df[rank == "S"]
  df <- df[!is.na(name) & name != ""]
  
  df[, .(
    sample = sample,
    species_raw = as.character(name),
    support_reads = as.numeric(taxon_reads),
    percent = as.numeric(percent)
  )]
}

get_counts_matrix <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) {
      GetAssayData(obj, assay = "RNA", slot = "counts")
    }
  )
}

# ============================================================
# 3. Read Kraken2 species reports
# ============================================================

kraken_species_reads_long <- rbindlist(
  lapply(kraken_report_files, read_kraken_species_report),
  fill = TRUE
)

fwrite(
  kraken_species_reads_long,
  file.path(plot_dir, "Fig_S7d_Kraken2_species_raw_long.csv")
)

# ============================================================
# 4. Harmonize Kraken2 species names
# ============================================================

cw_raw <- fread(crosswalk_file, check.names = FALSE)
cw_hm <- finalize_species_mgnify_logic(cw_raw)

cw_hm_df <- as.data.frame(cw_hm, check.names = FALSE)
names(cw_hm_df) <- tolower(names(cw_hm_df))

species_map <- cw_hm_df %>%
  as_tibble() %>%
  transmute(
    species_raw = as.character(mgnify_species),
    species_harmonized = as.character(species),
    species_reason = as.character(species_reason)
  ) %>%
  filter(!is.na(species_raw), species_raw != "") %>%
  mutate(
    species_harmonized = ifelse(
      is.na(species_harmonized) | species_harmonized == "",
      "Unknown sp.",
      species_harmonized
    ),
    species_reason = ifelse(is.na(species_reason), "", species_reason)
  ) %>%
  distinct(species_raw, .keep_all = TRUE)

fwrite(
  species_map,
  file.path(plot_dir, "Fig_S7d_species_harmonization_map.csv")
)

kraken_species_reads_hm_long <- kraken_species_reads_long %>%
  left_join(species_map, by = "species_raw") %>%
  mutate(
    species_harmonized = ifelse(
      is.na(species_harmonized) | species_harmonized == "",
      species_raw,
      species_harmonized
    )
  ) %>%
  group_by(sample, species_harmonized) %>%
  summarise(
    support_reads = sum(support_reads, na.rm = TRUE),
    percent = sum(percent, na.rm = TRUE),
    species_raw_merged = paste(unique(species_raw), collapse = " | "),
    .groups = "drop"
  )

fwrite(
  kraken_species_reads_hm_long,
  file.path(plot_dir, "Fig_S7d_Kraken2_species_harmonized_long.csv")
)

kraken_species_support <- kraken_species_reads_hm_long %>%
  group_by(species_harmonized) %>%
  summarise(
    kraken_support = mean(support_reads, na.rm = TRUE),
    kraken_median_support = median(support_reads, na.rm = TRUE),
    kraken_total_support = sum(support_reads, na.rm = TRUE),
    kraken_mean_percent = mean(percent, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(kraken_support))

fwrite(
  kraken_species_support,
  file.path(plot_dir, "Fig_S7d_Kraken2_species_support_summary.csv")
)

# ============================================================
# 5. Summarize smClassify species support
# ============================================================

obj <- readRDS(sce_file)

counts_mat <- get_counts_matrix(obj)
cell_total_umi <- Matrix::colSums(counts_mat)

cell_df <- data.table(
  species_harmonized = as.character(obj$species),
  total_umi = as.numeric(cell_total_umi)
)

cell_df <- cell_df[
  !is.na(species_harmonized) &
    species_harmonized != "" &
    !species_harmonized %in% c("Unknown", "Unknown sp.", "NA", "unassigned")
]

sm_species_support <- cell_df[, .(
  sm_support = mean(total_umi, na.rm = TRUE),
  sm_total_support = sum(total_umi, na.rm = TRUE),
  n_cells = .N
), by = species_harmonized][order(-sm_support)]

fwrite(
  sm_species_support,
  file.path(plot_dir, "Fig_S7d_smClassify_species_support_summary.csv")
)

# ============================================================
# 6. Merge Kraken2 and smClassify species support
# ============================================================

comp_df <- full_join(
  kraken_species_support %>%
    select(species_harmonized, kraken_support),
  sm_species_support %>%
    as_tibble() %>%
    select(species_harmonized, sm_support),
  by = "species_harmonized"
) %>%
  mutate(
    kraken_support = ifelse(is.na(kraken_support), 0, kraken_support),
    sm_support = ifelse(is.na(sm_support), 0, sm_support)
  )

fwrite(
  comp_df,
  file.path(plot_dir, "Fig_S7d_Kraken2_smClassify_species_support_summary.csv")
)

# ============================================================
# 7. Prepare source data
# ============================================================

top_n <- 20

top_species_union <- unique(c(
  comp_df %>%
    arrange(desc(sm_support)) %>%
    slice_head(n = top_n) %>%
    pull(species_harmonized),
  comp_df %>%
    arrange(desc(kraken_support)) %>%
    slice_head(n = top_n) %>%
    pull(species_harmonized)
))

species_order <- comp_df %>%
  filter(species_harmonized %in% top_species_union) %>%
  mutate(max_support = pmax(sm_support, kraken_support, na.rm = TRUE)) %>%
  arrange(max_support) %>%
  pull(species_harmonized)

plot_df <- comp_df %>%
  filter(species_harmonized %in% top_species_union) %>%
  pivot_longer(
    cols = c(sm_support, kraken_support),
    names_to = "Method",
    values_to = "Support"
  ) %>%
  mutate(
    Method = recode(
      Method,
      sm_support = "smClassify",
      kraken_support = "Metagenome Kraken2"
    ),
    Support_plot = ifelse(Support <= 0, NA, Support),
    species_harmonized = factor(species_harmonized, levels = species_order)
  )

fwrite(
  plot_df,
  file.path(plot_dir, "Fig_S7d_species_support_Kraken2_smClassify_SourceData.csv")
)

# ============================================================
# 8. Plot Supplementary Figure 7d
# ============================================================

p_s7d <- ggplot(
  plot_df,
  aes(x = Support_plot, y = species_harmonized, fill = Method)
) +
  geom_col(position = "dodge", width = 0.70, na.rm = TRUE) +
  scale_x_log10() +
  scale_fill_manual(values = c(
    "smClassify" = "#BCBD22",
    "Metagenome Kraken2" = "#F4A3A8"
  )) +
  labs(
    title = "Support for top species",
    x = "Support reads / UMI (log10)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = c(0.78, 0.18),
    legend.background = element_blank(),
    legend.key.size = unit(0.35, "cm")
  )

ggsave(
  file.path(plot_dir, "Fig_S7d_species_support_Kraken2_smClassify.pdf"),
  p_s7d,
  width = 4.2,
  height = 5.3
)

ggsave(
  file.path(plot_dir, "Fig_S7d_species_support_Kraken2_smClassify.png"),
  p_s7d,
  width = 4.2,
  height = 5.3,
  dpi = 300
)

print(p_s7d)

