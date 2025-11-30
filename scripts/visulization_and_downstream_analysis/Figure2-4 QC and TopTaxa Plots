# ============================================================
# smClassify-Mouse microbiome QC and Top-taxa composition figures
#
# This script generates:
#   - Figure 2F: Top-20 Family / Genus / Species composition tables
#   - Figure 2G: raw_umi_raw_gene.pdf (UMI / gene count distributions)
#
# Expected inputs (under project_dir/output):
#   - sce_bac_annotated.rds
#       Seurat object with species & taxonomy annotated.
#   - barcode_species_calls.csv
#       Per-barcode species calls (output of smClassify_mouse_species_assignment.R).
#
# Additional input (under project_dir/input):
#   - Raw_mouse_microbial_sce.rds
#       Original Seurat object used for raw UMI / gene QC.
#
# Outputs (written to project_dir/output):
#   - QC_report_nature_style.pdf    (margin, doublet, Top-20 species composition)
#   - raw_umi_raw_gene.pdf          (Figure 2G)
#   - comp_list (list of data.frames for Figure 2F: Family/Genus/Species)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# ============================================================
# Paths for GitHub version
# ============================================================
project_dir <- "path/to/project_root"  # <-- modify this for your environment

input_dir   <- file.path(project_dir, "input")
output_dir  <- file.path(project_dir, "output")

dir.create(input_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sce_raw_path   <- file.path(input_dir,  "Raw_mouse_microbial_sce.rds")
sce_annot_path <- file.path(output_dir, "sce_bac_annotated.rds")
calls_csv      <- file.path(output_dir, "barcode_species_calls.csv")
out_calls_csv  <- file.path(output_dir, "barcode_species_calls.csv")  # keep identical name

# ============================================================
# Load objects
# ============================================================
message("[plots] Loading Seurat objects and call table ...")
sce           <- readRDS(sce_raw_path)      # raw object (for raw UMI / genes, same as original)
sce_bac_annot <- readRDS(sce_annot_path)    # annotated object
res_calls_dt  <- fread(calls_csv, check.names = FALSE)

# For consistency with dplyr code below
res <- list(calls = res_calls_dt)

# ============================================================
# Common helpers / palettes (identical to original)
# ============================================================
`%||%` <- function(x, y) if (is.null(x)) y else x

theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      axis.line   = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.4, colour = "black"),
      axis.title  = element_text(size = base_size),
      axis.text   = element_text(size = base_size - 1, colour = "black"),
      plot.title  = element_text(face = "bold", hjust = 0, size = base_size + 1),
      legend.title= element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      plot.margin = margin(6,6,6,6)
    )
}

cell_type_cols <- c(
  "#1F77B4",  "#2CA02C","#FF7F0E","#6A5ACD", "#8C564B", "#E377C2", "#7F7F7F", 
  "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32",
  "#8B0000",  "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B", 
  "#20B2AA", "#FF1493", "#FF4500", "#32CD32", "#3E8E41", "#20B2AA"
)
okabe_ito <- c("#000000","#E69F00","#56B4E9","#009E73",
               "#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

# ============================================================
# -------- Figure 2F: Top-20 Family / Genus / Species --------
#   - Build composition tables per position × sample
#   - Parameters & thresholds identical to original:
#       * Family: Top 5
#       * Genus:  Top 10
#       * Species: Top 15
# ============================================================

if (!exists("cell_type_cols", inherits = FALSE)) {
  cell_type_cols <- c(
    "#1F77B4","#2CA02C","#FF7F0E","#6A5ACD","#8C564B","#E377C2","#7F7F7F",
    "#BCBD22","#17BECF","#F5A0A1","#C2B5D8","#FFB6C1","#3CB371","#9ACD32",
    "#8B0000","#FFD700","#DC143C","#228B22","#FF6347","#483D8B","#BDB76B",
    "#20B2AA","#FF1493","#FF4500","#32CD32","#3E8E41","#20B2AA"
  )
}

.meta_src <- sce_bac_annot  # in original: if (exists("sce_bac_annot")) ... else sce
stopifnot(!is.null(.meta_src))

meta_df2 <- .meta_src@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(
    sample   = if ("orig.ident" %in% names(.)) as.character(orig.ident) else "(batch)",
    position = if ("loc"        %in% names(.)) as.character(loc)        else "(all)"
  ) %>%
  dplyr::select(barcode, sample, position)

calls2 <- res$calls %>%
  left_join(meta_df2, by = "barcode") %>%
  mutate(
    Species = ifelse(is.na(species) | species == "", "Unknown", species),
    Genus   = ifelse(is.na(Tax_Genus)  | Tax_Genus  == "", "Unknown", Tax_Genus),
    Family  = ifelse(is.na(Tax_Family) | Tax_Family == "", "Unknown", Tax_Family)
  )

.build_comp_rank <- function(df_calls, rank_col = c("Family","Genus","Species"), topN = 20) {
  rank_col <- match.arg(rank_col)
  top_taxa <- df_calls %>%
    count(!!sym(rank_col), name = "n_tot") %>%
    arrange(desc(n_tot)) %>%
    slice_head(n = topN) %>%
    pull(!!sym(rank_col))
  
  comp <- df_calls %>%
    mutate(taxon = ifelse(.data[[rank_col]] %in% top_taxa, .data[[rank_col]], "Others")) %>%
    count(position, sample, taxon, name = "n") %>%
    group_by(position, sample) %>%
    mutate(Percent = round(100 * n / sum(n), 2)) %>%
    ungroup() %>%
    arrange(position, sample, desc(Percent))
  
  comp
}

comp_list <- list(
  Family  = .build_comp_rank(calls2, "Family",  topN = 5),
  Genus   = .build_comp_rank(calls2, "Genus",   topN = 10),
  Species = .build_comp_rank(calls2, "Species", topN = 15)
)

topN_map <- list(Family = 5, Genus = 10, Species = 15)

# 如果你有自己的绘图函数，可以在这里生成 Figure 2F：
# p_family  <- plot_rank_nice(comp_list$Family,  "Family")
# p_genus   <- plot_rank_nice(comp_list$Genus,   "Genus")
# p_species <- plot_rank_nice(comp_list$Species, "Species")
# pdf(file.path(output_dir, "Family_Genus_Species_Top20.pdf"), width = 6, height = 6)
# print(p_family)
# print(p_genus)
# print(p_species)
# dev.off()

message("[plots] Top-20 Family / Genus / Species composition tables (Figure 2F) built in `comp_list`.")

# ============================================================
# -------- Figure 2G: raw_umi_raw_gene & QC report -----------
#   - Raw UMI and gene counts (violin/box, log10)
#   - Assignment success / doublets / Top-20 species composition
#   - Outputs:
#       * QC_report_nature_style.pdf
#       * raw_umi_raw_gene.pdf
# ============================================================

M_raw <- Seurat::GetAssayData(sce, assay = "RNA", slot = "counts")
if (!inherits(M_raw, "dgCMatrix")) M_raw <- as(M_raw, "dgCMatrix")

raw_qc <- data.frame(
  barcode     = colnames(M_raw),
  total_umi   = Matrix::colSums(M_raw),
  gene_count  = Matrix::colSums(sign(M_raw)),
  orig.ident  = if ("orig.ident" %in% colnames(sce@meta.data)) sce$orig.ident else "(batch)",
  loc         = if ("loc" %in% colnames(sce@meta.data)) sce$loc else "(all)"
)

meta_df <- sce_bac_annot@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(
    orig.ident = if (!"orig.ident" %in% names(.)) "(batch)" else orig.ident,
    loc        = if (!"loc"        %in% names(.)) "(all)"   else loc
  )

res_calls <- res$calls %>%
  left_join(meta_df[, c("barcode","orig.ident","loc","rRNA_fraction")], by = "barcode") %>%
  mutate(
    pass = factor(pass, levels = c(FALSE, TRUE), labels = c("Fail","Pass")),
    SpeciesDoubletCandidate = factor(
      SpeciesDoubletCandidate,
      levels = c(FALSE, TRUE),
      labels = c("Singlet","Doublet")
    )
  )

# Plot 1: Raw total UMI per cell
p_raw_umi <- ggplot(raw_qc, aes(x = orig.ident, fill = orig.ident, y = total_umi)) +
  geom_violin(colour = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  scale_y_log10() +
  labs(title = "Rawdata: UMI counts per cell", y = "Total UMI (log10)", x = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = cell_type_cols)

# Plot 2: Raw detected genes per cell
p_raw_gene <- ggplot(raw_qc, aes(x = orig.ident, fill = orig.ident, y = gene_count)) +
  geom_violin(colour = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  scale_y_log10() +
  labs(title = "Rawdata: Gene counts per cell", y = "Detected genes (log10)", x = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = cell_type_cols)

# Plot 3: Species assignment success rate per sample
p_success <- ggplot(res_calls, aes(x = orig.ident, fill = pass)) +
  geom_bar(position = "fill", colour = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = okabe_ito[c(6,3)]) +
  labs(title = "Species assignment success rate", y = "Fraction", x = NULL, fill = "Call") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 4: Margin (Top1 - Top2) distribution
p_margin <- ggplot(res_calls, aes(x = margin12)) +
  geom_histogram(bins = 60, fill = "grey75", colour = "black") +
  geom_vline(xintercept = 0.2, linetype = "dashed") +
  labs(title = "Confidence margin (Top1 - Top2)", x = "Margin12", y = "Cells") +
  theme_nature()

# Plot 5: Putative doublet fraction
p_doublet <- ggplot(res_calls, aes(x = orig.ident, fill = SpeciesDoubletCandidate)) +
  geom_bar(position = "fill", colour = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = okabe_ito[c(2,7)]) +
  labs(title = "Putative doublet fraction", y = "Fraction", x = NULL, fill = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 6: Top-20 species composition per loc x sample
species_tbl <- table(res_calls$species)
species_tbl <- species_tbl[!is.na(names(species_tbl))]
top_k <- min(20, length(species_tbl))
top_species <- if (top_k > 0) names(sort(species_tbl, decreasing = TRUE))[1:top_k] else character(0)

comp_species <- res_calls %>%
  mutate(
    species_plot = ifelse(is.na(species), "Unknown", species),
    species_plot = ifelse(species_plot %in% top_species, species_plot, "Others")
  ) %>%
  group_by(loc, orig.ident, species_plot) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  group_by(loc, orig.ident) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

p_comp <- ggplot(comp_species, aes(x = orig.ident, y = frac, fill = species_plot)) +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~loc, scales = "free_x") +
  labs(title = "Top20 species composition", y = "Composition", x = NULL, fill = "Species") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# Export PDFs (Figure 2G + QC report)
# ============================================================
pdf(file.path(output_dir, "QC_report_nature_style.pdf"), width = 4, height = 3)
print(p_margin)
print(p_doublet)
print(p_comp)
dev.off()

pdf(file.path(output_dir, "raw_umi_raw_gene.pdf"), width = 7.2, height = 3)
print(p_raw_umi)
print(p_raw_gene)
dev.off()

# Keep same call-table output
fwrite(res$calls, out_calls_csv, bom = TRUE)
message("[plots] QC report saved: ", file.path(output_dir, "QC_report_nature_style.pdf"))
message("[plots] Species calls saved: ", out_calls_csv)
message("[plots] Figure 2F (comp_list) and Figure 2G (raw_umi_raw_gene.pdf) generation complete.")
