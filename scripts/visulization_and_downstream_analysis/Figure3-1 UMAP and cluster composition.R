# ============================================================

# This script performs downstream visualization for the
# mouse microbiome single-cell dataset, after species
# assignment and basic QC have been completed.
#
# Major steps:
#   A) Remove rRNA/tmRNA features based on genome annotations.
#   B) Align a GTF-wide annotation table to gene IDs.
#   C) Aggregate counts from locus-level features to gene symbols.
#   D) Rebuild a Seurat object with aggregated symbols and
#      carry over cell-level metadata.
#   E) Perform light QC (feature filter, per-sample downsampling,
#      basic UMI thresholding).
#   F) SCTransform, dimensional reduction, clustering, UMAP.
#   G) Add derived metadata (group, position, phenotype/region).
#   H) Generate a set of UMAP plots (by species, group, sample,
#      phenotype, region, cluster) and a bar plot of cluster
#      composition by group.
#
# Notes:
#   - All parameter values, thresholds, and file paths are kept
#     identical to your original script.
#   - Only comments / documentation have been expanded for clarity.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(future)
  library(tibble)
  library(paletteer)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)             # for grid.text, gpar
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(Seurat)
  library(scales)           # hue_pal, percent_format, etc.
})

# ----------------------------------------------------------------
# Global options
# ----------------------------------------------------------------
set.seed(1234)   # Fix random seed for reproducibility

work_dir <- "D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009"
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

# A 21-color palette: ggthemes::Classic_20 + black
color21 <- c(paletteer_d("ggthemes::Classic_20", n = 20), "black")

# ============================================================
# Helper functions
# ============================================================

# .blankish:
#   Returns TRUE for NA / empty / whitespace-only entries.
.blankish <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  !nzchar(trimws(x))
}

# .norm_id:
#   Standardize IDs by:
#     - Converting to uppercase
#     - Replacing "_" with "-" (for robust matching)
.norm_id  <- function(x) {
  toupper(gsub("_","-", as.character(x), fixed = TRUE))
}

# build_rrna_set:
#   Identify rRNA / tmRNA features using:
#     1) Explicit boolean columns is_rRNA / is_tmRNA if present.
#     2) Otherwise, regular expression matching on 'product'.
#   Returns a vector of normalized locus_tag IDs for rRNA / tmRNA.
build_rrna_set <- function(attr_df) {
  stopifnot(is.data.frame(attr_df))
  
  # Preferred path: explicit flags
  if ("is_rRNA" %in% names(attr_df) || "is_tmRNA" %in% names(attr_df)) {
    isR <- if ("is_rRNA" %in% names(attr_df)) attr_df$is_rRNA %in% TRUE else FALSE
    isT <- if ("is_tmRNA" %in% names(attr_df)) attr_df$is_tmRNA %in% TRUE else FALSE
    return(.norm_id(unique(na.omit(attr_df$locus_tag[isR | isT]))))
  }
  
  # Fallback: product-based inference
  prod <- as.character(attr_df$product)
  mask <- (!is.na(prod)) &
    (
      grepl("(?i)(ribosomal\\s*RNA|[0-9]+S\\s*rRNA|RNASE\\s*P\\s*RNA)", prod, perl = TRUE) |
      grepl("(?i)tmRNA|transfer-messenger", prod, perl = TRUE)
    ) &
    (!grepl("(?i)tRNA|riboswitch", prod, perl = TRUE))
  
  .norm_id(unique(na.omit(attr_df$locus_tag[mask])))
}

# align_gtf_to_genes:
#   Align a wide GTF annotation table (gtf_wide) to rownames of the
#   expression matrix (genes).
#   - Tries multiple key columns in order: locus_tag, ID, gene_id, Name.
#   - Fills a canonical 'gene_id' column based on normalized matching.
#   - Returns a unique table with columns: gene_id, symbol, ko, product
#     (limited to those present in the input table).
align_gtf_to_genes <- function(gtf_wide, genes) {
  setDT(gtf_wide)
  genes_norm <- .norm_id(genes)
  key_order <- c("locus_tag", "ID", "gene_id", "Name")
  
  gtf_wide[, gene_id := as.character(gene_id)]
  
  for (k in key_order) {
    to_fill <- .blankish(gtf_wide$gene_id)
    cand    <- .norm_id(gtf_wide[[k]])
    hit     <- to_fill & !.blankish(cand) & (cand %in% genes_norm)
    gtf_wide[hit, gene_id := cand[hit]]
  }
  
  keep_cols <- intersect(c("gene_id","symbol","ko","product"), names(gtf_wide))
  unique(gtf_wide[! .blankish(gene_id), ..keep_cols])
}

# fast_rowsum_by_group:
#   Fast sparse row aggregation: group features (rows) according to
#   a grouping vector gvec and sum the counts per group.
#   - M: dgCMatrix (genes x cells)
#   - gvec: character vector of group labels (length == nrow(M))
#   - min_total: drop original rows with total < min_total before grouping
#   - Returns a dgCMatrix with one row per unique group label.
fast_rowsum_by_group <- function(M, gvec, min_total = 1L, verbose = TRUE) {
  if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  stopifnot(length(gvec) == nrow(M))
  
  rs <- Matrix::rowSums(M)
  keep <- !is.na(gvec) & rs >= min_total
  M2 <- M[keep,,drop = FALSE]
  grp <- gvec[keep]
  
  f <- factor(grp, levels = unique(grp))
  G <- Matrix::sparseMatrix(
    i    = as.integer(f),
    j    = seq_along(f),
    x    = 1L,
    dims = c(nlevels(f), length(f))
  )
  
  out <- G %*% M2
  rownames(out) <- levels(f)
  
  if (verbose) {
    message(sprintf("Grouped sum done: %d groups × %d cells; nnz=%d",
                    nrow(out), ncol(out), Matrix::nnzero(out)))
  }
  out
}

# ============================================================
# Main workflow
# ============================================================

# ------------------------------------------------------------
# Step A: Remove rRNA/tmRNA genes from the Seurat object
# ------------------------------------------------------------
seu <- readRDS("sce_bac_annotated.rds")
attr_df <- fread("genome_feature_attributes_M1.csv", check.names = FALSE)

rrna_ids <- build_rrna_set(attr_df)
seu <- seu[!(rownames(seu) %in% rrna_ids), ]

# ------------------------------------------------------------
# Step B: Align GTF-wide gene annotation to expression matrix
# ------------------------------------------------------------
gtf_wide <- fread("gtf_gene_annotation.csv", check.names = FALSE)
genes_nr <- rownames(GetAssayData(seu, "RNA", "counts"))
sym_map  <- align_gtf_to_genes(gtf_wide, genes_nr)

# Drop entries with "hypothetical protein" (low-information annotation)
if ("product" %in% names(sym_map)) {
  sym_map <- sym_map[!(tolower(product) == "hypothetical protein"), ]
}

# ------------------------------------------------------------
# Step C: Aggregate counts by gene symbol
# ------------------------------------------------------------
# For non-MGYG symbols, strip trailing "-number" or "_number" to merge
# variants from the same locus/gene.
sym_map <- sym_map %>%
  mutate(symbol = ifelse(
    grepl("^MGYG", symbol),
    symbol,
    gsub("[-_][0-9]+$", "", symbol)
  ))

# Build mapping from original gene_id -> gene symbol, ordered to match genes_nr
gvec <- setNames(sym_map$symbol, sym_map$gene_id)[.norm_id(genes_nr)]

# Aggregate counts across genes that share the same symbol
M_symbol <- fast_rowsum_by_group(
  GetAssayData(seu, "RNA", "counts"),
  gvec,
  min_total = 5
)

# ------------------------------------------------------------
# Step D: Rebuild Seurat object and carry over metadata
# ------------------------------------------------------------
# We keep the original cell metadata (columns of meta.data)
# and attach it to the new aggregated count matrix.
seu <- CreateSeuratObject(
  counts    = M_symbol,
  meta.data = seu@meta.data[colnames(M_symbol), ]
)
Idents(seu) <- seu$species

# ------------------------------------------------------------
# Step E: Light QC filters
# ------------------------------------------------------------

# (1) Keep features (symbols) expressed in >= 5 cells
keep_features <- names(which(Matrix::rowSums(GetAssayData(seu, "RNA", "counts") > 0) >= 5))
seu <- subset(seu, features = keep_features)

# (2) Downsample to at most 7000 cells per sample (orig.ident)
cells_keep <- seu@meta.data %>%
  rownames_to_column("cell") %>%
  group_by(orig.ident) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(nrow(.x), 7000))) %>%
  pull(cell)
seu <- subset(seu, cells = cells_keep)

# (3) Basic UMI count filter: retain cells with 30 < nCount_RNA < 1000
seu <- subset(seu, subset = nCount_RNA > 30 & nCount_RNA < 1000)

# ------------------------------------------------------------
# Step F: SCTransform, PCA, neighbors, clustering, UMAP
# ------------------------------------------------------------
plan(sequential)
options(future.globals.maxSize = 32 * 1024^3)

seu <- SCTransform(
  seu,
  vst.flavor          = "v2",
  vars.to.regress     = "nCount_RNA",
  variable.features.n = 5000,
  conserve.memory     = TRUE,
  method = if (requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else "poisson"
)

seu <- RunPCA(seu, npcs = 50, assay = "SCT")
seu <- FindNeighbors(seu, dims = 1:6)
seu <- FindClusters(seu, resolution = 0.1)
seu <- RunUMAP(seu, dims = 1:6)

DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size  = .1,
  raster   = FALSE,
  cols     = color21
)

saveRDS(seu, "stepF_seu_SCT_clustered.rds")

# ============================================================
# Step 1. Add higher-level metadata (group / position / phenotype)
# ============================================================

meta <- seu@meta.data

# group:
#   Remove trailing "-digit" from orig.ident, so that
#   samples like "Cecum-WT-1" / "Cecum-WT-2" share the same group.
group_vec <- sub("-[0-9]$", "", meta$orig.ident)

# position:
#   Derive anatomical region (Cecum / Colon / Rectum) from orig.ident prefix.
position_vec <- case_when(
  grepl("^Cecum",  meta$orig.ident) ~ "Cecum",
  grepl("^Colon",  meta$orig.ident) ~ "Colon",
  grepl("^Rectum", meta$orig.ident) ~ "Rectum",
  TRUE ~ NA_character_
)

# phenotype:
#   Label any sample containing "WT" as WT; otherwise DB.
phenotype_vec <- ifelse(grepl("WT", meta$orig.ident), "WT", "DB")

# Assign new metadata columns
seu$group     <- group_vec
seu$position  <- position_vec
seu$phenotype <- phenotype_vec
seu$region    <- seu$position   # 'region' is a copy of 'position' for convenience

# Set factor levels explicitly for consistent plotting order
seu$phenotype <- factor(seu$phenotype, levels = c("WT", "DB"))
seu$position  <- factor(seu$position,  levels = c("Cecum", "Colon", "Rectum"))
seu$region    <- factor(seu$region,    levels = c("Cecum", "Colon", "Rectum"))
seu$group     <- factor(seu$group,     levels = unique(group_vec))

# Quick summary of new annotations
message("Summary after grouping:")
print(table(Phenotype = seu$phenotype, useNA = "ifany"))
print(table(Position  = seu$position,  useNA = "ifany"))
print(head(seu@meta.data[, c("orig.ident","group","phenotype","position","region")], 8))

# ============================================================
# Step 2. UMAP plotting utilities
# ============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# expand_palette:
#   Given a base color vector and desired length n,
#   expand or truncate to get a palette of length n.
#   Uses scales::hue_pal() if base_cols is empty or too short.
expand_palette <- function(base_cols, n) {
  if (is.null(base_cols) || length(base_cols) == 0) return(scales::hue_pal()(n))
  if (length(base_cols) >= n) return(unname(base_cols[seq_len(n)]))
  c(unname(base_cols), scales::hue_pal()(n - length(base_cols)))
}

# make_named_palette:
#   Build a named palette for given factor levels.
make_named_palette <- function(levels_vec, base_cols) {
  cols <- expand_palette(base_cols, length(levels_vec))
  names(cols) <- levels_vec
  cols
}

# .ensure_factor_levels:
#   Ensure that the metadata column 'group.by' is a factor.
.ensure_factor_levels <- function(obj, group.by) {
  stopifnot(inherits(obj, "Seurat"))
  if (!group.by %in% colnames(obj@meta.data)) {
    stop("`group.by` not found in metadata: ", group.by)
  }
  v <- obj[[group.by, drop = TRUE]]
  if (!is.factor(v)) v <- factor(as.character(v))
  obj@meta.data[[group.by]] <- v
  obj
}

# plot_umap:
#   General-purpose UMAP plotting wrapper with:
#   - flexible grouping (group.by)
#   - custom palettes
#   - consistent minimal theme suitable for figures
#   - automatic PNG/PDF export if filename is provided.
plot_umap <- function(
    obj, group.by, title = NULL,
    palette = NULL,
    filename = NULL,
    out_dir = ".",
    pt.size = 0.1, label = FALSE, raster = FALSE,
    w = 12, h = 8, dpi = 600
) {
  obj <- .ensure_factor_levels(obj, group.by)
  lv  <- levels(obj[[group.by, drop = TRUE]])
  
  if (!is.null(names(palette)) && all(lv %in% names(palette))) {
    pal_named <- palette[lv]
  } else {
    pal_named <- make_named_palette(lv, base_cols = palette %||% character())
  }
  
  p <- DimPlot(
    obj, reduction = "umap", group.by = group.by,
    pt.size = pt.size, raster = raster, cols = pal_named, label = label
  ) +
    theme_minimal(base_size = 30) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      legend.title = element_text(size = 35, face = "bold"),
      legend.text  = element_text(size = 33),
      plot.title   = element_text(size = 34, face = "bold", hjust = 0.5)
    ) +
    labs(title = title %||% paste("UMAP by", group.by))
  
  if (!is.null(filename)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    png_path <- file.path(out_dir, filename)
    pdf_path <- sub("\\.png$", ".pdf", png_path)
    
    ggsave(png_path, p, width = w, height = h, dpi = dpi, units = "in", device = "png")
    ggsave(pdf_path, p, width = w, height = h, units = "in", device = "pdf")
  }
  p
}

# plot_topK_species_umap:
#   Collapse rare species into "Others" and plot UMAP by
#   top-K most abundant species across samples.
plot_topK_species_umap <- function(
    seu,
    species_col = "species",
    sample_col  = "orig.ident",
    top_k = 20,
    others_label = "Others",
    palette_topk = c(paletteer::paletteer_d("ggthemes::Classic_20", n = 20), "black"),
    out_dir = ".", filename = "UMAP_Top20_Species.png",
    pt.size = 0.1, raster = FALSE, w = 12, h = 5, dpi = 600
) {
  stopifnot(species_col %in% colnames(seu@meta.data))
  stopifnot(sample_col  %in% colnames(seu@meta.data))
  
  top_names <- seu@meta.data %>%
    as_tibble() %>%
    count(.data[[species_col]], .data[[sample_col]], name = "n") %>%
    group_by(.data[[species_col]]) %>%
    summarise(total_count = sum(n), .groups = "drop") %>%
    arrange(desc(total_count)) %>%
    slice_head(n = top_k) %>%
    pull(!!sym(species_col))
  
  seu@meta.data <- seu@meta.data %>%
    mutate(topK_species = ifelse(.data[[species_col]] %in% top_names,
                                 .data[[species_col]], others_label))
  
  plot_umap(
    obj      = seu,
    group.by = "topK_species",
    title    = paste0("UMAP by Top", top_k, " Species"),
    palette  = palette_topk,
    filename = filename,
    out_dir  = out_dir,
    pt.size  = pt.size,
    raster   = raster,
    w        = w,
    h        = h,
    dpi      = dpi
  )
}

# ============================================================
# Step 3. UMAP plotting for multiple groupings
# ============================================================

# Pre-defined palettes used for different facets
nature_palette  <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
                     "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
nature_palette2 <- c("#4c78a8","#f58518","#54a24b","#e45756","#72b7b2",
                     "#b279a2","#ff9da6","#9d755d","#bab0ab","#edc948")

nature_palette3 <- c("#edc948","#4c78a8","#f58518","#54a24b","#e45756",
                     "#72b7b2","#b279a2","#ff9da6")

pal_pheno  <- c(WT = "#4DBBD5", T2DM = "#E64B35")
pal_region <- c(Cecum = "#1f77b4", Colon = "#2ca02c", Rectum = "#ff7f0e")

out_dir <- "fig_umap"

# ------------------------------------------------------------
# UMAP by Top20 species (collapse others)
# ------------------------------------------------------------
plot_topK_species_umap(
  seu,
  species_col = "species",
  sample_col  = "orig.ident",
  top_k       = 20,
  others_label = "Others",
  out_dir     = out_dir,
  filename    = "UMAP_Top20_Species.png",
  pt.size     = 0.1,
  raster      = FALSE,
  w           = 24,
  h           = 10,
  dpi         = 600
)

# Explicit Top20 species factor in metadata
sp_counts <- sort(table(seu$species), decreasing = TRUE)
n_top     <- min(20, length(sp_counts))
top_species <- names(sp_counts)[seq_len(n_top)]

seu$top20_species <- ifelse(
  seu$species %in% top_species,
  as.character(seu$species),
  "Others"
)
seu$top20_species <- factor(
  seu$top20_species,
  levels = c(top_species, "Others")
)
n_other <- sum(!(seu$species %in% top_species))
levels(seu$top20_species)[levels(seu$top20_species) == "Others"]
# Quick summary
table(seu$top20_species)

# ------------------------------------------------------------
# UMAP by Group
# ------------------------------------------------------------
plot_umap(
  seu,
  "group",
  "UMAP by Group",
  palette  = nature_palette3,
  filename = "UMAP_group.png",
  out_dir  = out_dir,
  w        = 12,
  h        = 9,
  dpi      = 600
)

# ------------------------------------------------------------
# UMAP by Sample (orig.ident) – palette 1
# ------------------------------------------------------------
plot_umap(
  seu,
  "orig.ident",
  "UMAP by Sample",
  palette  = nature_palette,
  filename = "UMAP_Sample.png",
  out_dir  = out_dir,
  w        = 12,
  h        = 9,
  dpi      = 600
)

# ------------------------------------------------------------
# UMAP by Phenotype (WT vs T2DM)
# ------------------------------------------------------------
plot_umap(
  seu,
  "phenotype",
  "UMAP by Phenotype",
  palette  = pal_pheno,
  filename = "UMAP_Phenotype.png",
  out_dir  = out_dir,
  w        = 12,
  h        = 9,
  dpi      = 600
)

# ------------------------------------------------------------
# UMAP by Region (Cecum / Colon / Rectum)
# ------------------------------------------------------------
plot_umap(
  seu,
  "position",
  "UMAP by Region",
  palette  = pal_region,
  filename = "UMAP_Region.png",
  out_dir  = out_dir,
  w        = 12,
  h        = 9,
  dpi      = 600
)

# ------------------------------------------------------------
# UMAP by Cluster (seurat_clusters)
# ------------------------------------------------------------
plot_umap(
  seu,
  "seurat_clusters",
  "UMAP by Cluster",
  palette  = nature_palette,
  filename = "UMAP_Clusters.png",
  out_dir  = out_dir,
  w        = 11.5,
  h        = 9,
  dpi      = 600
)

# ------------------------------------------------------------
# UMAP by Sample (orig.ident) – palette 2 (alternative look)
# ------------------------------------------------------------
plot_umap(
  seu,
  "orig.ident",
  "UMAP by Sample",
  palette  = nature_palette2,
  filename = "UMAP_Sample.png",
  out_dir  = out_dir,
  w        = 20,
  h        = 8,
  dpi      = 600
)

# ============================================================
# Cluster composition barplot: proportion of clusters per group
# ============================================================

seu@active.ident <- as.factor(seu$group)

# Proportion table: Group × Cluster
prop_df <- as.data.frame(prop.table(table(Idents(seu), seu$seurat_clusters)))
colnames(prop_df) <- c("Group", "Cluster", "Freq")
prop_df$Group <- factor(prop_df$Group, levels = unique(prop_df$Group))

# Stacked bar: cluster composition by group
p_bar <- ggplot(prop_df, aes(x = Cluster, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = expand_palette(nature_palette3, length(levels(prop_df$Group)))) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(size = 11, face = "bold"),
    plot.title    = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(
    x     = "Cluster",
    y     = "Proportion",
    title = "Cluster Composition by Group"
  )

p_bar
ggsave("Cluster_Composition.png", p_bar, w = 4.5, h = 3, dpi = 600)
