# ============================================================
# Figure 2 smClassify analysis

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
})

# ============================================================
base_dir   <- "Figure2"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

sce_path  <- file.path(input_dir, "Mouse_microbial_raw_seurat.rds")
attr_csv  <- file.path(input_dir, "Mouse_gene_annotation_rRNA_flags.csv")
cw_csv    <- file.path(input_dir, "Mouse_taxonomy_crosswalk_MGnify_RefSeq.csv")

work_dir  <- output_dir

qc_report_dir             <- file.path(output_dir, "QC_report")
raw_umi_raw_gene_dir      <- file.path(output_dir, "raw_umi_raw_gene")
taxonomic_composition_dir <- file.path(output_dir, "taxonomic_composition")

annotated_sce_rds  <- file.path(output_dir, "sce_bac_annotated.rds")
benchmark_calls_rds <- file.path(output_dir, "barcode_species_calls_benchmark.rds")

# QC_report outputs.
out_calls_csv      <- file.path(qc_report_dir, "barcode_species_calls.csv")
prefilter_qc_csv   <- file.path(qc_report_dir, "bacterial_prefilter_qc.csv")
plot_calls_csv     <- file.path(qc_report_dir, "species_calls_with_sample_metadata.csv")
qc_report_pdf      <- file.path(qc_report_dir, "QC_report_nature_style.pdf")
summary_dir        <- file.path(qc_report_dir, "smClassify_assignment_summary")

# raw_umi_raw_gene outputs.
raw_qc_csv         <- file.path(raw_umi_raw_gene_dir, "raw_cell_qc_metrics.csv")
raw_qc_pdf         <- file.path(raw_umi_raw_gene_dir, "raw_umi_raw_gene.pdf")

# taxonomic_composition outputs.
species_comp_csv   <- file.path(taxonomic_composition_dir, "Top20_species_composition_by_sample.csv")
family_comp_csv    <- file.path(taxonomic_composition_dir, "Family_composition_Top5.csv")
genus_comp_csv     <- file.path(taxonomic_composition_dir, "Genus_composition_Top22.csv")
tax_species_csv    <- file.path(taxonomic_composition_dir, "Species_composition_Top15.csv")
tax_comp_pdf       <- file.path(taxonomic_composition_dir, "Family_Genus_Species_TopN.pdf")
summary_composition_dir <- file.path(taxonomic_composition_dir, "smClassify_taxonomic_composition")

# Small wrappers for consistent table and multi-page PDF export.
write_csv_table <- function(x, file) {
  data.table::fwrite(as.data.table(x), file = file, bom = TRUE)
  message("[output] Table saved: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

save_pdf_report <- function(file, plots, width, height) {
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  invisible(lapply(plots, print))
  message("[output] Figure saved: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(raw_umi_raw_gene_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(taxonomic_composition_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(summary_composition_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# ============================================================

.blankish <- function(x){
  # Identify missing, empty, or whitespace-only strings.
  x <- as.character(x); x[is.na(x)] <- ""; !nzchar(trimws(x))
}
.norm_id  <- function(x){
  # Standardize feature identifiers by converting to upper case and replacing underscores with dashes.
  toupper(gsub("_","-", as.character(x), fixed = TRUE))
}
.get_mgyg <- function(x){
  # Extract MGYG genome identifiers from feature names.
  x2 <- toupper(as.character(x))
  m  <- regexpr("MGYG[_-]?[0-9]+", x2, perl = TRUE)
  ans <- ifelse(m > 0, regmatches(x2, m), NA_character_)
  gsub("MGYG[_-]([0-9]+)", "MGYG\\1", ans, perl = TRUE)  # Standardize to MGYG<digits>.
}
.norm_mgyg <- function(x){
  # Standardize MGYG identifiers by removing hyphen or underscore separators.
  toupper(gsub("^\\s*(MGYG)[_-]?(\\d+)\\s*$", "\\1\\2", as.character(x), perl = TRUE))
}

# Regular expressions for low-information species labels.
LOW_PAT_GENERAL <- "^\\s*([A-Za-z0-9_\\-]+)\\s+sp\\.?\\s*\\d*\\s*$"
GENUS_SP_UPPER  <- "^\\s*([A-Z0-9-]+)\\s+sp\\.?\\s*\\d*\\s*$"

# ============================================================
build_rrna_set <- function(attr_df) {
  if ("is_rRNA" %in% names(attr_df) || "is_tmRNA" %in% names(attr_df)) {
    isR <- if ("is_rRNA" %in% names(attr_df)) (attr_df$is_rRNA %in% TRUE) else FALSE
    isT <- if ("is_tmRNA" %in% names(attr_df)) (attr_df$is_tmRNA %in% TRUE) else FALSE
    ids <- unique(na.omit(as.character(attr_df$locus_tag[isR | isT])))
    return(.norm_id(ids))
  }
  # Annotation-based detection when explicit rRNA/tmRNA flags are unavailable.
  prod <- as.character(attr_df$product)
  mask <- (!is.na(prod)) &
    (grepl("(?i)(ribosomal\\s*RNA|\\b[0-9]+S\\s*rRNA\\b|RNASE\\s*P\\s*RNA)", prod, perl=TRUE) |
       grepl("(?i)\\btmRNA\\b|transfer-messenger\\s*RNA", prod, perl=TRUE)) &
    (!grepl("(?i)\\btRNA\\b|riboswitch", prod, perl=TRUE))
  .norm_id(unique(na.omit(as.character(attr_df$locus_tag[mask]))))
}


finalize_species_mgnify_logic <- function(cw) {
  stopifnot(is.data.frame(cw))
  cw <- as.data.table(cw)
  
  # Harmonize column names across crosswalk file versions.
  ln <- tolower(names(cw))
  ln <- sub("^refseq[ ._]?species$", "refseq.species", ln)
  names(cw) <- ln
  if (!"mgnify_species" %in% names(cw)) cw[, mgnify_species := NA_character_]
  if (!"refseq.species" %in% names(cw)) cw[, `refseq.species` := NA_character_]
  
  mg <- cw[["mgnify_species"]]
  rf <- cw[["refseq.species"]]
  mg_blank <- .blankish(mg)
  rf_blank <- .blankish(rf)
  
  # Canonical non-specific species-name patterns, for example "Bacteroides sp." and "Bacteroides sp1234".
  GENUS_SP_DOT <- "^\\s*([A-Za-z0-9_-]+)\\s+sp\\.$"       # Matches "Genus sp." labels.
  GENUS_SP_NUM <- "^\\s*([A-Za-z0-9_-]+)\\s+sp[0-9]+$"    # Matches "Genus sp1234" labels.
  mg_is_spdot <- !mg_blank & grepl(GENUS_SP_DOT, mg, perl = TRUE)
  mg_is_spnum <- !mg_blank & grepl(GENUS_SP_NUM, mg, perl = TRUE)
  rf_low <- !rf_blank & (grepl(GENUS_SP_DOT, rf, perl = TRUE) | grepl(GENUS_SP_NUM, rf, perl = TRUE))
  
  get_genus <- function(x, pat) sub(pat, "\\1", x, perl = TRUE)
  genus_allcaps <- function(genus) grepl("^[A-Z0-9-]+$", genus)
  
  # Default rule: retain the MGnify species label when present.
  chosen <- mg
  reason <- ifelse(mg_blank, "MGnify_blank", "MGnify_used")
  
  # Rule 1: use a specific RefSeq label when the MGnify label is missing.
  use_ref_blank <- mg_blank & !rf_blank & !rf_low
  chosen[use_ref_blank] <- rf[use_ref_blank]
  reason[use_ref_blank] <- "Use_RefSeq_due_to_MGnify_blank"
  
  # Rule 2: resolve "Genus sp." labels using RefSeq when possible; otherwise assign "Uncultured Genus".
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
  
  # Rule 3: resolve "Genus sp1234" labels using RefSeq when possible.
  # All-uppercase genus tokens are rewritten as "Uncultured Genus".
  if (any(mg_is_spnum, na.rm = TRUE)) {
    use_ref <- mg_is_spnum & !rf_blank & !rf_low
    chosen[use_ref] <- rf[use_ref]
    reason[use_ref] <- "Use_RefSeq_due_to_MGnify_spnum"
    
    still_sp <- mg_is_spnum & (rf_blank | rf_low | .blankish(rf))
    if (any(still_sp, na.rm = TRUE)) {
      genus_tok <- get_genus(mg[still_sp], GENUS_SP_NUM)
      is_allcaps <- genus_allcaps(genus_tok)
      chosen[still_sp] <- ifelse(is_allcaps, paste0("Uncultured ", genus_tok), mg[still_sp])
      reason[still_sp] <- ifelse(is_allcaps,
                                 "MGnify_spnum_allcaps_rewritten_to_Uncultured_Genus",
                                 "Keep_MGnify_spnum")
    }
  }
  
  # Final fallback for unresolved labels.
  chosen[.blankish(chosen)] <- "Unknown sp."
  
  cw[, Species := chosen]
  cw[, Species_reason := reason]
  cw
}


bacterial_prefilter <- function(
    sce, attr_df,
    assay = "RNA", slot = "counts",
    min_nonrRNA_count    = 100,
    min_nonrRNA_features = 30,
    max_rrna_frac        = 0.95,
    verbose = TRUE
){
  t0 <- proc.time()[3]
  stopifnot(inherits(sce, "Seurat"))
  if (verbose) message("[prefilter] Reading assay=", assay, " slot=", slot)
  M <- Seurat::GetAssayData(sce, assay = assay, slot = slot)
  if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  
  if (verbose) {
    message(sprintf("[prefilter] Raw matrix: %d genes x %d cells; density=%.4f",
                    nrow(M), ncol(M),
                    Matrix::nnzero(M) / as.double(nrow(M)) / as.double(ncol(M))))
  }
  
  rrna_ids <- build_rrna_set(attr_df)
  if (verbose) message(sprintf("[prefilter] rRNA/tmRNA set size: %d", length(rrna_ids)))
  
  rrna_mask <- rownames(M) %in% rrna_ids
  M_rrna <- if (any(rrna_mask)) M[rrna_mask,, drop = FALSE] else M[0,, drop = FALSE]
  M_non  <- if (any(rrna_mask)) M[!rrna_mask,, drop = FALSE] else M
  
  if (verbose) {
    message(sprintf("[prefilter] Non-rRNA matrix: %d x %d; nnz=%d",
                    nrow(M_non), ncol(M_non), Matrix::nnzero(M_non)))
  }
  
  total_umi      <- Matrix::colSums(M)
  nonr_umi       <- Matrix::colSums(M_non)
  nonr_features  <- Matrix::colSums(sign(M_non))
  rrna_fraction  <- if (nrow(M_rrna) > 0) Matrix::colSums(M_rrna) / pmax(1, total_umi) else rep(0, ncol(M))
  
  qc <- data.frame(
    total_umi     = as.numeric(total_umi),
    nonr_umi      = as.numeric(nonr_umi),
    nonr_features = as.numeric(nonr_features),
    rRNA_fraction = as.numeric(rrna_fraction),
    row.names = colnames(sce)
  )
  sce <- Seurat::AddMetaData(sce, metadata = qc)
  
  retain_cell <- (nonr_umi >= min_nonrRNA_count) &
    (nonr_features >= min_nonrRNA_features) &
    (rrna_fraction <= max_rrna_frac)
  
  sce$BacPrefilterRetained <- retain_cell
  sce_filt <- subset(sce, cells = colnames(sce)[retain_cell])
  
  kept <- sum(retain_cell); allc <- length(retain_cell)
  if (verbose) {
    message(sprintf("[prefilter] Retained: %d / %d cells (%.1f%%)", kept, allc, 100*kept/allc))
    message(sprintf("[prefilter] Elapsed %.2f s", proc.time()[3] - t0))
  }
  list(sce = sce_filt, qc = qc)
}

assign_species_per_barcode <- function(
    sce, attr_df, crosswalk,
    assay = "RNA", slot = "counts",
    min_nonrRNA_count = 100,
    min_support_counts = 50,
    min_frac = 0.6,
    min_margin = 0.2,
    write_tax_to_meta = TRUE,
    verbose = TRUE
){
  t0 <- proc.time()[3]
  stopifnot(inherits(sce, "Seurat"))
  if (verbose) message("[assign] Reading assay=", assay, " slot=", slot)
  
  M <- Seurat::GetAssayData(sce, assay = assay, slot = slot)
  if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  rownames(M) <- .norm_id(rownames(M))
  
  if (verbose) {
    message(sprintf("[assign] Raw matrix: %d x %d; nnz=%d",
                    nrow(M), ncol(M), Matrix::nnzero(M)))
  }
  
  total_umi <- Matrix::colSums(M)
  
  # Remove rRNA and tmRNA features before species assignment.
  if (verbose) message("[assign] Splitting rRNA vs. non-rRNA")
  rrna_ids <- build_rrna_set(attr_df)
  rrna_mask <- rownames(M) %in% rrna_ids
  M_rrna <- if (any(rrna_mask)) M[rrna_mask,,drop=FALSE] else M[0,,drop=FALSE]
  M_non  <- if (any(rrna_mask)) M[!rrna_mask,,drop=FALSE] else M
  if (nrow(M_non) == 0) stop("No usable non-rRNA features.")
  
  rrna_fraction <- if (nrow(M_rrna) > 0) Matrix::colSums(M_rrna) / pmax(1, total_umi) else rep(0, ncol(M))
  nonr_total <- Matrix::colSums(M_non)
  
  if (verbose) {
    message(sprintf("[assign] Non-rRNA matrix: %d x %d; nnz=%d",
                    nrow(M_non), ncol(M_non), Matrix::nnzero(M_non)))
  }
  
  # Map non-rRNA features to MGYG genome identifiers.
  if (verbose) message("[assign] Parsing MGYG ids from feature names")
  gene_mgyg <- .norm_mgyg(.get_mgyg(rownames(M_non)))
  
  # Prepare species and higher-rank taxonomy lookups.
  if (verbose) message("[assign] Normalizing crosswalk; mapping MGYG -> species")
  cw_low <- as.data.frame(crosswalk, check.names = FALSE)
  names(cw_low) <- tolower(names(cw_low))
  if (!"species" %in% names(cw_low)) {
    if ("mgnify_species" %in% names(cw_low)) cw_low$species <- cw_low$mgnify_species
    else if ("refseq.species" %in% names(cw_low)) cw_low$species <- cw_low$`refseq.species`
    else stop("Crosswalk lacks species column.")
  }
  if (!"genome" %in% names(cw_low)) stop("Crosswalk must include 'Genome'/'genome'.")
  cw_low$genome <- .norm_mgyg(cw_low$genome)
  
  mgyg2sp <- setNames(as.character(cw_low$species), cw_low$genome)
  species_for_gene <- unname(mgyg2sp[gene_mgyg])
  
  keep <- !is.na(species_for_gene) & species_for_gene != ""
  M_map <- M_non[keep,,drop=FALSE]
  sp_for_gene <- species_for_gene[keep]
  mgyg_non    <- gene_mgyg[keep]
  if (nrow(M_map) == 0) stop("No non-rRNA features can be mapped to species.")
  
  # Build lookup tables for higher taxonomic ranks.
  rank_keys <- c("Domain","Phylum","Class","Order","Family","Genus")
  rk_low <- tolower(rank_keys)
  need_cols <- c("genome", rk_low)
  if (!all(need_cols %in% names(cw_low))) {
    stop("Crosswalk must include: Genome + ", paste(rank_keys, collapse=", "))
  }
  mgyg2rank <- setNames(vector("list", length(rank_keys)), rank_keys)
  for (i in seq_along(rank_keys)) {
    rk <- rk_low[i]
    mgyg2rank[[i]] <- setNames(as.character(cw_low[[rk]]), cw_low$genome)
  }
  names(mgyg2rank) <- rank_keys
  
  # Aggregate feature-level counts into species-by-cell matrices using sparse matrix multiplication.
  if (verbose) message("[assign] Computing species-by-cell counts and support gene counts")
  t_mul <- proc.time()[3]
  sp_levels <- sort(unique(sp_for_gene))
  sp_index  <- match(sp_for_gene, sp_levels)
  P <- sparseMatrix(i = sp_index, j = seq_len(nrow(M_map)), x = 1,
                    dims = c(length(sp_levels), nrow(M_map)))
  S <- P %*% M_map           # species-by-cell count matrix
  nzM <- sign(M_map)
  S_gene_ct <- P %*% nzM     # species-by-cell supporting gene count matrix
  if (verbose) {
    message(sprintf("[assign] Sparse ops done in %.2f s | %d x %d; nnz=%d",
                    proc.time()[3]-t_mul, nrow(S), ncol(S), Matrix::nnzero(S)))
  }
  
  # Identify the top two species per cell using non-rRNA count fractions.
  if (verbose) message("[assign] Ranking Top1 / Top2 per cell")
  ST <- as(S, "dgTMatrix")
  ii <- ST@i + 1L
  jj <- ST@j + 1L
  vx <- ST@x
  frac <- vx / nonr_total[jj]
  dt <- data.table(j = jj, i = ii, frac = frac, cnt = vx)
  setkey(dt, j)
  setorder(dt, j, -frac)
  dt[, rk := seq_len(.N), by = j]
  top2 <- dt[rk <= 2L][, rk := NULL]
  
  top1 <- top2[, .SD[1], by = j]
  top2b <- top2[top2[, .I[.N >= 2L], by = j]$V1]
  top2b <- top2[j %in% top2b$j][, .SD[2], by = j]
  
  species_vec <- sp_levels
  species_top_vec    <- rep(NA_character_, ncol(M))
  species_second_vec <- rep(NA_character_, ncol(M))
  top_frac_vec       <- numeric(ncol(M))
  second_frac_vec    <- numeric(ncol(M))
  support_cnt_top    <- numeric(ncol(M))
  support_cnt_second <- numeric(ncol(M))
  
  species_top_vec[top1$j] <- species_vec[top1$i]
  top_frac_vec[top1$j]    <- top1$frac
  support_cnt_top[top1$j] <- top1$cnt
  if (nrow(top2b)) {
    species_second_vec[top2b$j] <- species_vec[top2b$i]
    second_frac_vec[top2b$j]    <- top2b$frac
    support_cnt_second[top2b$j] <- top2b$cnt
  }
  
  # Retrieve supporting gene counts for the top two species assignments.
  if (verbose) message("[assign] Getting Top1/Top2 supporting gene counts")
  SG <- as(S_gene_ct, "dgTMatrix")
  dtg <- data.table(j = SG@j+1L, i = SG@i+1L, g = as.integer(SG@x))
  setkey(dtg, j, i)
  top1_idx <- data.table(j = top1$j, i = top1$i)
  top1_idx[, g := dtg[top1_idx, on = .(j, i)]$g]
  top2_idx <- data.table(j = top2b$j, i = top2b$i)
  top2_idx[, g := dtg[top2_idx, on = .(j, i)]$g]
  support_genes_top    <- integer(ncol(M)); support_genes_top[top1_idx$j]    <- ifelse(is.na(top1_idx$g), 0L, top1_idx$g)
  support_genes_second <- integer(ncol(M)); support_genes_second[top2_idx$j] <- ifelse(is.na(top2_idx$g), 0L, top2_idx$g)
  
  # Assign higher taxonomic ranks by majority vote among genes supporting the top species.
  if (verbose) message("[assign] Majority-vote taxonomy (Top1)")
  nzM2 <- nzM; nzM2@x[] <- 1
  NZ <- as(nzM2, "dgTMatrix")
  dt_nz <- data.table(j = NZ@j+1L, gene_row = NZ@i+1L)
  dt_nz[, sp := sp_index[gene_row]]
  top1_for_join <- data.table(j = top1$j, sp = top1$i)
  setkey(dt_nz, j, sp); setkey(top1_for_join, j, sp)
  dt_top1_genes <- dt_nz[top1_for_join, nomatch = 0L]
  mgyg_vec <- mgyg_non
  rank_vals <- vector("list", length(rank_keys)); names(rank_vals) <- rank_keys
  for (rk in rank_keys) {
    m <- mgyg2rank[[rk]][mgyg_vec]
    rank_vals[[rk]] <- m
  }
  for (rk in rank_keys) {
    data.table::set(dt_top1_genes, j = rk, value = rank_vals[[rk]][dt_top1_genes$gene_row])
  }
  Mode1 <- function(v) {
    v <- v[!is.na(v) & nzchar(v)]
    if (!length(v)) return(NA_character_)
    tt <- table(v); names(tt)[which.max(tt)]
  }
  Tax_Domain <- Tax_Phylum <- Tax_Class <- Tax_Order <- Tax_Family <- Tax_Genus <- rep(NA_character_, ncol(M))
  if (nrow(dt_top1_genes)) {
    Tax_tbl <- dt_top1_genes[, .(
      Tax_Domain = Mode1(Domain),
      Tax_Phylum = Mode1(Phylum),
      Tax_Class  = Mode1(Class),
      Tax_Order  = Mode1(Order),
      Tax_Family = Mode1(Family),
      Tax_Genus  = Mode1(Genus)
    ), by = j]
    Tax_Domain[Tax_tbl$j] <- Tax_tbl$Tax_Domain
    Tax_Phylum[Tax_tbl$j] <- Tax_tbl$Tax_Phylum
    Tax_Class[Tax_tbl$j]  <- Tax_tbl$Tax_Class
    Tax_Order[Tax_tbl$j]  <- Tax_tbl$Tax_Order
    Tax_Family[Tax_tbl$j] <- Tax_tbl$Tax_Family
    Tax_Genus[Tax_tbl$j]  <- Tax_tbl$Tax_Genus
  }
  
  # Apply assignment-confidence and putative-doublet criteria.
  margin12 <- top_frac_vec - second_frac_vec
  high_confidence <- (!is.na(species_top_vec)) &
    (nonr_total >= min_nonrRNA_count) &
    (support_cnt_top >= min_support_counts) &
    (top_frac_vec >= min_frac) &
    (margin12 >= min_margin)
  confidence <- ifelse(high_confidence, "High", "Low")
  doublet_candidate <- (!is.na(species_second_vec)) &
    (nonr_total >= min_nonrRNA_count) &
    (support_cnt_top >= min_support_counts) &
    (top_frac_vec >= min_frac) &
    (margin12 < (2 * min_margin))
  
  # Build the per-barcode species-call table.
  calls <- data.table::data.table(
    barcode            = colnames(M),
    species            = species_top_vec,
    SpeciesSecond      = species_second_vec,
    nonr_total         = as.numeric(nonr_total),
    total_umi          = as.numeric(total_umi),
    rRNA_fraction      = as.numeric(rrna_fraction),
    top_frac           = top_frac_vec,
    SecondFrac         = second_frac_vec,
    margin12           = margin12,
    support_counts     = as.numeric(support_cnt_top),
    SecondCounts       = as.numeric(support_cnt_second),
    SupportGenes       = as.integer(support_genes_top),
    SecondGenes        = as.integer(support_genes_second),
    confidence         = confidence,
    SpeciesDoubletCandidate = doublet_candidate,
    Tax_Domain = Tax_Domain,
    Tax_Phylum = Tax_Phylum,
    Tax_Class  = Tax_Class,
    Tax_Order  = Tax_Order,
    Tax_Family = Tax_Family,
    Tax_Genus  = Tax_Genus
  )
  
  if (isTRUE(write_tax_to_meta)) {
    add_meta <- as.data.frame(calls[, .(species, SpeciesSecond, rRNA_fraction,
                                        Tax_Domain,Tax_Phylum,Tax_Class,Tax_Order,Tax_Family,Tax_Genus)],
                              stringsAsFactors = FALSE)
    rownames(add_meta) <- calls$barcode
    sce <<- Seurat::AddMetaData(sce, metadata = add_meta)
  }
  
  if (verbose) message(sprintf("[assign] Done in %.2f s", proc.time()[3] - t0))
  list(species_counts = S, calls = calls)
}

# ============================================================
# Main analysis workflow

t_all <- proc.time()[3]


# Load input Seurat object and gene annotation table.

sce <- readRDS(sce_path)
attr_df <- fread(attr_csv, check.names = FALSE)

# Step 1: prefilter cells using non-rRNA evidence.
pf <- bacterial_prefilter(
  sce, attr_df,
  min_nonrRNA_count = 100,
  min_nonrRNA_features = 30,
  max_rrna_frac = 0.95,
  verbose = TRUE
)
sce_bac <- pf$sce

# Step 2: read taxonomy crosswalk and reconcile species labels.
cw_raw  <- fread(cw_csv, check.names = FALSE)
if (!any(tolower(names(cw_raw)) == "genome")) stop("crosswalk must include 'Genome'")
names(cw_raw)[tolower(names(cw_raw)) == "genome"] <- "Genome"
cw <- finalize_species_mgnify_logic(cw_raw)

# Report mapping diagnostics for feature-to-taxonomy matching.
M_tmp <- Seurat::GetAssayData(sce_bac, assay = "RNA", slot = "counts")
if (!inherits(M_tmp, "dgCMatrix")) M_tmp <- as(M_tmp, "dgCMatrix")
cw_keys_norm <- .norm_mgyg(.get_mgyg(as.character(cw$genome)))
mgyg2sp <- setNames(
  as.character(cw$Species)[!is.na(cw_keys_norm)],
  cw_keys_norm[!is.na(cw_keys_norm)]
)
mg_from_rows <- .norm_mgyg(.get_mgyg(rownames(M_tmp)))

rrna_ids_now <- build_rrna_set(attr_df)

# Step 3a: assign species to individual cell barcodes for benchmarking.
res_benchmark <- assign_species_per_barcode(
  sce_bac, attr_df, cw,
  min_nonrRNA_count = 0,
  min_support_counts = 0,
  min_frac = 0.6,
  min_margin = 0.2,
  write_tax_to_meta = TRUE,
  verbose = TRUE
)
res_benchmark$calls

# Step 3b: assign species to individual cell barcodes for downstream clustering.
res_cluster <- assign_species_per_barcode(
  sce_bac, attr_df, cw,
  min_nonrRNA_count = 100,
  min_support_counts = 50,
  min_frac = 0.6,
  min_margin = 0.2,
  write_tax_to_meta = TRUE,
  verbose = TRUE
)
res_cluster$calls



# ============================================================
# Output block 1: QC_report

write_csv_table(res_benchmark$calls, out_calls_csv)
write_csv_table(pf$qc %>% tibble::rownames_to_column("barcode"), prefilter_qc_csv)

cols_to_write <- c(
  "species", "SpeciesSecond", "confidence", "SpeciesDoubletCandidate",
  "rRNA_fraction", "nonr_total", "total_umi",
  "top_frac", "SecondFrac", "margin12",
  "support_counts", "SecondCounts", "SupportGenes", "SecondGenes",
  "Tax_Domain", "Tax_Phylum", "Tax_Class", "Tax_Order", "Tax_Family", "Tax_Genus"
)

md_benchmark <- as.data.frame(res_benchmark$calls[, cols_to_write, with = FALSE])
rownames(md_benchmark) <- res_benchmark$calls$barcode

sce_bac_benchmark_annot <- Seurat::AddMetaData(sce_bac, metadata = md_benchmark)
Seurat::Idents(sce_bac_benchmark_annot) <- "species"

md_cluster <- as.data.frame(res_cluster$calls[, cols_to_write, with = FALSE])
rownames(md_cluster) <- res_cluster$calls$barcode

sce_bac_annot <- Seurat::AddMetaData(sce_bac, metadata = md_cluster)
Seurat::Idents(sce_bac_annot) <- "species"

saveRDS(sce_bac_benchmark_annot, file = benchmark_calls_rds)
saveRDS(sce_bac_annot, file = annotated_sce_rds)

# ============================================================
# Quality-control data tables
M_raw <- Seurat::GetAssayData(sce, assay = "RNA", slot = "counts")
if (!inherits(M_raw, "dgCMatrix")) M_raw <- as(M_raw, "dgCMatrix")

raw_qc <- data.frame(
  barcode     = colnames(M_raw),
  total_umi   = Matrix::colSums(M_raw),
  gene_count  = Matrix::colSums(sign(M_raw)),
  orig.ident  = if ("orig.ident" %in% colnames(sce@meta.data)) sce$orig.ident else "(batch)",
  loc         = if ("loc" %in% colnames(sce@meta.data)) sce$loc else "(all)"
)

meta_df <- sce_bac@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(
    orig.ident = if (!"orig.ident" %in% names(.)) "(batch)" else orig.ident,
    loc        = if (!"loc"        %in% names(.)) "(all)"   else loc
  )

res_calls <- res_benchmark$calls %>%
  left_join(meta_df[, c("barcode","orig.ident","loc","rRNA_fraction")], by = "barcode") %>%
  mutate(
    confidence = factor(confidence, levels = c("Low", "High")),
    SpeciesDoubletCandidate = factor(SpeciesDoubletCandidate, levels = c(FALSE, TRUE),
                                     labels = c("No", "Yes"))
  )

# ============================================================
# Quality-control figures

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

# Color palettes for phenotype and taxonomy visualizations.
cell_type_cols <- c("#1F77B4",  "#2CA02C","#FF7F0E","#6A5ACD", "#8C564B", "#E377C2", "#7F7F7F", 
                    "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32",
                    "#8B0000",  "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B", 
                    "#20B2AA", "#FF1493", "#FF4500", "#32CD32", "#3E8E41", "#20B2AA")
okabe_ito <- c("#000000","#E69F00","#56B4E9","#009E73",
               "#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

table(raw_qc$orig.ident)
sample_order <- c(
  "Cecum-WT-1",
  "Cecum-WT-2",
  "Cecum-DB-1",
  "Cecum-DB-2",
  "Colon-WT-1",
  "Colon-WT-2",
  "Colon-DB-1",
  "Colon-DB-2",
  "Rectum-WT-1",
  "Rectum-WT-2",
  "Rectum-DB-1",
  "Rectum-DB-2"
)

position_order <- c("Cecum", "Colon", "Rectum")

raw_qc$orig.ident <- factor(raw_qc$orig.ident, levels = sample_order)
raw_qc$loc <- factor(raw_qc$loc, levels = position_order)
res_calls$orig.ident <- factor(res_calls$orig.ident, levels = sample_order)
res_calls$loc <- factor(res_calls$loc, levels = position_order)
pal_pheno <- c(
  WT = "#4DBBD5",
  DB = "#E64B35"
)
raw_qc$pheno <- ifelse(grepl("WT", raw_qc$orig.ident), "WT", "DB")
# Figure: raw total UMI counts per cell.
p_raw_umi <- ggplot(raw_qc, aes(x = orig.ident, y = total_umi, fill = pheno)) +
  geom_violin(colour = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  scale_y_log10() +
  scale_fill_manual(values = pal_pheno) +
  labs(title = "Rawdata: UMI counts per cell", y = "Total UMI (log10)", x = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_raw_umi

# Figure: raw detected genes per cell.
p_raw_gene <- ggplot(raw_qc, aes(x = orig.ident, y = gene_count, fill = pheno)) +
  geom_violin(colour = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  scale_y_log10() +
  scale_fill_manual(values = pal_pheno) +
  labs(title = "Rawdata: Gene counts per cell", y = "Detected genes (log10)", x = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_raw_gene
# Figure: species-assignment confidence by sample.
p_confidence <- ggplot(res_calls, aes(x = orig.ident, fill = confidence)) +
  geom_bar(position = "fill", colour = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = okabe_ito[c(9,3)]) +
  labs(title = "Species assignment confidence by sample", y = "Fraction", x = NULL, fill = "Confidence") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_confidence

# Figure: distribution of the Top1-minus-Top2 species assignment margin.
# The dashed line indicates the min_margin threshold used for assignment confidence.
p_margin <- ggplot(res_calls, aes(x = margin12)) +
  geom_histogram(bins = 60, fill = "grey75", colour = "black") +
  geom_vline(xintercept = 0.2, linetype = "dashed") +
  labs(title = "Confidence margin (Top1 - Top2)", x = "Margin12", y = "Cells") +
  theme_nature()
p_margin 
# Figure: putative doublet fraction based on the species-assignment margin rule.
p_doublet <- ggplot(res_calls, aes(x = orig.ident, fill = SpeciesDoubletCandidate)) +
  geom_bar(position = "fill", colour = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = okabe_ito[c(2,7)]) +
  labs(title = "Fraction of cells with ambiguous species", 
       y = "Fraction", x = NULL, fill = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_doublet
# Figure: top-20 species composition by intestinal region and sample.
species_tbl <- table(res_calls$species)
species_tbl <- species_tbl[!is.na(names(species_tbl))]
top_k <- min(20, length(species_tbl))
top_species <- if (top_k > 0) names(sort(species_tbl, decreasing = TRUE))[1:top_k] else character(0)

comp_species <- res_calls %>%
  mutate(
    species_plot = ifelse(is.na(species), "Unknown", species),
    species_plot = ifelse(species_plot %in% top_species, species_plot, "Others"),
    orig.ident = factor(orig.ident, levels = sample_order),
    loc = factor(loc, levels = position_order)
  ) %>%
  group_by(loc, orig.ident, species_plot) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  group_by(loc, orig.ident) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  arrange(loc, orig.ident, desc(frac))

p_comp <- ggplot(comp_species, aes(x = orig.ident, y = frac, fill = species_plot)) +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~loc, scales = "free_x") +
  labs(title = "Top20 species composition", y = "Composition", x = NULL, fill = "Species") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


write_csv_table(res_calls, plot_calls_csv)
write_csv_table(comp_species, species_comp_csv)

save_pdf_report(
  file   = qc_report_pdf,
  plots  = list(p_confidence, p_margin, p_doublet, p_comp),
  width  = 4,
  height = 3
)

# ============================================================
# Output block 2: raw_umi_raw_gene

write_csv_table(raw_qc, raw_qc_csv)

save_pdf_report(
  file   = raw_qc_pdf,
  plots  = list(p_raw_umi, p_raw_gene),
  width  = 7.2,
  height = 3
)

# ============================================================
# Output block 3: taxonomic_composition
#
# Build and save family-, genus-, and species-level composition summaries.
# Non-top taxa are grouped as Others.
# ============================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(purrr)
})

# Null-coalescing helper used when optional plotting parameters are not provided.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Default color palette used when no external palette has been defined.
if (!exists("cell_type_cols", inherits = FALSE)) {
  cell_type_cols <- c(
    "#1F77B4","#2CA02C","#FF7F0E","#6A5ACD","#8C564B","#E377C2","#7F7F7F",
    "#BCBD22","#17BECF","#F5A0A1","#C2B5D8","#FFB6C1","#3CB371","#9ACD32",
    "#8B0000","#FFD700","#DC143C","#228B22","#FF6347","#483D8B","#BDB76B",
    "#20B2AA","#FF1493","#FF4500","#32CD32","#3E8E41","#20B2AA"
  )
}

# Extract sample and intestinal-region metadata from the annotated Seurat object when available.
.meta_src <- if (exists("sce_bac_annot")) sce_bac_annot else sce
stopifnot(!is.null(.meta_src))

meta_df2 <- .meta_src@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(
    sample   = if ("orig.ident" %in% names(.)) as.character(orig.ident) else "(batch)",
    position = if ("loc"        %in% names(.)) as.character(loc)        else "(all)"
  ) %>%
  dplyr::select(barcode, sample, position)

# Merge species calls with sample metadata and harmonize taxonomy-rank labels.
calls2 <- res_benchmark$calls %>%
  left_join(meta_df2, by = "barcode") %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    position = factor(position, levels = position_order),
    Species = ifelse(is.na(species) | species == "", "Unknown", species),
    Genus   = ifelse(is.na(Tax_Genus)  | Tax_Genus  == "", "Unknown", Tax_Genus),
    Family  = ifelse(is.na(Tax_Family) | Tax_Family == "", "Unknown", Tax_Family)
  )

# Compute Top-N composition for a selected taxonomic rank; all remaining taxa are grouped as Others.
.build_comp_rank <- function(df_calls, rank_col = c("Family", "Genus", "Species"), topN = 20) {
  rank_col <- match.arg(rank_col)
  
  df_rank <- df_calls %>%
    dplyr::mutate(
      rank_value = as.character(.data[[rank_col]]),
      rank_value = ifelse(is.na(rank_value) | rank_value == "", "Unknown", rank_value),
      sample = factor(sample, levels = sample_order),
      position = factor(position, levels = position_order)
    )
  
  top_taxa <- df_rank %>%
    dplyr::count(rank_value, name = "n_tot") %>%
    dplyr::arrange(dplyr::desc(n_tot)) %>%
    dplyr::slice_head(n = topN) %>%
    dplyr::pull(rank_value)
  
  comp <- df_rank %>%
    dplyr::mutate(
      taxon = ifelse(rank_value %in% top_taxa, rank_value, "Others")
    ) %>%
    dplyr::count(position, sample, taxon, name = "n") %>%
    dplyr::group_by(position, sample) %>%
    dplyr::mutate(Percent = round(100 * n / sum(n), 2)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(position, sample, dplyr::desc(Percent))
  
  comp
}

# Generate composition tables for family, genus, and species levels.


make_named_palette <- function(df, base_cols) {
  taxa <- levels(df$taxon)
  if (length(base_cols) < length(taxa)) {
    # Extend the palette automatically if the number of taxa exceeds the base palette length.
    extra <- scales::hue_pal()(length(taxa) - length(base_cols))
    cols <- c(base_cols, extra)
  } else {
    cols <- base_cols
  }
  names(cols) <- taxa
  cols
}

plot_rank_nice <- function(df, rank_name, base_cols = cell_type_cols,
                           sample_order = NULL, position_order = NULL,
                           show_pct_labels = FALSE) {
  df <- df %>%
    mutate(
      position = factor(position, levels = position_order %||% unique(position)),
      sample   = factor(sample,   levels = sample_order   %||% unique(sample)),
      taxon    = fct_relevel(taxon, "Others", after = Inf)
    )
  
  pal <- make_named_palette(df, base_cols)
  
  p <- ggplot(df, aes(sample, Percent, fill = taxon)) +
    geom_col(width = 0.85) +
    facet_wrap(~position, scales = "free_x") +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(
      title = paste0(rank_name, " composition (Top", topN_map[[rank_name]], " + Others)"),
      x = "Sample", y = "Percentage", fill = "taxon"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank()
    )
  
  if (show_pct_labels) {
    p <- p + geom_text(
      aes(label = ifelse(Percent >= 5, paste0(Percent, "%"), "")),
      position = position_stack(vjust = 0.5), size = 3, color = "white"
    )
  }
  p
}

comp_list <- list(
  Family  = .build_comp_rank(calls2, "Family",  topN = 5),
  Genus   = .build_comp_rank(calls2, "Genus",   topN = 10),
  Species = .build_comp_rank(calls2, "Species", topN = 15)
)
# Top-N values used in plot titles.
topN_map <- list(Family = 5, Genus = 10, Species = 15)


# Generate family-, genus-, and species-level composition plots.
p_family  <- plot_rank_nice(comp_list$Family,  "Family",
                            sample_order = sample_order, position_order = position_order)
p_genus   <- plot_rank_nice(comp_list$Genus,   "Genus",
                            sample_order = sample_order, position_order = position_order)
p_species <- plot_rank_nice(comp_list$Species, "Species",
                            sample_order = sample_order, position_order = position_order)
p_family
p_genus
p_species

write_csv_table(comp_list$Family,  family_comp_csv)
write_csv_table(comp_list$Genus,   genus_comp_csv)
write_csv_table(comp_list$Species, tax_species_csv)

save_pdf_report(
  file   = tax_comp_pdf,
  plots  = list(p_family, p_genus, p_species),
  width  = 6,
  height = 6
)
