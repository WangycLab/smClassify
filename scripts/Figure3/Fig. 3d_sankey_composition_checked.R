library(data.table)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(KEGGREST)
library(stringr)
library(readr)
library(tidyr)

set.seed(1)

base_dir <- "Figure3"
input_dir <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

fig3b_dir <- file.path(output_dir, "Figure3b_KEGG_enrichment")
cache_dir <- file.path(fig3b_dir, "cache")
enrich_dir <- file.path(fig3b_dir, "cluster_enrichment")
combined_dir <- file.path(enrich_dir, "combined")

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

gtf_file <- file.path(input_dir, "Figure3_mouse_Mgnify_genomes.gtf")
top_marker_csv <- file.path(input_dir, "Total cells_markers_top50.csv")
hierarchy_csv <- file.path(input_dir, "Figure3_KEGG_pathway_hierarchy_clean.csv")

padj_cut <- 0.25
data_min_count <- 1
top_n_plot <- 20

.norm_id <- function(x) {
  toupper(gsub("_", "-", as.character(x), fixed = TRUE))
}

.blankish <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  !nzchar(trimws(x))
}

get_attr_gtf <- function(x, key) {
  p <- paste0("(^|[;\\t\\s])", key, '\\s+"?([^";]+)"?')
  m <- regexpr(p, x, perl = TRUE)
  ans <- rep(NA_character_, length(x))
  hit <- m > 0
  ans[hit] <- sub(p, "\\2", regmatches(x, m)[hit], perl = TRUE)
  ans
}

split_kos <- function(s) {
  if (is.na(s) || !nzchar(s)) return(character(0))
  toks <- unlist(strsplit(s, "[,;\\s|]+"))
  toks <- gsub("^(ko:|KEGG:)", "", toks, ignore.case = TRUE)
  toks <- toupper(toks)
  toks[grepl("^K\\d{5}$", toks)]
}

.safe_dir <- function(x) {
  gsub('[\\\\/:*?"<>|]', "_", as.character(x))
}

read_markers <- function(DT, alias_long) {

  DT <- data.table::as.data.table(DT)
  names(DT) <- tolower(names(DT))

  if ("gene" %in% names(DT)) {
    DT[, gene_input := gene]
  } else if ("gene_id" %in% names(DT)) {
    DT[, gene_input := gene_id]
  } else {
    DT[, gene_input := symbol]
  }

  DT[, gene_input := .norm_id(gene_input)]
  DT[, cluster := as.character(cluster)]

  mapped <- merge(
    DT[, .(gene_input, cluster)],
    alias_long,
    by.x = "gene_input",
    by.y = "alias",
    all.x = TRUE
  )

  mapped[is.na(gene_id), gene_id := gene_input]
  mapped[, .(gene_id, cluster)]
}

get_interest_universe <- function(sub_df, KO_by_gene, gene_universe) {

  gids <- unique(na.omit(sub_df$gene_id))

  interest <- unique(KO_by_gene$ko[KO_by_gene$gene_id %in% gids])
  universe <- unique(KO_by_gene$ko[KO_by_gene$gene_id %in% gene_universe])

  list(interest = interest, universe = universe)
}

run_enrich <- function(interest, universe, TERM2GENE, TERM2NAME,
                       padj_cut = 0.25, min_count = 1) {

  space <- unique(TERM2GENE$gene)
  interest <- intersect(interest, space)
  universe <- intersect(universe, space)

  ek <- clusterProfiler::enricher(
    gene = interest,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    universe = universe,
    pAdjustMethod = "BH",
    minGSSize = 1,
    maxGSSize = 5000,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )

  df <- as.data.frame(ek)

  dplyr::filter(df, p.adjust <= padj_cut, Count >= min_count)
}

plot_dot <- function(df, title, png_file, pdf_file) {

  df$Description[is.na(df$Description)] <- df$ID[is.na(df$Description)]

  p <- ggplot(df, aes(
    x = p.adjust,
    y = reorder(Description, -p.adjust)
  )) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_x_reverse() +
    labs(title = title, x = "p.adjust", y = NULL) +
    theme_bw(base_size = 12)

  ggsave(png_file, p, width = 8, height = 6, dpi = 180)
  ggsave(pdf_file, p, width = 8, height = 6)
}

DT <- fread(
  gtf_file,
  sep = "\t",
  header = FALSE,
  select = c(3, 9),
  col.names = c("feature", "attributes")
)

A <- DT$attributes

lt <- get_attr_gtf(A, "locus_tag")
nm <- get_attr_gtf(A, "Name")
gn <- get_attr_gtf(A, "gene")
prd <- get_attr_gtf(A, "product")
keg <- get_attr_gtf(A, "KEGG")

gene_id <- .norm_id(lt)

symbol <- nm
symbol[.blankish(symbol)] <- gn[.blankish(symbol)]
symbol[.blankish(symbol)] <- lt[.blankish(symbol)]
symbol[.blankish(symbol)] <- prd[.blankish(symbol)]

symbol <- .norm_id(symbol)

gtf_map <- data.table(
  gene_id = gene_id,
  symbol = symbol,
  KEGG = keg
)

gtf_map <- gtf_map[!is.na(gene_id)]

gtf_map <- gtf_map[, {
  ks <- split_kos(KEGG)
  if (length(ks) == 0) list(ko = NA_character_) else list(ko = ks)
}, by = .(gene_id, symbol)]

gtf_map <- unique(gtf_map)

alias_long <- unique(rbind(
  data.table(alias = gtf_map$gene_id, gene_id = gtf_map$gene_id),
  data.table(alias = gtf_map$symbol, gene_id = gtf_map$gene_id)
))

KO_T2G <- na.omit(as.data.table(gtf_map)[, .(term = ko, gene = gene_id)])
gene_universe <- unique(KO_T2G$gene)

kegg_cache_ok <- all(file.exists(file.path(cache_dir,
  c("Figure3b_PATH2KO.rds",
    "Figure3b_MOD2KO.rds",
    "Figure3b_PATH_NM.rds",
    "Figure3b_MOD_NM.rds"))))

if (kegg_cache_ok) {
  PATH2KO <- readRDS(file.path(cache_dir, "Figure3b_PATH2KO.rds"))
  MOD2KO  <- readRDS(file.path(cache_dir, "Figure3b_MOD2KO.rds"))
  PATH_NM <- readRDS(file.path(cache_dir, "Figure3b_PATH_NM.rds"))
  MOD_NM  <- readRDS(file.path(cache_dir, "Figure3b_MOD_NM.rds"))
} else {

  lk_path <- KEGGREST::keggLink("pathway", "ko")
  lk_mod  <- KEGGREST::keggLink("module", "ko")

  PATH2KO <- data.table(
    term = sub("^ko:", "", names(lk_path)),
    gene = sub("^ko:", "", lk_path)
  )

  MOD2KO <- data.table(
    term = sub("^ko:", "", names(lk_mod)),
    gene = sub("^ko:", "", lk_mod)
  )

  PATH_NM <- data.frame(
    term = names(KEGGREST::keggList("pathway")),
    name = as.character(KEGGREST::keggList("pathway"))
  )

  MOD_NM <- data.frame(
    term = names(KEGGREST::keggList("module")),
    name = as.character(KEGGREST::keggList("module"))
  )

  saveRDS(PATH2KO, file.path(cache_dir, "Figure3b_PATH2KO.rds"))
  saveRDS(MOD2KO, file.path(cache_dir, "Figure3b_MOD2KO.rds"))
  saveRDS(PATH_NM, file.path(cache_dir, "Figure3b_PATH_NM.rds"))
  saveRDS(MOD_NM, file.path(cache_dir, "Figure3b_MOD_NM.rds"))
}

top_markers <- read.csv(top_marker_csv)
markers_df <- read_markers(top_markers, alias_long)

KO_by_gene <- unique(na.omit(as.data.table(gtf_map)[, .(gene_id, ko)]))
clusters <- sort(unique(markers_df$cluster))

all_path <- list()
all_mod <- list()

for (cl in clusters) {

  subm <- markers_df[cluster == cl]

  KU <- get_interest_universe(subm, KO_by_gene, gene_universe)

  cl_dir <- file.path(enrich_dir, .safe_dir(cl))
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)

  res_p <- run_enrich(KU$interest, KU$universe, PATH2KO, PATH_NM,
                      padj_cut, data_min_count)

  res_m <- run_enrich(KU$interest, KU$universe, MOD2KO, MOD_NM,
                      padj_cut, data_min_count)

  fwrite(res_p, file.path(cl_dir, "pathway.tsv"))
  fwrite(res_m, file.path(cl_dir, "module.tsv"))

  plot_dot(head(res_p[order(res_p$p.adjust)], top_n_plot),
           paste("Pathway", cl),
           file.path(cl_dir, "pathway.png"),
           file.path(cl_dir, "pathway.pdf"))

  plot_dot(head(res_m[order(res_m$p.adjust)], top_n_plot),
           paste("Module", cl),
           file.path(cl_dir, "module.png"),
           file.path(cl_dir, "module.pdf"))

  all_path[[cl]] <- res_p
  all_mod[[cl]] <- res_m
}

fwrite(rbindlist(all_path, fill = TRUE),
       file.path(combined_dir, "all_pathways.tsv"))

fwrite(rbindlist(all_mod, fill = TRUE),
       file.path(combined_dir, "all_modules.tsv"))
