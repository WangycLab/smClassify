# ============================================================
# Figure 3b KEGG enrichment analysis (SCI minimal version)
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(KEGGREST)
library(stringr)
library(readr)
library(tidyr)

set.seed(1)

# ============================================================
# 1. paths
# ============================================================

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
top_n_plot <- 20

# ============================================================
# 2. helpers (minimal)
# ============================================================

.norm_id <- function(x) {
  toupper(gsub("_", "-", as.character(x), fixed = TRUE))
}

.blank <- function(x) {
  is.na(x) | !nzchar(trimws(as.character(x)))
}

get_attr_gtf <- function(x, key) {
  pattern <- paste0(key, '\\s+"?([^";]+)"?')
  m <- regexpr(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  hit <- m > 0
  out[hit] <- regmatches(x, m)[hit] |>
    gsub(pattern, "\\1", x[hit])
  out
}

split_kos <- function(s) {
  toks <- unlist(strsplit(as.character(s), "[,;\\s|]+"))
  toks <- gsub("^(ko:|KEGG:)", "", toks)
  toks <- toupper(toks)
  toks[grepl("^K\\d{5}$", toks)]
}

safe_name <- function(x) gsub('[\\\\/:*?"<>|]', "_", x)

# ============================================================
# 3. read GTF + build mapping
# ============================================================

DT <- fread(gtf_file, sep = "\t", header = FALSE,
            select = c(3, 9),
            col.names = c("feature", "attributes"))

A <- DT$attributes

lt  <- get_attr_gtf(A, "locus_tag")
nm  <- get_attr_gtf(A, "Name")
gn  <- get_attr_gtf(A, "gene")
prd <- get_attr_gtf(A, "product")
keg <- get_attr_gtf(A, "KEGG")

gene_id <- .norm_id(lt)

symbol <- nm
symbol[.blank(symbol)] <- gn[.blank(symbol)]
symbol[.blank(symbol)] <- lt[.blank(symbol)]
symbol[.blank(symbol)] <- prd[.blank(symbol)]

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

# ============================================================
# 4. KEGG enrichment (NO stop, NO guards)
# ============================================================

run_enrich <- function(interest, universe, TERM2GENE, TERM2NAME) {

  interest <- intersect(interest, TERM2GENE$gene)
  universe <- intersect(universe, TERM2GENE$gene)

  ek <- enricher(
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

  as.data.frame(ek)
}

plot_dot <- function(df, title, out_png, out_pdf) {

  df$Description[is.na(df$Description)] <- df$ID[is.na(df$Description)]

  p <- ggplot(df, aes(
    x = p.adjust,
    y = reorder(Description, -p.adjust)
  )) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_x_reverse() +
    labs(title = title, x = "p.adjust", y = NULL) +
    theme_bw()

  ggsave(out_png, p, width = 8, height = 6, dpi = 180)
  ggsave(out_pdf, p, width = 8, height = 6)
}

# ============================================================
# 5. KEGG databases
# ============================================================

lk_path <- keggLink("pathway", "ko")
lk_mod  <- keggLink("module", "ko")

PATH2KO <- data.table(
  term = sub("^ko:", "", names(lk_path)),
  gene = sub("^ko:", "", lk_path)
)

MOD2KO <- data.table(
  term = sub("^ko:", "", names(lk_mod)),
  gene = sub("^ko:", "", lk_mod)
)

# ============================================================
# 6. markers
# ============================================================

markers <- fread(top_marker_csv)
names(markers) <- tolower(names(markers))

markers[, gene_input :=
          if ("gene" %in% names(markers)) gene else
            if ("gene_id" %in% names(markers)) gene_id else symbol]

markers[, gene_input := .norm_id(gene_input)]

markers <- merge(
  markers[, .(gene_input, cluster)],
  alias_long,
  by.x = "gene_input",
  by.y = "alias",
  all.x = TRUE
)

markers[is.na(gene_id), gene_id := gene_input]

clusters <- unique(markers$cluster)

# ============================================================
# 7. main loop (fully minimal)
# ============================================================

all_path <- list()
all_mod <- list()

for (cl in clusters) {

  sub <- markers[cluster == cl]

  interest <- unique(KO_T2G$term[KO_T2G$gene %in% sub$gene_id])
  universe <- unique(KO_T2G$term)

  res_p <- run_enrich(interest, universe, PATH2KO, NULL)
  res_m <- run_enrich(interest, universe, MOD2KO, NULL)

  cl_dir <- file.path(enrich_dir, safe_name(cl))
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(res_p, file.path(cl_dir, "pathway.tsv"))
  fwrite(res_m, file.path(cl_dir, "module.tsv"))

  top_p <- res_p %>% arrange(p.adjust) %>% head(top_n_plot)
  top_m <- res_m %>% arrange(p.adjust) %>% head(top_n_plot)

  plot_dot(top_p,
           paste("Pathway", cl),
           file.path(cl_dir, "pathway.png"),
           file.path(cl_dir, "pathway.pdf"))

  plot_dot(top_m,
           paste("Module", cl),
           file.path(cl_dir, "module.png"),
           file.path(cl_dir, "module.pdf"))

  all_path[[cl]] <- res_p
  all_mod[[cl]] <- res_m
}

# ============================================================
# 8. export
# ============================================================

fwrite(rbindlist(all_path, fill = TRUE),
       file.path(combined_dir, "all_pathways.tsv"))

fwrite(rbindlist(all_mod, fill = TRUE),
       file.path(combined_dir, "all_modules.tsv"))

# ============================================================
# 9. final plot
# ============================================================

df <- fread(file.path(combined_dir, "all_pathways.tsv"))

hier <- fread(hierarchy_csv)

df <- merge(df, hier, by.x = "ID", by.y = "PathwayID", all.x = TRUE)

df <- df[p.adjust < 0.05]

df_top <- df[, head(.SD[order(GeneRatio)], 20), by = cluster]

p <- ggplot(df_top, aes(cluster, Description, size = -log10(p.adjust))) +
  geom_point() +
  theme_bw()

ggsave(file.path(fig3b_dir, "KEGG_overview.png"), p, width = 10, height = 10)
ggsave(file.path(fig3b_dir, "KEGG_overview.pdf"), p, width = 10, height = 10)
