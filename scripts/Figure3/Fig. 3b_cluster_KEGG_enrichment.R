# ============================================================
# Figure 3b KEGG enrichment analysis of functional clusters
#
# This script performs KEGG pathway over-representation analysis

# Required input files:
#   Figure3/Input/Total cells_markers_top50.csv
#   Figure3/Input/Figure3_mouse_selected_genomes.gtf
#   Figure3/Input/Figure3_KEGG_pathway_hierarchy_clean.csv
#

# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(KEGGREST)
  library(stringr)
  library(readr)
  library(tidyr)
})

options(stringsAsFactors = FALSE)
set.seed(1)

# ============================================================
# 1. Input, output, and analysis parameters
# ============================================================
base_dir   <- "Figure3"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

fig3b_dir  <- file.path(output_dir, "Figure3b_KEGG_enrichment")
cache_dir  <- file.path(fig3b_dir, "cache")
enrich_dir <- file.path(fig3b_dir, "cluster_enrichment")
combined_dir <- file.path(enrich_dir, "combined")
log_dir <- file.path(fig3b_dir, "logs")

invisible(lapply(
  list(input_dir, output_dir, fig3b_dir, cache_dir, enrich_dir, combined_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# Input files. Keep these names consistent with the current Figure 3 workflow.
gtf_file <- file.path(input_dir, "Figure3_mouse_Mgnify_genomes.gtf")
top_marker_csv  <- file.path(input_dir, "Total cells_markers_top50.csv")
hierarchy_csv   <- file.path(input_dir, "Figure3_KEGG_pathway_hierarchy_clean.csv")


# Enrichment and plotting parameters. These are kept from the original script.
padj_cut   <- 0.25
data_min_count  <- 1
top_n_plot <- 20

stopifnot(file.exists(gtf_file))
stopifnot(file.exists(top_marker_csv))
stopifnot(file.exists(hierarchy_csv))

############################################################
# 1. Helper functions
############################################################

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
  if (any(hit)) {
    ans[hit] <- sub(p, "\\2", regmatches(x, m)[hit], perl = TRUE)
  }
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

  if (!"cluster" %in% names(DT)) {
    stop("Marker file must contain a cluster column. Current columns are: ",
         paste(names(DT), collapse = ", "))
  }

  if ("gene" %in% names(DT)) {
    DT[, gene_input := gene]
  } else if ("gene_id" %in% names(DT)) {
    DT[, gene_input := gene_id]
  } else if ("symbol" %in% names(DT)) {
    DT[, gene_input := symbol]
  } else {
    stop("Marker file must contain one of these columns: gene, gene_id, or symbol. Current columns are: ",
         paste(names(DT), collapse = ", "))
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
  if (length(universe) < 10) universe <- unique(KO_by_gene$ko)
  list(interest = interest, universe = universe)
}

run_enrich <- function(interest, universe, TERM2GENE, TERM2NAME,
                       padj_cut = 0.25, min_count = 1) {
  space <- unique(TERM2GENE$gene)
  interest <- intersect(interest, space)
  universe <- intersect(universe, space)

  if (!length(interest)) {
    return(list(all = data.frame(), keep = data.frame()))
  }

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

  if (is.null(ek)) {
    return(list(all = data.frame(), keep = data.frame()))
  }

  df <- as.data.frame(ek)
  keep <- dplyr::filter(df, p.adjust <= padj_cut, Count >= min_count)
  list(all = df, keep = keep)
}

plot_dot <- function(df, title, png_file, pdf_file) {
  if (!nrow(df)) return(invisible(NULL))
  df$Description[is.na(df$Description)] <- df$ID[is.na(df$Description)]

  p <- ggplot(df, aes(x = p.adjust, y = reorder(Description, -p.adjust))) +
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
  col.names = c("feature", "attributes"),
  quote = "",
  fill = TRUE,
  showProgress = TRUE
)

A <- as.character(DT$attributes)

lt  <- get_attr_gtf(A, "locus_tag")
nm  <- get_attr_gtf(A, "Name")
gn  <- get_attr_gtf(A, "gene")
prd <- get_attr_gtf(A, "product")

keg1 <- get_attr_gtf(A, "KEGG")
keg2 <- get_attr_gtf(A, "Dbxref")
keg3 <- ifelse(
  grepl("ko:K\\d{5}", A, ignore.case = TRUE),
  str_extract(A, "ko:K\\d{5}(?:[,;|\\s]ko:K\\d{5})*"),
  NA_character_
)
keg <- ifelse(!.blankish(keg1), keg1, ifelse(!.blankish(keg2), keg2, keg3))

gene_id <- .norm_id(lt)

sym_raw <- nm
sym_raw[.blankish(sym_raw)] <- gn[.blankish(sym_raw)]
sym_raw[.blankish(sym_raw)] <- lt[.blankish(sym_raw)]

prod_tok <- ifelse(!is.na(prd), sub("[,;].*$", "", prd), NA_character_)
prod_tok[nchar(prod_tok) > 40] <- NA_character_
sym_raw[.blankish(sym_raw)] <- prod_tok[.blankish(sym_raw)]

symbol <- .norm_id(sym_raw)

gtf_map_wide <- data.table(
  gene_id = gene_id,
  symbol = symbol,
  KEGG = keg
)
gtf_map_wide <- gtf_map_wide[!.blankish(gene_id)]

gtf_map <- gtf_map_wide[, {
  ks <- split_kos(KEGG)
  if (length(ks) == 0) list(ko = NA_character_) else list(ko = ks)
}, by = .(gene_id, symbol)]

setkeyv(gtf_map, c("gene_id", "symbol", "ko"))
gtf_map <- unique(gtf_map)

alias_long <- unique(rbind(
  data.table(alias = gtf_map$gene_id, gene_id = gtf_map$gene_id),
  data.table(alias = gtf_map$symbol,  gene_id = gtf_map$gene_id)
))
alias_long <- alias_long[!is.na(alias) & nzchar(alias)]

KO_T2G <- unique(na.omit(gtf_map[, .(term = ko, gene = gene_id)]))
gene_universe <- unique(KO_T2G$gene)

saveRDS(gtf_map,       file.path(cache_dir, "Figure3b_gtf_map.rds"))
saveRDS(alias_long,    file.path(cache_dir, "Figure3b_alias_long.rds"))
saveRDS(KO_T2G,        file.path(cache_dir, "Figure3b_KO_T2G.rds"))
saveRDS(gene_universe, file.path(cache_dir, "Figure3b_gene_universe.rds"))



kegg_cache_files <- c(
  "Figure3b_PATH2KO.rds",
  "Figure3b_MOD2KO.rds",
  "Figure3b_PATH_NM.rds",
  "Figure3b_MOD_NM.rds"
)
kegg_cache_ok <- all(file.exists(file.path(cache_dir, kegg_cache_files)))

if (kegg_cache_ok) {
  PATH2KO <- readRDS(file.path(cache_dir, "Figure3b_PATH2KO.rds"))
  MOD2KO  <- readRDS(file.path(cache_dir, "Figure3b_MOD2KO.rds"))
  PATH_NM <- readRDS(file.path(cache_dir, "Figure3b_PATH_NM.rds"))
  MOD_NM  <- readRDS(file.path(cache_dir, "Figure3b_MOD_NM.rds"))
} else {
  message("[Figure 3b] Downloading KEGG pathway/module mappings for the first run.")

  lk_path <- KEGGREST::keggLink("pathway", "ko")
  lk_mod  <- KEGGREST::keggLink("module",  "ko")

  PATH2KO <- if (length(lk_path)) {
    data.table(
      term = sub("^ko", "map", sub("^path:", "", unname(lk_path))),
      gene = sub("^ko:", "", names(lk_path))
    ) |> unique()
  } else {
    data.table(term = character(), gene = character())
  }

  MOD2KO <- if (length(lk_mod)) {
    data.table(
      term = sub("^(md:)?", "", unname(lk_mod)),
      gene = sub("^ko:", "", names(lk_mod))
    ) |> unique()
  } else {
    data.table(term = character(), gene = character())
  }

  lst_path <- KEGGREST::keggList("pathway")
  lst_mod  <- KEGGREST::keggList("module")

  PATH_NM <- data.frame(
    term = sub("^path:", "", names(lst_path)),
    name = as.character(lst_path)
  )
  MOD_NM <- data.frame(
    term = sub("^(md:)?", "", names(lst_mod)),
    name = as.character(lst_mod)
  )

  saveRDS(PATH2KO, file.path(cache_dir, "Figure3b_PATH2KO.rds"))
  saveRDS(MOD2KO,  file.path(cache_dir, "Figure3b_MOD2KO.rds"))
  saveRDS(PATH_NM, file.path(cache_dir, "Figure3b_PATH_NM.rds"))
  saveRDS(MOD_NM,  file.path(cache_dir, "Figure3b_MOD_NM.rds"))
}


############################################################
# 4. Read marker genes
############################################################

top_markers <- read.csv(top_marker_csv, check.names = FALSE)
markers_df <- read_markers(top_markers, alias_long)

KO_by_gene <- unique(na.omit(as.data.table(gtf_map)[, .(gene_id, ko)]))
clusters <- sort(unique(markers_df$cluster))


############################################################
# 5. Cluster-wise KEGG pathway/module enrichment
############################################################

all_path_all  <- list()
all_path_keep <- list()
all_mod_all   <- list()
all_mod_keep  <- list()

for (cl in clusters) {
  message("[Figure 3b] Processing cluster: ", cl)

  subm <- markers_df[cluster == cl, ]
  KU <- get_interest_universe(subm, KO_by_gene, gene_universe)

  cl_dir <- file.path(enrich_dir, paste0("cluster_", .safe_dir(cl)))
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)

  # KEGG pathway ORA.
  res_p <- run_enrich(
    KU$interest,
    KU$universe,
    PATH2KO,
    PATH_NM,
    padj_cut = padj_cut,
    min_count = data_min_count
  )

  fwrite(res_p$all,  file.path(cl_dir, "Figure3b_KEGG_pathway_results_all.tsv"),      sep = "\t", quote = FALSE, na = "")
  fwrite(res_p$keep, file.path(cl_dir, "Figure3b_KEGG_pathway_results_filtered.tsv"), sep = "\t", quote = FALSE, na = "")

  if (nrow(res_p$keep)) {
    dfp <- res_p$keep %>% arrange(p.adjust) %>% dplyr::slice_head(n = top_n_plot)
    plot_dot(
      dfp,
      paste("KEGG Pathway -", cl),
      file.path(cl_dir, "Figure3b_KEGG_pathway_dotplot.png"),
      file.path(cl_dir, "Figure3b_KEGG_pathway_dotplot.pdf")
    )
  }

  if (nrow(res_p$all))  res_p$all$cluster  <- cl
  if (nrow(res_p$keep)) res_p$keep$cluster <- cl

  all_path_all[[as.character(cl)]]  <- res_p$all
  all_path_keep[[as.character(cl)]] <- res_p$keep

  # KEGG module ORA.
  res_m <- run_enrich(
    KU$interest,
    KU$universe,
    MOD2KO,
    MOD_NM,
    padj_cut = padj_cut,
    min_count = data_min_count
  )

  fwrite(res_m$all,  file.path(cl_dir, "Figure3b_KEGG_module_results_all.tsv"),      sep = "\t", quote = FALSE, na = "")
  fwrite(res_m$keep, file.path(cl_dir, "Figure3b_KEGG_module_results_filtered.tsv"), sep = "\t", quote = FALSE, na = "")

  if (nrow(res_m$keep)) {
    dfm <- res_m$keep %>% arrange(p.adjust) %>% dplyr::slice_head(n = top_n_plot)
    plot_dot(
      dfm,
      paste("KEGG Module -", cl),
      file.path(cl_dir, "Figure3b_KEGG_module_dotplot.png"),
      file.path(cl_dir, "Figure3b_KEGG_module_dotplot.pdf")
    )
  }

  if (nrow(res_m$all))  res_m$all$cluster  <- cl
  if (nrow(res_m$keep)) res_m$keep$cluster <- cl

  all_mod_all[[as.character(cl)]]  <- res_m$all
  all_mod_keep[[as.character(cl)]] <- res_m$keep
}

############################################################
# 6. Export combined enrichment tables
############################################################

comb_dir <- file.path(enrich_dir, "combined")
dir.create(comb_dir, recursive = TRUE, showWarnings = FALSE)

PA_all  <- data.table::rbindlist(all_path_all,  fill = TRUE, use.names = TRUE)
PA_keep <- data.table::rbindlist(all_path_keep, fill = TRUE, use.names = TRUE)
MO_all  <- data.table::rbindlist(all_mod_all,   fill = TRUE, use.names = TRUE)
MO_keep <- data.table::rbindlist(all_mod_keep,  fill = TRUE, use.names = TRUE)

fwrite(PA_all,  file.path(comb_dir, "Figure3b_KEGG_pathway_results_all_clusters_all.tsv"),      sep = "\t", quote = FALSE, na = "")
fwrite(PA_keep, file.path(comb_dir, "Figure3b_KEGG_pathway_results_all_clusters_filtered.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(MO_all,  file.path(comb_dir, "Figure3b_KEGG_module_results_all_clusters_all.tsv"),       sep = "\t", quote = FALSE, na = "")
fwrite(MO_keep, file.path(comb_dir, "Figure3b_KEGG_module_results_all_clusters_filtered.tsv"),  sep = "\t", quote = FALSE, na = "")


df <- readr::read_tsv(
  file.path(comb_dir, "Figure3b_KEGG_pathway_results_all_clusters_all.tsv"),
  show_col_types = FALSE
)

hier <- read.csv(hierarchy_csv, check.names = FALSE)

# Merge pathway hierarchy and calculate plotting metrics.
df_anno <- df %>%
  mutate(
    ID = sub("^ko", "map", ID),
    Count = if (!"Count" %in% names(.)) as.numeric(sub("/.*$", "", GeneRatio)) else Count,
    GeneRatio_num = {
      sp <- strsplit(GeneRatio, "/", fixed = TRUE)
      as.numeric(vapply(sp, `[`, 1L, FUN.VALUE = character(1))) /
        as.numeric(vapply(sp, `[`, 2L, FUN.VALUE = character(1)))
    },
    p.adjust = suppressWarnings(as.numeric(p.adjust)),
    log10padj = ifelse(!is.na(p.adjust) & p.adjust > 0, -log10(p.adjust), NA_real_)
  ) %>%
  dplyr::left_join(
    hier %>% dplyr::select(PathwayID, Category, Subcategory),
    by = c("ID" = "PathwayID")
  ) %>%
  mutate(
    Category = tidyr::replace_na(Category, "Unclassified"),
    Subcategory = tidyr::replace_na(Subcategory, "Unclassified")
  )

# Original manual filtering rules.
remove_big_pathway <- TRUE

remove_categories <- c(
  "Human Diseases",
  "Drug Development"
)

remove_Description <- c(
  "Autophagy - animal",
  "Autophagy - yeast",
  "Mitophagy - animal",
  "Mitophagy - yeast",
  "Longevity regulating pathway - multiple species",
  "Cell cycle - Caulobacter",
  "Biofilm formation - Escherichia coli",
  "Biofilm formation - Vibrio cholerae",
  "Biofilm formation - Pseudomonas aeruginosa",
  "Drug metabolism - other enzymes",
  "ABC transporters",
  "Glycosphingolipid biosynthesis - globo and isoglobo series",
  "Streptomycin biosynthesis",
  "Photosynthesis",
  "Carbon fixation by Calvin cycle",
  "Glucosinolate biosynthesis",
  "Biosynthesis of various plant secondary metabolites",
  "Monobactam biosynthesis"
)

big_pathway_ids <- c(
  "ko01100", "map01100",
  "ko01110", "map01110",
  "ko01120", "map01120",
  "ko01200", "map01200",
  "ko01210", "map01210",
  "ko01212", "map01212",
  "ko01230", "map01230",
  "ko01232", "map01232",
  "ko01250", "map01250",
  "ko01240", "map01240",
  "ko01220", "map01220",
  "ko01310", "map01310",
  "ko01320", "map01320",
  "ko03010", "map03010"
)

norm_big_ids <- unique(sub("^ko", "map", big_pathway_ids))

df_filtered <- df_anno %>%
  dplyr::mutate(ID_norm = sub("^ko", "map", ID)) %>%
  dplyr::filter(!is.na(Description), !is.na(p.adjust)) %>%
  dplyr::mutate(
    Count = if (!"Count" %in% names(.)) as.numeric(sub("/.*$", "", GeneRatio)) else Count
  ) %>%
  { if (remove_big_pathway) dplyr::filter(., !(ID_norm %in% norm_big_ids)) else . } %>%
  dplyr::filter(!Category %in% remove_categories) %>%
  dplyr::filter(!Description %in% remove_Description) %>%
  mutate(
    GeneRatio_num = as.numeric(GeneRatio_num),
    log10p = -log10(p.adjust),
    cat_sub = paste(Subcategory, sep = " | ")
  ) %>%
  dplyr::filter(p.adjust < 0.05)

# Keep the original selection rule: each cluster keeps the 20 entries with
# the smallest GeneRatio_num.
df_top <- df_filtered %>%
  group_by(cluster) %>%
  slice_min(order_by = GeneRatio_num, n = 20) %>%
  ungroup()

cat_sub_levels <- unique(df_top$cat_sub)

kegg_sig <- df_top %>%
  mutate(
    cat_sub = factor(cat_sub, levels = cat_sub_levels),
    Description = factor(Description, levels = unique(Description)),
    cluster = as.factor(cluster)
  )

write.csv(kegg_sig, file.path(fig3b_dir, "Figure3b_kegg_sig.csv"), row.names = FALSE)
write.csv(kegg_sig, file.path(fig3b_dir, "Figure3b_KEGG_pathway_overview_source_data.csv"), row.names = FALSE)



############################################################
# 8. Draw Figure 3b
############################################################

p_fig3b <- ggplot(
  kegg_sig,
  aes(
    x = cluster,
    y = Description,
    size = log10p,
    fill = GeneRatio_num
  )
) +
  geom_point(shape = 21, colour = "black", stroke = 0.3) +
  scale_size(range = c(3, 7), name = expression(-log[10](p~value))) +
  scale_fill_gradient(low = "white", high = "firebrick4", name = "Gene Ratio") +
  facet_grid(
    rows = vars(cat_sub),
    cols = NULL,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    x = "Cluster",
    y = NULL,
    title = "KEGG Enrichment Pathways"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", colour = "grey30"),
    strip.text.y = element_text(face = "bold", angle = 0, hjust = 0),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.spacing.y = unit(3, "mm"),
    legend.position = "right"
  )

print(p_fig3b)

ggsave(
  file.path(fig3b_dir, "Figure3b_KEGG_pathway_overview.png"),
  plot = p_fig3b,
  width = 10.5,
  height = 11,
  dpi = 500
)

ggsave(
  file.path(fig3b_dir, "Figure3b_KEGG_pathway_overview.pdf"),
  plot = p_fig3b,
  width = 10.5,
  height = 11
)

