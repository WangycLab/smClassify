
suppressPackageStartupMessages({
  library(data.table)      
  library(dplyr)           
  library(ggplot2)          
  library(clusterProfiler) 
  library(KEGGREST)         
  library(stringr)         
})
options(stringsAsFactors = FALSE)
set.seed(1)
setwd("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009")
OUTDIR <- "D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009"
GTFFILE  <- "M1_selected_genome_fixed.gtf"

OUTDIR   <- get0("OUTDIR",   ifnotfound = file.path(getwd(), "ora_out"))
GTFFILE  <- get0("GTFFILE",  ifnotfound = file.path(getwd(), "example.gtf"))

# Diameters
PADJ_CUT   <- 0.25     # p.adjust 
MIN_COUNT  <- 1        # Min gene count
TOP_N_PLOT <- 20       # Top in dotplot

# Dirs
DIR_CACHE <- file.path(OUTDIR, "cache")
DIR_ENR   <- file.path(OUTDIR, "enrich_out")
DIR_LOG   <- file.path(OUTDIR, "logs")
invisible(lapply(list(OUTDIR, DIR_CACHE, DIR_ENR, DIR_LOG), dir.create, recursive = TRUE, showWarnings = FALSE))


.norm_id <- function(x) toupper(gsub("_","-", as.character(x), fixed = TRUE)) 
.blankish <- function(x){ x <- as.character(x); x[is.na(x)] <- ""; !nzchar(trimws(x)) }

# Get attributes name
get_attr_gtf <- function(x, key){
  p <- paste0("(^|[;\\t\\s])", key, '\\s+"?([^";]+)"?')
  m <- regexpr(p, x, perl = TRUE)
  ans <- rep(NA_character_, length(x)); hit <- m > 0
  if (any(hit)) ans[hit] <- sub(p, "\\2", regmatches(x, m)[hit], perl = TRUE)
  ans
}

# Split KEGG KO
split_kos <- function(s){
  if (is.na(s) || !nzchar(s)) return(character(0))
  toks <- unlist(strsplit(s, "[,;\\s|]+"))
  toks <- gsub("^(ko:|KEGG:)", "", toks, ignore.case = TRUE)
  toks <- toupper(toks)
  toks[grepl("^K\\d{5}$", toks)]
}

.safe_dir <- function(x){ gsub('[\\\\/:*?"<>|]', "_", as.character(x)) } # ��ȫ�ļ�����

# Step A: Get GTF reference

stopifnot(file.exists(GTFFILE))
DT <- fread(GTFFILE, sep = "\t", header = FALSE, select = c(3,9),
            col.names = c("feature","attributes"), quote = "", fill = TRUE, showProgress = TRUE)
A   <- as.character(DT$attributes)

lt  <- get_attr_gtf(A, "locus_tag")
nm  <- get_attr_gtf(A, "Name")
gn  <- get_attr_gtf(A, "gene")
prd <- get_attr_gtf(A, "product")

keg1 <- get_attr_gtf(A, "KEGG")
keg2 <- get_attr_gtf(A, "Dbxref")
keg3 <- ifelse(grepl("ko:K\\d{5}", A, ignore.case = TRUE),
               str_extract(A, "ko:K\\d{5}(?:[,;|\\s]ko:K\\d{5})*"), NA_character_)
keg <- ifelse(!.blankish(keg1), keg1, ifelse(!.blankish(keg2), keg2, keg3))

# Get gene symbol 
gene_id <- .norm_id(lt)
sym_raw <- nm
sym_raw[.blankish(sym_raw)] <- gn[.blankish(sym_raw)]
sym_raw[.blankish(sym_raw)] <- lt[.blankish(sym_raw)]
prod_tok <- ifelse(!is.na(prd), sub("[,;].*$","", prd), NA_character_)
prod_tok[nchar(prod_tok) > 40] <- NA_character_
sym_raw[.blankish(sym_raw)] <- prod_tok[.blankish(sym_raw)]
symbol <- .norm_id(sym_raw)


gtf_map_wide <- data.table(gene_id = gene_id, symbol = symbol, KEGG = keg)
gtf_map_wide <- gtf_map_wide[! .blankish(gene_id)]

gtf_map <- gtf_map_wide[, {
  ks <- split_kos(KEGG); if (length(ks) == 0) list(ko = NA_character_) else list(ko = ks)
}, by = .(gene_id, symbol)]
setkeyv(gtf_map, c("gene_id","symbol","ko"))
gtf_map <- unique(gtf_map)
alias_long <- unique(rbind(
  data.table(alias = gtf_map$gene_id, gene_id = gtf_map$gene_id),
  data.table(alias = gtf_map$symbol,  gene_id = gtf_map$gene_id)
))
alias_long <- alias_long[!is.na(alias) & nzchar(alias)]

# KO to gene
KO_T2G <- unique(na.omit(gtf_map[, .(term = ko, gene = gene_id)]))
gene_universe <- unique(KO_T2G$gene)

saveRDS(gtf_map,       file.path(DIR_CACHE, "gtf_map.rds"))
saveRDS(alias_long,    file.path(DIR_CACHE, "alias_long.rds"))
saveRDS(KO_T2G,        file.path(DIR_CACHE, "KO_T2G.rds"))
saveRDS(gene_universe, file.path(DIR_CACHE, "gene_universe.rds"))
message("Step A finished, GTF reference saved to ", DIR_CACHE)

## Step B: Read reference

# Read reference
gtf_map       <- readRDS(file.path(DIR_CACHE, "gtf_map.rds"))
alias_long    <- readRDS(file.path(DIR_CACHE, "alias_long.rds"))
KO_T2G        <- readRDS(file.path(DIR_CACHE, "KO_T2G.rds"))
gene_universe <- readRDS(file.path(DIR_CACHE, "gene_universe.rds"))

# KEGG Pathway/Module data
kegg_cache_ok <- all(file.exists(file.path(DIR_CACHE, c("PATH2KO.rds","MOD2KO.rds","PATH_NM.rds","MOD_NM.rds"))))
if (kegg_cache_ok){
  PATH2KO <- readRDS(file.path(DIR_CACHE, "PATH2KO.rds"))
  MOD2KO  <- readRDS(file.path(DIR_CACHE, "MOD2KO.rds"))
  PATH_NM <- readRDS(file.path(DIR_CACHE, "PATH_NM.rds"))
  MOD_NM  <- readRDS(file.path(DIR_CACHE, "MOD_NM.rds"))
} else {
  lk_path <- KEGGREST::keggLink("pathway", "ko")
  lk_mod  <- KEGGREST::keggLink("module",  "ko")
  PATH2KO <- if (length(lk_path)) data.table(
    term = sub("^path:", "", unname(lk_path)) |> sub("^ko", "map", x = _),
    gene = sub("^ko:", "", names(lk_path))
  ) |> unique() else data.table(term=character(), gene=character())
  MOD2KO <- if (length(lk_mod)) data.table(
    term = sub("^(md:)?", "", unname(lk_mod)),
    gene = sub("^ko:", "", names(lk_mod))
  ) |> unique() else data.table(term=character(), gene=character())
  lst_path <- KEGGREST::keggList("pathway")
  lst_mod  <- KEGGREST::keggList("module")
  PATH_NM  <- data.frame(term = sub("^path:", "", names(lst_path)), name = as.character(lst_path))
  MOD_NM   <- data.frame(term = sub("^(md:)?", "", names(lst_mod)),  name = as.character(lst_mod))
  saveRDS(PATH2KO, file.path(DIR_CACHE, "PATH2KO.rds"))
  saveRDS(MOD2KO,  file.path(DIR_CACHE, "MOD2KO.rds"))
  saveRDS(PATH_NM, file.path(DIR_CACHE, "PATH_NM.rds"))
  saveRDS(MOD_NM,  file.path(DIR_CACHE, "MOD_NM.rds"))
}
message(" KEGG referenced saved. ")


read_markers <- function(DT, alias_long) {
  DT <- data.table::as.data.table(DT)
  data.table::setnames(DT, tolower(names(DT)))
  DT[, gene := .norm_id(gene)]
  mapped <- merge(DT, alias_long, by.x = "gene", by.y = "alias", all.x = TRUE)
  mapped[is.na(gene_id), gene_id := gene]
  mapped[, .(gene_id, cluster)]
}


# Read top cell markers 
topN_markers<-read.csv("Total cells_markers_top50.csv")
topN_markers
markers_df <- read_markers(topN_markers, alias_long)

KO_by_gene <- unique(na.omit(as.data.table(gtf_map)[, .(gene_id, ko)]))

get_interest_universe <- function(sub_df, KO_by_gene, gene_universe){
  gids <- unique(na.omit(sub_df$gene_id))
  interest <- unique(KO_by_gene$ko[KO_by_gene$gene_id %in% gids])
  universe <- unique(KO_by_gene$ko[KO_by_gene$gene_id %in% gene_universe])
  if (length(universe) < 10) universe <- unique(KO_by_gene$ko)
  list(interest = interest, universe = universe)
}

run_enrich <- function(interest, universe, TERM2GENE, TERM2NAME, padj_cut=0.25, minCount=1){
  space <- unique(TERM2GENE$gene)
  interest <- intersect(interest, space)
  universe <- intersect(universe, space)
  if (!length(interest)) return(list(all=data.frame(), keep=data.frame()))
  ek <- clusterProfiler::enricher(
    gene = interest,
    TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
    universe = universe,
    pAdjustMethod = "BH",
    minGSSize = 1, maxGSSize = 5000,
    pvalueCutoff = 1, qvalueCutoff = 1
  )
  if (is.null(ek)) return(list(all=data.frame(), keep=data.frame()))
  df <- as.data.frame(ek)
  keep <- dplyr::filter(df, p.adjust <= padj_cut, Count >= minCount)
  list(all=df, keep=keep)
}

plot_dot <- function(df, title, png_file, pdf_file){
  if (!nrow(df)) return(invisible(NULL))
  df$Description[is.na(df$Description)] <- df$ID[is.na(df$Description)]
  p <- ggplot2::ggplot(df, aes(x = p.adjust, y = reorder(Description, -p.adjust))) +
    ggplot2::geom_point(aes(size = Count, color = p.adjust)) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(title = title, x = "p.adjust", y = NULL) +
    ggplot2::theme_bw(base_size = 12)
  ggplot2::ggsave(png_file, p, width = 8, height = 6, dpi = 180)
  ggplot2::ggsave(pdf_file, p, width = 8, height = 6)
}

## Enrichment analysis by each cluster
clusters <- sort(unique(markers_df$cluster))
clusters
all_path_all  <- list()
all_path_keep <- list()
all_mod_all   <- list()
all_mod_keep  <- list()


for (cl in clusters){
  message("Run cluster", cl)
  subm <- markers_df[cluster == cl, ]
  KU <- get_interest_universe(subm, KO_by_gene, gene_universe)
  
  cl_dir <- file.path(DIR_ENR, paste0("cluster_", .safe_dir(cl)))
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Pathway
  res_p <- run_enrich(KU$interest, KU$universe, PATH2KO, PATH_NM, PADJ_CUT, MIN_COUNT)
  fwrite(res_p$all,  file.path(cl_dir, "KEGG_pathway_results_all.tsv"), sep="\t", quote=FALSE, na="")
  fwrite(res_p$keep, file.path(cl_dir, "KEGG_pathway_results_filtered.tsv"), sep="\t", quote=FALSE, na="")
  if (nrow(res_p$keep)){
    dfp <- res_p$keep %>% arrange(p.adjust) %>% dplyr::slice_head(n = TOP_N_PLOT)
    plot_dot(dfp, paste("KEGG Pathway -", cl),
             file.path(cl_dir, "KEGG_pathway_dotplot.png"),
             file.path(cl_dir, "KEGG_pathway_dotplot.pdf"))
  }
  if (nrow(res_p$all))  res_p$all$cluster  <- cl
  if (nrow(res_p$keep)) res_p$keep$cluster <- cl
  
  all_path_all[[as.character(cl)]]  <- res_p$all
  all_path_keep[[as.character(cl)]] <- res_p$keep
  
  # Module
  res_m <- run_enrich(KU$interest, KU$universe, MOD2KO, MOD_NM, PADJ_CUT, MIN_COUNT)
  fwrite(res_m$all,  file.path(cl_dir, "KEGG_module_results_all.tsv"), sep="\t", quote=FALSE, na="")
  fwrite(res_m$keep, file.path(cl_dir, "KEGG_module_results_filtered.tsv"), sep="\t", quote=FALSE, na="")
  if (nrow(res_m$keep)){
    dfm <- res_m$keep %>% arrange(p.adjust) %>% dplyr::slice_head(n = TOP_N_PLOT)
    plot_dot(dfm, paste("KEGG Module -", cl),
             file.path(cl_dir, "KEGG_module_dotplot.png"),
             file.path(cl_dir, "KEGG_module_dotplot.pdf"))
  }
  if (nrow(res_m$all))  res_m$all$cluster  <- cl
  if (nrow(res_m$keep)) res_m$keep$cluster <- cl
  
  all_mod_all[[as.character(cl)]]  <- res_m$all
  all_mod_keep[[as.character(cl)]] <- res_m$keep
}

# Combine and output
comb_dir <- file.path(DIR_ENR, "combined")
dir.create(comb_dir, recursive = TRUE, showWarnings = FALSE)

PA_all  <- data.table::rbindlist(all_path_all,  fill = TRUE, use.names = TRUE)
PA_keep <- data.table::rbindlist(all_path_keep, fill = TRUE, use.names = TRUE)
MO_all  <- data.table::rbindlist(all_mod_all,   fill = TRUE, use.names = TRUE)
MO_keep <- data.table::rbindlist(all_mod_keep,  fill = TRUE, use.names = TRUE)

fwrite(PA_all,  file.path(comb_dir, "KEGG_pathway_results_all_clusters_all.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(PA_keep, file.path(comb_dir, "KEGG_pathway_results_all_clusters_filtered.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(MO_all,  file.path(comb_dir, "KEGG_module_results_all_clusters_all.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(MO_keep, file.path(comb_dir, "KEGG_module_results_all_clusters_filtered.tsv"), sep = "\t", quote = FALSE, na = "")

utils::capture.output(sessionInfo(), file = file.path(DIR_LOG, "sessionInfo.txt"))



library(tidyverse)
df <- readr::read_tsv("enrich_out/combined/KEGG_pathway_results_all_clusters_all.tsv", show_col_types = FALSE)

# Read kegg pathway hierarchy file
hier <- read.csv("kegg_pathway_hierarchy_clean_v4.csv")
# hier <- readr::read_csv("kegg_pathway_hierarchy_clean_v4.csv", show_col_types = FALSE)

df_anno <- df %>%
  mutate(
    ID = sub("^ko", "map", ID),
    Count = if (!"Count" %in% names(.)) as.numeric(sub("/.*$", "", GeneRatio)) else Count,
    GeneRatio_num = {
      sp <- strsplit(GeneRatio, "/", fixed = TRUE)
      as.numeric(vapply(sp, `[`, 1L, FUN.VALUE = character(1))) /
        as.numeric(vapply(sp, `[`, 2L, FUN.VALUE = character(1)))
    },
    p.adjust  = suppressWarnings(as.numeric(p.adjust)),
    log10padj = ifelse(!is.na(p.adjust) & p.adjust > 0, -log10(p.adjust), NA_real_)
  ) %>%
  dplyr::left_join(hier %>% dplyr::select(PathwayID, Category, Subcategory),
                   by = c("ID" = "PathwayID")) %>%
  mutate(
    Category    = tidyr::replace_na(Category,    "Unclassified"),
    Subcategory = tidyr::replace_na(Subcategory, "Unclassified")
  )


# Remove terms related to human diseases

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

# Remove broad pathways
big_pathway_ids <- c(
  "ko01100","map01100",
  "ko01110","map01110",
  "ko01120","map01120",
  "ko01200","map01200",
  "ko01210","map01210",
  "ko01212","map01212",
  "ko01230","map01230",
  "ko01232","map01232",
  "ko01250","map01250",
  "ko01240","map01240",
  "ko01220","map01220",
  "ko01310","map01310",
  "ko01320","map01320",
  "ko03010","map03010"
)


norm_big_ids <- unique(sub("^ko", "map", big_pathway_ids))

df_filtered <- df_anno %>%
  dplyr::mutate(ID_norm = sub("^ko", "map", ID)) %>%
  dplyr::filter(!is.na(Description), !is.na(p.adjust)) %>%
  dplyr::mutate(
    Count = if (!"Count" %in% names(.)) as.numeric(sub("/.*$", "", GeneRatio)) else Count
  ) %>%
  { if (remove_big_pathway) dplyr::filter(., !(ID_norm %in% norm_big_ids)) else . } %>%
  dplyr::filter(
    !Category %in% remove_categories,   
  )%>%
  dplyr::filter(
    !Description %in% remove_Description,   
  )%>%
  mutate(
    GeneRatio_num = as.numeric(GeneRatio_num),
    log10p = -log10(p.adjust),
    cat_sub = paste(Subcategory, sep = " | ")
  )%>%
  dplyr::filter(p.adjust <0.05)



df_top <- df_filtered %>%
  group_by(cluster) %>%
  slice_min(order_by = GeneRatio_num, n = 20) %>%
  ungroup()
df_filtered$GeneRatio_num

cat_sub_levels <- unique(df_top$cat_sub)

kegg_sig <- df_top %>%
  mutate(
    cat_sub     = factor(cat_sub, levels = cat_sub_levels),
    Description = factor(Description, levels = unique(Description))
  )
kegg_sig $GeneRatio_num


write.csv(kegg_sig, "kegg_sig.csv", row.names = FALSE)

# Plot
library(ggplot2)
library(ggh4x)

strip_colors <- c("#f0f0f0","#d9f0d3","#c6dbef","#fdd0a2") # �ɸ�

kegg_sig$cluster<-as.factor(kegg_sig$cluster)

p <- ggplot(kegg_sig, aes(x = cluster,
                          y = Description,
                          size = log10p,
                          fill = GeneRatio_num)) +
  geom_point(shape = 21, colour = "black", stroke = 0.3) +
  scale_size(range = c(3, 7), name = expression(-log[10](p~value))) +
  scale_fill_gradient(low = "white", high = "firebrick4",
                      name ="Gene Ratio" ) +
  ggh4x::facet_grid2(
    rows  = vars(cat_sub),
    cols  = NULL,
    scales = "free_y",
    space  = "free_y",
    strip = ggh4x::strip_themed(
      background_y = elem_list_rect(fill = strip_colors, colour = "grey30"),
      text_y       = elem_list_text(face = "bold", hjust = 0, angle = 0) # ������������
    )
  ) +
  labs(x = "Cluster", y = NULL, title = "KEGG Enrichment Pathways") +
  theme_bw(base_size = 12) +
  theme(
    strip.placement = "outside",
    axis.text.y     = element_text(size = 12, color = "black"),
    panel.spacing.y = unit(3, "mm"),
    legend.position = "right"
  )

print(p)

ggsave("DEGs_cluster_KEGGpathways_overview.png", plot = p, width = 10.5, height = 11, dpi = 500)
ggsave("DEGs_cluster_KEGGpathways_overview.pdf", plot = p, width = 10.5, height = 11)

