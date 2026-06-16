############################################################
## Figure 4j host enterocyte bile-acid signalling DotPlot

############################################################

# ============================================================
# 1. Load packages
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# ============================================================
# 2. Input and output paths
# ============================================================
base_dir   <- "Figure4"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

input_rds <- file.path(input_dir, "whole_final_annotation.rds")

output_pdf <- file.path(output_dir, "Fig4j_host_enterocyte_bile_acid_signalling_dotplot.pdf")
output_png <- file.path(output_dir, "Fig4j_host_enterocyte_bile_acid_signalling_dotplot.png")
source_csv <- file.path(output_dir, "Fig4j_host_enterocyte_bile_acid_signalling_dotplot_source_data.csv")
gene_check_csv <- file.path(output_dir, "Fig4j_host_enterocyte_bile_acid_signalling_gene_check.csv")
celltype_count_csv <- file.path(output_dir, "Fig4j_host_celltype_counts_after_initial_filter.csv")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1000)
options(future.globals.maxSize = 2500 * 1024^2)

# ============================================================
# 3. Prepare enterocyte object
# ============================================================
seu <- readRDS(input_rds)

seu <- subset(
  seu,
  subset = tissue %in% c("cecum", "colon") & mice == "SPF"
)

sample_q <- CreateSeuratObject(
  counts = GetAssayData(seu, assay = "RNA", layer = "counts"),
  meta.data = seu[[]],
  min.cells = 50,
  min.features = 50
)

celltype_counts <- table(sample_q$celltypes)

write.csv(
  data.frame(
    celltypes = names(celltype_counts),
    cell_number = as.integer(celltype_counts)
  ),
  celltype_count_csv,
  row.names = FALSE
)

valid_celltypes <- names(celltype_counts[celltype_counts >= 100])

sample_q <- subset(
  sample_q,
  subset = celltypes %in% valid_celltypes
)

sample_q <- SCTransform(
  sample_q,
  vst.flavor = "v2",
  vars.to.regress = "nCount_RNA",
  conserve.memory = TRUE
)

sample_q$tissue_celltype <- paste(
  sample_q$tissue,
  sample_q$celltypes,
  sep = "_"
)

Idents(sample_q) <- sample_q$tissue_celltype

ent <- subset(sample_q, subset = celltypes == "Enterocyte")

# ============================================================
# 4. Define Figure 4j genes and pathways
# ============================================================
pathway_levels <- c("ASBT-IBABP-OSTalpha/beta", "FXR-FGF15", "TGR5-GLP-1")

gene_pathway_df <- data.frame(
  gene = c(
    "Slc10a2", "Fabp6", "Slc51a", "Slc51b",
    "Nr1h4", "Fgf15", "Fgfr4", "Klb", "Cyp7a1", "Nr0b2",
    "Gpb1r", "Gcg", "Pcsk1", "Glp1r", "Ins1", "Ins2"
  ),
  pathway = c(
    rep(pathway_levels[1], 4),
    rep(pathway_levels[2], 6),
    rep(pathway_levels[3], 6)
  ),
  stringsAsFactors = FALSE
)

gene_check <- gene_pathway_df %>%
  mutate(Present_in_enterocyte_object = gene %in% rownames(ent))

write.csv(gene_check, gene_check_csv, row.names = FALSE)

filtered_df <- gene_pathway_df %>%
  filter(gene %in% rownames(ent))

if (nrow(filtered_df) == 0) {
  stop("None of the Figure 4j genes were found in the enterocyte object.")
}

# ============================================================
# 5. Generate DotPlot and source data
# ============================================================
dot_data <- DotPlot(
  ent,
  features = filtered_df$gene
)$data %>%
  left_join(filtered_df, by = c("features.plot" = "gene"))

dot_data$pathway <- factor(
  dot_data$pathway,
  levels = pathway_levels
)

write.csv(dot_data, source_csv, row.names = FALSE)

# ============================================================
# 6. Plot Figure 4j
# ============================================================
fig4j <- ggplot(dot_data, aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  facet_wrap(~pathway, scales = "free_x", nrow = 1) +
  scale_color_gradient(low = "lightgrey", high = "#B2182B") +
  scale_size(range = c(1, 6)) +
  RotatedAxis() +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank()
  )

ggsave(output_pdf, fig4j, width = 8, height = 2.6)
ggsave(output_png, fig4j, width = 8, height = 2.6, dpi = 300, bg = "white")

