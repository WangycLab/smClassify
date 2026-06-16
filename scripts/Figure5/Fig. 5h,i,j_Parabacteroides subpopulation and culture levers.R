############################################################
## Figure 5h,i,j Parabacteroides subpopulation and culture levers
##

# ============================================================
# 1. Load packages
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(pheatmap)
  library(reshape2)
  library(stringr)
})

# ============================================================
# 2. Input and output paths
# ============================================================
base_dir   <- "Figure5"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

annotated_seurat_rds <- "Figure3/Output/Figure3_mouse_bacterial_functional_cluster_annotated.rds"
parabac_type_csv     <- file.path(input_dir, "Parabac_type.csv")

dpi_out <- 300

# ============================================================
# 3. Load annotated Seurat object
# ============================================================
seu <- readRDS(annotated_seurat_rds)

# ============================================================
# 4. Subset Parabacteroides cells and recluster
# ============================================================
sample <- seu
sample@active.ident <- as.factor(sample$Tax_Genus)

Parabac <- subset(sample, idents = "Parabacteroides")

write.csv(
  as.data.frame(table(Parabac$species)),
  file = file.path(output_dir, "Figure5hij_Parabacteroides_species_cell_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(Parabac$group)),
  file = file.path(output_dir, "Figure5hij_Parabacteroides_group_cell_counts.csv"),
  row.names = FALSE
)

Parabac <- SCTransform(
  Parabac,
  verbose = FALSE,
  vars.to.regress = "nCount_RNA"
)

Parabac <- RunPCA(Parabac, verbose = FALSE)
Parabac <- FindNeighbors(Parabac, dims = 1:5)
Parabac <- FindClusters(Parabac, resolution = 0.15)
Parabac <- RunUMAP(Parabac, dims = 1:5)

saveRDS(
  Parabac,
  file = file.path(output_dir, "Figure5hij_Parabacteroides_reclustered.rds")
)

# ============================================================
# 5. Supplementary UMAPs before subtype annotation
# ============================================================
cluster_cols <- c(
  "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1",
  "#3CB371", "#9ACD32", "#8B0000", "#FFD700", "#DC143C",
  "#228B22", "#FF6347", "#483D8B", "#BDB76B", "#FF1493",
  "#FF4500", "#32CD32", "#3E8E41", "#20B2AA"
)

p_cluster <- DimPlot(
  Parabac,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "seurat_clusters"
) +
  scale_color_manual(values = cluster_cols) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )
p_cluster
ggsave(
  file.path(output_dir, "Supplementary_Parabacteroides_UMAP_by_cluster.png"),
  p_cluster,
  width = 4,
  height = 3,
  dpi = dpi_out
)

p_phenotype <- DimPlot(
  Parabac,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "phenotype"
) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )

ggsave(
  file.path(output_dir, "Supplementary_Parabacteroides_UMAP_by_phenotype.png"),
  p_phenotype,
  width = 4,
  height = 3,
  dpi = dpi_out
)

p_region <- DimPlot(
  Parabac,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.1,
  group.by = "region"
) +
  scale_color_manual(values = c(
    "Cecum" = "#1f77b4",
    "Colon" = "#2ca02c",
    "Rectum" = "#ff7f0e"
  )) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )

ggsave(
  file.path(output_dir, "Supplementary_Parabacteroides_UMAP_by_region.png"),
  p_region,
  width = 4.5,
  height = 3,
  dpi = dpi_out
)

# ============================================================
# 6. Figure 5h: UMAP colored by species
# ============================================================
species_cols <- c(
  "#BCBD22",
  "#228B22",
  "#7F7F7F"
)

p_5h <- DimPlot(
  Parabac,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  group.by = "species"
) +
  scale_color_manual(values = species_cols) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )

ggsave(
  file.path(output_dir, "Figure5h_Parabacteroides_UMAP_by_species.png"),
  p_5h,
  width = 6,
  height = 3,
  dpi = dpi_out
)

ggsave(
  file.path(output_dir, "Figure5h_Parabacteroides_UMAP_by_species.pdf"),
  p_5h,
  width = 6,
  height = 3
)

# ============================================================
# 7. Supplementary BCAA-related gene dotplot
# ============================================================
p_bcaa <- DotPlot(
  Parabac,
  features = c("mutB", "epi"),
  group.by = "group"
) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(face = "bold")
  ) +
  labs(title = "BCAA-related genes - Parabacteroides distasonis")

ggsave(
  file.path(output_dir, "Supplementary_Parabacteroides_BCAA_gene_dotplot.png"),
  p_bcaa,
  width = 5,
  height = 4,
  dpi = dpi_out
)

# ============================================================
# 8. Marker analysis and subtype annotation
# ============================================================
Idents(Parabac) <- Parabac$seurat_clusters

message("Finding markers ...")
markers <- FindAllMarkers(
  object = Parabac,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.15
)

write.csv(
  markers,
  file = file.path(output_dir, "Supplementary_Parabacteroides_markers.csv"),
  row.names = FALSE
)

topN <- markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 50) |>
  ungroup()

write.csv(
  topN,
  file = file.path(output_dir, "Supplementary_Parabacteroides_markers_top50.csv"),
  row.names = FALSE
)

parabac_type <- read.csv(parabac_type_csv)

new.cluster.ids <- as.character(parabac_type[, 2])

Parabac@active.ident <- as.factor(Parabac$seurat_clusters)
names(new.cluster.ids) <- levels(Parabac)

Parabac <- RenameIdents(Parabac, new.cluster.ids)
Parabac$cellsubtype <- Parabac@active.ident

write.csv(
  as.data.frame(table(Parabac$cellsubtype)),
  file = file.path(output_dir, "Figure5i_Parabacteroides_subtype_cell_counts.csv"),
  row.names = FALSE
)

saveRDS(
  Parabac,
  file = file.path(output_dir, "Figure5hi_Parabacteroides_subtype_annotated.rds")
)

# ============================================================
# 9. Figure 5i: UMAP colored by subtype
# ============================================================
subtype_cols <- c(
  "#6A5ACD",
  "#F5A0A1",
  "#3CB371",
  "#FFD700"
)

p_5i <- DimPlot(
  Parabac,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  group.by = "cellsubtype"
) +
  scale_color_manual(values = subtype_cols) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14)
  )

ggsave(
  file.path(output_dir, "Figure5i_Parabacteroides_UMAP_by_subtype.png"),
  p_5i,
  width = 7,
  height = 3,
  dpi = dpi_out
)

ggsave(
  file.path(output_dir, "Figure5i_Parabacteroides_UMAP_by_subtype.pdf"),
  p_5i,
  width = 7,
  height = 3
)

# ============================================================
# 10. Supplementary top marker-gene dotplot
# ============================================================
Idents(Parabac) <- "seurat_clusters"

top5_nonMGYG <- topN %>%
  filter(!grepl("^MGYG", gene)) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  ungroup()

write.csv(
  top5_nonMGYG,
  file = file.path(output_dir, "Supplementary_Parabacteroides_top5_nonMGYG_markers.csv"),
  row.names = FALSE
)

marker_genes <- unique(top5_nonMGYG$gene)

p_marker <- DotPlot(Parabac, features = marker_genes) +
  scale_color_gradient(low = "#2ECC71", high = "#8E44AD") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "Functional Subtypes", title = "DotPlot of Top Marker Genes") +
  scale_size_continuous(range = c(1, 8))

ggsave(
  file.path(output_dir, "Supplementary_Parabacteroides_top_marker_gene_dotplot.pdf"),
  p_marker,
  width = 12,
  height = 3
)

ggsave(
  file.path(output_dir, "Supplementary_Parabacteroides_top_marker_gene_dotplot.png"),
  p_marker,
  width = 12,
  height = 3,
  dpi = dpi_out
)

# ============================================================
# 11. Figure 5j: culture-related levers
# ============================================================

# 11.1 Culture module gene sets
culture_modules <- list(
  Inulin_Amylose = c("susC", "susD", "malQ", "glgA"),
  Mucin_NAG = c("nanA", "nagA", "nagB", "nagC"),
  Fucose = c("fucI", "fucO"),
  Xylose_Man_Gal = c("xylA", "xylB", "xylE", "manA", "galE", "galK"),
  Vitamin_B12 = c("cobB", "cobD", "cbiA", "cbiD"),
  Fumarate_acceptor = c("frdA", "frdB"),
  Bile_salts = c("cbh", "acrF", "tolC", "ompR"),
  Branched_SCFAs = c("ilvE", "leuA", "leuB", "leuC"),
  Vitamins_broad = c(
    "bioA", "bioB", "bioC", "bioD", "bioF",
    "thiC", "thiD", "thiE",
    "folB", "folC", "folD",
    "nadA", "nadB", "nadC", "nadD", "nadE",
    "menA", "menB", "menD", "menE", "menF"
  )
)

# 11.2 Gene-evidence keywords
gene_keywords <- list(
  Mucin_NAG = c("nagA", "nagB", "nanA"),
  Fucose = c("fucI", "fucP"),
  XMG = c("xylA", "xylB", "manA", "galE", "galK"),
  Prop_core = c("mutA", "mutB", "pccA", "pccB", "mcmA", "mcmB", "prpB", "prpC", "prpD"),
  B12_core = c("cobB", "cobD", "cbiA", "cbiD")
)

# 11.3 Lever-to-module weights
lever_module_weights <- list(
  "Mucin/NAG" = c(PUL_sugars = 1.0),
  "Fucose" = c(PUL_sugars = 1.0),
  "Xylose/Man/Gal" = c(PUL_sugars = 1.0),
  "Fumarate (e- acceptor)" = c(Propionate = 1.0),
  "Ethanolamine" = c(Ethanolamine_B12 = 1.0),
  "Vitamin B12" = c(Ethanolamine_B12 = 1.0),
  "Bile salts (tolerate)" = c(Bile_tolerance = 1.0)
)

# 11.4 Collapse into six module sets for AddModuleScore
module_sets <- list(
  PUL_sugars = unique(unlist(c(
    culture_modules$Inulin_Amylose,
    culture_modules$Mucin_NAG,
    culture_modules$Fucose,
    culture_modules$Xylose_Man_Gal
  ))),
  Ethanolamine_B12 = unique(unlist(c(
    culture_modules$Vitamin_B12
  ))),
  Propionate = culture_modules$Fumarate_acceptor,
  Bile_tolerance = culture_modules$Bile_salts,
  AA_bios_deg = culture_modules$Branched_SCFAs,
  Vitamins_cofactors = culture_modules$Vitamins_broad
)

# 11.5 Add module scores
DefaultAssay(Parabac) <- DefaultAssay(Parabac)

for (nm in names(module_sets)) {
  Parabac <- AddModuleScore(
    Parabac,
    features = list(module_sets[[nm]]),
    name = paste0(nm, "_score"),
    vst = FALSE,
    assay = DefaultAssay(Parabac),
    slot = "data",
    nbin = 6
  )
}

# 11.6 Aggregate module scores by species
scores_species <- Parabac@meta.data %>%
  group_by(species) %>%
  summarise(across(ends_with("_score1"), mean, na.rm = TRUE)) %>%
  arrange(species)

mat_mod <- as.matrix(scores_species[, -1])
rownames(mat_mod) <- scores_species$species
colnames(mat_mod) <- gsub("_score1$", "", colnames(mat_mod))

z_mod <- scale(mat_mod)

# 11.7 Gene-level evidence
avg_expr_by_species <- function(gene_set) {
  feats <- unique(unlist(gene_set))
  feats <- feats[feats %in% rownames(Parabac)]
  if (length(feats) == 0) return(NULL)
  
  avg <- AverageExpression(
    Parabac,
    features = feats,
    group.by = "species",
    assays = DefaultAssay(Parabac)
  )[[1]]
  
  rowMeans(avg, na.rm = TRUE)
}

gene_scores <- list()

for (n in names(gene_keywords)) {
  sc <- avg_expr_by_species(gene_keywords[[n]])
  if (!is.null(sc)) {
    names(sc) <- gsub("-", "_", names(sc))
    gene_scores[[n]] <- sc
  }
}

if (length(gene_scores) > 0) {
  gene_df <- do.call(cbind, lapply(gene_scores, function(x) {
    x <- x[rownames(z_mod)]
    return(x)
  }))
  colnames(gene_df) <- names(gene_scores)
  z_gene <- scale(gene_df)
} else {
  z_gene <- NULL
}

# 11.8 Combine module-score and gene-expression evidence into lever scores
lever_score <- matrix(0, nrow = nrow(z_mod), ncol = length(lever_module_weights))
rownames(lever_score) <- rownames(z_mod)
colnames(lever_score) <- names(lever_module_weights)

for (lv in names(lever_module_weights)) {
  w <- lever_module_weights[[lv]]
  mod_cols <- intersect(names(w), colnames(z_mod))
  
  mod_part <- 0
  if (length(mod_cols) > 0) {
    mod_part <- as.matrix(z_mod[, mod_cols, drop = FALSE]) %*%
      matrix(w[mod_cols], ncol = 1)
    mod_part <- as.numeric(mod_part)
  }
  
  gene_part <- rep(0, nrow(z_mod))
  if (!is.null(z_gene)) {
    if (lv == "Mucin/NAG" && "Mucin_NAG" %in% colnames(z_gene)) {
      gene_part <- z_gene[, "Mucin_NAG"]
    }
    if (lv == "Fucose" && "Fucose" %in% colnames(z_gene)) {
      gene_part <- z_gene[, "Fucose"]
    }
    if (lv == "Xylose/Man/Gal" && "XMG" %in% colnames(z_gene)) {
      gene_part <- z_gene[, "XMG"]
    }
    if (lv == "Fumarate (e- acceptor)" && "Prop_core" %in% colnames(z_gene)) {
      gene_part <- z_gene[, "Prop_core"]
    }
    if (lv %in% c("Ethanolamine", "Vitamin B12") && "B12_core" %in% colnames(z_gene)) {
      gene_part <- z_gene[, "B12_core"]
    }
  }
  
  gene_part[is.na(gene_part)] <- 0
  
  if (all(gene_part == 0)) {
    lever_score[, lv] <- mod_part
  } else {
    lever_score[, lv] <- 0.5 * mod_part + 0.5 * gene_part
  }
}

# Derived column: avoid bile is the negative of bile tolerance.
lever_score <- cbind(
  lever_score,
  "Avoid bile" = -lever_score[, "Bile salts (tolerate)"]
)

# 11.9 Display label mapping
lever_display_names <- c(
  "Mucin/NAG" = "Mucin/GlcNAc utilization",
  "Fucose" = "Fucose utilization",
  "Xylose/Man/Gal" = "Xylose/Mannose/Galactose use",
  "Fumarate (e- acceptor)" = "Fumarate as electron acceptor",
  "Ethanolamine" = "Ethanolamine use (B12-linked)",
  "Vitamin B12" = "B12 synthesis/requirement",
  "Bile salts (tolerate)" = "Bile salt tolerance",
  "Avoid bile" = "Avoid bile in medium"
)

module_display_names <- c(
  "PUL_sugars" = "PUL & simple-sugar use",
  "Ethanolamine_B12" = "Ethanolamine-B12 module",
  "Propionate" = "Propionate / fumarate reduction",
  "Bile_tolerance" = "Bile efflux & stress",
  "AA_bios_deg" = "Branched SCFAs (AA metabolism)",
  "Vitamins_cofactors" = "Vitamin/cofactor synthesis"
)

gene_layer_display <- c(
  "Mucin_NAG" = "Mucin/GlcNAc gene evidence",
  "Fucose" = "Fucose gene evidence",
  "XMG" = "Xylose/Mannose/Galactose gene evidence",
  "Prop_core" = "Propionate core gene evidence",
  "B12_core" = "B12 gene evidence"
)

colnames(lever_score) <- ifelse(
  colnames(lever_score) %in% names(lever_display_names),
  lever_display_names[colnames(lever_score)],
  colnames(lever_score)
)

if (exists("z_mod")) {
  colnames(z_mod) <- ifelse(
    colnames(z_mod) %in% names(module_display_names),
    module_display_names[colnames(z_mod)],
    colnames(z_mod)
  )
}

if (exists("z_gene")) {
  colnames(z_gene) <- ifelse(
    colnames(z_gene) %in% names(gene_layer_display),
    gene_layer_display[colnames(z_gene)],
    colnames(z_gene)
  )
}

# 11.10 Symbols and heatmap
to_symbol <- function(x, pos_thr = 0.5, neg_thr = -0.5) {
  ifelse(x > pos_thr, "+", ifelse(x < neg_thr, "-", "."))
}

sym_mat <- apply(lever_score, 2, to_symbol)

sym_long <- reshape2::melt(
  sym_mat,
  varnames = c("species", "lever"),
  value.name = "sym"
)

num_long <- reshape2::melt(
  lever_score,
  varnames = c("species", "lever"),
  value.name = "score"
) %>%
  left_join(sym_long, by = c("species", "lever"))

write.csv(
  num_long,
  file = file.path(output_dir, "Figure5j_Parabacteroides_culture_lever_heatmap_source_data.csv"),
  row.names = FALSE
)

p_5j <- ggplot(num_long, aes(x = lever, y = species, fill = score)) +
  geom_tile() +
  geom_text(aes(label = sym), size = 6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(
    title = "Culture levers for isolation (+ favorable / - unfavorable / . neutral)",
    x = "Culture levers (actionable medium conditions)",
    y = "Parabacteroides species"
  ) +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_5j 
ggsave(
  file.path(output_dir, "Figure5j_Parabacteroides_culture_lever_evidence_heatmap.pdf"),
  p_5j,
  width = 8,
  height = 3.5
)

ggsave(
  file.path(output_dir, "Figure5j_Parabacteroides_culture_lever_evidence_heatmap.png"),
  p_5j,
  width = 8,
  height = 3.5,
  dpi = dpi_out
)

write.csv(
  as.data.frame(lever_score),
  file = file.path(output_dir, "Figure5j_Parabacteroides_species_culture_scores.csv")
)

write.csv(
  as.data.frame(sym_mat),
  file = file.path(output_dir, "Figure5j_Parabacteroides_species_culture_symbols.csv")
)

saveRDS(
  Parabac,
  file = file.path(output_dir, "Figure5hij_Parabacteroides_final_object.rds")
)