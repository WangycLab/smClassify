# This script generates Figure 3g, 3h, 3i, 3j, 3k panel plots
# This script performs:
#   1) Subsetting the main Seurat object to cluster C0.
#   2) Reclustering and UMAP embedding of C0.
#   3) UMAP visualizations by:
#        - subcluster,
#        - group,
#        - phenotype (DB vs WT),
#        - location,
#        - species,
#        - family (Tax_Family).  
#   4) Module score calculation for functional gene sets
#      (protective / risk / SCFA-context modules).
#   5) Dot plots of module scores and individual genes
#   6) Cluster composition bar plots by group / phenotype / location
#
# ============================================================

suppressPackageStartupMessages({
  req_pkgs <- c(
    "Seurat", "tidyverse", "Matrix", "data.table", "ggrepel",
    "RColorBrewer", "gplots", "rtracklayer", "scales",
    "pheatmap", "ComplexHeatmap", "circlize", "paletteer"
  )
  for (pkg in req_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
  library(plot1cell)
})

# ------------------------------------------------------------
# 0. Project directory and input object
# ------------------------------------------------------------
project_dir <- "path/to/project_root"  # <-- modify to your own path
setwd(project_dir)

load("seu_symbol_SCT_clustered_annotated.RData")  # loads `seu`

# ------------------------------------------------------------
# 1. Subset main object to cluster C0 and recluster
# ------------------------------------------------------------

# Use original seurat_clusters as identities, then subset C0
seu@active.ident <- as.factor(seu$seurat_clusters)

glycan_seu <- subset(seu, idents = "0")

# SCTransform on cluster C0 (regress nCount_RNA)
glycan_seu <- SCTransform(
  glycan_seu,
  verbose        = FALSE,
  vars.to.regress = "nCount_RNA"
)

# PCA and elbow plot
glycan_seu <- RunPCA(glycan_seu, verbose = FALSE)
ElbowPlot(glycan_seu)

# Neighbors, clustering, and UMAP on top 4 PCs
glycan_seu <- FindNeighbors(glycan_seu, dims = 1:4)
glycan_seu <- FindClusters(glycan_seu, resolution = 0.15)
glycan_seu <- RunUMAP(glycan_seu, dims = 1:4)

# Basic palettes (kept consistent with original usage)
nature_palette2 <- c(
  "#3C5488", "#F39B7F", "#8491B4", "#91D1C2",
  "#F5A0A1", "#C2B5D8", "#FFB6C1", "#7E6148",
  "#E64B35", "#4DBBD5", "#00A087", "#B09C85"
)

# Same helper as in previous scripts (used in barplots)
expand_palette <- function(base_colors, n) {
  if (n <= length(base_colors)) base_colors[1:n] else colorRampPalette(base_colors)(n)
}

# A generic qualitative palette used in several plots
cell_type_cols <- c(
  "#1F77B4", "#2CA02C", "#FF7F0E", "#6A5ACD", "#8C564B",
  "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#F5A0A1",
  "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000",
  "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B",
  "#BDB76B", "#20B2AA", "#FF1493", "#FF4500", "#32CD32"
)

# Output directory for this figure panel (optional)
out_dir <- file.path(project_dir, "fig_Fig3_g_to_k")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------
# 2. UMAP plots – cluster, group, phenotype, location, species, family
#    (Figures 3g / 3h-style)
# ------------------------------------------------------------

# 2.1 UMAP by glycan subclusters
p1 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label     = FALSE,
  pt.size   = 0.5,
  group.by  = "seurat_clusters"
) +
  scale_color_manual(values = nature_palette2) +
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave(file.path(out_dir, "C0_glycan_UMAP_Clusters.png"), p1, width = 8, height = 7, dpi = 1200)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Clusters.pdf"),   p1, width = 8, height = 7)

# 2.2 UMAP by sample group
p2 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label     = FALSE,
  pt.size   = 1,
  group.by  = "group"
) +
  scale_color_manual(values = cell_type_cols) +
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave(file.path(out_dir, "C0_glycan_UMAP_Group.png"), p2, width = 8, height = 8, dpi = 1200)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Group.pdf"), p2, width = 8, height = 8)

# 2.3 UMAP by phenotype (DB vs WT)
p3 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label     = FALSE,
  pt.size   = 0.5,
  group.by  = "phenotype"
) +
  scale_color_manual(values = c("DB" = "#E64B35", "WT" = "#4DBBD5")) +
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave(file.path(out_dir, "C0_glycan_UMAP_Phenotype.png"), p3, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Phenotype.pdf"), p3, width = 8, height = 6)

# 2.4 UMAP by anatomical location (loc)
p4 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label     = FALSE,
  pt.size   = 0.5,
  group.by  = "loc"
) +
  scale_color_manual(values = c(
    "Cecum"  = "#1F77B4",
    "Colon"  = "#2CA02C",
    "Rectum" = "#FF7F0E"
  )) +
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave(file.path(out_dir, "C0_glycan_UMAP_Location.png"), p4, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Location.pdf"), p4, width = 8, height = 6)

# 2.5 UMAP by species
p5 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label     = FALSE,
  pt.size   = 0.5,
  group.by  = "species"
) +
  scale_color_manual(values = cell_type_cols) +
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave(file.path(out_dir, "C0_glycan_UMAP_Species.png"), p5, width = 9, height = 6, dpi = 300)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Species.pdf"), p5, width = 9, height = 6)

# 2.6 UMAP by family (Tax_Family)
p6 <- DimPlot(
  glycan_seu,
  reduction = "umap",
  label     = FALSE,
  pt.size   = 0.5,
  group.by  = "Tax_Family"
) +
  scale_color_manual(values = cell_type_cols) +
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 18)
  )

ggsave(file.path(out_dir, "C0_glycan_UMAP_Tax_Family.png"), p6, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Tax_Family.pdf"), p6, width = 8, height = 6)

# ------------------------------------------------------------
# 3. Functional genesets and module scores
#    (for dot plots – Figures 3i / 3j)
# ------------------------------------------------------------

genesets <- list(
  # A: Protective / beneficial

  # Polysaccharide Utilization Locus (PUL)
  Plant_PUL = c(
    "susC","susD","susE","susF",
    "GH13","GH32","GH28","GH10","GH11","CE1","CE4",
    "amyA","xynA","xynB","araF","araG","lacZ","malE","malF"
  ),
  # Butyrate core pathway
  Butyrate_core = c(
    "thl","hbd","crt","bcd","etfA","etfB","croR","bdhA","paaH","paaF","paaZ"
  ),
  # Alternative butyrate synthesis (ptb/buk)
  Butyrate_alt  = c("ptb","buk"),
  # Vitamin B12 and menaquinone synthesis
  Vitamin_B12_MK = c(
    "cobA","cobB","cobC","cobD","cobQ","cobS","cobT","cobU",
    "cbiA","cbiC","cbiD","cbiK","menA","menB","menC","ubiE"
  ),
  # Antioxidant defense
  Antioxidant = c(
    "katA","katE","sodA","sodB","sodC","ahpC","ahpF","trxA","trxB",
    "dps","perR","oxyR","msrA","msrB","tsaA","gor","grxA","tpx","bcp"
  ),

  # B: Risk / harmful

  # Mucin degradation
  Mucin_GH = c(
    "nanA","nanE","nanK","nanT","nagA","nagB","nagC","nagZ",
    "fucA","fucI","fucK","fucU","fucP",
    "GH33","GH29","GH95","GH20","sulfatase"
  ),
  # Secondary bile acid genes
  Bile_acid_secondary = c("baiA","baiCD","baiF"),
  Bile_Hydrolase       = c("cbh","bsh","bshA","bshB","bshC"),
  Efflux_Resistance    = c(
    "tolC","acrA","acrB","acrD","ompR","ompF","ompC",
    "emrA","emrB","norM","mdtK","mdtF","cusA","cusB","mexB"
  ),
  Bile_total = c(
    "cbh","bsh","bshA","bshB","bshC",
    "tolC","acrA","acrB","acrD","ompR","ompF","ompC",
    "emrA","emrB","norM","mdtK","mdtF","cusA","cusB","mexB"
  ),
  # LPS core synthesis
  LPS_Core = c(
    "lpxA","lpxB","lpxC","lpxD","lpxK","lpxL","lpxM","lpxP","kdsA","kdsB",
    "waaC","waaF","waaL","msbA","lptA","lptB","lptC"
  ),
  # LPS modification
  LPS_Mod_PEtN = c("eptA","pagP","lpxT","arnT","pmrA","pmrB","phoP","phoQ"),
  # Outer membrane lipid transport
  OM_Lipid_Transport = c(
    "mlaA","mlaB","mlaC","mlaD","mlaE","mlaF","lptD","lptE","vacJ"
  ),
  # Ethanolamine metabolism
  Ethanolamine = c(
    "eutA","eutB","eutC","eutD","eutE","eutG","eutH","eutJ","eutK","eutL",
    "eutM","eutN","eutP","eutQ","eutR","eutS","eutT","eutV","eutW","eutX"
  ),
  # Succinate production block
  Succinate_Block = c(
    "frdA","frdB","frdC","frdD","sdhA","sdhB","sdhC","sdhD",
    "mdh","fumA","fumB","fumC","sucA","sucB","sucC","sucD"
  ),
  # T6SS genes
  T6SS = c(
    "hcp","vgrG","tssA","tssB","tssC","tssD","tssE","tssF","tssG","tssH",
    "tssI","tssJ","tssK","tssL","tssM","tssN","clpV","dotU","icmF"
  ),
  # Flagella
  Flagella = c(
    "fliC","fliD","fliE","fliF","fliG","fliH","fliI","fliJ","fliK","fliL",
    "fliM","fliN","fliO","fliP","fliQ","fliR","fliS","fliT","flgA","flgB"
  ),
  # TMA/TMAO metabolism
  TMA_TMAO = c(
    "cutC","cutD","cntA","cntB","yeaW","yeaX",
    "caiA","caiB","caiC","caiD","caiE","caiF","caiT",
    "betT","torA","torZ","torC","dmsA","dmsB","dmsC"
  ),
  # Ethanol/acetaldehyde production
  Ethanol_Acetaldehyde = c(
    "pdc","adhE","adhA","adhB","adhP","adhC","aldA","aldB","aldH",
    "fucO","yqhD","yiaY","frmA","frmB","exaA","exaB","exaC","mdh"
  ),
  # Glycolipid load & lipid biosynthesis
  Sugar_Overload = c("ptsI","ptsH","lacC","lacE","lacF","lacG","galT"),
  Lipid_Biosynthesis = c("fabD","fabF"),

  # C: Contextual / double-edged

  # Propionate synthesis
  Propionate = c(
    "pccA","pccB","mutA","mutB","epi","prpE","prpC","prpD","prpB","acsA",
    "mmdA","lcdA","mmcA","prpF",
    "methylcitrate synthase","methylcitrate dehydratase"
  ),
  # Acetate synthesis
  Acetate = c(
    "pta","ackA","acs","poxB","yfiD","ldhA","pflB","pflA",
    "pykF","pykA","eno","pgk","gapA","tpiA","pgm"
  ),
  # Overall SCFA-related genes
  SCFA_total = c(
    "thl","hbd","crt","bcd","etfA","etfB","croR","bdhA","paaH","paaF","paaZ",
    "ptb","buk","atoA","atoD","fadE","etfQ","etfR",
    "pccA","pccB","mutA","mutB","epi","prpE","prpC","prpD","prpB","acsA",
    "mmdA","lcdA","mmcA","prpF",
    "pta","ackA","acs","poxB","yfiD","ldhA","pflB","pflA",
    "pykF","pykA","eno","pgk","gapA","tpiA","pgm"
  )
)

# Keep genes present in glycan_seu
all_genes <- rownames(glycan_seu)
genesets_clean <- lapply(genesets, function(g) intersect(g, all_genes))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]

# ------------------------------------------------------------
# 4. Remove old module score columns (if present) and recompute
# ------------------------------------------------------------

# Remove module-related columns by pattern
modules_to_remove <- grep(
  "Plant_PUL|Butyrate|Vitamin_B12|Antioxidant|Mucin_GH|Bile|LPS|Ethanolamine|Succinate|T6SS|Flagella|TMA_TMAO|Ethanol_Acetaldehyde|Propionate|Acetate|SCFA|Sugar_Overload|Lipid_Biosynthesis",
  colnames(glycan_seu@meta.data),
  value = TRUE
)

glycan_seu@meta.data <- glycan_seu@meta.data[
  , !(colnames(glycan_seu@meta.data) %in% modules_to_remove)
]

cols_to_remove <- c(
  "OM_Lipid_Transport9",
  "Efflux_Resistance8",
  "OM_Lipid_Transport12",
  "Plant_PUL1",
  "Butyrate_core2",
  "Butyrate_total3",
  "Vitamin_B12_MK4",
  "Antioxidant5",
  "Mucin_GH6",
  "Bile_Hydrolase7",
  "Bile_total9",
  "LPS_Core10",
  "LPS_Mod_PEtN11",
  "LPS_total13",
  "Ethanolamine14",
  "Succinate_Block15",
  "T6SS16",
  "Flagella17",
  "TMA_TMAO18",
  "Ethanol_Acetaldehyde19",
  "Sugar_Overload20",
  "Lipid_Biosynthesis21",
  "Propionate22",
  "Acetate23",
  "SCFA_total24",
  "Efflux_Resistance7"
)

glycan_seu@meta.data <- glycan_seu@meta.data[
  , !(colnames(glycan_seu@meta.data) %in% cols_to_remove)
]

# Recompute module scores with specified nbin and seed
glycan_seu <- AddModuleScore(
  glycan_seu,
  features = genesets_clean,
  name     = names(genesets_clean),
  nbin     = 24,
  seed     = 123
)

# ------------------------------------------------------------
# 5. Dot plots of module scores (Protective / Risk / SCFA)
#    -> Figures 3i / 3j-like
# ------------------------------------------------------------

# Protective module scores (selected)
p <- DotPlot(
  glycan_seu,
  features = c("Plant_PUL1", "Vitamin_B12_MK3", "Antioxidant4"),
  group.by = "seurat_clusters",
  scale    = TRUE
) +
  scale_color_gradient(low = "white", high = "red3") +
  coord_flip() +
  scale_size_continuous(range = c(4, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text   = element_text(face = "bold")
  ) +
  ggtitle("Protective")

ggsave(file.path(out_dir, "Dotplot_Protective.pdf"), p, width = 7,   height = 5.6)
ggsave(file.path(out_dir, "Dotplot_Protective.png"), p, width = 7,   height = 5.6, dpi = 600)

# Risk module scores
p <- DotPlot(
  glycan_seu,
  features = c(
    "LPS_Mod_PEtN11","LPS_Core10",
    "Ethanol_Acetaldehyde18","Sugar_Overload19",
    "Lipid_Biosynthesis20","Mucin_GH5"
  ),
  group.by = "seurat_clusters",
  scale    = TRUE
) +
  scale_color_gradient(low = "white", high = "#1F77B4") +
  coord_flip() +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text   = element_text(face = "bold")
  ) +
  ggtitle("Risk")

ggsave(file.path(out_dir, "Dotplot_Risk.pdf"), p, width = 7.5, height = 5.5)
ggsave(file.path(out_dir, "Dotplot_Risk.png"), p, width = 7.5, height = 5.5, dpi = 600)

# SCFA / contextual flux module scores
p <- DotPlot(
  glycan_seu,
  features = c("Acetate22","Propionate21","Butyrate_core2","SCFA_total23"),
  group.by = "seurat_clusters",
  scale    = TRUE
) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_flip() +
  scale_size_continuous(range = c(4, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text   = element_text(face = "bold")
  ) +
  ggtitle("SCFA_ContextFlux")

ggsave(file.path(out_dir, "Dotplot_SCFA_ContextFlux.pdf"), p, width = 6.5, height = 5.5)
ggsave(file.path(out_dir, "Dotplot_SCFA_ContextFlux.png"), p, width = 6.5, height = 5.5, dpi = 600)

# ------------------------------------------------------------
# 6. Dot plots of individual genes (protective / risk / SCFA)
#    (supporting detail for Figure 3g–k panel)
# ------------------------------------------------------------

# Protective genes
protective_genes <- unique(unlist(
  genesets[c("Plant_PUL","Butyrate_core","Butyrate_alt","Vitamin_B12_MK","Antioxidant")]
))
protective_genes_valid <- intersect(protective_genes, rownames(glycan_seu))

DotPlot(
  glycan_seu,
  features = protective_genes_valid,
  group.by = "seurat_clusters",
  scale    = TRUE
) +
  scale_color_gradient(low = "white", high = "red3") +
  coord_flip() +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text   = element_text(face = "bold")
  ) +
  ggtitle("Protective genes")

# Risk genes
risk_genes <- unique(unlist(
  genesets[c(
    "Mucin_GH","Bile_acid_secondary","Bile_Hydrolase",
    "Efflux_Resistance","LPS_Core","LPS_Mod_PEtN",
    "OM_Lipid_Transport","Ethanolamine","Succinate_Block",
    "T6SS","Flagella","TMA_TMAO",
    "Ethanol_Acetaldehyde","Sugar_Overload","Lipid_Biosynthesis"
  )]
))
risk_genes_valid <- intersect(risk_genes, rownames(glycan_seu))

DotPlot(
  glycan_seu,
  features = risk_genes_valid,
  group.by = "seurat_clusters",
  scale    = TRUE
) +
  scale_color_gradient(low = "white", high = "#1F77B4") +
  coord_flip() +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text   = element_text(face = "bold")
  ) +
  ggtitle("Risk genes")

# SCFA / contextual genes
scfa_genes <- unique(unlist(genesets[c("Propionate","Acetate","SCFA_total")]))
scfa_genes_valid <- intersect(scfa_genes, rownames(glycan_seu))

DotPlot(
  glycan_seu,
  features = scfa_genes_valid,
  group.by = "seurat_clusters",
  scale    = TRUE
) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_flip() +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text   = element_text(face = "bold")
  ) +
  ggtitle("SCFA genes")

# ------------------------------------------------------------
# 7. Cluster composition bar plots
#    (Figure 3k-style)
# ------------------------------------------------------------

# 7.1 Composition by group
glycan_seu@active.ident <- as.factor(glycan_seu$group)
prop_df <- as.data.frame(
  prop.table(table(glycan_seu$seurat_clusters, Idents(glycan_seu)))
)
colnames(prop_df) <- c("Group", "Cluster", "Freq")
prop_df$Group <- factor(prop_df$Group, levels = unique(prop_df$Group))

p_bar <- ggplot(prop_df, aes(x = Cluster, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = cell_type_cols) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title   = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(x = "Cluster", y = "Proportion", title = "Cluster Composition by Group")

p_bar
ggsave(file.path(out_dir, "C0_Bar_ClusterByGroup.png"), p_bar, width = 5, height = 4, dpi = 600)

# 7.2 Composition by phenotype (DB vs WT)
glycan_seu@active.ident <- as.factor(glycan_seu$phenotype)
prop_df <- as.data.frame(
  prop.table(table(Idents(glycan_seu), glycan_seu$seurat_clusters))
)
colnames(prop_df) <- c("Group", "Cluster", "Freq")
prop_df$Group <- factor(prop_df$Group, levels = unique(prop_df$Group))

p_bar <- ggplot(prop_df, aes(x = Cluster, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = expand_palette(nature_palette2, length(levels(prop_df$Group)))) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title   = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(x = "Cluster", y = "Proportion", title = "Cluster Composition by Phenotype")

p_bar
ggsave(file.path(out_dir, "C0_Bar_ClusterByPhenotype.png"), p_bar, width = 5, height = 4, dpi = 600)

# 7.3 Composition by location (loc)
glycan_seu@active.ident <- as.factor(glycan_seu$loc)
prop_df <- as.data.frame(
  prop.table(table(Idents(glycan_seu), glycan_seu$seurat_clusters))
)
colnames(prop_df) <- c("Group", "Cluster", "Freq")
prop_df$Group <- factor(prop_df$Group, levels = unique(prop_df$Group))

p_bar <- ggplot(prop_df, aes(x = Cluster, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "Cecum"  = "#1F77B4",
    "Colon"  = "#2CA02C",
    "Rectum" = "#FF7F0E"
  )) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title   = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(x = "Cluster", y = "Proportion", title = "Cluster Composition by Location")

p_bar
ggsave(file.path(out_dir, "C0_Bar_ClusterByLocation.png"), p_bar, width = 5, height = 4, dpi = 600)

# 7.4 Composition by group (alternative palette)
glycan_seu@active.ident <- as.factor(glycan_seu$group)
prop_df <- as.data.frame(
  prop.table(table(Idents(glycan_seu), glycan_seu$seurat_clusters))
)
colnames(prop_df) <- c("Group", "Cluster", "Freq")
prop_df$Group <- factor(prop_df$Group, levels = unique(prop_df$Group))

p_bar <- ggplot(prop_df, aes(x = Cluster, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "#1F77B4", "#17BECF", "#2CA02C", "#9ACD32", "#FF7F0E", "#F5A0A1"
  )) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title   = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(x = "Cluster", y = "Proportion", title = "Group Composition by Cluster")

p_bar
ggsave(file.path(out_dir, "C0_Bar_GroupByCluster.png"), p_bar, width = 5, height = 4, dpi = 600)

