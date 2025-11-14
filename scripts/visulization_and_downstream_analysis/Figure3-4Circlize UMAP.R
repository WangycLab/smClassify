# Circlize UMAP Visualization
# 1. Load packages
req_pkgs <- c(
  "Seurat", "tidyverse", "Matrix", "data.table", "ggrepel",
  "RColorBrewer", "gplots", "rtracklayer", "scales",
  "pheatmap", "ComplexHeatmap", "circlize", "paletteer", "plot1cell"
)
for (pkg in req_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}
library(plot1cell)
library(circlize)

# 2. Global parameters and colors
ROOT_DIR <- "D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009"
out_dir  <- ROOT_DIR
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

nature_palette <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85"
)
nature_palette2 <- c(
  "#3C5488", "#F39B7F", "#8491B4", "#91D1C2",
  "#F5A0A1", "#C2B5D8", "#FFB6C1", "#7E6148",
  "#E64B35", "#4DBBD5", "#00A087", "#B09C85"
)

expand_palette <- function(base_colors, n) {
  if (n <= length(base_colors)) base_colors[1:n] else colorRampPalette(base_colors)(n)
}

# Function for save plots
save_plot <- function(p, name, w=7, h=5, dpi=600) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width=w, height=h, dpi=dpi)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width=w, height=h)
}

# 3. Load Seurat object
set.seed(1234)
setwd(ROOT_DIR)

seu <- readRDS("stepF_seu_SCT_clustered.rds")

seu$group  <- sub("-[0-9]$", "", seu$orig.ident)
seu$position <- case_when(
  grepl("^Cecum",  seu$orig.ident) ~ "Cecum",
  grepl("^Colon",  seu$orig.ident) ~ "Colon",
  grepl("^Rectum", seu$orig.ident) ~ "Rectum",
  TRUE ~ NA_character_
)
seu$phenotype <- ifelse(grepl("WT", seu$orig.ident), "WT", "T2DM")
seu$region <- seu$position


anno_df <- read.csv("Cluster_functional_summary.csv", header=TRUE)
new.ids <- setNames(as.character(anno_df[[2]]), levels(seu))
seu@active.ident <- seu$seurat_clusters
seu <- RenameIdents(seu, new.ids)
seu$celltype <- seu@active.ident

save(seu, file="seu_symbol_SCT_clustered_annotated.RData")
load("seu_symbol_SCT_clustered_annotated.RData")

# 4. Prepare Circlize source data
seu@active.ident <- as.factor(seu$celltype)
circ_data <- prepare_circlize_data(seu, scale=0.7)

celltype_colors <- c("#4DBBD5","#E64B35","#00A087","#3C5488",
                     "#F39B7F","#91D1C2","#C2B5D8","#7E6148")

region_cols    <- expand_palette(nature_palette,  length(unique(seu$region)))
phenotype_cols <- expand_palette(brewer.pal(6, "Set2"), length(unique(seu$phenotype)))

# 5. Draw Circlize 
png(file.path(out_dir, "Circlize_Celltype.png"), width=8000, height=8000, res=1200)
plot_circlize(
  circ_data,
  do.label = F,
  pt.size = 0.4,
  col.use = celltype_colors,
  bg.color = "white",
  kde2d.n = 200,
  repel = TRUE,
  label.cex = 0.8
)
dev.off()

cairo_pdf(file.path(out_dir, "Circlize_Celltype.pdf"), width = 5, height = 5, onefile = FALSE, family = "sans")
op <- par(mar = rep(0, 4))
plot_circlize(
  circ_data,
  do.label = F,
  pt.size = 0.4,
  col.use = celltype_colors,
  bg.color = "white",
  kde2d.n = 200,
  repel = TRUE,
  label.cex = 1
)
add_track(circ_data, group = "region", colors = region_cols, track_num = 2)
add_track(circ_data, group = "phenotype", colors = phenotype_cols, track_num = 3)
par(op)
dev.off()


cairo_pdf(file.path(out_dir, "Circlize_Celltype.pdf"), width = 8, height = 8, onefile = FALSE, family = "sans")
op <- par(mar = rep(0, 4))
plot_circlize(
  circ_data,
  do.label = F,
  pt.size = 0.5,
  col.use = celltype_colors,
  bg.color = "white",
  kde2d.n = 200,
  repel = TRUE,
  label.cex = 1
)
par(op)
dev.off()



## Load packages for network
# install.packages(c("networkD3","dplyr","htmlwidgets","webshot2"))
library(networkD3)
library(dplyr)
library(htmlwidgets)
library(webshot2)     # Require Chrome/Chromium

## Load data
# Idents(sample) = Function£»sample$celltype = Cluster
seu@active.ident <- seu$celltype
data <- as.data.frame(prop.table(table(Idents(seu), seu$group)))
data <- as.data.frame(prop.table(table(Idents(seu), seu$group), margin = 1))
head(data)
df <- data.frame(Function = data$Var2,
                 Cluster  = data$Var1,
                 Freq     = data$Freq) 
head(df)
## Set nodes & links
nodes <- data.frame(name = c(unique(df$Function), unique(df$Cluster)))

links <- df |>
  mutate(IDsource = match(Function, nodes$name) - 1,
         IDtarget = match(Cluster,  nodes$name) - 1) |>
  dplyr::select(IDsource, IDtarget, Freq)


## Draw Sankey Plot
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value  = "Freq",   NodeID = "name",
                   nodeWidth = 20, fontSize = 30,
                   sinksRight = FALSE)
p


# Get metadata
meta_df <- seu@meta.data

# P1: Stacked bar plot group by phenotype
frac_df <- meta_df %>%
  group_by(celltype, phenotype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(freq = n / sum(n))

p1 <- ggplot(frac_df, aes(x = freq, y = fct_rev(factor(celltype)), fill = phenotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Fraction", y = NULL) +
  scale_fill_manual(values = c("T2DM" = "#E74C3C", "WT" = "#3498DB")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "right",
    panel.grid = element_blank()
  )

# P2: Stacked bar plot group by loc

frac_df <- meta_df %>%
  group_by(celltype, loc) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(freq = n / sum(n))

p1loc <- ggplot(frac_df, aes(x = freq, y = fct_rev(factor(celltype)), fill = loc)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Fraction", y = NULL) +
  scale_fill_manual(values = c("Cecum" = "#1f77b4",
                               "Colon" = "#2ca02c",
                               "Rectum" = "#ff7f0e")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "right",
    panel.grid = element_blank()
  )


# P2: Sum cell numbers
meta_df <- as.data.frame(seu@meta.data)

total_df <- meta_df %>%
  dplyr::count(celltype)

p2 <- ggplot(total_df, aes(x = n / 1000, y = fct_rev(factor(celltype)))) +
  geom_col(fill =  "#BCBD22", width = 0.7) +
  labs(x = "Cell Nums (x103)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    panel.grid = element_blank()
  )

# -------------------------------
library(patchwork)

# -------------------------------
p_combined <- p1 + p1loc + p2 + plot_layout(widths = c(1, 1, 0.8))
p_combined 
ggsave("seu_Cellsubtype_Phenotype_Stacked.pdf", p_combined, width = 8, height = 2.5, units = "in", dpi = 600)
ggsave("seu_Cellsubtype_Phenotype_Stacked.png", p_combined, width = 8, height = 2.5, units = "in", dpi = 600)
