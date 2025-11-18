library(networkD3)
library(dplyr)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(gplots)
library(rtracklayer)
library(scales)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
## Set work directory

load("seu_symbol_SCT_clustered_annotated.RData")

setwd("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009/Figure6")
nature_palette <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85"
)
cell_type_cols <- c("#1F77B4",  "#2CA02C","#FF7F0E","#6A5ACD", "#8C564B", "#E377C2", "#7F7F7F", 
                    "#BCBD22", "#17BECF", "#F5A0A1", "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32",
                    "#8B0000",  "#FFD700", "#DC143C", "#228B22", "#FF6347", "#483D8B", "#BDB76B", 
                    "#20B2AA", "#FF1493", "#FF4500", "#32CD32", "#3E8E41", "#20B2AA")

sample<-seu
sample@active.ident <- as.factor(sample$species)

Murigor <- subset(sample, idents = "Muribaculum gordoncarteri")


Murigor <- SCTransform(Murigor, verbose = FALSE,vars.to.regress = "nCount_RNA")
Murigor <- RunPCA(Murigor, verbose = FALSE)
ElbowPlot(Murigor)

Murigor <- FindNeighbors(Murigor, dims = 1:6)
Murigor <- FindClusters(Murigor, resolution = 0.1)
Murigor <- RunUMAP(Murigor, dims = 1:6)

DimPlot(Murigor, reduction = "umap", label = F, pt.size = 1,
        group.by = "seurat_clusters") +
  scale_color_manual(values = cell_type_cols)+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 

saveRDS(Murigor,"Murigor.Rds")
Murigor <-readRDS("Murigor.Rds")

p<-DimPlot(Murigor , reduction = "umap", label = F, pt.size = 0.1,
           group.by = "region") +
  theme_minimal() +
  scale_color_manual(values=c( "Cecum" = "#1f77b4",
                               "Colon" = "#2ca02c",
                               "Rectum" = "#ff7f0e"))+
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 14)) 
p
ggsave("Murigor _DimPlot_seurat_region.png", p, width =6, height = 4, dpi = 300)

p<-DimPlot(Murigor , reduction = "umap", label = F, pt.size = 0.1,
           group.by = "phenotype") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("Murigor _DimPlot_seurat_phenotype.png", p, width = 6, height = 4, dpi = 300)

nature_palette <- c(
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85","#E64B35", "#4DBBD5", "#00A087", "#3C5488"
)
p<-DimPlot(Murigor , reduction = "umap", label = F, pt.size = 0.1,
           group.by = "seurat_clusters") +
  theme_minimal() +
  scale_color_manual(values=nature_palette)+
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("Murigor _DimPlot_seurat_clusters.png", p, width = 6, height = 4, dpi = 300)



## Load Murigor object (use load("Murigor.RData") if previously saved)
Idents(Murigor) <- Murigor$seurat_clusters       # Ensure identity

message("Finding markers")
markers <- FindAllMarkers(
  object          = Murigor,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.15
) 
markers

topN <- markers %>% group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 50) %>% ungroup()
table(topN $cluster)
write.csv(topN,"Murigor_markers_top50.csv")
topN
top5_nonMGYG <- topN %>%
  dplyr::filter(!grepl("^MGYG", gene)) %>%      # Remove genes starting with MGYG
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%       # Take top 10 per cluster
  ungroup()

genes <- unique(top5_nonMGYG$gene)
genes


# Plot
p <- DotPlot(Murigor, features = genes) +
  scale_color_gradient(low = "#2ECC71", high = "#8E44AD") +
  scale_size_continuous(range = c(1, 8)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title   = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.grid   = element_blank()
  ) +
  labs(x = "", y = "Functional Subtypes", title = "DotPlot of Top Marker Genes")

# Show plot
print(p)

# Save plot
ggsave("Murigor_dotplot_topgenes.pdf", p, width = 14, height = 4)
ggsave("Murigor_dotplot_topgenes.png", p, width = 14, height = 4, dpi = 300)


a <- read.csv('Murib_type.csv')
a
levels(Murigor)
new.cluster.ids <- as.character(a[,2])
new.cluster.ids
levels(Murigor)
Murigor@active.ident <- as.factor(Murigor$seurat_clusters)
names(new.cluster.ids) <- levels(Murigor)
Murigor <- RenameIdents(Murigor, new.cluster.ids)
Murigor$cellsubtype <- Murigor@active.ident


cell_type_cols <- c(
 "#C2B5D8","#1F77B4",  "#BCBD22", "#17BECF",
  "#C2B5D8", "#FFB6C1", "#3CB371", "#9ACD32", "#8B0000", "#FFD700",
  "#DC143C", "#228B22", "#483D8B", "#BDB76B", "#20B2AA",
  "#FF1493", "#FF4500", "#32CD32", "#3E8E41")
p_Subtype <- DimPlot(Murigor, reduction = "umap", label = FALSE, pt.size = 0.1, 
                     group.by = "cellsubtype") +
  scale_color_manual(values = cell_type_cols) +
  ggtitle("Cell Subtype") + common_theme
p_Subtype
ggsave("Murigor_UMAP_Subtype.png", p_Subtype, width = 7, height = 3.5, dpi = dpi_out)
ggsave("Murigor_UMAP_Subtype.pdf", p_Subtype, width = 7, height = 3.5, device = cairo_pdf)

pal_group <- c(
  "Cecum-WT"   = "#6baed6",  
  "Cecum-DB"   = "#08519c",  
  "Colon-WT"   = "#74c476", 
  "Colon-DB"   = "#006d2c",  
  "Rectum-WT"  = "#fdae6b",  
  "Rectum-DB"  = "#e6550d"   
)

meta <- Murigor@meta.data

# 1. Compute fractions
prop_df <- meta %>%
  group_by(seurat_clusters, group) %>%
  summarise(n = n()) %>%
  group_by(group) %>%
  mutate(Fraction = n / sum(n))

# 2. Draw stacked bar chart
p<-ggplot(prop_df, aes(x = seurat_clusters, y = Fraction, fill = group)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() +   # Horizontal bar chart
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Fraction", fill = "Cell type") +
  theme_classic(base_size = 18) +
  scale_fill_manual(values = pal_group) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )
p
ggsave("Murigor_prop_df_seurat_clusters_group.pdf", p, w=6, h=4)


#==========monocle3=====================
library(Seurat)
library(monocle3) 
library(ggplot2)


Murigor@active.ident <- as.factor(Murigor$seurat_clusters)

seurat_obj<-Murigor
          
# Step 2: Extract count matrix
expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

# Step 3: Build Monocle3 metadata and feature data
cell_metadata <- seurat_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expr_matrix))
rownames(gene_annotation) <- rownames(expr_matrix)

# Step 4: Create cell_data_set object
cds <- new_cell_data_set(expr_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# Step 5: Pass Seurat UMAP embeddings to Monocle3
reducedDims(cds)$UMAP <- seurat_obj@reductions$umap@cell.embeddings

# Step 6: Add cluster info for visualization
clusters <- seurat_obj$seurat_clusters
colData(cds)$cluster <- as.character(clusters)

# Step 7: Learn trajectory and compute pseudotime
cds <- cluster_cells(cds, reduction_method = "UMAP")  # Cluster based on UMAP
cds <- learn_graph(cds, use_partition = TRUE)
cds <- order_cells(cds)  # Optionally choose root cells manually


# Step 8: Plot trajectory
pseudotime<-plot_cells(cds,
                       color_cells_by = "pseudotime",
                       show_trajectory_graph = TRUE,
                       label_groups_by_cluster = TRUE,cell_size = 0.5,alpha = 0.5,
                       label_leaves = TRUE,
                       label_branch_points = TRUE)
pseudotime
ggsave("Murigor_pseudotime.png", pseudotime, width =6, height = 4, dpi = 300)
ggsave("Murigor_pseudotime.pdf", pseudotime, width =6, height = 4)


# Find genes significantly changing along pseudotime
deg_pseudotime <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# Sort by pseudotime significance
deg_pseudotime <- deg_pseudotime[order(deg_pseudotime$q_value), ]
head(deg_pseudotime)


target_genes<- c(# Metabolism-related
  "gdhA", "glnA", "gadB", "gadC",
  # Fatty acid/Urea response
  "fabF",  "ureG",
  # ROS/NO response
  "sodB",  "rbr")


rowData(cds)$gene_short_name <- rownames(cds)

ptop_genes <- plot_genes_in_pseudotime(
  cds[target_genes, ],
  ncol = 4,
  min_expr = 0.5,
  color_cells_by = "pseudotime",cell_size = 1
) +
  theme(
    axis.text = element_text(size = 14),       # Text size
    axis.title = element_text(size = 16, face = "bold"), # Bold titles
    strip.text = element_text(size = 16, face = "bold"), # Facet labels
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

ptop_genes
ggsave("Murigor_ptop_T2DM_core_genes1.png", ptop_genes, width =14, height = 5, dpi = 300)
ggsave("Murigor_ptop_T2DM_core_genes1.pdf", ptop_genes, width =14, height = 5)


library(Seurat)
library(ggplot2)
library(dplyr)


pseudotime_vec <- pseudotime(cds)
head(pseudotime_vec)
# Build dataframe
pseudotime_df <- data.frame(cell = names(pseudotime_vec), pseudotime = pseudotime_vec)

# Add metadata
Murigor <- AddMetaData(Murigor, metadata = pseudotime_df)

meta_df <- Murigor@meta.data %>%
  dplyr::filter(!is.na(pseudotime))  # Remove missing values

# Plot density curves (ridgeline-like)
p_phenotype <- ggplot(meta_df, aes(x = pseudotime, fill = phenotype)) +
  geom_density(alpha = 0.6, color = NA) +
  facet_wrap(~phenotype, ncol = 1) +
  scale_fill_manual(values = c("T2DM" = "#E74C3C", "WT" = "#3498DB")) +
  labs(x = "Pseudotime", y = "Cell density") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.title = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )
p_phenotype
ggsave("Murigor_Nature_Phenotype_Pseudotime.pdf", p_phenotype, width = 4, height = 3.5, units = "in", dpi = 600)
ggsave("Murigor_Nature_Phenotype_Pseudotime.png", p_phenotype, width = 4, height = 3.5, units = "in", dpi = 600)


p_region <- ggplot(meta_df, aes(x = pseudotime, fill = region)) +
  geom_density(alpha = 0.6, color = NA) +
  facet_wrap(~region, ncol = 1) +
  scale_fill_manual(values = c("Cecum" = "#1F77B4", "Colon" = "#2CA02C", "Rectum" = "#FF7F0E")) +
  labs(x = "Pseudotime", y = "Cell density") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.title = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )

ggsave("Murigor_Nature_Region_Pseudotime.pdf", p_region, width = 4, height = 3.5, units = "in", dpi = 600)
ggsave("Murigor_Nature_Region_Pseudotime.png", p_region, width = 4, height = 3.5, units = "in", dpi = 600)



#=========== Heatmap ============================================
# Select significant genes (customizable)
sig_genes <- rownames(subset(deg_pseudotime, q_value < 1e-4))[1:100]
sig_genes
sig_genes <- intersect(sig_genes, rownames(seurat_obj[["SCT"]]))


sig_genes

log_expr <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")[sig_genes, ]
log_expr <- as.matrix(log_expr[, cell_order])
log_expr_scaled <- t(scale(t(log_expr)))


# Step 1: Add gene_id (use rownames as gene_id)
deg_pseudotime$gene_id <- rownames(deg_pseudotime)

# Step 3: Filter out MGYG and NA names; take top 100 by q_value
sig_genes <- deg_pseudotime %>%
  filter(!is.na(gene_id), !grepl("^MGYG", gene_id)) %>%
  arrange(q_value) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  slice_head(n = 100)
sig_genes$gene_id


sig_genes <- intersect(sig_genes$gene_id, rownames(seurat_obj[["SCT"]]))
sig_genes
# 1. Preprocess expression matrix
log_expr <- GetAssayData(seurat_obj, assay = "SCT", slot = "count")[sig_genes, ]
log_expr <- as.matrix(log_expr[, cell_order])
log_expr_scaled <- t(scale(t(log_expr)))


# 2. Hierarchical clustering to get modules (k = 4; customizable)
row_dist <- dist(log_expr_scaled)
row_clust <- hclust(row_dist, method = "complete")
gene_modules <- cutree(row_clust, k = 4)  # Modify k to change number of modules

# 3. Build annotation dataframe
row_anno <- data.frame(
  Module = factor(paste0("Module", gene_modules)),
  row.names = rownames(log_expr_scaled)
)
annotation_colors = list(
  Module = c(
    Module1 = "#66C2A5",
    Module2 = "#FC8D62",
    Module3 = "#8DA0CB",
    Module4 = "#E78AC3"
  )
)
# Customize module colors using Set2 palette
module_levels <- levels(row_anno$Module)
module_colors <- RColorBrewer::brewer.pal(length(module_levels), "Set2")
names(module_colors) <- module_levels

# Combine into list
anno_colors <- list(Module = module_colors)

# 4. Set heatmap colors
my_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
my_breaks <- seq(-2, 2, length.out = 101)


## -------- 1. Use k = 4 to define gene_modules --------
row_dist <- stats::dist(log_expr_scaled)
row_clust <- hclust(row_dist, method = "complete")

gene_modules <- cutree(row_clust, k = 4)

## -------- 2. Build gene_info and compute sd --------
gene_info <- data.frame(
  Name   = rownames(log_expr_scaled),
  Module = paste0("Module", gene_modules),
  stringsAsFactors = FALSE
)
gene_info$sd <- apply(log_expr_scaled, 1, sd)

## -------- 3. For each module, take top 5 genes by sd --------
top_genes_by_module <- gene_info %>%
  group_by(Module) %>%
  slice_max(order_by = sd, n =5, with_ties = FALSE) %>%
  ungroup()

## -------- 4. Prepare annotation_row: per-row Module --------
row_anno_all <- data.frame(
  Module = gene_info$Module,
  row.names = gene_info$Name
)

## -------- 5. Prepare labels_row: show names only for top genes --------
labels_vec <- ifelse(
  gene_info$Name %in% top_genes_by_module$Name,
  gene_info$Name,
  ""                     # No label
)

## -------- 6. Module colors ----------
module_levels  <- sort(unique(gene_info$Module))
module_colors  <- RColorBrewer::brewer.pal(length(module_levels), "Set2")
names(module_colors) <- module_levels
anno_colors <- list(Module = module_colors)


pheatmap(log_expr_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,           # Must be TRUE to show labels_row
         labels_row   = labels_vec,      # Show labels only for top genes
         show_colnames = FALSE,
         main  = "Pseudotime Heatmap (Top-labelled)",
         color = my_colors,
         breaks = my_breaks,
         annotation_row    = row_anno_all,
         annotation_colors = anno_colors)


png(filename = "Murigor_monocle3_pheatmap.png", width = 4, height =8, units = "in", res = 600)

# 5. Heatmap with module colors
pheatmap(log_expr_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Pseudotime Heatmap with Modules",
         color = my_colors,
         breaks = my_breaks,
         annotation_row = row_anno,
         annotation_colors = anno_colors)
dev.off()



png(filename = "Murigor_monocle3_pheatmap.pdf", width = 4, height =8)

# 5. Heatmap with module colors
pheatmap(log_expr_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Pseudotime Heatmap with Modules",
         color = my_colors,
         breaks = my_breaks,
         annotation_row = row_anno,
         annotation_colors = anno_colors)
dev.off()



urease_genes <- c( 
  # ---- Accessory/assembly ----
 "ureG", "ureF","ureE","ureD",
  
  # ---- Core subunits ----
  "ureC","ureB","ureA")

p<-DotPlot (seu,
         features = urease_genes  , # alpha-galactosidase,
         group.by = "top20_species",
         scale = TRUE) +
  scale_color_gradient(
    low =  "white", high = "red3") +
  coord_flip() +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Nitrogen metabolism")
p
ggsave("DotPlot_seu_Nitrogen metabolism_species.png", p, width = 10, height = 5, dpi = 300)
ggsave("DotPlot_seu_Nitrogen metabolism_species.pdf", p, width = 10, height = 5)


p<-DotPlot(seu,
         features = urease_genes , # alpha-galactosidase,
         group.by = "group",
         scale = TRUE) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3",
    midpoint = 0   # Midpoint for diverging color scale
  ) +
  coord_flip() +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Nitrogen metabolism")
p
ggsave("DotPlot_seu_Nitrogen metabolism_group.png", p, width = 5, height = 5, dpi = 300)
ggsave("DotPlot_seu_Nitrogen metabolism_group.pdf", p, width = 5, height = 5)



p<-DotPlot (Murigor,
            features = c( 
              # ---- Accessory/assembly ----
              "ureG", "ureF","ureE","ureD",
              # ---- Core subunits ----
              "ureC","ureB") , # alpha-galactosidase,
            group.by = "cellsubtype",
            scale = TRUE) +
  scale_color_gradient(
    low =  "white", high = "red3") +
  coord_flip() +
  scale_size_continuous(range = c(0, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold"),
      plot.margin  = margin(t = 25, r = 25, b = 25, l = 35, unit = "mm")
  ) +
  ggtitle("Nitrogen metabolism")
p
ggsave("DotPlot_Murigor_Nitrogen metabolism_cellsubtype.png", p, width = 7, height = 8, dpi = 300)
ggsave("DotPlot_Murigor_Nitrogen metabolism_cellsubtype.pdf", p, w=6, h=3)


p<-DotPlot(Murigor,
           features =  c( 
             "ureG", "ureF","ureE","ureD",
             "ureC","ureB","ureA"),
           group.by = "group",
           scale = TRUE) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3",
    midpoint = 0   # Midpoint for diverging color scale
  ) +
  coord_flip() +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Nitrogen metabolism")
p
ggsave("DotPlot_Murigor_Nitrogen metabolism_group.png", p, width = 6, height = 5, dpi = 300)
ggsave("DotPlot_Murigor_Nitrogen metabolism_group.pdf", p, w=5, h=5)



p<-DotPlot (Murigor,
            features =  c("susC",   
            "susD",   
            "pulA"),   # Pullulanase; degrades pullulan, starch/maltose  , # alpha-galactosidase,
            group.by = "cellsubtype",
            scale = TRUE) +
  scale_color_gradient(
    low = "white", high = "red3"
  ) +
  coord_flip() +
  scale_size_continuous(range = c(0, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold"),
    plot.margin  = margin(t = 25, r = 25, b = 25, l = 35, unit = "mm")
  ) +
  ggtitle("polysaccharide degradation")
p
ggsave("DotPlot_seu_polysaccharide degradation_cellsubtype.png", p, width = 7, height = 8, dpi = 300)
ggsave("DotPlot_seu_polysaccharide degradation_cellsubtype.pdf", p, width = 7, height = 8)


p<-DotPlot(Murigor,
           features =  c( # ---- Accessory/assembly ----
                          "ureG", "ureF","ureE","ureD",
                          # ---- Core subunits ----
                          "ureC","ureB",
                         "bepF",
                         "bepG",
                         "MGYG000413837-01569", 
                         "MGYG000413837-01743",
                         "pulA",
                         "susD", 
                         "susC") , # alpha-galactosidase,
           group.by = "group",
           scale = TRUE) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red3",
    midpoint = 0   # Midpoint for diverging color scale
  ) +
  coord_flip() +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Polysaccharide degradation")
p
ggsave("DotPlot_Murigor_Polysaccharide degradation-Urease_group.png", p, width = 5, height = 5, dpi = 300)
ggsave("DotPlot_Murigor_Polysaccharide degradation-Urease_group.pdf", p, w=6.5, h=5)

