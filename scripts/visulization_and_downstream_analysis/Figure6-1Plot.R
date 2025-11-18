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

load("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009/seu_symbol_SCT_clustered_annotated.RData")
setwd("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009/Figure6")
# Top 20 species
# Count and rank
sp_counts <- sort(table(seu$species), decreasing = TRUE)
n_top <- min(20, length(sp_counts))
top_species <- names(sp_counts)[seq_len(n_top)]
seu$top20_species <- ifelse(seu$species %in% top_species,
                            as.character(seu$species), "Others")
seu$top20_species <- factor(
  seu$top20_species,
  levels = c(top_species, "Others")
)
n_other <- sum(!(seu$species %in% top_species))
levels(seu$top20_species)[levels(seu$top20_species) == "Others"] 
# Quick check
table(seu$top20_species)

p<-FeaturePlot(seu, features=c("bepF"), cols=c("lightgrey","red"),pt.size = 0.1)+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("FeaturePlot_seu_bepF.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_bepF.png", p, width = 7.5, height = 6, dpi = 300)

p<-FeaturePlot(seu, features=c("bepG"), cols=c("lightgrey","red"),pt.size = 0.1)+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("FeaturePlot_seu_bepG.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_bepG.png", p, width = 7.5, height = 6, dpi = 300)


p<-FeaturePlot(seu, features=c("MGYG000413837-01743"), cols=c("lightgrey","red"),pt.size = 0.1)+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("FeaturePlot_seu_MGYG000413837-01743.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_MGYG000413837-01743.png", p, width = 7.5, height = 6, dpi = 300)

p<-FeaturePlot(seu, features=c("MGYG000413837-01568"), cols=c("lightgrey","red"),pt.size = 0.1)+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("FeaturePlot_seu_MGYG000413837-01568.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_MGYG000413837-01568.png", p, width = 7.5, height = 6, dpi = 300)



p<-FeaturePlot(seu, features=c("MGYG000413837-00889"), cols=c("lightgrey","red"),pt.size = 0.1)+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("FeaturePlot_seu_MGYG000413837-00889.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_MGYG000413837-00889.png", p, width = 7.5, height = 6, dpi = 300)



p<-FeaturePlot(seu, features=c("MGYG000413837-01569"), cols=c("lightgrey","red"),pt.size = 0.1)+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.title = element_text(size = 14)) 
p
ggsave("FeaturePlot_seu_MGYG000413837_01569.pdf", p, width = 7.5, height = 6)
ggsave("FeaturePlot_seu_MGYG000413837_01569.png", p, width = 7.5, height = 6, dpi = 300)


