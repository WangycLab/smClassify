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

setwd("D:/BaiduSyncdisk/post-doc/T2DMmicrobiome/mouse_colon/20251009")
load("seu_symbol_SCT_clustered_annotated.RData")


seu@active.ident <- as.factor(seu$seurat_clusters)
#Plot subcluster C0
glycan_seu <- subset(seu, idents = "0")
glycan_seu <- SCTransform(glycan_seu, verbose = FALSE,vars.to.regress = "nCount_RNA")
glycan_seu <- RunPCA(glycan_seu, verbose = FALSE)
ElbowPlot(glycan_seu)

glycan_seu <- FindNeighbors(glycan_seu, dims = 1:4)
glycan_seu <- FindClusters(glycan_seu, resolution = 0.15)
glycan_seu <- RunUMAP(glycan_seu, dims = 1:4)

nature_palette2 <- c(
  "#3C5488", "#F39B7F", "#8491B4", "#91D1C2",
  "#F5A0A1", "#C2B5D8", "#FFB6C1", "#7E6148",
  "#E64B35", "#4DBBD5", "#00A087", "#B09C85"
)
## Cluster
p1 <- DimPlot(glycan_seu, reduction = "umap", label = F, pt.size = 0.5,
              group.by = "seurat_clusters") +
  scale_color_manual(values = nature_palette2) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text  = element_text(size = 18),
        plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 18))
p1
ggsave( "C0_glycan_UMAP_Clusters.png", p1, width=8, height=7, dpi=1200)
ggsave("C0_glycan_UMAP_Clusters.pdf", p1, width=8, height=7, device="pdf")


## Group
p2 <- DimPlot(glycan_seu, reduction = "umap", label = F, pt.size = 1,
              group.by = "group") +
  scale_color_manual(values = cell_type_cols) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text  = element_text(size = 18),
        plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 18))

ggsave(file.path(out_dir, "C0_glycan_UMAP_Group.png"), p2, width=8, height=8, dpi=1200)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Group.pdf"), p2, width=8, height=8, device="pdf")


## Phenotype
p3 <- DimPlot(glycan_seu, reduction = "umap", label = F, pt.size = 0.5,
              group.by = "phenotype") +
  scale_color_manual(values = c("#E64B35","#4DBBD5")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text  = element_text(size = 18),
        plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 18))

ggsave(file.path(out_dir, "C0_glycan_UMAP_Phenotype.png"), p3, width=8, height=6, dpi=300)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Phenotype.pdf"), p3, width=8, height=6, device="pdf")


## Location
p4 <- DimPlot(glycan_seu, reduction = "umap", label = F, pt.size = 0.5,
              group.by = "loc") +
  scale_color_manual(values = c("#1F77B4", "#2CA02C", "#FF7F0E")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text  = element_text(size = 18),
        plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 18))

ggsave(file.path(out_dir, "C0_glycan_UMAP_Location.png"), p4, width=8, height=6, dpi=300)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Location.pdf"), p4, width=8, height=6, device="pdf")


## Top 20 species
p5 <- DimPlot(glycan_seu, reduction = "umap", label = F, pt.size = 0.5,
              group.by = "species") +
  scale_color_manual(values = cell_type_cols) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text  = element_text(size = 18),
        plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 18))

ggsave(file.path(out_dir, "C0_glycan_UMAP_Species.png"), p5, width=9, height=6, dpi=300)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Species.pdf"), p5, width=9, height=6, device="pdf")

glycan_seu$Tax_Family
## Top 20 species
p6 <- DimPlot(glycan_seu, reduction = "umap", label = F, pt.size = 0.5,
              group.by = "Tax_Family") +
  scale_color_manual(values = cell_type_cols) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text  = element_text(size = 18),
        plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 18))

ggsave(file.path(out_dir, "C0_glycan_UMAP_Tax_Family.png"), p6, width=8, height=6, dpi=300)
ggsave(file.path(out_dir, "C0_glycan_UMAP_Tax_Family.pdf"), p6, width=8, height=6, device="pdf")



genesets <- list(
  
  # A:Protective/bennificial
  # Polysaccharide Utilization Locus (PUL)
  Plant_PUL = c("susC","susD","susE","susF",
                "GH13","GH32","GH28","GH10","GH11","CE1","CE4",
                "amyA","xynA","xynB","araF","araG","lacZ","malE","malF"), 
  # Butyrate core pathway
  Butyrate_core = c("thl","hbd","crt","bcd","etfA","etfB","croR","bdhA","paaH","paaF","paaZ"), 
  # Alternative butyrate synthesis pathway: production of butyrate via phosphotransbutyrylase (ptb) and butyrate kinase (buk)
  Butyrate_alt  = c("ptb","buk"),
  # Vitamin B12 and menaquinone synthesis
  Vitamin_B12_MK = c("cobA","cobB","cobC","cobD","cobQ","cobS","cobT","cobU",
                     "cbiA","cbiC","cbiD","cbiK","menA","menB","menC","ubiE"), 
  # Antioxidant defense
  Antioxidant = c("katA","katE","sodA","sodB","sodC","ahpC","ahpF","trxA","trxB",
                  "dps","perR","oxyR","msrA","msrB","tsaA","gor","grxA","tpx","bcp"), 

  
  # B.Risk/Harmful
  # Mucin degradation
  Mucin_GH = c("nanA","nanE","nanK","nanT","nagA","nagB","nagC","nagZ",
               "fucA","fucI","fucK","fucU","fucP","GH33","GH29","GH95","GH20","sulfatase"), 
  # Bile acid related
  Bile_acid_secondary = c("baiA","baiCD","baiF"),
  Bile_Hydrolase = c("cbh","bsh","bshA","bshB","bshC"), 
  Efflux_Resistance = c("tolC","acrA","acrB","acrD","ompR","ompF","ompC",
                        "emrA","emrB","norM","mdtK","mdtF","cusA","cusB","mexB"), 
  # Overview of bile acid tolerance
  Bile_total = c("cbh","bsh","bshA","bshB","bshC",
                 "tolC","acrA","acrB","acrD","ompR","ompF","ompC",
                 "emrA","emrB","norM","mdtK","mdtF","cusA","cusB","mexB"), 
  # LPS core synthesis
  LPS_Core = c("lpxA","lpxB","lpxC","lpxD","lpxK","lpxL","lpxM","lpxP","kdsA","kdsB",
               "waaC","waaF","waaL","msbA","lptA","lptB","lptC"), 
  # LPS modification
  LPS_Mod_PEtN = c("eptA","pagP","lpxT","arnT","pmrA","pmrB","phoP","phoQ"), 
  # Outer membrane lipid transport
  OM_Lipid_Transport = c("mlaA","mlaB","mlaC","mlaD","mlaE","mlaF","lptD","lptE","vacJ"), 
  # Ethanolamine metabolism
  Ethanolamine = c("eutA","eutB","eutC","eutD","eutE","eutG","eutH","eutJ","eutK","eutL",
                   "eutM","eutN","eutP","eutQ","eutR","eutS","eutT","eutV","eutW","eutX"), 
  # Succinate production
  Succinate_Block = c("frdA","frdB","frdC","frdD","sdhA","sdhB","sdhC","sdhD",
                      "mdh","fumA","fumB","fumC","sucA","sucB","sucC","sucD"), 
  # Type VI secretion system
  T6SS = c("hcp","vgrG","tssA","tssB","tssC","tssD","tssE","tssF","tssG","tssH",
           "tssI","tssJ","tssK","tssL","tssM","tssN","clpV","dotU","icmF"), 
  # Flagellar system
  Flagella = c("fliC","fliD","fliE","fliF","fliG","fliH","fliI","fliJ","fliK","fliL",
               "fliM","fliN","fliO","fliP","fliQ","fliR","fliS","fliT","flgA","flgB"), 
  # TMA/TMAO metabolism
  TMA_TMAO = c("cutC","cutD","cntA","cntB","yeaW","yeaX",
               "caiA","caiB","caiC","caiD","caiE","caiF","caiT",
               "betT","torA","torZ","torC","dmsA","dmsB","dmsC"), 
  #  Ethanol/acetaldehyde production
  Ethanol_Acetaldehyde = c("pdc","adhE","adhA","adhB","adhP","adhC","aldA","aldB","aldH",
                           "fucO","yqhD","yiaY","frmA","frmB","exaA","exaB","exaC","mdh"), 
  # Glycolipid metabolic load & membrane lipid enhancement
  Sugar_Overload = c("ptsI","ptsH","lacC","lacE","lacF","lacG","galT"), 
  Lipid_Biosynthesis = c("fabD","fabF"), 

  
  # C. Contextual/Double-edged sword
  # Propionate synthesis
  Propionate = c("pccA","pccB","mutA","mutB","epi","prpE","prpC","prpD","prpB","acsA",
                 "mmdA","lcdA","mmcA","prpF","methylcitrate synthase","methylcitrate dehydratase"), 
  # Acetate synthesis
  Acetate = c("pta","ackA","acs","poxB","yfiD","ldhA","pflB","pflA",
              "pykF","pykA","eno","pgk","gapA","tpiA","pgm"), 
  # Overview of short-chain fatty acids
  SCFA_total = c("thl","hbd","crt","bcd","etfA","etfB","croR","bdhA","paaH","paaF","paaZ",
                 "ptb","buk","atoA","atoD","fadE","etfQ","etfR",
                 "pccA","pccB","mutA","mutB","epi","prpE","prpC","prpD","prpB","acsA",
                 "mmdA","lcdA","mmcA","prpF",
                 "pta","ackA","acs","poxB","yfiD","ldhA","pflB","pflA",
                 "pykF","pykA","eno","pgk","gapA","tpiA","pgm") ,
  Bile_acid_secondary = c("baiA","baiCD","baiF")
  
)



# Keep genes that exist in the Seurat object
all_genes <- rownames(glycan_seu)

genesets_clean <- lapply(genesets, function(g) intersect(g, all_genes))
genesets_clean <- genesets_clean[sapply(genesets_clean, length) > 0]
genesets_clean
# Remove columns whose names contain these modules
modules_to_remove <- grep("Plant_PUL|Butyrate|Vitamin_B12|Antioxidant|Mucin_GH|Bile|LPS|Ethanolamine|Succinate|T6SS|Flagella|TMA_TMAO|Ethanol_Acetaldehyde|Propionate|Acetate|SCFA|Sugar_Overload|Lipid_Biosynthesis",
                          colnames(glycan_seu@meta.data), value = TRUE)

glycan_seu@meta.data <- glycan_seu@meta.data[ , !(colnames(glycan_seu@meta.data) %in% modules_to_remove) ]
cols_to_remove <- c("OM_Lipid_Transport9",
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
                    "SCFA_total24","Efflux_Resistance7")

glycan_seu@meta.data <- glycan_seu@meta.data[ , !(colnames(glycan_seu@meta.data) %in% cols_to_remove) ]

# Calculate module score
glycan_seu <- AddModuleScore(
  glycan_seu,
  features = genesets_clean,
  name = names(genesets_clean),
  nbin = 24,    
  seed = 123
)
colnames(glycan_seu@meta.data)

p<-DotPlot(glycan_seu, features = c("Plant_PUL1" ,  "Vitamin_B12_MK3"  , "Antioxidant4"  ),
        group.by = "seurat_clusters",scale = T
) +
  scale_color_gradient(
    low =  "white", high = "red3" 
  ) +
  coord_flip() +
  scale_size_continuous(range = c(4, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Protective")
ggsave("Dotplot_Protective.pdf", p, width=7, height=5.6, device="pdf")


p<-DotPlot(glycan_seu, features = c(    "LPS_Mod_PEtN11" ,  "LPS_Core10" ,
                                     "Ethanol_Acetaldehyde18" ,        
                                     "Sugar_Overload19"  ,     
                                     "Lipid_Biosynthesis20"  ,  "Mucin_GH5"    
),
group.by = "seurat_clusters",scale = T
) +
  scale_color_gradient(
    low =  "white",  high = "#1F77B4"
  ) +
  coord_flip() +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Risk")
ggsave("Dotplot_Risk.pdf", p, width=7.5, height=5.5, device="pdf")

p<-DotPlot(glycan_seu, features = c("Acetate22"  ,"Propionate21" , "Butyrate_core2" ,            
                                 "SCFA_total23"      ),
        group.by = "seurat_clusters",scale = T
) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_flip() +
  scale_size_continuous(range = c(4, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("SCFA_ContextFlux")
ggsave("Dotplot_SCFA_ContextFlux.pdf", p, width=6.5, height=5.5, device="pdf")



# Protective
protective_genes <- unique(unlist(genesets[c("Plant_PUL","Butyrate_core","Butyrate_alt",
                                             "Vitamin_B12_MK","Antioxidant")]))
# Keep genes that exist in the Seurat object
protective_genes_valid <- intersect(protective_genes, rownames(glycan_seu))

DotPlot(glycan_seu, features = protective_genes_valid,
        group.by = "seurat_clusters", scale = TRUE) +
  scale_color_gradient(low = "white", high = "red3") +
  coord_flip() +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Protective genes")


#Risk
risk_genes <- unique(unlist(genesets[c("Mucin_GH","Bile_acid_secondary","Bile_Hydrolase",
                                       "Efflux_Resistance","LPS_Core","LPS_Mod_PEtN",
                                       "OM_Lipid_Transport","Ethanolamine","Succinate_Block",
                                       "T6SS","Flagella","TMA_TMAO",
                                       "Ethanol_Acetaldehyde","Sugar_Overload","Lipid_Biosynthesis")]))

risk_genes_valid <- intersect(risk_genes, rownames(glycan_seu))

DotPlot(glycan_seu, features = risk_genes_valid,
        group.by = "seurat_clusters", scale = TRUE) +
  scale_color_gradient(low = "white", high = "#1F77B4") +
  coord_flip() +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("Risk genes")


#SCFA/Contextual
scfa_genes <- unique(unlist(genesets[c("Propionate","Acetate","SCFA_total")]))
scfa_genes_valid <- intersect(scfa_genes, rownames(glycan_seu))

DotPlot(glycan_seu, features = scfa_genes_valid,
        group.by = "seurat_clusters", scale = TRUE) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_flip() +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(face = "bold")
  ) +
  ggtitle("SCFA genes")


#Plot barplot
glycan_seu@active.ident <- as.factor(glycan_seu$group)
prop_df <- as.data.frame(prop.table(table( glycan_seu$seurat_clusters,Idents(glycan_seu))))
colnames(prop_df) <- c("Group", "Cluster", "Freq")
prop_df$Group <- factor(prop_df$Group, levels=unique(prop_df$Group))

p_bar <- ggplot(prop_df, aes(x=Cluster, y=Freq, fill=Group)) +
  geom_bar(stat="identity", position="fill") +
  scale_y_continuous(labels=percent_format()) +
  scale_fill_manual(values= cell_type_cols) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(size=11, face="bold"),
        plot.title=element_text(size=13, face="bold", hjust=0.5)) +
  labs(x="Cluster", y="Proportion", title="Cluster Composition by Group")

p_bar


glycan_seu@active.ident <- as.factor(glycan_seu$phenotype)
prop_df <- as.data.frame(prop.table(table(Idents(glycan_seu), glycan_seu$seurat_clusters)))
colnames(prop_df) <- c("Group", "Cluster", "Freq")
prop_df$Group <- factor(prop_df$Group, levels=unique(prop_df$Group))

p_bar <- ggplot(prop_df, aes(x=Cluster, y=Freq, fill=Group)) +
  geom_bar(stat="identity", position="fill") +
  scale_y_continuous(labels=percent_format()) +
  scale_fill_manual(values=expand_palette(nature_palette2, length(levels(prop_df$Group)))) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(size=11, face="bold"),
        plot.title=element_text(size=13, face="bold", hjust=0.5)) +
  labs(x="Cluster", y="Proportion", title="Cluster Composition by Group")
p_bar


glycan_seu@active.ident <- as.factor(glycan_seu$loc)
prop_df <- as.data.frame(prop.table(table(Idents(glycan_seu), 
                                          glycan_seu$seurat_clusters)))
colnames(prop_df) <- c("Group", "Cluster", "Freq")
prop_df$Group <- factor(prop_df$Group, levels=unique(prop_df$Group))

p_bar <- ggplot(prop_df, aes(x=Cluster, y=Freq, fill=Group)) +
  geom_bar(stat="identity", position="fill") +
  scale_y_continuous(labels=percent_format()) +
  scale_fill_manual(values=c(
    "#1F77B4", "#2CA02C", "#FF7F0E"
  )) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(size=11, face="bold"),
        plot.title=element_text(size=13, face="bold", hjust=0.5)) +
  labs(x="Cluster", y="Proportion", title="Cluster Composition by Group")
p_bar


glycan_seu@active.ident <- as.factor(glycan_seu$group)
prop_df <- as.data.frame(prop.table(table(Idents(glycan_seu), 
                                          glycan_seu$seurat_clusters)))
colnames(prop_df) <- c("Group", "Cluster", "Freq")
prop_df$Group <- factor(prop_df$Group, levels=unique(prop_df$Group))

p_bar <- ggplot(prop_df, aes(x=Cluster, y=Freq, fill=Group)) +
  geom_bar(stat="identity", position="fill") +
  scale_y_continuous(labels=percent_format()) +
  scale_fill_manual(values=c(
    "#1F77B4", "#17BECF", "#2CA02C","#9ACD32", "#FF7F0E", "#F5A0A1"
  )) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(size=11, face="bold"),
        plot.title=element_text(size=13, face="bold", hjust=0.5)) +
  labs(x="Cluster", y="Proportion", title="Group Composition by clsuter")
p_bar
