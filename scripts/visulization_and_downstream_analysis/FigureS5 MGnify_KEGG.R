library(dplyr)
library(tidyverse)
data<-read.delim2("~/Mouse_gut/genomes-all_metadata.tsv",header = T)
data <- data %>%
  separate(Lineage, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(across(everything(), ~ sub("^[a-z]__", "", .)))




species_info<- data %>% distinct(Genome,p, s, .keep_all = FALSE)
phylum_species<-data %>% distinct(p, s, .keep_all = FALSE)


# Species corresponding to top 5 phyla
eggnog<-read.delim2("~/Mouse_gut/M1_eggNOG.tsv",header = T,check.names = F)
eggnog$species<-substr(eggnog$`#query`,1,13)
eggnog<-merge(eggnog,species_info,by.x="species",by.y="Genome")
eggnog_freq<-table(eggnog$p)%>%as.data.frame()%>%arrange(desc(Freq))%>%head(8)
df2<-species_info[species_info$p %in% eggnog_freq$Var1,]

eggnog<-eggnog[eggnog$species %in% df2$Genome ,]
eggnog$`#query`<-gsub("_","-",eggnog$`#query`)
eggnog1<-eggnog#[eggnog$`#query` %in% genes_in_sample ,]
#eggnog<-eggnog[eggnog$s %in% metadata$species_info ,]

# KEGG
KEGG_table <- eggnog1 %>%
  dplyr::select(KEGG_Pathway, `#query`,species) %>%
  filter(KEGG_Pathway != "-") %>%
  dplyr::mutate(KEGG_Pathway = strsplit(KEGG_Pathway, ',')) %>%
  tidyr::unnest()
KEGG_table<-KEGG_table[grepl("map", KEGG_table$KEGG_Pathway),]
KEGG_table<-merge(KEGG_table,species_info,by.x="species",by.y="Genome")
KEGG_table<-arrange(KEGG_table,KEGG_Pathway)
library(clusterProfiler)
koid<-ko2name(KEGG_table$KEGG_Pathway)
KEGG_table$Description<-koid$name
KEGG_table$s<-ifelse(KEGG_table$s=="",KEGG_table$species,KEGG_table$s)
KEGG_table<-subset(KEGG_table, Description!="" & s !="")

interested_pathway<-c("Metabolic pathways","Propanoate metabolism",
                      "Carbon metabolism","Pyruvate metabolism",
                      "Butanoate metabolism","Starch and sucrose metabolism",
                      "Biosynthesis of amino acids","Galactose metabolism",
                      "Valine, leucine and isoleucine biosynthesis",
                      "Falty acid metabolism","Fatty acid biosynthesis")
species_description_count <- KEGG_table %>% dplyr::count(KEGG_Pathway,Description,s,p)
species_description_count<-species_description_count[species_description_count$Description %in% interested_pathway,]
species_description_count$p<-factor(species_description_count$p,levels = eggnog_freq$Var1)
species_description_count<-arrange(species_description_count,p)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)  
library(dplyr)
library(RColorBrewer)
# Convert data format suitable for heatmap
data_matrix <- species_description_count %>%
  select(Description, s, n) %>%
  pivot_wider(names_from = s, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Description") %>%
  as.matrix()
zscore_matrix <- t(scale(t(data_matrix)))  # Transpose, standardize by column (original rows), then transpose back
log_matrix <- log10(data_matrix+1)

# Create annotation colors for X axis (s)
annotation_df <- species_description_count %>% distinct(s, p)  # Keep only unique s-p combinations
annotation_df <- column_to_rownames(annotation_df, var = "s")

# Generate color mapping (assign different colors to each `p`)
annotation_colors <- structure(
  brewer.pal(length(unique(annotation_df$p)), "Set3"), 
  names = as.character(unique(annotation_df$p))
)

# Create X axis annotation
col_anno <- columnAnnotation(
  Phylum = annotation_df[colnames(data_matrix), "p"],
  col = list(Phylum = annotation_colors),
  annotation_legend_param = list(title = "Phylum")
)

# Plot heatmap
Heatmap(log_matrix, name = "Log10 count",
        col = colorRamp2(c(1, max(log_matrix)), c("white", "red")),
        cluster_rows = TRUE, cluster_columns = F,show_column_names = F,
        top_annotation = col_anno,row_names_gp = gpar(fontsize = 10))  # Add X axis annotation

