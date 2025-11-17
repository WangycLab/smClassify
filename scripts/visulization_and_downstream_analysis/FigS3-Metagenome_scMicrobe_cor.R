library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(cowplot)
library(tidyverse)
library(stringr)
library(purrr)
# Plot metagenome genus composition

metagenome<-read.csv("~/Mouse_gut/taxonomy_Genus_abund_Sample.csv",header = T)
df<-metagenome

df_long <- df %>%
  filter(!str_detect(Genus, "unclassified")) %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "Abundance")

# Calculate top 20 
top20 <- df_long %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

# Label genus not in top 20 as Others
df_plot <- df_long %>%
  mutate(Genus = if_else(Genus %in% top20_species, Genus, "Others")) %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Normalize by sample
df_plot <- df_plot %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance))
species_order <- df_plot %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  pull(Genus)

# Set stack order: Others at the bottom
df_plot$Genus <- factor(
  df_plot$Genus,
  levels = c(setdiff(species_order, "Others"),"Others")
)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#57C3F3', "#68A180",'#476D87',"#CCE0F5",
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175', '#D6E7A3')
# Plot
df_plot$Sample<-gsub("MC","Cecum",df_plot$Sample)
ggplot(df_plot, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Relative abundance", x = NULL, fill = "Genus")+
  scale_fill_manual(values=my36colors)+  
  guides(fill = guide_legend(reverse = FALSE))  


# Plot metagenome species composion
df<-read.csv("~/Mouse_gut/taxonomy_Species_abund_Sample.csv",header = T)

df_long <- df %>%
  filter(!str_detect(Species, "unclassified")) %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")

# Calculate top 20 species
top20_species <- df_long %>%
  group_by(Species) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 20) %>%
  pull(Species)

# Label species not in top 20 as Others
df_plot <- df_long %>%
  mutate(Species = if_else(Species %in% top20_species, Species, "Others")) %>%
  group_by(Sample, Species) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Normalize by sample
df_plot <- df_plot %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance))
species_order <- df_plot %>%
  group_by(Species) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  pull(Species)

# Set stack order: Others at the bottom
df_plot$Species <- factor(
  df_plot$Species,
  levels = c(setdiff(species_order, "Others"),"Others")
)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#57C3F3', "#68A180",'#476D87',"#CCE0F5",
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175', '#D6E7A3')
# Plot
df_plot$Sample<-gsub("MC","Cecum",df_plot$Sample)
ggplot(df_plot, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Relative abundance", x = NULL, fill = "Species")+
  scale_fill_manual(values=my36colors)+  
  guides(fill = guide_legend(reverse = FALSE))  



# Correlation between genus in metagenome and single microbe data
# Genus abundance in metagenome
metagenome<-read.csv("~/Mouse_gut/taxonomy_Genus_abund_Group.csv",header = T)
metagenome[[1]] <- sub("^.{3}", "", metagenome[[1]])
metagenome$MC_WT<-as.numeric(metagenome$MC_WT)
metagenome$MC_DB<-as.numeric(metagenome$MC_DB)
metagenome[,c(2:3)] <- sweep(metagenome[,c(2:3)], 2, colSums(metagenome[,c(2:3)]), FUN = "/")
colnames(metagenome)<-c("Name","W","D")


# Genus abundance in single microbe
sample_q<-readRDS("~/fsdownload/sce_bac_annotated.rds")
sample_q <- sample_q[ , grepl("Cecum", sample_q@meta.data$orig.ident) ]
sample_q <- sample_q[ , grepl("DB", sample_q@meta.data$orig.ident) ]
genus_count_D<-table(sample_q$Tax_Genus) %>% as.data.frame()
sample_q<-readRDS("~/fsdownload/sce_bac_annotated.rds")
sample_q <- sample_q[ , grepl("Cecum", sample_q@meta.data$orig.ident) ]
sample_q <- sample_q[ , grepl("WT", sample_q@meta.data$orig.ident) ]
genus_count_W<-table(sample_q$Tax_Genus) %>% as.data.frame()

count_all<-merge(genus_count_D,genus_count_W,by="Var1",all=T)
count_all<-count_all[,c(1,3,2)]
colnames(count_all)<-c("Name","W","D")
count_all[is.na(count_all)]<-0
count_all[,c(2:3)] <- sweep(count_all[,c(2:3)], 2, colSums(count_all[,c(2:3)]), FUN = "/")
colnames(count_all)<-colnames(metagenome)
count_all$Name<-gsub("_[A-Z]","",count_all$Name)

# Unify genus info
mgnify_info<-read.delim2("~/Mouse_gut/genomes-all_metadata.tsv",header = T)
mgnify_info <- mgnify_info %>%
  separate(Lineage, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(across(everything(), ~ sub("^[a-z]__", "", .)))
mgnify_info$s<-ifelse(mgnify_info$s =="",mgnify_info$Genome,mgnify_info$s)
info<-mgnify_info[,c(18:20)]%>%distinct()
info<-subset(info,f!="")
info1<-info[info$g %in% setdiff(count_all$Name,metagenome$Name),]
metagenome[grepl(paste(info1$f, collapse = "|"), metagenome$Name), ]

update_word <- function(df1, df2, family = "Lachnospiraceae") {

  df_merged <- merge(df1, df2, by.x = "Name", by.y = "g", all.x = TRUE)
  #
  df_merged$Name <- ifelse(!is.na(df_merged$f) & df_merged$f == family & (df_merged$Name %in% setdiff(count_all$Name,metagenome$Name)),
                        paste0(family, "_unclassified"),
                        df_merged$Name)
  
  df1_new <- df_merged[, names(df1)]
  return(df1_new)
}
count_all1 <- update_word(count_all, info, family = "Lachnospiraceae")
count_all1 <- update_word(count_all1, info, family = "Muribaculaceae")
count_all1 <- update_word(count_all1, info, family = "Oscillospiraceae")
count_all1 <- aggregate(. ~ Name, data = count_all1, FUN = sum)

# Merge genus composition of metagenome and single microbe data
df<-merge(count_all1,metagenome,by="Name")
df[,c(2:5)] <- sweep(df[,c(2:5)], 2, colSums(df[,c(2:5)]), FUN = "/")

#df <- df %>%
#  mutate(Name = ifelse(if_any(everything(), is.na), "Others", Name))
#df[is.na(df)]<-0
#df <- aggregate(. ~ Name, data = df, FUN = sum)

rownames(df)<-df$Name
df<-df[,-1]
#df <- df[rowSums(df > 0.001) >= 1, ]

# Plot scatter and calculate correlation
p1<-ggscatter(df, 
          x = "W.x", 
          y = "W.y",
          color = "#1f78b4", shape = 20, size = 1, 
          add = "reg.line",  
          add.params = list(color = "#cc3333", fill = "lightgray", size = 0.4), 
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          cor.coef.size = 5,
          ggtheme = theme_minimal(),
          xlab = "Single microbe", 
          ylab = "Metagenome"
) +
  theme(
    panel.grid = element_blank(),        
    axis.line = element_line(),          
    axis.ticks = element_line()          
  ) 

p2<-ggscatter(df, 
              x = "D.x", 
              y = "D.y",
              color = "#1f78b4", shape = 20, size = 1, 
              add = "reg.line",  
              add.params = list(color = "#cc3333", fill = "lightgray", size = 0.4), 
              cor.coef = TRUE,
              cor.coeff.args = list(method = "pearson", label.sep = "\n"),
              cor.coef.size = 5,
              ggtheme = theme_minimal(),
              xlab = "Single microbe", 
              ylab = "Metagenome"
) +
  theme(
    panel.grid = element_blank(),        
    axis.line = element_line(),         
    axis.ticks = element_line()         
  ) 
