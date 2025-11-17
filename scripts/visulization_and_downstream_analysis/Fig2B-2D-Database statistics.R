# Load genome completeness data
library(cowplot)
library(tidyverse)
library(ggplot2)
mgnify_info<-read.delim2("~/Mouse_gut/genomes-all_metadata.tsv",header = T)
mgnify_info <- mgnify_info %>%
  separate(Lineage, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(across(everything(), ~ sub("^[a-z]__", "", .)))
mgbc_info<-read.delim2("~/fsdownload/MGBC_md_26640.tsv",header = T)
mgbc_info <- mgbc_info %>%
  separate(`Taxonomy..GTDB.r95.`, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(across(everything(), ~ sub("^[a-z]__", "", .)))

refseq_info<-read.delim2("~/Mouse_gut/ncbi_dataset.tsv",header = T)
refseq_info<-subset(refseq_info,CheckM.completeness!="")
ml<-data.frame(Completeness= mgnify_info$Completeness,
               Database= rep("Mgnify",nrow(mgnify_info))
)
ml2<-data.frame(Completeness= mgbc_info$Completeness,
                Database= rep("MGBC",nrow(mgbc_info))
)
ml3<-data.frame(Completeness= refseq_info$CheckM.completeness,
                Database= rep("Refseq Bacteria",nrow(refseq_info))
)
  
ml<-rbind(ml,ml2,ml3)
ml$Completeness=as.numeric(ml$Completeness)
ml$group <- "Completeness"
p2<-ggplot(ml, aes(x = Database, y = Completeness, fill = Database)) +
  geom_boxplot(position = position_dodge(width = 0.8),size = 0.5,outlier.size = 0.1,width=0.5) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Completeness") + theme(
         panel.background = element_blank(),     
         panel.grid.major = element_blank(),     
         panel.grid.minor = element_blank(),     
         axis.line = element_line(color = "black", size = 0.5),  
         axis.ticks = element_line(color = "black")  
       )+
  theme(strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.position = "none")+
  scale_fill_manual(values =  c("#1f78b4","#fb9a99","#a6cee3")) + facet_wrap(~group, scales = "free_y")
write.csv(ml,file = "Completeness_compare.csv",row.names = F)
# Species number in each database
tb<-read.csv("~/Mouse_gut/Database_summary.csv",header = T)

p1<-ggplot(tb, aes(x = Database, y = Count, fill = Database)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),width=0.5) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Count") + theme(
         panel.background = element_blank(),     
         panel.grid.major = element_blank(),     
         panel.grid.minor = element_blank(),     
         axis.line = element_line(color = "black", size = 0.5),  
         axis.ticks = element_line(color = "black")  
       )+
  theme(strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none") +
  scale_fill_manual(values =  c("#1f78b4","#fb9a99","#a6cee3")) +
  facet_wrap(~Group, scales = "free_y")

plot_grid(p1, p2, 
          nrow = 1,              
          rel_widths = c(3, 1))


# Classified rate in each database
cl<-read.csv("~/Mouse_gut/classified_20250321.csv",header = T)
cl$genome <- factor(cl$genome, levels = c("MGBC", "MGnify", "Refseq"))
cl<-cl[cl$level %in% c("classified","Genus","Species"),]

p3<-ggplot(cl, aes(x = genome, y = value, fill = genome)) +
  geom_boxplot(position = position_dodge(width = 0.8),size = 0.5,outlier.size = 0.1,width=0.5) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Classified Rate %") + theme(
         panel.background = element_blank(),     
         panel.grid.major = element_blank(),     
         panel.grid.minor = element_blank(),     
         axis.line = element_line(color = "black", size = 0.5),  
         axis.ticks = element_line(color = "black")  
       )+
  theme(strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")+  
  geom_jitter(aes(color = genome),size = 1.5,shape = 16) + 
  scale_fill_manual(values =  c("#1f78b4","#fb9a99","#a6cee3")) + facet_wrap(~level, scales = "free_y",nrow = 1)+ scale_color_manual(values =  c("#1f78b4","#fb9a99","#a6cee3"))


species<-read.csv("~/Mouse_gut/genome_sample_unique_species_count.csv",header = T)
genus<-read.csv("~/Mouse_gut/genome_sample_unique_genus_count.csv",header = T)
species$Genome <- factor(species$Genome, levels = c("MGBC", "MGnify", "Refseq"))
genus$Genome <- factor(genus$Genome, levels = c("MGBC", "MGnify", "Refseq"))
colnames(species)<-colnames(genus)
tb<-rbind(species,genus)
colnames(tb)[3]<-"Uniq Count"
tb$Level<-c(rep("Species",36),rep("Genus",36))

p4<-ggplot(tb, aes(x = Genome, y = `Uniq Count`, fill = Genome)) +
  geom_boxplot(position = position_dodge(width = 0.8),size = 0.5,outlier.size = 0.1,width=0.5) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Uniq Count") + theme(
         panel.background = element_blank(),     
         panel.grid.major = element_blank(),     
         panel.grid.minor = element_blank(),     
         axis.line = element_line(color = "black", size = 0.5),  
         axis.ticks = element_line(color = "black")  
       )+
  theme(strip.text = element_text(size = 10),legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))+  
  geom_jitter(aes(color = Genome),size = 1.5,shape = 16) + 
  scale_fill_manual(values =  c("#1f78b4","#fb9a99","#a6cee3"))+
  scale_color_manual(values =  c("#1f78b4","#fb9a99","#a6cee3"))+
   facet_wrap(~Level,nrow = 1,scales = "free_y")
plot_grid(p3, p4, 
          nrow = 1,                
          rel_widths = c(3, 2))
