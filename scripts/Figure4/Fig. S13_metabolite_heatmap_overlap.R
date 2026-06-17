############################################################
## Figure S13 metabolite heatmap and overlap analysis

############################################################

# ============================================================
# 1. Load packages
# ============================================================

library(tidyverse)
library(readxl)
library(limma)
library(pheatmap)
library(UpSetR)

# ============================================================
# 2. Input and output paths
# ============================================================
base_dir   <- "Figure4"
input_dir  <- file.path(base_dir, "Input")
output_dir <- file.path(base_dir, "Output")

serum_file   <- file.path(input_dir, "serum_norm_ms2_metabolites_intensity.xlsx")
content_file <- file.path(input_dir, "content_norm_ms2_metabolites_intensity.xlsx")

dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 3. Colours
# ============================================================
pal_pheno <- c(
  WT = "#4DBBD5",
  DB = "#E64B35"
)

pal_region <- c(
  Cecum  = "#1f77b4",
  Colon  = "#2ca02c",
  Rectum = "#ff7f0e"
)

heat_cols <- colorRampPalette(c("#2166ac","white","#b2182b"))(100)

ann_colors <- list(
  Tissue = pal_region,
  Genotype = pal_pheno
)

# ============================================================
# 4. Read metabolomics data
# ============================================================
serum   <- read_xlsx(serum_file)
content <- read_xlsx(content_file)

names(serum) <- str_replace_all(names(serum), "^Seurm", "Serum")

feature_cols <- c(
  "Type","ID","MZ","RT","MS2.name","MS2_score",
  "Formula","SuperClass","Class","Subclass",
  "HMDB","kegg","kegg_pathway","MS1.name"
)

serum_long <- serum %>%
  pivot_longer(-all_of(feature_cols),
               names_to="Sample",
               values_to="Intensity")

content_long <- content %>%
  pivot_longer(-all_of(feature_cols),
               names_to="Sample",
               values_to="Intensity")

all_long <- bind_rows(serum_long,content_long)

# ============================================================
# 5. Parse sample metadata
# ============================================================
meta <- all_long %>%
  distinct(Sample) %>%
  mutate(
    Tissue_raw = str_extract(Sample,"^[A-Za-z]+"),
    Genotype   = str_extract(Sample,"(WT|DB)"),
    Rep        = as.integer(str_extract(Sample,"\\d{1,2}$")),
    Tissue = case_when(
      Tissue_raw %in% c("MC","Cecum") ~ "Cecum",
      Tissue_raw %in% c("JC","Colon") ~ "Colon",
      Tissue_raw %in% c("ZC","Rectum") ~ "Rectum",
      Tissue_raw %in% c("Serum") ~ "Serum",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!str_detect(Sample,"^QC")) %>%
  filter(!is.na(Tissue),!is.na(Genotype)) %>%
  mutate(
    Tissue=factor(Tissue,levels=c("Cecum","Colon","Rectum","Serum")),
    Genotype=factor(Genotype,levels=c("WT","DB"))
  )

write_csv(meta,file.path(output_dir,"sample_metadata.csv"))

all_long <- all_long %>%
  filter(Sample %in% meta$Sample)

# ============================================================
# 6. Build metabolite intensity matrix
# ============================================================
expr_mat <- all_long %>%
  filter(!is.na(MS2.name)) %>%
  group_by(MS2.name,Sample) %>%
  summarise(Intensity=max(Intensity,na.rm=TRUE),.groups="drop") %>%
  pivot_wider(names_from=Sample,values_from=Intensity) %>%
  distinct(MS2.name,.keep_all=TRUE) %>%
  column_to_rownames("MS2.name") %>%
  as.matrix()

expr_log <- log2(expr_mat+1)

# ============================================================
# 7. Helper functions
# ============================================================
impute_median <- function(mat){
  t(apply(mat,1,function(x){
    if(all(is.na(x))) return(x)
    x[is.na(x)]<-median(x,na.rm=TRUE)
    x
  }))
}

clean_matrix <- function(mat){
  mat <- mat[rowSums(is.na(mat)) < ncol(mat),]
  mat <- impute_median(mat)
  mat <- mat[apply(mat,1,sd)>0,]
  mat
}

# ============================================================
# 8. Differential metabolite analysis
# ============================================================
tissues <- c("Serum","Cecum","Colon","Rectum")
diff_list <- list()

for(tis in tissues){
  
  cat("Processing:",tis,"\n")
  
  samp <- meta %>% filter(Tissue==tis) %>% pull(Sample)
  
  mat <- expr_log[,samp,drop=FALSE]
  mat <- clean_matrix(mat)
  
  geno <- meta$Genotype[match(colnames(mat),meta$Sample)]
  
  design <- model.matrix(~0+geno)
  colnames(design) <- c("WT","DB")
  
  fit <- lmFit(mat,design)
  cont <- makeContrasts(DB-WT,levels=design)
  
  fit2 <- contrasts.fit(fit,cont)
  fit2 <- eBayes(fit2)
  
  tab <- topTable(fit2,number=Inf) %>%
    rownames_to_column("Metabolite")
  
  write_csv(tab,file.path(output_dir,paste0(tis,"_DBvsWT_limma.csv")))
  
  sig <- tab %>%
    filter(P.Value<0.05,abs(logFC)>0.5) %>%
    pull(Metabolite)
  
  diff_list[[tis]] <- sig
}

# ============================================================
# 9. Heatmap annotation
# ============================================================
luminal_meta <- meta %>%
  filter(Tissue %in% c("Cecum","Colon","Rectum")) %>%
  arrange(Tissue,Genotype,Rep)

ann_col <- luminal_meta %>%
  select(Sample,Tissue,Genotype) %>%
  column_to_rownames("Sample")

samples <- luminal_meta$Sample

# ============================================================
# 10. Heatmap of all luminal metabolites
# ============================================================
mat_all <- expr_log[,samples]
mat_all <- clean_matrix(mat_all)
mat_all <- mat_all[,rownames(ann_col)]

pdf(file.path(output_dir,"Heatmap_AllMetabolites_T.pdf"),8,10)

pheatmap(
  mat_all,
  scale="row",
  color=heat_cols,
  annotation_col=ann_col,
  annotation_colors=ann_colors,
  cluster_cols=TRUE,
  show_rownames=FALSE,
  show_colnames=FALSE,
  border_color=NA,
  main="All luminal metabolites"
)

dev.off()

# ============================================================
# 11. Heatmap of significant luminal metabolites
# ============================================================
sig_luminal <- unique(unlist(diff_list[c("Cecum","Colon","Rectum")]))

mat_sig <- expr_log[rownames(expr_log)%in%sig_luminal,samples]
mat_sig <- clean_matrix(mat_sig)
mat_sig <- mat_sig[,rownames(ann_col)]

colnames(mat_sig)

pdf(file.path(output_dir,"Heatmap_SignificantMetabolites_T.pdf"),6,8)

pheatmap(
  mat_sig,
  scale="row",
  color=heat_cols,
  annotation_col=ann_col,
  annotation_colors=ann_colors,
  cluster_cols=TRUE,
  border_color=NA,
  show_rownames=FALSE,
  show_colnames=T,
  fontsize_row=6,
  main="Significant differential metabolites across intestinal regions"
)

dev.off()

# ============================================================
# 12. Overlap plot of differential luminal metabolites
# ============================================================
diff_list_lumen <- diff_list[c("Cecum","Colon","Rectum")]

up_in <- fromList(diff_list_lumen)

png(file.path(output_dir,"DifferentialMetabolites_UpSet.png"),
    width=2600,height=2200,res=350)

UpSetR::upset(
  up_in,
  sets=c("Cecum","Colon","Rectum"),
  order.by="freq",
  sets.bar.color="#4DBBD5",
  main.bar.color="#E64B35"
)

dev.off()

