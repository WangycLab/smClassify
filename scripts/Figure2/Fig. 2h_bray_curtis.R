# Fig2h Bray-Curtis distance


# 1. Load packages ---------------------------------------------------------

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggplot2)
library(ggpubr)
library(cowplot)


# 2. Load scMicrobe object and sample metadata -----------------------------

sample_q <- readRDS("~/Mouse_gut/sce_bac_annotated_withoutUMIGeneFilter.rds")
sample_q <- subset(sample_q, subset = pass == TRUE)

metadata <- read.csv("~/Mouse_gut/sample_metadata.csv", header = TRUE)
metadata <- arrange(metadata, batch, Segment, Phenotype, sample)


# 3. Shared helper code ----------------------------------------------------

calc_ratio_per_sample <- function(seurat_obj, tax_col){
  
  sample_ids <- unique(seurat_obj@meta.data$orig.ident)
  
  res_list <- list()
  
  for(sid in sample_ids){
    
    sub <- seurat_obj[, seurat_obj@meta.data$orig.ident == sid]
    
    tab <- table(sub[[tax_col]])
    df <- as.data.frame(tab)
    
    df <- df[df[[tax_col]] != "Unknown sp.",]
    
    df$ratio <- df$Freq / sum(df$Freq)
    
    colnames(df)[1] <- "taxon"
    colnames(df)[3] <- sid
    
    res_list[[sid]] <- df[,c("taxon", sid)]
  }
  
  res <- Reduce(function(x,y) full_join(x,y, by = "taxon"), res_list)
  res[is.na(res)] <- 0
  
  return(res)
}

calc_bray_rep_pairs <- function(species_mat, label, threshold = 0){
  
  df <- species_mat %>%
    rowwise() %>%
    mutate(mean_abundance = mean(c_across(-taxon))) %>%
    ungroup() %>%
    filter(mean_abundance > threshold) %>%
    select(-mean_abundance)
  
  mat <- df %>%
    column_to_rownames("taxon") %>%
    t() %>%
    as.data.frame()
  
  bray <- vegdist(mat, method = "bray")
  bray_mat <- as.matrix(bray)
  
  sample_names <- rownames(bray_mat)
  rep_group <- sapply(strsplit(sample_names, "-"), function(x){
    paste(x[1:2], collapse = "-")
  })
  names(rep_group) <- sample_names
  
  df_long <- as.data.frame(as.table(bray_mat))
  colnames(df_long) <- c("sample1", "sample2", "distance")
  
  df_long <- df_long %>%
    mutate(
      rep1 = rep_group[sample1],
      rep2 = rep_group[sample2]
    ) %>%
    filter(rep1 == rep2, sample1 != sample2) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(sample1, sample2)), collapse = "_")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE)
  
  df_long$method <- label
  df_long$threshold <- threshold
  
  return(df_long)
}


# 4. Genus-level Bray-Curtis distance -------------------------------------

base_dir <- "~/fsdownload/kraken_result/Mgnify_g/"
taxonomy_data_list <- list()
report_files <- list.files(base_dir, pattern = "_sc_taxonomy_G\\.report$", full.names = TRUE)

for (file in report_files) {
  sample_name <- basename(file)
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  data$Sample_Name <- sample_name
  taxonomy_data_list[[sample_name]] <- data
}

all_taxonomy_genus <- do.call(rbind, taxonomy_data_list)
rownames(all_taxonomy_genus) <- 1:nrow(all_taxonomy_genus)
all_taxonomy_genus$Sample_Name <- gsub("_sc_taxonomy_G.report", "", all_taxonomy_genus$Sample_Name)
all_taxonomy_genus <- merge(all_taxonomy_genus, metadata, by.x = "Sample_Name", by.y = "sample")
all_taxonomy_genus <- mutate(all_taxonomy_genus, BC = paste(barcode, Segment, Phenotype, batch, sep = "_"))

new_genus <- all_taxonomy_genus$name[match(sample_q@meta.data$BC, all_taxonomy_genus$BC)]
sample_q <- AddMetaData(sample_q, metadata = new_genus, col.name = "kraken2_genus")

sample_q$smClassify_genus <- sample_q$Tax_Genus
sample_q$kraken2_genus <- sample_q@meta.data$kraken2_genus

sm_genus_mat <- calc_ratio_per_sample(sample_q, "smClassify_genus")
kr_genus_mat <- calc_ratio_per_sample(sample_q, "kraken2_genus")

df_genus_sm <- calc_bray_rep_pairs(sm_genus_mat, "smClassify", 0)
df_genus_kr <- calc_bray_rep_pairs(kr_genus_mat, "Kraken2", 0)
df_genus <- bind_rows(df_genus_kr, df_genus_sm)
df_genus$method <- factor(df_genus$method, levels = c("Kraken2", "smClassify"))

df_genus_wide <- df_genus %>%
  select(pair, rep1, sample1, sample2, method, distance) %>%
  pivot_wider(names_from = method, values_from = distance)

wilcox.test(
  df_genus_wide$smClassify,
  df_genus_wide$Kraken2,
  paired = TRUE,
  exact = FALSE,
  correct = TRUE
)

write.csv(df_genus_wide, file = "bc_dist_compare0_G.csv", row.names = FALSE)

p_genus <- ggplot(df_genus, aes(x = method, y = distance, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8), size = 0.5, outlier.size = 0.1, width = 0.5) +
  theme_minimal() +
  labs(
    title = "Genus",
    x = "",
    y = "Bray-Curtis distance"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    strip.text = element_text(size = 10),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_point(aes(color = method), size = 1.5, shape = 16) +
  scale_fill_manual(values = c("Kraken2" = "#1f78b4", "smClassify" = "#fb9a99")) +
  scale_color_manual(values = c("Kraken2" = "#1f78b4", "smClassify" = "#fb9a99")) +
  stat_compare_means(
    method = "wilcox.test",
    paired = TRUE,
    method.args = list(exact = FALSE, correct = TRUE),
    label = "p.format",
    label.x = 1.3,
    label.y = max(df_genus$distance) * 1.1
  )

p_genus


# 5. Species-level Bray-Curtis distance -----------------------------------

base_dir <- "~/fsdownload/kraken_result/Mgnify_s/"
taxonomy_data_list <- list()
report_files <- list.files(base_dir, pattern = "_sc_taxonomy\\.report$", full.names = TRUE)

for (file in report_files) {
  sample_name <- basename(file)
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  data$Sample_Name <- sample_name
  taxonomy_data_list[[sample_name]] <- data
}

all_taxonomy_species <- do.call(rbind, taxonomy_data_list)
rownames(all_taxonomy_species) <- 1:nrow(all_taxonomy_species)
all_taxonomy_species$Sample_Name <- gsub("_sc_taxonomy.report", "", all_taxonomy_species$Sample_Name)
all_taxonomy_species <- merge(all_taxonomy_species, metadata, by.x = "Sample_Name", by.y = "sample")
all_taxonomy_species <- mutate(all_taxonomy_species, BC = paste(barcode, Segment, Phenotype, batch, sep = "_"))

new_species <- all_taxonomy_species$name[match(sample_q@meta.data$BC, all_taxonomy_species$BC)]
sample_q <- AddMetaData(sample_q, metadata = new_species, col.name = "kraken2_species")

sample_q$smClassify_species <- sample_q$species
sample_q$kraken2_species <- sample_q@meta.data$kraken2_species

sm_species_mat <- calc_ratio_per_sample(sample_q, "smClassify_species")
kr_species_mat <- calc_ratio_per_sample(sample_q, "kraken2_species")

df_species_sm <- calc_bray_rep_pairs(sm_species_mat, "smClassify", 0)
df_species_kr <- calc_bray_rep_pairs(kr_species_mat, "Kraken2", 0)
df_species <- bind_rows(df_species_kr, df_species_sm)
df_species$method <- factor(df_species$method, levels = c("Kraken2", "smClassify"))

df_species_wide <- df_species %>%
  select(pair, rep1, sample1, sample2, method, distance) %>%
  pivot_wider(names_from = method, values_from = distance)

wilcox.test(
  df_species_wide$smClassify,
  df_species_wide$Kraken2,
  paired = TRUE,
  exact = FALSE,
  correct = TRUE
)

write.csv(df_species_wide, file = "bc_dist_compare0_S.csv", row.names = FALSE)

p_species <- ggplot(df_species, aes(x = method, y = distance, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8), size = 0.5, outlier.size = 0.1, width = 0.5) +
  theme_minimal() +
  labs(
    title = "Species",
    x = "",
    y = "Bray-Curtis distance"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    strip.text = element_text(size = 10),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_point(aes(color = method), size = 1.5, shape = 16) +
  scale_fill_manual(values = c("Kraken2" = "#1f78b4", "smClassify" = "#fb9a99")) +
  scale_color_manual(values = c("Kraken2" = "#1f78b4", "smClassify" = "#fb9a99")) +
  stat_compare_means(
    method = "wilcox.test",
    paired = TRUE,
    method.args = list(exact = FALSE, correct = TRUE),
    label = "p.format",
    label.x = 1.3,
    label.y = max(df_species$distance) * 1.1
  )

p_species


# 6. Combined panel --------------------------------------------------------

fig2h <- plot_grid(
  p_genus,
  p_species,
  nrow = 1
)

fig2h
