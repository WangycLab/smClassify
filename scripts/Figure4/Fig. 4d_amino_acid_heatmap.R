############################################################
## Figure 4d amino-acid metabolite heatmap

# ============================================================
# 1. Load packages
# ============================================================
library(tidyverse)
library(readxl)
library(pheatmap)
library(ggpubr)

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
# 3. Amino-acid metabolite list
# ============================================================
amino_acids <- c(
  "Alanine","Aspartic acid","Glutamic acid","Glutamine","Isoleucine","Leucine",
  "Lysine","Methionine","Phenylalanine","Proline","Serine","Threonine",
  "Tryptophan","Tyrosine","Valine","Ornithine","Citrulline","Arginine"
)

# ============================================================
# 4. Read metabolomics data
# ============================================================
feature_cols <- c(
  "Type","ID","MZ","RT","MS2.name","MS2_score",
  "Formula","SuperClass","Class","Subclass",
  "HMDB","kegg","kegg_pathway","MS1.name"
)

serum <- read_xlsx(serum_file)
content <- read_xlsx(content_file)

names(serum) <- str_replace_all(names(serum), "^Seurm", "Serum")

serum_long <- serum |>
  pivot_longer(-all_of(feature_cols),
               names_to="Sample",
               values_to="Intensity")

content_long <- content |>
  pivot_longer(-all_of(feature_cols),
               names_to="Sample",
               values_to="Intensity")

dat <- bind_rows(serum_long, content_long)

dat <- dat %>%
  mutate(
    Metabolite = ifelse(is.na(MS2.name) | !nzchar(MS2.name), MS1.name, MS2.name),
    Treatment = case_when(
      str_detect(Sample, "WT") ~ "WT",
      str_detect(Sample, "DB") ~ "DB",
      TRUE ~ NA_character_
    ),
    Tissue = case_when(
      str_detect(Sample, "Serum") ~ "Serum",
      str_detect(Sample, "MC") ~ "Cecum",
      str_detect(Sample, "JC") ~ "Colon",
      str_detect(Sample, "ZC") ~ "Rectum",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Treatment)) %>%
  mutate(
    Metabolite = str_trim(Metabolite),
    Tissue = factor(Tissue, levels=c("Cecum","Colon","Rectum","Serum")),
    Treatment = factor(Treatment, levels=c("WT","DB")),
    Intensity = as.numeric(Intensity)
  )

# ============================================================
# 5. Build heatmap matrix
# ============================================================
make_heatmap_matrix <- function(df, metabolite_vec){
  
  desired_cols <- c(
    "Cecum_WT","Cecum_DB",
    "Colon_WT","Colon_DB",
    "Rectum_WT","Rectum_DB",
    "Serum_WT","Serum_DB"
  )
  
  mat_df <- df %>%
    filter(Metabolite %in% metabolite_vec) %>%
    filter(!is.na(Tissue), !is.na(Treatment)) %>%
    group_by(Metabolite, Tissue, Treatment) %>%
    summarise(mean_intensity = mean(Intensity, na.rm = TRUE), .groups="drop") %>%
    mutate(Group = paste(Tissue, Treatment, sep="_")) %>%
    select(Metabolite, Group, mean_intensity) %>%
    pivot_wider(names_from=Group, values_from=mean_intensity)
  
  mat_df <- tibble(Metabolite = metabolite_vec) %>%
    left_join(mat_df, by="Metabolite")
  
  for(cc in desired_cols){
    if(!cc %in% colnames(mat_df)){
      mat_df[[cc]] <- NA_real_
    }
  }
  
  mat_df <- mat_df %>%
    select(Metabolite, all_of(desired_cols))
  
  mat <- mat_df %>%
    as.data.frame()
  
  rownames(mat) <- mat$Metabolite
  mat$Metabolite <- NULL
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  keep <- apply(mat,1,function(x) sum(is.finite(x)) >= 2)
  mat <- mat[keep,,drop=FALSE]
  
  if(nrow(mat)==0) return(NULL)
  
  mat_z <- t(apply(mat,1,function(x){
    
    idx <- is.finite(x)
    
    if(sum(idx)<2){
      return(rep(NA,length(x)))
    }
    
    m <- mean(x[idx])
    s <- sd(x[idx])
    
    if(s==0){
      out <- rep(0,length(x))
      out[!idx] <- NA
      return(out)
    }
    
    out <- (x-m)/s
    out[!idx] <- NA
    out
  }))
  
  rownames(mat_z) <- rownames(mat)
  colnames(mat_z) <- colnames(mat)
  
  mat_z
}

# ============================================================
# 6. Export source data tables
# ============================================================
export_source_data <- function(df, metabolite_vec, prefix){
  
  desired_cols <- c(
    "Cecum_WT","Cecum_DB",
    "Colon_WT","Colon_DB",
    "Rectum_WT","Rectum_DB",
    "Serum_WT","Serum_DB"
  )
  
  source_long <- df %>%
    filter(Metabolite %in% metabolite_vec) %>%
    filter(!is.na(Tissue), !is.na(Treatment)) %>%
    select(
      Sample,
      Tissue,
      Treatment,
      Metabolite,
      Intensity,
      Type,
      ID,
      MZ,
      RT,
      MS2.name,
      MS2_score,
      Formula,
      SuperClass,
      Class,
      Subclass,
      HMDB,
      kegg,
      kegg_pathway,
      MS1.name
    ) %>%
    arrange(Metabolite, Tissue, Treatment, Sample)
  
  write_csv(
    source_long,
    file.path(output_dir, paste0(prefix, "_source_data_long.csv"))
  )
  
  source_mean <- source_long %>%
    group_by(Metabolite, Tissue, Treatment) %>%
    summarise(
      mean_intensity = mean(Intensity, na.rm = TRUE),
      n_samples = sum(!is.na(Intensity)),
      .groups = "drop"
    ) %>%
    mutate(Group = paste(Tissue, Treatment, sep="_")) %>%
    arrange(Metabolite, Tissue, Treatment)
  
  write_csv(
    source_mean,
    file.path(output_dir, paste0(prefix, "_source_data_group_mean_long.csv"))
  )
  
  mean_matrix_df <- source_mean %>%
    select(Metabolite, Group, mean_intensity) %>%
    pivot_wider(names_from=Group, values_from=mean_intensity)
  
  mean_matrix_df <- tibble(Metabolite = metabolite_vec) %>%
    left_join(mean_matrix_df, by="Metabolite")
  
  for(cc in desired_cols){
    if(!cc %in% colnames(mean_matrix_df)){
      mean_matrix_df[[cc]] <- NA_real_
    }
  }
  
  mean_matrix_df <- mean_matrix_df %>%
    select(Metabolite, all_of(desired_cols))
  
  write_csv(
    mean_matrix_df,
    file.path(output_dir, paste0(prefix, "_source_data_group_mean_matrix.csv"))
  )
  
  invisible(
    list(
      source_long = source_long,
      source_mean = source_mean,
      mean_matrix = mean_matrix_df
    )
  )
}

# ============================================================
# 7. Plot heatmap
# ============================================================
plot_heatmap <- function(metabolite_vec, prefix){
  
  mat_z <- make_heatmap_matrix(dat, metabolite_vec)
  
  if(is.null(mat_z)){
    message("No valid data for ", prefix)
    return(NULL)
  }
  
  ann_col <- data.frame(
    Tissue = factor(sub("_.*","",colnames(mat_z)),
                    levels=c("Cecum","Colon","Rectum","Serum")),
    Treatment = factor(sub(".*_","",colnames(mat_z)),
                       levels=c("WT","DB"))
  )
  
  rownames(ann_col) <- colnames(mat_z)
  
  ann_colors <- list(
    Treatment = c(
      WT = "#4DBBD5",DB = "#D94F45"
    ),
    Tissue = c(
      Cecum = "#1F77B4",
      Colon = "#2CA02C",
      Rectum = "#FF7F0E",
      Serum = "#9467BD"
    )
  )
  
  mat_cluster <- mat_z
  for(i in 1:nrow(mat_cluster)){
    row_vals <- mat_cluster[i,]
    idx <- is.finite(row_vals)
    if(any(idx)){
      row_mean <- mean(row_vals[idx])
      row_vals[!idx] <- row_mean
      mat_cluster[i,] <- row_vals
    }
  }
  
  row_hc <- hclust(dist(mat_cluster))
  
  render_check <- tibble(
    prefix = prefix,
    n_metabolites_plotted = nrow(mat_z),
    n_groups_plotted = ncol(mat_z),
    n_finite_values = sum(is.finite(mat_z)),
    pdf_file = file.path(output_dir,paste0(prefix,"_heatmap_F.pdf")),
    png_file = file.path(output_dir,paste0(prefix,"_heatmap_F.png"))
  )
  
  write_csv(
    render_check,
    file.path(output_dir, paste0(prefix, "_heatmap_render_check.csv"))
  )
  
  if(nrow(mat_z) == 0 || ncol(mat_z) == 0 || sum(is.finite(mat_z)) == 0){
    stop("No plottable values for ", prefix, ". Please check the source data and metabolite names.")
  }
  
  pheatmap::pheatmap(
    mat_z,
    cluster_rows=row_hc,
    cluster_cols=FALSE,
    annotation_col=ann_col,
    annotation_colors=ann_colors,
    scale="none",
    border_color=NA,
    fontsize_row=10,
    fontsize_col=10,
    main=prefix,
    na_col="grey90",
    filename=file.path(output_dir,paste0(prefix,"_heatmap_F.pdf")),
    width=6,
    height=max(5,0.35*nrow(mat_z)+2)
  )
  
  pheatmap::pheatmap(
    mat_z,
    cluster_rows=row_hc,
    cluster_cols=FALSE,
    annotation_col=ann_col,
    annotation_colors=ann_colors,
    scale="none",
    border_color=NA,
    fontsize_row=10,
    fontsize_col=10,
    main=prefix,
    na_col="grey90",
    filename=file.path(output_dir,paste0(prefix,"_heatmap_F.png")),
    width=6,
    height=max(5,0.35*nrow(mat_z)+2)
  )
  
  write.csv(
    mat_z,
    file.path(output_dir, paste0(prefix, "_heatmap_zscore_matrix.csv")),
    row.names = TRUE
  )
  
  invisible(mat_z)
}

# ============================================================
# 8. Check metabolite detection
# ============================================================
all_mets <- unique(dat$Metabolite)

cat("Total metabolites detected:", length(all_mets), "\n")
cat("Amino acids requested:", length(amino_acids), "\n")
cat("Amino acids detected:", sum(amino_acids %in% all_mets), "\n")

amino_acid_detection_summary <- tibble(
  Metabolite = amino_acids,
  Detected = amino_acids %in% all_mets
)

write_csv(
  amino_acid_detection_summary,
  file.path(output_dir, "AminoAcids_detection_summary.csv")
)

fig4d_summary <- tibble(
  total_metabolites_detected = length(all_mets),
  amino_acids_requested = length(amino_acids),
  amino_acids_detected = sum(amino_acids %in% all_mets)
)

write_csv(
  fig4d_summary,
  file.path(output_dir, "Fig4d_amino_acid_heatmap_summary.csv")
)

# ============================================================
# 9. Export source data and draw Figure 4d heatmap
# ============================================================
fig4d_source_data <- export_source_data(df = dat, metabolite_vec = amino_acids, prefix = "AminoAcids")
fig4d_heatmap_matrix <- plot_heatmap(amino_acids, "AminoAcids")

