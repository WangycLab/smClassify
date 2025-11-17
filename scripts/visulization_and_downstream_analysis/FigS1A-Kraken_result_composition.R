library(tidyverse)
library(dplyr)
refseq_mapping<-read.delim2("~/fsdownload/refseq_mapping.tsv",header = F)
mgnify_info<-read.delim2("~/Mouse_gut/genomes-all_metadata.tsv",header = T)
mgnify_info <- mgnify_info %>%
  separate(Lineage, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(across(everything(), ~ sub("^[a-z]__", "", .)))

refseq_mapping <- refseq_mapping %>%
  mutate(species = str_extract(V3, "^\\S+\\s+\\S+"))
refseq_mapping$genus <- sub(" .*", "", refseq_mapping$V3)

base_dir <- "~/fsdownload/kraken_result/Refseq_s/"  
# Initialize list to store data
taxonomy_data_list <- list()
report_files <- list.files(base_dir, pattern = "_sc_taxonomy\\.report$", full.names = TRUE)
report_files <- report_files[grepl("412", report_files)]
for (file in report_files) {
  sample_name <- basename((file))  # Get sample folder name
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,quote = "")  # Read data
  data$Sample_Name <- sample_name  # Add sample name column
  taxonomy_data_list[[sample_name]] <- data  # Store in list
}

if(T){
all_taxonomy_data <- do.call(rbind, taxonomy_data_list)
rownames(all_taxonomy_data)<-1:nrow(all_taxonomy_data)
all_taxonomy_data$Sample_Name<-gsub("_sc_taxonomy.report","",all_taxonomy_data$Sample_Name)
}
#anx
if(T){
  df <- read.table("~/fsdownload/refseq_in_mgnify.result", header = FALSE, 
                   sep = "\t", stringsAsFactors = FALSE)
# Extract GCF_xxx.y from column 1 (first two underscore parts)
  df$id1 <- sapply(df$V1, function(x) {
    parts <- strsplit(basename(x), "_")[[1]]
    paste(parts[1], parts[2], sep = "_")
  })
  
# Extract MGYG ID from column 2 (remove path and suffix)
  df$id2 <- sapply(df$V2, function(x) {
    sub(".fna$", "", basename(x))
  })
  df<-merge(df,refseq_mapping,by.x="id1",by.y="V1")
  mgnify_info <- unique(mgnify_info[, c("Species_rep","g" ,"s")])
  df<-merge(df,mgnify_info,by.x="id2",by.y="Species_rep")
  anx95<-subset(df,V3.x>95)
}

all_taxonomy_data <- all_taxonomy_data %>%
  mutate(selected_name = case_when(
    name %in% anx95$V3.y | name %in% anx95$species ~ name,  # If name1 is in the list, choose it
    name2 %in% anx95$V3.y| name2 %in% anx95$species ~ name2,  # Else if name2 is in the list, choose it
    TRUE ~ NA_character_            # Otherwise set to NA
  ),selected_puri = case_when(
    name %in% anx95$V3.y | name %in% anx95$species~ fraction_total_reads,  # If name1 is in the list, choose it
    name2 %in% anx95$V3.y | name2 %in% anx95$species ~ fraction_total_reads2,  # Else if name2 is in the list, choose it
    TRUE ~ NA_real_            # Otherwise set to NA
  )
  )%>%
  filter(!is.na(selected_name)) 

# Save retained barcode ID table

all_taxonomy_data1<-subset(all_taxonomy_data,selected_puri>0.1 & kraken_assigned_reads>50)
# Add rules for retaining names
anx95 <- anx95 %>%
  mutate(species = ifelse(
    V3.y %in% setdiff(all_taxonomy_data1$selected_name, anx95$species),
    V3.y,
    species
  ))

anx95<-anx95[,c(1,2,5,6,7,9,10,13)]
anxmax <- anx95 %>%
  group_by(species) %>%                    
  slice_max(order_by = V3.x, n = 1, with_ties = FALSE) %>%  
  ungroup()
all_taxonomy_data1<-merge(all_taxonomy_data1,anxmax,by.x="selected_name",
                          by.y = "species")
all_taxonomy_data1$selected_name<-all_taxonomy_data1$s

all_taxonomy_data1$barcode<-paste(all_taxonomy_data1$barcode,all_taxonomy_data1$Sample_Name,sep = "_")
all_taxonomy_data1<-arrange(all_taxonomy_data1,Sample_Name,barcode)
tb<-all_taxonomy_data1[,c(2,1)]
all_taxonomy_data1$selected_name <- gsub("_[A-Z]+(\\s|$)", "\\1", all_taxonomy_data1$selected_name)
write.table(tb,file = "~/Mouse_gut/barcode_species.txt",row.names = F,quote = F,sep = '\t')
if(T){
  all_taxonomy_data1$selected_name<-gsub("\\[","",all_taxonomy_data1$selected_name)
  all_taxonomy_data1$selected_name<-gsub("\\]","",all_taxonomy_data1$selected_name)
}


top_names <- all_taxonomy_data1 %>%
  count(selected_name, Sample_Name) %>%
  group_by(selected_name) %>%
  summarise(total_count = sum(n)) %>%
  arrange(desc(total_count))%>%
  slice_head(n = 15) %>%
  pull(selected_name)  # Extract top 15 `selected_name` list

all_taxonomy_data1$selected_g<-sub(" .*", "", all_taxonomy_data1$selected_name)
# Mark which are top15
all_taxonomy_data1 <- all_taxonomy_data1 %>%
  mutate(
    is_top = selected_name %in% top_names
  )

# Extract all genera of top15 (for genus-level matching)
top_genus <- all_taxonomy_data1 %>%
  filter(is_top) %>%
  pull(selected_g) %>%
  unique()

# Generate cleaned_name
all_taxonomy_data1 <- all_taxonomy_data1 %>%
  mutate(
    cleaned_name = case_when(
      selected_name %in% top_names ~ selected_name,
      selected_g %in% top_genus ~ paste0("Other_", selected_g),
      TRUE ~ "Other species"
    )
  )
# Classify names outside top 15 as "Other"
all_taxonomy_data1 <- all_taxonomy_data1 %>%
  mutate(selected_g = ifelse(selected_g %in% top_genus, selected_g, "Other genus"))
# Inner rings represent for genus
inner <- all_taxonomy_data1 %>%
  group_by(selected_g) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count)) %>%
  mutate(
    ring = 1,
    group = selected_g,
    genus_order = row_number()
  )
inner <- inner %>%
  mutate(order = ifelse(group == "Other genus", 999, row_number())) %>%
  arrange(order) %>%
  mutate(group = factor(group, levels = group))

# Outer ring: group by genus order, preserve intra-group order
outer <- all_taxonomy_data1 %>%
  group_by(selected_g, cleaned_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(inner %>% select(selected_g = group, genus_order = order), by = "selected_g") %>%
  mutate(
    ring = 2,
    group = cleaned_name,
    is_other = grepl("^Other_", cleaned_name)
  ) %>%
  arrange(genus_order, is_other, desc(count)) %>%
  mutate(group = factor(group, levels = unique(group)))

# Combine data
plot_data <- bind_rows(inner, outer)

# Set legend order: inner (genus) then outer (species), keep order consistent
group_levels <- c(inner$group, outer$group)
plot_data$group <- factor(plot_data$group, levels = (group_levels))

# Percent and label

plot_data <- plot_data %>%
  group_by(ring) %>%
  mutate(
    percent = count / sum(count),
    label = paste0(group, "\n", round(percent * 100,1), "%")
  )

# Plot
library(RColorBrewer)
genus_list <- unique(inner$group)
paired_colors <- brewer.pal(12, "Paired")
light_colors <- paired_colors[seq(1, 12, by = 2)]  # Odd positions: light
dark_colors  <- paired_colors[seq(2, 12, by = 2)]  # Even positions: dark
reordered_colors <- c(dark_colors, light_colors)
base_colors<-c(reordered_colors,"gray","white","darkgray","lightblue")
names(base_colors) <- genus_list

# Assign sub-colors to species
get_shades <- function(base_col, n) {
  (colorRampPalette(c("white", base_col))(n + 1)[-1])  # Dark to light, drop lightest white
}

# Build color mapping table
species_list <- outer %>%
  group_by(selected_g) %>%
  summarise(species = list(group), .groups = "drop")

color_mapping <- c()

for (i in seq_len(nrow(species_list))) {
  genus <- species_list$selected_g[i]
  species <- unlist(species_list$species[i])
  
  n_species <- length(species)
  shades <- get_shades(base_colors[genus], n_species)
  
  # Assign colors: first gets dark color
  names(shades) <- species
  color_mapping <- c(color_mapping, shades)
}

# Add inner ring colors
color_mapping <- c(color_mapping, base_colors)
ggplot(plot_data, aes(x = ring, y = count, fill = group)) +
  geom_bar(stat = "identity", width = 0.9, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.7) +  # Leave extra space for outer ring labels
  theme_void() +
  geom_text(
    data = subset(plot_data, ring == 1),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  geom_text(
    data = subset(plot_data, ring == 2),
    aes(label = label),
    position = position_stack(vjust = 0.5),  # Place on outer ring
    size = 2.8
  ) +
  scale_fill_manual(values = color_mapping) +
  guides(fill = guide_legend(title = "Taxa")) +
  theme(legend.position = "none")


##MGBC/Mgnify (can be adjusted by kraken result )
base_dir <- "~/fsdownload/kraken_result/MGBC_s/"  # Please change to your actual path for MGBC/Mgnify result
# Initialize list to store data
taxonomy_data_list <- list()
report_files <- list.files(base_dir, pattern = "_sc_taxonomy\\.report$", full.names = TRUE)
report_files <- report_files[grepl("412", report_files)]
for (file in report_files) {
  sample_name <- basename((file))  # Get sample folder name
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,quote = "")  # Read data
  data$Sample_Name <- sample_name  # Add sample name column
  taxonomy_data_list[[sample_name]] <- data  # Store in list
}

if(T){
  all_taxonomy_data <- do.call(rbind, taxonomy_data_list)
  rownames(all_taxonomy_data)<-1:nrow(all_taxonomy_data)
  all_taxonomy_data$Sample_Name<-gsub("_sc_taxonomy.report","",all_taxonomy_data$Sample_Name)
}

all_taxonomy_data1<-subset(all_taxonomy_data,fraction_total_reads>0.5 & kraken_assigned_reads>50)

top_names <- all_taxonomy_data1 %>%
  count(name, Sample_Name) %>%
  group_by(name) %>%
  summarise(total_count = sum(n)) %>%
  arrange(desc(total_count))%>%
  slice_head(n = 10) %>%
  pull(name)  # Extract top 10 `name` list

all_taxonomy_data1$selected_g<-sub(" .*", "", all_taxonomy_data1$name)
# Mark which are top10
all_taxonomy_data1 <- all_taxonomy_data1 %>%
  mutate(
    is_top = name %in% top_names
  )

# Extract all genera of top10 (for genus-level matching)
top_genus <- all_taxonomy_data1 %>%
  filter(is_top) %>%
  pull(selected_g) %>%
  unique()

# Generate cleaned_name
all_taxonomy_data1 <- all_taxonomy_data1 %>%
  mutate(
    cleaned_name = case_when(
      name %in% top_names ~ name,
      selected_g %in% top_genus ~ paste0("Other_", selected_g),
      TRUE ~ "Other species"
    )
  )
# Classify names outside top 10 as "Other"
all_taxonomy_data1 <- all_taxonomy_data1 %>%
  mutate(selected_g = ifelse(selected_g %in% top_genus, selected_g, "Other genus"))




# Assume df is your data frame (contains name, selected_g, selected_puri)

# Inner ring: compute total abundance by genus (or count)
inner <- all_taxonomy_data1 %>%
  group_by(selected_g) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count)) %>%
  mutate(
    ring = 1,
    group = selected_g,
    genus_order = row_number()
  )
inner <- inner %>%
  mutate(order = ifelse(group == "Other genus", 999, row_number())) %>%
  arrange(order) %>%
  mutate(group = factor(group, levels = group))

# Outer ring: group by genus order, preserve intra-group order
outer <- all_taxonomy_data1 %>%
  group_by(selected_g, cleaned_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(inner %>% select(selected_g = group, genus_order = order), by = "selected_g") %>%
  mutate(
    ring = 2,
    group = cleaned_name,
    is_other = grepl("^Other_", cleaned_name)
  ) %>%
  arrange(genus_order, is_other, desc(count)) %>%
  mutate(group = factor(group, levels = unique(group)))

# Combine data
plot_data <- bind_rows(inner, outer)

# Set legend order: inner (genus) then outer (species), keep order consistent
group_levels <- c(inner$group, outer$group)
plot_data$group <- factor(plot_data$group, levels = (group_levels))

# Percent and label

plot_data <- plot_data %>%
  group_by(ring) %>%
  mutate(
    percent = count / sum(count),
    label = paste0(group, "\n", round(percent * 100,1), "%")
  )

# Plot
library(RColorBrewer)
genus_list <- unique(inner$group)
paired_colors <- brewer.pal(12, "Paired")
light_colors <- paired_colors[seq(1, 12, by = 2)]  # Odd positions: light
dark_colors  <- paired_colors[seq(2, 12, by = 2)]  # Even positions: dark
reordered_colors <- c(dark_colors, light_colors)
base_colors<-c(reordered_colors,"gray","white","darkgray","lightblue")
names(base_colors) <- genus_list

# Set colors for species
get_shades <- function(color, n) {
  colorRampPalette(c("white", color))(n + 1)[-1]  # Drop the lightest
}
# Build color mapping
species_list <- outer %>%
  group_by(selected_g) %>%
  summarise(species = list(group), .groups = "drop")

color_mapping <- c()

for (i in seq_len(nrow(species_list))) {
  genus <- species_list$selected_g[i]
  species <- unlist(species_list$species[i])
  shades <- get_shades(base_colors[genus], length(species))
  names(shades) <- species
  color_mapping <- c(color_mapping, shades)
}

# Add inner ring color
color_mapping <- c(color_mapping, base_colors)
ggplot(plot_data, aes(x = ring, y = count, fill = group)) +
  geom_bar(stat = "identity", width = 0.9, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.7) +  # Leave extra space for outer ring labels
  theme_void() +
  geom_text(
    data = subset(plot_data, ring == 1),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  geom_text(
    data = subset(plot_data, ring == 2),
    aes(label = label),
    position = position_stack(vjust = 0.5),  # Place on outer ring
    size = 2.8
  ) +
  scale_fill_manual(values = color_mapping) +
  guides(fill = guide_legend(title = "Taxa")) +
  theme(legend.position = "none")

