eggnog<-read.delim2("~/Mouse_gut/M1_eggNOG.tsv",header = T,check.names = F)
count_non_dash <- function(df) {
  sapply(df, function(column) sum(column != "-"))
}
count_rows_non_dash <- function(df) {
  sum(apply(df, 1, function(row) any(row != "-")))
}
# Use the functions
count_non_dash(eggnog[c(7,10,13,19)])
anno_num <- count_rows_non_dash(eggnog[c(7,10,13,19)])

gtfinfo<-read.csv("~/Mouse_gut/M1_genome_all.fix_info.csv",header = T,stringsAsFactors = F)

# Create data frame for plotting
pie_data <- data.frame(
  category = c("Annotated genes", "Unannotated genes"),
  count = c(anno_num, nrow(gtfinfo) - anno_num)
)

# Plot pie chart
pie_data$percentage <- round(pie_data$count / sum(pie_data$count) * 100, 1)

# Customize legend labels
pie_data$label <- paste0(pie_data$category, " ", pie_data$count, " (", pie_data$percentage, "%)")

# Plot pie chart
ggplot(pie_data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1,color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c( "Annotated genes" = "#1F78B4","Unannotated genes"="lightgray")) +
  theme(legend.title = element_blank()) + 
  geom_text(aes(label = paste0(percentage, "%")), color = "black", size = 6, position = position_stack(vjust = 0.5),
) 

pie_data2 <- data.frame(
  category = c("Annotated genes", "Unannotated genes"),
  percentage = c("39.1","60.9")
)

ggplot(pie_data2, aes(x = "", y = percentage, fill = category)) +
  geom_bar(stat = "identity", width = 1,color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c( "Annotated genes" = "#1F78B4","Unannotated genes"="lightgray")) +
  theme(legend.title = element_blank()) + 
  geom_text(aes(label = paste0(percentage, "%")), color = "black", size = 6, position = position_stack(vjust = 0.5),
  ) 

DF = eggnog[c(7,13,10,19)]
DF[DF == "-"] <- NA
database_counts <- colSums(!is.na(DF))

# Create data frame for horizontal bar chart
database_data <- data.frame(
  database = names(database_counts),
  count = database_counts
)

# Sort the data
database_data <- database_data %>%
  arrange(desc(count))
database_data$database<-factor(database_data$database,levels = database_data$database)
# Plot horizontal bar chart
ggplot(database_data, aes(x = reorder(database, count), y = count, fill = database))+
  geom_bar(stat = "identity") +
  coord_flip() +scale_y_reverse() +ylim(0, max(database_data$count) * 1.2) +
  labs(title = "Annotated Genes per Database", x = "Database", y = "Count") +
  theme_minimal() + theme(axis.title.x = element_blank(), 
                          axis.text.x = element_blank(), 
                          axis.ticks.x = element_blank(),
                          legend.position = "none") +
  geom_text(aes(label = count), hjust = -0.1, size = 3)+
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = paletteer_d("ggthemes::Classic_20",n=4)) 
  
DF = eggnog[c(7,13,10,19)]
# Count annotated genes across combinations
DF[DF == "-"] <- NA
DF<-DF %>%
  mutate(across(everything(), ~ !is.na(.)))

combinations <- DF %>%
  summarise(
    DB1_DB2 = sum(DF[,1] & DF[,2], na.rm = TRUE),
    DB1_DB3 = sum(DF[,1] & DF[,3], na.rm = TRUE),
    DB1_DB4 = sum(DF[,1] & DF[,4], na.rm = TRUE),
    DB2_DB3 = sum(DF[,2] & DF[,3], na.rm = TRUE),
    DB2_DB4 = sum(DF[,2] & DF[,4], na.rm = TRUE),
    DB3_DB4 = sum(DF[,3] & DF[,4], na.rm = TRUE),
    DB1_DB2_DB3 = sum(DF[,1] & DF[,2] & DF[,3], na.rm = TRUE),
    DB1_DB2_DB4 = sum(DF[,1] & DF[,2] & DF[,4], na.rm = TRUE),
    DB1_DB3_DB4 = sum(DF[,1] & DF[,3] & DF[,4], na.rm = TRUE),
    DB2_DB3_DB4 = sum(DF[,2] & DF[,3] & DF[,4], na.rm = TRUE),
    DB1_DB2_DB3_DB4 = sum(DF[,1] & DF[,2] & DF[,3] & DF[,4], na.rm = TRUE)
  ) 
bubble_data <- data.frame(
  DB1_DB2 = c(T,T,F,F),
  DB1_DB3 = c(T,F,T,F),
  DB1_DB4 = c(T,F,F,T),
  DB2_DB3 = c(F,T,T,F),
  DB2_DB4 = c(F,T,F,T),
  DB3_DB4 = c(F,F,T,T),
  DB1_DB2_DB3 = c(T,T,T,F),
  DB1_DB2_DB4 = c(T,T,F,T),
  DB1_DB3_DB4 = c(T,F,T,T),
  DB2_DB3_DB4 = c(F,T,T,T),
  DB1_DB2_DB3_DB4 = c(T,T,T,T)
)
library(reshape2)
bubble_data$database <- colnames(DF)

ml<-melt(bubble_data,id="database")
ml$database<-factor(ml$database,levels = rev(bubble_data$database))
p2<-ggplot(ml, aes(x = variable, y = database, size = value, color = value)) +
  geom_point() +
  scale_size_manual(values = c(`FALSE` = 0, `TRUE` = 5)) +
  scale_color_manual(values = c(`FALSE` = "white", `TRUE` = "red")) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),axis.title = element_blank(), legend.position = "none")
p2
# Create data frame for vertical bar chart
library(tidyr)
combinations<-t(combinations)%>%as.data.frame()
combinations$combination<-rownames(combinations)
combinations$V1<-log10(combinations$V1)
combinations$combination<-factor(combinations$combination,levels = combinations$combination)
# Plot vertical bar chart
p1<-ggplot(combinations, aes(x = combination, y = V1)) +
  geom_bar(stat = "identity") +
  labs(title = "Gene Count by Annotation Combination", x = "", y = "log10 Count") +
  theme_minimal()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), legend.position = "none")
library(cowplot)
plot_grid(p1, p2, ncol = 1, align = "v",rel_heights = c(2,1))


