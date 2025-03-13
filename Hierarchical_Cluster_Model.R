####LOAD PACKAGES####
library(readr)
library(dplyr)
library(tidyverse)
library(klaR)
library(scatterplot3d)
library(reshape2)
library(ggplot2)
library(dendextend)
library(circlize)

####SET SEED####
set.seed(7895)

####LOAD THE DATA####

# Ensure the analysed data is tidy, with categorical factor variables. 
# Presence of a comorbidity should be encoded as '1', and absence as '0'.
# Split the data into 'case' and 'control' cohorts, based on the outcome variable.


####CASE HCA ANALYSIS####

# Cluster the patients.
distances_hospitalised <- dist(df_hospitalised, method = "euclidean")
clusterHosp <- hclust(distances_hospitalised, method = "ward.D2")
clusterGroups = cutree(clusterHosp, k = 4) 

# Initialise an empty data frame to store the results.
final_data_hospitalised <- data.frame()
for (col in names(df_hospitalised)) {
  col_means <- tapply(df_hospitalised[[col]], clusterGroups, mean)
  col_means_df <- data.frame(cluster = names(col_means), mean = col_means)
  colnames(col_means_df) <- c("cluster", col)
  if (nrow(final_data_hospitalised) == 0) {
    final_data_hospitalised <- col_means_df
  } else {
    final_data_hospitalised <- merge(final_data_hospitalised, col_means_df, by = "cluster")
  }
}

# Cluster comorbidities. 
df_circ_dendrogram <- df_hospitalised
df_circ_dendrogram_transposed <- t(df_circ_dendrogram)
distances_transposed <- dist(df_circ_dendrogram_transposed, method = "euclidean")
clusterHosp_transposed <- hclust(distances_transposed, method = "ward.D2")
dendrogram_sampled_transposed <- as.dendrogram(clusterHosp_transposed)

# Circular dendrogram.
circ_dendro_colors = c("olivedrab", "indianred", "goldenrod", "steelblue")
k = 4
dend_sampled_transposed <- dend_sampled_transposed %>% 
  color_branches(k=k, col = circ_dendro_colors) %>%
  set("branches_lwd", 2) %>%
  set("branches_lty", 1) %>%
  color_labels(k=k, col = circ_dendro_colors)
circlize_dendrogram(
  dend_sampled_transposed,
  facing = 'outside',
  labels = TRUE,
  labels_track_height = 0.410,
  labels_cex = 0.5,
  labels_just = c(1, 0.5),
  dend_track_height = 0.54,
  gap.degree = 1)

#==================================================================================#

####CONTROL HCA ANALYSIS####

# Cluster the patients.
distances_not_hospitalised <- dist(df_not_hospitalised, method = "euclidean")
clusterNotHosp <- hclust(distances_not_hospitalised, method = "ward.D2")
clusterGroups = cutree(clusterNotHosp, k = 5)

# Initialise an empty data frame to store the results.
final_data_not_hospitalised <- data.frame()
for (col in names(df_not_hospitalised)) {
  col_means <- tapply(df_not_hospitalised[[col]], clusterGroups, mean)
  col_means_df <- data.frame(cluster = names(col_means), mean = col_means)
  colnames(col_means_df) <- c("cluster", col)
  if (nrow(final_data_not_hospitalised) == 0) {
    final_data_not_hospitalised <- col_means_df
  } else {
    final_data_not_hospitalised <- merge(final_data_not_hospitalised, col_means_df, by = "cluster")
  }
}