
# Libraries ---------------------------------------------------------------
library(ggplot2)
library(reshape)
library(dplyr)
# Data --------------------------------------------------------------------
heatmap_data <- read.csv(here::here("data/multivariate_results/clustering_results_dataframe.csv"), row.names = 1)


# Preparing Data ----------------------------------------------------------

## For continuous data
heatmap_data_scaled <- heatmap_data
heatmap_data_scaled$CX3 <- as.numeric(scale(heatmap_data_scaled$CX3))
heatmap_data_scaled$CX4 <- as.numeric(scale(heatmap_data_scaled$CX4))
heatmap_data_scaled$CX10 <- as.numeric(scale(heatmap_data_scaled$CX10))
heatmap_data_scaled$CX14 <- as.numeric(scale(heatmap_data_scaled$CX14))
heatmap_data_scaled$age <- as.numeric(scale(heatmap_data_scaled$age))

heatmap_data_scaled_cn <- as.matrix(heatmap_data_scaled[, c('CX3', 'CX4', 'CX10', 'CX14', 'age')])
heatmap_data_scaled_cn_melt <- melt(heatmap_data_scaled_cn)
heatmap_data_scaled_cn_melt$X2 <- factor(heatmap_data_scaled_cn_melt$X2, levels = c('age', 'CX3', 'CX4', 'CX10', 'CX14'))

## For categorical data
heatmap_data_categorical_gene <- subset(heatmap_data_scaled, 
                                          select = -c(CX3, CX4, CX10, CX14,
                                                      kproto_clustering, kmedoids_clustering,
                                                      hier_clustering,
                                                      final_cluster_labels, age))
# heatmap_data_categorical_wgd <- subset(heatmap_data_scaled, 
#                                                select = c(WGD))
# 
# heatmap_data_categorical_condi <- subset(heatmap_data_scaled, 
#                                        select = c(Condition))


heatmap_data_categorical_all <- melt(as.matrix(heatmap_data_categorical_gene))
heatmap_data_categorical_all$value <- factor(heatmap_data_categorical_all$value, 
                                             levels = c(' 0', ' 1', ' 2', '-1', '-2','TRUE', 'FALSE', 'Resistant', 'Sensitive'))
heatmap_data_categorical_all$X2 <- factor(heatmap_data_categorical_all$X2,
                                          levels = c('Condition', 'WGD', 'ATM', 'BARD1', 'BRCA2',
                                                     'BRIP1', 'CHEK1', 'CHEK2', 'FAM175A',
                                                     'MRE11A', 'NBN', 'PALB2', 'RAD51C'))


# heatmap_data_scaled_categorical_melt_gene <- melt(as.matrix(heatmap_data_categorical_gene))
# heatmap_wgd_melt <- melt(as.matrix(heatmap_data_categorical_wgd))
# heatmap_condit_melt <- melt(as.matrix(heatmap_data_categorical_condi))

# combined_categorical_melt <- merge(heatmap_data_scaled_categorical_melt_gene,
#                                    heatmap_wgd_melt, by.x = 1, by.y = 1)
# 
# combined_categorical_melt2 <- merge(combined_categorical_melt,
#                                     heatmap_condit_melt, by.x = 1, by.y = 1)


# For clustering results



# Plotting heatmap for continuous  -----------------------------------------
heatmap_data_scaled_cn_melt %>%
  ggplot(aes(x = X1, y = X2, fill = value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

heatmap(heatmap_data_scaled_cn, scale = "none")


# Plotting heatmap for categorical ----------------------------------------

heatmap_data_categorical_all %>%
  ggplot(aes(x = X1, y = X2, fill = value)) + 
  geom_tile(width = 0.9, height = 0.9) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plotting heatmap for clustering results ---------------------------------


