
# Libraries ---------------------------------------------------------------
library(ggplot2)
library(reshape)
library(dplyr)
library(egg)
# Data --------------------------------------------------------------------
heatmap_data <- read.csv(here::here("data/multivariate_results/clustering_results_dataframe.csv"), row.names = 1)


# Preparing Data ----------------------------------------------------------
# Ordering data
## Sort by therapy outcome
levels_sample_Condition <- subset(heatmap_data, select = c(Condition))
levels_sample_Condition$id <- rownames(levels_sample_Condition)
rownames(levels_sample_Condition) <- NULL
sorted_levels_sample_Condition <- levels_sample_Condition[order(levels_sample_Condition$Condition), ]
sorted_levels_sample_Condition$id <- factor(sorted_levels_sample_Condition$id)

## Sort by CX14
levels_sample_cx14 <- subset(heatmap_data, select = c(CX14))
levels_sample_cx14$id <- rownames(levels_sample_cx14)
rownames(levels_sample_cx14) <- NULL
sorted_levels_sample_cx14<- levels_sample_cx14[order(levels_sample_cx14$CX14), ]
sorted_levels_sample_cx14$id <- factor(sorted_levels_sample_cx14$id)

## Sort by therapy outcome and CX14
levels_sample_Condition_cx14 <- subset(heatmap_data, select = c(Condition, CX14))
levels_sample_Condition_cx14$id <- rownames(levels_sample_Condition_cx14)
rownames(levels_sample_Condition_cx14) <- NULL

### Arranging by Condition and increasing values of CX14
sorted_levels_cx14_condition <-
  levels_sample_Condition_cx14 %>%
  arrange(Condition, CX14)

sorted_levels_cx14_condition$id <- factor(sorted_levels_cx14_condition$id)

## Order of data we want
# sample_order <- as.vector(sorted_levels_sample_Condition$id)
# sample_order <- as.vector(sorted_levels_sample_cx14$id)
sample_order <- as.vector(sorted_levels_cx14_condition$id)


## For continuous data
heatmap_data_scaled <- heatmap_data
heatmap_data_scaled$CX3 <- as.numeric(scale(heatmap_data_scaled$CX3))
heatmap_data_scaled$CX4 <- as.numeric(scale(heatmap_data_scaled$CX4))
heatmap_data_scaled$CX10 <- as.numeric(scale(heatmap_data_scaled$CX10))
heatmap_data_scaled$CX14 <- as.numeric(scale(heatmap_data_scaled$CX14))
heatmap_data_scaled$age <- as.numeric(scale(heatmap_data_scaled$age))
heatmap_data_scaled$WGD <- ifelse(heatmap_data_scaled$WGD == TRUE, "Yes", "No")

heatmap_data_scaled_cn <- as.matrix(heatmap_data_scaled[, c('CX3', 'CX4', 'CX10', 'CX14', 'age')])
heatmap_data_scaled_cn_melt <- melt(heatmap_data_scaled_cn)
heatmap_data_scaled_cn_melt$X2 <- factor(heatmap_data_scaled_cn_melt$X2, levels = c('age', 'CX3', 'CX4', 'CX10', 'CX14'))

### Ordering data
heatmap_data_scaled_cn_melt$X1 <- factor(heatmap_data_scaled_cn_melt$X1, levels = sample_order)

## For categorical data
heatmap_data_categorical_gene <- subset(heatmap_data_scaled, 
                                          select = -c(CX3, CX4, CX10, CX14,
                                                      kproto_clustering, kmedoids_clustering,
                                                      hier_clustering,
                                                      final_cluster_labels, Condition, age))

# heatmap_data_categorical_wgd <- subset(heatmap_data_scaled, 
#                                                select = c(WGD))
# 
# heatmap_data_categorical_condi <- subset(heatmap_data_scaled, 
#                                        select = c(Condition))


heatmap_data_categorical_all <- melt(as.matrix(heatmap_data_categorical_gene))
heatmap_data_categorical_all$value <- factor(heatmap_data_categorical_all$value, 
                                             levels = c(' 0', ' 1', ' 2', '-1', '-2', 'Yes', 'No'))
heatmap_data_categorical_all$X2 <- factor(heatmap_data_categorical_all$X2,
                                          levels = c('Condition', 'WGD', 'ATM', 'BARD1', 'BRCA2',
                                                     'BRIP1', 'CHEK1', 'CHEK2', 'FAM175A',
                                                     'MRE11A', 'NBN', 'PALB2', 'RAD51C'))

### Ordering data

heatmap_data_categorical_all$X1 <- factor(heatmap_data_categorical_all$X1, levels = sample_order)

# heatmap_data_scaled_categorical_melt_gene <- melt(as.matrix(heatmap_data_categorical_gene))
# heatmap_wgd_melt <- melt(as.matrix(heatmap_data_categorical_wgd))
# heatmap_condit_melt <- melt(as.matrix(heatmap_data_categorical_condi))

# combined_categorical_melt <- merge(heatmap_data_scaled_categorical_melt_gene,
#                                    heatmap_wgd_melt, by.x = 1, by.y = 1)
# 
# combined_categorical_melt2 <- merge(combined_categorical_melt,
#                                     heatmap_condit_melt, by.x = 1, by.y = 1)


# For clustering results
heatmap_data_clustering <- subset(heatmap_data_scaled, 
                           select = c( kproto_clustering, kmedoids_clustering,
                                       hier_clustering,
                                       final_cluster_labels, Condition))
heatmap_data_clustering_melt <- melt(as.matrix(heatmap_data_clustering))

### Setting sample order
heatmap_data_clustering_melt$X1 <- factor(heatmap_data_clustering_melt$X1, levels = sample_order)

### Setting clustering algorithms order
heatmap_data_clustering_melt$X2 <- factor(heatmap_data_clustering_melt$X2,
                                          levels = c("Condition", "final_cluster_labels",
                                                     "hier_clustering", "kmedoids_clustering",
                                                     "kproto_clustering"))


# Plotting heatmap for continuous  -----------------------------------------
mid <- 0

continuous_heatmap <-
heatmap_data_scaled_cn_melt %>%
  ggplot(aes(x = X1, y = X2, fill = value)) + 
  geom_tile(width = 0.9, height = 0.9) +
  # theme(axis.text.x = element_blank(), 
  #       axis.ticks.x = element_blank()) +
  scale_fill_gradient2(midpoint = mid, low = "blue", mid = "white", high = "red") +
  scale_x_discrete(labels = c(1:38)) +
  xlab('Patients') +
  ylab('Continuous predictors') + labs(fill = "Scaled level")

# heatmap(heatmap_data_scaled_cn, scale = "none")


# Plotting heatmap for categorical ----------------------------------------

categorical_heatmap <-
heatmap_data_categorical_all %>%
  ggplot(aes(x = X1, y = X2, fill = value)) + 
  geom_tile(width = 0.9, height = 0.9) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  ylab('Categorical predictors') + labs(fill = "Attributes")

# Plotting heatmap for clustering results ---------------------------------
clustering_results_heatmap <- 
heatmap_data_clustering_melt %>%
  ggplot(aes(x = X1, y = X2, fill = value)) +
  geom_tile(width = 0.9, height = 0.9) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  ylab('Clustering Algorithm') + labs(fill = "Outcome") +
  scale_y_discrete(labels = c("Therapy", "Ensemble", "Hierarchical", "K-Medoids", "K-prototype")) 

# Plotting combined heatmap -----------------------------------------------

combined_heatmap <-
ggarrange(clustering_results_heatmap, categorical_heatmap, continuous_heatmap)

# ggsave(here::here("graphs/analysis/heatmap_predictors_orderedBy_cx14_condition_with_patient_numbering.pdf"), combined_heatmap, height = 9, width = 16)

       