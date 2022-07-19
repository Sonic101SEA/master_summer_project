
# Libraries ---------------------------------------------------------------

library(glmnet)
library(klaR)
library(clustMixType)
library(cluster)
library(mclust)
library(factoextra)
library(ggplot2)
library(dplyr)
library(dbscan)
# Data --------------------------------------------------------------------

modelling_data <- read.csv(here::here("data/final_dataframe.csv"), row.names = 1)

## Remove patients with missing data
modelling_data_removed_na<- modelling_data[complete.cases(modelling_data), ]

## Removing variables that have low variability: CX6, 8, 13, 17 and mhBRCA1, mhCHEK2, mhPALB2
columns_to_remove <- c('CX6', 'CX8', 'CX13', 'CX17', 
                       'mhBRCA1', 'mhCHEK2', 'mhPALB2')
final_modelling_data <- modelling_data_removed_na[, !(names(modelling_data_removed_na) %in% columns_to_remove)]

## Factorising columns
columns_to_factorise <- c('BARD1', 'FAM175A', 'NBN', 'MRE11A', 'ATM', 'CHEK1', 
                          'BRCA2', 'PALB2', 'RAD51D', 'BRCA1', 
                          'RAD51C', 'BRIP1', 'CHEK2')
final_modelling_data[columns_to_factorise] <- lapply(final_modelling_data[columns_to_factorise], 
                                                     factor, levels = c("0", "-2", "-1", "1", "2"))
final_modelling_data$mhBRCA2 <- factor(final_modelling_data$mhBRCA2)
# final_modelling_data$Condition <- ifelse(final_modelling_data$Condition == "Resistant", 1, 0) # Resistant = 1, Sensitive = 0
# final_modelling_data$WGD <- ifelse(final_modelling_data$WGD == "TRUE", 1, 0) # TRUE = 1, FALSE = 0

## Extracting labels
labels <- subset(final_modelling_data, select = c('tumour_specimen_aliquot_id', 
                                                  'Condition'))

## Setting row names then removing id columns
# rownames(final_modelling_data) <- final_modelling_data$tumour_specimen_aliquot_id
# final_modelling_data <- subset(final_modelling_data, 
#                                select = -c(TCGA_id, tumour_specimen_aliquot_id))

# Functions ---------------------------------------------------------------

chooseBestModel <- function(x) {
  tabulatedOutcomes <- table(x)
  sortedOutcomes <- sort(tabulatedOutcomes, decreasing=TRUE)
  mostCommonLabel <- names(sortedOutcomes)[1]
  mostCommonLabel
}

# Multivariate Analysis - Unsupervised Clustering -------------------------
set.seed(1)
# Data Preparation
## Selecting variables that have more than 0.25 p-value
select_categorical <- c('ATM', 'NBN', 'CHEK1')

select_mixed <- c('ATM', 'NBN', 'CHEK1', 'age', 'FAM175A', 'CX14',
                  'PALB2', 'WGD', 'MRE11A', 'CX4', 'CX10', 'BRIP1',
                  'CHEK2', 'CX3', 'RAD51C', 'BARD1', 'BRCA2')

categorical_modelling_data <- subset(final_modelling_data, # Subset categorical data
                                select =  select_categorical)

mixed_modelling_data <- subset(final_modelling_data,
                               select = select_mixed)

## Scaling numeric variables
# mixed_modelling_data$age <- as.numeric(scale(mixed_modelling_data$age))
# mixed_modelling_data$CX14 <- as.numeric(scale(mixed_modelling_data$CX14))
# mixed_modelling_data$CX4 <- as.numeric(scale(mixed_modelling_data$CX4))
# mixed_modelling_data$CX3 <- as.numeric(scale(mixed_modelling_data$CX3))
# mixed_modelling_data$CX10 <- as.numeric(scale(mixed_modelling_data$CX10))

# Clustering

## K-prototype clustering
kproto_results <- kproto(mixed_modelling_data, k = 2, iter.max = 10)

## Gowers distance then use k-medoids clustering
gowers_distances <- daisy(mixed_modelling_data, metric = "gower", 
                          type = list(asymm = c(8))) # Setting WGD as asymmetric binary

kmedoids_results <- pam(gowers_distances, k = 2)

## Hierarchical clustering using gowers distance
hcluster_results <- hclust(gowers_distances, method = "complete")
cut_2_clusters <- cutree(hcluster_results, k = 2) # Cutting tree to obtain 2 clusters

# ## DBSCAN Clustering
# dbscan_results <- dbscan(gowers_distances, eps = 0.3)


## kmodes clustering - only for categorical variables
# kmodes_results <- kmodes(categorical_modelling_data, 2, 5)

## Checking the clusters results with adjusted Rand Index
adjustedRandIndex(kproto_results$cluster, cut_2_clusters)

# Combining clustering results into a dataframe to see their cluster labels from each algorithm
# kproto_results$cluster
# kmedoids_results$clustering
# cut_2_clusters

variables_with_clustering <- mixed_modelling_data
variables_with_clustering$kproto_clustering <- factor(kproto_results$cluster)
variables_with_clustering$kmedoids_clustering <- factor(kmedoids_results$clustering)
variables_with_clustering$hier_clustering <- factor(cut_2_clusters)
variables_with_clustering$Condition <- final_modelling_data$Condition
# variables_with_clustering$dbscan_clustering <- factor(dbscan_results$cluster)


# Checking which cluster is sensitive or resistant
kproto_glm <- 
  glm(kproto_clustering ~ CX3, variables_with_clustering, 
      family = binomial(link = "logit"))
summary(kproto_glm)
  ## Cluster 2 is resistance, cluster 1 is sensitivity

kmedoids_glm <-
  glm(kmedoids_clustering ~ CX3, variables_with_clustering, 
      family = binomial(link = "logit"))
summary(kmedoids_glm)
  ## Cluster 2 is resistance, cluster 1 is sensitivity

hier_glm <-
  glm(hier_clustering ~ CX3, variables_with_clustering, 
      family = binomial(link = "logit"))
summary(hier_glm)
  ## Cluster 2 is sensitivty, cluster 1 is resistance

# dbscan_glm <-
#   glm(dbscan_clustering ~ CX3, variables_with_clustering, 
#       family = binomial(link = "logit"))
# summary(dbscan_glm)
#   ## Cluster 1 is resistance, cluster 0 is sensitivity

# Changing the labels for the cluster names
variables_with_clustering_labelled <- variables_with_clustering

## For kproto
variables_with_clustering_labelled$kproto_clustering <- ifelse(variables_with_clustering$kproto_clustering == 2,
                                                      "Resistant", "Sensitive")

## For kmedoids
variables_with_clustering_labelled$kmedoids_clustering <- ifelse(variables_with_clustering$kmedoids_clustering == 2,
                                                               "Resistant", "Sensitive")

## For hierarchical clustering
variables_with_clustering_labelled$hier_clustering <- ifelse(variables_with_clustering$hier_clustering == 2,
                                                                 "Sensitive", "Resistant")

# Computing final cluster labels from 3 algorithms
variables_with_clustering_labelled$final_cluster_labels <- NA

variables_with_clustering_labelled$final_cluster_labels <-
  apply(variables_with_clustering_labelled[c('kproto_clustering',
                                            'kmedoids_clustering',
                                            'hier_clustering')],
        1, chooseBestModel)

# Computing accuracy of cluster results compared to ground truth labels

## Kprototype vs ground truth
sum(variables_with_clustering_labelled$kproto_clustering == variables_with_clustering_labelled$Condition) /
  nrow(variables_with_clustering_labelled) * 100

## Kmedoids vs ground truth
sum(variables_with_clustering_labelled$kmedoids_clustering== variables_with_clustering_labelled$Condition) /
  nrow(variables_with_clustering_labelled) * 100

## Hierarchical clustering vs ground truth
sum(variables_with_clustering_labelled$hier_clustering == variables_with_clustering_labelled$Condition) /
  nrow(variables_with_clustering_labelled) * 100

## Final cluster after ensemble of 3 algorithms
sum(variables_with_clustering_labelled$final_cluster_labels == variables_with_clustering_labelled$Condition) /
  nrow(variables_with_clustering_labelled) * 100

# Plotting the clusters
## Barplots
kproto_clustering_labelled <-
  variables_with_clustering %>%
  ggplot(aes(y = kproto_clustering)) +
  geom_bar(aes(fill = Condition))

kmedoids_clustering_labelled <-
  variables_with_clustering %>%
  ggplot(aes(y = kmedoids_clustering)) +
  geom_bar(aes(fill = Condition))

hier_clustering_labelled <-
  variables_with_clustering %>%
    ggplot(aes(y = cut_2_clusters)) +
    geom_bar(aes(fill = Condition))

## 2D plot
### Kmedoids results
clusplot(kmedoids_results, color = TRUE,
         shade = TRUE, labels = 2, line = 0)

fviz_cluster(kmedoids_results, kmedoids_results$clustering)

# Lasso
# lasso_model <- cv.glmnet(x = data.matrix(final_modelling_data[4:ncol(final_modelling_data)]), y = final_modelling_data$Condition, family = "binomial", alpha = 1, nfolds = 10)
