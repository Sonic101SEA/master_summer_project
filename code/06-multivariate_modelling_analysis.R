
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
library(caret)
library(dendextend)
library(Rtsne)
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

draw_confusion_matrix <- function(cm, title_cm = "") {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(title_cm, cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'Resistant', cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'Sensitive', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'Resistant', cex=1.2, srt=90)
  text(140, 335, 'Sensitive', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(50, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(50, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
}  
# Multivariate Analysis - Conducting Unsupervised Clustering --------------
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
# variables_with_clustering$dbscan_clustering <- as.numeric(dbscan_results$cluster)


# Checking which cluster is sensitive or resistant
kproto_glm <- 
  glm(kproto_clustering ~ WGD, variables_with_clustering, 
      family = binomial(link = "logit"))
summary(kproto_glm)
  ## Cluster 2 is resistance, cluster 1 is sensitivity

kmedoids_glm <-
  glm(kmedoids_clustering ~ WGD, variables_with_clustering, 
      family = binomial(link = "logit"))
summary(kmedoids_glm)
  ## Cluster 2 is resistance, cluster 1 is sensitivity

hier_glm <-
  glm(hier_clustering ~ WGD, variables_with_clustering, 
      family = binomial(link = "logit"))
summary(hier_glm)
  ## Cluster 2 is resistance, cluster 1 is sensitivity

# dbscan_glm <-
#   glm(dbscan_clustering ~ CX3, variables_with_clustering,
#       family = binomial(link = "logit"))
# summary(dbscan_glm)
#   ## Cluster 1 is resistant, cluster 0 is sensitivity

# ## Discriminant analysis
# ### First looking at how the grouping compares between clusters
# table(hier = variables_with_clustering$hier_clustering, kmedoids = variables_with_clustering$kmedoids_clustering)
# 
# lda_clustering <- lda(Condition ~ hier_clustering + kmedoids_clustering, variables_with_clustering)
# table(variables_with_clustering$Condition, predict(lda_clustering)$class)


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
                                                                 "Resistant", "Sensitive")

# Checkpoint saving dataframe after clustering
# variables_with_clustering_labelled_ids <- variables_with_clustering_labelled
# rownames(variables_with_clustering_labelled_ids) <- labels$tumour_specimen_aliquot_id
# write.csv(variables_with_clustering_labelled_ids, here::here("data/multivariate_results/clustering_results_dataframe.csv"))

# ## For DBSCAN
# variables_with_clustering_labelled$dbscan_clustering <- ifelse(variables_with_clustering$hier_clustering == 1,
#                                                              "Resistant", "Sensitive")

# Computing final cluster labels from 3 algorithms
variables_with_clustering_labelled$final_cluster_labels <- NA

variables_with_clustering_labelled$final_cluster_labels <-
  apply(variables_with_clustering_labelled[c('kproto_clustering',
                                            'kmedoids_clustering',
                                            'hier_clustering')],
        1, chooseBestModel)


# Multivariate analysis - Visualisation of Results -----------------------
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

# ## DBSCAN clustering vs ground truth
# sum(variables_with_clustering_labelled$dbscan_clustering == variables_with_clustering_labelled$Condition) /
#   nrow(variables_with_clustering_labelled) * 100

## Final cluster after ensemble of 3 algorithms
sum(variables_with_clustering_labelled$final_cluster_labels == variables_with_clustering_labelled$Condition) /
  nrow(variables_with_clustering_labelled) * 100

## Producing confusion matrix based on the prediction compared to ground truth
### Normal table to see values
table(Pred = variables_with_clustering_labelled$final_cluster_labels,
      Actual = variables_with_clustering_labelled$Condition)

### More detailed results from confusion matrix (includes sensitivity etc.)
conf_matrix_kproto <- confusionMatrix(factor(variables_with_clustering_labelled$kproto_clustering),
                                      reference = variables_with_clustering_labelled$Condition,
                                      mode = "everything")

conf_matrix_kmedoids <- confusionMatrix(factor(variables_with_clustering_labelled$kmedoids_clustering),
                                        reference = variables_with_clustering_labelled$Condition,
                                        mode = "everything")

conf_matrix_hier <- confusionMatrix(factor(variables_with_clustering_labelled$hier_clustering),
                                    reference = variables_with_clustering_labelled$Condition,
                                    mode = "everything")

conf_matrix_ensemble <- confusionMatrix(factor(variables_with_clustering_labelled$final_cluster_labels),
                                          reference = variables_with_clustering_labelled$Condition,
                                          mode = "everything")

### Plotting confusion matrix
#pdf(here::here("graphs/analysis/confusion_matrix_ensemble.pdf"))
draw_confusion_matrix(conf_matrix_ensemble, "Ensemble Confusion Matrix")
#dev.off()

#pdf(here::here("graphs/analysis/confusion_matrix_kproto.pdf"))
draw_confusion_matrix(conf_matrix_kproto, "K-prototypes Clustering Confusion Matrix")
#dev.off()

#pdf(here::here("graphs/analysis/confusion_matrix_kmedoids.pdf"))
draw_confusion_matrix(conf_matrix_kmedoids, "K-Medoids Clustering with Gower's Distance Confusion Matrix")
#dev.off()

#pdf(here::here("graphs/analysis/confusion_matrix_hier.pdf"))
draw_confusion_matrix(conf_matrix_hier, "Hierarchical Clustering with Gower's Distance Confusion Matrix")
#dev.off()


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

# dbscan_clustering_labelled <-
#   variables_with_clustering %>%
#     ggplot(aes(y = dbscan_clustering)) +
#     geom_bar(aes(fill = Condition))

## Visualisation of clusters
### Kmedoids results
tsne_obj <- Rtsne(gowers_distances, is_distance = TRUE, perplexity = 12) # Converting to T-SNE object

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(kmedoids_results$clustering),
         name = labels$Condition)
  
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))

clusplot(kmedoids_results, color = TRUE, line = 0)


### Hierarchical clustering results
dendrogram_object <- as.dendrogram(hcluster_results) # Creating dendrogram object

# nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
#                 cex = 0.7, col = "blue")
labels(dendrogram_object) <- labels$Condition[order.dendrogram(dendrogram_object)]
# dendrogram_object <- set(dendrogram_object, "labels_cex", 0.2)
dendrogram_object <- color_branches(dendrogram_object, k = 2, col = as.numeric(unique(labels$Condition)))
dendrogram_object <- dendrogram_object %>%
  set("labels_colors", as.numeric(labels$Condition),
                         order_value = TRUE) %>%
  set("labels_cex", 0.7)

#pdf(here::here("graphs/analysis/hierarchical_clustering_dendrogram_cluster_coloured.pdf"))
plot(dendrogram_object, 
     main = "Hierarchical Clustering", ylab = "Height")
#dev.off()

#pdf(here::here("graphs/analysis/hierarchical_clustering_dendrogram_labelled.pdf"))
plot(hcluster_results, labels = labels$tumour_specimen_aliquot_id, cex = 0.3,
     main = "Hierarchical Clustering", xlab = "Height", hang = -1)
#dev.off()


