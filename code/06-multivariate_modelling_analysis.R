
# Libraries ---------------------------------------------------------------

library(glmnet)
library(klaR)
library(clustMixType)
library(cluster)
# Data --------------------------------------------------------------------

modelling_data <- read.csv(here::here("data/final_dataframe.csv"), row.names = 1)

## Remove patients with missing data
modelling_data_removed_na<- modelling_data[complete.cases(modelling_data), ]

## Removing variables that have low variability: CX6, 8, 13, 17 and mhBRCA1, mhCHEK2, mhPALB2
columns_to_remove <- c('CX6', 'CX8', 'CX13', 'CX17', 'mhBRCA1', 'mhCHEK2', 'mhPALB2')
final_modelling_data <- modelling_data_removed_na[, !(names(modelling_data_removed_na) %in% columns_to_remove)]

## Factorising columns
columns_to_factorise <- c('BARD1', 'FAM175A', 'NBN', 'MRE11A', 'ATM', 'CHEK1', 'BRCA2', 'PALB2', 'RAD51D', 'BRCA1', 'RAD51C', 'BRIP1', 'CHEK2')
final_modelling_data[columns_to_factorise] <- lapply(final_modelling_data[columns_to_factorise], factor, levels = c("0", "-2", "-1", "1", "2"))
final_modelling_data$mhBRCA2 <- factor(final_modelling_data$mhBRCA2)
final_modelling_data$Condition <- ifelse(final_modelling_data$Condition == "Resistant", 1, 0) # Resistant = 1, Sensitive = 0
# final_modelling_data$WGD <- ifelse(final_modelling_data$WGD == "TRUE", 1, 0) # TRUE = 1, FALSE = 0

## Extracting labels
labels <- subset(final_modelling_data, select = c('tumour_specimen_aliquot_id', 
                                                  'Condition'))

## Setting row names then removing id columns
# rownames(final_modelling_data) <- final_modelling_data$tumour_specimen_aliquot_id
# final_modelling_data <- subset(final_modelling_data, 
#                                select = -c(TCGA_id, tumour_specimen_aliquot_id))


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

# Clustering

## kmodes clustering - only for categorical variables

kmodes_results <- kmodes(categorical_modelling_data, 2, 5)

## K-prototype clustering
kproto_results <- kproto(mixed_modelling_data, k = 2, iter.max = 10)

## Gowers distance then use k-medoids clustering
gowers_distances <- daisy(mixed_modelling_data, metric = "gower", 
                          type = list(asymm = c(8))) # Setting WGD as asymmetric binary

kmedoids_results <- pam(gowers_distances, k = 2)

# Lasso
lasso_model <- cv.glmnet(x = data.matrix(final_modelling_data[4:ncol(final_modelling_data)]), y = final_modelling_data$Condition, family = "binomial", alpha = 1, nfolds = 10)
