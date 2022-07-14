
# Libraries ---------------------------------------------------------------

library(glmnet)

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
final_modelling_data$WGD <- ifelse(final_modelling_data$WGD == "TRUE", 1, 0) # TRUE = 1, FALSE = 0

# Multivariate Analysis ---------------------------------------------------

lasso_model <- cv.glmnet(x = final_modelling_data[4:17], y = final_modelling_data$Condition, family = "binomial", alpha = 1, nfolds = 10)
