
# Libraries ---------------------------------------------------------------
library(purrr)
library(broom)
# Data --------------------------------------------------------------------
modelling_data <- read.csv(here::here("data/final_dataframe.csv"), 
                           row.names = 1)

## Remove patients with missing data
modelling_data_removed_na<- modelling_data[complete.cases(modelling_data), ]

## Removing variables that have low variability: CX6, 8, 13, 17 and mhBRCA1, mhCHEK2, mhPALB2
columns_to_remove <- c('CX6', 'CX8', 'CX13', 'CX17', 
                       'mhBRCA1', 'mhCHEK2', 'mhPALB2')
final_modelling_data <- modelling_data_removed_na[, !(names(modelling_data_removed_na) %in% columns_to_remove)]

## Factorising columns
columns_to_factorise <- c('BARD1', 'FAM175A', 'NBN', 'MRE11A', 
                          'ATM', 'CHEK1', 'BRCA2', 'PALB2', 'RAD51D', 
                          'BRCA1', 'RAD51C', 'BRIP1', 'CHEK2')
final_modelling_data[columns_to_factorise] <- lapply(final_modelling_data[columns_to_factorise], 
                                                     factor, levels = c("0", "-2", "-1", "1", "2"))
final_modelling_data$mhBRCA2 <- factor(final_modelling_data$mhBRCA2)
final_modelling_data$Condition <- ifelse(final_modelling_data$Condition == "Resistant", 1, 0) # Resistant = 1, Sensitive = 0
final_modelling_data$WGD <- ifelse(final_modelling_data$WGD == "TRUE", 1, 0) # TRUE = 1, FALSE = 0
# Functions ---------------------------------------------------------------

logistic_regression <- function(variable, dataset){
  # To conduct logistic regressions based on the variable names given to it for final_modelling_data
  model <- glm(as.formula("Condition ~ " %+% variable), data = dataset, family = binomial(link = logit))
}

# Function to paste variable name 
"%+%" <- function(x,y) paste(x, y, sep = "")

# Logistic Regression Analysis ----------------------------------------------------------------
# Model for one variable
model_test <- glm(Condition ~ BRCA2, family = binomial(link = 'logit'), data = final_modelling_data)

# Model loop through all the variables
predictors <- colnames(final_modelling_data[4:ncol(final_modelling_data)])
model_multiple <- lapply(predictors, logistic_regression, final_modelling_data) # Applying multiple logistic regressions on each variable

## Extracting p-values and coefficients
results <- map_df(model_multiple, tidy) # Extracting results directly from the model, not the summary

## Changing coefficients to odds ratio
results$estimate_odds <- exp(results$estimate)

## Output model results
# write.csv(results, here::here("data/univariate_results/univariate_results.csv"))

# Model for all variables in one model (Do not use)
model_all <- glm(reformulate(paste(predictors, sep = ""), "Condition"), 
                 family = binomial(link = 'logit'), data = final_modelling_data)


# Linear Regression Analysis ----------------------------------------------
# Analysis to see why certain effects are so large

## CX 4

## CX 10




