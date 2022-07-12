
# Libraries ---------------------------------------------------------------


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
# Functions ---------------------------------------------------------------

logistic_regression <- function(variable){
  # To conduct logistic regressions based on the variable names given to it for final_modelling_data
  model <- glm(as.formula("Condition ~ " %+% variable), data = final_modelling_data, family = binomial(link = logit))
  model_summary <- list(summary(model))
}

# Function to paste variable name 
"%+%" <- function(x,y) paste(x, y, sep = "")

# Analysis ----------------------------------------------------------------
# Model for one variable
model_test <- glm(Condition ~ BRCA2, family = binomial(link = 'logit'), data = final_modelling_data)

# Model loop through all the variables
predictors <- colnames(final_modelling_data[4:ncol(final_modelling_data)])
lapply(predictors, logistic_regression)
