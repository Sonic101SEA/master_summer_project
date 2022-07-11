
# Libraries ---------------------------------------------------------------


# Data --------------------------------------------------------------------
modelling_data <- read.csv(here::here("data/final_dataframe.csv"), row.names = 1)

## Remove patients with missing data
modelling_data_removed_na<- modelling_data[complete.cases(modelling_data), ]

## Removing variables that have low variability: CX6, 8, 13, 17 and mhBRCA1, mhCHEK2, mhPALB2
columns_to_remove <- c('CX6', 'CX8', 'CX13', 'CX17', 'mhBRCA1', 'mhCHEK2', 'mhPALB2')
final_modelling_data <- modelling_data_removed_na[, !(names(modelling_data_removed_na) %in% columns_to_remove)]


# Functions ---------------------------------------------------------------

logistic_regression <- function(outcome, variable){
  model <- glm(outcome ~ variable, family = binomial(link = 'logit'))
 # summary(model)$coefficients[2]
  
}


# Analysis ----------------------------------------------------------------
set.seed(1)

model_test <- glm(Condition ~ CX1, family = binomial(link = 'logit'), data = final_modelling_data)

for (i in 4:ncol(final_modelling_data)){
  model <- glm(final_modelling_data[, 'Condition'] ~ final_modelling_data[, i], family = binomial(link = 'logit'))
}
