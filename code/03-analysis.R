
# Data --------------------------------------------------------------------

analysis_data <- read.csv(here::here("data/final_dataframe.csv"), row.names = 1)

## Remove patients with missing data
analysis_data_na_removed <- analysis_data[complete.cases(analysis_data), ]

# Summary Statistics -----------------------------------------------------
## Distribution of sensitive and resistant patients
barplot(table(analysis_data_na_removed$Condition))

## Age distribution
summary(analysis_data_na_removed$age)
plot(density(analysis_data_na_removed$age))
