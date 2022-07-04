
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

## CN 1
# density_cn <- apply(analysis_data_na_removed[, 5:21], 2, density)
# 
# plot(NA, xlim=range(sapply(density_cn, "[", "x")), ylim=range(sapply(density_cn, "[", "y")))
# mapply(lines, density_cn, col=1:length(density_cn))
# legend("topright", legend=names(density_cn), fill=1:length(density_cn))
