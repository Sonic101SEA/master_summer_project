
# Libraries ---------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyverse)

# Data --------------------------------------------------------------------

analysis_data <- read.csv(here::here("data/final_dataframe.csv"), row.names = 1)

## Remove patients with missing data
analysis_data_na_removed <- analysis_data[complete.cases(analysis_data), ]

## Converting WGD to 1 and 0
# analysis_data_na_removed$WGD <- as.numeric(analysis_data_na_removed$WGD)

# Summary Statistics -----------------------------------------------------
## Distribution of sensitive and resistant patients
barplot(table(analysis_data_na_removed$Condition),
        main = "Distribution of resistant and sensitive patients to platinum therapy",
        ylab = "Frequency")

## Age distribution
summary(analysis_data_na_removed$age)
hist(analysis_data_na_removed$age, breaks = 10)
plot(density(analysis_data_na_removed$age))

## CN
analysis_data_na_removed %>%
  pivot_longer(5:21) %>%
  ggplot(aes(x = name, y = value, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~ name, scales = "free") +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")


ggplot(analysis_data_na_removed, aes(x = CX1))

# density_cn <- apply(analysis_data_na_removed[, 5:21], 2, density)
# 
# plot(NA, xlim=range(sapply(density_cn, "[", "x")), ylim=range(sapply(density_cn, "[", "y")))
# mapply(lines, density_cn, col=1:length(density_cn))
# legend("topright", legend=names(density_cn), fill=1:length(density_cn))

## WGD
analysis_data_na_removed %>%
  count(Condition, WGD) %>%
  ggplot(aes(x = Condition, y = n, fill = WGD)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  geom_text(aes(label=n), position = position_dodge(width = 1), vjust = -1.0) +
  theme_classic() + 
  labs(title = "Distribution of WGD among resistant and sensitive groups", x = "Condition", y = "No. of patients", fill="WGD event")

