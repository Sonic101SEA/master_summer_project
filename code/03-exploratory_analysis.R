
# Libraries ---------------------------------------------------------------

library(ggplot2)
library(tidyverse)

# Data --------------------------------------------------------------------

analysis_data <- read.csv(here::here("data/final_dataframe.csv"), row.names = 1)

## Remove patients with missing data
analysis_data_na_removed <- analysis_data[complete.cases(analysis_data), ]

## Converting WGD to 1 and 0
# analysis_data_na_removed$WGD <- as.numeric(analysis_data_na_removed$WGD)

# Summary Statistics -----------------------------------------------------
## Distribution of sensitive and resistant patients
plat_outcome_plot <-
barplot(table(analysis_data_na_removed$Condition),
        main = "Distribution of resistant and sensitive patients to platinum therapy",
        ylab = "Frequency")

## Age distribution
summary(analysis_data_na_removed$age)
hist(analysis_data_na_removed$age, breaks = 10)
plot(density(analysis_data_na_removed$age))

## CN
CN_plot <-
analysis_data_na_removed %>%
  pivot_longer(5:21) %>%
  ggplot(aes(x = name, y = value, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~ factor(name, levels = c("CX1", "CX2", "CX3", "CX4", "CX5", 
                                       "CX6", "CX7", "CX8", "CX9", "CX10", 
                                       "CX11", "CX12", "CX13", "CX14", "CX15",
                                       "CX16", "CX17")), scales = "free") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Boxplots of copy number signature activity split between the resistant and sensitive groups", x = "", y = "Signature activity")


# density_cn <- apply(analysis_data_na_removed[, 5:21], 2, density)
# 
# plot(NA, xlim=range(sapply(density_cn, "[", "x")), ylim=range(sapply(density_cn, "[", "y")))
# mapply(lines, density_cn, col=1:length(density_cn))
# legend("topright", legend=names(density_cn), fill=1:length(density_cn))

## WGD
wgd_plot <-
analysis_data_na_removed %>%
  count(Condition, WGD) %>%
  ggplot(aes(x = Condition, y = n, fill = WGD)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  geom_text(aes(label=n), position = position_dodge(width = 1), vjust = -1.0) +
  theme_classic() + 
  labs(title = "Distribution of WGD among resistant and sensitive groups", x = "Condition", y = "No. of patients", fill= "WGD event")

## Gene level calls

### Plotting by condition
gene_level_calls_plot_by_condition <-
  analysis_data_na_removed %>%
  pivot_longer(22:34) %>%
  ggplot(aes(x = value, fill = Condition)) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~ name, scales = "free") +
  labs(title = "Distribution of gene level calls among resistant and sensitive groups for genes of interest", 
       x = "Gene copy changes", y = "No. of patients", fill = "Gene level calls")

### Plotting by the level calls
gene_level_calls_plot_by_calls <-
analysis_data_na_removed %>%
  pivot_longer(22:34) %>%
  ggplot(aes(x = Condition, fill = factor(value, levels = c("-2", "-1", "0", "1", "2")))) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~ name, scales = "free") +
  labs(title = "Distribution of gene level calls among resistant and sensitive groups for genes of interest", 
       x = "Gene copy changes", y = "No. of patients", fill = "Gene level calls")

## Mutations in genes from SNV data
