
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
  analysis_data_na_removed %>%
  count(Condition) %>%
  ggplot(aes(x = Condition, y = n)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label=n), position = position_dodge(width = 1), vjust = -1.0) +
  labs(title = "Distribution of resistant and sensitive patients",
       x = "Condition", y = "No. of patients")

ggsave(here::here("graphs/plat_outcome_distribution.png"), plat_outcome_plot)

## Age distribution
summary(analysis_data_na_removed$age)
hist(analysis_data_na_removed$age, breaks = 10)
plot(density(analysis_data_na_removed$age))

## CN

### Non-stratified

CN_non_stratified_plot <-
  
  analysis_data_na_removed %>%
  pivot_longer(5:21) %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  facet_wrap(~ factor(name, levels = c("CX1", "CX2", "CX3", "CX4", "CX5",
                                       "CX6", "CX7", "CX8", "CX9", "CX10",
                                       "CX11", "CX12", "CX13", "CX14", "CX15",
                                       "CX16", "CX17")), scales = "free") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(title = "Boxplots of copy number signature activity", 
       x = "", y = "Signature activity")

#### Output of plot
ggsave(here::here("graphs/cn_activity_non_stratified_distribution.png"), CN_non_stratified_plot)

### Stratified by therapy outcome
CN_stratified_plot <-

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
  labs(title = "Boxplots of copy number signature activity split between the resistant and sensitive groups", 
       x = "", y = "Signature activity")

#### Output of plot
ggsave(here::here("graphs/cn_activity_stratified_distribution.png"), CN_stratified_plot)

# density_cn <- apply(analysis_data_na_removed[, 5:21], 2, density)
# 
# plot(NA, xlim=range(sapply(density_cn, "[", "x")), ylim=range(sapply(density_cn, "[", "y")))
# mapply(lines, density_cn, col=1:length(density_cn))
# legend("topright", legend=names(density_cn), fill=1:length(density_cn))

## WGD

### Plotting by non-stratified
wgd_plot_non_strat <-
  
  analysis_data_na_removed %>%
  count(WGD) %>%
  ggplot(aes(x = WGD, y = n)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  geom_text(aes(label=n), position = position_dodge(width = 1), vjust = -1.0) +
  theme_classic() + 
  scale_x_discrete(labels = c("No", "Yes")) +
  labs(title = "Distribution of WGD events", 
       x = "Presence of WGD", y = "No. of patients")

#### Output plot
ggsave(here::here("graphs/wgd_distribution_non_stratified.png"), wgd_plot_non_strat)

### Plotting by condition stratified
wgd_plot_strat <-

  analysis_data_na_removed %>%
  count(Condition, WGD) %>%
  ggplot(aes(x = Condition, y = n, fill = WGD)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  geom_text(aes(label=n), position = position_dodge(width = 1), vjust = -1.0) +
  scale_fill_discrete(labels = c("No", "Yes")) + 
  labs(title = "Distribution of WGD among resistant and sensitive groups", 
        x = "Condition", y = "No. of patients", fill= "WGD event")

#### Output plot
ggsave(here::here("graphs/wgd_distribution_stratified.png"), wgd_plot_strat)

## Gene level calls

### Plotting by condition non-stratified

gene_level_calls_plot_non_strat <-
  
  analysis_data_na_removed %>%
  pivot_longer(22:34) %>%
  ggplot(aes(x = value)) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~ name, scales = "free") +
  labs(title = "Distribution of gene level calls among patients for genes of interest", 
       x = "Gene copy changes", y = "No. of patients", fill = "Gene level calls")

#### Output plot
ggsave(here::here("graphs/gene_level_calls_distribution_non_stratified.png"), gene_level_calls_plot_non_strat)

### Plotting by condition stratified
gene_level_calls_plot_strat <-

  analysis_data_na_removed %>%
  pivot_longer(22:34) %>%
  ggplot(aes(x = value, fill = Condition)) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~ name, scales = "free") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Distribution of gene level calls among resistant and sensitive groups for genes of interest", 
        x = "Gene copy changes", y = "No. of patients", fill = "Gene level calls")

#### Output plot
ggsave(here::here("graphs/gene_level_calls_distribution_stratified.png"), gene_level_calls_plot_strat)

### Plotting by the level calls stratified by therapy outcome
gene_level_calls_plot_by_calls <-

  analysis_data_na_removed %>%
  pivot_longer(22:34) %>%
  ggplot(aes(x = Condition, fill = factor(value, levels = c("-2", "-1", "0", "1", "2")))) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~ name, scales = "free") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Distribution of gene level calls among resistant and sensitive groups for genes of interest", 
        x = "Gene copy changes", y = "No. of patients", fill = "Gene level calls")

#### Output plot
ggsave(here::here("graphs/gene_level_calls_distribution_stratified_by_level_calls.png"), gene_level_calls_plot_by_calls)

## Mutations in genes from SNV data

mutations_genes_snv_plot_strat <-

  analysis_data_na_removed %>%
  pivot_longer(36:39) %>%
  ggplot(aes(x = factor(value), fill = Condition)) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~name, scales = "free") +
  scale_x_discrete(labels = c("No", "Yes")) +
  labs(title = "Distribution of occurrence of mutations in genes of interest from SNV data stratified by condition",
       x = "Presence of mutation", y = "No. of patients")

#### Output plot
ggsave(here::here("graphs/mutations_genes_snv_stratified.png"), mutations_genes_snv_plot_strat)

### Without splitting to sensitive and resistant groups
mutations_genes_snv_plot_non_strat <-
  
  analysis_data_na_removed %>%
  pivot_longer(36:39) %>%
  ggplot(aes(x = factor(value))) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~name, scales = "free") +
  scale_x_discrete(labels = c("No", "Yes")) +
  labs(title = "Distribution of occurrence of mutations in genes of interest from SNV data",
       x = "Presence of mutation", y = "No. of patients")

#### Output plot
ggsave(here::here("graphs/mutations_genes_snv_non_stratified.png"), mutations_genes_snv_plot_non_strat)


