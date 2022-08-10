
# Libraries ---------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(table1)
library(tableone)
library(ggpubr)
# Data --------------------------------------------------------------------

analysis_data <- read.csv(here::here("data/final_dataframe.csv"), row.names = 1)

## Remove patients with missing data
analysis_data_na_removed <- analysis_data[complete.cases(analysis_data), ]

## Converting WGD to 1 and 0
# analysis_data_na_removed$WGD <- as.numeric(analysis_data_na_removed$WGD)


# Table One Summary -------------------------------------------------------
## Using Table1
table1_dataframe <- analysis_data_na_removed

levels_in_data_snv <- c(0, 1)
labels_for_data_snv <- c("Mutation not present", "Mutation present")

levels_in_data_gene_level <- c(-2, -1, 0, 1, 2)
labels_for_data_gene_level <- c("Two copy loss", "One copy loss", "No change", "One copy gain", "Two copy gain")

changing_labels <- function(data, levels_input, labels_input){
  data = factor(data, levels = levels_input, labels = labels_input)
}

table1_dataframe[c("mhBRCA1", "mhBRCA2", "mhCHEK2", "mhPALB2")] <- lapply(table1_dataframe[c("mhBRCA1", "mhBRCA2", "mhCHEK2", "mhPALB2")], 
                                                                          changing_labels, levels_in_data_snv, labels_for_data_snv)
table1_dataframe[c("BARD1", "FAM175A", "NBN", "MRE11A", 
                   "ATM", "CHEK1", "BRCA2", "PALB2", "RAD51D", 
                   "BRCA1", "RAD51C", "BRIP1", "CHEK2")] <- lapply(table1_dataframe[c("BARD1", "FAM175A", "NBN", "MRE11A", 
                                                                                      "ATM", "CHEK1", "BRCA2", "PALB2", "RAD51D", 
                                                                                      "BRCA1", "RAD51C", "BRIP1", "CHEK2")], 
                                                                   changing_labels, levels_in_data_gene_level, labels_for_data_gene_level) 

label(table1_dataframe$age) <- "Age"

table1_descriptive_statistics <- 
  table1(~ age + CX1 + CX2 + CX3 + CX4 + CX5 + CX6 + CX7 + CX8 + CX9 + CX10 + CX11 + CX12 + CX13 +
         CX14 + CX15 + CX16 + CX17 + BARD1 + FAM175A + NBN + MRE11A + ATM + CHEK1 + BRCA2 + PALB2 + RAD51D + 
         BRCA1 + RAD51C + BRIP1 + CHEK2 + WGD + mhBRCA1 + mhBRCA2 + mhCHEK2 + mhPALB2 | Condition, data = table1_dataframe, overall = "Total")


## Using tableone
## Vector of variables to summarise
# vars_summarise <- c("Condition", "age", "CX1", "CX2", "CX3", "CX4", "CX5", "CX6", "CX7", "CX8", "CX9", "CX10", "CX11",
#                     "CX12", "CX13", "CX14", "CX15", "CX16", "CX17", "WGD", "BARD1", "FAM175A", "NBN", "MRE11A", "ATM", "CHEK1", 
#                     "BRCA2", "PALB2", "RAD51D", "BRCA1", "RAD51C", "BRIP1", "CHEK2",
#                     "mhBRCA1", "mhBRCA2", "mhCHEK2", "mhPALB2")
# 
# ## Vector of categorical variables that need factoring
# 
# vars_cat <- c("BARD1", "FAM175A", "NBN", "MRE11A", "ATM", "CHEK1", 
#               "BRCA2", "PALB2", "RAD51D", "BRCA1", "RAD51C", "BRIP1", "CHEK2",
#               "mhBRCA1", "mhBRCA2", "mhCHEK2", "mhPALB2")
# 
# CreateTableOne(vars = vars_summarise, factorVars = vars_cat, data = analysis_data_na_removed, strata = "Condition")

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

ggsave(here::here("graphs/exploratory/plat_outcome_distribution.pdf"), plat_outcome_plot)

## Age distribution

analysis_data_na_removed_age_grouped <-
  analysis_data_na_removed %>%
  mutate(Age_class = cut(age, breaks = seq(40, 85, by = 5)))

# analysis_data_na_removed %>%
#   ggplot(aes(x = age)) +
#   geom_boxplot(width = 0.01) +
#   xlim(25, 100) +
#   labs(x = "Age")
# 
# summary(analysis_data_na_removed$age)
# hist(analysis_data_na_removed$age, breaks = 10)
# plot(density(analysis_data_na_removed$age))

age_labels <- c("41-45", "46-50", '51-55', '56-60', '61-65', '66-70', '71-75', '76-80', '80 and above')

age_distribution <-
analysis_data_na_removed_age_grouped %>%
  ggplot(aes(x = Age_class)) +
  geom_bar() +
  scale_x_discrete(labels = age_labels) +
  labs(x = "Age Group", y = "No. of Patients")

ggsave(here::here("graphs/exploratory/age_distribution.pdf"), age_distribution)

## CN

### Non-stratified same scale

CN_non_stratified_plot_same_scale <-
  
  analysis_data_na_removed %>%
  pivot_longer(5:21) %>%
  mutate(name_reorder = factor(name, levels = c("CX1", "CX2", "CX3", "CX4", "CX5",
                                                "CX6", "CX7", "CX8", "CX9", "CX10",
                                                "CX11", "CX12", "CX13", "CX14", "CX15",
                                                "CX16", "CX17"))) %>%
  ggplot(aes(x = reorder(name_reorder, value, median), y = value)) + # Reordering by median values
  geom_boxplot() +
  labs(x = "", y = "Signature activity", colour = "Copy Number")

ggsave(here::here("graphs/exploratory/cn_activity_non_stratified_distribution_same_scale.pdf"), CN_non_stratified_plot_same_scale)

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
  labs(x = "", y = "Signature activity")

#### Output of plot
ggsave(here::here("graphs/exploratory/cn_activity_non_stratified_distribution.pdf"), CN_non_stratified_plot)

cn_signatures_non_strat_figure <- ggarrange(CN_non_stratified_plot, CN_non_stratified_plot_same_scale,
                                  labels = c("A", "B"),
                                  ncol = 2, nrow = 1)

ggsave(here::here("graphs/exploratory/cn_non_stratified_figure.pdf"), cn_signatures_non_strat_figure, limitsize = TRUE)
### Stratified by therapy outcome on same scale

CN_stratified_plot_same_scale <-
  
  analysis_data_na_removed %>%
  pivot_longer(5:21) %>%
  ggplot(aes(x = reorder(name, value, median), y = value, fill = Condition)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Copy number signature activity split between resistant and sensitive groups on same scale", 
       x = "", y = "Signature activity")

ggsave(here::here("graphs/exploratory/cn_activity_stratified_distribution_same_scale.pdf"), CN_stratified_plot_same_scale)

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
  labs(title = "Boxplots of copy number signature activity split between resistant and sensitive groups", 
       x = "", y = "Signature activity")

#### Output of plot
ggsave(here::here("graphs/exploratory/cn_activity_stratified_distribution.pdf"), CN_stratified_plot)

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
ggsave(here::here("graphs/exploratory/wgd_distribution_non_stratified.pdf"), wgd_plot_non_strat)

### Plotting by condition stratified
wgd_plot_strat <-

  analysis_data_na_removed %>%
  count(Condition, WGD) %>%
  ggplot(aes(x = WGD, y = n, fill = Condition)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  geom_text(aes(label=n), position = position_dodge(width = 1), vjust = -1.0) +
  scale_x_discrete(labels = c("No", "Yes")) +
  labs(title = "Distribution of WGD among resistant and sensitive groups",
        x = "Presence of WGD", y = "No. of patients", fill= "Condition")

#### Output plot
ggsave(here::here("graphs/exploratory/wgd_distribution_stratified.pdf"), wgd_plot_strat)

## Gene level calls

### Plotting by condition non-stratified

gene_level_calls_plot_non_strat <-
  
  analysis_data_na_removed %>%
  pivot_longer(22:34) %>%
  ggplot(aes(x = factor(value))) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~ name, scales = "free") +
  scale_x_discrete(drop = FALSE) +
  labs(title = "Distribution of gene level calls among patients for genes of interest", 
       x = "Gene copy changes", y = "No. of patients", fill = "Gene level calls")

#### Output plot
ggsave(here::here("graphs/exploratory/gene_level_calls_distribution_non_stratified.pdf"), gene_level_calls_plot_non_strat)

### Plotting by condition stratified
gene_level_calls_plot_strat <-

  analysis_data_na_removed %>%
  pivot_longer(22:34) %>%
  ggplot(aes(x = factor(value), fill = Condition)) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~ name, scales = "free") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Distribution of gene level calls among resistant and sensitive groups for genes of interest", 
        x = "Gene copy changes", y = "No. of patients", fill = "Gene level calls")

#### Output plot
ggsave(here::here("graphs/exploratory/gene_level_calls_distribution_stratified.pdf"), gene_level_calls_plot_strat)

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
ggsave(here::here("graphs/exploratory/gene_level_calls_distribution_stratified_by_level_calls.pdf"), gene_level_calls_plot_by_calls)

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
ggsave(here::here("graphs/exploratory/mutations_genes_snv_stratified.pdf"), mutations_genes_snv_plot_strat)

### Without splitting to sensitive and resistant groups
mutations_genes_snv_plot_non_strat <-
  
  analysis_data_na_removed %>%
  pivot_longer(36:39) %>%
  ggplot(aes(x = factor(value))) +
  geom_bar(position = position_dodge(preserve = "single")) +
  facet_wrap(~name, scales = "free") +
  scale_x_discrete(labels = c("No", "Yes")) +
  labs(title = "Distribution of occurrence of mutations in genes of interest from SNV and Indel data",
       x = "Presence of mutation", y = "No. of patients")

#### Output plot
ggsave(here::here("graphs/exploratory/mutations_genes_snv_non_stratified.pdf"), mutations_genes_snv_plot_non_strat)


