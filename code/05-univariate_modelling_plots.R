
# Libraries ---------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(ggrepel)

# Data --------------------------------------------------------------------

univariate_results <- read.csv(here::here("data/univariate_results/univariate_results.csv"), row.names = 1, stringsAsFactors = FALSE)

# Processing before plotting ----------------------------------------------
# For 0.05
## Removing Intercept rows
univariate_results_rm_interecept <- subset(univariate_results, subset = univariate_results$term != "(Intercept)")

## Creating new column for effect direction
univariate_results_rm_interecept$effect_direction <- "Insignificant"

## if estimate_odds > 1 and pvalue  < 0.05, set as "positive"
univariate_results_rm_interecept$effect_direction[univariate_results_rm_interecept$estimate_odds > 1 &
                                                    univariate_results_rm_interecept$p.value < 0.05] <- "Positive"

## if estimate_odds < 1 and pvalue < 0.05, set as "negative"
univariate_results_rm_interecept$effect_direction[univariate_results_rm_interecept$estimate_odds < 1 &
                                                    univariate_results_rm_interecept$p.value < 0.05] <- "Negative"

## Setting labels for significant values
univariate_results_rm_interecept$label_name <- NA
univariate_results_rm_interecept$label_name[univariate_results_rm_interecept$effect_direction != "Insignificant"] <-
  univariate_results_rm_interecept$term[univariate_results_rm_interecept$effect_direction != "Insignificant"]

# For 0.10
univariate_results_10percent <- univariate_results_rm_interecept
## Creating new column for effect direction
univariate_results_10percent$effect_direction <- "Insignificant"

## if estimate_odds > 1 and pvalue  < 0.05, set as "positive"
univariate_results_10percent$effect_direction[univariate_results_10percent$estimate_odds > 1 &
                                                    univariate_results_10percent$p.value < 0.10] <- "Positive"

## if estimate_odds < 1 and pvalue < 0.05, set as "negative"
univariate_results_10percent$effect_direction[univariate_results_10percent$estimate_odds < 1 &
                                                    univariate_results_10percent$p.value < 0.10] <- "Negative"

## Setting labels for significant values
univariate_results_10percent$label_name <- NA
univariate_results_10percent$label_name[univariate_results_10percent$effect_direction != "Insignificant"] <-
  univariate_results_10percent$term[univariate_results_10percent$effect_direction != "Insignificant"]

# For 0.15
univariate_results_15percent <- univariate_results_rm_interecept
## Creating new column for effect direction
univariate_results_15percent$effect_direction <- "Insignificant"

## if estimate_odds > 1 and pvalue  < 0.05, set as "positive"
univariate_results_15percent$effect_direction[univariate_results_15percent$estimate_odds > 1 &
                                                univariate_results_15percent$p.value < 0.15] <- "Positive"

## if estimate_odds < 1 and pvalue < 0.05, set as "negative"
univariate_results_15percent$effect_direction[univariate_results_15percent$estimate_odds < 1 &
                                                univariate_results_15percent$p.value < 0.15] <- "Negative"

## Setting labels for significant values
univariate_results_15percent$label_name <- NA
univariate_results_15percent$label_name[univariate_results_15percent$effect_direction != "Insignificant"] <-
  univariate_results_15percent$term[univariate_results_15percent$effect_direction != "Insignificant"]

# Plotting ----------------------------------------------------------------
# Plotting a volcano plot, with positive variables labelled for 5% significance
univariate_results_rm_interecept %>%
  ggplot(aes(x = estimate, y = -log10(p.value), col = effect_direction, label = label_name)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  coord_cartesian(xlim = c(-30, 30)) +
  geom_text_repel()

# For 10% significance
univariate_results_10percent %>%
  ggplot(aes(x = estimate, y = -log10(p.value), col = effect_direction, label = label_name)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.10), col = "red") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  coord_cartesian(xlim = c(-30, 30)) +
  geom_text_repel()

# For 15% significance
univariate_results_15percent %>%
  ggplot(aes(x = estimate, y = -log10(p.value), col = effect_direction, label = label_name)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.15), col = "red") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  coord_cartesian(xlim = c(-30, 30)) +
  geom_text_repel()
