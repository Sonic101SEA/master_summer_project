
# Libraries ---------------------------------------------------------------
library(ggplot2)
library(reshape)
library(dplyr)
# Data --------------------------------------------------------------------
heatmap_data <- read.csv(here::here("data/multivariate_results/clustering_results_dataframe.csv"), row.names = 1)


# Preparing Data ----------------------------------------------------------

heatmap_data_scaled <- heatmap_data
heatmap_data_scaled$CX3 <- as.numeric(scale(heatmap_data_scaled$CX3))
heatmap_data_scaled$CX4 <- as.numeric(scale(heatmap_data_scaled$CX4))
heatmap_data_scaled$CX10 <- as.numeric(scale(heatmap_data_scaled$CX10))
heatmap_data_scaled$CX14 <- as.numeric(scale(heatmap_data_scaled$CX14))

heatmap_data_scaled_cn <- as.matrix(heatmap_data_scaled[, c('CX3', 'CX4', 'CX10', 'CX14')])
heatmap_data_scaled_cn_melt <- melt(heatmap_data_scaled_cn)
heatmap_data_scaled_cn_melt$X2 <- factor(heatmap_data_scaled_cn_melt$X2, levels = c('CX3', 'CX4', 'CX10', 'CX14'))

# Plotting heatmap --------------------------------------------------------
heatmap_data_scaled_cn_melt %>%
  ggplot(aes(x = X1, y = X2, fill = value)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

heatmap(heatmap_data_scaled_cn, scale = "none")
