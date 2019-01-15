library(tidyverse)
library(pheatmap)
library(gtools)
source('figure_functions.R')
##################Heatmap of omics##############################

omics <- read_csv('../data/tidy_omics.csv')

hm <- average_and_summarise_omics(omics) %>% 
  unite(temp,molecule,time_point.x) %>%
  spread(key = temp, value = avg_fold_change) %>%
  na.omit() %>%
  as.data.frame()

hm[sapply(hm,is.infinite)] <- 0

hm <- hm[,mixedorder(colnames(hm))]

pretty_col <- colorRampPalette(c('blue','darkgrey','yellow'))
pheatmap(as.matrix(hm[,-1]), cluster_cols = F, color = pretty_col(20), breaks = seq(-10,10,length.out = 20))
