library(tidyverse)
library(pheatmap)
library(gtools)
source('figure_functions.R')

ribosomes <- read_delim('../data/GO_ribosome_biogenesis.tsv', delim = '\t')
omics <- read_csv('../data/tidy_omics.csv')

hm <- average_and_summarise_omics(omics) %>%
  inner_join(ribosomes, by = c('ID' = 'Systematic ID')) %>%
  unite(temp,molecule,time_point.x) %>%
  spread(key = temp, value = avg_fold_change) %>%
  na.omit() %>%
  as.data.frame()

hm[sapply(hm,is.infinite)] <- 0
hm[sapply(hm,is.nan)] <- 0


hm <- hm[,mixedorder(colnames(hm))]

pretty_col <- colorRampPalette(c('blue','darkgrey','yellow'))
cls <- pheatmap(as.matrix(hm[,c(-1:-3)]), cluster_cols = F, color = pretty_col(20), breaks = seq(-5,5,length.out = 20), cutree_rows = 4,
                labels_row = hm[,1]) 
t <- cutree(cls$tree_row, k = 4)
