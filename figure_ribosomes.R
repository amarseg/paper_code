library(tidyverse)
source('figure_functions.R')

ribosomes <- read_delim('../data/GO_ribosome_biogenesis.tsv', delim = '\t')
omics <- read_csv('../data/tidy_omics.csv')

normalised_omics <- average_and_summarise_omics(omics) %>%
  inner_join(ribosomes, by = c('ID' = 'Systematic ID')) %>%
  unite(temp,molecule,time_point.x) %>%
  spread(key = temp, value = avg_fold_change) %>%
  na.omit() %>%
  as.data.frame()

