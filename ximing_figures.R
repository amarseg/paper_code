library(tidyverse)
library(pheatmap)
library(gtools)
source('figure_functions.R')
library(clusterProfiler)
library(ggrepel)

isx_data <- read_csv('../data/isx_data.csv') %>%
  mutate(time = time - 1)

int_genes <- c('taf2','taf3','taf5','taf11','taf12','spt7','wbp8','rpb1','rpb2','rpb9')

gene_ids <- read_tsv('ftp://ftp.pombase.org/pombe/names_and_identifiers/gene_IDs_names_products.tsv',
                     col_names = c('PomBaseID','second_id','name','chromosome','description',
                                   'UniProt','Type','Synonyms'))

omics <- read_csv('../data/tidy_omics.csv') %>%
  average_and_summarise_omics()

int_omics <- omics[[2]] %>%
  left_join(gene_ids, by = c('ID' = 'PomBaseID')) %>%
  left_join(isx_data, by = c('replicate' = 'rep', 'time_point.x' = 'time')) %>%
  filter(name %in% int_genes)

ggplot(subset(int_omics, molecule == 'Protein'), aes(x = Length_Erode.M03..4., y = read_number.x, group = replicate, colour = time_point.x)) +
  facet_wrap(~name, scales = 'free') +
  geom_line() +
  geom_point() 

ggsave('xi-ming_figure_read_number4.pdf')

ggplot(subset(int_omics, molecule == 'Protein'), aes(x = Length_Erode.M03..4., y = log2foldchange, group = replicate, colour = replicate)) +
  facet_wrap(~name, scales = 'free') +
  geom_line() +
  geom_point() 
ggsave('xi-ming_figure_fold_change5.pdf')
