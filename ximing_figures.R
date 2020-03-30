library(tidyverse)
library(pheatmap)
library(gtools)
source('figure_functions.R')
library(clusterProfiler)
library(ggrepel)

isx_data <- read_csv('../data/isx_data.csv') %>%
  mutate(time = time - 1)

int_genes <- c('clr1','clr2','clr3','mit1','ccq1','chp2','set2')

gene_ids <- read_tsv('ftp://ftp.pombase.org/pombe/names_and_identifiers/gene_IDs_names_products.tsv',
                     col_names = c('PomBaseID','second_id','name','chromosome','description',
                                   'UniProt','Type','Synonyms'))

omics <- read_csv('../data/tidy_omics.csv') %>%
  average_and_summarise_omics()

int_omics <- omics[[2]] %>%
  left_join(gene_ids, by = c('ID' = 'PomBaseID')) %>%
  left_join(isx_data, by = c('replicate' = 'rep', 'time_point.x' = 'time')) %>%
  filter(name %in% int_genes)

ggplot(subset(int_omics, molecule == 'RNA'), aes(x = Length_Erode.M03..4., y = read_number.x, group = replicate, colour = time_point.x)) +
  facet_wrap(~name, scales = 'free') +
  geom_line() +
  geom_point() 

ggsave('SHREC_complex_reads.pdf')

ggplot(subset(int_omics, molecule == 'RNA'), aes(x = Length_Erode.M03..4., y = log2foldchange, group = replicate, colour = replicate)) +
  facet_wrap(~name, scales = 'free') +
  geom_line() +
  geom_point() 
ggsave('SHREC_complex_fold_change.pdf')
