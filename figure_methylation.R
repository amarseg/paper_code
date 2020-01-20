library(tidyverse)
source('figure_functions.R')

high_transcripts <- read.delim('Z:/Pers_Amalia/chip_ncRNA/methil_genes.txt')

omics <- read_csv('../data/tidy_omics.csv')

normalised_omics <- average_and_summarise_omics(omics)

high <- normalised_omics[[1]] %>%
  filter(ID %in% high_transcripts[,1] & molecule == 'RNA')


high[is.nan(high$avg_fold_change),]$avg_fold_change <- NA
high[is.infinite(high$avg_fold_change),]$avg_fold_change <- NA

high$chromosome <- 'Unknown' 
high[str_detect(high$ID, 'SPA'),]$chromosome <- 'I'
high[str_detect(high$ID, 'SPB'),]$chromosome <- 'II'
high[str_detect(high$ID, 'SPC'),]$chromosome <- 'III'

ggplot(high, aes(x = time_point.x, y = avg_fold_change, colour = chromosome)) +
  geom_point( ) +
  geom_line(aes(group = ID), color = 'grey') +
  stat_summary( fun.y=mean, geom="line", colour="black", size = 1.25, linetype = 'dashed') +
  theme_bw()
ggsave('../figure_ouput/methylated_transcripts.pdf')


new_methyl_list <- read_csv('../data/sam_list.csv')

high <- normalised_omics[[1]] %>%
  filter(molecule == 'RNA') %>%
  inner_join(new_methyl_list, by = c('ID' = 'id'))

high[is.nan(high$avg_fold_change),]$avg_fold_change <- NA
high[is.infinite(high$avg_fold_change),]$avg_fold_change <- NAi

ggplot(high, aes(x = time_point.x, y = avg_fold_change, colour = cluster)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_line(aes(group = ID), color = 'grey', size = 1) +
  geom_point(size = 2.5) +
  stat_summary( fun.y=mean, geom="line", colour="black", size = 1.25, linetype = 'dashed') +
  theme_bw() +
  facet_wrap(~cluster) +
  scale_colour_brewer(palette = 'Dark2')

ggsave('../figure_ouput/methyl_figure_for_sam.pdf')

ggplot(filter(high, cluster == 'Chromosome II cluster'), aes(x = time_point.x, y = avg_fold_change, colour = cluster)) +
  geom_point( ) +
  geom_line(aes(group = ID), color = 'grey') +
  stat_summary( fun.y=mean, geom="line", colour="black", size = 1.25, linetype = 'dashed') +
  theme_bw()

