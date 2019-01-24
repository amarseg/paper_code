library(tidyverse)
source('figure_functions.R')

high_transcripts <- read.delim('Z:/Pers_Amalia/chip_ncRNA/methil_genes.txt')

omics <- read_csv('../data/tidy_omics.csv')

normalised_omics <- average_and_summarise_omics(omics)

high <- normalised_omics %>%
  filter(ID %in% high_transcripts[,1] & molecule == 'RNA')


high[is.nan(high$avg_fold_change),]$avg_fold_change <- NA
high[is.infinite(high$avg_fold_change),]$avg_fold_change <- NA

high$chromosome <- 'Unknown' 
high[str_detect(high$ID, 'AC'),]$chromosome <- 'I'
high[str_detect(high$ID, 'BC'),]$chromosome <- 'II'
high[str_detect(high$ID, 'CC'),]$chromosome <- 'III'

ggplot(high, aes(x = time_point.x, y = avg_fold_change, colour = chromosome)) +
  geom_point( ) +
  geom_line(aes(group = ID), color = 'grey') +
  stat_summary( fun.y=mean, geom="line", colour="black", size = 1.25, linetype = 'dashed') +
  theme_bw()
ggsave('../figure_ouput/methylated_transcripts.pdf')
