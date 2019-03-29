library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ggt)
#Load data#
all_data <- read_csv('../data/screening_data/live_cell/output_data_cleaned.csv') %>%
  rename('ID' = `Systematic ID`) %>%
  #filter(Time < 4) %>%
  group_by(ID, Time) %>%
  summarise(avg_area = median(AreaShape_Area))

tp0 <- filter(all_data, Time == 0)

max_size <- all_data %>%
  group_by(ID) %>%
  summarise(max_area = max(avg_area)) %>%
  left_join(tp0, by = 'ID', suffix = c('_max','_tp0')) %>%
  mutate(max_ratio = max_area/avg_area)

hit_list <- read_csv('../data/screening_data/output_rep1/summary_rep1_pval_0.15.csv') %>%
  filter(hits == 'hit')

second_list <- read_csv('../data/screening_data/output_rep2/statistics_rep2_pvalue_0.15.csv') %>%
  filter(hits == 'hit')

all_hits <- inner_join(hit_list, second_list, by = 'Systematic ID', suffix = c('_rep1','_rep2')) %>%
  select(`Systematic ID`,z_score_rep1, z_score_rep2)  %>%
  left_join(max_size, by = c('Systematic ID' = 'ID')) 

int <- select(all_hits, z_score_rep1, z_score_rep2, max_ratio, `Systematic ID`)
int$max_ratio <- log2(int$max_ratio)

col <- colorRampPalette(c('blue','gray','red'))
cl <- pheatmap(int[,-4], cluster_cols = F, labels_row = int$`Systematic ID`, cutree_rows = 10, breaks= seq(-3.5,3.5,1),
               color = col(7))


#######Figure for 1st screen##########
hit_list <- read_csv('../data/screening_data/output_rep1/summary_rep1_pval_0.15.csv')
second_list <- read_csv('../data/screening_data/output_rep2/statistics_rep2_pvalue_0.15.csv')

all_hits <- inner_join(hit_list, second_list, by = 'Systematic ID', suffix = c('_rep1','_rep2')) %>%
  mutate(hit_compilation = paste(hits_rep1, hits_rep2, sep = '-'))

ggplot(all_hits, aes(x = mean_area_rep1, y = mean_area_rep2, colour = hit_compilation, shape = hit_compilation)) +
  geom_point(size = 2) +
  theme_classic() +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  xlab('Mean Area Replicate 1 [px^2]') +
  ylab('Mean Area Replicate 2 [px^2]')

ggsave('../figure_ouput/screening_1.pdf')


