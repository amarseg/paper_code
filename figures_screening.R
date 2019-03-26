library(tidyverse)
library(pheatmap)
#Load data#
all_data <- read_csv('live_cell_imaging_cp/output_data_cleaned.csv') %>%
  rename('ID' = `Systematic ID`) %>%
  #filter(Time < 4) %>%
  group_by(ID, Time) %>%
  summarise(avg_area = median(AreaShape_Area))

annot_data <- read_csv('final_hit_list.csv')

tp0 <- filter(all_data, Time == 0)

max_size <- all_data %>%
  group_by(ID) %>%
  summarise(max_area = max(avg_area)) %>%
  left_join(annot_data, by = c('ID' = 'Systematic ID')) %>%
  left_join(tp0, by = 'ID', suffix = c('_max','_tp0')) %>%
  mutate(max_ratio = max_area/avg_area)


hit_list <- read_csv('output_rep1/summary_rep1_pval_0.15.csv') %>%
  filter(hits == 'hit')
second_list <- read_csv('output_rep2/statistics_rep2_pvalue_0.15.csv') %>%
  filter(hits == 'hit')

all_hits <- inner_join(hit_list, second_list, by = 'Systematic ID', suffix = c('_rep1','_rep2')) %>%
  select(`Systematic ID`,z_score_rep1, z_score_rep2)  %>%
  left_join(max_size, by = c('Systematic ID' = 'ID')) %>%
  write_csv('summary_all_hits')

int <- select(all_hits, z_score_rep1, z_score_rep2, max_ratio, `Systematic ID`)
int$max_ratio <- log2(int$max_ratio)

col <- colorRampPalette(c('blue','gray','red'))
cl <- pheatmap(int[,-4], cluster_cols = F, labels_row = int$`Systematic ID`, cutree_rows = 10, breaks= seq(-3.5,3.5,1),
               color = col(7))

hcl <- data.frame(cluster = cutree(cl$tree_row, k = 10), id = int[,4]) %>%
  write_csv('live_cell_imaging_cp/result_clustering.csv')


########Figure for 1st screen##########


