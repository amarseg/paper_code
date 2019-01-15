library('tidyverse')
library('ggpubr')
library('wesanderson')
##########Kinetics of cell growth#############

files <- list.files(path = '../data/movie_data/')
data_path <- '../data/movie_data'

data_live_imaging <- data_frame(filename = files) %>% # create a data frame# holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(data_path, .))) # a new data column
  )  %>%
  unnest()


ggplot(data_live_imaging, aes(x = t, y = l)) +
  geom_point( color = 'grey', alpha = 0.8) +
  geom_line(aes(group = filename), color = 'grey', alpha = 0.8) +
  xlab('Time exposed to 1NM-PP1 [min]') +
  ylab('Length [um]') +
  geom_smooth() +
  theme_light()

ggsave('../figure_ouput/live_imaging_plot.pdf')

####################Boxplot showing hits############

screening_data_1 <- read_csv('../data/screening_data/output_rep1/summary_rep1_pval.csv')
screening_data_2 <- read_csv('../data/screening_data/output_rep2/statistics_rep2_pvalue.csv')
hits <- read_csv('../data/screening_data/z_score_hits.csv') %>%
  filter(!is.na(size))

common_data <- inner_join(screening_data_1, screening_data_2, by = 'Systematic ID', suffix = c('_rep1','_rep2')) %>%
  add_column(hit = ifelse(.$`Systematic ID` %in% hits$`Systematic ID`, 'hit','not hit'))


ggplot(common_data, aes(x = mean_area_rep1, y = mean_area_rep2, colour = hit)) +
  geom_point() +
  xlab('Mean Area Replicate 1 [px^2]') +
  ylab('Mean Area Replicate 2 [px^2]') +
  theme_light() +
  scale_color_brewer(palette = 'Accent')

ggsave('../figure_ouput/Screening_hits.pdf')
