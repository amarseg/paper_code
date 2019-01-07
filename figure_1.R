library('tidyverse')
library('ggpubr')
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
