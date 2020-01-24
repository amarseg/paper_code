library('tidyverse')
library('ggpubr')
library('wesanderson')
library('cowplot')
##########Kinetics of cell growth#############

gene_names <- read_tsv('../data/sysID2product.tsv',
                       col_names = c('Systematic_ID','Name','Other_Names','Description'),
                       skip = 1) %>%
  mutate(Name = case_when(is.na(Name) ~ Systematic_ID,
         TRUE ~ Name))

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
hits <- read_csv('../data/summary_all_hits.csv') 

common_data <- inner_join(screening_data_1, screening_data_2, by = 'Systematic ID', suffix = c('_rep1','_rep2')) %>%
  add_column(hit = ifelse(.$`Systematic ID` %in% hits$`Systematic ID`, 'hit','not hit'))


ggplot(common_data, aes(x = mean_area_rep1, y = mean_area_rep2, colour = hit)) +
  geom_point() +
  xlab('Mean Area Replicate 1 [px^2]') +
  ylab('Mean Area Replicate 2 [px^2]') +
  theme_light() +
  scale_color_brewer(palette = 'Accent')

ggsave('../figure_ouput/Screening_hits.pdf')


#######################Areas boxplot################################
areas_1 <- read_csv('../data/screening_data/output_rep1/cell_areas.csv',
                    col_names = c('n','AreaShape_Area','Systematic ID','Metadata_Plate_Name','Well'),
                    skip = 1) %>%
  add_column(replicate = 'rep1') %>%
  select(-n)  
areas_2 <- read_csv('../data/screening_data/output_rep2/areas_rep2.csv') %>%
  add_column(replicate = 'rep2')


wt_areas <- read_csv('../data/wild_type_areas.csv') %>%
  select('AreaShape_Area','Metadata_Plate_Name','Well') %>%
  add_column(Systematic_ID = 'wt')

all_areas <- bind_rows(areas_1, areas_2) %>%
  rename(Systematic_ID = `Systematic ID`) %>%
  bind_rows(wt_areas) %>%
  filter(Systematic_ID == 'wt'| Systematic_ID %in% hits$`Systematic ID`) %>%
  left_join(gene_names, by = 'Systematic_ID') %>%
  mutate(Name = case_when(is.na(Name) ~ 'wt',
                             TRUE ~ Name))



ggplot(all_areas, aes(x = reorder(Name, AreaShape_Area, FUN = 'median'), y = AreaShape_Area)) +
  geom_boxplot() +
  geom_boxplot(data = filter(all_areas, Name == 'wt'), colour = 'red') +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  ylab('Area [Number of pixels]') +
  xlab('') +
  scale_y_continuous(limits = c(0,2000))

ggsave('boxplot_hits_primary.pdf',
       width = 210,
       height  = 297/3,
       units = 'mm')
  
  
  
