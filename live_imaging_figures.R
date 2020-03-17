library(tidyverse)
library(sicegar)
library(nlme)
library(broom)
library(ggpubr)
library(cowplot)

a4_width = 210
a4_height = 297

gene_names <- read_tsv('../data/sysID2product.tsv', skip = 1,
                       col_names = c('Systematic_ID','Name','Other_names','Description'))


lineage_data <-  read_csv('../data/live_imaging_data/new_lineage_cell_data.csv') %>%
  filter(cell_n != 0) %>%
  mutate(cell_n = paste0(Systematic_ID,';', cell_n)) %>%
  group_by(Systematic_ID, Metadata_Well, cell_n) %>%
  mutate(number_of_timepoints = n()) %>%
  ungroup()

##############################Fit logistic growth model using nls, plot results###############

filtered_lineage <- lineage_data %>%
  filter(number_of_timepoints <= 20 & number_of_timepoints > 5)

lm_models <- lmList(AreaShape_Area ~ Metadata_Time|cell_n, 
                    data = filtered_lineage)

sig_models <- nlsList(AreaShape_Area ~ SSlogis(Metadata_Time, 
                                               Asym, 
                                               xmid, 
                                               scal)|cell_n,
                      data = filtered_lineage)

model_coeff_lm <- lapply(lm_models, tidy) %>%
  bind_rows(.id = 'id') %>%
  separate(id, into = c('Systematic_ID','cell_n'),
           sep = ';') %>%
  left_join(gene_names, by = 'Systematic_ID') 



model_coeff_sig <- sapply(sig_models, tidy) %>%
  bind_rows(.id = 'id') %>%
  separate(id, into = c('Systematic_ID','cell_n'),
           sep = ';') 

max_values <- model_coeff_sig %>%
  filter(term == 'Asym') %>%
  mutate(p_adj = p.adjust(p.value)) 


growth_rates <- model_coeff_sig %>%
  filter(term == 'scal') %>%
  mutate(growth_rate = -1/estimate, p_adj = p.adjust(p.value)) %>%
  write_csv('growth_rate.csv')

od_50 <- model_coeff_sig %>%
  filter(term == 'xmid') %>%
  mutate(p_adj = p.adjust(p.value))

ggplot(filter(max_values, p_adj < 0.05),
       aes(x =reorder(Systematic_ID, estimate, FUN = 'median', na.rm = T), y = estimate)) +
  geom_boxplot()+ 
  scale_y_continuous(limits = c(250,8000)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


max_values %>%
  group_by(Systematic_ID) %>%
  summarise(median_est = median(estimate, na.rm = T)) %>%
  arrange(desc(median_est)) %>%
  write_csv('max_values.csv')


ggplot(growth_rates, 
       aes(x =reorder(Systematic_ID, growth_rate, FUN = 'median', na.rm = T), y = growth_rate)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_params <- full_join(max_values, growth_rates, by = c('Systematic_ID','cell_n'),
                        suffix = c('_max_val','_gr')) %>%
  full_join(od_50,by = c('Systematic_ID','cell_n'), 
            suffix = c('_od50','') ) %>%
  write_csv('all_model_parameters.csv')

ggplot(all_params, aes(x = estimate_max_val, y = estimate_gr)) +
  geom_point()

ggplot(all_params, aes(x = estimate, y = estimate_gr)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

ggplot(all_params, aes(x = estimate, y = estimate_max_val)) +
  geom_point()

data_with_params <- filtered_lineage %>%
  inner_join()


slopes <- model_coeff_lm %>%
  filter(term == 'Metadata_Time') %>%
  inner_join(all_params, by = c('Systematic_ID','cell_n'),
             suffix = c('_slope','')) %>%
  select_at(vars(-contains('term'))) %>%
  mutate(p_adj = p.adjust(p.value)) %>%
  write_csv('summary_both_models.csv')

ggplot(slopes, aes(x = estimate_slope, y =estimate_gr)) +
  geom_point()

g


slopes_for_output <- model_coeff_lm %>%
  filter(term == 'Metadata_Time') %>%
  mutate(p_adj = p.adjust(p.value)) %>%
  filter(p_adj < 0.05) %>%
  group_by(Systematic_ID) %>%
  summarise(median_slope = median(estimate), sd_slope = sd(estimate, na.rm = T)) %>%
  arrange(desc(median_slope)) %>%
  View()

slopes_p_val_filtered <- model_coeff_lm %>%
  filter(term == 'Metadata_Time') %>%
  mutate(p_adj = p.adjust(p.value)) 

p1<-ggplot(filter(slopes_p_val_filtered ,p_adj <0.05), 
           aes(x = reorder(Systematic_ID, estimate, FUN = 'median', na.rm = T), y = estimate))+
  geom_boxplot(outlier.size=0.25,lwd=0.25) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 6),
        axis.text = element_text(size = 4)) +
  xlab('')+
  ylab('Linear elongation rate')


max_sizes <- filtered_lineage %>%
  group_by(Systematic_ID, cell_n) %>%
  summarise(max_size = max(AreaShape_Area))

p2 <- ggplot(max_sizes, 
             aes(x = reorder(Systematic_ID, max_size, FUN = 'median', na.rm = T), y = max_size)) +
  geom_boxplot(outlier.size=0.25,lwd=0.25) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 6),
        axis.text = element_text(size = 4)) +
  xlab('')+
  ylab('Maximum size [px]')

bop<-plot_grid(p1, p2,
               ncol = 1)

ggsave(plot = bop,
       filename = 'boxplots_models.pdf',
       width = 210/2,
       units = 'mm')

###Plot trajectories and models of nas6 compared to wild type as an example##############
int_ids <- c('SPAC23H3.05c','SPAC6C3.08','wt')

subset_lineage <- filtered_lineage %>%
  filter(Systematic_ID %in% int_ids) %>%
  separate(cell_n, into = c('Systematic_ID','cell_n'),
           sep = ';')
  

ggplot(subset_lineage, aes(x = Metadata_Time, y = AreaShape_Area, colour = Systematic_ID,
                         group =  interaction(Systematic_ID, cell_n))) +
  geom_smooth(method = 'lm', aes(group = Systematic_ID)) +
  #geom_line(alpha = 0.25) +
  theme_pubr() +
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 4))

ggsave(filename = 'trace_example.pdf',
       height = a4_height/5,
       width = a4_width/3,
       units = 'mm')


#######################Variability analysis#############################



filtered_data <- read_csv('../data/live_imaging_data/new_lineage_cell_data.csv') %>%
  filter(cell_n != 0) %>%
  mutate(cell_n = paste0(Systematic_ID, cell_n)) %>%   
  group_by(Systematic_ID, Metadata_Well, cell_n) %>%
  mutate(number_of_timepoints = n()) %>%
  ungroup()


gr_threshold <- 50

proportion_data <- read_csv('../data/live_imaging_data/summary_both_models.csv') %>%
  group_by(Systematic_ID, well) %>%
  summarise(n_growing = sum(estimate_slope > gr_threshold),
            n_not_growing = sum(estimate_slope <= gr_threshold)) %>%
  mutate(proportion = n_growing/(n_growing + n_not_growing)) %>%
  arrange(proportion)


ggplot(proportion_data , aes(x = reorder(Systematic_ID, proportion)
                             , y = proportion)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


cv_data <- read_csv('analysis_live_imaging/output_data/summary_both_models.csv') %>%
  group_by(Systematic_ID) %>%
  summarise(mean_gr = mean(estimate_slope),
            sd_gr = sd(estimate_slope)) %>%
  mutate(cv_gr = sd_gr/mean_gr) 

ggplot(cv_data , aes(x = reorder(Systematic_ID, cv_gr)
                     , y = cv_gr)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


################################Batch effects for wt###############################################

only_wt <- filtered_data %>%
  filter(Systematic_ID == 'wt') %>%
  unite("new_id", Systematic_ID, Metadata_Row, cell_n, Metadata_Well,  sep = ';')

filtered_lineage <- only_wt %>%
  filter(number_of_timepoints <= 20 & number_of_timepoints > 5)

lm_models <- lmList(AreaShape_Area ~ Metadata_Time|new_id, 
                    data = filtered_lineage)


model_coeff_lm <- lapply(lm_models, tidy) %>%
  bind_rows(.id = 'id') %>%
  separate(id, into = c('id','batch','cell_n', 'well'),
           sep = ';')  

model_coeff_lm %>%
  filter(term == 'Metadata_Time') %>%
  ggplot(aes(x = interaction(well, batch), y = estimate, fill = batch)) +
  geom_boxplot() 

########################All data with batch effect#################

all_data_new_id <- filtered_data %>%
  unite("new_id", Systematic_ID, Metadata_Row, cell_n, Metadata_Well,  sep = ';') %>%
  filter(number_of_timepoints <= 20 & number_of_timepoints > 5)

lm_models <- lmList(AreaShape_Area ~ Metadata_Time|new_id, 
                    data = all_data_new_id)


model_coeff_lm <- lapply(lm_models, tidy) %>%
  bind_rows(.id = 'id') %>%
  separate(id, into = c('id','batch','cell_n', 'well'),
           sep = ';')  

model_coeff_lm %>%
  filter(term == 'Metadata_Time') %>%
  ggplot(aes(x =reorder(interaction(id, well), estimate, FUN = 'median', na.rm = T), y = estimate, fill = batch)) +
  geom_boxplot() +
  geom_boxplot(data = filter(model_coeff_lm, id == 'wt' & term == 'Metadata_Time'), aes(fill = 'black')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~batch, scales = 'free_x')


model_coeff_lm %>%
  filter(term == 'Metadata_Time') %>%
  ggplot(aes(x =batch, y = estimate, fill = batch)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  stat_compare_means(ref.group = 'row_1', label = 'p.signif')



gr_threshold <- 50

proportion_data <- model_coeff_lm %>%
  filter(term == 'Metadata_Time') %>%
  group_by(id, well, batch) %>%
  summarise(n_growing = sum(estimate > gr_threshold),
            n_not_growing = sum(estimate <= gr_threshold)) %>%
  mutate(proportion = n_growing/(n_growing + n_not_growing)) %>%
  arrange(proportion)


ggplot(proportion_data , aes(x = reorder(interaction(id, batch, well), proportion)
                             , y = proportion)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


cv_data <- read_csv('analysis_live_imaging/output_data/summary_both_models.csv') %>%
  group_by(Systematic_ID) %>%
  summarise(mean_gr = mean(estimate_slope),
            sd_gr = sd(estimate_slope)) %>%
  mutate(cv_gr = sd_gr/mean_gr) 

ggplot(cv_data , aes(x = reorder(Systematic_ID, cv_gr)
                     , y = cv_gr)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
########################PCA plot##################
tidy_pca_nest <- all_data_new_id %>%
  select(new_id, AreaShape_Area, Metadata_Time) %>%
  spread(value = 'AreaShape_Area', key = Metadata_Time) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames('new_id') %>%
  as.data.frame()

pcs <- prcomp(tidy_pca_nest) %>%
  as.data.frame()


plot(pcs)
pca_components <- broom::tidy(pcs, matrix = "samples") %>%
  spread(PC, value, sep = "") %>%
  separate(row, into = c('id','batch','cell_n', 'well'),
           sep = ';')  


ggplot(pca_components, aes(PC1, PC2, colour = batch)) +
  geom_point(alpha = 0.2) +
  theme_bw() 

####################Max size###############
all_data_new_id %>%
  group_by(new_id)  %>%
  filter(AreaShape_Area == max(AreaShape_Area)) %>%
  rename(max_time = Metadata_Time, max_area = AreaShape_Area) %>%
  separate(new_id, into = c('id','batch','cell_n', 'well'),
           sep = ';') %>%
  ggplot(aes(x = reorder(interaction(id, batch, well), max_area, FUN = 'median'), y = max_area)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

