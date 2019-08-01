library(tidyverse)
library(pheatmap)
library(gtools)
source('figure_functions.R')
library(ggrepel)

load_essential_list <- function(){
  essential <- read_tsv('../data/essential_genes.tsv') %>%
    add_column(phenotype = 'essential')
  non_ess <- read_tsv('../data/non_essential.tsv') %>%
    add_column(phenotype = 'non-essential')
  
  pheno_table <- bind_rows(essential, non_ess)
  return(pheno_table)
}


gene_lists <- load_gene_lists()


haplo <- read_csv('../data/haploinsufficient/tidy_kim.csv')
ribosomes <- read_delim('../data/GO_ribosome_biogenesis.tsv', delim = '\t')
omics <- read_csv('../data/tidy_omics.csv')


##########Heatmap ribosome biogenesis#################

hm <- average_and_summarise_omics(omics)[[1]] %>%
  inner_join(ribosomes, by = c('ID' = 'Systematic ID')) %>%
  unite(temp,molecule,time_point.x) %>%
  spread(key = temp, value = avg_fold_change) %>%
  na.omit() %>%
  as.data.frame() 

hm <- hm[,grep(colnames(hm), pattern = 'RNA|ID')]

hm[sapply(hm,is.infinite)] <- 0
hm[sapply(hm,is.nan)] <- 0


hm <- hm[,mixedorder(colnames(hm))]

pretty_col <- colorRampPalette(c('blue','darkgrey','yellow'))
cls <- pheatmap(as.matrix(hm[,c(-1:-3)]), cluster_cols = F, color = pretty_col(20), breaks = seq(-5,5,length.out = 20), kmeans_k = 3,
                labels_row = hm[,1]) 

k_means_cls <- data.frame(IDs = hm[,1], k_means_cluster = cls$kmeans$cluster)

#Plot proteins in k_means

summary_prot <- average_and_summarise_omics(omics)[[2]] %>%
  filter(molecule == 'Protein') %>%
  inner_join(k_means_cls, by = c('ID'= 'IDs')) %>%
  transform(k_means_cluster = as.factor(k_means_cluster))

ggplot(summary_prot, aes(x = time_point.x, y = log2foldchange)) +
  geom_point(alpha = 0.8, color = 'grey') +
  stat_summary(fun.y = 'mean', geom = 'line', aes(group = ID)) +
  stat_summary(fun.y = 'mean', geom = 'line', color = 'red') +
  facet_grid(rows = vars(k_means_cluster)) +
  ylim(c(-5,5))

##################Plotting haplo data in ribosomes#################

avg_omics <- average_and_summarise_omics(omics)[[2]]

haplo_ribo <- inner_join(avg_omics, haplo, by ='ID') %>%
  inner_join(ribosomes, by = c('ID' = 'Systematic ID'))

ggplot(haplo_ribo, aes(x = time_point.x, y = log2foldchange, label = ID)) +
  geom_point(alpha = 0.8, color = 'grey') +
  stat_summary(fun.y = 'mean', geom = 'line', aes(group = ID)) +
  stat_summary(fun.y = 'mean', geom = 'line', color = 'red') +
  facet_grid(rows = vars(behaviour)) +
  ylim(c(-5,5)) 


################Paralogue ribosomes##################################
ribo_table <- read_csv('../data/RP.table.280911.csv')
ribo_for_hm <- data.frame(row.names = ribo_table$Sp_name, ribo_type = ribo_table$Sp_family)

hm <- average_and_summarise_omics(omics)[[1]] %>%
  inner_join(ribo_table, by = c('ID' = 'Sp_name')) %>%
  unite(temp,molecule,time_point.x) %>%
  spread(key = temp, value = avg_fold_change) %>%
  na.omit() %>%
  as.data.frame() 

hm <- hm[,grep(colnames(hm), pattern = 'RNA|ID')]

hm[sapply(hm,is.infinite)] <- 0
hm[sapply(hm,is.nan)] <- 0


hm <- hm[,mixedorder(colnames(hm))]

# pretty_col <- colorRampPalette(c('blue','darkgrey','yellow'))
# cls <- pheatmap(as.matrix(hm[,c(-1:-3)]), cluster_cols = F, color = pretty_col(20), breaks = seq(-5,5,length.out = 20),
#                 labels_row = hm[,1], annotation_row = ribo_for_hm) 
# 
# k_means_cls <- data.frame(IDs = hm[,1], k_means_cluster = cls$kmeans$cluster)

phenotypes <- load_essential_list()

thing <- average_and_summarise_omics(omics)[[1]] %>%
  inner_join(ribo_table, by = c('ID' = 'Sp_name')) %>%
  inner_join(phenotypes, by = c('ID' = 'Systematic ID')) %>%
  left_join(haplo, by = 'ID')

ggplot(thing, aes(x = time_point.x, y = avg_fold_change, colour = phenotype)) +
  geom_point() +
  facet_wrap(~Sp_family)

ggplot(thing, aes(x = time_point.x, y = avg_fold_change, colour = behaviour)) +
  geom_point() +
  facet_wrap(~Sp_family)


ggplot(thing, aes(x = time_point.x, y = avg_fold_change, colour = phenotype)) +
  #geom_point() +
  stat_summary(fun.y = 'mean', geom = 'line', aes(group = phenotype)) +
  #ylim(c(-5,5)) +
  #facet_wrap(~phenotype) +
  geom_smooth()

ggplot(thing, aes(x = time_point.x, y = avg_fold_change, colour = behaviour)) +
  #geom_point() +
  stat_summary(fun.y = 'mean', geom = 'line', aes(group = phenotype)) +
  #facet_wrap(~behaviour) +
  geom_smooth()

only_de <- thing %>%
  add_column(de = if_else(.$ID %in% gene_lists$`Systematic ID`, 'DE','Not DE'))

ggplot(only_de , aes(x = time_point.x, y = avg_fold_change, colour = behaviour, shape = de)) +
  geom_point() +
  facet_wrap(~Sp_family)
