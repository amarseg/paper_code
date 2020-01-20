library(tidyverse)
library(pheatmap)
library(gtools)
source('figure_functions.R')
library(clusterProfiler)
library(ggrepel)
set.seed(502)
##################Heatmap of omics##############################
gene_lists <- load_gene_lists()
go_db <- load_go()

omics <- read_csv('../data/tidy_omics.csv')

hm <- average_and_summarise_omics(omics)[[1]] %>% 
  unite(temp,molecule,time_point.x) %>%
  spread(key = temp, value = avg_fold_change) %>%
  na.omit() %>%
  as.data.frame() %>%
  filter(ID %in% gene_lists$`Systematic ID`)

hm[sapply(hm,is.infinite)] <- 0

hm <- hm[,mixedorder(colnames(hm))]

pretty_col <- colorRampPalette(c('blue','darkgrey','yellow'))
heatmap_significant_molecules<- pheatmap(as.matrix(hm[,-1]), cluster_cols = F, color = pretty_col(10), breaks = seq(-5,5,length.out = 10),
         clustering_method = 'ward.D2',
         cutree_rows = 4)

cls <- tibble(cluster = cutree(heatmap_significant_molecules$tree_row, k = 4),
              Systematic_ID = hm$ID)


cl_go_enrichment <- clusterProfiler::compareCluster(data = cls, 
                                     Systematic_ID ~ cluster, 
                                     fun = 'enricher', 
                                     TERM2GENE = go_db$term2gene, 
                                     TERM2NAME = go_db$term2name)

cl_go_enrichment_df <- as.data.frame(cl_go_enrichment)
dotplot(cl_go_enrichment, showCategory = 50)

####Plot genes belonging to GO category#####
cl_go_enrichment_df %>%
  mutate(geneID = str_split(geneID, pattern = '/')) %>%
  unnest() %>%
  inner_join(average_and_summarise_omics(omics)[[2]], by = c('geneID' = 'ID')) %>%
  mutate(broad_description = case_when(str_detect(Description ,'ribo|translation') ~ 'Translation',
                                       str_detect(Description , 'ubiquitin|endopeptidase') ~ 'Protein degradation',
                                       str_detect(Description, 'glycolisis|tricarboxylic') ~ 'Metabolism',
                                       str_detect(Description, 'post|anaphase') ~ 'Cell Cycle',
                                       TRUE ~ Description)) %>%
  mutate(broad_description = str_to_sentence(broad_description)) %>%
  group_by(cluster, broad_description, molecule, time_point.x, replicate) %>%
  summarise(median_fold_change = median(log2foldchange)) %>%
  filter(cluster!= 4) %>%
  ggplot(aes(x = time_point.x , y = median_fold_change, group = broad_description, colour = broad_description)) +
  geom_point(alpha = 0.5) +
  stat_summary(fun.y = 'median', geom = 'line', size = 1.2) +
  facet_wrap(~cluster*molecule, scales = 'free') +
  theme_bw() +
  scale_color_brewer(palette = 'Dark2')

ggsave('pseudo_GO.pdf')
#################Enriquecimiento usando GO#########################


sum_omics <- average_and_summarise_omics(omics)[[1]] %>%
  na.omit()

plot <- sum_omics %>%
  inner_join(gene_lists, by = c('ID' = 'Systematic ID')) %>%
  separate(type, into = c('Direction','Molecule'), sep = ' ',remove = F)

ggplot(plot, aes(x = time_point.x, y = avg_fold_change, fill = Molecule, color = Molecule)) +
  geom_boxplot(aes(group = interaction(time_point.x, molecule), fill = Molecule, color = Molecule)) +
  stat_summary(fun.y=median, geom="line") +
  facet_grid(~Direction)

##############################clustered transcripts##########################

de_transcripts <- filter(gene_lists, type == 'Up RNA' | type == 'Down RNA')

de_trans <- average_and_summarise_omics(omics) %>%
  filter(molecule == 'RNA' & ID %in% de_transcripts$`Systematic ID`) %>%
  spread(key = time_point.x, value = avg_fold_change) %>%
  select(-molecule)

de_trans[sapply(de_trans,is.infinite)] <- 0
de_trans[sapply(de_trans,is.nan)] <- 0


transcript_clusters <- pheatmap(de_trans[,-1], cluster_cols = F, color = pretty_col(20), breaks = seq(-10,10,length.out = 20), kmeans_k =  5)
cl <- data.frame(ID = de_trans$ID, cluster = transcript_clusters$kmeans$cluster)

t <- clusterProfiler::compareCluster(data = cl, ID ~ cluster, fun = 'enricher', TERM2GENE = go_db$term2gene, TERM2NAME = go_db$term2name)
enrich_df <- as.data.frame(t)
dotplot(t, showCategory = 50)

ggplot(enrich_df, aes(y = -log10(p.adjust), x = Count, colour = cluster, label = Description)) +
  geom_point() +
  geom_label_repel()

proteomics <- average_and_summarise_omics(omics)[[1]] %>%
  inner_join(cl, by = 'ID')

ggplot(proteomics, aes(x = time_point.x, y = avg_fold_change, colour = molecule)) +
  stat_summary(fun.y=median, geom="line")+
  geom_boxplot(aes(group = interaction(time_point.x, molecule))) +
  facet_wrap(~cluster) +
  theme_bw()
ggsave('../output/clusters_in_protein.pdf')

ord_trans <- de_trans[order(cl$cluster),]
n_cl <- as.vector(table(cl$cluster))
t <- n_cl
for(i in 1:length(n_cl))
{
  t[i] <- sum(n_cl[1:i])
}
pheatmap(ord_trans[,-1], cluster_cols = F, cluster_rows = F, color = pretty_col(20), breaks = seq(-10,10,length.out = 20), gaps_row = t,
         labels_row = NULL)

#######################Get protein clusters and do heatmaps############################
prot_cl <- proteomics %>%
  filter(molecule == 'Protein')

for(i in unique(prot_cl$cluster))
{
  sub_data <- prot_cl %>%
    filter(cluster == i) %>%
    spread(key = time_point.x, value = avg_fold_change)
  file_name = paste0('Proteomics_cluster_',i,'.pdf')
  
  pdf(file_name)
  pheatmap(sub_data[,c(-1:-3)], labels_row = sub_data$ID, cluster_cols = F,color = pretty_col(10), breaks = seq(-5,5,length.out = 10), main = i)
  dev.off() 
}
