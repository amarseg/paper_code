average_and_summarise_omics <- function(tidy_omics)
{
  require(tidyverse)
  tp0 <- dplyr::filter(tidy_omics,time_point == 0)
  
  normalised_omics <- tidy_omics %>%
    left_join(tp0, by = c('replicate','ID','molecule')) %>%
    mutate(read_number.y = read_number.x  / read_number.y) %>%
    mutate(log2foldchange = log2(read_number.y))  %>%
    transform(time_point.x = as.numeric(time_point.x))
  
  avg_omics <- normalised_omics %>%
    group_by(time_point.x, ID, molecule) %>%
    summarise(avg_fold_change = mean(log2foldchange, na.rm = T)) %>%
    return()
}

load_gene_lists <- function()
{
  require(tidyverse)
  
  up_rna <- read_tsv('../data/up_rna.tsv') %>%
    add_column( type = 'Up RNA')
  
  down_rna <- read_tsv('../data/down_rna.tsv') %>%
    add_column( type = 'Down RNA')
  
  up_prot <- read_tsv('../data/up_prot.tsv') %>%
    add_column( type = 'Up prot')
  
  down_prot <- read_tsv('../data/down_prot.tsv') %>%
    add_column( type = 'Down prot')
  
  omics_lists <- bind_rows(up_rna, down_rna, up_prot, down_prot) %>%
    return()
  
}
