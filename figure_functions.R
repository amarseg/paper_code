average_and_summarise_omics <- function(tidy_omics)
{
  require(tidyverse)
  tidy_omics <- tidy_omics %>%
    mutate(read_number = read_number + 0.00001)
  
  tp0 <- dplyr::filter(tidy_omics,time_point == 0)
  
  normalised_omics <- tidy_omics %>%
    left_join(tp0, by = c('replicate','ID','molecule')) %>%
    mutate(read_number.y = read_number.x  / read_number.y) %>%
    mutate(log2foldchange = log2(read_number.y))  %>%
    transform(time_point.x = as.numeric(time_point.x))
  
  avg_omics <- normalised_omics %>%
    group_by(time_point.x, ID, molecule) %>%
    summarise(avg_fold_change = mean(log2foldchange, na.rm = T))
  
  return(list(avg_omics, normalised_omics))
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

load_go <- function(){
  require(AnnotationDbi)
  require(GO.db)
  
  go_data <- read_delim('ftp://ftp.pombase.org/pombe/annotations/Gene_ontology/gene_association.pombase.gz',
                        skip = 44, delim = '\t', col_names = F)
  term2gene <- go_data[,c(5,2)]
  
  term2name <- AnnotationDbi::select(GO.db, keys = keys(GO.db), columns = c('GOID','TERM'))
  
  go_database <- list(term2gene = term2gene, term2name = term2name)
  return(go_database)
}

load_omics_data <- function()
{
  require(tidyverse)
  require(DESeq2)
  
  #RNa procesing
  rna <- read_csv('rna_seq.csv') 
  size_factors <- rna %>%
    dplyr::select(-ID) %>%
    estimateSizeFactorsForMatrix()
  
  rna[,-1] <- sweep(rna[,-1],MARGIN=2,FUN="/",STATS=size_factors) %>%
    as.data.frame()
  
  final_rna <- rna %>%
    gather(key = sample_name, value  = read_number, -ID) %>%
    separate(sample_name, into = c('not needed','time_point','replicate'), sep = '_') %>%
    dplyr::select(-`not needed`) %>%
    mutate(time_point = str_extract(time_point, pattern = '[:digit:]{1,2}')) %>%
    mutate(replicate = str_extract(replicate, pattern = '[:digit:]{1}')) %>%
    add_column(molecule = 'RNA')
  
  
  #Protein procesing
  prot <- read_delim('SQ_Results_PROTEIN.tsv', delim = '\t') %>%
    dplyr::select(proteinName, b018p004AM_T0_01:b018p004AM_T11_G3_02) %>%
    filter(str_detect(proteinName, pattern = 'SP')) %>%
    gather(key = sample_name, value = read_number, -proteinName) %>%
    separate(sample_name, into = c('no','time_point','replicate','technical_replicate'), sep = '_') %>%
    separate(proteinName, into = c('ID','gene_name','chr','desc'), sep = '\\|') %>%
    mutate(time_point = str_extract(time_point, pattern = '[:digit:]{1,2}')) %>%
    mutate(replicate = str_extract(replicate, pattern = '[:digit:]{1}')) %>%
    dplyr::select(-(gene_name:desc), -no) %>%
    add_column(molecule = 'Protein')
  
  t <- prot[which(is.na(prot$technical_replicate)),]
  t$technical_replicate <- t$replicate  
  t$replicate <-1
  
  prot[which(is.na(prot$technical_replicate)),] <- t
  
  mean_prot <- prot %>%
    group_by(ID, time_point, replicate, molecule) %>%
    summarise_at('read_number', mean)
  
  omics_dataset <- bind_rows(final_rna, mean_prot) %>%
    return()
  
  write_csv(omics_dataset, 'tidy_omics.csv')
}