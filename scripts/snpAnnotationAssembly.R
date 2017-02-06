
#### Libraries
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(data.table)
library(stringr)
library(dtplyr)
select = dplyr::select

load("~/plateletExpressionModelling/data/snpAnnotationAssembly.RData")

for(i in 1:22){
  # chrdat = snpAnnotations %>% filter(chr == paste0('chr', i))
  # chrdat %>% mutate(rsIDwithTab = paste0(RsID, '\t')) %>% 
  #   select(rsIDwithTab) %>% 
  #   write_tsv(path = paste0('data/rsIdsByChromosome/chr', i, '_rsIds.txt'))
  dbSNPdati = read_tsv(paste0('/mnt/bigData2/resources/dbSNP_146_GRCh37p13/byChromosome/', i, '.vcf'),
                       col_names = FALSE)
  names(dbSNPdati) = c('CHROM', 
                       'POS', 
                       'RsID',
                       'REF', 
                       'ALT',
                       'QUAL', 
                       'FILTER',
                       'INFO')
  
  #TERRIBLE CODE INCOMING!
  if (i == 1){
    snpAnnotations %<>% left_join(dbSNPdati, by = 'RsID')
  } else{
    dbSNPdati %<>% filter(RsID %in% snpAnnotations$RsID)
    snpAnnotations[match(dbSNPdati$RsID, snpAnnotations$RsID), c('CHROM', 
                                                                 'POS', 
                                                                 'RsID',
                                                                 'REF', 
                                                                 'ALT',
                                                                 'QUAL', 
                                                                 'FILTER',
                                                                 'INFO')] = dbSNPdati
  }
}

snpAnnotations %>% 
  mutate(VariantID = str_c(CHROM, POS, REF, ALT, 'hg19', sep = '_')) %>% 
  select(CHROM, POS, VariantID, REF, ALT, RsID) %>% 
  na.omit %>% #only cuts out 100,513 / 3,602,58
  write_tsv(path = '~/plateletExpressionModelling/data/snpAnnotationFile.txt')