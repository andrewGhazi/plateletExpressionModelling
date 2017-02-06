# Let's try to run PrediXcan (both the PRAX and GTEx models) to try to find
# associations between predicted gene expression and phenotypes.
# Following along the instructions given here: https://github.com/hakyimlab/PrediXcan/tree/master/Software

library(dplyr)
library(readr)
library(magrittr)
library(biomaRt)
library(stringr)


outDirPath = '~/plateletExpressionModelling/outputs/phenoAssocOutputs/ '
dosageDir = '~/plateletExpressionModelling/data/genotypeByChr/ ' #Needs to point to a directory containing genotype file 

system(paste0('python ',
              '/usr/local/src/PrediXcan-master/Software/PrediXcan.py --predict ',
              '--output_dir ', outDirPath,
              '--dosages ', dosageDir,
              '--dosages_prefix praxGenotype.chr ',
              '--weights ', '~/plateletExpressionModelling/data/GTEx-V6p-HapMap-2016-09-08/TW_Whole_Blood_0.5.db ',
              '--samples ', '~/plateletExpressionModelling/data/genotypeByChr/samples.txt'))

# Association tests require you to be in the software's directory
setwd('/usr/local/src/PrediXcan-master/Software/')
phen = phenotypes %>% select(ADP:vWF) %>% names

for (i in 1:length(phen)) {
  system(paste0('mkdir ',
                '~/plateletExpressionModelling/outputs/phenoAssocOutputs/', 
                phen[i]))
}

for (i in 1:length(phen)) {
  system(paste0('python ',
                '/usr/local/src/PrediXcan-master/Software/PrediXcan.py --assoc ',
                '--output_dir ', outDirPath %>% gsub(' ', '', .), phen[i], '/ ',
                '--pred_exp ', '~/plateletExpressionModelling/outputs/phenoAssocOutputs/predicted_expression.txt ',
                '--pheno ', '~/plateletExpressionModelling/data/phenotypeFile.tsv ', 
                '--linear ',
                '--pheno_name ', phen[i]))
}

for (i in 1:length(phen)) {
  if (i == 1) {
    results = read_delim(paste0(outDirPath %>% gsub(' ', '', .), phen[i], '/association.txt'), 
                         skip = 1,
                         delim = ' ',
                         col_names = c('gene', 'beta', 't', 'p', 'SE')) %>% 
      mutate(phenotype = phen[i])
  } else {
    tmp = read_delim(paste0(outDirPath %>% gsub(' ', '', .), phen[i], '/association.txt'), 
                         skip = 1,
                         delim = ' ',
                         col_names = c('gene', 'beta', 't', 'p', 'SE')) %>% 
      mutate(phenotype = phen[i])
    
    results %<>% rbind(tmp)
  }
}

results %<>% mutate(q = p.adjust(p, method = 'fdr'),
                    bonf.p = p.adjust(p, method = 'bonferroni'),
                    gene = str_extract(gene, '[A-Z0-9]+(?=\\.)'))
results %>% filter(q < .05) %>% arrange(q)

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
resGenes = getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'), 
                 filters = 'ensembl_gene_id', 
                 values = (results %>% filter(q < .05) %>% .$gene), 
                 mart = mt) %>% 
  as.tbl %>% 
  rename(gene = ensembl_gene_id)

sigResults = results %>% 
  filter(q < .05) %>% 
  left_join(resGenes, by = 'gene') %>% 
  arrange(q) %>%
  select(gene, hgnc_symbol, phenotype, beta, t, SE, p, q, bonf.p) %T>%
  write_tsv('~/plateletExpressionModelling/outputs/phenoAssocOutputs/gtexModelPraxPhenotypeAssociations.tsv')

