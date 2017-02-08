#Now that we've compared a prax model to GTEx and DGN, let's try comparing the
#true RNA expression of a holdout set of PRAX patients with their predicted
#expression based on a predixcan model trained on the remaining patients.

#We might as well cross-validate phenotype associations in the midst of this I guess.


#### Libraries
library(dplyr)
library(magrittr)
library(ggplot2)
library(readr)
library(tidyr)
library(data.table)
library(parallel); select = dplyr::select; select_ = dplyr::select_; contains = dplyr::contains

#### Main script
setwd('~/plateletExpressionModelling/')

# Read in data
subjects = system(paste0('head -1 ', 
                         '/mnt/labhome/simonlm/projects/PRAX/Papers/eQTLpaper/MatrixEQTL/data/genotype.txt'), 
                  intern = TRUE) %>% 
  strsplit('\t') %>% 
  unlist %>% 
  .[-1]

genotype = read_tsv('data/genotypeFile.tsv') #These are both from compileInputs.R
exprData = read_tsv('data/expressionFile.tsv')
kgpToRs = read_delim('~/plateletExpressionModelling/data/InfiniumOmni5-4v1-2_A1_b144_rsids.txt',
                     col_names = TRUE,
                     delim = '\t')
names(kgpToRs)[1] = 'id'

snpAnnotations = read_tsv('~/plateletExpressionModelling/data/snpAnnotationFile.txt')

reverseDosage = function(dosage){
  2 - dosage
}

# Set CV parameters
nfold = 7 # 7-fold because 154 %% 7 = 0

set.seed(123456)
folds = split(sample(subjects), 1:nfold) %>% 
  lapply(sort)

# Run CV
mclapply(2:nfold, function(foldnum){
  setwd('~/plateletExpressionModelling/')
  
  #### Create genotype and expression files consisting of just the data of the training set
  genotype %>% 
    select_(.dots = c('VariantID', folds[-foldnum] %>% unlist %>% as.vector)) %>% #as.vector removes the names
    write_tsv(path = paste0('intermediates/', foldnum, '/trainGenotypeFold', foldnum, '.tsv')) #mutate_each(funs(reverseDosage), vars = -VariantID) %>% #We reverse the dosage because PrediXcan is expecting dosage of the minor allele # Took this out as I fixed it in compileInputs.R
  
  rm(genotype)
  gc()
  
  testGenotypes = fread('/mnt/labhome/simonlm/projects/PRAX/Papers/eQTLpaper/MatrixEQTL/data/genotype.txt',
                    header = TRUE) %>% 
    as.tbl %>% 
    select_(.dots = c("id", folds[[foldnum]])) %>% 
    left_join(kgpToRs, by = 'id', copy = TRUE) %>% #I don't exactly know what the copy = TRUE argument does but it throws and error otherwise
    left_join(snpAnnotations, by = 'RsID', copy = TRUE) %>% 
    na.omit %>% 
    mutate_each(funs(reverseDosage), vars = contains('X')) %>% 
    mutate(MAF = .1) %>% # Just faking this for now
    select(CHROM, RsID, POS, REF, ALT, MAF, contains('X'))
  
  testGenotypes %>% names %>% 
    grep('X', ., value = TRUE) %>% 
    data_frame(FID = ., 
               IID = .) %>%
    write_tsv(paste0('~/plateletExpressionModelling/intermediates/', foldnum, '/test/samples.txt'),
              col_names = FALSE)

  for (i in 1:22) {
    tmp = testGenotypes %>% filter(CHROM == i)
    write_tsv(tmp, paste0("~/plateletExpressionModelling/intermediates/", foldnum, "/test/testGenotypeFold", foldnum, '.chr', unique(tmp$CHROM), ".txt"),
              col_names = FALSE)
    system(paste0('gzip ',
                  '~/plateletExpressionModelling/intermediates/', foldnum, '/test/testGenotypeFold', foldnum, '.chr', i, '.txt'))
    #cat(i)
  }
  rm(testGenotypes)
  gc()
  
  exprData %>% 
    select_(.dots = c('ensembl_id', folds[-foldnum] %>% unlist %>% as.vector)) %>% 
    write_tsv(path = paste0('intermediates/', foldnum, '/trainExprFold', foldnum, '.tsv'))
  
  exprData %>% 
    select_(.dots = c('ensembl_id', folds[foldnum] %>% unlist %>% as.vector)) %>% 
    write_tsv(path = paste0('intermediates/', foldnum, '/testExprFold', foldnum, '.tsv'))
  
  #### Run predixcan model training pre-processing
  setwd('~/plateletExpressionModelling/prediXcanModeller/scripts/')
  
  system(paste0('python split_genotype_by_chr.py ', #split training genotype file by chr
                '~/plateletExpressionModelling/intermediates/', foldnum, '/trainGenotypeFold', foldnum, '.tsv ',
                '~/plateletExpressionModelling/intermediates/', foldnum, '/trainGenotypeFold', foldnum))
  
  system(paste0('Rscript ',
                'expr_to_transposed_RDS.R ', # convert expression file to RDS
                '~/plateletExpressionModelling/intermediates/', foldnum, '/trainExprFold', foldnum, '.tsv ',
                '~/plateletExpressionModelling/intermediates/', foldnum, '/trainExprFold', foldnum, '.RDS'))
  
  system(paste0('python ', #create meta data file
                '~/plateletExpressionModelling/prediXcanModeller/scripts/create_meta_data.py ',
                '--geno ', '~/plateletExpressionModelling/intermediates/', foldnum, '/trainGenotypeFold', foldnum, '.tsv ',
                '--expr ', '~/plateletExpressionModelling/intermediates/', foldnum, '/trainExprFold', foldnum, '.tsv ',
                '--snpset ', 'PRAXOmni5M ',
                '--alpha ', '0.5 ',
                '--n_k_folds ', '10 ',
                '--rsid_label ', 'unknown ',
                '--window ', '1000000 ',
                '--out_prefix ', '~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/allMetaData/praxFold', foldnum))
  
  
  
  #### Create a predixcan model on the training set by chromosome
  for (chrNum in 1:22) {
    study = 'prax '
    exprRDS = paste0('/mnt/labhome/andrew/plateletExpressionModelling/intermediates/', foldnum, '/trainExprFold', foldnum, '.RDS ')
    geno = paste0('/mnt/labhome/andrew/plateletExpressionModelling/intermediates/', foldnum, '/trainGenotypeFold', foldnum, '.chr', chrNum, '.txt ')
    gene_annot = '/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/intermediate/annotations/gene_annotation/praxGeneAnno2.RDS '
    snp_annot = paste0('/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/intermediate/annotations/snp_annotation/snp_annot.chr', chrNum, '.RDS ')
    n_k_folds = '10 '
    alpha = '.5 '
    out_dir = paste0('/mnt/labhome/andrew/plateletExpressionModelling/outputs/model_by_chr_fold/', foldnum, '/ ')
    chrom = paste0(chrNum %>% as.character(), ' ')
    snpset = 'PRAXOmni5M '
    window = '1000000 '
    
    system(paste0('Rscript ',
                  '~/plateletExpressionModelling/scripts/create_model_edit.R ', 
                  study,
                  exprRDS,
                  geno,
                  gene_annot,
                  snp_annot,
                  n_k_folds,
                  alpha,
                  out_dir, 
                  chrom,
                  snpset,
                  window))
    
    
  }
  
  #### Do post-processing to combine chr models
  study = 'prax '
  all_results_file = paste0('/mnt/labhome/andrew/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/allResults/praxFold', foldnum, 'Results.txt ')
  all_betas_file = paste0('/mnt/labhome/andrew/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/allBetas/praxFold', foldnum, 'Betas.txt ')
  all_logs_file = paste0('/mnt/labhome/andrew/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/allLogs/praxFold', foldnum, 'Logs.txt ')
  all_covariances_file = paste0('/mnt/labhome/andrew/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/allCovariances/praxFold', foldnum, 'Covariances.txt ')
  all_meta_data_file = paste0('~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/allMetaData/praxFold', foldnum, '.allMetaData.txt')
  alpha = '0.5 '
  snpset = 'PRAXOmni5M '
  
  setwd("/mnt/labhome/andrew/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/")
  
  system(paste0('sh /mnt/labhome/andrew/plateletExpressionModelling/scripts/make_all_results_edit.sh ',
                study,
                all_results_file,
                alpha,
                snpset,
                foldnum))
  
  system(paste0('sh /mnt/labhome/andrew/plateletExpressionModelling/scripts/make_all_betas_edit.sh ',
                study,
                all_betas_file,
                alpha,
                snpset,
                foldnum))
  
  system(paste0('sh /mnt/labhome/andrew/plateletExpressionModelling/scripts/make_all_logs_edit.sh ',
                study,
                all_logs_file,
                foldnum))
  
  system(paste0('rm ~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/allCovariances/praxFold', foldnum, 'Covariances.txt.gz'))
  system(paste0('sh /mnt/labhome/andrew/plateletExpressionModelling/scripts/make_all_covariances_edit.sh ',
                study,
                all_covariances_file,
                alpha,
                snpset,
                foldnum))
  
  # create SQLite DB from models
  system(paste0('python ',
                '~/plateletExpressionModelling/prediXcanModeller/scripts/make_sqlite_db.py ', 
                '--output ', '~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/dbs/praxFold', foldnum, '.db ',
                '--betas ', all_betas_file,
                '--results ', all_results_file,
                '--construction ', all_logs_file, 
                '--meta ', all_meta_data_file))
  
  # Filter by significance
  system(paste0('Rscript ',
                '~/plateletExpressionModelling/prediXcanModeller/scripts/filter_on_significance.R ',
                '~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/dbs/praxFold', foldnum, '.db ',
                '~/plateletExpressionModelling/prediXcanModeller/data/intermediate/annotations/gene_annotation/praxGeneAnno2.RDS ',
                '~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/dbs/praxFilteredFold', foldnum, '.db'))
  
  #### Run predixcan with the created model to predict expression in the test set
  
  outDirPath = paste0('~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/exprPredictions/', foldnum, '/ ')
  dosageDir = paste0('~/plateletExpressionModelling/intermediates/', foldnum, '/test/ ') #Needs to point to a directory containing genotype files
  
  # toBeGzipped = list.files(paste0('~/plateletExpressionModelling/intermediates/', foldnum, '/test/'),
  #                          pattern = '*.txt')
  # 
  # for (i in 1:length(toBeGzipped)) {
  #   system(paste0('gzip ',
  #                 '~/plateletExpressionModelling/intermediates/', foldnum, '/test/',
  #                 toBeGzipped[i]))
  # }
  
  setwd('/usr/local/src/PrediXcan-master/Software/')
  system(paste0('python ',
                '/usr/local/src/PrediXcan-master/Software/PrediXcan.py --predict ',
                '--output_dir ', outDirPath,
                '--dosages ', dosageDir,
                '--dosages_prefix testGenotypeFold', foldnum, '.chr ',
                '--weights ', '~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/dbs/praxFilteredFold', foldnum, '.db ',
                '--samples ', '~/plateletExpressionModelling/intermediates/', foldnum, '/test/samples.txt'))
  
  #### 2/5/17 We'll write this part later because it's Sunday and I want to chill
  return(paste0('Fold ', foldnum, ' done at ', Sys.time()))
  # See how well the model predicts expression on the test set. Log this information.
  
  # Look for associations with phenotypes in the test set. Log this information.
  
  # save(list = c('exprPredictionModel',
  #               'exprPredictionError',
  #               'phenotypeAssociations'), 
  #      file = paste0('~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/resultsFold', foldnum, '.RData'))
  
}, mc.cores = 3)



##### Performance on the test set looked disappointing. Let's see how it performs on the training sets.

customFun  = function(DF) { #Thanks http://stackoverflow.com/questions/41233173/how-can-i-write-dplyr-groups-to-separate-files
  write_tsv(DF %>% select(-CHROM), paste0("~/plateletExpressionModelling/data/genotypeByChr/praxGenotype.chr", unique(DF$CHROM), ".tsv"))
  return(DF)
}

kgpToRs = read_delim('~/plateletExpressionModelling/data/InfiniumOmni5-4v1-2_A1_b144_rsids.txt',
                     col_names = TRUE,
                     delim = '\t')
names(kgpToRs)[1] = 'id'

snpAnnotations = read_tsv('~/plateletExpressionModelling/data/snpAnnotationFile.txt')

genotypes = fread('/mnt/labhome/simonlm/projects/PRAX/Papers/eQTLpaper/MatrixEQTL/data/genotype.txt',
                  header = TRUE) %>% 
  as.tbl %>% 
  left_join(kgpToRs, by = 'id', copy = TRUE) %>% #I don't exactly know what the copy = TRUE argument does but it throws and error otherwise
  left_join(snpAnnotations, by = 'RsID', copy = TRUE) %>% 
  na.omit %>% 
  mutate_each(funs(reverseDosage), vars = contains('X')) %>% 
  mutate(MAF = .1) %>% # Just faking this for now
  select(CHROM, RsID, POS, REF, ALT, MAF, contains('X')) 

for (foldnum in 1:7) {
  outDirPath = paste0('~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/exprPredictionsFull/', foldnum, '/ ')
  dosageDir = paste0('~/plateletExpressionModelling/data/genotypeByChr/ ') #Needs to point to a directory containing genotype files
  
  setwd('/usr/local/src/PrediXcan-master/Software/')
  system(paste0('python ',
                '/usr/local/src/PrediXcan-master/Software/PrediXcan.py --predict ',
                '--output_dir ', outDirPath,
                '--dosages ', dosageDir,
                '--dosages_prefix praxGenotype.chr ',
                '--weights ', '~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/dbs/praxFilteredFold', foldnum, '.db ',
                '--samples ', '~/plateletExpressionModelling/data/genotypeByChr/samples.txt'))
}

exprData = read_tsv('data/expressionFile.tsv')
for (i in 1:7) {
  #read in the predicted expression for this fold
  predExpr = read_tsv(paste0('~/plateletExpressionModelling/outputs/crossValidatedExpressionAndPhenotypes/exprPredictionsFull/', i, '/predicted_expression.txt')) %>% 
    select(-FID) %>% 
    filter(IID %in% (folds[-i] %>% unlist %>% as.vector)) %>% 
    mutate(type = 'PrediXcan')
  
  rawExpr = exprData %>% 
    filter(ensembl_id %in% names(predExpr)) %>% 
    select_(.dots = c('ensembl_id', folds[-i] %>% unlist %>% as.vector))
  
  genes = rawExpr$ensembl_id
  
  realExpr = rawExpr %>% 
    select(-ensembl_id) %>% 
    t %>% 
    as.data.frame %>% 
    as.tbl
  
  names(realExpr) = genes
  
  realExpr$IID = names(rawExpr)[-1]
  
  realExpr %<>% mutate(type = 'PRAX') %>% select(IID, type, contains('ENSG'))
  
  foldDat = rbind(predExpr, realExpr) %>% 
    select(IID, type, contains('ENSG')) %>% 
    gather(gene, expression, -IID, -type) %>% 
    spread(type, expression)
  
  foldDat %>% ggplot(aes(PrediXcan, PRAX)) + geom_point()
  
  if (i == 1) {
    totDat = foldDat
  } else {
    totDat %<>% rbind(foldDat)
  }
}
#totDat %>% ggplot(aes(PrediXcan, PRAX)) + geom_point()

foldDat %>% 
  split(.$gene) %>% 
  map(~ lm(PRAX ~ PrediXcan, data = .)) %>% 
  map(summary) %>% 
  map_dbl('r.squared')


