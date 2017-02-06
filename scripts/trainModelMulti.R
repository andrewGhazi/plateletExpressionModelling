#Now we train the predixcan model. Apparently it can be slow so I'm going to use mclapply

##### Libraries
library(dplyr)
library(magrittr)
library(parallel)
select = dplyr::select

setwd("~/plateletExpressionModelling/prediXcanModeller/scripts/")

#From the github tutorial: https://github.com/hakyimlab/PredictDBPipeline/wiki/Detailed_Description
# For a single tissue, and single chromosome, use the command:
#   
#   Rscript ./create_model.R $study $expr_RDS $geno $gene_annot $snp_annot \
# $n_k_folds $alpha $out_dir $chrom $snpset $window
# where:
#   
#   - $study is the identifier you are using for the sqlite database
# you will eventually build, e.g. the tissue name for GTEx.
# - $expr_RDS is the path to the expression RDS file for that tissue
# - $geno is preprocessed genotype file for that chromosome
# - $gene_annot is the preprocessed gene annotation file.
# - $snp_annot is the preprocess snp annotation file for that 
# chromosome
# - $n_k_folds is the number of folds for cross-validation
# - $alpha is the mixing parameter for glmnet
# - $out_dir is the output directory (for our purposes this should
#                                     be data/intermediate/model_by_chr/)
# - $chrom is the chromosome number
# - $snpset is the name of the snpset you used
# - $window is the size of the search window to look for snps outside
# of a gene which may impact expression of that gene

trainPrediXcan = mclapply(1:22, function(chr){
  
  #Just putting all the parts here to make the paste in the system call easier to read/change. These include spaces on the end!
  study = 'prax '
  exprRDS = '/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/intermediate/expression_phenotypes/expressionFile.RDS '
  geno = paste0('/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/intermediate/genotypes/praxGenotype.chr', chr, '.txt ')
  gene_annot = '/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/intermediate/annotations/gene_annotation/praxGeneAnno.RDS '
  snp_annot = paste0('/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/intermediate/annotations/snp_annotation/snp_annot.chr', chr, '.RDS ')
  n_k_folds = '10 '
  alpha = '.5 '
  out_dir = '/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/intermediate/model_by_chr/ '
  chrom = paste0(chr %>% as.character(), ' ')
  snpset = 'PRAXOmni5M '
  window = '1000000 '
  
  system(paste0('Rscript create_model_edit.R ', 
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
  return(paste0('chr', chr, ' done'))
}, mc.cores = 12)

#### Combining results across chr
study = 'prax '
all_results_file = '/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/output/allResults/praxResults.txt '
all_betas_file = '/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/output/allBetas/praxBetas.txt '
all_logs_file = '/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/output/allLogs/praxLogs.txt '
all_covariances_file = '/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/output/allCovariances/praxCovariances.txt '
alpha = '0.5 '
snpset = 'PRAXOmni5M'

setwd("/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/output/")

system(paste0('sh /mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/scripts/make_all_results.sh ',
              study,
              all_results_file,
              alpha,
              snpset))

system(paste0('sh /mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/scripts/make_all_betas.sh ',
              study,
              all_betas_file,
              alpha,
              snpset))

system(paste0('sh /mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/scripts/make_all_logs.sh ',
              study,
              all_logs_file))

system(paste0('sh /mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/scripts/make_all_covariances.sh ',
              study,
              all_covariances_file,
              alpha,
              snpset))

#### Create output database
setwd("/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/")

system(paste0('python scripts/make_sqlite_db.py ',
              '--output /mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/output/dbs/firstPraxModel.db ',
              '--results ',
              all_results_file,
              '--construction ',
              all_logs_file,
              '--betas ',
              all_betas_file,
              '--meta ',
              '/mnt/labhome/andrew/plateletExpressionModelling/prediXcanModeller/data/output/allMetaData/firstAttempt.allMetaData.txt'))


#### Filter on FDR significance and 
system(paste0('Rscript scripts/filter_on_significance.R ',
              '~/plateletExpressionModelling/prediXcanModeller/data/output/dbs/firstPraxModel.db ',
              '~/plateletExpressionModelling/prediXcanModeller/data/intermediate/annotations/gene_annotation/praxGeneAnno2.RDS ',
              '~/plateletExpressionModelling/prediXcanModeller/data/output/dbs/firstPraxModelSignificanceFiltered.db'))
