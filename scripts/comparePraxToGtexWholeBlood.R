#Now that I've trained a PrediXcan model on the PRAX data, let's just try to
#compare against the PrediXcan people's model that they trained on the GTEx data
#-- specifically in the whole blood tissue

#### Libraries
library(RSQLite)
library(dplyr)
library(magrittr)

#### First let's get the GTEx data from predixcan's website
setwd('~/plateletExpressionModelling/data/')
# system('wget https://s3.amazonaws.com/predictdb/GTEx-V6p-HapMap-2016-09-08.tar.gz')
# system('gunzip GTEx-V6p-HapMap-2016-09-08.tar.gz')
# system('tar -xvf GTEx-V6p-HapMap-2016-09-08.tar') #could have just done tar -zxvf to start, doh

praxDB = dbConnect(drv = SQLite(), dbname = '~/plateletExpressionModelling/prediXcanModeller/data/output/dbs/firstPraxModelSignificanceFiltered.db')
gtexWholeBloodDB = dbConnect(SQLite(), dbname = '~/plateletExpressionModelling/data/GTEx-V6p-HapMap-2016-09-08/TW_Whole_Blood_0.5.db')

#options(max.print = 30)

# dbReadTable(praxDB, 'construction')
# dbReadTable(praxDB, 'extra')
# dbReadTable(praxDB, 'sample_info')
# dbReadTable(praxDB, 'weights')

praxModel = praxDB %>% 
  dbReadTable('weights') %>% 
  as.tbl() %>% 
  mutate(dataSource = 'PRAX')

gtexModel = gtexWholeBloodDB %>% 
  dbReadTable('weights') %>% 
  as.tbl %>% 
  mutate(dataSource = 'GTEx',
         gene = str_extract(gene, '[A-Z0-9]+(?=\\.)')) #regex to remove transcript numbers

bothModel = praxModel %>% rbind(gtexModel)

commonGenes = bothModel %>% 
  group_by(gene) %>% 
  summarise(numSnpsPRAX = sum(dataSource == 'PRAX'),
            numSnpsGTEx = sum(dataSource == 'GTEx')) %>% 
  filter(numSnpsPRAX > 0, numSnpsGTEx > 0) %>% 
  .$gene

commonModel = bothModel %>% 
  filter(gene %in% commonGenes)

genes = unique(commonModel$gene)

for (i in 1:length(genes)) { #for each gene that is modeled by both Prax and GTEx
  geneModel = commonModel %>% filter(gene == genes[i])
  commonSnps = geneModel %>% .$rsid %>% table
  
  if (!any(commonSnps > 1)) {next}
  
  
  if (i == 1) {
    modelComparison = geneModel %>% filter(rsid %in% names(commonSnps[commonSnps > 1])) %>% spread(dataSource, weight)
  }else{
    modelComparison %<>% rbind( (geneModel %>% filter(rsid %in% names(commonSnps[commonSnps > 1])) %>% spread(dataSource, weight)))
  }
}

modelComparison %>% 
  ggplot(aes(GTEx, PRAX)) + 
  geom_point() + 
  ggtitle(paste0('PRAX vs GTEx whole blood PrediXcan models - snp weight \ncomparison of ', 
                 nrow(modelComparison), 
                 ' snps used in both models to predict expression of ', 
                 modelComparison %>% .$gene %>% unique %>% length, 
                 ' genes')) +
  xlab('GTEx model SNP weights')  +
  ylab('PRAX model SNP weights') + 
  annotate('text', x = -.75, y = -.75, label = paste0('R^2 = ', summary(lm(PRAX ~GTEx, data = modelComparison))$r.squared %>% format(digits = 4)))
ggsave('~/plateletExpressionModelling/outputs/plots/preliminaryPRAX_GTExComparison.png')

#### Same thing for predictDB's DGN model
setwd('~/plateletExpressionModelling/data/')
# system('wget https://s3.amazonaws.com/predictdb/DGN-HapMap-2015.tar.gz')
# system(paste('tar -zxvf', 'DGN-HapMap-2015.tar.gz'))

praxDB = dbConnect(drv = SQLite(), dbname = '~/plateletExpressionModelling/prediXcanModeller/data/output/dbs/firstPraxModelSignificanceFiltered.db')
dgnDB = dbConnect(SQLite(), dbname = '~/plateletExpressionModelling/data/DGN-HapMap-2015/DGN-WB_0.5.db')

praxModel = praxDB %>% 
  dbReadTable('weights') %>% 
  as.tbl() %>% 
  mutate(dataSource = 'PRAX')

dgnModel = dgnDB %>% 
  dbReadTable('weights') %>% 
  as.tbl %>% 
  mutate(dataSource = 'dgn',
         gene = str_extract(gene, '[A-Z0-9]+(?=\\.)')) #regex to remove transcript numbers

bothModel = praxModel %>% rbind(dgnModel)

commonGenes = bothModel %>% 
  group_by(gene) %>% 
  summarise(numSnpsPRAX = sum(dataSource == 'PRAX'),
            numSnpsdgn = sum(dataSource == 'dgn')) %>% 
  filter(numSnpsPRAX > 0, numSnpsdgn > 0) %>% 
  .$gene

commonModel = bothModel %>% 
  filter(gene %in% commonGenes)

genes = unique(commonModel$gene)

for (i in 1:length(genes)) { #for each gene that is modeled by both Prax and dgn
  geneModel = commonModel %>% filter(gene == genes[i])
  commonSnps = geneModel %>% .$rsid %>% table
  
  if (!any(commonSnps > 1)) {next}
  
  
  if (i == 1) {
    PRAX_DGN_modelComparison = geneModel %>% filter(rsid %in% names(commonSnps[commonSnps > 1])) %>% spread(dataSource, weight)
  }else{
    PRAX_DGN_modelComparison %<>% rbind( (geneModel %>% filter(rsid %in% names(commonSnps[commonSnps > 1])) %>% spread(dataSource, weight)))
  }
}

PRAX_DGN_modelComparison %>% 
  ggplot(aes(dgn, PRAX)) + 
  geom_point() + 
  ggtitle(paste0('PRAX vs DGN PrediXcan models - snp weight comparison of ', 
                 nrow(PRAX_DGN_modelComparison), ' snps \nused in both models to predict expression of ', 
                 PRAX_DGN_modelComparison %>% .$gene %>% unique %>% length, 
                 ' genes')) +
  xlab('DGN model SNP weights')  +
  ylab('PRAX model SNP weights') + 
  annotate('text', x = -.75, y = -.75, label = paste0('R^2 = ', summary(lm(PRAX ~ dgn, data = PRAX_DGN_modelComparison))$r.squared %>% format(digits = 4)))
ggsave('~/plateletExpressionModelling/outputs/plots/preliminaryPRAX_DGNComparison.png')

#### Compare DGN & GTEx as a standard for how well things should overlap across different datasets

bothModel = gtexModel %>% rbind(dgnModel)

commonGenes = bothModel %>% 
  group_by(gene) %>% 
  summarise(numSnpsGTEx = sum(dataSource == 'GTEx'),
            numSnpsdgn = sum(dataSource == 'dgn')) %>% 
  filter(numSnpsGTEx > 0, numSnpsdgn > 0) %>% 
  .$gene

commonModel = bothModel %>% 
  filter(gene %in% commonGenes)

genes = unique(commonModel$gene)

for (i in 1:length(genes)) { #for each gene that is modeled by both GTEx and dgn
  geneModel = commonModel %>% filter(gene == genes[i])
  commonSnps = geneModel %>% .$rsid %>% table
  
  if (!any(commonSnps > 1)) {next}
  
  
  if (i == 1 || !exists('GTEx_DGN_modelComparison')) {
    GTEx_DGN_modelComparison = geneModel %>% filter(rsid %in% names(commonSnps[commonSnps > 1])) %>% spread(dataSource, weight)
  }else{
    GTEx_DGN_modelComparison %<>% rbind( (geneModel %>% filter(rsid %in% names(commonSnps[commonSnps > 1])) %>% spread(dataSource, weight)))
  }
}

p = GTEx_DGN_modelComparison %>% 
  ggplot(aes(dgn, GTEx)) + 
  geom_point(color = rgb(0,0,0,.2), stroke = 0) + 
  ggtitle(paste0('GTEx vs DGN PrediXcan models - snp weight comparison of ', nrow(GTEx_DGN_modelComparison), ' snps \nused in both models to predict expression of ', 
                 GTEx_DGN_modelComparison %>% .$gene %>% unique %>% length, 
                 ' genes')) +
  xlab('DGN model SNP weights')  +
  ylab('GTEx model SNP weights') + 
  annotate('text', x = 1.5, y = -.75, label = paste0('R^2 = ', summary(lm(GTEx ~ dgn, data = GTEx_DGN_modelComparison))$r.squared %>% format(digits = 4)))
ggsave('~/plateletExpressionModelling/outputs/plots/preliminary_GTEx_DGN_Comparison.png', p) #ggMarginal(p)

save(list = c('modelComparison', 'GTEx_DGN_modelComparison', 'PRAX_DGN_modelComparison'), file = '~/plateletExpressionModelling/outputs/PRAX_DGN_GTEx_modelComparison.RData')


