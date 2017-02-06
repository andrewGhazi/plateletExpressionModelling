#Let's start trying to get everything prepared for the platelet PrediXcan model building
#I'm following along from here: https://github.com/hakyimlab/PredictDBPipeline/wiki/Detailed_Description

# We need the following inputs:
# The start and end position of the gene
# All snps located within 1Mb base pairs of that gene (1Mb upstream of trascription start site and 1Mb downstream of transcription termination site).
# Expression quantification data for the gene in the specific tissue for the GTEx individuals
# Genotypes for all these individuals on the snp set

#Except we're going to be using our cohort instead of the GTEx data, hopefully that isn't a huge pain in the ass

##### Libraries
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(data.table)
library(stringr)
library(dtplyr)
library(RSQLite)
library(biomaRt); select = dplyr::select; contains = dplyr::contains

setwd("~/plateletExpressionModelling")

##### SNP annotation file - A tab delimited file, with a header row, containing the fields: chromosome, position, VariantID, RefAllele, AlternativeAllele, original RSID number, more recent RSID number, Num_alt_per_site.
# At a minimum, there must be a header row, with chromosome in the first column, position in the second, VariantID in the third (form: chr_pos_refAllele_altAllele_build), ref Allele in the fourth, alt Allele in the fifth, and rsid in the seventh.
snpAnnotations = read_tsv('/mnt/labhome/simonlm/projects/PRAX/Papers/eQTLpaper/MatrixEQTL/data/snploc.txt')

#Let's convert the kgp snps to rs identifiers using this method http://rstudio-pubs-static.s3.amazonaws.com/9686_5d78b74c4b9743cfb5ac63aba8d309b5.html
#kgpAnnotations %>% mutate(pos1 = pos + 1) %>% select(chr, pos, pos1, snp) %>% write_tsv('inputs/kgpToRsUCSCTableInput.txt', col_names = FALSE)
setwd('data/')
#system('wget http://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanomni5-4/v1-2/infinium-omni5-4-v1-2-a1-b144-rsids.zip',
#       ignore.stdout = TRUE)
#system('unzip infinium-omni5-4-v1-2-a1-b144-rsids.zip')
kgpToRs = read_delim('InfiniumOmni5-4v1-2_A1_b144_rsids.txt',
                     col_names = TRUE,
                     delim = '\t')
names(kgpToRs)[1] = 'snp'
snpAnnotations %<>% inner_join(kgpToRs, by = 'snp') %>% 
  filter(!grepl(',', RsID))

#Now we need to get the ref/alt alleles
#mutate(VariantID = paste0(chr, pos, refAllele, altAllele, build, sep = '_'))

write_tsv(snpAnnotations, 
          'data/snpAnnotationIntermediate.txt',
          col_names = TRUE)
save(snpAnnotations, 
     file = 'data/snpAnnotationIntermediate.RData')
write_tsv(snpAnnotations %>% select(RsID),
          path = 'data/snpRsIds.txt',
          col_names = FALSE)

# 
# for(i in 1:22){
#   # chrdat = snpAnnotations %>% filter(chr == paste0('chr', i))
#   # chrdat %>% mutate(rsIDwithTab = paste0(RsID, '\t')) %>% 
#   #   select(rsIDwithTab) %>% 
#   #   write_tsv(path = paste0('data/rsIdsByChromosome/chr', i, '_rsIds.txt'))
#   dbSNPdati = read_tsv(paste0('/mnt/bigData2/resources/dbSNP_146_GRCh37p13/byChromosome/', i, '.vcf'),
#                        col_names = FALSE)
#   names(dbSNPdati) = c('CHROM', 
#                        'POS', 
#                        'RsID',
#                        'REF', 
#                        'ALT',
#                        'QUAL', 
#                        'FILTER',
#                        'INFO')
#   
#   #TERRIBLE CODE INCOMING!
#   if (i == 1){
#     snpAnnotations %<>% left_join(dbSNPdati, by = 'RsID')
#   } else{
#     dbSNPdati %<>% filter(RsID %in% snpAnnotations$RsID)
#     snpAnnotations[match(dbSNPdati$RsID, snpAnnotations$RsID), c('CHROM', 
#                                                                  'POS', 
#                                                                  'RsID',
#                                                                  'REF', 
#                                                                  'ALT',
#                                                                  'QUAL', 
#                                                                  'FILTER',
#                                                                  'INFO')] = dbSNPdati
#   }
# }
#Input file written by scripts/snpAnnotationAssembly.R

# snpAnnotations %>% 
#   mutate(VariantID = str_c(CHROM, POS, REF, ALT, 'hg19', sep = '_')) %>% 
#   select(CHROM, POS, VariantID, REF, ALT, RsID)

##### genotype file - tab delimited, has a header row, columns are variantID then individual IDs. Value = dosage of second allele listed in variantID
# A tab delimited file, with a header row, first column being the variantID of the SNP, and all remaining columns being individual ids. The values indicate the dosage for the second allele listed in the variantID.
kgpToRs %<>% rename(id = snp)

reverseDosage = function(dosage){
  2 - dosage
}
snpAnnotations = read_tsv('snpAnnotationFile.txt')
setwd("~/plateletExpressionModelling")

genotypes = fread('/mnt/labhome/simonlm/projects/PRAX/Papers/eQTLpaper/MatrixEQTL/data/genotype.txt',
                  header = TRUE) %>% 
  as.tbl %>% 
  left_join(kgpToRs, by = 'id', copy = TRUE) %>% #I don't exactly know what the copy = TRUE argument does but it throws and error otherwise
  left_join(snpAnnotations, by = 'RsID', copy = TRUE) %>% 
  na.omit %>% 
  select(VariantID, X001.1:X163.189) %>% #TODO do 2 - values because it needs to be the dosage of the minor allele
  mutate_each(funs(reverseDosage), vars = -VariantID) %>% #We reverse the dosage because PrediXcan is expecting dosage of the minor allele
  write_tsv('data/genotypeFile.tsv', col_names = TRUE)

# Compile genotype files for running PrediXcan
genotypes %>% select(contains('X')) %>% names %>% data_frame(FID = ., IID = .) %>% write_tsv('~/plateletExpressionModelling/data/genotypeByChr/samples.txt',
                                                                                             col_names = FALSE)

customFun  = function(DF) { #Thanks http://stackoverflow.com/questions/41233173/how-can-i-write-dplyr-groups-to-separate-files
  write_tsv(DF %>% select(-CHROM), paste0("~/plateletExpressionModelling/data/genotypeByChr/praxGenotype.chr", unique(DF$CHROM), ".tsv"))
  return(DF)
}

genotypes = fread('/mnt/labhome/simonlm/projects/PRAX/Papers/eQTLpaper/MatrixEQTL/data/genotype.txt',
                  header = TRUE) %>% 
  as.tbl %>% 
  left_join(kgpToRs, by = 'id', copy = TRUE) %>% #I don't exactly know what the copy = TRUE argument does but it throws and error otherwise
  left_join(snpAnnotations, by = 'RsID', copy = TRUE) %>% 
  na.omit %>% 
  mutate_each(funs(reverseDosage), vars = contains('X')) %>% 
  mutate(MAF = .1) %>% # Just faking this for now
  select(CHROM, RsID, POS, REF, ALT, MAF, contains('X')) %>% 
  filter(RsID %in% gtexSnps)

setwd('~/plateletExpressionModelling/data/genotypeByChr/')
totSize = 0
for (i in 1:22) {
  tmp = genotypes %>% filter(CHROM == i)
  write_tsv(tmp, paste0("~/plateletExpressionModelling/data/genotypeByChr/praxGenotype.chr", unique(tmp$CHROM), ".tsv"),
            col_names = FALSE)
  totSize = object.size(tmp) + totSize
  system(paste0('gzip ',
                'praxGenotype.chr', i, '.tsv'))
  cat(i)
}


# andrew@smirnov:~/plateletExpressionModelling/prediXcanModeller/scripts$ python split_genotype_by_chr.py ~/plateletExpressionModelling/data/genotypeFile.tsv ~/plateletExpressionModelling/data/genotypeByChr/

#genotypes$VariantID = snpAnnotations$VariantID[match()]

##### expression files - tab delimited, header row, first column = ensembl ID, remaining = individual IDs. Values = expression levels
# A tab delimited file, with a header row, first column being the ensembl ID of the gene, and all remaining columns being individual ids. The values indicate the expresion levels of the gene.

exprData = fread('/mnt/labhome/simonlm/projects/PRAX/Papers/eQTLpaper/MatrixEQTL/data/ExprCorrected.txt',
                 skip = 1,
                 header = FALSE) %>% as.tbl
names(exprData) = c('hgnc_symbol', system('head -1 /mnt/labhome/simonlm/projects/PRAX/Papers/eQTLpaper/MatrixEQTL/data/ExprCorrected.txt', intern = TRUE) %>% strsplit('\t') %>% unlist)

#get the ensembl ids necessary
symbToID = getBM(c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'hgnc_symbol', 
                 values = exprData$geneName, 
                 useMart('ensembl', dataset = 'hsapiens_gene_ensembl')) %>% 
  as.tbl

exprData %<>% left_join(symbToID, by = 'hgnc_symbol', copy = TRUE) %>% 
  rename(ensembl_id = ensembl_gene_id) %>% 
  select(ensembl_id, X001.1:X163.189) %>% 
  na.omit %T>% 
  write_tsv('data/expressionFile.tsv')

##### Gene annotation file, GTF format with gene_id (ensembl ID) and gene_name attributes
#I pulled this from ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
setwd('data')

#system('wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz')
#system('gunzip Homo_sapiens.GRCh37.75.gtf.gz')

extractGeneID = function(x){
  x %>% 
    str_split(';') %>% 
    unlist %>% 
    grep('gene_id', ., value = TRUE) %>% 
    gsub('gene_id ', '', .) %>% 
    gsub('\"', '', .)
}



system('python scripts/parse_gtf_edit.py ~/plateletExpressionModelling/data/Homo_sapiens.GRCh37.75.gtf data/intermediate/annotations/gene_annotation/praxGeneAnno2.parsed.txt')
#andrew@smirnov:~/plateletExpressionModelling/prediXcanModeller/scripts$ Rscript geno_annot_to_RDS.R ../data/intermediate/annotations/gene_annotation/praxGeneAnno2.parsed.txt ../data/intermediate/annotations/gene_annotation/praxGeneAnno2.RDS



# tail -n +2 expressionFile.tsv | cut -f 1 > praxEnsemblIds.txt
# grep -f praxEnsemblIds.txt geneAnnotationIntermediate.tsv > praxGenes.gtf #I manually added 'seqname' into the pattern file to get the header

#### Covariate data for the expression
# A tab delimited text file,

#Changes I had to make to their scripts :/
# 1. Edit how type is handled in their gtf parser
# 2. line 51 in split_snp_annot_by_chr was indexed improperly to element 6 instead of 5
# 3. Not an edit but in the second part of snp annotation preprocessing you need to use Rscript and put 'snp_annot.chr' as the suffix of the input 

#python create_meta_data.py --geno=../data/input/genotypes/genotypeFile.tsv --expr=../data/input/expression_phenotypes/expressionFile.tsv --out_prefix=../data/output/allMetaData/firstAttempt --snpset=PRAXOmni5M --rsid_label=RsID --window=1000000

# fullGTF = fread('data/Homo_sapiens.GRCh37.75.gtf', 
#                 skip = 5, 
#                 header = FALSE,
#                 col.names = c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'),
#                 sep = '\t') %>% 
#   as.tbl %>% 
#   filter(grepl('gene_id', attribute), 
#          grepl('gene_name', attribute)) %>% 
#   filter(feature == 'gene') %T>% 
#   write.table('data/geneAnnotationIntermediate.gtf',
#               quote = FALSE,
#               sep = '\t',
#               row.names = FALSE)

#### Phenotype file - used for PrediXcan association detection
phenotypes = read_tsv('/mnt/bigData2/andrew/PRAX/data/phenotypes.txt', skip = 1,
                      col_names = FALSE)
phenNames = c('phenotype', genotypes %>% head %>% select(contains('X')) %>% names)
names(phenotypes) = phenNames

phenotypes %<>% 
  gather(IID, score, -phenotype) %>% 
  spread(phenotype, score) %>% 
  mutate(FID = IID) %>% 
  select(FID, IID, ADP:vWF) %T>% 
  write_tsv('~/plateletExpressionModelling/data/phenotypeFile.tsv')
  
#### Gene list - list of chromosome / gene pairs that PrediXcan needs apparently even though it says it's optional
setwd('~/plateletExpressionModelling/data/')
gtexWholeBloodDB = dbConnect(SQLite(), dbname = '~/plateletExpressionModelling/data/GTEx-V6p-HapMap-2016-09-08/TW_Whole_Blood_0.5.db')

genes = gtexWholeBloodDB %>% 
  dbReadTable('weights') %>%
  as.tbl %>% 
  select(gene) %>% 
  mutate(gene = str_extract(gene, '[A-Z0-9]+(?=\\.)')) %>% 
  .$gene %>% 
  unique

mt = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
g = getBM(attributes = c('chromosome_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', values = genes, mart = mt)
write_tsv(g, '~/plateletExpressionModelling/data/geneList.tsv')

