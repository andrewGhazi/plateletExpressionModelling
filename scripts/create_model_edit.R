#This is Andrew Ghazi's edited version of create_model.R that's changed to avoid some of the inbuilt assumptions about file names / locations
#Also they wrote %&% as their own version of paste0 in the file the source, lol

argv <- commandArgs(trailingOnly = TRUE)
source("~/plateletExpressionModelling/scripts/GTEx_Tissue_Wide_CV_elasticNet_edit.R")

tis <- argv[1]
chrom <- as.numeric(argv[9])
alpha <- as.numeric(argv[7])
window <- as.numeric(argv[11])


data_dir <- "~/plateletExpressionModelling/prediXcanModeller/data/intermediate/"

#expression_RDS <- data_dir %&% "expression_phenotypes/" %&% tis %&% "_Analysis.expr.RDS"
expression_RDS <- argv[2]
#geno_file <- data_dir %&% "genotypes/" %&% tis %&% "_chr" %&% chrom %&% "_Analysis.snps.biallelic.txt"
geno_file <- argv[3]
#gene_annot_RDS <- data_dir %&% "annotations/gene_annotation/gencode.v19.genes.v6p.patched_contigs.parsed.RDS"
gene_annot_RDS <- argv[4]

#snp_annot_RDS <- data_dir %&% "annotations/snp_annotation/GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.chr" %&% chrom %&% ".RDS"
snp_annot_RDS <- argv[5]
n_k_folds <- 10
#out_dir <- data_dir %&% "model_by_chr/"
out_dir <- argv[8]
snpset <- argv[10] # "1KG_snps"

TW_CV_model(expression_RDS, geno_file, gene_annot_RDS, snp_annot_RDS,
    n_k_folds, alpha, out_dir, tis, chrom, snpset, window)
