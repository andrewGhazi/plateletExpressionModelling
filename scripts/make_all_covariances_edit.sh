#!/bin/bash
# Script to concatenate covariance fiels split by chromosome together

tissue=$1
allCovariances=$2
alpha=$3
snpset=$4
foldnum=$5
echo "GENE RSID1 RSID2 VALUE" > $allCovariances
for covfile in $(ls /mnt/labhome/andrew/plateletExpressionModelling/outputs/model_by_chr_fold/${foldnum}/${tissue}_chr*_${snpset}_alpha_${alpha}_covariances.txt); do
	cat $covfile >> $allCovariances
done
gzip $allCovariances

