#!/bin/bash
# Script to concatenate result files split by chromosome together.

# Andrew Ghazi's edit to change the ls of the directory 

tissue=$1
allResults=$2
alpha=$3
snpset=$4
foldnum=$5
i=0
for resultsfile in $(ls /mnt/labhome/andrew/plateletExpressionModelling/outputs/model_by_chr_fold/${foldnum}/working_TW_${tissue}_exp_10-foldCV_elasticNet_alpha${alpha}_${snpset}_chr*); do
        if [ $i -eq 0 ] ; then
                head -n 1 $resultsfile > $allResults
                i=1
        fi
        tail -n +2 $resultsfile >> $allResults
done
