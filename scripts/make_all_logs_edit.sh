#!/bin/bash
# Script to concatenate log files split by chromosome together.

tissue=$1
allLogs=$2
i=0
foldnum=$3
for logfile in $(ls /mnt/labhome/andrew/plateletExpressionModelling/outputs/model_by_chr_fold/${foldnum}/${tissue}_chr*log.txt); do
        if [ $i -eq 0 ] ; then
                head -n 1 $logfile > $allLogs
                i=1
        fi
        tail -n +2 $logfile >> $allLogs
done
