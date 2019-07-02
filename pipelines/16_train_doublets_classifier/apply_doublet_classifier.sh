#!/bin/bash

#$ -cwd
#$ -N apply_doublets_SVM
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=100G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript apply_doublet_classifier.R $1

echo "End on `date`"
