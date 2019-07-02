#!/bin/bash

#$ -cwd
#$ -N train_doublets_SVM
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=300G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript train_doublets_SVM.R $1

echo "End on `date`"
