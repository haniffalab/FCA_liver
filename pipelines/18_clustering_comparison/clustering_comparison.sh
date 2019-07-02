#!/bin/bash

#$ -cwd
#$ -N clustering_comparison
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=600G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript clustering_comparison.R $1

echo "End on `date`"
