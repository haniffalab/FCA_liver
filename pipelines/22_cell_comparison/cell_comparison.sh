#!/bin/bash

#$ -cwd
#$ -N cell_comparison
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=400G
#$ -pe smp 3

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript cell_comparison.R $1

echo "End on `date`"
