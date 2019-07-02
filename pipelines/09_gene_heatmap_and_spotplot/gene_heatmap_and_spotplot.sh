#!/bin/bash

#$ -cwd
#$ -N gene_heatmap_and_spotplot
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=200G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript gene_heatmap_and_spotplot.R $1

echo "End on `date`"
