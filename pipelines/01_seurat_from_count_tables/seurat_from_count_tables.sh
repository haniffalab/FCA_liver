#!/bin/bash

#$ -cwd
#$ -N seurat_from_count_tables
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=300G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript seurat_from_count_tables.R $1

echo "End on `date`"
