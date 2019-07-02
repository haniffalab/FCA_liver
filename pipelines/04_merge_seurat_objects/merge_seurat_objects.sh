#!/bin/bash

#$ -cwd
#$ -N merge_seurat_objects
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=200G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript merge_seurat_objects.R $1

echo "End on `date`"
