#!/bin/bash

#$ -cwd
#$ -N split_seurat_by_category
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=400G

Rscript split_seurat_by_category.R

echo "End on `date`"
