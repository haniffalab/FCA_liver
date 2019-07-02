#!/bin/bash

#$ -cwd
#$ -N compute_DEG
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=200G

Rscript compute_DEG.R 

echo "End on `date`"
