#!/bin/bash

#$ -cwd
#$ -N write_data
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=400G

Rscript write_data.R

echo "End on `date`"
