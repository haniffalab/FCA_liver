#!/bin/bash

#$ -cwd
#$ -N prepare_input
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=400G

Rscript prepare_input.R

echo "End on `date`"
