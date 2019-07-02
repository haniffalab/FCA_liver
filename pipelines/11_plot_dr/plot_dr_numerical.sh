#!/bin/bash

#$ -cwd
#$ -N plot_dr
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=200G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript plot_dr_numerical.R $1 

echo "End on `date`"
