#!/bin/bash

#$ -cwd
#$ -N update_annotation
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=200G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript update_annotation.R $1

echo "End on `date`"
