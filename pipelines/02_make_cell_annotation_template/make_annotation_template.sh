#!/bin/bash

#$ -cwd
#$ -N make_annotation_template
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=200G
#$ -pe smp 6

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript make_annotation_template.R $1

echo "End on `date`"
