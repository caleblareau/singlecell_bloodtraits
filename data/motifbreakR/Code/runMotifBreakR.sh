#!/bin/bash

#$ -cwd
#$ -q long
#$ -N mbreaker.allvariants
#$ -cwd
#$ -t 1-16
#$ -tc 100
#$ -l h_vmem=50g 

# Software
source /broad/software/scripts/useuse
reuse R-3.3

# Variables
filelist="/broad/sankaranlab/ebao/MotifbreakR/mbreakerfiles_PP0.001.txt"

# Get sample from file list
file=$(cat $filelist | awk -v ln=$SGE_TASK_ID "NR==ln")
trait=$(cat /broad/sankaranlab/ebao/FINEMAP/traits.txt | awk -v ln=$SGE_TASK_ID "NR==ln")

# file=$(cat $filelist | awk -v ln=17 "NR==ln")
# trait=$(cat /broad/sankaranlab/ebao/FINEMAP/traits.txt | awk -v ln=1 "NR==ln")


# Call Rscript to run motifbreakR
Rscript /broad/sankaranlab/ebao/MotifbreakR/motifbreakR.R $file $trait
#Rscript motifbreakR.R $file 

