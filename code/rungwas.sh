#!/bin/bash
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --partition=DPB
#SBATCH --job-name=aspen_phenol$1
#SBATCH --output=aspen_phenol$1.log

gemma -bfile ../data/aspen_phenology -n $1  -lmm 2 -k ../data/aspen_phenology.cXX.txt -c ../data/covariates.tsv -o  ../resutls/aspen_phenol$1

