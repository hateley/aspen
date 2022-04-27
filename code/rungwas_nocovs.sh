#!/bin/bash
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --partition=DPB
#SBATCH --job-name=aspen_phenol_nofilter$1
#SBATCH --output=aspen_phenol_nocovs.log

gemma -bfile ../data/aspen_phenology -miss 1.0 -maf 0.0 -n $1 -lmm 2 -k ../data/aspen_phenology.cXX.txt -o aspen_phenol_no_covs$1

