#!/bin/bash

# Shannon Hateley, PhD
# December 2021
# submit all aspen phenology gwas

for i in `seq 1 20`;
do
 echo item: $i
 sbatch ./rungwas_nocovs.sh $i
done


