#!/bin/bash

# Shannon Hateley, PhD
# December 2021
# submit all aspen phenology gwas

for i in `seq 1 16`;
do
 echo item: $i
 sbatch ./rungwas.sh $i
done


