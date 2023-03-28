# aspen heritability

Collab with Ben Blonder at UC Berkeley to calculate heritability of phenology traits in Colorado, USA aspen samples.
Performed 2022
by: Shannon Hateley, PhD

Some absolute paths to reference data may need to be changed as I moved the aspen analysis from my personal directory into a shared aspen project directory prior to leaving my postdoc. Other than that, things *should* work fine. A lot of the analysis (and the initial heritability results writeup) in this directory are from the initial dataset (data_for_heritability.csv) that did not include the additional variables of snowmelt, max temp, etc. The results that were actually published in the paper are from the dataset that includes these variables (data_for_heritability_updated.csv). There is some code related to GEMMA, but not many files related to GEMMA. This is because I ended up not using GEMMA in the analysis after playing around and considering my preferred method of assessing heritability and desiding to use r to regress out covariates and then GCTA to calculate heritability. First I thought to use GEMMA to just model everything together but I am more familiar with GCTA for heritability and trust its output, plus this allowed me to play around with covariates in a more interpretable manner. You can pretty much disregard the GEMMA stuff.
 
Basically the project consists of the following:


1) data
data_for_heritability_updated.csv
aspen250 filtered genotype calls subset to samples in the phenology study and keeping genotype calls from 1 sample of sequencing data per tree only. The resulting files are called aspen_phenology.X. The steps leading up to this are the plink1/2 files and the list of smaples to keep and remove. The aspen250 samples are from previous work in in Blonder/Moi labs with the radseq data.
the gcta subfolder just has some files generated from on the genotype call files that are used in analyses like the grm and eigen values.

2) code
OLD_aspen_heritability.R old r analysis file from initial analysis. This did not go into the final paper
aspen_heritability_updated_covariates.R the real r analysis file for the paper plus a few messy notes and extra analyses / qc stuff
misc.txt is scripts and stuff for running software and bash commands
greml files are the fun scripts for running gcta on all the subset samples

3) aspenDir
links back to the main aspen directory. was useful when I was working out of my personal directory

4) results
basically everything generated from either the old or new R scripts or gcta. Some are results plots and some are intermediate files like residuals output from R that feed into gcta. the up to date file in the residuals and stuff are the jan_2023 updated stuff and the ones not in the update file are the original analysis. A little messy but you get the idea. a lot of plots are qc or extra analyses stuff that didn't make it into the paper. The names should be pretty self explanatory but you can also search for them in the code to see what they are if you need details. I can't guarantee that all plots will have code but I can guarantee that all important and published plots have code. There is also an "albi" folder that has some plots looking at kinship - not super important.

5) main directory files
this readme
a markdown file presenting initial analysis in a pretty-ish way

software: 
I used GCTA version gcta_v1.94.0Beta_linux_kernel_3_x86_64

