# aspen heritability

Collab with Ben Blonder at UC Berkeley to calculate heritability of phenology traits in Colorado, USA aspen samples.
Performed 2022
by: Shannon Hateley, PhD

Some absolute paths to reference data may need to be changed as I moved the aspen analysis from my personal directory into a shared aspen project directory prior to leaving my postdoc. Other than that, things *should* work fine. A lot of the analysis (and the initial heritability results writeup) in this directory are from the initial dataset (data_for_heritability.csv) that did not include the additional variables of snowmelt, max temp, etc. The results that were actually published in the paper are from the dataset that includes these variables (data_for_heritability_updated.csv). There is some code related to GEMMA, but no files related to GEMMA. This is because I ended up not using GEMMA in the analysis after playing around and considering my preferred method of assessing heritability and desiding to use r to regress out covariates and then GCTA to calculate heritability. First I thought to use GEMMA to just model everything together but I am more familiar with GCTA for heritability and trust its output, plus this allowed me to play around with covariates in a more interpretable manner. You can disregard the GEMMA stuff.
 
Basically the project consists of the following:
1) Data




