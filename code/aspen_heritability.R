# Author: Dr. Shannon Hateley
# Date: Fall 2021
# Project: Aspen Phenology with Dr. Benjamin Blonder

###############################
# Heritability of traits
# 
# 1) remove confounders
#    phenology = f(elve, slope, ...) + e resid + cytotype
# 2) map residuals
#    e phenology resids = u + bx + Z(pop struc) + e
#

# phenology = u + bx + pcs + cov. + .. e 
#
#
# going to use all the covariates in GEMMA, run each phenotype separately, 
# also take averages of the mutliple year observations




###############################


# set up environment
library(tidyverse)
library(caret) #makes one-hot encoding easy
library(reshape2) #for melt
library(cowplot)
theme_set(theme_cowplot())

setwd("~/safedata/slhateley/aspen/")



###############################
# load data
phenology <- read_csv("data/data_for_heritability.csv")
names <- colnames(phenology)

site_code <- phenology[1]
vars <- phenology[-1] #data to one-hot encode
tmp <- dummyVars("~.", data=vars, fullRank=F) #set fullRank=F for the covariate analaysis; later fullRank=T for outputting covariates, since we want to avoid fully correlated covariates 
tmpdf <- data.frame(predict(tmp, newdata = vars)) #put results from one-hot encoding into dataframe again
pheno_coded <- cbind(site_code, tmpdf) # add back the ids

rm(tmp, tmpdf)


###############################
# data checks


#######
# check correlations between phenology variables
# use fullRank=F on cavariates
#######


var_cor <- cor(pheno_coded[-1], use="pairwise.complete.obs")

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
  
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
var_cor <- reorder_cormat(var_cor)
upper_tri <- get_upper_tri(var_cor)
# Melt the correlation matrix
melted_var_cor <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_var_cor, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)
ggsave("results/phenology_correlation.pdf", ggheatmap)

#check for normality of phenology
phenology.continuous <- phenology %>%
  select(Site_Code, x, y, Elevation, Cos.aspect, Slope, Summer.Insolation, 
         DBH.mean, fraction_aspen, names[grep("pheno", names)])

shapiro <- sapply(phenology.continuous[-1], shapiro.test) %>% as.tibble() # shapiro wilkes nomality test
res1 <- shapiro[1,] %>%
  pivot_longer(colnames(shapiro), names_to = "variable", values_to = "stat")
res2 <- shapiro[2,] %>%
  pivot_longer(colnames(shapiro), names_to = "variable", values_to = "pval")
res_shapiro <- cbind(res1,res2[,2])
res_shapiro$stat <- res_shapiro$stat %>% gsub("^c.*= ","",.) %>% gsub(")","",.)
res_shapiro$pval <- unlist(res_shapiro$pval)
write_tsv(res_shapiro, "results/shapiro_wilk.tsv")

pheno.gathered.continuous <- phenology.continuous %>%
  gather(key = "variable", value = "value",
         -Site_Code)

ggplot(pheno.gathered.continuous, aes(x = value)) +
  geom_histogram() +
  facet_wrap(~variable, scales = 'free')
ggsave("results/phenology_histogram.pdf")

### Check for outliers ###

phenology.continuous[sapply(phenology.continuous, is.null)] <- NA
pca <- prcomp(na.omit(phenology.continuous[-1]), scale. = TRUE, rank. =10)
U <- pca$x
pc_plot <- qplot(U[, 1], U[, 2], xlab = 'PC 1', ylab = "PC 2", main = "aspen phenology PCA") + coord_equal()
print(pc_plot)
ggsave("results/pc1pc2.pdf", pc_plot)
apply(U, 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) )) #> 6sd from mean
# integer(0) no outliers are 6sd away from mean (so not that bad, but tis is somewhat arbitrary)

###############################
# Prep files for LMM input
#############

# change one-hot encoding to fullRank=T to remove fully-correlated columns
tmp <- dummyVars("~.", data=vars, fullRank=T) # fullRank=T for outputting covariates, since we want to avoid fully correlated covariates 
tmpdf <- data.frame(predict(tmp, newdata = vars)) #put results from one-hot encoding into dataframe again
pheno_coded <- cbind(site_code, tmpdf) # add back the ids

# take averages of the annual data and add to list of phenotypes

pheno_coded$pheno.mean.OGI <- rowMeans(pheno_coded[, names[grep("OGI", names)]], na.rm=T)
pheno_coded$pheno.mean.OGMn <- rowMeans(pheno_coded[, names[grep("OGMn", names)]], na.rm=T)
pheno_coded$pheno.mean.EVImax <- rowMeans(pheno_coded[, names[grep("EVImax", names)]], na.rm=T)
pheno_coded$pheno.mean.GSL <- rowMeans(pheno_coded[, names[grep("GSL", names)]], na.rm=T)
pheno_coded$intercept <- 1 # needed for GEMMA covariate file format

write_tsv(as.tibble(pheno_coded$Site_Code), "data/aspen.ids", col_names = FALSE) # output initial set of sample ids

# subset to samples with genotypes
geno_fam <- read_table2("data/plink_cleanup_2.fam", col_names = F)
pheno_coded <- pheno_coded %>% filter(Site_Code %in% geno_fam$X2)


#############
# output environmental variables to run in GEMMA LMM
#############

write_tsv(select(pheno_coded, names[grep("pheno", names)]), "data/phenotypes.tsv", col_names = F) # output phenology phenos
write_tsv(select(pheno_coded, names[grep("pheno", names)])[0,], "data/phenotypes.labels")
write_tsv(select(pheno_coded, intercept, colnames(pheno_coded[2:17])), "data/covariates.tsv", col_names = F) # output covariates
write_tsv(pheno_coded[0, c(39,2:17)], "data/covariates.labels") # output covariate names


'''
notes

We will need:
  the genomic data, which we have already:
    aspen250.bed/bim/fam/cXX.txt
    aspen250 is better quality but less snps than aspen10
    (really this stuff should be re-called now that there is a reference)
    also there are replicates of the trees so we could pool that sequence and call together for better coverage.
  the matched phenology value for each of the trees (503 trees)
   data_for_heritability.csv
  any potential confounding factor that could affect phenology: 
    slope, altitude, soil, etc. (this will capture variation in the environment across trees so we can account for it).
    also data_for_heritability.csv
    

___________________

full draft of the phenology manuscript:
  https://docs.google.com/document/d/1-fTVokEaYnh7hoWIJSzfgbuHKWDSDc6mB8tqwCGzvi4/edit#

Would appreciate any input by the 1st week of November if it is feasible.

___________________

Hi Moi and Shannon,

Great!! I attach the information you need. Note we have 4 years of phenology data so there are multiple entries for each phenology 
variable. I am guessing we would need to treat this as phenotypic plasticity, maybe by including a random intercept for year in the
model you are going to run. I also attach a large set of covariates that we are including in the main text models as well. 
The phenological variables of interest are: OGI (greenup date), OGMn (greendown date), GSL (growing season length), EVImax (peak 
greenness). 
I only included spatial covariates that are time invariant: 
"Summer.Insolation" "DBH.mean" "Canopy_openness" "Soil.type" "Cytotype" "geneticSexID" "Rock_Unit". 
We could include other time dependent covariates too (like snowmelt date), like in the main text models, but I am not sure if this 
is a good idea for a heritability analysis or not.
I also included info on cytotype and sex; not sure if these should be included, or if the model should be run separately for each
cytotype.
I also included a column "fraction_aspen" which indicates the % of each Sentinel-2 pixel (from which the phenology data are derived) 
is aspen forest. I think you should repeat the analysis using 100% of the data and then also thresholding only fraction_aspen>0.5 to 
ensure the phenology data is of high quality in the analysis.

Let me know if you have any questions.

Ben


___________________

Statistical analyses – Question 2
We extracted phenological data over each plot for which we had RADseq genomic data available. Values for greenup date, 
greendown date, and growing season length were paired when phenology algorithm detected an annual cycle and potentially
inter-annual gap-filling was used to resolve groups (i.e. a more liberal criterion than the above spatial analysis). 
Values for maximum grenness were always included. For all the above data, values were removed from the analysis when remotely
sensed aspen cover at 1 m resolution within a 30 m HLS pixel above the plot was less than 50%. A total of XX/503 plots were 
retained for the analysis. 

To ensure that variance was appropriately assigned to genetic sources, we also included several spatially variable but 
temporally constant variables for each plot, all of which were available from ground-based surveys reported in (Blonder et 
al. 2021). These variables included elevation, slope, cosine aspect, summer insolation, mean diameter of breast height of trees
in the plot, canopy openness, regolith type, and rock unit type. Information on cytotype and sex was also included.

We then carried out a genome-wide association (GWA) analysis for each of the phenological statistics. We developed linear mixed
models: y = ßX + Q + u + ε (Yu et al. 2006), in which we map the phenotype (y) variation onto genomic variation (X, coded as 
0, 1, 2 for each allele) and environmental variation (Q). Given the high-dimensionality of the genetic variance, we control 
for co-variation between predictors (linkage disequilibrium and population structure) via a random individual factor informed 
by the K kinship matrix across individuals: u ~ Normal(µ=0, σ=KVa); where Va represents the overall variance in the trait 
explained by genetics. The remaining variance in the trait is captured by the error ε~Normal(µ=0, σ=Ve). We then calculate
the broad-sense heritability (h2) of each phenology metric as the variance in y explained by variance in X.

##################

I also had to remove the replicates from the dataset, and used the aspen10 fam file to do this, since it still had the correct names.
If I recall genotypes using the reference genome now available, I should also combine sequence from these replicates so we have more data.

Also not all of these samples are ones that we have genotype data for. Only 488 samples:
Of the 489 samples with genotypes, 1 does not have phenology data: RGBXO
Of the 503 samples with phenology data, 14 do not have geno data


