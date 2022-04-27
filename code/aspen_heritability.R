# Author: Dr. Shannon Hateley
# Date: Fall 2021 - Spring 2022
# Project: Aspen Phenology with Dr. Benjamin Blonder

###############################
# Heritability of traits
# 
# 1) remove confounders
#    phenology = f(elve, slope, ...) + e resid + cytotype (or model cytotype separately)
# 2) map residuals
#    e phenology resids = u + bx + Z(pop struc) + e (or use kin.blup)
#
# phenology = u + bx + pcs + cov. + .. e 
#
# GEMMA alternative method:
# use all the covariates in GEMMA, run each phenotype separately, 
# also take averages of the multiple year observations

# sets to run heritability on:
# diploids, triploids, all
# all aspen coverage, greater than 20%


###############################


# set up environment
library(tidyverse)
library(caret) #makes one-hot encoding easy
library(reshape2) #for melt
library(stringr)
library(cowplot)
library(ggpubr)
library(grid)
theme_set(theme_cowplot())

setwd("~/safedata/slhateley/aspen/")


##############################################################
# load data
##############################################################

phenology <- read_csv("data/data_for_heritability.csv")
names <- colnames(phenology)

site_code <- phenology[1]
vars <- phenology[-1] #data to one-hot encode
pheno_coded <- cbind(site_code, data.frame(predict(
                    dummyVars("~.", data=vars, fullRank=F), newdata = vars))
                    #put results from one-hot encoding into dataframe
                    #set fullRank=F for the covariate analaysis;
                    #later fullRank=T for outputting covariates, since we want to avoid fully correlated covariates
                     )


##############################################################
# data checks and exploration
##############################################################

####### check correlations between phenology variables #######
# use fullRank=F on covariates

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
ggheatmap <- ggplot(melted_var_cor, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + # minimal theme
  theme(axis.text.x = element_text(angle = 60, vjust = 1, 
                                   size = 6, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  coord_fixed()

# Print the heatmap
print(ggheatmap)
ggsave("results/phenology_correlation.png", ggheatmap)
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
  facet_wrap(~variable, scales = 'free') +
  theme(strip.text.x = element_text(size=6), axis.text=element_text(size=6))
ggsave("results/phenology_histogram.png")

### Check for outliers ###

phenology.continuous[sapply(phenology.continuous, is.null)] <- NA
pca <- prcomp(na.omit(phenology.continuous[-1]), scale. = TRUE, rank. =10)
U <- pca$x
pc_plot <- qplot(U[, 1], U[, 2], xlab = 'PC 1', ylab = "PC 2", main = "aspen phenology PCA") + coord_equal()
print(pc_plot)
ggsave("results/phenologies_pc1pc2.png", pc_plot)
apply(U, 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) )) #> 6sd from mean
# integer(0) no outliers are 6sd away from mean (so not that bad, but tis is somewhat arbitrary)

###########
# some more exploration

#pheno_vals <- apply(select(pheno_coded, all_of(phenos)), 2, function(x) {length(which(!is.na(x)))})
#na_covariates <- apply(pheno_coded, 2, function(x) {length(which(is.na(x)))})

#look at relatedness matrix pca
kinship <- read_tsv("data/GEMMA_files/aspen_phenology.cXX.txt", col_names = F)

kin_pca <- prcomp(kinship, scale. = TRUE, rank. =10)
U <- kin_pca$x
Udf <- tibble(PC1=U[,1],PC2=U[,2])
pc_plot <- qplot(U[, 1], U[, 2], xlab = 'PC 1', ylab = "PC 2", main = "aspen kinship PCA") + coord_equal()
print(pc_plot)

ggplot(data=Udf, aes(x=PC1, y=PC2)) +
  geom_point() +
  labs(x="PC 1", y="PC 2")

ggsave("results/kinship_pc1pc2.png", pc_plot)
summ <- summary(kin_pca)
summ # output percent variance explained. 1 and 2 are each ~10%
#this looks a lot like the genetic distance by geographic distance plots that Moi was getting during the hackathon




#plot xy coordinates
ggplot(data=phenology, aes(x=x, y=y)) +
  geom_point()

# look at relationship between predictors and phenologies

plot_cov2phenol_relationship <- function(long_df){
  df_list = list()
  for (p in c("GSL", "OGI", "OGMn", "EVImax")){
    print(p)
    df <- long_df %>% select(c(1:18, p))
    df.gathered <- df %>%
      as_tibble() %>%
      gather(key = "variable", value = "value",
             -Site_Code, -p)
    df_list[[p]] <- df.gathered
  }
  return(df_list)
}

plot_dfs <- plot_cov2phenol_relationship(long_coded_df)


plot <- ggplot(plot_dfs[["GSL"]], aes(x = value, y = GSL)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~variable, scales = 'free') +
  geom_smooth(method = "lm", se = FALSE, color='red') +
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label", size=3) +
  xlab("covariate relation to GSL")
  
plot
ggsave("results/cov_phen_scatterplots/GSL.pdf", plot, height = 6, width=10)

plot <- ggplot(plot_dfs[["OGI"]], aes(x = value, y = OGI)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~variable, scales = 'free') +
  geom_smooth(method = "lm", se = FALSE, color='red') +
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label", size=3) +
  xlab("covariate relation to OGI")

plot
ggsave("results/cov_phen_scatterplots/OGI.pdf", plot, height = 6, width=10)

plot <- ggplot(plot_dfs[["OGMn"]], aes(x = value, y = OGMn)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~variable, scales = 'free') +
  geom_smooth(method = "lm", se = FALSE, color='red') +
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label", size=3) +
  xlab("covariate relation to OGMn")
plot
ggsave("results/cov_phen_scatterplots/OGMn.pdf", plot, height = 6, width=10)

plot <- ggplot(plot_dfs[["EVImax"]], aes(x = value, y = EVImax)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~variable, scales = 'free') +
  geom_smooth(method = "lm", se = FALSE, color='red') +
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label", size =3) +
  xlab("covariate relation to EVImax")
plot
ggsave("results/cov_phen_scatterplots/EVImax.pdf", plot, height = 6, width=10)

##############################################################
# prepping fastq files for grenepipe 
# because I wanted to map to the genome instead of the rad-seq scaffolds
# but this was taking forever and I don't think it will help much so skipping this step
# also using this to pull out ploidy info
##############################################################

#phenology <- phenology %>% filter(Site_Code %in% geno_fam$X2)
#grenepipe_samples <- read_tsv("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/aspen/grenepipe/aspen/samples.tsv")
#gp_df <- grenepipe_samples %>%
#mutate(Site_Code = str_split(sample, "_", simplify = T)[,1])

ploidy_est <- read_csv("aspenDir/estploidy_pp98.csv", skip=1, 
                       col_names= c("var1", "id", "trip", "dip", "est_ploidy"))
ploidy_est <- ploidy_est %>% mutate(Site_Code = str_split(id, "_", simplify = T)[,1])

ploidy_est <- ploidy_est %>% group_by(Site_Code)  %>% 
  mutate(est_mismatch = if(n_distinct(est_ploidy)==1) 'F' else 'T')


tmp <- ploidy_est %>% group_by(Site_Code) %>%
  summarize(n_distinct(id))

#questionable ploidy = GAULT, ~GJBOV, JCSAT, KLVF

ploidy <- phenology %>% select(Site_Code, Cytotype) %>%
  full_join(., ploidy_est[-1], by = "Site_Code")

ploidy <- ploidy %>% mutate(ploidy_group = case_when(
  is.na(Cytotype) & est_mismatch == 'F' ~ est_ploidy,
  is.na(Cytotype) & est_mismatch == 'T' ~ 'Ambiguous',
  Cytotype == est_ploidy ~ Cytotype,
  est_ploidy == "Ambiguous" ~ Cytotype,
  TRUE ~ Cytotype))

#make a diploid and a triploid list

diploid_samples <- grenepipe_samples %>% 
  filter(sample %in% filter(ploidy, ploidy_group=="Diploid")$id) %>% arrange(sample)
write_tsv(diploid_samples, "variant_calling/diploids/samples.diploid.tsv")

triploid_samples <- grenepipe_samples %>% 
  filter(sample %in% filter(ploidy, ploidy_group!="Diploid")$id) %>% arrange(sample)
write_tsv(triploid_samples, "variant_calling/triploids/samples.trip_or_ambiguous.tsv")



##############################################################
# linear regression of the covariates
##############################################################

# make sure one-hot encoding removes a level (full rank)
pheno_coded <- cbind(site_code, data.frame(predict(
  dummyVars("~.", data=vars, fullRank=T), newdata = vars))
  )

#fill in missing sex with mean (0.5)
pheno_coded  <- pheno_coded %>% mutate(geneticSexIDM = ifelse(is.na(geneticSexIDM), 0.5, geneticSexIDM))


#pick the covariates I want to regress out
covariates <- c("Elevation", "Cos.aspect", "Slope", "Summer.Insolation", "DBH.mean", "Canopy_openness", "Soil.typeSoil", "Soil.typeTalus", "geneticSexIDM", "Rock_UnitQuaternary.deposit", "Rock_UnitQuaternary.talus...rock.glacier", "Rock_UnitSandstone.or.conglomerate.or.siltstone", "Rock_UnitShale.or.limestone")

# add year as a covariate (treat as fixed effect since we only have 3 years and that's too small to estimate independent slopes with random effect)
# make a gathered df for each pheno, then join them together matching on year and id name

gathered_df <- list()

for (pheno in c("GSL", "OGI", "EVImax", "OGMn")){
  
  pheno_columns <- c(colnames(pheno_coded[grep(pheno, colnames(pheno_coded))]))
  
  gathered_df[[pheno]] <- gather(data = pheno_coded, key = year, value = measurement, all_of(pheno_columns))
  gathered_df[[pheno]] <- gathered_df[[pheno]] %>%
    mutate(year = paste('y_', substr(year, 7, 10), sep='')) %>%
    select(Site_Code, year, measurement)
  names(gathered_df[[pheno]])[names(gathered_df[[pheno]]) == "measurement"] <- pheno
}

gathered_df <- gathered_df %>% reduce(left_join, by=c('Site_Code', 'year'))

pheno_coded_stub <- pheno_coded %>% select(colnames(pheno_coded[grep("pheno", colnames(pheno_coded), invert = T)]))
long_coded_df <- left_join(pheno_coded_stub, gathered_df, on="Site_Code")

# dummy encode the years and make full rank
long_coded_df <- long_coded_df %>% mutate(value=1) %>% spread(year, value, fill=0)
long_coded_df <- long_coded_df %>% select(-'y_2016')

#add ploidy to later split data into diploid and triploid
long_coded_df <- ploidy %>% select(Site_Code, ploidy_group) %>%
  distinct() %>%
  right_join(long_coded_df, on='Site_Code')



# split into all the sets I want to do analysis on
df_subsets <- list()

#divide into ploidies
df_subsets[['coverany_all']] <- long_coded_df
df_subsets[['coverany_diploid']] <- long_coded_df %>% filter(ploidy_group=='Diploid')
df_subsets[['coverany_nodips']] <- long_coded_df %>% filter(ploidy_group!='Diploid') #include ambiguous and unknown since most are triploid
df_subsets[['coverany_triploid']] <- long_coded_df %>% filter(ploidy_group=='Triploid') 

#remove samples with aspen cover less than 20%
df_subsets[['cover20_all']] <- long_coded_df %>% filter(fraction_aspen >=.2)
df_subsets[['cover20_diploid']] <- long_coded_df %>% filter(ploidy_group=='Diploid' & fraction_aspen >=.2)
df_subsets[['cover20_nodips']] <- long_coded_df %>% filter(ploidy_group!='Diploid' & fraction_aspen >=.2)
df_subsets[['cover20_triploid']] <- long_coded_df %>% filter(ploidy_group=='Triploid' & fraction_aspen >=.2)

#remove samples with aspen cover less than 25%
df_subsets[['cover25_all']] <- long_coded_df %>% filter(fraction_aspen >=.25)
df_subsets[['cover25_diploid']] <- long_coded_df %>% filter(ploidy_group=='Diploid' & fraction_aspen >=.25)
df_subsets[['cover25_nodips']] <- long_coded_df %>% filter(ploidy_group!='Diploid' & fraction_aspen >=.25)
df_subsets[['cover25_triploid']] <- long_coded_df %>% filter(ploidy_group=='Triploid' & fraction_aspen >=.25)

#remove samples with aspen cover less than 50%
df_subsets[['cover50_all']] <- long_coded_df %>% filter(fraction_aspen >=.50)
df_subsets[['cover50_diploid']] <- long_coded_df %>% filter(ploidy_group=='Diploid' & fraction_aspen >=.50)
df_subsets[['cover50_nodips']] <- long_coded_df %>% filter(ploidy_group!='Diploid' & fraction_aspen >=.50)
df_subsets[['cover50_triploid']] <- long_coded_df %>% filter(ploidy_group=='Triploid' & fraction_aspen >=.50)



## assign the regressors

regressors <- covariates %>% paste(collapse=' + ')
regressors <- paste(c(regressors,'y_2017 + y_2018 + y_2019'), collapse = " + ")


## do the regression
## running on all samples regardless of available genotype to benefit from larger dataset
lm_func <- function(df){
  
  lm_df <- list()
  for (pheno in c("GSL", "OGI", "EVImax", "OGMn")){
    lm_df[[pheno]] <- lm(paste(pheno,"~", regressors), data = df, na.action = na.exclude)
  }
  return(lm_df)
}

res <- lapply(df_subsets, function(x) lm_func(x))



#check for linear relationship between independent and dependent variables
#(so phenology and covariates)
#and homoscedasticity
#as well as leverage for any given outlier sample


# write out the lm summary and plots for each condition and subset
for (i in seq(1,16)){
  subset <- names(res)[i]
  print(subset)
  for (j in seq(1,4)){
    pheno <- names(res[[i]][j])
    print(pheno)
    pdf(sprintf("results/covariate_residuals/regression_plots/%s_%s.pdf", subset, pheno))
    par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
    plot(res[[i]][[j]])
    par(mfrow=c(1,1)) # Change back to 1 x 1
    dev.off()
  }
}



# write out the residuals for each condition and subset
for (i in seq(1,16)){
  subset <- names(res)[i]
  print(subset)
  for (j in seq(1,4)){
    pheno <- names(res[[i]][j])
    print(pheno)
    residuals <- res[[i]][[j]]$residuals
    inds <- as.numeric(names(residuals))
    sitecode <- df_subsets[[i]]$Site_Code
    phenval <- df_subsets[[i]][[pheno]]
    sitecodes2use <- sitecode[inds]
    phenval2use <- phenval[inds]
    resid_df <- tibble(fam_id = 0,
                       ind_id = sitecodes2use,
                       residuals =residuals,
                       phenval = phenval2use)
    write_tsv(resid_df, 
              sprintf("results/covariate_residuals/%s_%s.phen", subset, pheno),
              col_names = F)
  }
}



# look at the GRM
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

grm <- ReadGRMBin("data/gcta/aspen_phenology")
# I actually went with .005 and removed the other files
#grm_01 <- ReadGRMBin("data/gcta/aspen_phenology_01")
#grm_005 <- ReadGRMBin("data/gcta/aspen_phenology_005")


############
# read in results from gcta greml
gcta_files <- list.files(path = "results/gcta_greml/", pattern = "*.hsq")

gcta_df <- do.call(rbind, lapply(gcta_files, function(x) cbind(
  read_tsv(paste("results/gcta_greml/", x, sep="")), stub=x)))

gcta_df <- gcta_df %>% separate(stub, into=c("cover", "ploidy", "rest"), sep="_")
gcta_df <- gcta_df %>% separate(rest, into=c("phenology", "val_used"), sep="[.]")

val_used_resid <- gcta_df %>% 
  filter(val_used == 'resids' & Source == 'V(G)/Vp' & ploidy != 'nodips')
val_used_full <- gcta_df %>% 
  filter(val_used == 'pheno' & Source == 'V(G)/Vp' & ploidy != 'nodips')

plot_resid <- ggplot(val_used_resid, aes(x=ploidy, y=Variance, color=cover)) +
  geom_point(alpha=0.5, size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=(Variance-SE), ymax=(Variance+SE)),
                    position=position_dodge(width=0.5)) +
  facet_grid(~phenology, switch="both") +
  theme(axis.text.x = element_text(angle = 75, hjust=1)) + 
  labs(y="Vg/Vp", x='covariate residuals') +
  coord_cartesian(ylim = c(0, 1))

plot_full <- ggplot(val_used_full, aes(x=ploidy, y=Variance, color=cover)) +
  geom_point(alpha=0.5, size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=(Variance-SE), ymax=(Variance+SE)), 
                position=position_dodge(width=0.5)) +
  facet_grid(~phenology, switch="both") +
  theme(axis.text.x = element_text(angle = 75, hjust=1)) + 
  labs(y="Vg/Vp", x='full phenology') +
  coord_cartesian(ylim = c(0, 1))

plot <- ggarrange(plot_resid, plot_full, ncol=2, common.legend=T, legend="bottom")

annotate_figure(plot, top=textGrob("Heritability Estimates", 
                                   gp=gpar(fontsize=20, font=8)))


ggsave("results/gcta_greml/plots/heritability.pdf", plot,
       width=10, height=7, units="in")

#look at genotypic and phenotypic variance

val_used_resid <- gcta_df %>% 
  filter(val_used == 'resids' & Source %in% (c("V(G)", "V(E)", "Vp")) & ploidy != 'nodips')
val_used_full <- gcta_df %>% 
  filter(val_used == 'phenos' & SSource %in% (c("V(G)", "V(E)", "Vp")) & ploidy != 'nodips')

plot_resid <- ggplot(val_used_resid, aes(x=ploidy, y=Variance, shape=cover, color=Source)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(ymin=(Variance-SE), ymax=(Variance+SE))) +
  facet_grid(~phenology, switch="both") +
  theme(axis.text.x = element_text(angle = 75, hjust=1)) + 
  labs(y="Variance", x='covariate residuals')

plot_full <- ggplot(val_used_full, aes(x=ploidy, y=Variance, shape=cover, color=Source)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(ymin=(Variance-SE), ymax=(Variance+SE))) +
  facet_grid(~phenology, switch="both") +
  theme(axis.text.x = element_text(angle = 75, hjust=1)) + 
  labs(y="Variance", x='full phenology')

require(gridExtra)
plot_variances <- grid.arrange(plot_resid, plot_full, ncol=2,
                                  top = "Contributing Variances")

ggsave("results/gcta_greml/plots/variances.png", plot_variances)


then

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

'''

##############################################################
# GEMMA prep
##############################################################

# Prep files for LMM input

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

write_tsv(as_tibble(pheno_coded$Site_Code), "data/aspen.ids", col_names = FALSE) # output initial set of sample ids

# subset to samples with genotypes
geno_fam <- read_table2("data/plink_cleanup_2.fam", col_names = F)
pheno_coded <- pheno_coded %>% filter(Site_Code %in% geno_fam$X2)

#fill in missing sex and ploidy values with mean (0.5)
pheno_coded  <- pheno_coded %>% mutate(geneticSexIDM = ifelse(is.na(geneticSexIDM), 0.5, geneticSexIDM),
                                       CytotypeTriploid = ifelse(is.na(CytotypeTriploid), 0.5, CytotypeTriploid))

###### output environmental variables to run in GEMMA LMM ######

phenos <- colnames(pheno_coded[grep("pheno", colnames(pheno_coded))])

write_tsv(select(pheno_coded, all_of(phenos)), "data/phenotypes.tsv", col_names = F) # output phenology phenos
write_tsv(select(pheno_coded, all_of(phenos))[0,], "data/phenotypes.labels")
write_tsv(select(pheno_coded, intercept, colnames(pheno_coded[2:17])), "data/covariates.tsv", col_names = F) # output covariates
write_tsv(pheno_coded[0, c(39,2:17)], "data/covariates.labels") # output covariate names



