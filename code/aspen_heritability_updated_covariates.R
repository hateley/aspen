# Author: Dr. Shannon Hateley
# Date: Fall 2021 - Spring 2022
# Project: Aspen Phenology with Dr. Benjamin Blonder

###############################
# Heritability of traits
# 
# 1) remove confounders
#    phenology = f(elve, slope, snowmelt...) + e resid (model cytotype separately)
# 2) map residuals
#    e phenology resids = u + bx + Z(pop struc) + e (using TCGA)
#
# phenology = u + bx + pcs + cov. + .. e 
#
# sets to run heritability on:
# diploids, triploids, all
# all aspen coverage, greater than 25%, 50%


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


##############################################################
# load data
##############################################################

phenology <- read_csv("data/data_for_heritability_updated.csv")
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

var_cor <- cor(pheno_coded[-1], use="pairwise.complete.obs", method="pearson")

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
ggsave("results/phenology_correlation_snowmelt.png", ggheatmap)
ggsave("results/phenology_correlation_snowmelt.pdf", ggheatmap)


#plot xy coordinates
ggplot(data=phenology, aes(x=x, y=y)) +
  geom_point()


##############################################################
# pull out ploidy info
##############################################################

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


new_covars <- c("tmax", "snowmelt", "sm_runs", "sm_q01")
new_covars <- lapply(new_covars, function(x) colnames(pheno_coded[grep(
  x, colnames(pheno_coded))])) %>% unlist()

  
  
covariates <- c("Elevation", "Cos.aspect", "Slope", "Summer.Insolation", "DBH.mean", 
                "Canopy_openness", "Soil.typeSoil", "Soil.typeTalus", "geneticSexIDM",
                "Rock_UnitQuaternary.deposit", "Rock_UnitQuaternary.talus...rock.glacier",
                "Rock_UnitSandstone.or.conglomerate.or.siltstone", "Rock_UnitShale.or.limestone",
                new_covars)

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

# fix the yearly values to be relative to phenol year

#tmax is pretty much the same for everything:
#> pheno_coded %>% select(contains('tmax')) %>% apply(2, max)
#tmax_q99.X2012 tmax_q99.X2013 tmax_q99.X2014 tmax_q99.X2015 tmax_q99.X2016 
#29.97137       29.97811       29.98526       29.87108       29.94518 
#tmax_q99.X2017 tmax_q99.X2018 tmax_q99.X2019 
#29.96096       30.05194       29.94654 
#> pheno_coded %>% select(contains('tmax')) %>% apply(2, min)
#tmax_q99.X2012 tmax_q99.X2013 tmax_q99.X2014 tmax_q99.X2015 tmax_q99.X2016 
#29.56789       29.56454       29.42282       29.41036       29.55413 
#tmax_q99.X2017 tmax_q99.X2018 tmax_q99.X2019 
#29.47145       29.65711       29.50122 

# if year = 2019, 2019 measurements = current, 2018 y-1, 2017 y-2, y-3
long_coded_df <- long_coded_df %>% mutate(snowmelt_t0 = case_when(year=="y_2019" ~ snowmelt.X2019,
                                                   year=="y_2018" ~ snowmelt.X2018,
                                                   year=="y_2017" ~ snowmelt.X2017,
                                                   year=="y_2016" ~ snowmelt.X2016)) %>%
  mutate(snowmelt_t1 = case_when(year=="y_2019" ~ snowmelt.X2018,
                                  year=="y_2018" ~ snowmelt.X2017,
                                  year=="y_2017" ~ snowmelt.X2016,
                                  year=="y_2016" ~ snowmelt.X2015)) %>%
  mutate(snowmelt_t2 = case_when(year=="y_2019" ~ snowmelt.X2017,
                                 year=="y_2018" ~ snowmelt.X2016,
                                 year=="y_2017" ~ snowmelt.X2015,
                                 year=="y_2016" ~ snowmelt.X2014)) %>%
  mutate(snowmelt_t3 = case_when(year=="y_2019" ~ snowmelt.X2016,
                                 year=="y_2018" ~ snowmelt.X2015,
                                 year=="y_2017" ~ snowmelt.X2014,
                                 year=="y_2016" ~ snowmelt.X2013))

long_coded_df <- long_coded_df %>% mutate(sm_runs_med_dur_t0 = case_when(year=="y_2019" ~ sm_runs_med_dur.X2019,
                                                        year=="y_2018" ~ sm_runs_med_dur.X2018,
                                                        year=="y_2017" ~ sm_runs_med_dur.X2017,
                                                        year=="y_2016" ~ sm_runs_med_dur.X2016)) %>%
  mutate(sm_runs_med_dur_t1 = case_when(year=="y_2019" ~ sm_runs_med_dur.X2018,
                                 year=="y_2018" ~ sm_runs_med_dur.X2017,
                                 year=="y_2017" ~ sm_runs_med_dur.X2016,
                                 year=="y_2016" ~ sm_runs_med_dur.X2015)) %>%
  mutate(sm_runs_med_dur_t2 = case_when(year=="y_2019" ~ sm_runs_med_dur.X2017,
                                 year=="y_2018" ~ sm_runs_med_dur.X2016,
                                 year=="y_2017" ~ sm_runs_med_dur.X2015,
                                 year=="y_2016" ~ sm_runs_med_dur.X2014)) %>%
  mutate(sm_runs_med_dur_t3 = case_when(year=="y_2019" ~ sm_runs_med_dur.X2016,
                                 year=="y_2018" ~ sm_runs_med_dur.X2015,
                                 year=="y_2017" ~ sm_runs_med_dur.X2014,
                                 year=="y_2016" ~ sm_runs_med_dur.X2013))

long_coded_df <- long_coded_df %>% mutate(sm_q01_t0 = case_when(year=="y_2019" ~ sm_q01.X2019,
                                                                         year=="y_2018" ~ sm_q01.X2018,
                                                                         year=="y_2017" ~ sm_q01.X2017,
                                                                         year=="y_2016" ~ sm_q01.X2016)) %>%
  mutate(sm_q01_t1 = case_when(year=="y_2019" ~ sm_q01.X2018,
                                        year=="y_2018" ~ sm_q01.X2017,
                                        year=="y_2017" ~ sm_q01.X2016,
                                        year=="y_2016" ~ sm_q01.X2015)) %>%
  mutate(sm_q01_t2 = case_when(year=="y_2019" ~ sm_q01.X2017,
                                        year=="y_2018" ~ sm_q01.X2016,
                                        year=="y_2017" ~ sm_q01.X2015,
                                        year=="y_2016" ~ sm_q01.X2014)) %>%
  mutate(sm_q01_t3 = case_when(year=="y_2019" ~ sm_q01.X2016,
                                        year=="y_2018" ~ sm_q01.X2015,
                                        year=="y_2017" ~ sm_q01.X2014,
                                        year=="y_2016" ~ sm_q01.X2013)) 


long_coded_df <- long_coded_df %>% mutate(tmax_q99_t0 = case_when(year=="y_2019" ~ tmax_q99.X2019,
                                                                year=="y_2018" ~ tmax_q99.X2018,
                                                                year=="y_2017" ~ tmax_q99.X2017,
                                                                year=="y_2016" ~ tmax_q99.X2016)) %>%
  mutate(tmax_q99_t1 = case_when(year=="y_2019" ~ tmax_q99.X2018,
                               year=="y_2018" ~ tmax_q99.X2017,
                               year=="y_2017" ~ tmax_q99.X2016,
                               year=="y_2016" ~ tmax_q99.X2015)) %>%
  mutate(tmax_q99_t2 = case_when(year=="y_2019" ~ tmax_q99.X2017,
                               year=="y_2018" ~ tmax_q99.X2016,
                               year=="y_2017" ~ tmax_q99.X2015,
                               year=="y_2016" ~ tmax_q99.X2014)) %>%
  mutate(tmax_q99_t3 = case_when(year=="y_2019" ~ tmax_q99.X2016,
                               year=="y_2018" ~ tmax_q99.X2015,
                               year=="y_2017" ~ tmax_q99.X2014,
                               year=="y_2016" ~ tmax_q99.X2013))


long_coded_df <- long_coded_df %>% select(-all_of(new_covars))

# dummy encode the years and make full rank
long_coded_df <- long_coded_df %>% mutate(value=1) %>% spread(year, value, fill=0)
long_coded_df <- long_coded_df %>% select(-'y_2016')



#################################################


# look at relationship between predictors and phenologies

plot_cov2phenol_relationship <- function(long_df){
  df_list = list()
  for (p in c("GSL", "OGI", "OGMn", "EVImax")){
    print(p)
    df <- long_df %>% select(c(1:18, p, 23:42))
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
ggsave("results/cov_phen_scatterplots/GSL_updated.pdf", plot, height = 6, width=10)

plot <- ggplot(plot_dfs[["OGI"]], aes(x = value, y = OGI)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~variable, scales = 'free') +
  geom_smooth(method = "lm", se = FALSE, color='red') +
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label", size=3) +
  xlab("covariate relation to OGI")

plot
ggsave("results/cov_phen_scatterplots/OGI_updated.pdf", plot, height = 6, width=10)

plot <- ggplot(plot_dfs[["OGMn"]], aes(x = value, y = OGMn)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~variable, scales = 'free') +
  geom_smooth(method = "lm", se = FALSE, color='red') +
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label", size=3) +
  xlab("covariate relation to OGMn")
plot
ggsave("results/cov_phen_scatterplots/OGMn_updated.pdf", plot, height = 6, width=10)

plot <- ggplot(plot_dfs[["EVImax"]], aes(x = value, y = EVImax)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~variable, scales = 'free') +
  geom_smooth(method = "lm", se = FALSE, color='red') +
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label", size =3) +
  xlab("covariate relation to EVImax")
plot
ggsave("results/cov_phen_scatterplots/EVImax_updated.pdf", plot, height = 6, width=10)


#################################################

#plot xy coordinates with temp
p <- new_covars[grep("tmax", new_covars)]
print(p)
df <- pheno_coded %>% select(c(1:3, p))
tmp <- df %>% gather(key=year, val=max_temp, p, -Site_Code, -x, -y) 

plot <- ggplot(tmp, aes(x=x, y=y, color=max_temp)) +
  geom_point(alpha=0.8) +
  facet_wrap(~year) +
  xlab("temperature by coordinate")
plot
ggsave("results/temp_by_coordinates.pdf", plot, height = 6, width=10)

#plot xy coordinates with snowmelt
p <- new_covars[grep("snowmelt", new_covars)]
print(p)
df <- pheno_coded %>% select(c(1:3, p))
tmp <- df %>% gather(key=year, val=melt_date, p, -Site_Code, -x, -y) 

plot <- ggplot(tmp, aes(x=x, y=y, color=melt_date)) +
  geom_point(alpha=0.8) +
  facet_wrap(~year) +
  xlab("snowmelt date by coordinate")
plot
ggsave("results/melt_date_by_coordinates.pdf", plot, height = 6, width=10)


#################################################

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

#remove samples with aspen cover less than 25%
df_subsets[['cover25_all']] <- long_coded_df %>% filter(aspen_cover >=.25)
df_subsets[['cover25_diploid']] <- long_coded_df %>% filter(ploidy_group=='Diploid' & aspen_cover >=.25)
df_subsets[['cover25_nodips']] <- long_coded_df %>% filter(ploidy_group!='Diploid' & aspen_cover >=.25)
df_subsets[['cover25_triploid']] <- long_coded_df %>% filter(ploidy_group=='Triploid' & aspen_cover >=.25)

#remove samples with aspen cover less than 50%
df_subsets[['cover50_all']] <- long_coded_df %>% filter(aspen_cover >=.50)
df_subsets[['cover50_diploid']] <- long_coded_df %>% filter(ploidy_group=='Diploid' & aspen_cover >=.50)
df_subsets[['cover50_nodips']] <- long_coded_df %>% filter(ploidy_group!='Diploid' & aspen_cover >=.50)
df_subsets[['cover50_triploid']] <- long_coded_df %>% filter(ploidy_group=='Triploid' & aspen_cover >=.50)



## assign the regressors

regressors <- colnames(long_coded_df[c(5:12, 14:18, 24:42)]) %>% paste(collapse=' + ')

####################################################################################
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



# write out the lm summary and plots for each condition and subset
for (i in seq(1,12)){
  subset <- names(res)[i]
  print(subset)
  for (j in seq(1,4)){
    pheno <- names(res[[i]][j])
    print(pheno)
    pdf(sprintf("results/covariate_residuals/regression_plots_updated/%s_%s.pdf", subset, pheno))
    par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
    plot(res[[i]][[j]])
    par(mfrow=c(1,1)) # Change back to 1 x 1
    dev.off()
  }
}



# write out the residuals for each condition and subset
for (i in seq(1,12)){
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
              sprintf("results/covariate_residuals/updated/%s_%s.phen", subset, pheno),
              col_names = F)
  }
}


####################################################################################
# RUN GCTA GREML, then do this
####################################################################################


# read in results from gcta greml
gcta_files <- list.files(path = "results/gcta_greml/updated/", pattern = "*.hsq")

gcta_df <- do.call(rbind, lapply(gcta_files, function(x) cbind(
  read_tsv(paste("results/gcta_greml/updated/", x, sep="")), stub=x)))

gcta_df <- gcta_df %>% separate(stub, into=c("cover", "ploidy", "rest"), sep="_")
gcta_df <- gcta_df %>% separate(rest, into=c("phenology", "val_used"), sep="[.]")
gcta_df <- gcta_df %>% filter(cover != "coverany")

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
  theme_minimal_hgrid() +
  labs(y="Vg/Vp", x='covariate residuals') +
  coord_cartesian(ylim = c(0, 1))

plot_full <- ggplot(val_used_full, aes(x=ploidy, y=Variance, color=cover)) +
  geom_point(alpha=0.5, size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=(Variance-SE), ymax=(Variance+SE)), 
                position=position_dodge(width=0.5)) +
  facet_grid(~phenology, switch="both") +
  theme(axis.text.x = element_text(angle = 75, hjust=1)) + 
  theme_minimal_hgrid() +
  labs(y="Vg/Vp", x='full phenology') +
  coord_cartesian(ylim = c(0, 1))

plot <- ggarrange(plot_resid, plot_full, ncol=2, common.legend=T, legend="bottom")

annotate_figure(plot, top=textGrob("Heritability Estimates", 
                                   gp=gpar(fontsize=20, font=8)))

ggsave("results/gcta_greml/plots/updated_heritability.pdf", plot,
       width=10, height=7, units="in")

#####
# get counts for paper

tmp2 <- long_coded_df %>% filter(
  aspen_cover >= 0.25) %>%
  distinct(Site_Code, .keep_all=T)



#################################################################
# output dataframe for Ben to make publication plot
#################################################################

paper_df <- val_used_resid %>%
  filter(cover %in% c("cover25", "cover50")) %>%
  select(-val_used)
write_tsv(paper_df, "results/gcta_greml/updated/paper_data_updated.tsv")


# set up environment
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

paper_df <- read_tsv("results/gcta_greml/updated/paper_data_updated.tsv")

plot_resid <- ggplot(paper_df, aes(x=ploidy, y=Variance, color=cover)) +
  geom_point(alpha=0.5, size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=(Variance-SE), ymax=(Variance+SE)),
                position=position_dodge(width=0.5)) +
  facet_grid(~phenology, switch="both") +
  theme(axis.text.x = element_text(angle = 75, hjust=1)) + 
  labs(y="Vg/Vp", x='covariate residuals', title="Heritability Estimates") +
  coord_cartesian(ylim = c(0, 1))

plot_resid

ggsave("results/gcta_greml/plots/heritability_25and50_updated.pdf", plot_resid,
       width=10, height=7, units="in")



########################

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



