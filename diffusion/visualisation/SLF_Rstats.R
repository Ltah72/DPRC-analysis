#Perform statistics on tract of interest (TOI) from diffusion analysis. This script will perform analysis for group 
#differences in the TOI which was selected from the user in the previous pipelines, done mainly from the software
#programme, MRtrix3. The main outputs are the text files of the fixel-based analysis (FBA) metrics from the TOI.m pipleine. 
#For this script, I have chosen to analyse the superior longitudinal fasciculus (SLF) between the 5 participant groups of 
#interest - Controls (C), subjective cognitive decline (SCD), amnestic mild cognitive impairment (aMCI), multiple-domain
#mild cognitive impairment (mMCI), and Alzheimer's disease (AD). 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 8/12/20

#install packages/open libraries
pacman::p_load(dplyr, ggplot2, psych, car, multcomp, lsr, BayesFactor)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#first read in the covariates group data file: 
setwd("V:/PartInfo")
covariates_data <- read.csv("covariates-participants-lined-up_update.csv")
#convert the group data to a factor variable
covariates_data$Group <- as.factor(covariates_data$Group)


#navigate to the correct pathway: 
setwd("H:/ltah262/NECTAR_data/LENORE/derivatives/F0/diff_data/test/metric_results_SLF")


#read in text file of data
SLF_data <- cbind.data.frame(read.table("FD_SLF_whole_TOI.txt", header = T), read.table("FC_log_SLF_whole_TOI.txt", header = T), read.table("FDC_SLF_whole_TOI.txt", header = T), read.table("FD_SLF_L_TOI.txt", header = T), read.table("FC_log_SLF_L_TOI.txt", header = T), read.table("FDC_SLF_L_TOI.txt", header = T), read.table("FD_SLF_R_TOI.txt", header = T), read.table("FC_log_SLF_R_TOI.txt", header = T), read.table("FDC_SLF_R_TOI.txt", header = T))
#rename columns with each FBA metric
colnames(SLF_data) <- c("mn_FD_SLF", "md_FD_SLF", "std_FD_SLF", "std_rv_FD_SLF", "min_FD_SLF", "max_FD_SLF", "count_FD_SLF", "mn_FC_SLF", "md_FC_SLF", "std_FC_SLF", "std_rv_FC_SLF", "min_FC_SLF", "max_FC_SLF", "count_FC_SLF", "mn_FDC_SLF", "md_FDC_SLF", "std_FDC_SLF", "std_rv_FDC_SLF", "min_FDC_SLF", "max_FDC_SLF", "count_FDC_SLF", "mn_FD_SLF_L", "md_FD_SLF_L", "std_FD_SLF_L", "std_rv_FD_SLF_L", "min_FD_SLF_L", "max_FD_SLF_L", "count_FD_SLF_L", "mn_FC_SLF_L", "md_FC_SLF_L", "std_FC_SLF_L", "std_rv_FC_SLF_L", "min_FC_SLF_L", "max_FC_SLF_L", "count_FC_SLF_L", "mn_FDC_SLF_L", "md_FDC_SLF_L", "std_FDC_SLF_L", "std_rv_FDC_SLF_L", "min_FDC_SLF_L", "max_FDC_SLF_L", "count_FDC_SLF_L", "mn_FD_SLF_R", "md_FD_SLF_R", "std_FD_SLF_R", "std_rv_FD_SLF_R", "min_FD_SLF_R", "max_FD_SLF_R", "count_FD_SLF_R", "mn_FC_SLF_R", "md_FC_SLF_R", "std_FC_SLF_R", "std_rv_FC_SLF_R", "min_FC_SLF_R", "max_FC_SLF_R", "count_FC_SLF_R", "mn_FDC_SLF_R", "md_FDC_SLF_R", "std_FDC_SLF_R", "std_rv_FDC_SLF_R", "min_FDC_SLF_R", "max_FDC_SLF_R", "count_FDC_SLF_R")
#add in Group classification column for each participant - data from another worksheet
SLF_data$Group <- covariates_data$Group
#add in covariates (clinical site) from the covariates_data dataframe to the SLF dataframe
SLF_data$ClinSite_name <- covariates_data$Clinical_site 
SLF_data$Age <- covariates_data$Age 
SLF_data$Sex <- covariates_data$Sex


#look at descriptive stats of the whole SLF FBA metrics between groups
SLF_FD_descrip <- describeBy(SLF_data$mn_FD_SLF, SLF_data$Group)
SLF_FC_descrip <- describeBy(SLF_data$mn_FC_SLF, SLF_data$Group)
SLF_FDC_descrip <- describeBy(SLF_data$mn_FDC_SLF, SLF_data$Group)

#look at descriptive stats of the Left SLF FBA metrics between groups
SLF_FD_L_descrip <- describeBy(SLF_data$mn_FD_SLF_L, SLF_data$Group)
SLF_FC_L_descrip <- describeBy(SLF_data$mn_FC_SLF_L, SLF_data$Group)
SLF_FDC_L_descrip <- describeBy(SLF_data$mn_FDC_SLF_L, SLF_data$Group)

#look at descriptive stats of the Right SLF FBA metrics between groups
SLF_FD_R_descrip <- describeBy(SLF_data$mn_FD_SLF_R, SLF_data$Group)
SLF_FC_R_descrip <- describeBy(SLF_data$mn_FC_SLF_R, SLF_data$Group)
SLF_FDC_R_descrip <- describeBy(SLF_data$mn_FDC_SLF_R, SLF_data$Group)


#run ANOVA to see if there are significant differences between groups
#for FD
#mean
SLF_FD_mod <- lm(mn_FD_SLF ~ Group, data = SLF_data)
SLF_FD_mod_L <- lm(mn_FD_SLF_L ~ Group, data = SLF_data)
SLF_FD_mod_R <- lm(mn_FD_SLF_R ~ Group, data = SLF_data)

#include the covariate of age (run an ANCOVA) in model
SLF_FD_age_mod <- lm(mn_FD_SLF ~ Group + Age, data = SLF_data)
SLF_FD_age_mod_L <- lm(mn_FD_SLF_L ~ Group + Age, data = SLF_data)
SLF_FD_age_mod_R <- lm(mn_FD_SLF_R ~ Group + Age, data = SLF_data)

#include the covariate of clinical site (run an ANCOVA) in model
SLF_FD_clinsite_mod <- lm(mn_FD_SLF ~ Group + ClinSite_name, data = SLF_data)
SLF_FD_clinsite_mod_L <- lm(mn_FD_SLF_L ~ Group + ClinSite_name, data = SLF_data)
SLF_FD_clinsite_mod_R <- lm(mn_FD_SLF_R ~ Group + ClinSite_name, data = SLF_data)

#include the covariate of all 3 variables (age, gender, clinical site) (run an ANCOVA) in model
SLF_FD_3covar_mod <- lm(mn_FD_SLF ~ Group + Age + Sex + ClinSite_name, data = SLF_data)
SLF_FD_3covar_mod_L <- lm(mn_FD_SLF_L ~ Group + Age + Sex + ClinSite_name, data = SLF_data)
SLF_FD_3covar_mod_R <- lm(mn_FD_SLF_R ~ Group + Age + Sex + ClinSite_name, data = SLF_data)


#median
#SLF_FD_mod <- lm(md_FD ~ Group, data = SLF_data)
#run ANOVA
anova(SLF_FD_mod)
anova(SLF_FD_mod_L)
anova(SLF_FD_mod_R)
#run Bayesian ANOVA
anovaBF(mn_FD_SLF ~ Group, data = SLF_data) 
anovaBF(mn_FD_SLF_L ~ Group, data = SLF_data) 
anovaBF(mn_FD_SLF_R ~ Group, data = SLF_data) 
#run ANCOVA
anova(SLF_FD_age_mod)
anova(SLF_FD_age_mod_L)
anova(SLF_FD_age_mod_R)
anova(SLF_FD_clinsite_mod)
anova(SLF_FD_clinsite_mod_L)
anova(SLF_FD_clinsite_mod_R)
anova(SLF_FD_3covar_mod)
anova(SLF_FD_3covar_mod_L)
anova(SLF_FD_3covar_mod_R)
#run Bayesian ANCOVA
#anovaBF(mn_FD_SLF ~ Group + Age, data = SLF_data) 
#anovaBF(mn_FD_SLF_L ~ Group, data = SLF_data) 
#anovaBF(mn_FD_SLF_R ~ Group, data = SLF_data) 
#to get multiple regression outputs
summary(SLF_FD_mod)
summary(SLF_FD_mod_L)
summary(SLF_FD_mod_R)
#for ANCOVA
summary(SLF_FD_age_mod)
summary(SLF_FD_age_mod_L)
summary(SLF_FD_age_mod_R)
summary(SLF_FD_clinsite_mod)
summary(SLF_FD_clinsite_mod_L)
summary(SLF_FD_clinsite_mod_R)
#run pairwise comparisons, given that the F-test was significant. 
post_hoc_SLF_FD_mod <- glht(SLF_FD_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_mod)
confint(post_hoc_SLF_FD_mod)

post_hoc_SLF_FD_mod_L <- glht(SLF_FD_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_mod_L)
confint(post_hoc_SLF_FD_mod_L)

post_hoc_SLF_FD_mod_R <- glht(SLF_FD_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_mod_R)
confint(post_hoc_SLF_FD_mod_R)
#for ANCOVA - age
post_hoc_SLF_FD_age_mod <- glht(SLF_FD_age_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_age_mod)
confint(post_hoc_SLF_FD_age_mod)

post_hoc_SLF_FD_age_mod_L <- glht(SLF_FD_age_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_age_mod_L)
confint(post_hoc_SLF_FD_age_mod_L)

post_hoc_SLF_FD_age_mod_R <- glht(SLF_FD_age_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_age_mod_R)
confint(post_hoc_SLF_FD_age_mod_R)

#for ANCOVA - clinical site
post_hoc_SLF_FD_clinsite_mod <- glht(SLF_FD_clinsite_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_clinsite_mod)
confint(post_hoc_SLF_FD_clinsite_mod)

post_hoc_SLF_FD_clinsite_mod_L <- glht(SLF_FD_clinsite_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_clinsite_mod_L)
confint(post_hoc_SLF_FD_clinsite_mod_L)

post_hoc_SLF_FD_clinsite_mod_R <- glht(SLF_FD_clinsite_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_clinsite_mod_R)
confint(post_hoc_SLF_FD_clinsite_mod_R)


#for ANCOVA - all 3 covariates (age, sex, clinical site)
post_hoc_SLF_FD_3covar_mod <- glht(SLF_FD_3covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_3covar_mod)
confint(post_hoc_SLF_FD_3covar_mod)

post_hoc_SLF_FD_3covar_mod_L <- glht(SLF_FD_3covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_3covar_mod_L)
confint(post_hoc_SLF_FD_3covar_mod_L)

post_hoc_SLF_FD_3covar_mod_R <- glht(SLF_FD_3covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_3covar_mod_R)
confint(post_hoc_SLF_FD_3covar_mod_R)

#calculate the effect size (eta-squared)
etaSquared(SLF_FD_mod)
etaSquared(SLF_FD_mod_L)
etaSquared(SLF_FD_mod_R)
#for ANCOVA
etaSquared(SLF_FD_clinsite_mod)
etaSquared(SLF_FD_clinsite_mod_L)
etaSquared(SLF_FD_clinsite_mod_R)
etaSquared(SLF_FD_age_mod)
etaSquared(SLF_FD_age_mod_L)
etaSquared(SLF_FD_age_mod_R)


#conduct power analysis for the whole FD
SLF_FD_group_means <- c(SLF_FD_descrip$`1`$mean, SLF_FD_descrip$`2`$mean, SLF_FD_descrip$`3`$mean, SLF_FD_descrip$`4`$mean, SLF_FD_descrip$`5`$mean)
power_SLF_FD_n <- power.anova.test(groups = length(SLF_FD_group_means), between.var = anova(SLF_FD_mod)$`Sum Sq`[1], within.var = anova(SLF_FD_mod)$`Sum Sq`[2], power = .8, sig.level = 0.05)
power_SLF_FD_power<- power.anova.test(groups = length(SLF_FD_group_means), between.var = anova(SLF_FD_mod)$`Sum Sq`[1], within.var = anova(SLF_FD_mod)$`Sum Sq`[2], n = 41, sig.level = 0.05)
#with age covariate
power_SLF_age_FD_n <- power.anova.test(groups = 5, between.var = anova(SLF_FD_age_mod)$`Sum Sq`[1], within.var = anova(SLF_FD_age_mod)$`Sum Sq`[3], power = .8, sig.level = 0.05)
power_SLF_age_FD_power<- power.anova.test(groups = 5, between.var = anova(SLF_FD_age_mod)$`Sum Sq`[1], within.var = anova(SLF_FD_age_mod)$`Sum Sq`[3], n = 41, sig.level = 0.05)

#conduct power analysis for the left FD
SLF_FD_L_group_means <- c(SLF_FD_L_descrip$`1`$mean, SLF_FD_L_descrip$`2`$mean, SLF_FD_L_descrip$`3`$mean, SLF_FD_L_descrip$`4`$mean, SLF_FD_L_descrip$`5`$mean)
power_SLF_FD_L_n <- power.anova.test(groups = length(SLF_FD_L_group_means), between.var = anova(SLF_FD_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FD_mod_L)$`Sum Sq`[2], power = .8, sig.level = 0.05)
power_SLF_FD_L_power<- power.anova.test(groups = length(SLF_FD_L_group_means), between.var = anova(SLF_FD_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FD_mod_L)$`Sum Sq`[2], n = 41, sig.level = 0.05)
#with age covariate
power_SLF_age_FD_L_n <- power.anova.test(groups = 5, between.var = anova(SLF_FD_age_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FD_age_mod_L)$`Sum Sq`[3], power = .8, sig.level = 0.05)
power_SLF_age_FD_L_power<- power.anova.test(groups = 5, between.var = anova(SLF_FD_age_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FD_age_mod_L)$`Sum Sq`[3], n = 41, sig.level = 0.05)

#conduct power analysis for the right FD
SLF_FD_R_group_means <- c(SLF_FD_R_descrip$`1`$mean, SLF_FD_R_descrip$`2`$mean, SLF_FD_R_descrip$`3`$mean, SLF_FD_R_descrip$`4`$mean, SLF_FD_R_descrip$`5`$mean)
power_SLF_FD_R_n <- power.anova.test(groups = length(SLF_FD_R_group_means), between.var = anova(SLF_FD_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FD_mod_R)$`Sum Sq`[2], power = .8, sig.level = 0.05)
power_SLF_FD_R_power<- power.anova.test(groups = length(SLF_FD_R_group_means), between.var = anova(SLF_FD_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FD_mod_R)$`Sum Sq`[2], n = 41, sig.level = 0.05)
#with age covariate
power_SLF_age_FD_R_n <- power.anova.test(groups = 5, between.var = anova(SLF_FD_age_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FD_age_mod_R)$`Sum Sq`[3], power = .8, sig.level = 0.05)
power_SLF_age_FD_R_power<- power.anova.test(groups = 5, between.var = anova(SLF_FD_age_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FD_age_mod_R)$`Sum Sq`[3], n = 41, sig.level = 0.05)




#for FC
#mean
SLF_FC_mod <- lm(mn_FC_SLF ~ Group, data = SLF_data)
SLF_FC_mod_L <- lm(mn_FC_SLF_L ~ Group, data = SLF_data)
SLF_FC_mod_R <- lm(mn_FC_SLF_R ~ Group, data = SLF_data)
#include the covariate of clinical site (run an ANCOVA) in model
SLF_FC_clinsite_mod <- lm(mn_FC_SLF ~ Group + ClinSite_covar, data = SLF_data)
SLF_FC_clinsite_mod_L <- lm(mn_FC_SLF_L ~ Group + ClinSite_covar, data = SLF_data)
SLF_FC_clinsite_mod_R <- lm(mn_FC_SLF_R ~ Group + ClinSite_covar, data = SLF_data)

#median
#SLF_FC_mod <- lm(md_FC ~ Group, data = SLF_data)
anova(SLF_FC_mod)
anova(SLF_FC_mod_L)
anova(SLF_FC_mod_R)
#run Bayesian ANOVA
anovaBF(mn_FC_SLF ~ Group, data = SLF_data) 
anovaBF(mn_FC_SLF_L ~ Group, data = SLF_data) 
anovaBF(mn_FC_SLF_R ~ Group, data = SLF_data) 
#for ANCOVA
anova(SLF_FC_clinsite_mod)
anova(SLF_FC_clinsite_mod_L)
anova(SLF_FC_clinsite_mod_R)

#for FDC
#mean
SLF_FDC_mod <- lm(mn_FDC_SLF ~ Group, data = SLF_data)
SLF_FDC_mod_L <- lm(mn_FDC_SLF_L ~ Group, data = SLF_data)
SLF_FDC_mod_R <- lm(mn_FDC_SLF_R ~ Group, data = SLF_data)
#include the covariate of clinical site (run an ANCOVA) in model
SLF_FDC_clinsite_mod <- lm(mn_FDC_SLF ~ Group + ClinSite_covar, data = SLF_data)
SLF_FDC_clinsite_mod_L <- lm(mn_FDC_SLF_L ~ Group + ClinSite_covar, data = SLF_data)
SLF_FDC_clinsite_mod_R <- lm(mn_FDC_SLF_R ~ Group + ClinSite_covar, data = SLF_data)


#median
#SLF_FDC_mod <- lm(md_FDC ~ Group, data = SLF_data)
anova(SLF_FDC_mod)
anova(SLF_FDC_mod_L)
anova(SLF_FDC_mod_R)
#run Bayesian ANOVA
anovaBF(mn_FDC_SLF ~ Group, data = SLF_data) 
anovaBF(mn_FDC_SLF_L ~ Group, data = SLF_data) 
anovaBF(mn_FDC_SLF_R ~ Group, data = SLF_data) 
#for ANCOVA
anova(SLF_FDC_clinsite_mod)
anova(SLF_FDC_clinsite_mod_L)
anova(SLF_FDC_clinsite_mod_R)


#plot data
#For FD - mean 
#whole SLF FD (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FD_SLF)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)

#whole SLF FD (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FD_SLF, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#left SLF FD (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FD_SLF_L)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)

#left SLF FD (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FD_SLF_L, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF_L, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#right SLF FD (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FD_SLF_R)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)

#right SLF FD (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FD_SLF_R, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF_R, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()


#FD - median
# ggplot(SLF_data, aes(x = Group, y = md_FD)) + 
#     geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
#     xlab("Group_status") + 
#     ylab("FD") +
#     scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none")

#FC - mean, whole SLF FC (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FC_SLF)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)

#whole SLF FC (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FC_SLF, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FC_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#left SLF FC (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FC_SLF_L)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)

#left SLF FC (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FC_SLF_L, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FC_SLF_L, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#right SLF FC (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FC_SLF_R)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)


#right SLF FC (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FC_SLF_R, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FC_SLF_R, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#FC - median
# ggplot(SLF_data, aes(x = Group, y = md_FC)) + 
#     geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
#     xlab("Group_status") + 
#     ylab("FC") +
#     scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none")

#for FDC - mean 
#whole SLF FDC (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)


#whole SLF FDC (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FDC_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#left SLF FDC (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF_L)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)

#left SLF FDC (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF_L, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FDC_SLF_L, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#right SLF FDC (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF_R)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")+
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)

#right SLF FDC (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF_R, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FDC_SLF_R, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#FDC - median
# ggplot(SLF_data, aes(x = Group, y = md_FDC)) + 
#     geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
#     xlab("Group_status") + 
#     ylab("FDC") +
#     scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none")




# Analysis by clinical site as the independent variable - test for clinical site. 

#convert clinsite names to a factor variable
SLF_data$ClinSite_name <- as.factor(SLF_data$ClinSite_name)

#run ANOVA to see if there are significant differences between groups
#for FD
clinsite_SLF_FD_mod <- lm(mn_FD_SLF ~ ClinSite_name, data = SLF_data)
clinsite_SLF_FD_mod_L <- lm(mn_FD_SLF_L ~ ClinSite_name, data = SLF_data)
clinsite_SLF_FD_mod_R <- lm(mn_FD_SLF_R ~ ClinSite_name, data = SLF_data)
anova(clinsite_SLF_FD_mod)
anova(clinsite_SLF_FD_mod_L)
anova(clinsite_SLF_FD_mod_R)
#for FC
clinsite_SLF_FC_mod <- lm(mn_FC_SLF ~ ClinSite_name, data = SLF_data)
clinsite_SLF_FC_mod_L <- lm(mn_FC_SLF_L ~ ClinSite_name, data = SLF_data)
clinsite_SLF_FC_mod_R <- lm(mn_FC_SLF_R ~ ClinSite_name, data = SLF_data)
anova(clinsite_SLF_FC_mod)
anova(clinsite_SLF_FC_mod_L)
anova(clinsite_SLF_FC_mod_R)
#for FDC
clinsite_SLF_FDC_mod <- lm(mn_FDC_SLF ~ ClinSite_name, data = SLF_data)
clinsite_SLF_FDC_mod_L <- lm(mn_FDC_SLF_L ~ ClinSite_name, data = SLF_data)
clinsite_SLF_FDC_mod_R <- lm(mn_FDC_SLF_R ~ ClinSite_name, data = SLF_data)
anova(clinsite_SLF_FDC_mod)
anova(clinsite_SLF_FDC_mod_L)
anova(clinsite_SLF_FDC_mod_R)

#run pairwise comparisons, given that the F-test was significant. 
#for FD
post_hoc_clinsite_SLF_FD_mod <- glht(clinsite_SLF_FD_mod, linfct = mcp(ClinSite_name = "Tukey"))
summary(post_hoc_clinsite_SLF_FD_mod)
confint(post_hoc_clinsite_SLF_FD_mod)
post_hoc_clinsite_SLF_FD_mod_L <- glht(clinsite_SLF_FD_mod_L, linfct = mcp(ClinSite_name = "Tukey"))
summary(post_hoc_clinsite_SLF_FD_mod_L)
confint(post_hoc_clinsite_SLF_FD_mod_L)
post_hoc_clinsite_SLF_FD_mod_R <- glht(clinsite_SLF_FD_mod_R, linfct = mcp(ClinSite_name = "Tukey"))
summary(post_hoc_clinsite_SLF_FD_mod_R)
confint(post_hoc_clinsite_SLF_FD_mod_R)
#for FC
post_hoc_clinsite_SLF_FC_mod_L <- glht(clinsite_SLF_FC_mod_L, linfct = mcp(ClinSite_name = "Tukey"))
summary(post_hoc_clinsite_SLF_FC_mod_L)
confint(post_hoc_clinsite_SLF_FC_mod_L)
#for FDC
post_hoc_clinsite_SLF_FDC_mod_R <- glht(clinsite_SLF_FDC_mod_R, linfct = mcp(ClinSite_name = "Tukey"))
summary(post_hoc_clinsite_SLF_FDC_mod_R)
confint(post_hoc_clinsite_SLF_FDC_mod_R)

#plot data
#whole SLF FD (raincloud plot)
ggplot(SLF_data, aes(x = ClinSite_name, y = mn_FD_SLF, fill = ClinSite_name)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF, color = ClinSite_name), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = ClinSite_name)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = ClinSite_name)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()




# Analysis for 3 groups (Older adults (C(1), SCD(2)) vs. MCI (aMCI(3), mMCI(4)) vs. AD(5))

#add in group data and combine the groups from 5 groups into 3 groups
SLF_data$Three_Groups <- covariates_data$Group

SLF_data$Three_Groups[SLF_data$Three_Groups == "2"] <- "1" 
SLF_data$Three_Groups[SLF_data$Three_Groups == "4"] <- "3" 


#look at descriptive stats of the whole SLF FBA metrics between groups
SLF_FD_3groups_descrip <- describeBy(SLF_data$mn_FD_SLF, SLF_data$Three_Groups)
SLF_FC_3groups_descrip <- describeBy(SLF_data$mn_FC_SLF, SLF_data$Three_Groups)
SLF_FDC_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF, SLF_data$Three_Groups)

#look at descriptive stats of the Left SLF FBA metrics between groups
SLF_FD_L_3groups_descrip <- describeBy(SLF_data$mn_FD_SLF_L, SLF_data$Three_Groups)
SLF_FC_L_3groups_descrip <- describeBy(SLF_data$mn_FC_SLF_L, SLF_data$Three_Groups)
SLF_FDC_L_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF_L, SLF_data$Three_Groups)


#look at descriptive stats of the Right SLF FBA metrics between groups
SLF_FD_R_3groups_descrip <- describeBy(SLF_data$mn_FD_SLF_R, SLF_data$Three_Groups)
SLF_FC_R_3groups_descrip <- describeBy(SLF_data$mn_FC_SLF_R, SLF_data$Three_Groups)
SLF_FDC_R_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF_R, SLF_data$Three_Groups)


#run ANOVA to see if there are significant differences between groups
#for FD
#mean
SLF_FD_3groups_mod <- lm(mn_FD_SLF ~ Three_Groups, data = SLF_data)
SLF_FD_3groups_mod_L <- lm(mn_FD_SLF_L ~ Three_Groups, data = SLF_data)
SLF_FD_3groups_mod_R <- lm(mn_FD_SLF_R ~ Three_Groups, data = SLF_data)
#run ANOVA
anova(SLF_FD_3groups_mod)
anova(SLF_FD_3groups_mod_L)
anova(SLF_FD_3groups_mod_R)
#run Bayesian ANOVA
anovaBF(mn_FD_SLF ~ Three_Groups, data = SLF_data) 
anovaBF(mn_FD_SLF_L ~ Three_Groups, data = SLF_data) 
anovaBF(mn_FD_SLF_R ~ Three_Groups, data = SLF_data) 
#calculate the effect size (eta-squared)
etaSquared(SLF_FD_3groups_mod)
etaSquared(SLF_FD_3groups_mod_L)
etaSquared(SLF_FD_3groups_mod_R)

#for FC
#mean
SLF_FC_3groups_mod <- lm(mn_FC_SLF ~ Three_Groups, data = SLF_data)
SLF_FC_3groups_mod_L <- lm(mn_FC_SLF_L ~ Three_Groups, data = SLF_data)
SLF_FC_3groups_mod_R <- lm(mn_FC_SLF_R ~ Three_Groups, data = SLF_data)
#run ANOVA
anova(SLF_FC_3groups_mod)
anova(SLF_FC_3groups_mod_L)
anova(SLF_FC_3groups_mod_R)
#run Bayesian ANOVA
anovaBF(mn_FC_SLF ~ Three_Groups, data = SLF_data) 
anovaBF(mn_FC_SLF_L ~ Three_Groups, data = SLF_data) 
anovaBF(mn_FC_SLF_R ~ Three_Groups, data = SLF_data) 

#for FDC
#mean
SLF_FDC_3groups_mod <- lm(mn_FDC_SLF ~ Three_Groups, data = SLF_data)
SLF_FDC_3groups_mod_L <- lm(mn_FDC_SLF_L ~ Three_Groups, data = SLF_data)
SLF_FDC_3groups_mod_R <- lm(mn_FDC_SLF_R ~ Three_Groups, data = SLF_data)
#run ANOVA
anova(SLF_FDC_3groups_mod)
anova(SLF_FDC_3groups_mod_L)
anova(SLF_FDC_3groups_mod_R)
#run Bayesian ANOVA
anovaBF(mn_FDC_SLF ~ Three_Groups, data = SLF_data) 
anovaBF(mn_FDC_SLF_L ~ Three_Groups, data = SLF_data) 
anovaBF(mn_FDC_SLF_R ~ Three_Groups, data = SLF_data) 

#run pairwise comparisons, given that the F-test was significant. 
post_hoc_SLF_FD_3groups_mod <- glht(SLF_FD_3groups_mod, linfct = mcp(Three_Groups = "Tukey"))
summary(post_hoc_SLF_FD_3groups_mod)
confint(post_hoc_SLF_FD_3groups_mod)

post_hoc_SLF_FD_3groups_mod_L <- glht(SLF_FD_3groups_mod_L, linfct = mcp(Three_Groups = "Tukey"))
summary(post_hoc_SLF_FD_3groups_mod_L)
confint(post_hoc_SLF_FD_3groups_mod_L)

post_hoc_SLF_FD_3groups_mod_R <- glht(SLF_FD_3groups_mod_R, linfct = mcp(Three_Groups = "Tukey"))
summary(post_hoc_SLF_FD_3groups_mod_R)
confint(post_hoc_SLF_FD_3groups_mod_R)

#conduct power analysis for the whole FD
SLF_FD_3groups_means <- c(SLF_FD_3groups_descrip$`1`$mean, SLF_FD_3groups_descrip$`2`$mean, SLF_FD_3groups_descrip$`3`$mean, SLF_FD_3groups_descrip$`4`$mean, SLF_FD_3groups_descrip$`5`$mean)
power_SLF_FD_3groups_n <- power.anova.test(groups = length(SLF_FD_3groups_means), between.var = anova(SLF_FD_3groups_mod)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod)$`Sum Sq`[2], power = .8, sig.level = 0.05)
power_SLF_FD_3groups_power<- power.anova.test(groups = length(SLF_FD_3groups_means), between.var = anova(SLF_FD_3groups_mod)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod)$`Sum Sq`[2], n = 41, sig.level = 0.05)
#conduct power analysis for the left FD
SLF_FD_L_3groups_means <- c(SLF_FD_L_3groups_descrip$`1`$mean, SLF_FD_L_3groups_descrip$`2`$mean, SLF_FD_L_3groups_descrip$`3`$mean, SLF_FD_L_3groups_descrip$`4`$mean, SLF_FD_L_3groups_descrip$`5`$mean)
power_SLF_FD_L_3groups_n <- power.anova.test(groups = length(SLF_FD_L_3groups_means), between.var = anova(SLF_FD_3groups_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod_L)$`Sum Sq`[2], power = .8, sig.level = 0.05)
power_SLF_FD_L_3groups_power<- power.anova.test(groups = length(SLF_FD_L_3groups_means), between.var = anova(SLF_FD_3groups_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod_L)$`Sum Sq`[2], n = 41, sig.level = 0.05)
#conduct power analysis for the right FD
SLF_FD_R_3groups_means <- c(SLF_FD_R_3groups_descrip$`1`$mean, SLF_FD_R_3groups_descrip$`2`$mean, SLF_FD_R_3groups_descrip$`3`$mean, SLF_FD_R_3groups_descrip$`4`$mean, SLF_FD_R_3groups_descrip$`5`$mean)
power_SLF_FD_R_3groups_n <- power.anova.test(groups = length(SLF_FD_R_3groups_means), between.var = anova(SLF_FD_3groups_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod_R)$`Sum Sq`[2], power = .8, sig.level = 0.05)
power_SLF_FD_R_3groups_power<- power.anova.test(groups = length(SLF_FD_R_3groups_means), between.var = anova(SLF_FD_3groups_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod_R)$`Sum Sq`[2], n = 41, sig.level = 0.05)



#plot data
#for FD - mean
#whole SLF FD (violin plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)

#whole SLF FD (raincloud plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF, fill = Three_Groups)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF, color = Three_Groups), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#left SLF FD (violin plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF_L)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)

#left SLF FD (raincloud plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF_L, fill = Three_Groups)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF_L, color = Three_Groups), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#right SLF FD (violin plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF_R)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)

#right SLF FD (raincloud plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF_R, fill = Three_Groups)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF_R, color = Three_Groups), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#for FC - mean
#whole SLF FC (violin plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FC_SLF)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)

#left SLF FC (violin plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FC_SLF_L)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)

#right SLF FC (violin plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FC_SLF_R)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)

#for FDC - mean
#whole SLF FDC (violin plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)

#Left SLF FDC (violin plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF_L)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)

#Right SLF FDC (violin plot)
ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF_R)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)





