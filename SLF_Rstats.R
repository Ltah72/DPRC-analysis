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
pacman::p_load(dplyr, ggplot2, psych, car, multcomp)


#first read in the covariates group data file: 
setwd("V:/PartInfo")
covariates_data <- read.csv("covariates-participants-lined-up.csv")
#convert the group data to a factor variable
covariates_data$Group <- as.factor(covariates_data$Group)


#navigate to the correct pathway: 
setwd("V:/NECTAR_data/LENORE/derivatives/diff_data/test/metric_results")


#read in text file of data
SLF_data <- cbind.data.frame(read.table("FD_SLF_TOI.txt", header = T), read.table("FC_log_SLF_TOI.txt", header = T), read.table("FDC_SLF_TOI.txt", header = T), read.table("FD_SLF_L_TOI.txt", header = T), read.table("FC_log_SLF_L_TOI.txt", header = T), read.table("FDC_SLF_L_TOI.txt", header = T), read.table("FD_SLF_R_TOI.txt", header = T), read.table("FC_log_SLF_R_TOI.txt", header = T), read.table("FDC_SLF_R_TOI.txt", header = T))
#rename columns with each FBA metric
colnames(SLF_data) <- c("mn_FD_SLF", "md_FD_SLF", "std_FD_SLF", "std_rv_FD_SLF", "min_FD_SLF", "max_FD_SLF", "count_FD_SLF", "mn_FC_SLF", "md_FC_SLF", "std_FC_SLF", "std_rv_FC_SLF", "min_FC_SLF", "max_FC_SLF", "count_FC_SLF", "mn_FDC_SLF", "md_FDC_SLF", "std_FDC_SLF", "std_rv_FDC_SLF", "min_FDC_SLF", "max_FDC_SLF", "count_FDC_SLF", "mn_FD_SLF_L", "md_FD_SLF_L", "std_FD_SLF_L", "std_rv_FD_SLF_L", "min_FD_SLF_L", "max_FD_SLF_L", "count_FD_SLF_L", "mn_FC_SLF_L", "md_FC_SLF_L", "std_FC_SLF_L", "std_rv_FC_SLF_L", "min_FC_SLF_L", "max_FC_SLF_L", "count_FC_SLF_L", "mn_FDC_SLF_L", "md_FDC_SLF_L", "std_FDC_SLF_L", "std_rv_FDC_SLF_L", "min_FDC_SLF_L", "max_FDC_SLF_L", "count_FDC_SLF_L", "mn_FD_SLF_R", "md_FD_SLF_R", "std_FD_SLF_R", "std_rv_FD_SLF_R", "min_FD_SLF_R", "max_FD_SLF_R", "count_FD_SLF_R", "mn_FC_SLF_R", "md_FC_SLF_R", "std_FC_SLF_R", "std_rv_FC_SLF_R", "min_FC_SLF_R", "max_FC_SLF_R", "count_FC_SLF_R", "mn_FDC_SLF_R", "md_FDC_SLF_R", "std_FDC_SLF_R", "std_rv_FDC_SLF_R", "min_FDC_SLF_R", "max_FDC_SLF_R", "count_FDC_SLF_R")
#add in Group classification column for each participant - data from another worksheet
SLF_data$Group <- covariates_data$Group


#look at descriptives of the SLF FBA metrices between groups
SLF_FD_descrip <- describeBy(SLF_data$mn_FD_SLF, SLF_data$Group)
SLF_FC_descrip <- describeBy(SLF_data$mn_FC_SLF, SLF_data$Group)
SLF_FDC_descrip <- describeBy(SLF_data$mn_FDC_SLF, SLF_data$Group)


#run ANOVA to see if there are significant differences between groups
#for FD
#mean
SLF_FD_mod <- lm(mn_FD ~ Group, data = SLF_data)
#median
SLF_FD_mod <- lm(md_FD ~ Group, data = SLF_data)
#run ANOVA
anova(SLF_FD_mod)
#to get multiple regression outputs
summary(SLF_FD_mod)
#run pairwise comparisons, given that the F-test was significant. 
post_hoc_SLF_FD_mod <- glht(SLF_FD_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_mod)
confint(post_hoc_SLF_FD_mod)

#for FC
#mean
SLF_FC_mod <- lm(mn_FC ~ Group, data = SLF_data)
#median
SLF_FC_mod <- lm(md_FC ~ Group, data = SLF_data)
anova(SLF_FC_mod)

#for FDC
#mean
SLF_FDC_mod <- lm(mn_FDC ~ Group, data = SLF_data)
#median
SLF_FDC_mod <- lm(md_FDC ~ Group, data = SLF_data)
anova(SLF_FDC_mod)


#plot data
#FD - mean
ggplot(SLF_data, aes(x = Group, y = mn_FD)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group_status") + 
    ylab("FD") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")
#FD - median
ggplot(SLF_data, aes(x = Group, y = md_FD)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group_status") + 
    ylab("FD") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")

#FC - mean
ggplot(SLF_data, aes(x = Group, y = mn_FC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group_status") + 
    ylab("FC") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")
#FC - median
ggplot(SLF_data, aes(x = Group, y = md_FC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group_status") + 
    ylab("FC") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")

#FDC - mean
ggplot(SLF_data, aes(x = Group, y = mn_FDC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group_status") + 
    ylab("FDC") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")
#FDC - median
ggplot(SLF_data, aes(x = Group, y = md_FDC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group_status") + 
    ylab("FDC") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")





















