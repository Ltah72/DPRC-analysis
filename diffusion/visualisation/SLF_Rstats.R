#Perform statistics on tract of interest (TOI) from diffusion analysis. This script will perform analysis for group 
#differences in the TOI which was selected from the user in the previous pipelines, done mainly from the software
#programme, MRtrix3. The main outputs are the text files of the fixel-based analysis (FBA) metrics from the TOI.m pipeline. 
#For this script, I have chosen to analyse the superior longitudinal fasciculus (SLF) between the 5 participant groups of 
#interest - Controls (C), subjective cognitive decline (SCD), amnestic mild cognitive impairment (aMCI), multiple-domain
#mild cognitive impairment (mMCI), and Alzheimer's disease (AD). 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 8/12/20


#------------------------------Setting up--------------------------------------#
#install packages/open libraries
pacman::p_load(dplyr, ggplot2, psych, car, multcomp, lsr, tidyr, BayesFactor, tidyverse, ppcor, nlme, effectsize, rstatix, sjstats, purrr, corrplot)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#first read in the covariates group data file: 
#setwd('/yourpathway/')
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')


####-----------------------Cross-sectional analysis-------------------------------####



DPRC_neuropsych_data <- read.csv("cross-sectional_DPRC_neuropsych_data_lined_up_valid_participants.csv")
#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'

#convert variables
DPRC_neuropsych_data$ParticipantID <- as.factor(DPRC_neuropsych_data$ParticipantID)
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Sex_binary <- as.factor(DPRC_neuropsych_data$Sex_binary)
DPRC_neuropsych_data$Sex<- as.factor(DPRC_neuropsych_data$Sex)

#navigate to the correct pathway which contains the SLF metric text files: 
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/template/TOI')

#read in text file of data
SLF_data <- cbind.data.frame(read.table("FD_SLF_whole_TOI.txt", header = T), 
                             read.table("FC_log_SLF_whole_TOI.txt", header = T), 
                             read.table("FDC_SLF_whole_TOI.txt", header = T), 
                             read.table("FD_SLF_L_TOI.txt", header = T), 
                             read.table("FC_log_SLF_L_TOI.txt", header = T), 
                             read.table("FDC_SLF_L_TOI.txt", header = T), 
                             read.table("FD_SLF_R_TOI.txt", header = T), 
                             read.table("FC_log_SLF_R_TOI.txt", header = T), 
                             read.table("FDC_SLF_R_TOI.txt", header = T), 
                             read.table("FD_SLF1_L_TOI.txt", header = T), 
                             read.table("FC_log_SLF1_L_TOI.txt", header = T), 
                             read.table("FDC_SLF1_L_TOI.txt", header = T), 
                             read.table("FD_SLF2_L_TOI.txt", header = T), 
                             read.table("FC_log_SLF2_L_TOI.txt", header = T), 
                             read.table("FDC_SLF2_L_TOI.txt", header = T), 
                             read.table("FD_SLF3_L_TOI.txt", header = T), 
                             read.table("FC_log_SLF3_L_TOI.txt", header = T), 
                             read.table("FDC_SLF3_L_TOI.txt", header = T), 
                             read.table("FD_SLF1_R_TOI.txt", header = T), 
                             read.table("FC_log_SLF1_R_TOI.txt", header = T), 
                             read.table("FDC_SLF1_R_TOI.txt", header = T), 
                             read.table("FD_SLF2_R_TOI.txt", header = T), 
                             read.table("FC_log_SLF2_R_TOI.txt", header = T), 
                             read.table("FDC_SLF2_R_TOI.txt", header = T), 
                             read.table("FD_SLF3_R_TOI.txt", header = T), 
                             read.table("FC_log_SLF3_R_TOI.txt", header = T), 
                             read.table("FDC_SLF3_R_TOI.txt", header = T))
#rename columns with each FBA metric
colnames(SLF_data) <- c("mn_FD_SLF", "md_FD_SLF", "std_FD_SLF", "std_rv_FD_SLF", "min_FD_SLF", "max_FD_SLF", "count_FD_SLF", 
                        "mn_FC_SLF", "md_FC_SLF", "std_FC_SLF", "std_rv_FC_SLF", "min_FC_SLF", "max_FC_SLF", "count_FC_SLF", 
                        "mn_FDC_SLF", "md_FDC_SLF", "std_FDC_SLF", "std_rv_FDC_SLF", "min_FDC_SLF", "max_FDC_SLF", "count_FDC_SLF", 
                        "mn_FD_SLF_L", "md_FD_SLF_L", "std_FD_SLF_L", "std_rv_FD_SLF_L", "min_FD_SLF_L", "max_FD_SLF_L", 
                        "count_FD_SLF_L", "mn_FC_SLF_L", "md_FC_SLF_L", "std_FC_SLF_L", "std_rv_FC_SLF_L", "min_FC_SLF_L", 
                        "max_FC_SLF_L", "count_FC_SLF_L", "mn_FDC_SLF_L", "md_FDC_SLF_L", "std_FDC_SLF_L", "std_rv_FDC_SLF_L", 
                        "min_FDC_SLF_L", "max_FDC_SLF_L", "count_FDC_SLF_L", "mn_FD_SLF_R", "md_FD_SLF_R", "std_FD_SLF_R", 
                        "std_rv_FD_SLF_R", "min_FD_SLF_R", "max_FD_SLF_R", "count_FD_SLF_R","mn_FC_SLF_R", "md_FC_SLF_R", 
                        "std_FC_SLF_R", "std_rv_FC_SLF_R", "min_FC_SLF_R", "max_FC_SLF_R", "count_FC_SLF_R", "mn_FDC_SLF_R", 
                        "md_FDC_SLF_R", "std_FDC_SLF_R", "std_rv_FDC_SLF_R", "min_FDC_SLF_R", "max_FDC_SLF_R", "count_FDC_SLF_R", 
                        "mn_FD_SLF1_L", "md_FD_SLF1_L", "std_FD_SLF1_L", "std_rv_FD_SLF1_L", "min_FD_SLF1_L", "max_FD_SLF1_L", 
                        "count_FD_SLF1_L", "mn_FC_SLF1_L", "md_FC_SLF1_L", "std_FC_SLF1_L", "std_rv_FC_SLF1_L", "min_FC_SLF1_L", 
                        "max_FC_SLF1_L", "count_FC_SLF1_L", "mn_FDC_SLF1_L", "md_FDC_SLF1_L", "std_FDC_SLF1_L", "std_rv_FDC_SLF1_L", 
                        "min_FDC_SLF1_L", "max_FDC_SLF1_L", "count_FDC_SLF1_L", "mn_FD_SLF2_L", "md_FD_SLF2_L", "std_FD_SLF2_L", 
                        "std_rv_FD_SLF2_L", "min_FD_SLF2_L", "max_FD_SLF2_L", "count_FD_SLF2_L", "mn_FC_SLF2_L", "md_FC_SLF2_L", 
                        "std_FC_SLF2_L", "std_rv_FC_SLF2_L", "min_FC_SLF2_L", "max_FC_SLF2_L", "count_FC_SLF2_L", "mn_FDC_SLF2_L", 
                        "md_FDC_SLF2_L", "std_FDC_SLF2_L", "std_rv_FDC_SLF2_L", "min_FDC_SLF2_L", "max_FDC_SLF2_L", "count_FDC_SLF2_L",
                        "mn_FD_SLF3_L", "md_FD_SLF3_L", "std_FD_SLF3_L", "std_rv_FD_SLF3_L", "min_FD_SLF3_L", "max_FD_SLF3_L", 
                        "count_FD_SLF3_L", "mn_FC_SLF3_L", "md_FC_SLF3_L", "std_FC_SLF3_L", "std_rv_FC_SLF3_L", "min_Fc_SLF3_L", 
                        "max_FC_SLF3_L", "count_FC_SLF3_L", "mn_FDC_SLF3_L", "md_FDC_SLF3_L", "std_FDC_SLF3_L", "std_rv_FDC_SLF3_L", 
                        "min_FDC_SLF3_L", "max_FDC_SLF3_L", "count_FDC_SLF3_L", "mn_FD_SLF1_R", "md_FD_SLF1_R", "std_FD_SLF1_R", 
                        "std_rv_FD_SLF1_R", "min_FD_SLF1_R", "max_FD_SLF1_R", "count_FD_SLF1_R", "mn_FC_SLF1_R", "md_FC_SLF1_R", 
                        "std_FC_SLF1_R", "std_rv_FC_SLF1_R", "min_FC_SLF1_R", "max_FC_SLF1_R", "count_FC_SLF1_R", "mn_FDC_SLF1_R", 
                        "md_FDC_SLF1_R", "std_FDC_SLF1_R", "std_rv_FDC_SLF1_R", "min_FDC_SLF1_R", "max_FDC_SLF1_R", "count_FDC_SLF1_R", 
                        "mn_FD_SLF2_R", "md_FD_SLF2_R", "std_FD_SLF2_R", "std_rv_FD_SLF2_R", "min_FD_SLF2_R", "max_FD_SLF2_R", 
                        "count_FD_SLF2_R","mn_FC_SLF2_R", "md_FC_SLF2_R", "std_FC_SLF2_R", "std_rv_FC_SLF2_R", "min_FC_SLF2_R", 
                        "max_FC_SLF2_R", "count_FC_SLF2_R", "mn_FDC_SLF2_R", "md_FDC_SLF2_R", "std_FDC_SLF2_R", "std_rv_FDC_SLF2_R", 
                        "min_FDC_SLF2_R", "max_FDC_SLF2_R", "count_FDC_SLF2_R", "mn_FD_SLF3_R", "md_FD_SLF3_R", "std_FD_SLF3_R", 
                        "std_rv_FD_SLF3_R", "min_FD_SLF3_R", "max_FD_SLF3_R", "count_FD_SLF3_R","mn_FC_SLF3_R", "md_FC_SLF3_R", 
                        "std_FC_SLF3_R", "std_rv_FC_SLF3_R", "min_FC_SLF3_R", "max_FC_SLF3_R", "count_FC_SLF3_R", "mn_FDC_SLF3_R",
                        "md_FDC_SLF3_R", "std_FDC_SLF3_R", "std_rv_FDC_SLF3_R", "min_FDC_SLF3_R", "max_FDC_SLF3_R", "count_FDC_SLF3_R")
#add in Participant ID  - data from another worksheet
SLF_data$ParticipantID <- DPRC_neuropsych_data$ParticipantID
#add in Group classification column for each participant - data from another worksheet
SLF_data$Group <- DPRC_neuropsych_data$Group
#add in covariates (clinical site) from the covariates_data dataframe to the SLF dataframe
SLF_data$ClinSite_name <- DPRC_neuropsych_data$Clinical_site 
SLF_data$Age <- DPRC_neuropsych_data$Age 
SLF_data$Sex <- DPRC_neuropsych_data$Sex
SLF_data$Sex_binary <- DPRC_neuropsych_data$Sex_binary
#add in neuropsych variables
neuropsych_test_names <- c("TrailsA.Raw","TrailsA.Z","TrailsB.Raw","TrailsB.Z","ColorNaming.Raw","ColorNaming.Z","WordReading.Raw","WordReading.Z","Inhibition.Raw","Inhibition.Z","LetFluency.Raw","LetFluency.Z","CatFluency.Raw","CatFluency.Z","Switching.Raw","Switching.z","HayBTime1.Raw","HayBTime1.z","HayBTime2.Raw","HayBTime2.z","HayBCatA.Raw","HayBCatA.z","HayBCatB.Raw","HayBCatB.z")
SLF_data[neuropsych_test_names] <- DPRC_neuropsych_data[neuropsych_test_names]
#add in trend Group variable
Trend_Group <- as.numeric(DPRC_neuropsych_data$Group)

####--------------------------Descriptives----------------------------------####
#look at descriptive stats of the whole SLF FBA metrics between groups
SLF_FD_descrip <- describeBy(SLF_data$mn_FD_SLF, SLF_data$Group)
SLF_FC_descrip <- describeBy(SLF_data$mn_FC_SLF, SLF_data$Group)
SLF_FDC_descrip <- describeBy(SLF_data$mn_FDC_SLF, SLF_data$Group)
mean(SLF_data$mn_FD_SLF)
sd(SLF_data$mn_FD_SLF)
mean(SLF_data$mn_FC_SLF)
sd(SLF_data$mn_FC_SLF)
mean(SLF_data$mn_FDC_SLF)
sd(SLF_data$mn_FDC_SLF)
#look at descriptive stats of the Left SLF FBA metrics between groups
SLF_FD_L_descrip <- describeBy(SLF_data$mn_FD_SLF_L, SLF_data$Group)
SLF_FC_L_descrip <- describeBy(SLF_data$mn_FC_SLF_L, SLF_data$Group)
SLF_FDC_L_descrip <- describeBy(SLF_data$mn_FDC_SLF_L, SLF_data$Group)
mean(SLF_data$mn_FD_SLF_L)
sd(SLF_data$mn_FD_SLF_L)
mean(SLF_data$mn_FC_SLF_L)
sd(SLF_data$mn_FC_SLF_L)
mean(SLF_data$mn_FDC_SLF_L)
sd(SLF_data$mn_FDC_SLF_L)
#look at descriptive stats of the Right SLF FBA metrics between groups
SLF_FD_R_descrip <- describeBy(SLF_data$mn_FD_SLF_R, SLF_data$Group)
SLF_FC_R_descrip <- describeBy(SLF_data$mn_FC_SLF_R, SLF_data$Group)
SLF_FDC_R_descrip <- describeBy(SLF_data$mn_FDC_SLF_R, SLF_data$Group)
mean(SLF_data$mn_FD_SLF_R)
sd(SLF_data$mn_FD_SLF_R)
mean(SLF_data$mn_FC_SLF_R)
sd(SLF_data$mn_FC_SLF_R)
mean(SLF_data$mn_FDC_SLF_R)
sd(SLF_data$mn_FDC_SLF_R)
#look at descriptive stats of the Left SLF1 FBA metrics between groups
SLF1_FD_L_descrip <- describeBy(SLF_data$mn_FD_SLF1_L, SLF_data$Group)
SLF1_FC_L_descrip <- describeBy(SLF_data$mn_FC_SLF1_L, SLF_data$Group)
SLF1_FDC_L_descrip <- describeBy(SLF_data$mn_FDC_SLF1_L, SLF_data$Group)
mean(SLF_data$mn_FD_SLF1_L)
sd(SLF_data$mn_FD_SLF1_L)
mean(SLF_data$mn_FC_SLF1_L)
sd(SLF_data$mn_FC_SLF1_L)
mean(SLF_data$mn_FDC_SLF1_L)
sd(SLF_data$mn_FDC_SLF1_L)
#look at descriptive stats of the Left SLF2 FBA metrics between groups
SLF2_FD_L_descrip <- describeBy(SLF_data$mn_FD_SLF2_L, SLF_data$Group)
SLF2_FC_L_descrip <- describeBy(SLF_data$mn_FC_SLF2_L, SLF_data$Group)
SLF2_FDC_L_descrip <- describeBy(SLF_data$mn_FDC_SLF2_L, SLF_data$Group)
mean(SLF_data$mn_FD_SLF2_L)
sd(SLF_data$mn_FD_SLF2_L)
mean(SLF_data$mn_FC_SLF2_L)
sd(SLF_data$mn_FC_SLF2_L)
mean(SLF_data$mn_FDC_SLF2_L)
sd(SLF_data$mn_FDC_SLF2_L)
#look at descriptive stats of the Left SLF2 FBA metrics between groups
SLF3_FD_L_descrip <- describeBy(SLF_data$mn_FD_SLF3_L, SLF_data$Group)
SLF3_FC_L_descrip <- describeBy(SLF_data$mn_FC_SLF3_L, SLF_data$Group)
SLF3_FDC_L_descrip <- describeBy(SLF_data$mn_FDC_SLF3_L, SLF_data$Group)
mean(SLF_data$mn_FD_SLF3_L)
sd(SLF_data$mn_FD_SLF3_L)
mean(SLF_data$mn_FC_SLF3_L)
sd(SLF_data$mn_FC_SLF3_L)
mean(SLF_data$mn_FDC_SLF3_L)
sd(SLF_data$mn_FDC_SLF3_L)
#look at descriptive stats of the Right SLF1 FBA metrics between groups
SLF1_FD_R_descrip <- describeBy(SLF_data$mn_FD_SLF1_R, SLF_data$Group)
SLF1_FC_R_descrip <- describeBy(SLF_data$mn_FC_SLF1_R, SLF_data$Group)
SLF1_FDC_R_descrip <- describeBy(SLF_data$mn_FDC_SLF1_R, SLF_data$Group)
mean(SLF_data$mn_FD_SLF1_R)
sd(SLF_data$mn_FD_SLF1_R)
mean(SLF_data$mn_FC_SLF1_R)
sd(SLF_data$mn_FC_SLF1_R)
mean(SLF_data$mn_FDC_SLF1_R)
sd(SLF_data$mn_FDC_SLF1_R)
#look at descriptive stats of the Right SLF2 FBA metrics between groups
SLF2_FD_R_descrip <- describeBy(SLF_data$mn_FD_SLF2_R, SLF_data$Group)
SLF2_FC_R_descrip <- describeBy(SLF_data$mn_FC_SLF2_R, SLF_data$Group)
SLF2_FDC_R_descrip <- describeBy(SLF_data$mn_FDC_SLF2_R, SLF_data$Group)
mean(SLF_data$mn_FD_SLF2_R)
sd(SLF_data$mn_FD_SLF2_R)
mean(SLF_data$mn_FC_SLF2_R)
sd(SLF_data$mn_FC_SLF2_R)
mean(SLF_data$mn_FDC_SLF2_R)
sd(SLF_data$mn_FDC_SLF2_R)
#look at descriptive stats of the Right SLF2 FBA metrics between groups
SLF3_FD_R_descrip <- describeBy(SLF_data$mn_FD_SLF3_R, SLF_data$Group)
SLF3_FC_R_descrip <- describeBy(SLF_data$mn_FC_SLF3_R, SLF_data$Group)
SLF3_FDC_R_descrip <- describeBy(SLF_data$mn_FDC_SLF3_R, SLF_data$Group)
mean(SLF_data$mn_FD_SLF3_R)
sd(SLF_data$mn_FD_SLF3_R)
mean(SLF_data$mn_FC_SLF3_R)
sd(SLF_data$mn_FC_SLF3_R)
mean(SLF_data$mn_FDC_SLF3_R)
sd(SLF_data$mn_FDC_SLF3_R)



###------------------ ANOVA & ANCOVA (age + sex) ---------------------------####

#run ANOVA to see if there are significant differences between groups
#for FD
#mean
SLF_FD_mod <- lm(mn_FD_SLF ~ Group, data = SLF_data)
SLF_FD_mod_L <- lm(mn_FD_SLF_L ~ Group, data = SLF_data)
SLF_FD_mod_R <- lm(mn_FD_SLF_R ~ Group, data = SLF_data)
SLF1_FD_mod_L <- lm(mn_FD_SLF1_L ~ Group, data = SLF_data)
SLF2_FD_mod_L <- lm(mn_FD_SLF2_L ~ Group, data = SLF_data)
SLF3_FD_mod_L <- lm(mn_FD_SLF3_L ~ Group, data = SLF_data)
SLF1_FD_mod_R <- lm(mn_FD_SLF1_R ~ Group, data = SLF_data)
SLF2_FD_mod_R <- lm(mn_FD_SLF2_R ~ Group, data = SLF_data)
SLF3_FD_mod_R <- lm(mn_FD_SLF3_R ~ Group, data = SLF_data)

#include the covariate of age (run an ANCOVA) in model
# SLF_FD_age_mod <- lm(mn_FD_SLF ~ Group + Age, data = SLF_data)
# SLF_FD_age_mod_L <- lm(mn_FD_SLF_L ~ Group + Age, data = SLF_data)
# SLF_FD_age_mod_R <- lm(mn_FD_SLF_R ~ Group + Age, data = SLF_data)

#include the covariate of clinical site (run an ANCOVA) in model
# SLF_FD_clinsite_mod <- lm(mn_FD_SLF ~ Group + ClinSite_name, data = SLF_data)
# SLF_FD_clinsite_mod_L <- lm(mn_FD_SLF_L ~ Group + ClinSite_name, data = SLF_data)
# SLF_FD_clinsite_mod_R <- lm(mn_FD_SLF_R ~ Group + ClinSite_name, data = SLF_data)

#include the covariate of 2 variables (age + sex) (run an ANCOVA) in model
SLF_FD_2covar_mod <- lm(mn_FD_SLF ~ Group + Age + Sex, data = SLF_data)
SLF_FD_2covar_mod_L <- lm(mn_FD_SLF_L ~ Group + Age + Sex, data = SLF_data)
SLF_FD_2covar_mod_R <- lm(mn_FD_SLF_R ~ Group + Age + Sex, data = SLF_data)
SLF1_FD_2covar_mod_L <- lm(mn_FD_SLF1_L ~ Group + Age + Sex, data = SLF_data)
SLF2_FD_2covar_mod_L <- lm(mn_FD_SLF2_L ~ Group + Age + Sex, data = SLF_data)
SLF3_FD_2covar_mod_L <- lm(mn_FD_SLF3_L ~ Group + Age + Sex, data = SLF_data)
SLF1_FD_2covar_mod_R <- lm(mn_FD_SLF1_R ~ Group + Age + Sex, data = SLF_data)
SLF2_FD_2covar_mod_R <- lm(mn_FD_SLF2_R ~ Group + Age + Sex, data = SLF_data)
SLF3_FD_2covar_mod_R <- lm(mn_FD_SLF3_R ~ Group + Age + Sex, data = SLF_data)

#test to see if the covariates and the treatment variable are independent from one another (assumption). 
#for whole FD
# ancova_sex_FD_check <- aov(mn_FD_SLF ~ Sex, data = SLF_data) 
# summary(ancova_sex_FD_check)
# ancova_age_FD_check <- aov(mn_FD_SLF ~ Age, data = SLF_data) 
# summary(ancova_age_FD_check)
# #for whole FC
# ancova_sex_FC_check <- aov(mn_FC_SLF ~ Sex, data = SLF_data) 
# summary(ancova_sex_FC_check)
# ancova_age_FC_check <- aov(mn_FC_SLF ~ Age, data = SLF_data) 
# summary(ancova_age_FC_check)
# #for whole FDC
# ancova_sex_FDC_check <- aov(mn_FDC_SLF ~ Sex, data = SLF_data) 
# summary(ancova_sex_FDC_check)
# ancova_age_FDC_check <- aov(mn_FDC_SLF ~ Age, data = SLF_data) 
# summary(ancova_age_FDC_check)

#median
#SLF_FD_mod <- lm(md_FD ~ Group, data = SLF_data)
#run ANOVA
anova(SLF_FD_mod)
anova(SLF_FD_mod_L)
anova(SLF_FD_mod_R)
anova(SLF1_FD_mod_L)
anova(SLF2_FD_mod_L)
anova(SLF3_FD_mod_L)
anova(SLF1_FD_mod_R)
anova(SLF2_FD_mod_R)
anova(SLF3_FD_mod_R)

#run Bayesian ANOVA
anovaBF(mn_FD_SLF ~ Group, data = SLF_data) 
anovaBF(mn_FD_SLF_L ~ Group, data = SLF_data) 
anovaBF(mn_FD_SLF_R ~ Group, data = SLF_data) 
#run ANCOVA
# anova(SLF_FD_age_mod)
# anova(SLF_FD_age_mod_L)
# anova(SLF_FD_age_mod_R)
# anova(SLF_FD_clinsite_mod)
# anova(SLF_FD_clinsite_mod_L)
# anova(SLF_FD_clinsite_mod_R)
#ANCOVA - for age and sex
anova(SLF_FD_2covar_mod)
anova(SLF_FD_2covar_mod_L)
anova(SLF_FD_2covar_mod_R)
anova(SLF1_FD_2covar_mod_L)
anova(SLF2_FD_2covar_mod_L)
anova(SLF3_FD_2covar_mod_L)
anova(SLF1_FD_2covar_mod_R)
anova(SLF2_FD_2covar_mod_R)
anova(SLF3_FD_2covar_mod_R)
#run Bayesian ANCOVA
#anovaBF(mn_FD_SLF ~ Group + Age, data = SLF_data) 
#anovaBF(mn_FD_SLF_L ~ Group, data = SLF_data) 
#anovaBF(mn_FD_SLF_R ~ Group, data = SLF_data) 
#to get multiple regression outputs
summary(SLF_FD_mod)
summary(SLF_FD_mod_L)
summary(SLF_FD_mod_R)

#Linear Trend Analysis in Linear Regression
SLF_FD_LinTrend_mod <- lm(mn_FD_SLF ~ Trend_Group + Group, data = SLF_data)
anova(SLF_FD_LinTrend_mod)
summary(SLF_FD_LinTrend_mod) 
SLF_L_FD_LinTrend_mod <- lm(mn_FD_SLF_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF_L_FD_LinTrend_mod)
summary(SLF_L_FD_LinTrend_mod) 
SLF_R_FD_LinTrend_mod <- lm(mn_FD_SLF_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF_R_FD_LinTrend_mod)
summary(SLF_R_FD_LinTrend_mod) 
SLF1_L_FD_LinTrend_mod <- lm(mn_FD_SLF1_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF1_L_FD_LinTrend_mod)
summary(SLF1_L_FD_LinTrend_mod) 
SLF2_L_FD_LinTrend_mod <- lm(mn_FD_SLF2_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF2_L_FD_LinTrend_mod)
summary(SLF2_L_FD_LinTrend_mod) 
SLF3_L_FD_LinTrend_mod <- lm(mn_FD_SLF3_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF3_L_FD_LinTrend_mod)
summary(SLF3_L_FD_LinTrend_mod)
SLF1_R_FD_LinTrend_mod <- lm(mn_FD_SLF1_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF1_R_FD_LinTrend_mod)
summary(SLF1_R_FD_LinTrend_mod) 
SLF2_R_FD_LinTrend_mod <- lm(mn_FD_SLF2_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF2_R_FD_LinTrend_mod)
summary(SLF2_R_FD_LinTrend_mod) 
SLF3_R_FD_LinTrend_mod <- lm(mn_FD_SLF3_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF3_R_FD_LinTrend_mod)
summary(SLF3_R_FD_LinTrend_mod) 

#Linear Trend Analysis in Linear Regression w/ covariates age + sex
SLF_FD_LinTrend_2covar_mod <- lm(mn_FD_SLF ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF_FD_LinTrend_2covar_mod)
summary(SLF_FD_LinTrend_2covar_mod) 
SLF_L_FD_LinTrend_2covar_mod <- lm(mn_FD_SLF_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF_L_FD_LinTrend_2covar_mod)
summary(SLF_L_FD_LinTrend_2covar_mod) 
SLF_R_FD_LinTrend_2covar_mod <- lm(mn_FD_SLF_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF_R_FD_LinTrend_2covar_mod)
summary(SLF_R_FD_LinTrend_2covar_mod) 
SLF1_L_FD_LinTrend_2covar_mod <- lm(mn_FD_SLF1_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF1_L_FD_LinTrend_2covar_mod)
summary(SLF1_L_FD_LinTrend_2covar_mod) 
SLF2_L_FD_LinTrend_2covar_mod <- lm(mn_FD_SLF2_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF2_L_FD_LinTrend_2covar_mod)
summary(SLF2_L_FD_LinTrend_2covar_mod) 
SLF3_L_FD_LinTrend_2covar_mod <- lm(mn_FD_SLF3_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF3_L_FD_LinTrend_2covar_mod)
summary(SLF3_L_FD_LinTrend_2covar_mod)
SLF1_R_FD_LinTrend_2covar_mod <- lm(mn_FD_SLF1_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF1_R_FD_LinTrend_2covar_mod)
summary(SLF1_R_FD_LinTrend_2covar_mod) 
SLF2_R_FD_LinTrend_2covar_mod <- lm(mn_FD_SLF2_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF2_R_FD_LinTrend_2covar_mod)
summary(SLF2_R_FD_LinTrend_2covar_mod) 
SLF3_R_FD_LinTrend_2covar_mod <- lm(mn_FD_SLF3_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF3_R_FD_LinTrend_2covar_mod)
summary(SLF3_R_FD_LinTrend_2covar_mod) 

#for ANCOVA
# summary(SLF_FD_age_mod)
# summary(SLF_FD_age_mod_L)
# summary(SLF_FD_age_mod_R)
# summary(SLF_FD_clinsite_mod)
# summary(SLF_FD_clinsite_mod_L)
# summary(SLF_FD_clinsite_mod_R)
#run pairwise comparisons, given that the F-test was significant. 
#whole SLF
post_hoc_SLF_FD_mod <- glht(SLF_FD_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_mod)
confint(post_hoc_SLF_FD_mod)
#Left and Right SLF whole
post_hoc_SLF_FD_mod_L <- glht(SLF_FD_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_mod_L)
confint(post_hoc_SLF_FD_mod_L)
post_hoc_SLF_FD_mod_R <- glht(SLF_FD_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_mod_R)
confint(post_hoc_SLF_FD_mod_R)
#Left SLF 1, 2, and 3
post_hoc_SLF1_FD_mod_L <- glht(SLF1_FD_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF1_FD_mod_L)
confint(post_hoc_SLF1_FD_mod_L)
post_hoc_SLF2_FD_mod_L <- glht(SLF2_FD_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_FD_mod_L)
confint(post_hoc_SLF2_FD_mod_L)
post_hoc_SLF3_FD_mod_L <- glht(SLF3_FD_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_FD_mod_L)
confint(post_hoc_SLF3_FD_mod_L)
#Right SLF 1, 2 and 3
post_hoc_SLF1_FD_mod_R <- glht(SLF1_FD_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF1_FD_mod_R)
confint(post_hoc_SLF1_FD_mod_R)
post_hoc_SLF2_FD_mod_R <- glht(SLF2_FD_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_FD_mod_R)
confint(post_hoc_SLF2_FD_mod_R)
post_hoc_SLF3_FD_mod_R <- glht(SLF3_FD_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_FD_mod_R)
confint(post_hoc_SLF3_FD_mod_R)

#post hoc for covariate (age & sex)
post_hoc_SLF_FD_ancova_mod <- glht(SLF_FD_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_ancova_mod)
confint(post_hoc_SLF_FD_ancova_mod)

post_hoc_SLF_L_FD_ancova_mod <- glht(SLF_FD_2covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_L_FD_ancova_mod)
confint(post_hoc_SLF_L_FD_ancova_mod)

post_hoc_SLF_R_FD_ancova_mod <- glht(SLF_FD_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_R_FD_ancova_mod)
confint(post_hoc_SLF_R_FD_ancova_mod)

post_hoc_SLF2_L_FD_ancova_mod <- glht(SLF2_FD_2covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_L_FD_ancova_mod)
confint(post_hoc_SLF2_L_FD_ancova_mod)

post_hoc_SLF3_L_FD_ancova_mod <- glht(SLF3_FD_2covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_L_FD_ancova_mod)
confint(post_hoc_SLF3_L_FD_ancova_mod)

post_hoc_SLF1_R_FD_ancova_mod <- glht(SLF1_FD_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF1_R_FD_ancova_mod)
confint(post_hoc_SLF1_R_FD_ancova_mod)

post_hoc_SLF2_R_FD_ancova_mod <- glht(SLF2_FD_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_R_FD_ancova_mod)
confint(post_hoc_SLF2_R_FD_ancova_mod)

post_hoc_SLF3_R_FD_ancova_mod <- glht(SLF3_FD_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_R_FD_ancova_mod)
confint(post_hoc_SLF3_R_FD_ancova_mod)

#plot 95% Confidence Interval
SLF_FD_95CI_data <- data.frame(SLF_group_number = c('1','1','1','2','2','3','3','3','4','4','4','5','5','5','5','6','7','7','7','8','8','8','8'),
                                      SLF_type = c('Whole_SLF', 'Whole_SLF','Whole_SLF','Left_SLF','Left_SLF', 'Right_SLF', 'Right_SLF', 'Right_SLF','Left_SLF2','Left_SLF2','Left_SLF2','Left_SLF3','Left_SLF3','Left_SLF3','Left_SLF3','Right_SLF1','Right_SLF2','Right_SLF2','Right_SLF2','Right_SLF3', 'Right_SLF3', 'Right_SLF3', 'Right_SLF3'),
                                      Group_contrast = c('C > aMCI', 'C > mMCI', 'C > AD', 'C > mMCI', 'C > AD','C > aMCI','C > mMCI','C > AD','C > mMCI','C > AD','SCD > AD','C > SCD','C > aMCI','C > mMCI','C > AD','C > AD','C > aMCI','C > mMCI','C > AD','C > SCD','C > aMCI','C > mMCI','C > AD'),
                                      estimate_diff = c(0.0152899,0.0191437,0.0252296,0.0179369,0.0249809,0.0165314,0.0202112,0.0254493,0.020646,0.032745,0.022266,0.0149789,0.0155556,0.0198357,0.0251135,0.0282489,0.019994,0.0245799,0.0309866,0.0155239,0.0163904,0.0214478,0.0194152),
                                      lower = c(0.029121,0.033129,0.0416139,0.0325738,0.0421285,0.0309934,0.0348344,0.0425811,0.040536,0.056048,0.043348,0.0280131,0.0288062,0.0332341,0.0408103,0.0490901,0.0394673,0.0442704,0.0540548,0.0293862,0.0304829,0.0356975,0.0361093), 
                                      upper = c(0.0014589,0.0051583,0.0088453,0.0033001,0.0078333,0.0020695,0.0055879,0.0083176,0.000755,0.009442,0.001183,0.0019447,0.0023049,0.0064372,0.0094167,0.0074076,0.0005208,0.0048893,0.0079185,0.0016616,0.0022978,0.0071981,0.0027211))  
     
#plot data
ggplot(SLF_FD_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Whole SLF", "2" = "Left SLF", "3" = "Right SLF", "4" = "Right SLF 1", "5" = "Left SLF 2", "6" = "Left SLF 3", "7" = "Right SLF 2", "8" = "Right SLF 3"))+
    theme_classic()

#only include the left and right SLF divisions (not combined SLF)
SLF_CI_FD_plot <- SLF_FD_95CI_data[9:23,]
#add in grouping variable order
SLF_CI_FD_plot$group_order_var <- factor(SLF_CI_FD_plot$group_order_var, levels = c("C > SCD", "C > SCD", "C > aMCI", "C > aMCI", "C > aMCI", "C > mMCI", "C > mMCI", "C > mMCI", "C > mMCI", "C > AD", "C > AD", "C > AD", "C > AD", "C > AD", "SCD > AD"))
SLF_CI_FD_plot$group_order_var <- factor(SLF_CI_FD_plot$group_order_var, levels = c("C > SCD", "C > aMCI", "C > mMCI", "C > AD", "SCD > AD"))

SLF_CI_FD_plot$group_order_var <- as.factor(c("C > SCD", "C > SCD", "C > aMCI", "C > aMCI", "C > aMCI", "C > mMCI", "C > mMCI", "C > mMCI", "C > mMCI", "C > AD", "C > AD", "C > AD", "C > AD", "C > AD", "SCD > AD"))

#plot data (just with the right and left SLF divisions)
ggplot(SLF_CI_FD_plot, aes(x=SLF_group_number, y=estimate_diff, group=group_order_var, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("6" = "Right SLF 1", "4" = "Left SLF 2", "5" = "Left SLF 3", "7" = "Right SLF 2", "8" = "Right SLF 3"))+
    labs(colour="Group Contrast")+ 
    scale_color_manual(values = c("#FF00FF", "#33CC66", "#33CCFF", "#999900", "#CCCCCC"))+
    theme_classic()+
    coord_flip()
    
#plot 95% Confidence Interval for ANCOVA (age & sex)
#for FD
#Create dataset: 
SLF_FD_95CI_2covar_data <- data.frame(SLF_group_number = c('1','1','3','3','5','5','5','7','8','8'),
                               SLF_type = c('Whole_SLF','Whole_SLF', 'Right_SLF', 'Right_SLF','Left_SLF3','Left_SLF3','Left_SLF3','Right_SLF2','Right_SLF3','Right_SLF3'),
                               Group_contrast = c('CvmMCI', 'CvAD', 'CvmMCI', 'CvAD','CvSCD','CvmMCI','CvAD','CvmMCI','CvSCD','CvmMCI'),
                               estimate_diff = c(0.0149701,0.0171247,0.016138,0.017475,0.013165,0.015781,0.0174047,0.0201535,0.0142024,0.019025),
                               lower = c(0.0286742,0.0336528,0.030519,0.03482,0.0257797,0.0289482,0.0332852,0.0397022,0.02792,0.0333434), 
                               upper = c(0.0012659,0.0005967,0.001757,0.000131,0.0005503,0.0026138,0.0015243,0.0006048,0.0004848,0.0047066))  
#plot data
ggplot(SLF_FD_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Whole SLF","3" = "Right SLF","5" = "Left SLF 3", "7" = "Right SLF 2", "8" = "Right SLF 3"))+
    theme_classic() 
    


# #for ANCOVA 
# post_hoc_SLF_FD_age_mod <- glht(SLF_FD_age_mod, linfct = mcp(Group = "Tukey"))
# summary(post_hoc_SLF_FD_age_mod)
# confint(post_hoc_SLF_FD_age_mod)
# 
# post_hoc_SLF_FD_age_mod_L <- glht(SLF_FD_age_mod_L, linfct = mcp(Group = "Tukey"))
# summary(post_hoc_SLF_FD_age_mod_L)
# confint(post_hoc_SLF_FD_age_mod_L)
# 
# post_hoc_SLF_FD_age_mod_R <- glht(SLF_FD_age_mod_R, linfct = mcp(Group = "Tukey"))
# summary(post_hoc_SLF_FD_age_mod_R)
# confint(post_hoc_SLF_FD_age_mod_R)
# 
# #for ANCOVA - clinical site
# post_hoc_SLF_FD_clinsite_mod <- glht(SLF_FD_clinsite_mod, linfct = mcp(Group = "Tukey"))
# summary(post_hoc_SLF_FD_clinsite_mod)
# confint(post_hoc_SLF_FD_clinsite_mod)
# 
# post_hoc_SLF_FD_clinsite_mod_L <- glht(SLF_FD_clinsite_mod_L, linfct = mcp(Group = "Tukey"))
# summary(post_hoc_SLF_FD_clinsite_mod_L)
# confint(post_hoc_SLF_FD_clinsite_mod_L)
# 
# post_hoc_SLF_FD_clinsite_mod_R <- glht(SLF_FD_clinsite_mod_R, linfct = mcp(Group = "Tukey"))
# summary(post_hoc_SLF_FD_clinsite_mod_R)
# confint(post_hoc_SLF_FD_clinsite_mod_R)
# 
# 
# #for ANCOVA - all 3 covariates (age, sex, clinical site)
# post_hoc_SLF_FD_3covar_mod <- glht(SLF_FD_3covar_mod, linfct = mcp(Group = "Tukey"))
# summary(post_hoc_SLF_FD_3covar_mod)
# confint(post_hoc_SLF_FD_3covar_mod)
# 
# post_hoc_SLF_FD_3covar_mod_L <- glht(SLF_FD_3covar_mod_L, linfct = mcp(Group = "Tukey"))
# summary(post_hoc_SLF_FD_3covar_mod_L)
# confint(post_hoc_SLF_FD_3covar_mod_L)
# 
# post_hoc_SLF_FD_3covar_mod_R <- glht(SLF_FD_3covar_mod_R, linfct = mcp(Group = "Tukey"))
# summary(post_hoc_SLF_FD_3covar_mod_R)
# confint(post_hoc_SLF_FD_3covar_mod_R)

#calculate the effect size (eta-squared)
etaSquared(SLF_FD_mod)
etaSquared(SLF_FD_mod_L)
etaSquared(SLF_FD_mod_R)
etaSquared(SLF1_FD_mod_L)
etaSquared(SLF2_FD_mod_L)
etaSquared(SLF3_FD_mod_L)
etaSquared(SLF1_FD_mod_R)
etaSquared(SLF2_FD_mod_R)
etaSquared(SLF3_FD_mod_R)
#for ANCOVA
#calculate the effect size (eta-squared)
etaSquared(SLF_FD_2covar_mod)
etaSquared(SLF_FD_2covar_mod_L)
etaSquared(SLF_FD_2covar_mod_R)
etaSquared(SLF1_FD_2covar_mod_L)
etaSquared(SLF2_FD_2covar_mod_L)
etaSquared(SLF3_FD_2covar_mod_L)
etaSquared(SLF1_FD_2covar_mod_R)
etaSquared(SLF2_FD_2covar_mod_R)
etaSquared(SLF3_FD_2covar_mod_R)

#effect size for sig. post hoc tests (Cohen's d)
#Whole SLF FD
    #for C vs. aMCI 
    DPRC_neuropsych_data_CvaMCI_wholeSLF <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 3)
    DPRC_neuropsych_data_CvaMCI_wholeSLF$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_wholeSLF$Group)
    cohensD(mn_FD_SLF ~ Group, data = DPRC_neuropsych_data_CvaMCI_wholeSLF) #this looks like Hedges' g? 
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_wholeSLF <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_wholeSLF$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_wholeSLF$Group)
    cohensD(mn_FD_SLF ~ Group, data = DPRC_neuropsych_data_CvmMCI_wholeSLF) #this looks like Hedges' g? 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_wholeSLF <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_wholeSLF$Group <- droplevels(DPRC_neuropsych_data_CvAD_wholeSLF$Group)
    cohensD(mn_FD_SLF ~ Group, data = DPRC_neuropsych_data_CvAD_wholeSLF) #this looks like Hedges' g? 
#Left SLF FD
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF_L$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF_L$Group)
    cohensD(mn_FD_SLF_L ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF_L) #this looks like Hedges' g? 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF_L$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF_L$Group)
    cohensD(mn_FD_SLF_L ~ Group, data = DPRC_neuropsych_data_CvAD_SLF_L) #this looks like Hedges' g? 
#Right SLF FD
    #for C vs. aMCI 
    DPRC_neuropsych_data_CvaMCI_SLF_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 3)
    DPRC_neuropsych_data_CvaMCI_SLF_R$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_SLF_R$Group)
    cohensD(mn_FD_SLF_R ~ Group, data = DPRC_neuropsych_data_CvaMCI_SLF_R) #this looks like Hedges' g? 
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF_R$Group)
    cohensD(mn_FD_SLF_R ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF_R) #this looks like Hedges' g? 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF_R$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF_R$Group)
    cohensD(mn_FD_SLF_R ~ Group, data = DPRC_neuropsych_data_CvAD_SLF_R) #this looks like Hedges' g? 
#Left SLF2 FD
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF2_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF2_L$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF2_L$Group)
    cohensD(mn_FD_SLF2_L ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF2_L) #this looks like Hedges' g? 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF2_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF2_L$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF2_L$Group)
    cohensD(mn_FD_SLF2_L ~ Group, data = DPRC_neuropsych_data_CvAD_SLF2_L) #this looks like Hedges' g? 
    #for SCD vs. AD
    DPRC_neuropsych_data_SCDvAD_SLF2_L <- subset(SLF_data, SLF_data$Group == 2 | SLF_data$Group == 5)
    DPRC_neuropsych_data_SCDvAD_SLF2_L$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_SLF2_L$Group)
    cohensD(mn_FD_SLF2_L ~ Group, data = DPRC_neuropsych_data_SCDvAD_SLF2_L) #this looks like Hedges' g? 
#Left SLF3 FD
    #for C vs. SCD 
    DPRC_neuropsych_data_CvSCD_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
    DPRC_neuropsych_data_CvSCD_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_L$Group)
    cohensD(mn_FD_SLF3_L ~ Group, data = DPRC_neuropsych_data_CvSCD_SLF3_L) #this looks like Hedges' g? 
    #for C vs. aMCI 
    DPRC_neuropsych_data_CvaMCI_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 3)
    DPRC_neuropsych_data_CvaMCI_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_SLF3_L$Group)
    cohensD(mn_FD_SLF3_L ~ Group, data = DPRC_neuropsych_data_CvaMCI_SLF3_L) #this looks like Hedges' g? 
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_L$Group)
    cohensD(mn_FD_SLF3_L ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF3_L) #this looks like Hedges' g? 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF3_L$Group)
    cohensD(mn_FD_SLF3_L ~ Group, data = DPRC_neuropsych_data_CvAD_SLF3_L) #this looks like Hedges' g? 
#Right SLF1 FD
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF1_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF1_R$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF1_R$Group)
    cohensD(mn_FD_SLF1_R ~ Group, data = DPRC_neuropsych_data_CvAD_SLF1_R) #this looks like Hedges' g? 
#Right SLF2 FD
    #for C vs. aMCI 
    DPRC_neuropsych_data_CvaMCI_SLF2_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 3)
    DPRC_neuropsych_data_CvaMCI_SLF2_R$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_SLF2_R$Group)
    cohensD(mn_FD_SLF2_R ~ Group, data = DPRC_neuropsych_data_CvaMCI_SLF2_R) #this looks like Hedges' g? 
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF2_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF2_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF2_R$Group)
    cohensD(mn_FD_SLF2_R ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF2_R) #this looks like Hedges' g? 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF2_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF2_R$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF2_R$Group)
    cohensD(mn_FD_SLF2_R ~ Group, data = DPRC_neuropsych_data_CvAD_SLF2_R) #this looks like Hedges' g? 
#Right SLF3 FD
    #for C vs. SCD 
    DPRC_neuropsych_data_CvSCD_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
    DPRC_neuropsych_data_CvSCD_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_R$Group)
    cohensD(mn_FD_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvSCD_SLF3_R) #this looks like Hedges' g? 
    #for C vs. aMCI 
    DPRC_neuropsych_data_CvaMCI_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 3)
    DPRC_neuropsych_data_CvaMCI_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_SLF3_R$Group)
    cohensD(mn_FD_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvaMCI_SLF3_R) #this looks like Hedges' g? 
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_R$Group)
    cohensD(mn_FD_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF3_R) #this looks like Hedges' g? 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF3_R$Group)
    cohensD(mn_FD_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvAD_SLF3_R) #this looks like Hedges' g? 
    
#For ANCOVA - effect size for sig. post hoc tests (Cohen's d) (using a.tes function from the compute.es package)
#Whole SLF FD
t_value_effect_size <- summary(post_hoc_SLF_FD_ancova_mod) 
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_wholeSLF <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_wholeSLF$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_wholeSLF$Group)
    #cohensD(mn_FD_SLF ~ Group + as.factor(Age) + Sex, data = DPRC_neuropsych_data_CvmMCI_wholeSLF) #nope, don't work like this w/ covariates...so see the solution below 
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_wholeSLF, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_wholeSLF)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_wholeSLF <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_wholeSLF$Group <- droplevels(DPRC_neuropsych_data_CvAD_wholeSLF$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_wholeSLF, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_wholeSLF)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#Left SLF FD - no sig f/us
#Right SLF FD
t_value_effect_size <- summary(post_hoc_SLF_R_FD_ancova_mod) 
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF_R, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF_R ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_SLF_R)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF_R$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_SLF_R, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF_R ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_SLF_R)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#Left SLF2 FD - no sig. f/us
#Left SLF3 FD
t_value_effect_size <- summary(post_hoc_SLF3_L_FD_ancova_mod) 
    #for C vs. SCD
    DPRC_neuropsych_data_CvSCD_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
    DPRC_neuropsych_data_CvSCD_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_L$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF3_L, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age + Sex, data = DPRC_neuropsych_data_CvSCD_SLF3_L)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_L$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF3_L, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_SLF3_L)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF3_L$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_SLF3_L, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_SLF3_L)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#Right SLF2 FD
t_value_effect_size <- summary(post_hoc_SLF2_R_FD_ancova_mod) 
   #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF2_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF2_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF2_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF2_R, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF2_R ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_SLF2_R)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#Right SLF3 FD
t_value_effect_size <- summary(post_hoc_SLF3_R_FD_ancova_mod) 
    #for C vs. SCD 
    DPRC_neuropsych_data_CvSCD_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
    DPRC_neuropsych_data_CvSCD_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_SLF3_R, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF3_R ~ Age + Sex, data = DPRC_neuropsych_data_CvSCD_SLF3_R)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_R$Group)
    cohensD(mn_FD_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF3_R) #this looks like Hedges' g? 
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF3_R, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF3_R ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_SLF3_R)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
    
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
SLF1_FC_mod_L <- lm(mn_FC_SLF1_L ~ Group, data = SLF_data)
SLF2_FC_mod_L <- lm(mn_FC_SLF2_L ~ Group, data = SLF_data)
SLF3_FC_mod_L <- lm(mn_FC_SLF3_L ~ Group, data = SLF_data)
SLF1_FC_mod_R <- lm(mn_FC_SLF1_R ~ Group, data = SLF_data)
SLF2_FC_mod_R <- lm(mn_FC_SLF2_R ~ Group, data = SLF_data)
SLF3_FC_mod_R <- lm(mn_FC_SLF3_R ~ Group, data = SLF_data)

#include the covariate of clinical site (run an ANCOVA) in model
SLF_FC_clinsite_mod <- lm(mn_FC_SLF ~ Group + ClinSite_covar, data = SLF_data)
SLF_FC_clinsite_mod_L <- lm(mn_FC_SLF_L ~ Group + ClinSite_covar, data = SLF_data)
SLF_FC_clinsite_mod_R <- lm(mn_FC_SLF_R ~ Group + ClinSite_covar, data = SLF_data)

#include the covariate of 2 variables (age + sex) (run an ANCOVA) in model
SLF_FC_2covar_mod <- lm(mn_FC_SLF ~ Group + Age + Sex, data = SLF_data)
SLF_FC_2covar_mod_L <- lm(mn_FC_SLF_L ~ Group + Age + Sex, data = SLF_data)
SLF_FC_2covar_mod_R <- lm(mn_FC_SLF_R ~ Group + Age + Sex, data = SLF_data)
SLF1_FC_2covar_mod_L <- lm(mn_FC_SLF1_L ~ Group + Age + Sex, data = SLF_data)
SLF2_FC_2covar_mod_L <- lm(mn_FC_SLF2_L ~ Group + Age + Sex, data = SLF_data)
SLF3_FC_2covar_mod_L <- lm(mn_FC_SLF3_L ~ Group + Age + Sex, data = SLF_data)
SLF1_FC_2covar_mod_R <- lm(mn_FC_SLF1_R ~ Group + Age + Sex, data = SLF_data)
SLF2_FC_2covar_mod_R <- lm(mn_FC_SLF2_R ~ Group + Age + Sex, data = SLF_data)
SLF3_FC_2covar_mod_R <- lm(mn_FC_SLF3_R ~ Group + Age + Sex, data = SLF_data)

#median
#SLF_FC_mod <- lm(md_FC ~ Group, data = SLF_data)

anova(SLF_FC_mod)
anova(SLF_FC_mod_L)
anova(SLF_FC_mod_R)
anova(SLF1_FC_mod_L)
anova(SLF2_FC_mod_L)
anova(SLF3_FC_mod_L)
anova(SLF1_FC_mod_R)
anova(SLF2_FC_mod_R)
anova(SLF3_FC_mod_R)
#run Bayesian ANOVA
anovaBF(mn_FC_SLF ~ Group, data = SLF_data) 
anovaBF(mn_FC_SLF_L ~ Group, data = SLF_data) 
anovaBF(mn_FC_SLF_R ~ Group, data = SLF_data) 
#for ANCOVA
anova(SLF_FC_clinsite_mod)
anova(SLF_FC_clinsite_mod_L)
anova(SLF_FC_clinsite_mod_R)
#ANCOVA - for age and sex
anova(SLF_FC_2covar_mod)
anova(SLF_FC_2covar_mod_L)
anova(SLF_FC_2covar_mod_R)
anova(SLF1_FC_2covar_mod_L)
anova(SLF2_FC_2covar_mod_L)
anova(SLF3_FC_2covar_mod_L)
anova(SLF1_FC_2covar_mod_R)
anova(SLF2_FC_2covar_mod_R)
anova(SLF3_FC_2covar_mod_R)

#Linear Trend Analysis in Linear Regression
SLF_FC_LinTrend_mod <- lm(mn_FC_SLF ~ Trend_Group + Group, data = SLF_data)
anova(SLF_FC_LinTrend_mod)
summary(SLF_FC_LinTrend_mod) 
SLF_L_FC_LinTrend_mod <- lm(mn_FC_SLF_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF_L_FC_LinTrend_mod)
summary(SLF_L_FC_LinTrend_mod) 
SLF_R_FC_LinTrend_mod <- lm(mn_FC_SLF_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF_R_FC_LinTrend_mod)
summary(SLF_R_FC_LinTrend_mod) 
SLF1_L_FC_LinTrend_mod <- lm(mn_FC_SLF1_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF1_L_FC_LinTrend_mod)
summary(SLF1_L_FC_LinTrend_mod) 
SLF2_L_FC_LinTrend_mod <- lm(mn_FC_SLF2_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF2_L_FC_LinTrend_mod)
summary(SLF2_L_FC_LinTrend_mod) 
SLF3_L_FC_LinTrend_mod <- lm(mn_FC_SLF3_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF3_L_FC_LinTrend_mod)
summary(SLF3_L_FC_LinTrend_mod)
SLF1_R_FC_LinTrend_mod <- lm(mn_FC_SLF1_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF1_R_FC_LinTrend_mod)
summary(SLF1_R_FC_LinTrend_mod) 
SLF2_R_FC_LinTrend_mod <- lm(mn_FC_SLF2_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF2_R_FC_LinTrend_mod)
summary(SLF2_R_FC_LinTrend_mod) 
SLF3_R_FC_LinTrend_mod <- lm(mn_FC_SLF3_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF3_R_FC_LinTrend_mod)
summary(SLF3_R_FC_LinTrend_mod) 

#Linear Trend Analysis in Linear Regression w/ covariates age + sex
SLF_FC_LinTrend_2covar_mod <- lm(mn_FC_SLF ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF_FC_LinTrend_2covar_mod)
summary(SLF_FC_LinTrend_2covar_mod) 
SLF_L_FC_LinTrend_2covar_mod <- lm(mn_FC_SLF_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF_L_FC_LinTrend_2covar_mod)
summary(SLF_L_FC_LinTrend_2covar_mod) 
SLF_R_FC_LinTrend_2covar_mod <- lm(mn_FC_SLF_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF_R_FC_LinTrend_2covar_mod)
summary(SLF_R_FC_LinTrend_2covar_mod) 
SLF1_L_FC_LinTrend_2covar_mod <- lm(mn_FC_SLF1_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF1_L_FC_LinTrend_2covar_mod)
summary(SLF1_L_FC_LinTrend_2covar_mod) 
SLF2_L_FC_LinTrend_2covar_mod <- lm(mn_FC_SLF2_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF2_L_FC_LinTrend_2covar_mod)
summary(SLF2_L_FC_LinTrend_2covar_mod) 
SLF3_L_FC_LinTrend_2covar_mod <- lm(mn_FC_SLF3_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF3_L_FC_LinTrend_2covar_mod)
summary(SLF3_L_FC_LinTrend_2covar_mod)
SLF1_R_FC_LinTrend_2covar_mod <- lm(mn_FC_SLF1_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF1_R_FC_LinTrend_2covar_mod)
summary(SLF1_R_FC_LinTrend_2covar_mod) 
SLF2_R_FC_LinTrend_2covar_mod <- lm(mn_FC_SLF2_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF2_R_FC_LinTrend_2covar_mod)
summary(SLF2_R_FC_LinTrend_2covar_mod) 
SLF3_R_FC_LinTrend_2covar_mod <- lm(mn_FC_SLF3_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF3_R_FC_LinTrend_2covar_mod)
summary(SLF3_R_FC_LinTrend_2covar_mod) 

#for FDC
#mean
SLF_FDC_mod <- lm(mn_FDC_SLF ~ Group, data = SLF_data)
SLF_FDC_mod_L <- lm(mn_FDC_SLF_L ~ Group, data = SLF_data)
SLF_FDC_mod_R <- lm(mn_FDC_SLF_R ~ Group, data = SLF_data)
SLF1_FDC_mod_L <- lm(mn_FDC_SLF1_L ~ Group, data = SLF_data)
SLF2_FDC_mod_L <- lm(mn_FDC_SLF2_L ~ Group, data = SLF_data)
SLF3_FDC_mod_L <- lm(mn_FDC_SLF3_L ~ Group, data = SLF_data)
SLF1_FDC_mod_R <- lm(mn_FDC_SLF1_R ~ Group, data = SLF_data)
SLF2_FDC_mod_R <- lm(mn_FDC_SLF2_R ~ Group, data = SLF_data)
SLF3_FDC_mod_R <- lm(mn_FDC_SLF3_R ~ Group, data = SLF_data)
#include the covariate of clinical site (run an ANCOVA) in model
SLF_FDC_clinsite_mod <- lm(mn_FDC_SLF ~ Group + ClinSite_covar, data = SLF_data)
SLF_FDC_clinsite_mod_L <- lm(mn_FDC_SLF_L ~ Group + ClinSite_covar, data = SLF_data)
SLF_FDC_clinsite_mod_R <- lm(mn_FDC_SLF_R ~ Group + ClinSite_covar, data = SLF_data)
#include the covariate of 2 variables (age + sex) (run an ANCOVA) in model
SLF_FDC_2covar_mod <- lm(mn_FDC_SLF ~ Group + Age + Sex, data = SLF_data)
SLF_FDC_2covar_mod_L <- lm(mn_FDC_SLF_L ~ Group + Age + Sex, data = SLF_data)
SLF_FDC_2covar_mod_R <- lm(mn_FDC_SLF_R ~ Group + Age + Sex, data = SLF_data)
SLF1_FDC_2covar_mod_L <- lm(mn_FDC_SLF1_L ~ Group + Age + Sex, data = SLF_data)
SLF2_FDC_2covar_mod_L <- lm(mn_FDC_SLF2_L ~ Group + Age + Sex, data = SLF_data)
SLF3_FDC_2covar_mod_L <- lm(mn_FDC_SLF3_L ~ Group + Age + Sex, data = SLF_data)
SLF1_FDC_2covar_mod_R <- lm(mn_FDC_SLF1_R ~ Group + Age + Sex, data = SLF_data)
SLF2_FDC_2covar_mod_R <- lm(mn_FDC_SLF2_R ~ Group + Age + Sex, data = SLF_data)
SLF3_FDC_2covar_mod_R <- lm(mn_FDC_SLF3_R ~ Group + Age + Sex, data = SLF_data)

#median
#SLF_FDC_mod <- lm(md_FDC ~ Group, data = SLF_data)
anova(SLF_FDC_mod)
anova(SLF_FDC_mod_L)
anova(SLF_FDC_mod_R)
anova(SLF1_FDC_mod_L)
anova(SLF2_FDC_mod_L)
anova(SLF3_FDC_mod_L)
anova(SLF1_FDC_mod_R)
anova(SLF2_FDC_mod_R)
anova(SLF3_FDC_mod_R)
#run Bayesian ANOVA
# anovaBF(mn_FDC_SLF ~ Group, data = SLF_data) 
# anovaBF(mn_FDC_SLF_L ~ Group, data = SLF_data) 
# anovaBF(mn_FDC_SLF_R ~ Group, data = SLF_data) 
#for ANCOVA
# anova(SLF_FDC_clinsite_mod)
# anova(SLF_FDC_clinsite_mod_L)
# anova(SLF_FDC_clinsite_mod_R)
#ANCOVA - for age and sex
anova(SLF_FDC_2covar_mod)
anova(SLF_FDC_2covar_mod_L)
anova(SLF_FDC_2covar_mod_R)
anova(SLF1_FDC_2covar_mod_L)
anova(SLF2_FDC_2covar_mod_L)
anova(SLF3_FDC_2covar_mod_L)
anova(SLF1_FDC_2covar_mod_R)
anova(SLF2_FDC_2covar_mod_R)
anova(SLF3_FDC_2covar_mod_R)

#Linear Trend Analysis in Linear Regression
SLF_FDC_LinTrend_mod <- lm(mn_FDC_SLF ~ Trend_Group + Group, data = SLF_data)
anova(SLF_FDC_LinTrend_mod)
summary(SLF_FDC_LinTrend_mod) 
SLF_L_FDC_LinTrend_mod <- lm(mn_FDC_SLF_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF_L_FDC_LinTrend_mod)
summary(SLF_L_FDC_LinTrend_mod) 
SLF_R_FDC_LinTrend_mod <- lm(mn_FDC_SLF_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF_R_FDC_LinTrend_mod)
summary(SLF_R_FDC_LinTrend_mod) 
SLF1_L_FDC_LinTrend_mod <- lm(mn_FDC_SLF1_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF1_L_FDC_LinTrend_mod)
summary(SLF1_L_FDC_LinTrend_mod) 
SLF2_L_FDC_LinTrend_mod <- lm(mn_FDC_SLF2_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF2_L_FDC_LinTrend_mod)
summary(SLF2_L_FDC_LinTrend_mod) 
SLF3_L_FDC_LinTrend_mod <- lm(mn_FDC_SLF3_L ~ Trend_Group + Group, data = SLF_data)
anova(SLF3_L_FDC_LinTrend_mod)
summary(SLF3_L_FDC_LinTrend_mod)
SLF1_R_FDC_LinTrend_mod <- lm(mn_FDC_SLF1_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF1_R_FDC_LinTrend_mod)
summary(SLF1_R_FDC_LinTrend_mod) 
SLF2_R_FDC_LinTrend_mod <- lm(mn_FDC_SLF2_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF2_R_FDC_LinTrend_mod)
summary(SLF2_R_FDC_LinTrend_mod) 
SLF3_R_FDC_LinTrend_mod <- lm(mn_FDC_SLF3_R ~ Trend_Group + Group, data = SLF_data)
anova(SLF3_R_FDC_LinTrend_mod)
summary(SLF3_R_FDC_LinTrend_mod) 

#Linear Trend Analysis in Linear Regression w/ covariates age + sex
SLF_FDC_LinTrend_2covar_mod <- lm(mn_FDC_SLF ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF_FDC_LinTrend_2covar_mod)
summary(SLF_FDC_LinTrend_2covar_mod) 
SLF_L_FDC_LinTrend_2covar_mod <- lm(mn_FDC_SLF_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF_L_FDC_LinTrend_2covar_mod)
summary(SLF_L_FDC_LinTrend_2covar_mod) 
SLF_R_FDC_LinTrend_2covar_mod <- lm(mn_FDC_SLF_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF_R_FDC_LinTrend_2covar_mod)
summary(SLF_R_FDC_LinTrend_2covar_mod) 
SLF1_L_FDC_LinTrend_2covar_mod <- lm(mn_FDC_SLF1_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF1_L_FDC_LinTrend_2covar_mod)
summary(SLF1_L_FDC_LinTrend_2covar_mod) 
SLF2_L_FDC_LinTrend_2covar_mod <- lm(mn_FDC_SLF2_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF2_L_FDC_LinTrend_2covar_mod)
summary(SLF2_L_FDC_LinTrend_2covar_mod) 
SLF3_L_FDC_LinTrend_2covar_mod <- lm(mn_FDC_SLF3_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF3_L_FDC_LinTrend_2covar_mod)
summary(SLF3_L_FDC_LinTrend_2covar_mod)
SLF1_R_FDC_LinTrend_2covar_mod <- lm(mn_FDC_SLF1_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF1_R_FDC_LinTrend_2covar_mod)
summary(SLF1_R_FDC_LinTrend_2covar_mod) 
SLF2_R_FDC_LinTrend_2covar_mod <- lm(mn_FDC_SLF2_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF2_R_FDC_LinTrend_2covar_mod)
summary(SLF2_R_FDC_LinTrend_2covar_mod) 
SLF3_R_FDC_LinTrend_2covar_mod <- lm(mn_FDC_SLF3_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF3_R_FDC_LinTrend_2covar_mod)
summary(SLF3_R_FDC_LinTrend_2covar_mod) 


#run pairwise comparisons, given that the F-test was significant. 
#whole SLF
post_hoc_SLF_FDC_mod <- glht(SLF_FDC_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FDC_mod)
confint(post_hoc_SLF_FDC_mod)
#Left and Right SLF whole
post_hoc_SLF_FDC_mod_L <- glht(SLF_FDC_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FDC_mod_L)
confint(post_hoc_SLF_FDC_mod_L)
post_hoc_SLF_FDC_mod_R <- glht(SLF_FDC_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FDC_mod_R)
confint(post_hoc_SLF_FDC_mod_R)
#Left SLF 1, 2, and 3
post_hoc_SLF1_FDC_mod_L <- glht(SLF1_FDC_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF1_FDC_mod_L)
confint(post_hoc_SLF1_FDC_mod_L)
post_hoc_SLF2_FDC_mod_L <- glht(SLF2_FDC_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_FDC_mod_L)
confint(post_hoc_SLF2_FDC_mod_L)
post_hoc_SLF3_FDC_mod_L <- glht(SLF3_FDC_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_FDC_mod_L)
confint(post_hoc_SLF3_FDC_mod_L)
#Right SLF 1, 2 and 3
post_hoc_SLF1_FDC_mod_R <- glht(SLF1_FDC_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF1_FDC_mod_R)
confint(post_hoc_SLF1_FDC_mod_R)
post_hoc_SLF2_FDC_mod_R <- glht(SLF2_FDC_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_FDC_mod_R)
confint(post_hoc_SLF3_FDC_mod_R)
post_hoc_SLF3_FDC_mod_R <- glht(SLF3_FDC_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_FDC_mod_R)
confint(post_hoc_SLF3_FDC_mod_R)

#post hoc for covariate (age & sex)
post_hoc_SLF2_L_FDC_ancova_mod <- glht(SLF2_FDC_2covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_L_FDC_ancova_mod) # no sig. f/u
confint(post_hoc_SLF2_L_FDC_ancova_mod)
post_hoc_SLF3_R_FDC_ancova_mod <- glht(SLF3_FDC_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_R_FDC_ancova_mod)
confint(post_hoc_SLF3_R_FDC_ancova_mod)


#plot 95% Confidence Interval
#for FDC
#Create dataset: 
SLF_FDC_95CI_data <- data.frame(SLF_group_number = c('1'),
                                SLF_type = c('Right_SLF3'),
                                Group_contrast = c('C > SCD'),
                                estimate_diff = c(0.0332702),
                                lower = c(0.0636956), 
                                upper = c(0.0028449)) 
#plot data
ggplot(SLF_FDC_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Right SLF 3"))+
    theme_classic() 

#plot data (just with the right and left SLF divisions)
ggplot(SLF_FDC_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Right SLF 3"))+
    labs(colour="Group Contrast")+ 
    scale_color_manual(values = "#999900")+
    theme_classic()+
    coord_flip()


#plot 95% Confidence Interval - with covariates (age + sex)
#for FDC
#Create dataset: 
SLF_FDC_95CI_2covar_data <- data.frame(SLF_group_number = c('1','1'),
                                SLF_type = c('Right_SLF3','Right_SLF3'),
                                Group_contrast = c('CvSCD', 'CvaMCI'),
                                estimate_diff = c(0.0345616,0.0317537),
                                lower = c(0.0647873,0.0630815), 
                                upper = c(0.0043359,0.0004259)) 
#plot data
ggplot(SLF_FDC_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Right SLF 3"))+
    theme_classic() 

#calculate the effect size (eta-squared)
etaSquared(SLF_FDC_mod)
etaSquared(SLF_FDC_mod_L)
etaSquared(SLF_FDC_mod_R)
etaSquared(SLF1_FDC_mod_L)
etaSquared(SLF2_FDC_mod_L)
etaSquared(SLF3_FDC_mod_L)
etaSquared(SLF1_FDC_mod_R)
etaSquared(SLF2_FDC_mod_R)
etaSquared(SLF3_FDC_mod_R)
#for ANCOVA
#calculate the effect size (eta-squared)
etaSquared(SLF_FDC_2covar_mod)
etaSquared(SLF_FDC_2covar_mod_L)
etaSquared(SLF_FDC_2covar_mod_R)
etaSquared(SLF1_FDC_2covar_mod_L)
etaSquared(SLF2_FDC_2covar_mod_L)
etaSquared(SLF3_FDC_2covar_mod_L)
etaSquared(SLF1_FDC_2covar_mod_R)
etaSquared(SLF2_FDC_2covar_mod_R)
etaSquared(SLF3_FDC_2covar_mod_R)

#effect size for sig. post hoc tests (Cohen's d)
#Left SLF2 FDC 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF2_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF2_L$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF2_L$Group)
    cohensD(mn_FDC_SLF2_L ~ Group, data = DPRC_neuropsych_data_CvAD_SLF2_L) #this looks like Hedges' g? 
#Right SLF3 FDC
    #for C vs. SCD
    DPRC_neuropsych_data_CvSCD_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
    DPRC_neuropsych_data_CvSCD_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_R$Group)
    cohensD(mn_FDC_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvSCD_SLF3_R) #this looks like Hedges' g? 
    #for C vs. aMCI (sig. for 2 covar)
    DPRC_neuropsych_data_CvaMCI_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 3)
    DPRC_neuropsych_data_CvaMCI_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_SLF3_R$Group)
    cohensD(mn_FDC_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvaMCI_SLF3_R) #this looks like Hedges' g? 
    #for C vs. mMCI (sig. for 2 covar)
    DPRC_neuropsych_data_CvmMCI_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_R$Group)
    cohensD(mn_FDC_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF3_R) #this looks like Hedges' g? 


    
#effect size for sig. post hoc tests (with sex + age as covariates)
#Left SLF2 FDC - no sig. f/us
#Right SLF3 FD
    t_value_effect_size <- summary(post_hoc_SLF3_R_FDC_ancova_mod) 
    #for C vs. SCD 
    DPRC_neuropsych_data_CvSCD_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
    DPRC_neuropsych_data_CvSCD_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_SLF3_R, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FDC_SLF3_R ~ Age + Sex, data = DPRC_neuropsych_data_CvSCD_SLF3_R)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_R$Group)
    cohensD(mn_FD_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF3_R) #this looks like Hedges' g? 
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF3_R, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FDC_SLF3_R ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_SLF3_R)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
    
    

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

#Combine all sig. plots on one plot
#put data into long format
# SLF_data_FD <- dplyr::select(SLF_data, 
#                                ParticipantID,
#                                Group,
#                                mn_FD_SLF,
#                                mn_FD_SLF_L,
#                                mn_FD_SLF_R,
#                                mn_FD_SLF1_L,
#                                mn_FD_SLF2_L,
#                                mn_FD_SLF3_L,
#                                mn_FD_SLF1_R,
#                                mn_FD_SLF2_R,
#                                mn_FD_SLF3_R)
# 
# SLF_data_FD_long <- gather(SLF_data_FD, 
#                            "SLF_type",
#                            "FD_metric",
#                            mn_FD_SLF,
#                            mn_FD_SLF_L,
#                            mn_FD_SLF_R,
#                            mn_FD_SLF1_L,
#                            mn_FD_SLF2_L,
#                            mn_FD_SLF3_L,
#                            mn_FD_SLF1_R,
#                            mn_FD_SLF2_R,
#                            mn_FD_SLF3_R)
# #All tracts SLF FD (raincloud plot)
# ggplot(SLF_data_FD_long, aes(x = SLF_type, y = FD_metric, fill = Group)) + 
#     geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_point(aes(y = FD_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
#     guides(colour=FALSE)+
#     guides(fill = guide_legend(override.aes = list(shape = NA)))+
#     xlab("SLF Tract") + 
#     ylab("Fibre Density (FD)") +
#     scale_x_discrete(labels = c("mn_FD_SLF" = "Whole SLF", "mn_FD_SLF_L" = "Left SLF", "mn_FD_SLF_R" = "Right SLF", "mn_FD_SLF1_L" = "Left SLF1","mn_FD_SLF2_L" = "Left SLF2","mn_FD_SLF3_L" = "Left SLF3","mn_FD_SLF1_R" = "Right SLF1","mn_FD_SLF2_R" = "Right SLF2","mn_FD_SLF3_R" = "Right SLF3")) + 
#     scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
#     theme_classic()

#just with left and right SLF 1, 2, 3
SLF_data_FD <- dplyr::select(SLF_data, 
                             ParticipantID,
                             Group,
                             mn_FD_SLF1_L,
                             mn_FD_SLF2_L,
                             mn_FD_SLF3_L,
                             mn_FD_SLF1_R,
                             mn_FD_SLF2_R,
                             mn_FD_SLF3_R)
SLF_data_FD_long <- gather(SLF_data_FD, 
                           "SLF_type",
                           "FD_metric",
                           mn_FD_SLF1_L,
                           mn_FD_SLF2_L,
                           mn_FD_SLF3_L,
                           mn_FD_SLF1_R,
                           mn_FD_SLF2_R,
                           mn_FD_SLF3_R)
#All tracts SLF FD (raincloud plot)
ggplot(SLF_data_FD_long, aes(x = SLF_type, y = FD_metric, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = FD_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    guides(colour=FALSE)+
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF Tract") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("mn_FD_SLF1_L" = "Left SLF1","mn_FD_SLF2_L" = "Left SLF2","mn_FD_SLF3_L" = "Left SLF3","mn_FD_SLF1_R" = "Right SLF1","mn_FD_SLF2_R" = "Right SLF2","mn_FD_SLF3_R" = "Right SLF3")) + 
    scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
    theme_classic()

#put data into long format
SLF_data_FDC <- dplyr::select(SLF_data, 
                             ParticipantID,
                             Group,
                             mn_FDC_SLF,
                             mn_FDC_SLF_L,
                             mn_FDC_SLF_R,
                             mn_FDC_SLF1_L,
                             mn_FDC_SLF2_L,
                             mn_FDC_SLF3_L,
                             mn_FDC_SLF1_R,
                             mn_FDC_SLF2_R,
                             mn_FDC_SLF3_R)

SLF_data_FDC_long <- gather(SLF_data_FDC, 
                           "SLF_type",
                           "FDC_metric",
                           mn_FDC_SLF,
                           mn_FDC_SLF_L,
                           mn_FDC_SLF_R,
                           mn_FDC_SLF1_L,
                           mn_FDC_SLF2_L,
                           mn_FDC_SLF3_L,
                           mn_FDC_SLF1_R,
                           mn_FDC_SLF2_R,
                           mn_FDC_SLF3_R)
#All tracts SLF FDC (raincloud plot)
ggplot(SLF_data_FDC_long, aes(x = SLF_type, y = FDC_metric, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = FDC_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    guides(colour=FALSE)+
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF Tract") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("mn_FDC_SLF" = "Whole SLF", "mn_FDC_SLF_L" = "Left SLF", "mn_FDC_SLF_R" = "Right SLF", "mn_FDC_SLF1_L" = "Left SLF1","mn_FDC_SLF2_L" = "Left SLF2","mn_FDC_SLF3_L" = "Left SLF3","mn_FDC_SLF1_R" = "Right SLF1","mn_FDC_SLF2_R" = "Right SLF2","mn_FDC_SLF3_R" = "Right SLF3")) + 
    scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
    theme_classic()

#just with left and right SLF 1, 2, 3
SLF_data_FDC <- dplyr::select(SLF_data, 
                              ParticipantID,
                              Group,
                              mn_FDC_SLF1_L,
                              mn_FDC_SLF2_L,
                              mn_FDC_SLF3_L,
                              mn_FDC_SLF1_R,
                              mn_FDC_SLF2_R,
                              mn_FDC_SLF3_R)

SLF_data_FDC_long <- gather(SLF_data_FDC, 
                            "SLF_type",
                            "FDC_metric",
                            mn_FDC_SLF1_L,
                            mn_FDC_SLF2_L,
                            mn_FDC_SLF3_L,
                            mn_FDC_SLF1_R,
                            mn_FDC_SLF2_R,
                            mn_FDC_SLF3_R)
#All tracts SLF FDC (raincloud plot)
ggplot(SLF_data_FDC_long, aes(x = SLF_type, y = FDC_metric, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = FDC_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    guides(colour=FALSE)+
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF Tract") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("mn_FDC_SLF1_L" = "Left SLF1","mn_FDC_SLF2_L" = "Left SLF2","mn_FDC_SLF3_L" = "Left SLF3","mn_FDC_SLF1_R" = "Right SLF1","mn_FDC_SLF2_R" = "Right SLF2","mn_FDC_SLF3_R" = "Right SLF3")) + 
    scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
    theme_classic()


###--------------- Correlation tests(SLF vs. neuropsych) -------------------####
#Look at correlation between SLF metrics and executive function neuropsych tests
#for FD:
#whole SLF-
cor.test(SLF_data$mn_FD_SLF, SLF_data$TrailsA.Raw)
#by group
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#scatterplot 
ggplot(SLF_data, aes(x = TrailsA.Raw, y = mn_FD_SLF))+
    geom_point(size = .5) +
    geom_smooth(method=lm)+
    xlab("Trail Making Test A") + 
    ylab("Whole SLF Tract Fibre Density (FD)") +
    theme_classic()
#scatterplot by group 
ggplot(SLF_data, aes(x = TrailsA.Raw, y = mn_FD_SLF, colour=Group))+
    geom_point(size = .5) +
    geom_smooth(method=lm, se=FALSE)+
    xlab("Trail Making Test A") + 
    ylab("Whole SLF Tract Fibre Density (FD)") +
    theme_classic()

#Whole SLF-
cor.test(SLF_data$mn_FD_SLF, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FD_SLF, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FD_SLF, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FD_SLF, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FD_SLF, SLF_data$CatFluency.Raw)
cor.test(SLF_data$mn_FD_SLF, SLF_data$Switching.Raw)
cor.test(SLF_data$mn_FD_SLF, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FD_SLF, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF-
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$TrailsA.Raw)
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$CatFluency.Raw)
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$Switching.Raw)
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$HayBTime1.Raw) 
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF_L, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF-
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$TrailsA.Raw)
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$CatFluency.Raw)
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF1-
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF1_L, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF2-
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$TrailsA.Raw)
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$CatFluency.Raw)
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$Switching.Raw)
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$HayBTime1.Raw) 
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF2_L, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF3-
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$TrailsA.Raw)
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$CatFluency.Raw)
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$Switching.Raw)
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF3_L, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF1-
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$TrailsA.Raw)
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$HayBCatA.Raw) 
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF2-
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$TrailsA.Raw)
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$CatFluency.Raw)
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$Switching.Raw)
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF2_R, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF3-
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$TrailsA.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$HayBTime2.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$HayBCatB.Raw) #not sig.

#for FC
cor.test(SLF_data$mn_FC_SLF, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF-
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF_L, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF-
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF_R, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF1-
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF1_L, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF2-
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF2_L, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF3-
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF3_L, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF1-
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF2-
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF2_R, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF3-
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$TrailsA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$HayBTime2.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FC_SLF3_R, SLF_data$HayBCatB.Raw) #not sig.

#for FDC
cor.test(SLF_data$mn_FDC_SLF, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FDC_SLF, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FDC_SLF, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FDC_SLF, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FDC_SLF, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF-
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$Switching.Raw)
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$HayBTime1.Raw)#not sig. 
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF-
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF_R, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF1-
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF1_L, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF2-
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$CatFluency.Raw)
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$Switching.Raw)
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$HayBCatB.Raw) #not sig.
#Left SLF3-
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF1-
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF2-
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$TrailsB.Raw)
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$ColorNaming.Raw)
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$HayBTime2.Raw)
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF3-
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$TrailsA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$Inhibition.Raw)
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$HayBTime2.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_R, SLF_data$HayBCatB.Raw) #not sig.

#visualise a correlation matrix plot
#create dataset with only the relevant variables (neuropsych variables and SLF metrics)
SLF_data_vars <- SLF_data[,c("TrailsA.Raw","TrailsB.Raw","ColorNaming.Raw","WordReading.Raw",
                             "Inhibition.Raw","LetFluency.Raw","CatFluency.Raw","Switching.Raw",
                             "HayBTime1.Raw","HayBTime2.Raw","HayBCatA.Raw","HayBCatB.Raw",
                             "mn_FD_SLF","mn_FC_SLF","mn_FDC_SLF","mn_FD_SLF_L","mn_FC_SLF_L",
                             "mn_FDC_SLF_L","mn_FD_SLF_R","mn_FC_SLF_R","mn_FDC_SLF_R","mn_FD_SLF1_L",
                             "mn_FC_SLF1_L","mn_FDC_SLF1_L","mn_FD_SLF2_L","mn_FC_SLF2_L",
                             "mn_FDC_SLF2_L","mn_FD_SLF3_L","mn_FC_SLF3_L","mn_FDC_SLF3_L","mn_FD_SLF1_R",
                             "mn_FC_SLF1_R","mn_FDC_SLF1_R","mn_FD_SLF2_R","mn_FC_SLF2_R",
                             "mn_FDC_SLF2_R","mn_FD_SLF3_R","mn_FC_SLF3_R","mn_FDC_SLF3_R")]
#need to drop NAs (n = 199)
SLF_data_vars_noNAS <- na.omit(SLF_data_vars)
#create a correlation matrix
cor_matrix <- cor(SLF_data_vars_noNAS)
#create a heatmap of the correlation matrix
#whole heatmap correlation plot
corrplot(cor_matrix, method = "color", tl.col = "black")
#plot just the interested variables - need to adjust dataset
cor_matrix_vars_interest <- cor_matrix[13:39,1:12]
#change correlations to all be in the same direction (for the Fluency tasks)
cor_matrix_vars_interest[,6:8] <- (cor_matrix_vars_interest[,6:8])*-1
#reorder and rename some columns
cor_matrix_vars_interest <- cor_matrix_vars_interest[,c(1,3,4,9,11,12,10,8,5,2,6,7)]
colnames(cor_matrix_vars_interest) <- c("TMT-A","Colour Naming","Word Reading","HayTime1","HayCatAError","HayCatBError","HayTime2","Category Switching","Inhibition","TMT-B","Letter Fluency","Category Fluency")
rownames(cor_matrix_vars_interest) <- c("Whole SLF FD","Whole SLF FC","Whole SLF FDC","Left SLF FD","Left SLF FC","Left SLF FDC","Right SLF FD","Right SLF FC","Right SLF FDC","Left SLF1 FD","Left SLF1 FC","Left SLF1 FDC","Left SLF2 FD","Left SLF2 FC","Left SLF2 FDC","Left SLF3 FD","Left SLF3 FC","Left SLF3 FDC","Right SLF1 FD","Right SLF1 FC","Right SLF1 FDC","Right SLF2 FD","Right SLF2 FC","Right SLF2 FDC","Right SLF3 FD","Right SLF3 FC","Right SLF3 FDC")
#plot your correlation plot!
corrplot(cor_matrix_vars_interest, method = "ellipse", tl.col = "black", is.corr=FALSE)

#visualise with p-values:
testRes<- cor.mtest(SLF_data_vars_noNAS, conf.level=0.95)
#whole heatmap plot
corrplot(cor_matrix, p.mat = testRes$p, sig.level = 0.05)
#plot just the interested variables - need to adjust dataset
testRes_pvalues_vars_interest <- testRes$p[13:39,1:12]
#reorder and rename some columns for the p-value matrix
testRes_pvalues_vars_interest <- testRes_pvalues_vars_interest[,c(1,3,4,9,11,12,10,8,5,2,6,7)]
colnames(testRes_pvalues_vars_interest) <- c("TMT-A","Colour Naming","Word Reading","HayTime1","HayCatAError","HayCatBError","HayTime2","Category Switching","Inhibition","TMT-B","Letter Fluency","Category Fluency")
rownames(testRes_pvalues_vars_interest) <- c("Whole SLF FD","Whole SLF FC","Whole SLF FDC","Left SLF FD","Left SLF FC","Left SLF FDC","Right SLF FD","Right SLF FC","Right SLF FDC","Left SLF1 FD","Left SLF1 FC","Left SLF1 FDC","Left SLF2 FD","Left SLF2 FC","Left SLF2 FDC","Left SLF3 FD","Left SLF3 FC","Left SLF3 FDC","Right SLF1 FD","Right SLF1 FC","Right SLF1 FDC","Right SLF2 FD","Right SLF2 FC","Right SLF2 FDC","Right SLF3 FD","Right SLF3 FC","Right SLF3 FDC")
#plot correlation plot with p-values!
corrplot(cor_matrix_vars_interest, p.mat = testRes_pvalues_vars_interest, sig.level = 0.05, method = "ellipse", tl.col = "black", is.corr=FALSE)


#partial correlation version: 
pcor_matrix <- pcor(SLF_data_vars_noNAS)
#whole heatmap correlation plot
corrplot(pcor_matrix$estimate, method = "color", tl.col = "black")
#plot just the interested variables - need to adjust dataset
pcor_matrix_vars_interest <- pcor_matrix$estimate[13:39,1:12]
#change correlations to all be in the same direction (for the Fluency tasks)
pcor_matrix_vars_interest[,6:8] <- (pcor_matrix_vars_interest[,6:8])*-1
#reorder and rename some columns
pcor_matrix_vars_interest <- pcor_matrix_vars_interest[,c(1,3,4,9,11,12,10,8,5,2,6,7)]
colnames(pcor_matrix_vars_interest) <- c("TMT-A","Colour Naming","Word Reading","HayTime1","HayCatAError","HayCatBError","HayTime2","Category Switching","Inhibition","TMT-B","Letter Fluency","Category Fluency")
rownames(pcor_matrix_vars_interest) <- c("Whole SLF FD","Whole SLF FC","Whole SLF FDC","Left SLF FD","Left SLF FC","Left SLF FDC","Right SLF FD","Right SLF FC","Right SLF FDC","Left SLF1 FD","Left SLF1 FC","Left SLF1 FDC","Left SLF2 FD","Left SLF2 FC","Left SLF2 FDC","Left SLF3 FD","Left SLF3 FC","Left SLF3 FDC","Right SLF1 FD","Right SLF1 FC","Right SLF1 FDC","Right SLF2 FD","Right SLF2 FC","Right SLF2 FDC","Right SLF3 FD","Right SLF3 FC","Right SLF3 FDC")
#plot your correlation plot!
corrplot(pcor_matrix_vars_interest, method = "ellipse", tl.col = "black", is.corr=FALSE)

#visualise with p-values:
#testRes<- cor.mtest(SLF_data_vars_noNAS, conf.level=0.95)
#whole heatmap plot
corrplot(pcor_matrix$estimate, p.mat = pcor_matrix$p.value, sig.level = 0.05)
#plot just the interested variables - need to adjust dataset
testRes_pcor_pvalues_vars_interest <- pcor_matrix$p.value[13:39,1:12]
#reorder and rename some columns for the p-value matrix
testRes_pcor_pvalues_vars_interest <- testRes_pcor_pvalues_vars_interest[,c(1,3,4,9,11,12,10,8,5,2,6,7)]
colnames(testRes_pcor_pvalues_vars_interest) <- c("TMT-A","Colour Naming","Word Reading","HayTime1","HayCatAError","HayCatBError","HayTime2","Category Switching","Inhibition","TMT-B","Letter Fluency","Category Fluency")
rownames(testRes_pcor_pvalues_vars_interest) <- c("Whole SLF FD","Whole SLF FC","Whole SLF FDC","Left SLF FD","Left SLF FC","Left SLF FDC","Right SLF FD","Right SLF FC","Right SLF FDC","Left SLF1 FD","Left SLF1 FC","Left SLF1 FDC","Left SLF2 FD","Left SLF2 FC","Left SLF2 FDC","Left SLF3 FD","Left SLF3 FC","Left SLF3 FDC","Right SLF1 FD","Right SLF1 FC","Right SLF1 FDC","Right SLF2 FD","Right SLF2 FC","Right SLF2 FDC","Right SLF3 FD","Right SLF3 FC","Right SLF3 FDC")
#plot correlation plot with p-values!
corrplot(pcor_matrix_vars_interest, p.mat = testRes_pcor_pvalues_vars_interest, sig.level = 0.05, method = "ellipse", tl.col = "black", is.corr=FALSE)


#Look at correlation for just 12 variables total (3 SLF tracts - whole, left, and right; 3 neurospsych 
#composite scores - processing speed, inhibition, generation)
#use z-scores of the neuropsych tests (to be on the same scale) to combine them 
#put proc speed variable z-scores onto a new dataset
proc_speed_zscores_data <- dplyr::select(DPRC_neuropsych_data, 
                                         ParticipantID,
                                         Group,
                                         Age,
                                         Sex,
                                         Sex_binary,
                                         TrailsA.Z, 
                                         ColorNaming.Z, 
                                         WordReading.Z, 
                                         HayBTime1.z)
#convert sex_binary to numeric class (for partial correlation)
proc_speed_zscores_data$Sex_binary<-as.numeric(proc_speed_zscores_data$Sex_binary)
#add in SLF data to this
proc_speed_zscores_data[,10:18] <- SLF_data[,c(1,8,15,22,29,36,43,50,57)] 
#remove NAs
proc_speed_zscores_data_noNAS <- na.omit(proc_speed_zscores_data)
#take the average z-score per every participant across the processing speed variables
proc_speed_average_zscore<-rowMeans(proc_speed_zscores_data_noNAS[,6:9])
#add this average proc speed zscore to the noNas dataframe
proc_speed_zscores_data_noNAS['proc_speed_average_zscore'] <- proc_speed_average_zscore
#measure the correlation for processing speed
#Whole SLF vs processing speed -
cor.test(proc_speed_zscores_data_noNAS$mn_FD_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
#Left SLF vs processing speed
cor.test(proc_speed_zscores_data_noNAS$mn_FD_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
#Right SLF vs processing speed
cor.test(proc_speed_zscores_data_noNAS$mn_FD_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)

#put inhibition variable z-scores onto a new dataset
inhibition_zscores_data <- dplyr::select(DPRC_neuropsych_data, 
                                         ParticipantID,
                                         Group,
                                         Age,
                                         Sex,
                                         Sex_binary,
                                         TrailsB.Z, 
                                         Inhibition.Z, 
                                         Switching.z,
                                         HayBTime2.z)
#convert sex_binary to numeric class (for partial correlation)
inhibition_zscores_data$Sex_binary<-as.numeric(inhibition_zscores_data$Sex_binary)
#add in SLF data to this
inhibition_zscores_data[,10:18] <- SLF_data[,c(1,8,15,22,29,36,43,50,57)] 
#remove NAs
inhibition_zscores_data_noNAS <- na.omit(inhibition_zscores_data)
#take the average z-score per every participant across the processing speed variables
inhibition_average_zscore<-rowMeans(inhibition_zscores_data_noNAS[,6:9])
#add this average proc speed zscore to the noNas dataframe
inhibition_zscores_data_noNAS['inhibition_average_zscore'] <- inhibition_average_zscore
#measure the correlation for inhibition
#Whole SLF vs inhibition -
cor.test(inhibition_zscores_data_noNAS$mn_FD_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore)
#Left SLF vs inhibition
cor.test(inhibition_zscores_data_noNAS$mn_FD_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore)
#Right SLF vs inhibition
cor.test(inhibition_zscores_data_noNAS$mn_FD_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FC_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore)

#put generation variable z-scores onto a new dataset
generation_zscores_data <- dplyr::select(DPRC_neuropsych_data, 
                                         ParticipantID,
                                         Group,
                                         Age,
                                         Sex,
                                         Sex_binary,
                                         LetFluency.Z,
                                         CatFluency.Z)
#convert sex_binary to numeric class (for partial correlation)
generation_zscores_data$Sex_binary<-as.numeric(generation_zscores_data$Sex_binary)
#add in SLF data to this
generation_zscores_data[,8:16] <- SLF_data[,c(1,8,15,22,29,36,43,50,57)] 
#remove NAs
generation_zscores_data_noNAS <- na.omit(generation_zscores_data)
#take the average z-score per every participant across the processing speed variables
generation_average_zscore<-rowMeans(generation_zscores_data_noNAS[,6:7])
#add this average proc speed zscore to the noNas dataframe
generation_zscores_data_noNAS['generation_average_zscore'] <- generation_average_zscore
#measure the correlation for generation
#Whole SLF vs generation -
cor.test(generation_zscores_data_noNAS$mn_FD_SLF, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FC_SLF, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FDC_SLF, generation_zscores_data_noNAS$generation_average_zscore)
#Left SLF vs generation
cor.test(generation_zscores_data_noNAS$mn_FD_SLF_L, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FDC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore)
#Right SLF vs generation
cor.test(generation_zscores_data_noNAS$mn_FD_SLF_R, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FC_SLF_R, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FDC_SLF_R, generation_zscores_data_noNAS$generation_average_zscore)

#examine partial correlation (age and sex as covariates) for this:
#Processing speed
#Whole SLF vs processing speed -
pcor.test(proc_speed_zscores_data_noNAS$mn_FD_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Left SLF vs processing speed -
pcor.test(proc_speed_zscores_data_noNAS$mn_FD_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Right SLF vs processing speed -
pcor.test(proc_speed_zscores_data_noNAS$mn_FD_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Inhibition
#Whole SLF vs inhibition -
pcor.test(inhibition_zscores_data_noNAS$mn_FD_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Left SLF vs inhibition -
pcor.test(inhibition_zscores_data_noNAS$mn_FD_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Right SLF vs inhibition -
pcor.test(inhibition_zscores_data_noNAS$mn_FD_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FC_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Generation
#Whole SLF vs generation -
pcor.test(generation_zscores_data_noNAS$mn_FD_SLF, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FC_SLF, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FDC_SLF, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Left SLF vs generation -
pcor.test(generation_zscores_data_noNAS$mn_FD_SLF_L, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FDC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Right SLF vs generation -
pcor.test(generation_zscores_data_noNAS$mn_FD_SLF_R, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FC_SLF_R, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FDC_SLF_R, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])

#Look at correlations per group
#for whole FD
#TMT-A
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#TMT-B
cor.test(formula = ~ mn_FD_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#ColorNaming
cor.test(formula = ~ mn_FD_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FD_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FD_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FD_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FD_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FD_SLF + Switching.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FD_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FD_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FD_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FD_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FD_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for Left FD
#TMT-A
cor.test(formula = ~ mn_FD_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#TMT-B
cor.test(formula = ~ mn_FD_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#ColorNaming
cor.test(formula = ~ mn_FD_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FD_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FD_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FD_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FD_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FD_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FD_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FD_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FD_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FD_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FD_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FD_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for right SLF FD
#TMT-A
cor.test(formula = ~ mn_FD_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FD_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#ColorNaming
cor.test(formula = ~ mn_FD_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FD_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FD_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FD_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FD_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FD_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FD_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FD_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FD_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FD_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FD_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for left SLF1 FD
#TMT-A
cor.test(formula = ~ mn_FD_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FD_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#ColorNaming
cor.test(formula = ~ mn_FD_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FD_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FD_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#LetFluency
cor.test(formula = ~ mn_FD_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FD_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FD_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FD_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FD_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for left SLF2 FD
#TMT-A
cor.test(formula = ~ mn_FD_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#TMT-B
cor.test(formula = ~ mn_FD_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FD_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FD_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FD_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FD_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FD_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FD_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FD_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FD_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("4")) #sig corr with mMCI group
cor.test(formula = ~ mn_FD_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for left SLF3 FD
#TMT-A
cor.test(formula = ~ mn_FD_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FD_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FD_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FD_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FD_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FD_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FD_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FD_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#HayBCatA
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for right SLF1 FD
#TMT-A
cor.test(formula = ~ mn_FD_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FD_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#ColorNaming
cor.test(formula = ~ mn_FD_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FD_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FD_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("5"))
#LetFluency
cor.test(formula = ~ mn_FD_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FD_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FD_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FD_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FD_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FD_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#HayBCatB
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for right SLF2 FD
#TMT-A
cor.test(formula = ~ mn_FD_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FD_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FD_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FD_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FD_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FD_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FD_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FD_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FD_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FD_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FD_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for right SLF3 FD
#TMT-A
cor.test(formula = ~ mn_FD_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FD_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FD_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FD_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FD_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FD_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FD_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FD_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FD_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FD_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FD_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FD_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 

#for whole FC
#TMT-A
cor.test(formula = ~ mn_FC_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FC_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FC_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FC_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FC_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FC_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FC_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FC_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FC_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("3")) #sig. corr with aMCI group
cor.test(formula = ~ mn_FC_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) #sig. corr with AD group
#Switching
cor.test(formula = ~ mn_FC_SLF + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FC_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FC_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FC_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FC_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for Left FC
#TMT-A
cor.test(formula = ~ mn_FC_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FC_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("3")) #sig. corr with aMCI group
cor.test(formula = ~ mn_FC_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FC_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FC_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FC_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FC_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("3")) #sig. corr with aMCI group
cor.test(formula = ~ mn_FC_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FC_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("3")) #sig. corr with aMCI group
cor.test(formula = ~ mn_FC_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#Switching
cor.test(formula = ~ mn_FC_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FC_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FC_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FC_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FC_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for right SLF FC
#TMT-A
cor.test(formula = ~ mn_FC_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FC_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FC_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FC_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FC_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FC_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FC_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FC_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FC_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FC_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#Switching
cor.test(formula = ~ mn_FC_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FC_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FC_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FC_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FC_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for left SLF1 FC
#TMT-A
cor.test(formula = ~ mn_FC_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FC_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FC_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FC_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FC_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FC_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FC_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FC_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("4")) #sig corr with mMCI group
cor.test(formula = ~ mn_FC_SLF1_L + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FC_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for left SLF2 FC
#TMT-A
cor.test(formula = ~ mn_FC_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#TMT-B
cor.test(formula = ~ mn_FC_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FC_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FC_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FC_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FC_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FC_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FC_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FC_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FC_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FC_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FC_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FC_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FC_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("4")) #sig corr with mMCI group
cor.test(formula = ~ mn_FC_SLF2_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for left SLF3 FC
#TMT-A
cor.test(formula = ~ mn_FC_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FC_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FC_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FC_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FC_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FC_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FC_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FC_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FC_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FC_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FC_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#HayBCatA
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_L + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for right SLF1 FC
#TMT-A
cor.test(formula = ~ mn_FC_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FC_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#ColorNaming
cor.test(formula = ~ mn_FC_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FC_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FC_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + Inhibition.Raw, data = SLF_data, subset = Group == c("5"))
#LetFluency
cor.test(formula = ~ mn_FC_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FC_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FC_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FC_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("5")) #sig corr with AD group
#HayBCatB
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF1_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for right SLF2 FC
#TMT-A
cor.test(formula = ~ mn_FC_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FC_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("1")) #sig corr with C group
cor.test(formula = ~ mn_FC_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FC_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FC_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FC_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FC_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FC_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FC_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF2_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 
#for right SLF3 FC
#TMT-A
cor.test(formula = ~ mn_FC_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + TrailsA.Raw, data = SLF_data, subset = Group == c("5")) 
#TMT-B
cor.test(formula = ~ mn_FC_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + TrailsB.Raw, data = SLF_data, subset = Group == c("5")) 
#ColorNaming
cor.test(formula = ~ mn_FC_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + ColorNaming.Raw, data = SLF_data, subset = Group == c("5")) 
#WordReading
cor.test(formula = ~ mn_FC_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + WordReading.Raw, data = SLF_data, subset = Group == c("5")) 
#Inhibition
cor.test(formula = ~ mn_FC_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + Inhibition.Raw, data = SLF_data, subset = Group == c("5")) 
#LetFluency
cor.test(formula = ~ mn_FC_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + LetFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#CatFluency
cor.test(formula = ~ mn_FC_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + CatFluency.Raw, data = SLF_data, subset = Group == c("5")) 
#Switching
cor.test(formula = ~ mn_FC_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("1")) 
cor.test(formula = ~ mn_FC_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("3")) #sig corr with aMCI group
cor.test(formula = ~ mn_FC_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + Switching.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime1
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("2")) 
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime1.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBTime2
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("2")) #sig corr with SCD group
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBTime2.Raw, data = SLF_data, subset = Group == c("5")) 
#HayBCatA
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatA.Raw, data = SLF_data, subset = Group == c("5"))
#HayBCatB
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("1"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("2"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("3"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("4"))
cor.test(formula = ~ mn_FC_SLF3_R + HayBCatB.Raw, data = SLF_data, subset = Group == c("5")) 


#Run partial correlation (controlling for covariates, age and sex) - note that
#missing values are not allowed.
#for TrailsA:
SLF_data_TrailsA <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "TrailsA.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                               "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                               "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                               "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                               "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                               "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_TrailsA <- SLF_data_TrailsA %>% drop_na(TrailsA.Raw)
SLF_data_TrailsA$Sex_binary <- as.numeric(SLF_data_TrailsA$Sex_binary)
#for FD
pcor.test(SLF_data_TrailsA$mn_FD_SLF, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FD_SLF_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FD_SLF_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FD_SLF1_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FD_SLF2_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FD_SLF3_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FD_SLF1_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FD_SLF2_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FD_SLF3_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
#for FC
pcor.test(SLF_data_TrailsA$mn_FC_SLF, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FC_SLF_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FC_SLF_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FC_SLF1_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FC_SLF2_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FC_SLF3_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FC_SLF1_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FC_SLF2_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FC_SLF3_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_TrailsA$mn_FDC_SLF, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FDC_SLF_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FDC_SLF_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FDC_SLF1_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FDC_SLF2_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FDC_SLF3_L, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FDC_SLF1_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FDC_SLF2_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsA$mn_FDC_SLF3_R, SLF_data_TrailsA$TrailsA.Raw, SLF_data_TrailsA[,c("Age", "Sex_binary")])
#for TrailsB:
SLF_data_TrailsB <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "TrailsB.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                               "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                               "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                               "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                               "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                               "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_TrailsB <- SLF_data_TrailsB %>% drop_na(TrailsB.Raw)
SLF_data_TrailsB$Sex_binary <- as.numeric(SLF_data_TrailsB$Sex_binary)
#for FD
pcor.test(SLF_data_TrailsB$mn_FD_SLF, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_TrailsB$mn_FD_SLF_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_TrailsB$mn_FD_SLF_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_TrailsB$mn_FD_SLF1_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FD_SLF2_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_TrailsB$mn_FD_SLF3_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_TrailsB$mn_FD_SLF1_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_TrailsB$mn_FD_SLF2_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_TrailsB$mn_FD_SLF3_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_TrailsB$mn_FC_SLF, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FC_SLF_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FC_SLF_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FC_SLF1_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FC_SLF2_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FC_SLF3_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FC_SLF1_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FC_SLF2_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FC_SLF3_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_TrailsB$mn_FDC_SLF, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FDC_SLF_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FDC_SLF_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FDC_SLF1_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FDC_SLF2_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FDC_SLF3_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FDC_SLF1_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FDC_SLF2_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FDC_SLF3_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
#for ColorNaming:
SLF_data_ColorNaming <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "ColorNaming.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                               "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                               "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                               "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                               "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                               "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_ColorNaming <- SLF_data_ColorNaming %>% drop_na(ColorNaming.Raw)
SLF_data_ColorNaming$Sex_binary <- as.numeric(SLF_data_ColorNaming$Sex_binary)
#for FD
pcor.test(SLF_data_ColorNaming$mn_FD_SLF, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_ColorNaming$mn_FD_SLF_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_ColorNaming$mn_FD_SLF_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_ColorNaming$mn_FD_SLF1_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FD_SLF2_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FD_SLF3_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_ColorNaming$mn_FD_SLF1_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FD_SLF2_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_ColorNaming$mn_FD_SLF3_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_ColorNaming$mn_FC_SLF, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FC_SLF_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FC_SLF_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FC_SLF1_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FC_SLF2_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FC_SLF3_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FC_SLF1_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FC_SLF2_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FC_SLF3_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF1_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF2_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF3_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF1_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF2_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF3_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
#for WordReading:
SLF_data_WordReading <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "WordReading.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                   "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                   "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                   "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                   "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                   "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_WordReading <- SLF_data_WordReading %>% drop_na(WordReading.Raw)
SLF_data_WordReading$Sex_binary <- as.numeric(SLF_data_WordReading$Sex_binary)
#for FD
pcor.test(SLF_data_WordReading$mn_FD_SLF, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_WordReading$mn_FD_SLF_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_WordReading$mn_FD_SLF_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_WordReading$mn_FD_SLF1_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FD_SLF2_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FD_SLF3_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_WordReading$mn_FD_SLF1_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FD_SLF2_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_WordReading$mn_FD_SLF3_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_WordReading$mn_FC_SLF, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FC_SLF_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FC_SLF_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FC_SLF1_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FC_SLF2_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FC_SLF3_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FC_SLF1_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FC_SLF2_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FC_SLF3_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_WordReading$mn_FDC_SLF, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FDC_SLF_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FDC_SLF_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FDC_SLF1_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FDC_SLF2_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FDC_SLF3_L, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FDC_SLF1_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FDC_SLF2_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_WordReading$mn_FDC_SLF3_R, SLF_data_WordReading$WordReading.Raw, SLF_data_WordReading[,c("Age", "Sex_binary")])
#for Inhibition:
SLF_data_Inhibition <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "Inhibition.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                  "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                  "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                  "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                  "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                  "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_Inhibition <- SLF_data_Inhibition %>% drop_na(Inhibition.Raw)
SLF_data_Inhibition$Sex_binary <- as.numeric(SLF_data_Inhibition$Sex_binary)
#for FD
pcor.test(SLF_data_Inhibition$mn_FD_SLF, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Inhibition$mn_FD_SLF_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Inhibition$mn_FD_SLF_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Inhibition$mn_FD_SLF1_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FD_SLF2_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Inhibition$mn_FD_SLF3_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Inhibition$mn_FD_SLF1_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FD_SLF2_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_Inhibition$mn_FD_SLF3_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
#for FC
pcor.test(SLF_data_Inhibition$mn_FC_SLF, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FC_SLF_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FC_SLF_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FC_SLF1_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FC_SLF2_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FC_SLF3_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FC_SLF1_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FC_SLF2_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FC_SLF3_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_Inhibition$mn_FDC_SLF, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Inhibition$mn_FDC_SLF_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Inhibition$mn_FDC_SLF_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Inhibition$mn_FDC_SLF1_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FDC_SLF2_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Inhibition$mn_FDC_SLF3_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FDC_SLF1_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FDC_SLF2_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_Inhibition$mn_FDC_SLF3_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
#for LetFluency:
SLF_data_LetFluency <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "LetFluency.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                  "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                  "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                  "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                  "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                  "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_LetFluency <- SLF_data_LetFluency %>% drop_na(LetFluency.Raw)
SLF_data_LetFluency$Sex_binary <- as.numeric(SLF_data_LetFluency$Sex_binary)
#for FD
pcor.test(SLF_data_LetFluency$mn_FD_SLF, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_LetFluency$mn_FD_SLF_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_LetFluency$mn_FD_SLF_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_LetFluency$mn_FD_SLF1_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FD_SLF2_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_LetFluency$mn_FD_SLF3_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_LetFluency$mn_FD_SLF1_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FD_SLF2_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_LetFluency$mn_FD_SLF3_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_LetFluency$mn_FC_SLF, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FC_SLF_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FC_SLF_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FC_SLF1_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FC_SLF2_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FC_SLF3_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FC_SLF1_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FC_SLF2_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FC_SLF3_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_LetFluency$mn_FDC_SLF, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FDC_SLF_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FDC_SLF_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FDC_SLF1_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FDC_SLF2_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])#sig
pcor.test(SLF_data_LetFluency$mn_FDC_SLF3_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FDC_SLF1_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FDC_SLF2_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FDC_SLF3_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
#for CatFluency:
SLF_data_CatFluency <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "CatFluency.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                  "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                  "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                  "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                  "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                  "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_CatFluency <- SLF_data_CatFluency %>% drop_na(CatFluency.Raw)
SLF_data_CatFluency$Sex_binary <- as.numeric(SLF_data_CatFluency$Sex_binary)
#for FD
pcor.test(SLF_data_CatFluency$mn_FD_SLF, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_CatFluency$mn_FD_SLF_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_CatFluency$mn_FD_SLF_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_CatFluency$mn_FD_SLF1_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FD_SLF2_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_CatFluency$mn_FD_SLF3_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_CatFluency$mn_FD_SLF1_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FD_SLF2_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_CatFluency$mn_FD_SLF3_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_CatFluency$mn_FC_SLF, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FC_SLF_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FC_SLF_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FC_SLF1_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FC_SLF2_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FC_SLF3_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FC_SLF1_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FC_SLF2_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FC_SLF3_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_CatFluency$mn_FDC_SLF, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FDC_SLF_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FDC_SLF_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FDC_SLF1_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FDC_SLF2_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FDC_SLF3_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FDC_SLF1_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FDC_SLF2_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FDC_SLF3_R, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
#for Switching:
SLF_data_Switching <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "Switching.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                  "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                  "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                  "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                  "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                  "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_Switching <- SLF_data_Switching %>% drop_na(Switching.Raw)
SLF_data_Switching$Sex_binary <- as.numeric(SLF_data_Switching$Sex_binary)
#for FD
pcor.test(SLF_data_Switching$mn_FD_SLF, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_Switching$mn_FD_SLF_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_Switching$mn_FD_SLF_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_Switching$mn_FD_SLF1_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FD_SLF2_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Switching$mn_FD_SLF3_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_Switching$mn_FD_SLF1_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FD_SLF2_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_Switching$mn_FD_SLF3_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_Switching$mn_FC_SLF, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FC_SLF_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FC_SLF_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FC_SLF1_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FC_SLF2_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FC_SLF3_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FC_SLF1_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FC_SLF2_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FC_SLF3_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_Switching$mn_FDC_SLF, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FDC_SLF_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FDC_SLF_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FDC_SLF1_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FDC_SLF2_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_Switching$mn_FDC_SLF3_L, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FDC_SLF1_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FDC_SLF2_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Switching$mn_FDC_SLF3_R, SLF_data_Switching$Switching.Raw, SLF_data_Switching[,c("Age", "Sex_binary")])
#for HayBTime1:
SLF_data_HayBTime1 <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "HayBTime1.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                 "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                 "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                 "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                 "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                 "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_HayBTime1 <- SLF_data_HayBTime1 %>% drop_na(HayBTime1.Raw)
SLF_data_HayBTime1$Sex_binary <- as.numeric(SLF_data_HayBTime1$Sex_binary)
#for FD
pcor.test(SLF_data_HayBTime1$mn_FD_SLF, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime1$mn_FD_SLF_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime1$mn_FD_SLF_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime1$mn_FD_SLF1_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FD_SLF2_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime1$mn_FD_SLF3_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime1$mn_FD_SLF1_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FD_SLF2_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime1$mn_FD_SLF3_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_HayBTime1$mn_FC_SLF, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FC_SLF_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FC_SLF_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FC_SLF1_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FC_SLF2_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FC_SLF3_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FC_SLF1_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FC_SLF2_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FC_SLF3_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_HayBTime1$mn_FDC_SLF, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FDC_SLF_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FDC_SLF_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FDC_SLF1_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FDC_SLF2_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime1$mn_FDC_SLF3_L, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FDC_SLF1_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FDC_SLF2_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime1$mn_FDC_SLF3_R, SLF_data_HayBTime1$HayBTime1.Raw, SLF_data_HayBTime1[,c("Age", "Sex_binary")])
#for HayBTime2:
SLF_data_HayBTime2 <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "HayBTime2.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                 "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                 "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                 "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                 "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                 "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_HayBTime2 <- SLF_data_HayBTime2 %>% drop_na(HayBTime2.Raw)
SLF_data_HayBTime2$Sex_binary <- as.numeric(SLF_data_HayBTime2$Sex_binary)
#for FD
pcor.test(SLF_data_HayBTime2$mn_FD_SLF, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime2$mn_FD_SLF_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime2$mn_FD_SLF_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime2$mn_FD_SLF1_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FD_SLF2_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_HayBTime2$mn_FD_SLF3_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime2$mn_FD_SLF1_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FD_SLF2_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBTime2$mn_FD_SLF3_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_HayBTime2$mn_FC_SLF, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FC_SLF_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FC_SLF_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FC_SLF1_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_HayBTime2$mn_FC_SLF2_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FC_SLF3_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FC_SLF1_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])#sig
pcor.test(SLF_data_HayBTime2$mn_FC_SLF2_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FC_SLF3_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_HayBTime2$mn_FDC_SLF, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FDC_SLF_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FDC_SLF_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FDC_SLF1_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FDC_SLF2_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")]) #sig 
pcor.test(SLF_data_HayBTime2$mn_FDC_SLF3_L, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FDC_SLF1_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FDC_SLF2_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBTime2$mn_FDC_SLF3_R, SLF_data_HayBTime2$HayBTime2.Raw, SLF_data_HayBTime2[,c("Age", "Sex_binary")])
#for HayBCatA:
SLF_data_HayBCatA <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "HayBCatA.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                 "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                 "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                 "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                 "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                 "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_HayBCatA <- SLF_data_HayBCatA %>% drop_na(HayBCatA.Raw)
SLF_data_HayBCatA$Sex_binary <- as.numeric(SLF_data_HayBCatA$Sex_binary)
#for FD
pcor.test(SLF_data_HayBCatA$mn_FD_SLF, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatA$mn_FD_SLF_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatA$mn_FD_SLF_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatA$mn_FD_SLF1_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FD_SLF2_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatA$mn_FD_SLF3_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatA$mn_FD_SLF1_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FD_SLF2_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatA$mn_FD_SLF3_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_HayBCatA$mn_FC_SLF, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FC_SLF_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FC_SLF_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FC_SLF1_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FC_SLF2_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FC_SLF3_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FC_SLF1_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FC_SLF2_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FC_SLF3_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_HayBCatA$mn_FDC_SLF, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FDC_SLF_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FDC_SLF_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FDC_SLF1_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FDC_SLF2_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatA$mn_FDC_SLF3_L, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FDC_SLF1_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FDC_SLF2_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatA$mn_FDC_SLF3_R, SLF_data_HayBCatA$HayBCatA.Raw, SLF_data_HayBCatA[,c("Age", "Sex_binary")])
#for HayBCatB:
SLF_data_HayBCatB <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "HayBCatB.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_HayBCatB <- SLF_data_HayBCatB %>% drop_na(HayBCatB.Raw)
SLF_data_HayBCatB$Sex_binary <- as.numeric(SLF_data_HayBCatB$Sex_binary)
#for FD
pcor.test(SLF_data_HayBCatB$mn_FD_SLF, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatB$mn_FD_SLF_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatB$mn_FD_SLF_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatB$mn_FD_SLF1_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FD_SLF2_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatB$mn_FD_SLF3_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatB$mn_FD_SLF1_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FD_SLF2_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatB$mn_FD_SLF3_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_HayBCatB$mn_FC_SLF, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FC_SLF_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FC_SLF_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FC_SLF1_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FC_SLF2_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FC_SLF3_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FC_SLF1_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FC_SLF2_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FC_SLF3_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_HayBCatB$mn_FDC_SLF, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FDC_SLF_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FDC_SLF_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FDC_SLF1_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FDC_SLF2_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_HayBCatB$mn_FDC_SLF3_L, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FDC_SLF1_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FDC_SLF2_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_HayBCatB$mn_FDC_SLF3_R, SLF_data_HayBCatB$HayBCatB.Raw, SLF_data_HayBCatB[,c("Age", "Sex_binary")])

#map out all correlations - but need to drop NAs...
SLF_data_TrailsA <- SLF_data[c("Group", "Age", "Sex", "Sex_binary", "TrailsA.Raw", "mn_FD_SLF")]
SLF_data_TrailsA <- SLF_data_TrailsA %>% drop_na(TrailsA.Raw)
SLF_data_TrailsA$Sex_binary <- as.numeric(SLF_data_TrailsA$Sex_binary)

gp.cor <- SLF_data_TrailsA %>%
    split(.$Group) %>%  
    map(~corr.test(x = .x %>% dplyr::select(mn_FD_SLF, TrailsA.Raw),
                   use = "everything",
                   method = "pearson",
                   adjust = "none",
                   alpha = 0.05,
                   ci = TRUE, minlength = 5)
    )

map(gp.cor, ~.x$r)



# Analysis by clinical site as the independent variable - test for clinical site. 

#convert clinsite names to a factor variable
# SLF_data$ClinSite_name <- as.factor(SLF_data$ClinSite_name)
# 
# #run ANOVA to see if there are significant differences between groups
# #for FD
# clinsite_SLF_FD_mod <- lm(mn_FD_SLF ~ ClinSite_name, data = SLF_data)
# clinsite_SLF_FD_mod_L <- lm(mn_FD_SLF_L ~ ClinSite_name, data = SLF_data)
# clinsite_SLF_FD_mod_R <- lm(mn_FD_SLF_R ~ ClinSite_name, data = SLF_data)
# anova(clinsite_SLF_FD_mod)
# anova(clinsite_SLF_FD_mod_L)
# anova(clinsite_SLF_FD_mod_R)
# #for FC
# clinsite_SLF_FC_mod <- lm(mn_FC_SLF ~ ClinSite_name, data = SLF_data)
# clinsite_SLF_FC_mod_L <- lm(mn_FC_SLF_L ~ ClinSite_name, data = SLF_data)
# clinsite_SLF_FC_mod_R <- lm(mn_FC_SLF_R ~ ClinSite_name, data = SLF_data)
# anova(clinsite_SLF_FC_mod)
# anova(clinsite_SLF_FC_mod_L)
# anova(clinsite_SLF_FC_mod_R)
# #for FDC
# clinsite_SLF_FDC_mod <- lm(mn_FDC_SLF ~ ClinSite_name, data = SLF_data)
# clinsite_SLF_FDC_mod_L <- lm(mn_FDC_SLF_L ~ ClinSite_name, data = SLF_data)
# clinsite_SLF_FDC_mod_R <- lm(mn_FDC_SLF_R ~ ClinSite_name, data = SLF_data)
# anova(clinsite_SLF_FDC_mod)
# anova(clinsite_SLF_FDC_mod_L)
# anova(clinsite_SLF_FDC_mod_R)

# #run pairwise comparisons, given that the F-test was significant. 
# #for FD
# post_hoc_clinsite_SLF_FD_mod <- glht(clinsite_SLF_FD_mod, linfct = mcp(ClinSite_name = "Tukey"))
# summary(post_hoc_clinsite_SLF_FD_mod)
# confint(post_hoc_clinsite_SLF_FD_mod)
# post_hoc_clinsite_SLF_FD_mod_L <- glht(clinsite_SLF_FD_mod_L, linfct = mcp(ClinSite_name = "Tukey"))
# summary(post_hoc_clinsite_SLF_FD_mod_L)
# confint(post_hoc_clinsite_SLF_FD_mod_L)
# post_hoc_clinsite_SLF_FD_mod_R <- glht(clinsite_SLF_FD_mod_R, linfct = mcp(ClinSite_name = "Tukey"))
# summary(post_hoc_clinsite_SLF_FD_mod_R)
# confint(post_hoc_clinsite_SLF_FD_mod_R)
# #for FC
# post_hoc_clinsite_SLF_FC_mod_L <- glht(clinsite_SLF_FC_mod_L, linfct = mcp(ClinSite_name = "Tukey"))
# summary(post_hoc_clinsite_SLF_FC_mod_L)
# confint(post_hoc_clinsite_SLF_FC_mod_L)
# #for FDC
# post_hoc_clinsite_SLF_FDC_mod_R <- glht(clinsite_SLF_FDC_mod_R, linfct = mcp(ClinSite_name = "Tukey"))
# summary(post_hoc_clinsite_SLF_FDC_mod_R)
# confint(post_hoc_clinsite_SLF_FDC_mod_R)
# 
# #plot data
# #whole SLF FD (raincloud plot)
# ggplot(SLF_data, aes(x = ClinSite_name, y = mn_FD_SLF, fill = ClinSite_name)) + 
#     geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_point(aes(y = mn_FD_SLF, color = ClinSite_name), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = ClinSite_name)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = ClinSite_name)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FD)") +
#     theme_classic() +
#     theme(legend.position = "none") +
#     coord_flip()




# # Analysis for 3 groups (Older adults (C(1), SCD(2)) vs. MCI (aMCI(3), mMCI(4)) vs. AD(5))
# 
# #add in group data and combine the groups from 5 groups into 3 groups
# SLF_data$Three_Groups <- covariates_data$Group
# 
# SLF_data$Three_Groups[SLF_data$Three_Groups == "2"] <- "1" 
# SLF_data$Three_Groups[SLF_data$Three_Groups == "4"] <- "3" 
# 
# 
# #look at descriptive stats of the whole SLF FBA metrics between groups
# SLF_FD_3groups_descrip <- describeBy(SLF_data$mn_FD_SLF, SLF_data$Three_Groups)
# SLF_FC_3groups_descrip <- describeBy(SLF_data$mn_FC_SLF, SLF_data$Three_Groups)
# SLF_FDC_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF, SLF_data$Three_Groups)
# 
# #look at descriptive stats of the Left SLF FBA metrics between groups
# SLF_FD_L_3groups_descrip <- describeBy(SLF_data$mn_FD_SLF_L, SLF_data$Three_Groups)
# SLF_FC_L_3groups_descrip <- describeBy(SLF_data$mn_FC_SLF_L, SLF_data$Three_Groups)
# SLF_FDC_L_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF_L, SLF_data$Three_Groups)
# 
# 
# #look at descriptive stats of the Right SLF FBA metrics between groups
# SLF_FD_R_3groups_descrip <- describeBy(SLF_data$mn_FD_SLF_R, SLF_data$Three_Groups)
# SLF_FC_R_3groups_descrip <- describeBy(SLF_data$mn_FC_SLF_R, SLF_data$Three_Groups)
# SLF_FDC_R_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF_R, SLF_data$Three_Groups)
# 
# 
# #run ANOVA to see if there are significant differences between groups
# #for FD
# #mean
# SLF_FD_3groups_mod <- lm(mn_FD_SLF ~ Three_Groups, data = SLF_data)
# SLF_FD_3groups_mod_L <- lm(mn_FD_SLF_L ~ Three_Groups, data = SLF_data)
# SLF_FD_3groups_mod_R <- lm(mn_FD_SLF_R ~ Three_Groups, data = SLF_data)
# #run ANOVA
# anova(SLF_FD_3groups_mod)
# anova(SLF_FD_3groups_mod_L)
# anova(SLF_FD_3groups_mod_R)
# #run Bayesian ANOVA
# anovaBF(mn_FD_SLF ~ Three_Groups, data = SLF_data) 
# anovaBF(mn_FD_SLF_L ~ Three_Groups, data = SLF_data) 
# anovaBF(mn_FD_SLF_R ~ Three_Groups, data = SLF_data) 
# #calculate the effect size (eta-squared)
# etaSquared(SLF_FD_3groups_mod)
# etaSquared(SLF_FD_3groups_mod_L)
# etaSquared(SLF_FD_3groups_mod_R)
# 
# #for FC
# #mean
# SLF_FC_3groups_mod <- lm(mn_FC_SLF ~ Three_Groups, data = SLF_data)
# SLF_FC_3groups_mod_L <- lm(mn_FC_SLF_L ~ Three_Groups, data = SLF_data)
# SLF_FC_3groups_mod_R <- lm(mn_FC_SLF_R ~ Three_Groups, data = SLF_data)
# #run ANOVA
# anova(SLF_FC_3groups_mod)
# anova(SLF_FC_3groups_mod_L)
# anova(SLF_FC_3groups_mod_R)
# #run Bayesian ANOVA
# anovaBF(mn_FC_SLF ~ Three_Groups, data = SLF_data) 
# anovaBF(mn_FC_SLF_L ~ Three_Groups, data = SLF_data) 
# anovaBF(mn_FC_SLF_R ~ Three_Groups, data = SLF_data) 
# 
# #for FDC
# #mean
# SLF_FDC_3groups_mod <- lm(mn_FDC_SLF ~ Three_Groups, data = SLF_data)
# SLF_FDC_3groups_mod_L <- lm(mn_FDC_SLF_L ~ Three_Groups, data = SLF_data)
# SLF_FDC_3groups_mod_R <- lm(mn_FDC_SLF_R ~ Three_Groups, data = SLF_data)
# #run ANOVA
# anova(SLF_FDC_3groups_mod)
# anova(SLF_FDC_3groups_mod_L)
# anova(SLF_FDC_3groups_mod_R)
# #run Bayesian ANOVA
# anovaBF(mn_FDC_SLF ~ Three_Groups, data = SLF_data) 
# anovaBF(mn_FDC_SLF_L ~ Three_Groups, data = SLF_data) 
# anovaBF(mn_FDC_SLF_R ~ Three_Groups, data = SLF_data) 
# 
# #run pairwise comparisons, given that the F-test was significant. 
# post_hoc_SLF_FD_3groups_mod <- glht(SLF_FD_3groups_mod, linfct = mcp(Three_Groups = "Tukey"))
# summary(post_hoc_SLF_FD_3groups_mod)
# confint(post_hoc_SLF_FD_3groups_mod)
# 
# post_hoc_SLF_FD_3groups_mod_L <- glht(SLF_FD_3groups_mod_L, linfct = mcp(Three_Groups = "Tukey"))
# summary(post_hoc_SLF_FD_3groups_mod_L)
# confint(post_hoc_SLF_FD_3groups_mod_L)
# 
# post_hoc_SLF_FD_3groups_mod_R <- glht(SLF_FD_3groups_mod_R, linfct = mcp(Three_Groups = "Tukey"))
# summary(post_hoc_SLF_FD_3groups_mod_R)
# confint(post_hoc_SLF_FD_3groups_mod_R)
# 
# #conduct power analysis for the whole FD
# SLF_FD_3groups_means <- c(SLF_FD_3groups_descrip$`1`$mean, SLF_FD_3groups_descrip$`2`$mean, SLF_FD_3groups_descrip$`3`$mean, SLF_FD_3groups_descrip$`4`$mean, SLF_FD_3groups_descrip$`5`$mean)
# power_SLF_FD_3groups_n <- power.anova.test(groups = length(SLF_FD_3groups_means), between.var = anova(SLF_FD_3groups_mod)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod)$`Sum Sq`[2], power = .8, sig.level = 0.05)
# power_SLF_FD_3groups_power<- power.anova.test(groups = length(SLF_FD_3groups_means), between.var = anova(SLF_FD_3groups_mod)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod)$`Sum Sq`[2], n = 41, sig.level = 0.05)
# #conduct power analysis for the left FD
# SLF_FD_L_3groups_means <- c(SLF_FD_L_3groups_descrip$`1`$mean, SLF_FD_L_3groups_descrip$`2`$mean, SLF_FD_L_3groups_descrip$`3`$mean, SLF_FD_L_3groups_descrip$`4`$mean, SLF_FD_L_3groups_descrip$`5`$mean)
# power_SLF_FD_L_3groups_n <- power.anova.test(groups = length(SLF_FD_L_3groups_means), between.var = anova(SLF_FD_3groups_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod_L)$`Sum Sq`[2], power = .8, sig.level = 0.05)
# power_SLF_FD_L_3groups_power<- power.anova.test(groups = length(SLF_FD_L_3groups_means), between.var = anova(SLF_FD_3groups_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod_L)$`Sum Sq`[2], n = 41, sig.level = 0.05)
# #conduct power analysis for the right FD
# SLF_FD_R_3groups_means <- c(SLF_FD_R_3groups_descrip$`1`$mean, SLF_FD_R_3groups_descrip$`2`$mean, SLF_FD_R_3groups_descrip$`3`$mean, SLF_FD_R_3groups_descrip$`4`$mean, SLF_FD_R_3groups_descrip$`5`$mean)
# power_SLF_FD_R_3groups_n <- power.anova.test(groups = length(SLF_FD_R_3groups_means), between.var = anova(SLF_FD_3groups_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod_R)$`Sum Sq`[2], power = .8, sig.level = 0.05)
# power_SLF_FD_R_3groups_power<- power.anova.test(groups = length(SLF_FD_R_3groups_means), between.var = anova(SLF_FD_3groups_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FD_3groups_mod_R)$`Sum Sq`[2], n = 41, sig.level = 0.05)
# 
# 
# 
# #plot data
# #for FD - mean
# #whole SLF FD (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FD)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #whole SLF FD (raincloud plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF, fill = Three_Groups)) + 
#     geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_point(aes(y = mn_FD_SLF, color = Three_Groups), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FD)") +
#     scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     coord_flip()
# 
# #left SLF FD (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF_L)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FD)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #left SLF FD (raincloud plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF_L, fill = Three_Groups)) + 
#     geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_point(aes(y = mn_FD_SLF_L, color = Three_Groups), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FD)") +
#     scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     coord_flip()
# 
# #right SLF FD (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF_R)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FD)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #right SLF FD (raincloud plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FD_SLF_R, fill = Three_Groups)) + 
#     geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_point(aes(y = mn_FD_SLF_R, color = Three_Groups), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FD)") +
#     scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     coord_flip()
# 
# #for FC - mean
# #whole SLF FC (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FC_SLF)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Cross-section (FC)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #left SLF FC (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FC_SLF_L)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Cross-section (FC)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #right SLF FC (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FC_SLF_R)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Cross-section (FC)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #for FDC - mean
# #whole SLF FDC (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density Cross-section (FDC)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #Left SLF FDC (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF_L)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density Cross-section (FDC)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #Right SLF FDC (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF_R)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density Cross-section (FDC)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)








####------------------ Longitudinal analysis (F0 vs. F2) -------------------####

DPRC_neuropsych_data <- read.csv("longitudinal_DPRC_neuropsych_data_lined_up_valid_participants.csv")
#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'

#Add in longitudinal values for participants
Individual_number <- c(1:124, 1:124)
DPRC_neuropsych_data$Individual_number <- Individual_number

#convert variables
DPRC_neuropsych_data$ParticipantID <- as.factor(DPRC_neuropsych_data$ParticipantID)
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Sex_binary <- as.factor(DPRC_neuropsych_data$Sex_binary)
DPRC_neuropsych_data$Sex <- as.factor(DPRC_neuropsych_data$Sex)
DPRC_neuropsych_data$Timepoint <- as.factor(DPRC_neuropsych_data$Timepoint)
DPRC_neuropsych_data$Individual_number <- as.factor(DPRC_neuropsych_data$Individual_number)


#navigate to the correct pathway which contains the SLF metric text files: 
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/template/TOI')

#read in text file of data
SLF_data <- cbind.data.frame(read.table("FD_SLF_whole_TOI.txt", header = T), 
                             read.table("FC_log_SLF_whole_TOI.txt", header = T), 
                             read.table("FDC_SLF_whole_TOI.txt", header = T), 
                             read.table("FD_SLF_L_TOI.txt", header = T), 
                             read.table("FC_log_SLF_L_TOI.txt", header = T), 
                             read.table("FDC_SLF_L_TOI.txt", header = T), 
                             read.table("FD_SLF_R_TOI.txt", header = T), 
                             read.table("FC_log_SLF_R_TOI.txt", header = T), 
                             read.table("FDC_SLF_R_TOI.txt", header = T), 
                             read.table("FD_SLF1_L_TOI.txt", header = T), 
                             read.table("FC_log_SLF1_L_TOI.txt", header = T), 
                             read.table("FDC_SLF1_L_TOI.txt", header = T), 
                             read.table("FD_SLF2_L_TOI.txt", header = T), 
                             read.table("FC_log_SLF2_L_TOI.txt", header = T), 
                             read.table("FDC_SLF2_L_TOI.txt", header = T), 
                             read.table("FD_SLF3_L_TOI.txt", header = T), 
                             read.table("FC_log_SLF3_L_TOI.txt", header = T), 
                             read.table("FDC_SLF3_L_TOI.txt", header = T), 
                             read.table("FD_SLF1_R_TOI.txt", header = T), 
                             read.table("FC_log_SLF1_R_TOI.txt", header = T), 
                             read.table("FDC_SLF1_R_TOI.txt", header = T), 
                             read.table("FD_SLF2_R_TOI.txt", header = T), 
                             read.table("FC_log_SLF2_R_TOI.txt", header = T), 
                             read.table("FDC_SLF2_R_TOI.txt", header = T), 
                             read.table("FD_SLF3_R_TOI.txt", header = T), 
                             read.table("FC_log_SLF3_R_TOI.txt", header = T), 
                             read.table("FDC_SLF3_R_TOI.txt", header = T))
#rename columns with each FBA metric
colnames(SLF_data) <- c("mn_FD_SLF", "md_FD_SLF", "std_FD_SLF", "std_rv_FD_SLF", "min_FD_SLF", "max_FD_SLF", "count_FD_SLF", 
                        "mn_FC_SLF", "md_FC_SLF", "std_FC_SLF", "std_rv_FC_SLF", "min_FC_SLF", "max_FC_SLF", "count_FC_SLF", 
                        "mn_FDC_SLF", "md_FDC_SLF", "std_FDC_SLF", "std_rv_FDC_SLF", "min_FDC_SLF", "max_FDC_SLF", "count_FDC_SLF", 
                        "mn_FD_SLF_L", "md_FD_SLF_L", "std_FD_SLF_L", "std_rv_FD_SLF_L", "min_FD_SLF_L", "max_FD_SLF_L", 
                        "count_FD_SLF_L", "mn_FC_SLF_L", "md_FC_SLF_L", "std_FC_SLF_L", "std_rv_FC_SLF_L", "min_FC_SLF_L", 
                        "max_FC_SLF_L", "count_FC_SLF_L", "mn_FDC_SLF_L", "md_FDC_SLF_L", "std_FDC_SLF_L", "std_rv_FDC_SLF_L", 
                        "min_FDC_SLF_L", "max_FDC_SLF_L", "count_FDC_SLF_L", "mn_FD_SLF_R", "md_FD_SLF_R", "std_FD_SLF_R", 
                        "std_rv_FD_SLF_R", "min_FD_SLF_R", "max_FD_SLF_R", "count_FD_SLF_R","mn_FC_SLF_R", "md_FC_SLF_R", 
                        "std_FC_SLF_R", "std_rv_FC_SLF_R", "min_FC_SLF_R", "max_FC_SLF_R", "count_FC_SLF_R", "mn_FDC_SLF_R", 
                        "md_FDC_SLF_R", "std_FDC_SLF_R", "std_rv_FDC_SLF_R", "min_FDC_SLF_R", "max_FDC_SLF_R", "count_FDC_SLF_R", 
                        "mn_FD_SLF1_L", "md_FD_SLF1_L", "std_FD_SLF1_L", "std_rv_FD_SLF1_L", "min_FD_SLF1_L", "max_FD_SLF1_L", 
                        "count_FD_SLF1_L", "mn_FC_SLF1_L", "md_FC_SLF1_L", "std_FC_SLF1_L", "std_rv_FC_SLF1_L", "min_FC_SLF1_L", 
                        "max_FC_SLF1_L", "count_FC_SLF1_L", "mn_FDC_SLF1_L", "md_FDC_SLF1_L", "std_FDC_SLF1_L", "std_rv_FDC_SLF1_L", 
                        "min_FDC_SLF1_L", "max_FDC_SLF1_L", "count_FDC_SLF1_L", "mn_FD_SLF2_L", "md_FD_SLF2_L", "std_FD_SLF2_L", 
                        "std_rv_FD_SLF2_L", "min_FD_SLF2_L", "max_FD_SLF2_L", "count_FD_SLF2_L", "mn_FC_SLF2_L", "md_FC_SLF2_L", 
                        "std_FC_SLF2_L", "std_rv_FC_SLF2_L", "min_FC_SLF2_L", "max_FC_SLF2_L", "count_FC_SLF2_L", "mn_FDC_SLF2_L", 
                        "md_FDC_SLF2_L", "std_FDC_SLF2_L", "std_rv_FDC_SLF2_L", "min_FDC_SLF2_L", "max_FDC_SLF2_L", "count_FDC_SLF2_L",
                        "mn_FD_SLF3_L", "md_FD_SLF3_L", "std_FD_SLF3_L", "std_rv_FD_SLF3_L", "min_FD_SLF3_L", "max_FD_SLF3_L", 
                        "count_FD_SLF3_L", "mn_FC_SLF3_L", "md_FC_SLF3_L", "std_FC_SLF3_L", "std_rv_FC_SLF3_L", "min_Fc_SLF3_L", 
                        "max_FC_SLF3_L", "count_FC_SLF3_L", "mn_FDC_SLF3_L", "md_FDC_SLF3_L", "std_FDC_SLF3_L", "std_rv_FDC_SLF3_L", 
                        "min_FDC_SLF3_L", "max_FDC_SLF3_L", "count_FDC_SLF3_L", "mn_FD_SLF1_R", "md_FD_SLF1_R", "std_FD_SLF1_R", 
                        "std_rv_FD_SLF1_R", "min_FD_SLF1_R", "max_FD_SLF1_R", "count_FD_SLF1_R", "mn_FC_SLF1_R", "md_FC_SLF1_R", 
                        "std_FC_SLF1_R", "std_rv_FC_SLF1_R", "min_FC_SLF1_R", "max_FC_SLF1_R", "count_FC_SLF1_R", "mn_FDC_SLF1_R", 
                        "md_FDC_SLF1_R", "std_FDC_SLF1_R", "std_rv_FDC_SLF1_R", "min_FDC_SLF1_R", "max_FDC_SLF1_R", "count_FDC_SLF1_R", 
                        "mn_FD_SLF2_R", "md_FD_SLF2_R", "std_FD_SLF2_R", "std_rv_FD_SLF2_R", "min_FD_SLF2_R", "max_FD_SLF2_R", 
                        "count_FD_SLF2_R","mn_FC_SLF2_R", "md_FC_SLF2_R", "std_FC_SLF2_R", "std_rv_FC_SLF2_R", "min_FC_SLF2_R", 
                        "max_FC_SLF2_R", "count_FC_SLF2_R", "mn_FDC_SLF2_R", "md_FDC_SLF2_R", "std_FDC_SLF2_R", "std_rv_FDC_SLF2_R", 
                        "min_FDC_SLF2_R", "max_FDC_SLF2_R", "count_FDC_SLF2_R", "mn_FD_SLF3_R", "md_FD_SLF3_R", "std_FD_SLF3_R", 
                        "std_rv_FD_SLF3_R", "min_FD_SLF3_R", "max_FD_SLF3_R", "count_FD_SLF3_R","mn_FC_SLF3_R", "md_FC_SLF3_R", 
                        "std_FC_SLF3_R", "std_rv_FC_SLF3_R", "min_FC_SLF3_R", "max_FC_SLF3_R", "count_FC_SLF3_R", "mn_FDC_SLF3_R",
                        "md_FDC_SLF3_R", "std_FDC_SLF3_R", "std_rv_FDC_SLF3_R", "min_FDC_SLF3_R", "max_FDC_SLF3_R", "count_FDC_SLF3_R")

#add in Timepoint (alternating pattern for SLF analysis)
Timepoint <- rep(c('F0','F2'), 124)
SLF_data$Timepoint<-Timepoint
#add in Participant ID in alternating order for SLF analysis
ParticipantID_SLF <- sort(as.character(DPRC_neuropsych_data$ParticipantID))
SLF_data$ParticipantID_SLF<-ParticipantID_SLF
#sort dataframe by Time Point + PT ID
SLF_data <- SLF_data[order(SLF_data$Timepoint),]
#rename Participant ID 
names(SLF_data)[names(SLF_data)=='ParticipantID_SLF'] <- 'ParticipantID'
#convert to factor variables
SLF_data$ParticipantID <- as.factor(SLF_data$ParticipantID)
SLF_data$Timepoint <- as.factor(SLF_data$Timepoint)
#add in Individual number
SLF_data$Individual_number <- DPRC_neuropsych_data$Individual_number
#add in Group classification column for each participant - data from another worksheet
SLF_data$Group <- DPRC_neuropsych_data$Group
#add in covariates (clinical site) from the covariates_data dataframe to the SLF dataframe
SLF_data$ClinSite_name <- DPRC_neuropsych_data$Clinical_site 
SLF_data$Age <- DPRC_neuropsych_data$Age 
SLF_data$Sex <- DPRC_neuropsych_data$Sex
SLF_data$Sex_binary <- DPRC_neuropsych_data$Sex_binary
#add in neuropsych variables
neuropsych_test_names <- c("TrailsA.Raw","TrailsA.Z","TrailsB.Raw","TrailsB.Z","ColorNaming.Raw","ColorNaming.Z","WordReading.Raw","WordReading.Z","Inhibition.Raw","Inhibition.Z","LetFluency.Raw","LetFluency.Z","CatFluency.Raw","CatFluency.Z","Switching.Raw","Switching.z","HayBTime1.Raw","HayBTime1.z","HayBTime2.Raw","HayBTime2.z","HayBCatA.Raw","HayBCatA.z","HayBCatB.Raw","HayBCatB.z")
SLF_data[neuropsych_test_names] <- DPRC_neuropsych_data[neuropsych_test_names]


###----------------------------Descriptives----------------------------------####
#look at descriptive stats of the whole SLF FBA metrics between groups for F0 
#and F2 time points
    SLF_FD_longit_descrip <- describeBy(SLF_data$mn_FD_SLF, list(SLF_data$Group, SLF_data$Timepoint))
    SLF_FC_longit_descrip <- describeBy(SLF_data$mn_FC_SLF, list(SLF_data$Group, SLF_data$Timepoint))
    SLF_FDC_longit_descrip <- describeBy(SLF_data$mn_FDC_SLF, list(SLF_data$Group, SLF_data$Timepoint))
    #find mean & SD from total sample in F0 and F2 timepoints:
    #F0 - FD
    F0_SLF_FD <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF_FD <- na.omit(F0_SLF_FD$mn_FD_SLF)
    mean(noNAsF0_SLF_FD)
    sd(noNAsF0_SLF_FD)
    #FC
    F0_SLF_FC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF_FC <- na.omit(F0_SLF_FC$mn_FC_SLF)
    mean(noNAsF0_SLF_FC)
    sd(noNAsF0_SLF_FC)
    #FDC
    F0_SLF_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF_FDC <- na.omit(F0_SLF_FDC$mn_FDC_SLF)
    mean(noNAsF0_SLF_FDC)
    sd(noNAsF0_SLF_FDC)
    #F2 - FD
    F2_SLF_FD <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF_FD <- na.omit(F2_SLF_FD$mn_FD_SLF)
    mean(noNAsF2_SLF_FD)
    sd(noNAsF2_SLF_FD)
    #FC
    F2_SLF_FC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF_FC <- na.omit(F2_SLF_FC$mn_FC_SLF)
    mean(noNAsF2_SLF_FC)
    sd(noNAsF2_SLF_FC)
    #FDC
    F2_SLF_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF_FDC <- na.omit(F2_SLF_FDC$mn_FDC_SLF)
    mean(noNAsF2_SLF_FDC)
    sd(noNAsF2_SLF_FDC)

#For Left SLF FBA metrics
    SLF_L_FD_longit_descrip <- describeBy(SLF_data$mn_FD_SLF_L, list(SLF_data$Group, SLF_data$Timepoint))
    SLF_L_FC_longit_descrip <- describeBy(SLF_data$mn_FC_SLF_L, list(SLF_data$Group, SLF_data$Timepoint))
    SLF_L_FDC_longit_descrip <- describeBy(SLF_data$mn_FDC_SLF_L, list(SLF_data$Group, SLF_data$Timepoint))
    #find mean & SD from total sample in F0 and F2 timepoints:
    #F0 - FD
    F0_SLF_L_FD <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF_L_FD <- na.omit(F0_SLF_L_FD$mn_FD_SLF_L)
    mean(noNAsF0_SLF_L_FD)
    sd(noNAsF0_SLF_L_FD)
    #FC
    F0_SLF_L_FC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF_L_FC <- na.omit(F0_SLF_L_FC$mn_FC_SLF_L)
    mean(noNAsF0_SLF_L_FC)
    sd(noNAsF0_SLF_L_FC)
    #FDC
    F0_SLF_L_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF_L_FDC <- na.omit(F0_SLF_L_FDC$mn_FDC_SLF_L)
    mean(noNAsF0_SLF_L_FDC)
    sd(noNAsF0_SLF_L_FDC)
    #F2 - FD
    F2_SLF_L_FD <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF_L_FD <- na.omit(F2_SLF_L_FD$mn_FD_SLF_L)
    mean(noNAsF2_SLF_L_FD)
    sd(noNAsF2_SLF_L_FD)
    #FC
    F2_SLF_L_FC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF_L_FC <- na.omit(F2_SLF_L_FC$mn_FC_SLF_L)
    mean(noNAsF2_SLF_L_FC)
    sd(noNAsF2_SLF_L_FC)
    #FDC
    F2_SLF_L_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF_L_FDC <- na.omit(F2_SLF_L_FDC$mn_FDC_SLF_L)
    mean(noNAsF2_SLF_L_FDC)
    sd(noNAsF2_SLF_L_FDC)

#For Right SLF FBA metrics
    SLF_R_FD_longit_descrip <- describeBy(SLF_data$mn_FD_SLF_R, list(SLF_data$Group, SLF_data$Timepoint))
    SLF_R_FC_longit_descrip <- describeBy(SLF_data$mn_FC_SLF_R, list(SLF_data$Group, SLF_data$Timepoint))
    SLF_R_FDC_longit_descrip <- describeBy(SLF_data$mn_FDC_SLF_R, list(SLF_data$Group, SLF_data$Timepoint))
    #find mean & SD from total sample in F0 and F2 timepoints:
    #F0 - FD
    F0_SLF_R_FD <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF_R_FD <- na.omit(F0_SLF_R_FD$mn_FD_SLF_R)
    mean(noNAsF0_SLF_R_FD)
    sd(noNAsF0_SLF_R_FD)
    #FC
    F0_SLF_R_FC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF_R_FC <- na.omit(F0_SLF_R_FC$mn_FC_SLF_R)
    mean(noNAsF0_SLF_R_FC)
    sd(noNAsF0_SLF_R_FC)
    #FDC
    F0_SLF_R_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF_R_FDC <- na.omit(F0_SLF_R_FDC$mn_FDC_SLF_R)
    mean(noNAsF0_SLF_R_FDC)
    sd(noNAsF0_SLF_R_FDC)
    #F2 - FD
    F2_SLF_R_FD <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF_R_FD <- na.omit(F2_SLF_R_FD$mn_FD_SLF_R)
    mean(noNAsF2_SLF_R_FD)
    sd(noNAsF2_SLF_R_FD)
    #FC
    F2_SLF_R_FC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF_R_FC <- na.omit(F2_SLF_R_FC$mn_FC_SLF_R)
    mean(noNAsF2_SLF_R_FC)
    sd(noNAsF2_SLF_R_FC)
    #FDC
    F2_SLF_R_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF_R_FDC <- na.omit(F2_SLF_R_FDC$mn_FDC_SLF_R)
    mean(noNAsF2_SLF_R_FDC)
    sd(noNAsF2_SLF_R_FDC)
#For Left SLF 1 FBA metrics
    SLF1_L_FD_longit_descrip <- describeBy(SLF_data$mn_FD_SLF1_L, list(SLF_data$Group, SLF_data$Timepoint))
    SLF1_L_FC_longit_descrip <- describeBy(SLF_data$mn_FC_SLF1_L, list(SLF_data$Group, SLF_data$Timepoint))
    SLF1_L_FDC_longit_descrip <- describeBy(SLF_data$mn_FDC_SLF1_L, list(SLF_data$Group, SLF_data$Timepoint))
    #find mean & SD from total sample in F0 and F2 timepoints:
    #F0 - FD
    F0_SLF1_L_FD <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF1_L_FD <- na.omit(F0_SLF1_L_FD$mn_FD_SLF1_L)
    mean(noNAsF0_SLF1_L_FD)
    sd(noNAsF0_SLF1_L_FD)
    #FC
    F0_SLF1_L_FC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF1_L_FC <- na.omit(F0_SLF1_L_FC$mn_FC_SLF1_L)
    mean(noNAsF0_SLF1_L_FC)
    sd(noNAsF0_SLF1_L_FC)
    #FDC
    F0_SLF1_L_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF1_L_FDC <- na.omit(F0_SLF1_L_FDC$mn_FDC_SLF1_L)
    mean(noNAsF0_SLF1_L_FDC)
    sd(noNAsF0_SLF1_L_FDC)
    #F2 - FD
    F2_SLF1_L_FD <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF1_L_FD <- na.omit(F2_SLF1_L_FD$mn_FD_SLF1_L)
    mean(noNAsF2_SLF1_L_FD)
    sd(noNAsF2_SLF1_L_FD)
    #FC
    F2_SLF1_L_FC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF1_L_FC <- na.omit(F2_SLF1_L_FC$mn_FC_SLF1_L)
    mean(noNAsF2_SLF1_L_FC)
    sd(noNAsF2_SLF1_L_FC)
    #FDC
    F2_SLF1_L_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF1_L_FDC <- na.omit(F2_SLF1_L_FDC$mn_FDC_SLF1_L)
    mean(noNAsF2_SLF1_L_FDC)
    sd(noNAsF2_SLF1_L_FDC)
#For Left SLF 2 FBA metrics
    SLF2_L_FD_longit_descrip <- describeBy(SLF_data$mn_FD_SLF2_L, list(SLF_data$Group, SLF_data$Timepoint))
    SLF2_L_FC_longit_descrip <- describeBy(SLF_data$mn_FC_SLF2_L, list(SLF_data$Group, SLF_data$Timepoint))
    SLF2_L_FDC_longit_descrip <- describeBy(SLF_data$mn_FDC_SLF2_L, list(SLF_data$Group, SLF_data$Timepoint))
    #find mean & SD from total sample in F0 and F2 timepoints:
    #F0 - FD
    F0_SLF2_L_FD <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF2_L_FD <- na.omit(F0_SLF2_L_FD$mn_FD_SLF2_L)
    mean(noNAsF0_SLF2_L_FD)
    sd(noNAsF0_SLF2_L_FD)
    #FC
    F0_SLF2_L_FC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF2_L_FC <- na.omit(F0_SLF2_L_FC$mn_FC_SLF2_L)
    mean(noNAsF0_SLF2_L_FC)
    sd(noNAsF0_SLF2_L_FC)
    #FDC
    F0_SLF2_L_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF2_L_FDC <- na.omit(F0_SLF2_L_FDC$mn_FDC_SLF2_L)
    mean(noNAsF0_SLF2_L_FDC)
    sd(noNAsF0_SLF2_L_FDC)
    #F2 - FD
    F2_SLF2_L_FD <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF2_L_FD <- na.omit(F2_SLF2_L_FD$mn_FD_SLF2_L)
    mean(noNAsF2_SLF2_L_FD)
    sd(noNAsF2_SLF2_L_FD)
    #FC
    F2_SLF2_L_FC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF2_L_FC <- na.omit(F2_SLF2_L_FC$mn_FC_SLF2_L)
    mean(noNAsF2_SLF2_L_FC)
    sd(noNAsF2_SLF2_L_FC)
    #FDC
    F2_SLF2_L_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF2_L_FDC <- na.omit(F2_SLF2_L_FDC$mn_FDC_SLF2_L)
    mean(noNAsF2_SLF2_L_FDC)
    sd(noNAsF2_SLF2_L_FDC)
#For Left SLF 3 FBA metrics
    SLF3_L_FD_longit_descrip <- describeBy(SLF_data$mn_FD_SLF3_L, list(SLF_data$Group, SLF_data$Timepoint))
    SLF3_L_FC_longit_descrip <- describeBy(SLF_data$mn_FC_SLF3_L, list(SLF_data$Group, SLF_data$Timepoint))
    SLF3_L_FDC_longit_descrip <- describeBy(SLF_data$mn_FDC_SLF3_L, list(SLF_data$Group, SLF_data$Timepoint))
    #find mean & SD from total sample in F0 and F2 timepoints:
    #F0 - FD
    F0_SLF3_L_FD <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF3_L_FD <- na.omit(F0_SLF3_L_FD$mn_FD_SLF3_L)
    mean(noNAsF0_SLF3_L_FD)
    sd(noNAsF0_SLF3_L_FD)
    #FC
    F0_SLF3_L_FC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF3_L_FC <- na.omit(F0_SLF3_L_FC$mn_FC_SLF3_L)
    mean(noNAsF0_SLF3_L_FC)
    sd(noNAsF0_SLF3_L_FC)
    #FDC
    F0_SLF3_L_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF3_L_FDC <- na.omit(F0_SLF3_L_FDC$mn_FDC_SLF3_L)
    mean(noNAsF0_SLF3_L_FDC)
    sd(noNAsF0_SLF3_L_FDC)
    #F2 - FD
    F2_SLF3_L_FD <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF3_L_FD <- na.omit(F2_SLF3_L_FD$mn_FD_SLF3_L)
    mean(noNAsF2_SLF3_L_FD)
    sd(noNAsF2_SLF3_L_FD)
    #FC
    F2_SLF3_L_FC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF3_L_FC <- na.omit(F2_SLF3_L_FC$mn_FC_SLF3_L)
    mean(noNAsF2_SLF3_L_FC)
    sd(noNAsF2_SLF3_L_FC)
    #FDC
    F2_SLF3_L_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF3_L_FDC <- na.omit(F2_SLF3_L_FDC$mn_FDC_SLF3_L)
    mean(noNAsF2_SLF3_L_FDC)
    sd(noNAsF2_SLF3_L_FDC)
#For Right SLF 1 FBA metrics
    SLF1_R_FD_longit_descrip <- describeBy(SLF_data$mn_FD_SLF1_R, list(SLF_data$Group, SLF_data$Timepoint))
    SLF1_R_FC_longit_descrip <- describeBy(SLF_data$mn_FC_SLF1_R, list(SLF_data$Group, SLF_data$Timepoint))
    SLF1_R_FDC_longit_descrip <- describeBy(SLF_data$mn_FDC_SLF1_R, list(SLF_data$Group, SLF_data$Timepoint))
    #find mean & SD from total sample in F0 and F2 timepoints:
    #F0 - FD
    F0_SLF1_R_FD <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF1_R_FD <- na.omit(F0_SLF1_R_FD$mn_FD_SLF1_R)
    mean(noNAsF0_SLF1_R_FD)
    sd(noNAsF0_SLF1_R_FD)
    #FC
    F0_SLF1_R_FC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF1_R_FC <- na.omit(F0_SLF1_R_FC$mn_FC_SLF1_R)
    mean(noNAsF0_SLF1_R_FC)
    sd(noNAsF0_SLF1_R_FC)
    #FDC
    F0_SLF1_R_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF1_R_FDC <- na.omit(F0_SLF1_R_FDC$mn_FDC_SLF1_R)
    mean(noNAsF0_SLF1_R_FDC)
    sd(noNAsF0_SLF1_R_FDC)
    #F2 - FD
    F2_SLF1_R_FD <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF1_R_FD <- na.omit(F2_SLF1_R_FD$mn_FD_SLF1_R)
    mean(noNAsF2_SLF1_R_FD)
    sd(noNAsF2_SLF1_R_FD)
    #FC
    F2_SLF1_R_FC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF1_R_FC <- na.omit(F2_SLF1_R_FC$mn_FC_SLF1_R)
    mean(noNAsF2_SLF1_R_FC)
    sd(noNAsF2_SLF1_R_FC)
    #FDC
    F2_SLF1_R_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF1_R_FDC <- na.omit(F2_SLF1_R_FDC$mn_FDC_SLF1_R)
    mean(noNAsF2_SLF1_R_FDC)
    sd(noNAsF2_SLF1_R_FDC)
#For Right SLF 2 FBA metrics
    SLF2_R_FD_longit_descrip <- describeBy(SLF_data$mn_FD_SLF2_R, list(SLF_data$Group, SLF_data$Timepoint))
    SLF2_R_FC_longit_descrip <- describeBy(SLF_data$mn_FC_SLF2_R, list(SLF_data$Group, SLF_data$Timepoint))
    SLF2_R_FDC_longit_descrip <- describeBy(SLF_data$mn_FDC_SLF2_R, list(SLF_data$Group, SLF_data$Timepoint))
    #find mean & SD from total sample in F0 and F2 timepoints:
    #F0 - FD
    F0_SLF2_R_FD <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF2_R_FD <- na.omit(F0_SLF2_R_FD$mn_FD_SLF2_R)
    mean(noNAsF0_SLF2_R_FD)
    sd(noNAsF0_SLF2_R_FD)
    #FC
    F0_SLF2_R_FC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF2_R_FC <- na.omit(F0_SLF2_R_FC$mn_FC_SLF2_R)
    mean(noNAsF0_SLF2_R_FC)
    sd(noNAsF0_SLF2_R_FC)
    #FDC
    F0_SLF2_R_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF2_R_FDC <- na.omit(F0_SLF2_R_FDC$mn_FDC_SLF2_R)
    mean(noNAsF0_SLF2_R_FDC)
    sd(noNAsF0_SLF2_R_FDC)
    #F2 - FD
    F2_SLF2_R_FD <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF2_R_FD <- na.omit(F2_SLF2_R_FD$mn_FD_SLF2_R)
    mean(noNAsF2_SLF2_R_FD)
    sd(noNAsF2_SLF2_R_FD)
    #FC
    F2_SLF2_R_FC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF2_R_FC <- na.omit(F2_SLF2_R_FC$mn_FC_SLF2_R)
    mean(noNAsF2_SLF2_R_FC)
    sd(noNAsF2_SLF2_R_FC)
    #FDC
    F2_SLF2_R_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF2_R_FDC <- na.omit(F2_SLF2_R_FDC$mn_FDC_SLF2_R)
    mean(noNAsF2_SLF2_R_FDC)
    sd(noNAsF2_SLF2_R_FDC)
#For Right SLF 3 FBA metrics
    SLF3_R_FD_longit_descrip <- describeBy(SLF_data$mn_FD_SLF3_R, list(SLF_data$Group, SLF_data$Timepoint))
    SLF3_R_FC_longit_descrip <- describeBy(SLF_data$mn_FC_SLF3_R, list(SLF_data$Group, SLF_data$Timepoint))
    SLF3_R_FDC_longit_descrip <- describeBy(SLF_data$mn_FDC_SLF3_R, list(SLF_data$Group, SLF_data$Timepoint))
    #find mean & SD from total sample in F0 and F2 timepoints:
    #F0 - FD
    F0_SLF3_R_FD <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF3_R_FD <- na.omit(F0_SLF3_R_FD$mn_FD_SLF3_R)
    mean(noNAsF0_SLF3_R_FD)
    sd(noNAsF0_SLF3_R_FD)
    #FC
    F0_SLF3_R_FC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF3_R_FC <- na.omit(F0_SLF3_R_FC$mn_FC_SLF3_R)
    mean(noNAsF0_SLF3_R_FC)
    sd(noNAsF0_SLF3_R_FC)
    #FDC
    F0_SLF3_R_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F0",]
    noNAsF0_SLF3_R_FDC <- na.omit(F0_SLF3_R_FDC$mn_FDC_SLF3_R)
    mean(noNAsF0_SLF3_R_FDC)
    sd(noNAsF0_SLF3_R_FDC)
    #F2 - FD
    F2_SLF3_R_FD <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF3_R_FD <- na.omit(F2_SLF3_R_FD$mn_FD_SLF3_R)
    mean(noNAsF2_SLF3_R_FD)
    sd(noNAsF2_SLF3_R_FD)
    #FC
    F2_SLF3_R_FC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF3_R_FC <- na.omit(F2_SLF3_R_FC$mn_FC_SLF3_R)
    mean(noNAsF2_SLF3_R_FC)
    sd(noNAsF2_SLF3_R_FC)
    #FDC
    F2_SLF3_R_FDC <- SLF_data[SLF_data[, "Timepoint"] == "F2",]
    noNAsF2_SLF3_R_FDC <- na.omit(F2_SLF3_R_FDC$mn_FDC_SLF3_R)
    mean(noNAsF2_SLF3_R_FDC)
    sd(noNAsF2_SLF3_R_FDC)
    
####--------------------- ANOVA & ANCOVA (age & sex) ------------------------####
    
#run ANOVA to see if there are significant differences between groups
#run mixed design, 2 x 5 ANOVA for whole FBA
#for FD
aov_SLF_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF_FD_mod)
#aov_SLF_FD_mod <- aov(mn_FD_SLF ~ Group*Timepoint + Error(ParticipantID/Timepoint), data=SLF_data)
#summary(aov_SLF_FD_mod)
#aov_SLF_FD_mod <- Anova(lm(mn_FD_SLF ~ Group*Timepoint, data=SLF_data), type = "III") #use this anova test to account for unbalanced designs/sample sizes
#aov_SLF_FD_mod
#for FC
aov_SLF_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF_FC_mod)
#FDC
aov_SLF_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF_FDC_mod)
#for left SLF
aov_SLF_L_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF_L_FD_mod)
aov_SLF_L_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF_L_FC_mod)
aov_SLF_L_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF_L_FDC_mod)
#for right SLF
aov_SLF_R_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF_R_FD_mod)
aov_SLF_R_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF_R_FC_mod)
aov_SLF_R_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF_R_FDC_mod)
#for left SLF 1
aov_SLF1_L_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF1_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_L_FD_mod)
aov_SLF1_L_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF1_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_L_FC_mod)
aov_SLF1_L_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF1_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_L_FDC_mod)
#for left SLF 2
aov_SLF2_L_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF2_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_FD_mod) 
aov_SLF2_L_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF2_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_FC_mod) 
aov_SLF2_L_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF2_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_FDC_mod)  
#for left SLF 3
aov_SLF3_L_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF3_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_FD_mod) 
aov_SLF3_L_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF3_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_FC_mod) 
aov_SLF3_L_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF3_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_FDC_mod) 
#for right SLF 1
aov_SLF1_R_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF1_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_R_FD_mod) 
aov_SLF1_R_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF1_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_R_FC_mod)
aov_SLF1_R_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF1_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_R_FDC_mod)
#for right SLF 2
aov_SLF2_R_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF2_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_R_FD_mod)
aov_SLF2_R_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF2_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_R_FC_mod)
aov_SLF2_R_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF2_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_R_FDC_mod)
#for right SLF 3
aov_SLF3_R_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF3_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_R_FD_mod)
aov_SLF3_R_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF3_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_R_FC_mod)
aov_SLF3_R_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF3_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_R_FDC_mod)

#follow up tests with sig. main effects
#Whole SLF
#Whole SLF FD by Group
aov(mn_FD_SLF ~ Group, data = SLF_data) %>% tukey_hsd()
#effect size for Group
SLF_data%>%cohens_d(mn_FD_SLF ~ Group)
#Whole SLF FC by Timepoint
aov(mn_FC_SLF ~ Timepoint, data = SLF_data) %>% tukey_hsd()
#Left SLF
#Left SLF FD by Group
aov(mn_FD_SLF_L ~ Group, data = SLF_data) %>% tukey_hsd()
#effect size for Group
SLF_data%>%cohens_d(mn_FD_SLF_L ~ Group)
#Left SLF FC by Timepoint
aov(mn_FC_SLF_L ~ Timepoint, data = SLF_data) %>% tukey_hsd()
#Right SLF
#Right SLF FD by Group
aov(mn_FD_SLF_R ~ Group, data = SLF_data) %>% tukey_hsd()
#effect size for Group
SLF_data%>%cohens_d(mn_FD_SLF_R ~ Group)
#Right SLF FC by Timepoint
aov(mn_FC_SLF_R ~ Timepoint, data = SLF_data) %>% tukey_hsd()
#Left SLF1
#Left SLF1 FC by Timepoint
aov(mn_FC_SLF1_L ~ Timepoint, data = SLF_data) %>% tukey_hsd()
#Left SLF1 FDC by Timepoint
aov(mn_FDC_SLF1_L ~ Timepoint, data = SLF_data) %>% tukey_hsd()
#Left SLF2
#Left SLF2 FD by Group
aov(mn_FD_SLF2_L ~ Group, data = SLF_data) %>% tukey_hsd()
#effect size for Group
SLF_data%>%cohens_d(mn_FD_SLF2_L ~ Group)
#Left SLF2 FC by Timepoint
aov(mn_FC_SLF2_L ~ Timepoint, data = SLF_data) %>% tukey_hsd()
#Right SLF1
#Right SLF1 FD by Group
aov(mn_FD_SLF1_R ~ Group, data = SLF_data) %>% tukey_hsd()
#effect size for Group
SLF_data%>%cohens_d(mn_FD_SLF1_R ~ Group)
#Right SLF FC by Timepoint
aov(mn_FC_SLF1_R ~ Timepoint, data = SLF_data) %>% tukey_hsd()
#Right SLF1 FDC by Timepoint
aov(mn_FDC_SLF1_L ~ Timepoint, data = SLF_data) %>% tukey_hsd()
#Right SLF2
#Right SLF2 FD by Group
aov(mn_FD_SLF2_R ~ Group, data = SLF_data) %>% tukey_hsd()
#effect size for Group
SLF_data%>%cohens_d(mn_FD_SLF2_R ~ Group)
#Right SLF3
#Right SLF3 FD by Group
aov(mn_FD_SLF3_R ~ Group, data = SLF_data) %>% tukey_hsd()
#effect size for Group
SLF_data%>%cohens_d(mn_FD_SLF3_R ~ Group)


# #view interaction plot - by Group
# SLF_data %>%
#     group_by(Group,Timepoint) %>%
#     summarise(s_mean=mean(mn_FD_SLF1_R)) %>%
#     ggplot(aes(y=s_mean,x=Group,colour=Timepoint,group=Timepoint))+
#     geom_point()+geom_line()+
#     ylab("SLF1 Right Fibre Cross-section (FC)") +
#     theme_classic()
# #view interaction plot - by Timepoint
# SLF_data %>%
#     group_by(Group,Timepoint) %>%
#     summarise(s_mean=mean(mn_FD_SLF1_R)) %>%
#     ggplot(aes(y=s_mean,x=Timepoint,colour=Group,group=Group))+
#     geom_point()+geom_line()+
#     ylab("SLF1 Right  Fibre Cross-section (FC)") +
#     theme_classic()
# #Left SLF 2 FD
# #non-sig. interaction - test by Time Point
# #Test by Time Point
# aov(mn_FD_SLF2_L ~ Timepoint, data = SLF_data) %>% tukey_hsd()
# # SLF_data %>%
# #     pairwise_t_test(
# #         mn_FD_SLF2_L ~ Timepoint, paired = TRUE,
# #         p.adjust.method = "fdr"
# #     )
# #effect size for Time Point
# SLF_data%>%cohens_d(mn_FD_SLF2_L~Timepoint, paired=TRUE)
# #Left SLF 2 FDC
# #non-sig. interaction - test by Time Point
# aov(mn_FDC_SLF2_L ~ Timepoint, data = SLF_data) %>% tukey_hsd()
# #effect size for Time Point
# SLF_data%>%cohens_d(mn_FDC_SLF2_L~Timepoint, paired=TRUE)
# 
# #follow up for sig. interaction 
# #Right SLF 1
# #interaction post hoc follow-up
# #view interaction plot - by Group
# SLF_data %>%
#     group_by(Group,Timepoint) %>%
#     summarise(s_mean=mean(mn_FD_SLF1_R)) %>%
#     ggplot(aes(y=s_mean,x=Group,colour=Timepoint,group=Timepoint))+
#     geom_point()+geom_line()+
#     ylab("SLF1 Right Fibre Density (FD)") +
#     theme_classic()
# #view interaction plot - by Timepoint
# SLF_data %>%
#     group_by(Group,Timepoint) %>%
#     summarise(s_mean=mean(mn_FD_SLF1_R)) %>%
#     ggplot(aes(y=s_mean,x=Timepoint,colour=Group,group=Group))+
#     geom_point()+geom_line()+
#     ylab("SLF1 Right Fibre Density (FD)") +
#     theme_classic()
# #Sig. interaction
#Simple main effect w/ interaction - for Group
# posthoc_ME_Group_SLF1_R_FD<- SLF_data %>%
#     group_by(Timepoint) %>%
#     anova_test(dv=mn_FD_SLF1_R,wid=Individual_number,between=Group) %>%
#     adjust_pvalue(method="fdr")
# posthoc_ME_Group_SLF1_R_FD
# #Pairwise comparison between groups levels
# posthoc_pairwise_Group_SLF1_R_FD <- SLF_data %>%
#     group_by(Timepoint) %>%
#     pairwise_t_test(mn_FD_SLF1_R ~ Group, p.adjust.method = "fdr")
# posthoc_pairwise_Group_SLF1_R_FD
# #Simple main effect w/ ineraction - for Timepoint
# posthoc_ME_Timepoint_SLF1_R_FD <- SLF_data %>%
#     group_by(Group) %>%
#     anova_test(dv=mn_FD_SLF1_R,wid=Individual_number,within=Timepoint,effect.size = "pes") %>%
#     get_anova_table() %>%
#     adjust_pvalue(method="fdr")
# posthoc_ME_Timepoint_SLF1_R_FD
# #Pairwise comparison between groups levels
# posthoc_pairwise_Timepoint_SLF1_R_FD <- SLF_data %>%
#     group_by(Group) %>%
#     pairwise_t_test(
#     mn_FD_SLF1_R ~ Timepoint, paired = TRUE, 
#     p.adjust.method = "fdr") %>%
#     select(-df, -statistic, -p) # Remove details
# posthoc_pairwise_Timepoint_SLF1_R_FD
# #effect size for interaction (SCD Group and Time Point)
# SLF_data_SCD_SLF1_R_FD_long <- subset(SLF_data, DPRC_neuropsych_data$Group == 2)
# SLF_data_SCD_SLF1_R_FD_long$Group <- droplevels(SLF_data_SCD_SLF1_R_FD_long$Group)
# SLF_data_SCD_SLF1_R_FD_long%>%cohens_d(mn_FD_SLF1_R~Timepoint, paired=TRUE)

# #Interaction - Tukey correction
# TukeyHSD(aov(mn_FD_SLF1_R~Group+Timepoint+Group:Timepoint,data=SLF_data))
# #another option:
# aov(mn_FD_SLF1_R ~ Group*Timepoint, data = SLF_data) %>% tukey_hsd()

# #CI for Interaction
# #Right SLF 1 FD
# post_hoc_aov_SLF1_R_FD_mod <- lme(mn_FD_SLF1_R ~ Group*Timepoint, random = ~1 | ParticipantID/Timepoint, data=SLF_data)
# summary(glht(post_hoc_aov_SLF1_R_FD_mod, linfct=mcp(Timepoint="Tukey")))
# nlme::intervals(post_hoc_aov_SLF1_R_FD_mod, level=0.95)

#plot 95% Confidence Interval (FD) - differences only found in Group
SLF_FD_95CI_data <- data.frame(SLF_group_number = c('1','1','1','2','2','2','3','3','3','3','4','4','4','4','5','5','5','5','6','6','6','6','7','7','7','7'),
                               SLF_type = c('Whole_SLF','Whole_SLF','Whole_SLF','Left_SLF','Left_SLF','Left_SLF','Right_SLF','Right_SLF','Right_SLF','Right_SLF','Left_SLF2','Left_SLF2','Left_SLF2','Left_SLF2','Right_SLF1','Right_SLF1','Right_SLF1','Right_SLF1','Right_SLF2','Right_SLF2','Right_SLF2','Right_SLF2','Right_SLF3','Right_SLF3','Right_SLF3','Right_SLF3'),
                               Group_contrast = c('C>SCD','C>aMCI','C>AD','C>SCD','C>aMCI','C>AD','C>SCD','C>aMCI','C>mMCI','C>AD','C>AD','SCD>AD','aMCI>AD','mMCI>AD','C>SCD','C>aMCI','C>mMCI','C>AD','C>SCD','C>aMCI','C>AD','mMCI>AD','C>SCD','C>aMCI','C>mMCI','C>AD'),
                               estimate_diff = c(-0.0145,-0.0130,-0.0251,-0.0111,-0.0117,-0.0230,-0.0180,-0.0144,-0.0129,-0.0272,-0.0335,-0.0254,-0.0205,-0.0224,-0.0180,-0.0141,-0.0157,-0.0305,-0.0209,-0.0171,-0.0324,-0.0207,-0.0196,-0.0160,-0.0138,-0.0234)*-1,
                               lower = c(-0.0248,-0.0240,-0.0395,-0.0217,-0.0229,-0.0379,-0.0290,-0.0260,-0.0255,-0.0425,-0.0526,-0.0430,-0.0387,-0.0416,-0.0304,-0.0273,-0.0300,-0.0478,-0.0346,-0.0316,-0.0515,-0.0399,-0.0313,-0.0283,-0.0272,-0.0396)*-1, 
                               upper = c(-0.00406, -0.00201,-0.0106,-0.000400,-0.000377,-0.00817, -0.00707,-0.00284,-0.000361,-0.0120,-0.0145,-0.00786,-0.00234,-0.00318,-0.00555,-0.000948,-0.00142,-0.0132,-0.00717,-0.00263,-0.0134,-0.00143,-0.00799,-0.00370,-0.000422,-0.00722)*-1)  
#plot data
ggplot(SLF_FD_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Whole SLF", "2" = "Left SLF", "3" = "Right SLF", "4" = "Left SLF 2", "5" = "Right SLF 1", "6" = "Right SLF 2", "7" = "Right SLF 3"))+
    theme_classic() 

#Run ANCOVA With covariates (age & sex)
#whole SLF
#does not work w/ covariate, Age
#SLF_data$Age <- as.factor(SLF_data$Age)
#aov_SLF_FD_2covar_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF, wid=Individual_number, between=Group, within=Timepoint, covariate=c(Sex,Age), effect.size ="pes")
#get_anova_table(aov_SLF_FD_2covar_mod) 
aov_SLF_FD_2covar_mod<- aov(mn_FD_SLF ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
#aov_SLF_FD_2covar_mod<- anova_test(mn_FD_SLF ~ Age + Sex + Group*Timepoint, data=SLF_data, wid=Individual_number, effect.size = "pes") #another way
summary(aov_SLF_FD_2covar_mod)
aov_SLF_FC_2covar_mod<- aov(mn_FC_SLF ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF_FC_2covar_mod)
aov_SLF_FDC_2covar_mod<- aov(mn_FDC_SLF ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF_FDC_2covar_mod)
#Left SLF
aov_SLF_L_FD_2covar_mod<- aov(mn_FD_SLF_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF_L_FD_2covar_mod)
aov_SLF_L_FC_2covar_mod<- aov(mn_FC_SLF_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF_L_FC_2covar_mod)
aov_SLF_L_FDC_2covar_mod<- aov(mn_FDC_SLF_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF_L_FDC_2covar_mod)
#Right SLF
aov_SLF_R_FD_2covar_mod<- aov(mn_FD_SLF_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF_R_FD_2covar_mod)
aov_SLF_R_FC_2covar_mod<- aov(mn_FC_SLF_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF_R_FC_2covar_mod)
aov_SLF_R_FDC_2covar_mod<- aov(mn_FDC_SLF_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF_R_FDC_2covar_mod)
#Left SLF 1
aov_SLF1_L_FD_2covar_mod<- aov(mn_FD_SLF1_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF1_L_FD_2covar_mod)
aov_SLF1_L_FC_2covar_mod<- aov(mn_FC_SLF1_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF1_L_FC_2covar_mod)
aov_SLF1_L_FDC_2covar_mod<- aov(mn_FDC_SLF1_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF1_L_FDC_2covar_mod)
#Left SLF 2
aov_SLF2_L_FD_2covar_mod<- aov(mn_FD_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_L_FD_2covar_mod) 
aov_SLF2_L_FC_2covar_mod<- aov(mn_FC_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_L_FC_2covar_mod)
aov_SLF2_L_FDC_2covar_mod<- aov(mn_FDC_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_L_FDC_2covar_mod) 
#Left SLF 3
aov_SLF3_L_FD_2covar_mod<- aov(mn_FD_SLF3_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF3_L_FD_2covar_mod)
aov_SLF3_L_FC_2covar_mod<- aov(mn_FC_SLF3_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF3_L_FC_2covar_mod)
aov_SLF3_L_FDC_2covar_mod<- aov(mn_FDC_SLF3_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF3_L_FDC_2covar_mod)
#Right SLF 1
aov_SLF1_R_FD_2covar_mod<- aov(mn_FD_SLF1_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF1_R_FD_2covar_mod)
aov_SLF1_R_FC_2covar_mod<- aov(mn_FC_SLF1_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF1_R_FC_2covar_mod)
aov_SLF1_R_FDC_2covar_mod<- aov(mn_FDC_SLF1_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF1_R_FDC_2covar_mod)
#Right SLF 2
aov_SLF2_R_FD_2covar_mod<- aov(mn_FD_SLF2_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_R_FD_2covar_mod) 
aov_SLF2_R_FC_2covar_mod<- aov(mn_FC_SLF2_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_R_FC_2covar_mod)
aov_SLF2_R_FDC_2covar_mod<- aov(mn_FDC_SLF2_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_R_FDC_2covar_mod) 
#Right SLF 3
aov_SLF3_R_FD_2covar_mod<- aov(mn_FD_SLF3_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF3_R_FD_2covar_mod)
aov_SLF3_R_FC_2covar_mod<- aov(mn_FC_SLF3_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF3_R_FC_2covar_mod)
aov_SLF3_R_FDC_2covar_mod<- aov(mn_FDC_SLF3_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF3_R_FDC_2covar_mod)

#effect size w/ covariates
sjstats::eta_sq(aov_SLF_FD_2covar_mod)
sjstats::eta_sq(aov_SLF_FC_2covar_mod)
sjstats::eta_sq(aov_SLF_FDC_2covar_mod)
sjstats::eta_sq(aov_SLF_L_FD_2covar_mod)
sjstats::eta_sq(aov_SLF_L_FC_2covar_mod)
sjstats::eta_sq(aov_SLF_L_FDC_2covar_mod)
sjstats::eta_sq(aov_SLF_R_FD_2covar_mod)
sjstats::eta_sq(aov_SLF_R_FC_2covar_mod)
sjstats::eta_sq(aov_SLF_R_FDC_2covar_mod)
sjstats::eta_sq(aov_SLF1_L_FD_2covar_mod)
sjstats::eta_sq(aov_SLF1_L_FC_2covar_mod)
sjstats::eta_sq(aov_SLF1_L_FDC_2covar_mod)
sjstats::eta_sq(aov_SLF2_L_FD_2covar_mod)
sjstats::eta_sq(aov_SLF2_L_FC_2covar_mod)
sjstats::eta_sq(aov_SLF2_L_FDC_2covar_mod)
sjstats::eta_sq(aov_SLF3_L_FD_2covar_mod)
sjstats::eta_sq(aov_SLF3_L_FC_2covar_mod)
sjstats::eta_sq(aov_SLF3_L_FDC_2covar_mod)
sjstats::eta_sq(aov_SLF1_R_FD_2covar_mod)
sjstats::eta_sq(aov_SLF1_R_FC_2covar_mod)
sjstats::eta_sq(aov_SLF1_R_FDC_2covar_mod)
sjstats::eta_sq(aov_SLF2_R_FD_2covar_mod)
sjstats::eta_sq(aov_SLF2_R_FC_2covar_mod)
sjstats::eta_sq(aov_SLF2_R_FDC_2covar_mod)
sjstats::eta_sq(aov_SLF3_R_FD_2covar_mod)
sjstats::eta_sq(aov_SLF3_R_FC_2covar_mod)
sjstats::eta_sq(aov_SLF3_R_FDC_2covar_mod)

#follow up tests with sig. main effects and interactions
#Whole SLF FD
post_hoc_aov_SLF_FD_2covar_mod <- lme(mn_FD_SLF ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
#post_hoc_aov_SLF_FD_2covar_mod <- lm(mn_FD_SLF ~ Group*Timepoint+Age+Sex, data=SLF_data)
summary(glht(post_hoc_aov_SLF_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#confint(post_hoc_aov_SLF_FD_2covar_mod)
#nlme::intervals(post_hoc_aov_SLF_FD_2covar_mod, level=0.95)
#Whole SLF FC
post_hoc_aov_SLF_FC_2covar_mod <- lme(mn_FC_SLF ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF_FC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Left SLF FD
post_hoc_aov_SLF_L_FD_2covar_mod <- lme(mn_FD_SLF_L ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF_L_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#Left SLF FC
post_hoc_aov_SLF_L_FC_2covar_mod <- lme(mn_FC_SLF_L ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF_L_FC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Right SLF FD
post_hoc_aov_SLF_R_FD_2covar_mod <- lme(mn_FD_SLF_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#Right SLF FC
post_hoc_aov_SLF_R_FC_2covar_mod <- lme(mn_FC_SLF_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF_R_FC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Left SLF1 FC
post_hoc_aov_SLF1_L_FC_2covar_mod <- lme(mn_FC_SLF1_L ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_L_FC_2covar_mod, linfct=mcp(Group="Tukey")))
summary(glht(post_hoc_aov_SLF1_L_FC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Left SLF1 FDC
post_hoc_aov_SLF1_L_FDC_2covar_mod <- lme(mn_FDC_SLF1_L ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_L_FDC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Left SLF 2 FD
post_hoc_aov_SLF2_L_FD_2covar_mod <- lme(mn_FD_SLF2_L ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_L_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#Left SLF 2 FC
post_hoc_aov_SLF2_L_FC_2covar_mod <- lme(mn_FC_SLF2_L ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_L_FC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Right SLF1 FD
post_hoc_aov_SLF1_R_FD_2covar_mod <- lme(mn_FD_SLF1_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#Right SLF1 FC
post_hoc_aov_SLF1_R_FC_2covar_mod <- lme(mn_FC_SLF1_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_R_FC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Right SLF1 FDC
post_hoc_aov_SLF1_R_FDC_2covar_mod <- lme(mn_FDC_SLF1_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_R_FDC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Right SLF2 FD
post_hoc_aov_SLF2_R_FD_2covar_mod <- lme(mn_FD_SLF2_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#Right SLF2 FC
post_hoc_aov_SLF2_R_FC_2covar_mod <- lme(mn_FC_SLF2_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_R_FC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Right SLF3 FD
post_hoc_aov_SLF3_R_FD_2covar_mod <- lme(mn_FD_SLF3_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF3_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#Right SLF3 FDC
post_hoc_aov_SLF3_R_FDC_2covar_mod <- lme(mn_FDC_SLF3_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF3_R_FDC_2covar_mod, linfct=mcp(Group="Tukey")))

# #follow up for sig. interaction 
# #Right SLF 1
# #interaction post hoc follow-up
# #view interaction plot
# SLF_data %>%
#     group_by(Group,Timepoint) %>%
#     summarise(s_mean=mean(mn_FD_SLF1_R)) %>%
#     ggplot(aes(y=s_mean,x=Group,colour=Timepoint,group=Timepoint))+
#     geom_point()+geom_line()+
#     ylab("SLF1 Right Fibre Density (FD)") +
#     theme_classic()
# #view interaction on the variable for 'Group' - no post hoc correction applied? 
# SLF_data %>% filter(Group=="1") %>%
#     aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
#     summary()
# SLF_data %>% filter(Group=="2") %>% #significant interaction for SCD group
#     aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
#     summary()
# SLF_data %>% filter(Group=="3") %>% 
#     aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
#     summary()
# SLF_data %>% filter(Group=="4") %>% 
#     aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
#     summary()
# SLF_data %>% filter(Group=="5") %>% 
#     aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
#     summary()
# #view interaction on the variable for 'Timepoint' - no post hoc correction applied? 
# SLF_data %>% filter(Timepoint=="F0") %>%
#     aov(mn_FD_SLF1_R~Group+Sex+Age,data=.)%>%
#     summary()
# SLF_data %>% filter(Timepoint=="F2") %>%
#     aov(mn_FD_SLF1_R~Group+Sex+Age,data=.)%>%
#     summary()
# 
# #Interaction - #Tukey test
# TukeyHSD(aov(mn_FD_SLF1_R~Group*Timepoint+as.factor(Age)+Sex, data = SLF_data)) 
# #another option
# aov(mn_FD_SLF1_R ~ Group*Timepoint+Sex+as.factor(Age), data = SLF_data)  %>% tukey_hsd()
# 
# #CI for Interaction w/ covariates
# #Right SLF 1 FD
# post_hoc_aov_SLF1_R_FD_2covar_mod <- lme(mn_FD_SLF1_R ~ Group*Timepoint + Age + Sex, random = ~1 | ParticipantID/Timepoint, data=SLF_data)
# summary(glht(post_hoc_aov_SLF1_R_FD_2covar_mod, linfct=mcp(Timepoint="Tukey")))
# nlme::intervals(post_hoc_aov_SLF1_R_FD_2covar_mod, level=0.95) #doesn't work w/ covariates in the model.
# 
# #plot 95% Confidence Interval w/ covariates (FD)
# SLF_FD_95CI_2covar_data <- data.frame(SLF_group_number = c('1'),
#                                SLF_type = c('Left_SLF2'),
#                                Group_contrast = c('F0vF2'),
#                                estimate_diff = c(-0.0058891960),
#                                lower = c(-0.0234024688), 
#                                upper = c(0.0116240768))  
# #plot data
# ggplot(SLF_FD_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
#     geom_point()+
#     geom_errorbar(aes(ymin=lower, ymax=upper))+
#     xlab("SLF Tract") + 
#     ylab("95% Confidence Interval")+
#     scale_x_discrete(labels = c("1" = "Left SLF 2"))+
#     theme_classic() 
# 
# #plot 95% Confidence Interval w/ covariates (FDC)
# SLF_FDC_95CI_2covar_data <- data.frame(SLF_group_number = c('1'),
#                                 SLF_type = c('Left_SLF2'),
#                                 Group_contrast = c('F0vF2'),
#                                 estimate_diff = c(-0.0077689181),
#                                 lower = c(-0.040474001), 
#                                 upper = c(0.0249361650))  
# #plot data
# ggplot(SLF_FDC_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
#     geom_point()+
#     geom_errorbar(aes(ymin=lower, ymax=upper))+
#     xlab("SLF Tract") + 
#     ylab("95% Confidence Interval")+
#     scale_x_discrete(labels = c("1" = "Left SLF 2"))+
#     theme_classic() 

#Combine all sig. plots on one plot for FD
#put data into long format
SLF_data_FD <- dplyr::select(SLF_data, 
                             ParticipantID,
                             Group,
                             Timepoint,
                             mn_FD_SLF,
                             mn_FD_SLF_L,
                             mn_FD_SLF_R,
                             mn_FD_SLF1_L,
                             mn_FD_SLF2_L,
                             mn_FD_SLF3_L,
                             mn_FD_SLF1_R,
                             mn_FD_SLF2_R,
                             mn_FD_SLF3_R)
SLF_data_FD_long <- gather(SLF_data_FD, 
                           "SLF_type",
                           "FD_metric",
                           mn_FD_SLF,
                           mn_FD_SLF_L,
                           mn_FD_SLF_R,
                           mn_FD_SLF1_L,
                           mn_FD_SLF2_L,
                           mn_FD_SLF3_L,
                           mn_FD_SLF1_R,
                           mn_FD_SLF2_R,
                           mn_FD_SLF3_R)
#All tracts SLF FD (raincloud plot)
ggplot(SLF_data_FD_long, aes(x = SLF_type, y = FD_metric, fill = Timepoint)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = FD_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF tract") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("mn_FD_SLF" = "Whole SLF", "mn_FD_SLF_L" = "Left SLF", "mn_FD_SLF_R" = "Right SLF", "mn_FD_SLF1_L" = "Left SLF1","mn_FD_SLF2_L" = "Left SLF2","mn_FD_SLF3_L" = "Left SLF3","mn_FD_SLF1_R" = "Right SLF1","mn_FD_SLF2_R" = "Right SLF2","mn_FD_SLF3_R" = "Right SLF3")) + 
    scale_colour_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
    theme_classic()

#All tracts SLF FD (raincloud plot - highlight the group differences only)
ggplot(SLF_data_FD_long, aes(x = SLF_type, y = FD_metric, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = FD_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    guides(colour=FALSE)+
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF tract") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("mn_FD_SLF" = "Whole SLF", "mn_FD_SLF_L" = "Left SLF", "mn_FD_SLF_R" = "Right SLF", "mn_FD_SLF1_L" = "Left SLF1","mn_FD_SLF2_L" = "Left SLF2","mn_FD_SLF3_L" = "Left SLF3","mn_FD_SLF1_R" = "Right SLF1","mn_FD_SLF2_R" = "Right SLF2","mn_FD_SLF3_R" = "Right SLF3")) + 
    scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
    theme_classic()

#Combine all sig. plots on one plot for FC
#put data into long format
SLF_data_FC <- dplyr::select(SLF_data, 
                             ParticipantID,
                             Group,
                             Timepoint,
                             mn_FC_SLF,
                             mn_FC_SLF_L,
                             mn_FC_SLF_R,
                             mn_FC_SLF1_L,
                             mn_FC_SLF2_L,
                             mn_FC_SLF3_L,
                             mn_FC_SLF1_R,
                             mn_FC_SLF2_R,
                             mn_FC_SLF3_R)
SLF_data_FC_long <- gather(SLF_data_FC, 
                           "SLF_type",
                           "FC_metric",
                           mn_FC_SLF,
                           mn_FC_SLF_L,
                           mn_FC_SLF_R,
                           mn_FC_SLF1_L,
                           mn_FC_SLF2_L,
                           mn_FC_SLF3_L,
                           mn_FC_SLF1_R,
                           mn_FC_SLF2_R,
                           mn_FC_SLF3_R)
#All tracts SLF FC (raincloud plot)
ggplot(SLF_data_FC_long, aes(x = SLF_type, y = FC_metric, fill = Timepoint)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = FC_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF tract") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("mn_FC_SLF" = "Whole SLF", "mn_FC_SLF_L" = "Left SLF", "mn_FC_SLF_R" = "Right SLF", "mn_FC_SLF1_L" = "Left SLF1","mn_FC_SLF2_L" = "Left SLF2","mn_FC_SLF3_L" = "Left SLF3","mn_FC_SLF1_R" = "Right SLF1","mn_FC_SLF2_R" = "Right SLF2","mn_FC_SLF3_R" = "Right SLF3")) + 
    scale_colour_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
    theme_classic()

# #plot individual changes from F0 to F2 for FC metric per each SLF tract
# ggplot(SLF_data_FC_long, aes(x = SLF_type, y = FC_metric, fill=Timepoint)) + 
#     #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_boxplot(aes(x = SLF_type, y = FC_metric, colour=Timepoint), width = 0.1, fill = "white", outlier.size = 1) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
#     geom_point(aes(colour = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_line(aes(group = ParticipantID, colour=ParticipantID)) +
#     theme(legend.position = "none") + 
#     #guides(fill = guide_legend(override.aes = list(shape = NA)))+
#     xlab("SLF tract") + 
#     ylab("Fibre Cross-section (FC)") +
#     scale_x_discrete(labels = c("mn_FC_SLF" = "Whole SLF", "mn_FC_SLF_L" = "Left SLF", "mn_FC_SLF_R" = "Right SLF", "mn_FC_SLF1_L" = "Left SLF1","mn_FC_SLF2_L" = "Left SLF2","mn_FC_SLF3_L" = "Left SLF3","mn_FC_SLF1_R" = "Right SLF1","mn_FC_SLF2_R" = "Right SLF2","mn_FC_SLF3_R" = "Right SLF3")) + 
#     #scale_colour_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
#     theme_classic() +
# 
# #for all SLF tracts - look at FC    
# #make a new dataset for this
# SLF_FC_all <- SLF_data_FC_long %>% group_by(ParticipantID) %>% summarise_at(vars(FC_metric),funs(mean(.,na.rm=TRUE)))
# SLF_FC_all['Timepoint'] <- SLF_data$Timepoint 
# SLF_FC_all['Group'] <- SLF_data$Group 
# 
# ggplot(SLF_FC_all, aes(x = Timepoint, y = FC_metric, colour=Timepoint, group=ParticipantID)) + 
#     geom_point(aes(colour=Timepoint),position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_line()+
#     theme(legend.position = "none") + 
#     geom_boxplot(aes(colour=Timepoint), width = 0.1, fill = "white", outlier.size = 1) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
#     #geom_point(aes(colour = ParticipantID)) +
#     #theme(legend.position = "none")
#     #guides(fill = guide_legend(override.aes = list(shape = NA)))+
#     xlab("SLF tract") + 
#     ylab("Fibre Cross-section (FC)") +
#     scale_x_discrete(labels = c("mn_FC_SLF" = "Whole SLF", "mn_FC_SLF_L" = "Left SLF", "mn_FC_SLF_R" = "Right SLF", "mn_FC_SLF1_L" = "Left SLF1","mn_FC_SLF2_L" = "Left SLF2","mn_FC_SLF3_L" = "Left SLF3","mn_FC_SLF1_R" = "Right SLF1","mn_FC_SLF2_R" = "Right SLF2","mn_FC_SLF3_R" = "Right SLF3")) + 
#     #scale_colour_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
#     theme_classic()    
    

#Combine all sig. plots on one plot for FDC
#put data into long format
SLF_data_FDC <- dplyr::select(SLF_data, 
                             ParticipantID,
                             Group,
                             Timepoint,
                             mn_FDC_SLF,
                             mn_FDC_SLF_L,
                             mn_FDC_SLF_R,
                             mn_FDC_SLF1_L,
                             mn_FDC_SLF2_L,
                             mn_FDC_SLF3_L,
                             mn_FDC_SLF1_R,
                             mn_FDC_SLF2_R,
                             mn_FDC_SLF3_R)
SLF_data_FDC_long <- gather(SLF_data_FDC, 
                           "SLF_type",
                           "FDC_metric",
                           mn_FDC_SLF,
                           mn_FDC_SLF_L,
                           mn_FDC_SLF_R,
                           mn_FDC_SLF1_L,
                           mn_FDC_SLF2_L,
                           mn_FDC_SLF3_L,
                           mn_FDC_SLF1_R,
                           mn_FDC_SLF2_R,
                           mn_FDC_SLF3_R)
#All tracts SLF FDC (raincloud plot)
ggplot(SLF_data_FDC_long, aes(x = SLF_type, y = FDC_metric, fill = Timepoint)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = FDC_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF tract") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("mn_FDC_SLF" = "Whole SLF", "mn_FDC_SLF_L" = "Left SLF", "mn_FDC_SLF_R" = "Right SLF", "mn_FDC_SLF1_L" = "Left SLF1","mn_FDC_SLF2_L" = "Left SLF2","mn_FDC_SLF3_L" = "Left SLF3","mn_FDC_SLF1_R" = "Right SLF1","mn_FDC_SLF2_R" = "Right SLF2","mn_FDC_SLF3_R" = "Right SLF3")) + 
    scale_colour_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
    theme_classic()


#Run correlation differences and linear regression comparison between SLF metrics and neuropsych test scores
#prep dataframe
SLF_data_F0 <- SLF_data[1:124,]
SLF_data_F2 <- SLF_data[125:248,]
SLF_data_subtracted <- SLF_data_F2[, c(1:189, 197:220)] - SLF_data_F0[, c(1:189, 197:220)]
#add in Participant ID 
SLF_data_subtracted$ParticipantID <- SLF_data_F0$ParticipantID
#add in Group classification column for each participant - data from another worksheet
SLF_data_subtracted$Group <-SLF_data_F0$Group
#add in Timepoint (F0 vs. F2)
SLF_data_subtracted$Timepoint <- SLF_data_F0$Timepoint
#add in covariates (clinical site) from the covariates_data dataframe to the SLF dataframe
SLF_data_subtracted$ClinSite_name <- SLF_data_F0$Clinical_site 
SLF_data_subtracted$Age <-SLF_data_F0$Age 
SLF_data_subtracted$Sex <- SLF_data_F0$Sex
SLF_data_subtracted$Sex_binary <- SLF_data_F0$Sex_binary

####------------------------ Correlation -----------------------------------####

#Run correlation test on the difference scores (F2 - F0)  of SLF metrics vs. neurospsych tests
#Whole SLF-
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$WordReading.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF-
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF_L, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF-
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$HayBTime2.Raw) #sig.
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF_R, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF1-
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$HayBTime2.Raw) #sig.
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF1_L, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF2-
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF2_L, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF3-
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$WordReading.Raw) #sig.
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF3_L, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF1-
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF1_R, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF2-
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF2_R, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF3-
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$TrailsA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$HayBTime2.Raw) #sig. 
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FD_SLF3_R, SLF_data_subtracted$HayBCatB.Raw) 

#for FC
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF-
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF_L, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF-
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF_R, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF1-
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF1_L, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF2-
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF2_L, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF3-
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF3_L, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF1-
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF1_R, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF2-
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF2_R, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF3-
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$TrailsA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$HayBTime2.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FC_SLF3_R, SLF_data_subtracted$HayBCatB.Raw) 

#for FDC - Whole FDC
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$HayBTime2.Raw) #sig.
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF-
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$HayBTime2.Raw) #sig.
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF_L, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF-
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$HayBTime2.Raw) #sig.
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF_R, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF1-
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$HayBTime2.Raw) #sig.
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF1_L, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF2-
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF2_L, SLF_data_subtracted$HayBCatB.Raw) 
#Left SLF3-
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF3_L, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF1-
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF1_R, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF2-
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$TrailsA.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$HayBTime2.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF2_R, SLF_data_subtracted$HayBCatB.Raw) 
#Right SLF3-
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$TrailsA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$TrailsB.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$ColorNaming.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$WordReading.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$Inhibition.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$LetFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$CatFluency.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$Switching.Raw)
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$HayBTime1.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$HayBTime2.Raw) #sig. 
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$HayBCatA.Raw) 
cor.test(SLF_data_subtracted$mn_FDC_SLF3_R, SLF_data_subtracted$HayBCatB.Raw) 



#Run partial correlation (controlling for covariates, age and sex) - note that
#missing values are not allowed.
#for TrailsA:
SLF_data_subtracted_TrailsA <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "TrailsA.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                               "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                               "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                               "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                               "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                               "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_TrailsA <- SLF_data_subtracted_TrailsA %>% drop_na(TrailsA.Raw)
SLF_data_subtracted_TrailsA$Sex_binary <- as.numeric(SLF_data_subtracted_TrailsA$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_TrailsA$mn_FD_SLF, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FD_SLF_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FD_SLF_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FD_SLF1_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FD_SLF2_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FD_SLF3_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FD_SLF1_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FD_SLF2_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FD_SLF3_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
#for FC
pcor.test(SLF_data_subtracted_TrailsA$mn_FC_SLF, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FC_SLF_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FC_SLF_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FC_SLF1_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FC_SLF2_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FC_SLF3_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FC_SLF1_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FC_SLF2_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FC_SLF3_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_TrailsA$mn_FDC_SLF, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FDC_SLF_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FDC_SLF_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FDC_SLF1_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FDC_SLF2_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FDC_SLF3_L, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FDC_SLF1_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FDC_SLF2_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsA$mn_FDC_SLF3_R, SLF_data_subtracted_TrailsA$TrailsA.Raw, SLF_data_subtracted_TrailsA[,c("Age", "Sex_binary")])
#for TrailsB:
SLF_data_subtracted_TrailsB <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "TrailsB.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                               "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                               "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                               "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                               "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                               "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_TrailsB <- SLF_data_subtracted_TrailsB %>% drop_na(TrailsB.Raw)
SLF_data_subtracted_TrailsB$Sex_binary <- as.numeric(SLF_data_subtracted_TrailsB$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_TrailsB$mn_FD_SLF, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_TrailsB$mn_FD_SLF_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_TrailsB$mn_FD_SLF_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_TrailsB$mn_FD_SLF1_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FD_SLF2_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_TrailsB$mn_FD_SLF3_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_TrailsB$mn_FD_SLF1_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_TrailsB$mn_FD_SLF2_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_TrailsB$mn_FD_SLF3_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_TrailsB$mn_FC_SLF, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FC_SLF_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FC_SLF_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FC_SLF1_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FC_SLF2_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FC_SLF3_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FC_SLF1_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FC_SLF2_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FC_SLF3_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_TrailsB$mn_FDC_SLF, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FDC_SLF_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FDC_SLF_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FDC_SLF1_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FDC_SLF2_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FDC_SLF3_L, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FDC_SLF1_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FDC_SLF2_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_TrailsB$mn_FDC_SLF3_R, SLF_data_subtracted_TrailsB$TrailsB.Raw, SLF_data_subtracted_TrailsB[,c("Age", "Sex_binary")])
#for ColorNaming:
SLF_data_subtracted_ColorNaming <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "ColorNaming.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                   "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                   "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                   "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                   "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                   "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_ColorNaming <- SLF_data_subtracted_ColorNaming %>% drop_na(ColorNaming.Raw)
SLF_data_subtracted_ColorNaming$Sex_binary <- as.numeric(SLF_data_subtracted_ColorNaming$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_ColorNaming$mn_FD_SLF, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_ColorNaming$mn_FD_SLF_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_ColorNaming$mn_FD_SLF_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_ColorNaming$mn_FD_SLF1_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FD_SLF2_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FD_SLF3_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_ColorNaming$mn_FD_SLF1_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FD_SLF2_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_ColorNaming$mn_FD_SLF3_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_ColorNaming$mn_FC_SLF, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FC_SLF_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FC_SLF_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FC_SLF1_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FC_SLF2_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FC_SLF3_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FC_SLF1_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FC_SLF2_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FC_SLF3_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_ColorNaming$mn_FDC_SLF, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FDC_SLF_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FDC_SLF_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FDC_SLF1_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FDC_SLF2_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_ColorNaming$mn_FDC_SLF3_L, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FDC_SLF1_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FDC_SLF2_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_ColorNaming$mn_FDC_SLF3_R, SLF_data_subtracted_ColorNaming$ColorNaming.Raw, SLF_data_subtracted_ColorNaming[,c("Age", "Sex_binary")])
#for WordReading:
SLF_data_subtracted_WordReading <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "WordReading.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                   "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                   "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                   "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                   "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                   "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_WordReading <- SLF_data_subtracted_WordReading %>% drop_na(WordReading.Raw)
SLF_data_subtracted_WordReading$Sex_binary <- as.numeric(SLF_data_subtracted_WordReading$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_WordReading$mn_FD_SLF, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_WordReading$mn_FD_SLF_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_WordReading$mn_FD_SLF_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_WordReading$mn_FD_SLF1_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FD_SLF2_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FD_SLF3_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_subtracted_WordReading$mn_FD_SLF1_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FD_SLF2_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_WordReading$mn_FD_SLF3_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_WordReading$mn_FC_SLF, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FC_SLF_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FC_SLF_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FC_SLF1_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FC_SLF2_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FC_SLF3_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FC_SLF1_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FC_SLF2_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FC_SLF3_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_WordReading$mn_FDC_SLF, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FDC_SLF_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FDC_SLF_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FDC_SLF1_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FDC_SLF2_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FDC_SLF3_L, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FDC_SLF1_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FDC_SLF2_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_WordReading$mn_FDC_SLF3_R, SLF_data_subtracted_WordReading$WordReading.Raw, SLF_data_subtracted_WordReading[,c("Age", "Sex_binary")])
#for Inhibition:
SLF_data_subtracted_Inhibition <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "Inhibition.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                  "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                  "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                  "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                  "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                  "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_Inhibition <- SLF_data_subtracted_Inhibition %>% drop_na(Inhibition.Raw)
SLF_data_subtracted_Inhibition$Sex_binary <- as.numeric(SLF_data_subtracted_Inhibition$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_Inhibition$mn_FD_SLF, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FD_SLF_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FD_SLF_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FD_SLF1_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FD_SLF2_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FD_SLF3_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FD_SLF1_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FD_SLF2_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FD_SLF3_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_Inhibition$mn_FC_SLF, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FC_SLF_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FC_SLF_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FC_SLF1_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FC_SLF2_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FC_SLF3_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FC_SLF1_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FC_SLF2_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FC_SLF3_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_Inhibition$mn_FDC_SLF, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FDC_SLF_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FDC_SLF_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FDC_SLF1_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FDC_SLF2_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FDC_SLF3_L, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FDC_SLF1_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Inhibition$mn_FDC_SLF2_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Inhibition$mn_FDC_SLF3_R, SLF_data_subtracted_Inhibition$Inhibition.Raw, SLF_data_subtracted_Inhibition[,c("Age", "Sex_binary")]) 
#for LetFluency:
SLF_data_subtracted_LetFluency <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "LetFluency.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                  "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                  "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                  "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                  "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                  "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_LetFluency <- SLF_data_subtracted_LetFluency %>% drop_na(LetFluency.Raw)
SLF_data_subtracted_LetFluency$Sex_binary <- as.numeric(SLF_data_subtracted_LetFluency$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_LetFluency$mn_FD_SLF, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_LetFluency$mn_FD_SLF_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_LetFluency$mn_FD_SLF_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_LetFluency$mn_FD_SLF1_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FD_SLF2_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_LetFluency$mn_FD_SLF3_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_LetFluency$mn_FD_SLF1_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FD_SLF2_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_LetFluency$mn_FD_SLF3_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_LetFluency$mn_FC_SLF, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FC_SLF_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FC_SLF_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FC_SLF1_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FC_SLF2_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FC_SLF3_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FC_SLF1_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FC_SLF2_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FC_SLF3_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_LetFluency$mn_FDC_SLF, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FDC_SLF_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FDC_SLF_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FDC_SLF1_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FDC_SLF2_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FDC_SLF3_L, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FDC_SLF1_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FDC_SLF2_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_LetFluency$mn_FDC_SLF3_R, SLF_data_subtracted_LetFluency$LetFluency.Raw, SLF_data_subtracted_LetFluency[,c("Age", "Sex_binary")])
#for CatFluency:
SLF_data_subtracted_CatFluency <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "CatFluency.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                  "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                  "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                  "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                  "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                  "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_CatFluency <- SLF_data_subtracted_CatFluency %>% drop_na(CatFluency.Raw)
SLF_data_subtracted_CatFluency$Sex_binary <- as.numeric(SLF_data_subtracted_CatFluency$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_CatFluency$mn_FD_SLF, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_CatFluency$mn_FD_SLF_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_CatFluency$mn_FD_SLF_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_CatFluency$mn_FD_SLF1_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FD_SLF2_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_CatFluency$mn_FD_SLF3_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_CatFluency$mn_FD_SLF1_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FD_SLF2_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_CatFluency$mn_FD_SLF3_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_CatFluency$mn_FC_SLF, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FC_SLF_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FC_SLF_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FC_SLF1_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FC_SLF2_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FC_SLF3_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FC_SLF1_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FC_SLF2_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FC_SLF3_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_CatFluency$mn_FDC_SLF, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FDC_SLF_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FDC_SLF_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FDC_SLF1_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FDC_SLF2_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FDC_SLF3_L, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FDC_SLF1_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FDC_SLF2_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_CatFluency$mn_FDC_SLF3_R, SLF_data_subtracted_CatFluency$CatFluency.Raw, SLF_data_subtracted_CatFluency[,c("Age", "Sex_binary")])
#for Switching:
SLF_data_subtracted_Switching <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "Switching.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                 "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                 "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                 "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                 "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                 "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_Switching <- SLF_data_subtracted_Switching %>% drop_na(Switching.Raw)
SLF_data_subtracted_Switching$Sex_binary <- as.numeric(SLF_data_subtracted_Switching$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_Switching$mn_FD_SLF, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Switching$mn_FD_SLF_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Switching$mn_FD_SLF_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Switching$mn_FD_SLF1_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FD_SLF2_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Switching$mn_FD_SLF3_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Switching$mn_FD_SLF1_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FD_SLF2_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Switching$mn_FD_SLF3_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_Switching$mn_FC_SLF, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FC_SLF_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FC_SLF_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FC_SLF1_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FC_SLF2_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FC_SLF3_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FC_SLF1_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FC_SLF2_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FC_SLF3_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_Switching$mn_FDC_SLF, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FDC_SLF_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FDC_SLF_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FDC_SLF1_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FDC_SLF2_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_Switching$mn_FDC_SLF3_L, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FDC_SLF1_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FDC_SLF2_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_Switching$mn_FDC_SLF3_R, SLF_data_subtracted_Switching$Switching.Raw, SLF_data_subtracted_Switching[,c("Age", "Sex_binary")])
#for HayBTime1:
SLF_data_subtracted_HayBTime1 <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "HayBTime1.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                 "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                 "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                 "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                 "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                 "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_HayBTime1 <- SLF_data_subtracted_HayBTime1 %>% drop_na(HayBTime1.Raw)
SLF_data_subtracted_HayBTime1$Sex_binary <- as.numeric(SLF_data_subtracted_HayBTime1$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_HayBTime1$mn_FD_SLF, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime1$mn_FD_SLF_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime1$mn_FD_SLF_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime1$mn_FD_SLF1_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FD_SLF2_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime1$mn_FD_SLF3_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime1$mn_FD_SLF1_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FD_SLF2_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime1$mn_FD_SLF3_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_HayBTime1$mn_FC_SLF, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FC_SLF_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FC_SLF_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FC_SLF1_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FC_SLF2_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FC_SLF3_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FC_SLF1_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FC_SLF2_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FC_SLF3_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_HayBTime1$mn_FDC_SLF, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FDC_SLF_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FDC_SLF_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FDC_SLF1_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FDC_SLF2_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime1$mn_FDC_SLF3_L, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FDC_SLF1_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FDC_SLF2_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime1$mn_FDC_SLF3_R, SLF_data_subtracted_HayBTime1$HayBTime1.Raw, SLF_data_subtracted_HayBTime1[,c("Age", "Sex_binary")])
#for HayBTime2:
SLF_data_subtracted_HayBTime2 <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "HayBTime2.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                 "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                 "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                 "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                 "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                 "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_HayBTime2 <- SLF_data_subtracted_HayBTime2 %>% drop_na(HayBTime2.Raw)
SLF_data_subtracted_HayBTime2$Sex_binary <- as.numeric(SLF_data_subtracted_HayBTime2$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_HayBTime2$mn_FD_SLF, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_subtracted_HayBTime2$mn_FD_SLF_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime2$mn_FD_SLF_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_subtracted_HayBTime2$mn_FD_SLF1_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_subtracted_HayBTime2$mn_FD_SLF2_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime2$mn_FD_SLF3_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime2$mn_FD_SLF1_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FD_SLF2_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime2$mn_FD_SLF3_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) #sig
#for FC
pcor.test(SLF_data_subtracted_HayBTime2$mn_FC_SLF, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FC_SLF_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FC_SLF_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FC_SLF1_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBTime2$mn_FC_SLF2_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FC_SLF3_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FC_SLF1_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FC_SLF2_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FC_SLF3_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_HayBTime2$mn_FDC_SLF, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_subtracted_HayBTime2$mn_FDC_SLF_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_subtracted_HayBTime2$mn_FDC_SLF_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_subtracted_HayBTime2$mn_FDC_SLF1_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_subtracted_HayBTime2$mn_FDC_SLF2_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])  
pcor.test(SLF_data_subtracted_HayBTime2$mn_FDC_SLF3_L, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FDC_SLF1_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FDC_SLF2_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBTime2$mn_FDC_SLF3_R, SLF_data_subtracted_HayBTime2$HayBTime2.Raw, SLF_data_subtracted_HayBTime2[,c("Age", "Sex_binary")]) #sig
#for HayBCatA:
SLF_data_subtracted_HayBCatA <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "HayBCatA.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_HayBCatA <- SLF_data_subtracted_HayBCatA %>% drop_na(HayBCatA.Raw)
SLF_data_subtracted_HayBCatA$Sex_binary <- as.numeric(SLF_data_subtracted_HayBCatA$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_HayBCatA$mn_FD_SLF, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatA$mn_FD_SLF_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatA$mn_FD_SLF_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatA$mn_FD_SLF1_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FD_SLF2_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatA$mn_FD_SLF3_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatA$mn_FD_SLF1_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FD_SLF2_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatA$mn_FD_SLF3_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_HayBCatA$mn_FC_SLF, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FC_SLF_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FC_SLF_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FC_SLF1_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FC_SLF2_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FC_SLF3_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FC_SLF1_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FC_SLF2_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FC_SLF3_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_HayBCatA$mn_FDC_SLF, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FDC_SLF_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FDC_SLF_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FDC_SLF1_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FDC_SLF2_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatA$mn_FDC_SLF3_L, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FDC_SLF1_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FDC_SLF2_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatA$mn_FDC_SLF3_R, SLF_data_subtracted_HayBCatA$HayBCatA.Raw, SLF_data_subtracted_HayBCatA[,c("Age", "Sex_binary")])
#for HayBCatB:
SLF_data_subtracted_HayBCatB <- SLF_data_subtracted[c("Group", "Age", "Sex", "Sex_binary", "HayBCatB.Raw", "mn_FD_SLF", "mn_FD_SLF_L", 
                                "mn_FD_SLF_R", "mn_FD_SLF1_L", "mn_FD_SLF2_L" ,"mn_FD_SLF3_L", "mn_FD_SLF1_R", 
                                "mn_FD_SLF2_R", "mn_FD_SLF3_R", "mn_FC_SLF", "mn_FC_SLF_L", "mn_FC_SLF_R", 
                                "mn_FC_SLF1_L", "mn_FC_SLF2_L" ,"mn_FC_SLF3_L", "mn_FC_SLF1_R", "mn_FC_SLF2_R", 
                                "mn_FC_SLF3_R","mn_FDC_SLF", "mn_FDC_SLF_L", "mn_FDC_SLF_R", "mn_FDC_SLF1_L", 
                                "mn_FDC_SLF2_L","mn_FDC_SLF3_L", "mn_FDC_SLF1_R", "mn_FDC_SLF2_R", "mn_FDC_SLF3_R")]
SLF_data_subtracted_HayBCatB <- SLF_data_subtracted_HayBCatB %>% drop_na(HayBCatB.Raw)
SLF_data_subtracted_HayBCatB$Sex_binary <- as.numeric(SLF_data_subtracted_HayBCatB$Sex_binary)
#for FD
pcor.test(SLF_data_subtracted_HayBCatB$mn_FD_SLF, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatB$mn_FD_SLF_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatB$mn_FD_SLF_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatB$mn_FD_SLF1_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FD_SLF2_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatB$mn_FD_SLF3_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatB$mn_FD_SLF1_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FD_SLF2_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatB$mn_FD_SLF3_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")]) 
#for FC
pcor.test(SLF_data_subtracted_HayBCatB$mn_FC_SLF, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FC_SLF_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FC_SLF_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FC_SLF1_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FC_SLF2_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FC_SLF3_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FC_SLF1_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FC_SLF2_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FC_SLF3_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
#for FDC
pcor.test(SLF_data_subtracted_HayBCatB$mn_FDC_SLF, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FDC_SLF_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FDC_SLF_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FDC_SLF1_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FDC_SLF2_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_subtracted_HayBCatB$mn_FDC_SLF3_L, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FDC_SLF1_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FDC_SLF2_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_subtracted_HayBCatB$mn_FDC_SLF3_R, SLF_data_subtracted_HayBCatB$HayBCatB.Raw, SLF_data_subtracted_HayBCatB[,c("Age", "Sex_binary")])



#####---------- Simple linear regression (F0 to F2 prediction) -------------####
#1. Hayling's Sentence
#HayTime1 vs. SLF FD whole
SLF_FD_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF FC whole
SLF_FC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF FDC whole
SLF_FDC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF L FD 
SLF_L_FD_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF L FC 
SLF_L_FC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF L FDC 
SLF_L_FDC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF R FD 
SLF_R_FD_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF R FC 
SLF_R_FC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF R FDC 
SLF_R_FDC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 1 L FD 
SLF1_L_FD_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 1 L FC 
SLF1_L_FC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_Hay1Raw_mod) #sig
#HayTime1 vs. SLF 1 L FDC 
SLF1_L_FDC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 2 L FD 
SLF2_L_FD_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 2 L FC
SLF2_L_FC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_Hay1Raw_mod) #sig
#HayTime1 vs. SLF 2 L FDC 
SLF2_L_FDC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 3 L FD 
SLF3_L_FD_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 3 L FC 
SLF3_L_FC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 3 L FDC 
SLF3_L_FDC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 1 R FD 
SLF1_R_FD_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_Hay1Raw_mod) #sig
#HayTime1 vs. SLF 1 R FC 
SLF1_R_FC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 1 R FDC 
SLF1_R_FDC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 2 R FD 
SLF2_R_FD_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 2 R FC 
SLF2_R_FC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 2 R FDC 
SLF2_R_FDC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 3 R FD 
SLF3_R_FD_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 3 R FC 
SLF3_R_FC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_Hay1Raw_mod) #not sig
#HayTime1 vs. SLF 3 R FDC 
SLF3_R_FDC_Hay1Raw_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_Hay1Raw_mod) #not sig

#HayTime2 vs. SLF FD whole
SLF_FD_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF FC whole
SLF_FC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF FDC whole
SLF_FDC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF L FD 
SLF_L_FD_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF L FC 
SLF_L_FC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF L FDC 
SLF_L_FDC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF R FD 
SLF_R_FD_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF R FC 
SLF_R_FC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF R FDC 
SLF_R_FDC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 1 L FD 
SLF1_L_FD_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 1 L FC 
SLF1_L_FC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 1 L FDC 
SLF1_L_FDC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 2 L FD 
SLF2_L_FD_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 2 L FC
SLF2_L_FC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 2 L FDC 
SLF2_L_FDC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 3 L FD 
SLF3_L_FD_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 3 L FC 
SLF3_L_FC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 3 L FDC 
SLF3_L_FDC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 1 R FD 
SLF1_R_FD_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_Hay2Raw_mod) #sig
#HayTime2 vs. SLF 1 R FC 
SLF1_R_FC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 1 R FDC 
SLF1_R_FDC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 2 R FD 
SLF2_R_FD_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 2 R FC 
SLF2_R_FC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 2 R FDC 
SLF2_R_FDC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 3 R FD 
SLF3_R_FD_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 3 R FC 
SLF3_R_FC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_Hay2Raw_mod) #not sig
#HayTime2 vs. SLF 3 R FDC 
SLF3_R_FDC_Hay2Raw_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_Hay2Raw_mod) #not sig

#HayCatAError vs. SLF FD whole
SLF_FD_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF FC whole
SLF_FC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF FDC whole
SLF_FDC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF L FD 
SLF_L_FD_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF L FC 
SLF_L_FC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF L FDC 
SLF_L_FDC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF R FD 
SLF_R_FD_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF R FC 
SLF_R_FC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF R FDC 
SLF_R_FDC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 1 L FD 
SLF1_L_FD_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 1 L FC 
SLF1_L_FC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 1 L FDC 
SLF1_L_FDC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 2 L FD 
SLF2_L_FD_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 2 L FC
SLF2_L_FC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 2 L FDC 
SLF2_L_FDC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 3 L FD 
SLF3_L_FD_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 3 L FC 
SLF3_L_FC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 3 L FDC 
SLF3_L_FDC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 1 R FD 
SLF1_R_FD_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 1 R FC 
SLF1_R_FC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 1 R FDC 
SLF1_R_FDC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 2 R FD 
SLF2_R_FD_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 2 R FC 
SLF2_R_FC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 2 R FDC 
SLF2_R_FDC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 3 R FD 
SLF3_R_FD_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 3 R FC 
SLF3_R_FC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_HayCatARaw_mod) #not sig
#HayCatAError vs. SLF 3 R FDC 
SLF3_R_FDC_HayCatARaw_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_HayCatARaw_mod) #not sig

#HayCatBError vs. SLF FD whole
SLF_FD_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF FC whole
SLF_FC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF FDC whole
SLF_FDC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF L FD 
SLF_L_FD_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF L FC 
SLF_L_FC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF L FDC 
SLF_L_FDC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF R FD 
SLF_R_FD_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF R FC 
SLF_R_FC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF R FDC 
SLF_R_FDC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 1 L FD 
SLF1_L_FD_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 1 L FC 
SLF1_L_FC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 1 L FDC 
SLF1_L_FDC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 2 L FD 
SLF2_L_FD_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 2 L FC
SLF2_L_FC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 2 L FDC 
SLF2_L_FDC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 3 L FD 
SLF3_L_FD_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 3 L FC 
SLF3_L_FC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 3 L FDC 
SLF3_L_FDC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 1 R FD 
SLF1_R_FD_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_HayCatBRaw_mod) #sig
#HayCatBError vs. SLF 1 R FC 
SLF1_R_FC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_HayCatBRaw_mod) #sig
#HayCatBError vs. SLF 1 R FDC 
SLF1_R_FDC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_HayCatBRaw_mod) #sig
#HayCatBError vs. SLF 2 R FD 
SLF2_R_FD_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 2 R FC 
SLF2_R_FC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 2 R FDC 
SLF2_R_FDC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 3 R FD 
SLF3_R_FD_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 3 R FC 
SLF3_R_FC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_HayCatBRaw_mod) #not sig
#HayCatBError vs. SLF 3 R FDC 
SLF3_R_FDC_HayCatBRaw_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_HayCatBRaw_mod) #not sig

#2. Stroop Test
#ColourNaming vs. SLF FD whole
SLF_FD_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF FC whole
SLF_FC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF FDC whole
SLF_FDC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF L FD 
SLF_L_FD_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF L FC 
SLF_L_FC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_ColourNamingRaw_mod) #sig
#ColourNaming vs. SLF L FDC 
SLF_L_FDC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF R FD 
SLF_R_FD_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF R FC 
SLF_R_FC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF R FDC 
SLF_R_FDC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 1 L FD 
SLF1_L_FD_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 1 L FC 
SLF1_L_FC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_ColourNamingRaw_mod) #sig
#ColourNaming vs. SLF 1 L FDC 
SLF1_L_FDC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 2 L FD 
SLF2_L_FD_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 2 L FC
SLF2_L_FC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 2 L FDC 
SLF2_L_FDC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 3 L FD 
SLF3_L_FD_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 3 L FC 
SLF3_L_FC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 3 L FDC 
SLF3_L_FDC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 1 R FD 
SLF1_R_FD_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 1 R FC 
SLF1_R_FC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 1 R FDC 
SLF1_R_FDC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 2 R FD 
SLF2_R_FD_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 2 R FC 
SLF2_R_FC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 2 R FDC 
SLF2_R_FDC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 3 R FD 
SLF3_R_FD_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 3 R FC 
SLF3_R_FC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_ColourNamingRaw_mod) #not sig
#ColourNaming vs. SLF 3 R FDC 
SLF3_R_FDC_ColourNamingRaw_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_ColourNamingRaw_mod) #not sig

#WordReading vs. SLF FD whole
SLF_FD_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_WordReadingRaw_mod) #not sig
#WordReading vs. SLF FC whole
SLF_FC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF FDC whole
SLF_FDC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF L FD 
SLF_L_FD_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_WordReadingRaw_mod) #not sig
#WordReading vs. SLF L FC 
SLF_L_FC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF L FDC 
SLF_L_FDC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF R FD 
SLF_R_FD_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_WordReadingRaw_mod) #not sig
#WordReading vs. SLF R FC 
SLF_R_FC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF R FDC 
SLF_R_FDC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 1 L FD 
SLF1_L_FD_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 1 L FC 
SLF1_L_FC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 1 L FDC 
SLF1_L_FDC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 2 L FD 
SLF2_L_FD_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 2 L FC
SLF2_L_FC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 2 L FDC 
SLF2_L_FDC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 3 L FD 
SLF3_L_FD_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 3 L FC 
SLF3_L_FC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 3 L FDC 
SLF3_L_FDC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 1 R FD 
SLF1_R_FD_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 1 R FC 
SLF1_R_FC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 1 R FDC 
SLF1_R_FDC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_WordReadingRaw_mod) #sig
#WordReading vs. SLF 2 R FD 
SLF2_R_FD_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 2 R FC 
SLF2_R_FC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 2 R FDC 
SLF2_R_FDC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 3 R FD 
SLF3_R_FD_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 3 R FC 
SLF3_R_FC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_WordReadingRaw_mod) #not sig
#WordReading vs. SLF 3 R FDC 
SLF3_R_FDC_WordReadingRaw_mod <- lm(WordReading.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_WordReadingRaw_mod) #not sig

#Inhibition vs. SLF FD whole
SLF_FD_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF FC whole
SLF_FC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF FDC whole
SLF_FDC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF L FD 
SLF_L_FD_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF L FC 
SLF_L_FC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF L FDC 
SLF_L_FDC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF R FD 
SLF_R_FD_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF R FC 
SLF_R_FC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF R FDC 
SLF_R_FDC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 1 L FD 
SLF1_L_FD_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 1 L FC 
SLF1_L_FC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 1 L FDC 
SLF1_L_FDC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 2 L FD 
SLF2_L_FD_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 2 L FC
SLF2_L_FC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 2 L FDC 
SLF2_L_FDC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 3 L FD 
SLF3_L_FD_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 3 L FC 
SLF3_L_FC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 3 L FDC 
SLF3_L_FDC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 1 R FD 
SLF1_R_FD_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 1 R FC 
SLF1_R_FC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 1 R FDC 
SLF1_R_FDC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 2 R FD 
SLF2_R_FD_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 2 R FC 
SLF2_R_FC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 2 R FDC 
SLF2_R_FDC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 3 R FD 
SLF3_R_FD_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 3 R FC 
SLF3_R_FC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_InhibitionRaw_mod) #not sig
#Inhibition vs. SLF 3 R FDC 
SLF3_R_FDC_InhibitionRaw_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_InhibitionRaw_mod) #not sig

#3. Trail Making Test (TMT)
#TrailsA vs. SLF FD whole
SLF_FD_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_TrailsARaw_mod) #not sig
#TrailsA vs. SLF FC whole
SLF_FC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF FDC whole
SLF_FDC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF L FD 
SLF_L_FD_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_TrailsARaw_mod) #not sig
#TrailsA vs. SLF L FC 
SLF_L_FC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF L FDC 
SLF_L_FDC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF R FD 
SLF_R_FD_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_TrailsARaw_mod) #not sig
#TrailsA vs. SLF R FC 
SLF_R_FC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF R FDC 
SLF_R_FDC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 1 L FD 
SLF1_L_FD_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 1 L FC 
SLF1_L_FC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 1 L FDC 
SLF1_L_FDC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 2 L FD 
SLF2_L_FD_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 2 L FC
SLF2_L_FC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 2 L FDC 
SLF2_L_FDC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 3 L FD 
SLF3_L_FD_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 3 L FC 
SLF3_L_FC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 3 L FDC 
SLF3_L_FDC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 1 R FD 
SLF1_R_FD_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 1 R FC 
SLF1_R_FC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 1 R FDC 
SLF1_R_FDC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 2 R FD 
SLF2_R_FD_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_TrailsARaw_mod) #sig
#TrailsA vs. SLF 2 R FC 
SLF2_R_FC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 2 R FDC 
SLF2_R_FDC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 3 R FD 
SLF3_R_FD_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 3 R FC 
SLF3_R_FC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_TrailsARaw_mod) #not sig
#TrailsA vs. SLF 3 R FDC 
SLF3_R_FDC_TrailsARaw_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_TrailsARaw_mod) #not sig

#TrailsB vs. SLF FD whole
SLF_FD_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF FC whole
SLF_FC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF FDC whole
SLF_FDC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF L FD 
SLF_L_FD_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF L FC 
SLF_L_FC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF L FDC 
SLF_L_FDC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF R FD 
SLF_R_FD_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF R FC 
SLF_R_FC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF R FDC 
SLF_R_FDC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 1 L FD 
SLF1_L_FD_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 1 L FC 
SLF1_L_FC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 1 L FDC 
SLF1_L_FDC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 2 L FD 
SLF2_L_FD_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 2 L FC
SLF2_L_FC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 2 L FDC 
SLF2_L_FDC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 3 L FD 
SLF3_L_FD_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 3 L FC 
SLF3_L_FC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 3 L FDC 
SLF3_L_FDC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 1 R FD 
SLF1_R_FD_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 1 R FC 
SLF1_R_FC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 1 R FDC 
SLF1_R_FDC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 2 R FD 
SLF2_R_FD_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 2 R FC 
SLF2_R_FC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 2 R FDC 
SLF2_R_FDC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 3 R FD 
SLF3_R_FD_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 3 R FC 
SLF3_R_FC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_TrailsBRaw_mod) #not sig
#TrailsB vs. SLF 3 R FDC 
SLF3_R_FDC_TrailsBRaw_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_TrailsBRaw_mod) #not sig

#LetFluency vs. SLF FD whole
SLF_FD_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF FC whole
SLF_FC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF FDC whole
SLF_FDC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF L FD 
SLF_L_FD_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF L FC 
SLF_L_FC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF L FDC 
SLF_L_FDC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF R FD 
SLF_R_FD_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF R FC 
SLF_R_FC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF R FDC 
SLF_R_FDC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 1 L FD 
SLF1_L_FD_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 1 L FC 
SLF1_L_FC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 1 L FDC 
SLF1_L_FDC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 2 L FD 
SLF2_L_FD_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 2 L FC
SLF2_L_FC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 2 L FDC 
SLF2_L_FDC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 3 L FD 
SLF3_L_FD_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_LetFluencyRaw_mod) #sig
#LetFluency vs. SLF 3 L FC 
SLF3_L_FC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 3 L FDC 
SLF3_L_FDC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_LetFluencyRaw_mod) #sig
#LetFluency vs. SLF 1 R FD 
SLF1_R_FD_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 1 R FC 
SLF1_R_FC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 1 R FDC 
SLF1_R_FDC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 2 R FD 
SLF2_R_FD_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 2 R FC 
SLF2_R_FC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 2 R FDC 
SLF2_R_FDC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 3 R FD 
SLF3_R_FD_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 3 R FC 
SLF3_R_FC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_LetFluencyRaw_mod) #not sig
#LetFluency vs. SLF 3 R FDC 
SLF3_R_FDC_LetFluencyRaw_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_LetFluencyRaw_mod) #not sig

#CatFluency vs. SLF FD whole
SLF_FD_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF FC whole
SLF_FC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF FDC whole
SLF_FDC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF L FD 
SLF_L_FD_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF L FC 
SLF_L_FC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF L FDC 
SLF_L_FDC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF R FD 
SLF_R_FD_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF R FC 
SLF_R_FC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF R FDC 
SLF_R_FDC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 1 L FD 
SLF1_L_FD_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 1 L FC 
SLF1_L_FC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 1 L FDC 
SLF1_L_FDC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 2 L FD 
SLF2_L_FD_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 2 L FC
SLF2_L_FC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 2 L FDC 
SLF2_L_FDC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 3 L FD 
SLF3_L_FD_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 3 L FC 
SLF3_L_FC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 3 L FDC 
SLF3_L_FDC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 1 R FD 
SLF1_R_FD_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 1 R FC 
SLF1_R_FC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 1 R FDC 
SLF1_R_FDC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 2 R FD 
SLF2_R_FD_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 2 R FC 
SLF2_R_FC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 2 R FDC 
SLF2_R_FDC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 3 R FD 
SLF3_R_FD_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 3 R FC 
SLF3_R_FC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_CatFluencyRaw_mod) #not sig
#CatFluency vs. SLF 3 R FDC 
SLF3_R_FDC_CatFluencyRaw_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_CatFluencyRaw_mod) #not sig

#Switching vs. SLF FD whole
SLF_FD_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FD_SLF, data = SLF_data_subtracted)
summary(SLF_FD_SwitchingRaw_mod) #not sig
#Switching vs. SLF FC whole
SLF_FC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FC_SLF, data = SLF_data_subtracted)
summary(SLF_FC_SwitchingRaw_mod) #not sig
#Switching vs. SLF FDC whole
SLF_FDC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FDC_SLF, data = SLF_data_subtracted)
summary(SLF_FDC_SwitchingRaw_mod) #not sig
#Switching vs. SLF L FD 
SLF_L_FD_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FD_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FD_SwitchingRaw_mod) #not sig
#Switching vs. SLF L FC 
SLF_L_FC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FC_SwitchingRaw_mod) #not sig
#Switching vs. SLF L FDC 
SLF_L_FDC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FDC_SLF_L, data = SLF_data_subtracted)
summary(SLF_L_FDC_SwitchingRaw_mod) #not sig
#Switching vs. SLF R FD 
SLF_R_FD_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FD_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FD_SwitchingRaw_mod) #not sig
#Switching vs. SLF R FC 
SLF_R_FC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FC_SwitchingRaw_mod) #not sig
#Switching vs. SLF R FDC 
SLF_R_FDC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FDC_SLF_R, data = SLF_data_subtracted)
summary(SLF_R_FDC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 1 L FD 
SLF1_L_FD_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FD_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FD_SwitchingRaw_mod) #not sig
#Switching vs. SLF 1 L FC 
SLF1_L_FC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 1 L FDC 
SLF1_L_FDC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FDC_SLF1_L, data = SLF_data_subtracted)
summary(SLF1_L_FDC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 2 L FD 
SLF2_L_FD_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FD_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FD_SwitchingRaw_mod) #not sig
#Switching vs. SLF 2 L FC
SLF2_L_FC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 2 L FDC 
SLF2_L_FDC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FDC_SLF2_L, data = SLF_data_subtracted)
summary(SLF2_L_FDC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 3 L FD 
SLF3_L_FD_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FD_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FD_SwitchingRaw_mod) #not sig
#Switching vs. SLF 3 L FC 
SLF3_L_FC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 3 L FDC 
SLF3_L_FDC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FDC_SLF3_L, data = SLF_data_subtracted)
summary(SLF3_L_FDC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 1 R FD 
SLF1_R_FD_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FD_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FD_SwitchingRaw_mod) #not sig
#Switching vs. SLF 1 R FC 
SLF1_R_FC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 1 R FDC 
SLF1_R_FDC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FDC_SLF1_R, data = SLF_data_subtracted)
summary(SLF1_R_FDC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 2 R FD 
SLF2_R_FD_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FD_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FD_SwitchingRaw_mod) #not sig
#Switching vs. SLF 2 R FC 
SLF2_R_FC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 2 R FDC 
SLF2_R_FDC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FDC_SLF2_R, data = SLF_data_subtracted)
summary(SLF2_R_FDC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 3 R FD 
SLF3_R_FD_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FD_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FD_SwitchingRaw_mod) #not sig
#Switching vs. SLF 3 R FC 
SLF3_R_FC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FC_SwitchingRaw_mod) #not sig
#Switching vs. SLF 3 R FDC 
SLF3_R_FDC_SwitchingRaw_mod <- lm(Switching.Raw ~ mn_FDC_SLF3_R, data = SLF_data_subtracted)
summary(SLF3_R_FDC_SwitchingRaw_mod) #not sig


#Run Linear Regression w/ covariates (age & sex) ------------------------------#
#1. Hayling's Sentence
#HayTime1 vs. SLF FD whole
SLF_FD_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF + Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF FC whole
SLF_FC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF FDC whole
SLF_FDC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF L FD 
SLF_L_FD_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF L FC 
SLF_L_FC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_Hay1Raw_2covar_mod) #sig
#HayTime1 vs. SLF L FDC 
SLF_L_FDC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF R FD 
SLF_R_FD_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF R FC 
SLF_R_FC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF R FDC 
SLF_R_FDC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 1 L FD 
SLF1_L_FD_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 1 L FC 
SLF1_L_FC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 1 L FDC 
SLF1_L_FDC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 2 L FD 
SLF2_L_FD_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 2 L FC
SLF2_L_FC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_Hay1Raw_2covar_mod) #sig
#HayTime1 vs. SLF 2 L FDC 
SLF2_L_FDC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 3 L FD 
SLF3_L_FD_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 3 L FC 
SLF3_L_FC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 3 L FDC 
SLF3_L_FDC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 1 R FD 
SLF1_R_FD_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_Hay1Raw_2covar_mod) #sig
#HayTime1 vs. SLF 1 R FC 
SLF1_R_FC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 1 R FDC 
SLF1_R_FDC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 2 R FD 
SLF2_R_FD_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 2 R FC 
SLF2_R_FC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 2 R FDC 
SLF2_R_FDC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 3 R FD 
SLF3_R_FD_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 3 R FC 
SLF3_R_FC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_Hay1Raw_2covar_mod) #not sig
#HayTime1 vs. SLF 3 R FDC 
SLF3_R_FDC_Hay1Raw_2covar_mod <- lm(HayBTime1.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_Hay1Raw_2covar_mod) #not sig

#HayTime2 vs. SLF FD whole
SLF_FD_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF FC whole
SLF_FC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF FDC whole
SLF_FDC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF L FD 
SLF_L_FD_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF L FC 
SLF_L_FC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF L FDC 
SLF_L_FDC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF R FD 
SLF_R_FD_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF R FC 
SLF_R_FC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF R FDC 
SLF_R_FDC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 1 L FD 
SLF1_L_FD_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 1 L FC 
SLF1_L_FC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 1 L FDC 
SLF1_L_FDC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 2 L FD 
SLF2_L_FD_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 2 L FC
SLF2_L_FC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 2 L FDC 
SLF2_L_FDC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 3 L FD 
SLF3_L_FD_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 3 L FC 
SLF3_L_FC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 3 L FDC 
SLF3_L_FDC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 1 R FD 
SLF1_R_FD_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_Hay2Raw_2covar_mod) #sig
#HayTime2 vs. SLF 1 R FC 
SLF1_R_FC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 1 R FDC 
SLF1_R_FDC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 2 R FD 
SLF2_R_FD_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 2 R FC 
SLF2_R_FC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 2 R FDC 
SLF2_R_FDC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 3 R FD 
SLF3_R_FD_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 3 R FC 
SLF3_R_FC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_Hay2Raw_2covar_mod) #not sig
#HayTime2 vs. SLF 3 R FDC 
SLF3_R_FDC_Hay2Raw_2covar_mod <- lm(HayBTime2.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_Hay2Raw_2covar_mod) #not sig

#HayCatAError vs. SLF FD whole
SLF_FD_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF FC whole
SLF_FC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF FDC whole
SLF_FDC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF L FD 
SLF_L_FD_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF L FC 
SLF_L_FC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF L FDC 
SLF_L_FDC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF R FD 
SLF_R_FD_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF R FC 
SLF_R_FC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF R FDC 
SLF_R_FDC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 1 L FD 
SLF1_L_FD_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 1 L FC 
SLF1_L_FC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 1 L FDC 
SLF1_L_FDC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 2 L FD 
SLF2_L_FD_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 2 L FC
SLF2_L_FC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 2 L FDC 
SLF2_L_FDC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 3 L FD 
SLF3_L_FD_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 3 L FC 
SLF3_L_FC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 3 L FDC 
SLF3_L_FDC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 1 R FD 
SLF1_R_FD_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 1 R FC 
SLF1_R_FC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 1 R FDC 
SLF1_R_FDC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 2 R FD 
SLF2_R_FD_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 2 R FC 
SLF2_R_FC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 2 R FDC 
SLF2_R_FDC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 3 R FD 
SLF3_R_FD_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 3 R FC 
SLF3_R_FC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_HayCatARaw_2covar_mod) #not sig
#HayCatAError vs. SLF 3 R FDC 
SLF3_R_FDC_HayCatARaw_2covar_mod <- lm(HayBCatA.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_HayCatARaw_2covar_mod) #not sig

#HayCatBError vs. SLF FD whole
SLF_FD_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF FC whole
SLF_FC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF FDC whole
SLF_FDC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF L FD 
SLF_L_FD_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF L FC 
SLF_L_FC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF L FDC 
SLF_L_FDC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF R FD 
SLF_R_FD_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF R FC 
SLF_R_FC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF R FDC 
SLF_R_FDC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 1 L FD 
SLF1_L_FD_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 1 L FC 
SLF1_L_FC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 1 L FDC 
SLF1_L_FDC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 2 L FD 
SLF2_L_FD_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 2 L FC
SLF2_L_FC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 2 L FDC 
SLF2_L_FDC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 3 L FD 
SLF3_L_FD_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 3 L FC 
SLF3_L_FC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 3 L FDC 
SLF3_L_FDC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 1 R FD 
SLF1_R_FD_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_HayCatBRaw_2covar_mod) #sig
#HayCatBError vs. SLF 1 R FC 
SLF1_R_FC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_HayCatBRaw_2covar_mod) #sig
#HayCatBError vs. SLF 1 R FDC 
SLF1_R_FDC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_HayCatBRaw_2covar_mod) #sig
#HayCatBError vs. SLF 2 R FD 
SLF2_R_FD_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 2 R FC 
SLF2_R_FC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 2 R FDC 
SLF2_R_FDC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 3 R FD 
SLF3_R_FD_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_HayCatBRaw_2covar_mod) #not sig
#HayCatBError vs. SLF 3 R FC 
SLF3_R_FC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_HayCatBRaw_2covar_mod) #sig
#HayCatBError vs. SLF 3 R FDC 
SLF3_R_FDC_HayCatBRaw_2covar_mod <- lm(HayBCatB.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_HayCatBRaw_2covar_mod) #not sig

#2. Stroop Test
#ColourNaming vs. SLF FD whole
SLF_FD_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF FC whole
SLF_FC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF FDC whole
SLF_FDC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF L FD 
SLF_L_FD_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF L FC 
SLF_L_FC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_ColourNamingRaw_2covar_mod) #sig
#ColourNaming vs. SLF L FDC 
SLF_L_FDC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF R FD 
SLF_R_FD_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF R FC 
SLF_R_FC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF R FDC 
SLF_R_FDC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 1 L FD 
SLF1_L_FD_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 1 L FC 
SLF1_L_FC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_ColourNamingRaw_2covar_mod) #sig
#ColourNaming vs. SLF 1 L FDC 
SLF1_L_FDC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 2 L FD 
SLF2_L_FD_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 2 L FC
SLF2_L_FC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 2 L FDC 
SLF2_L_FDC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 3 L FD 
SLF3_L_FD_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 3 L FC 
SLF3_L_FC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 3 L FDC 
SLF3_L_FDC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 1 R FD 
SLF1_R_FD_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 1 R FC 
SLF1_R_FC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 1 R FDC 
SLF1_R_FDC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 2 R FD 
SLF2_R_FD_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 2 R FC 
SLF2_R_FC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 2 R FDC 
SLF2_R_FDC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 3 R FD 
SLF3_R_FD_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 3 R FC 
SLF3_R_FC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_ColourNamingRaw_2covar_mod) #not sig
#ColourNaming vs. SLF 3 R FDC 
SLF3_R_FDC_ColourNamingRaw_2covar_mod <- lm(ColorNaming.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_ColourNamingRaw_2covar_mod) #not sig

#WordReading vs. SLF FD whole
SLF_FD_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF FC whole
SLF_FC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF FDC whole
SLF_FDC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF L FD 
SLF_L_FD_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF L FC 
SLF_L_FC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF L FDC 
SLF_L_FDC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF R FD 
SLF_R_FD_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF R FC 
SLF_R_FC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF R FDC 
SLF_R_FDC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 1 L FD 
SLF1_L_FD_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 1 L FC 
SLF1_L_FC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 1 L FDC 
SLF1_L_FDC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 2 L FD 
SLF2_L_FD_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 2 L FC
SLF2_L_FC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 2 L FDC 
SLF2_L_FDC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 3 L FD 
SLF3_L_FD_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 3 L FC 
SLF3_L_FC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 3 L FDC 
SLF3_L_FDC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 1 R FD 
SLF1_R_FD_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 1 R FC 
SLF1_R_FC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_WordReadingRaw_2covar_mod) #sig
#WordReading vs. SLF 1 R FDC 
SLF1_R_FDC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_WordReadingRaw_2covar_mod) #sig
#WordReading vs. SLF 2 R FD 
SLF2_R_FD_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 2 R FC 
SLF2_R_FC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 2 R FDC 
SLF2_R_FDC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 3 R FD 
SLF3_R_FD_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 3 R FC 
SLF3_R_FC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_WordReadingRaw_2covar_mod) #not sig
#WordReading vs. SLF 3 R FDC 
SLF3_R_FDC_WordReadingRaw_2covar_mod <- lm(WordReading.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_WordReadingRaw_2covar_mod) #not sig

#Inhibition vs. SLF FD whole
SLF_FD_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF FC whole
SLF_FC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF FDC whole
SLF_FDC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF L FD 
SLF_L_FD_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF L FC 
SLF_L_FC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF L FDC 
SLF_L_FDC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF R FD 
SLF_R_FD_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF R FC 
SLF_R_FC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF R FDC 
SLF_R_FDC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 1 L FD 
SLF1_L_FD_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 1 L FC 
SLF1_L_FC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 1 L FDC 
SLF1_L_FDC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 2 L FD 
SLF2_L_FD_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 2 L FC
SLF2_L_FC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 2 L FDC 
SLF2_L_FDC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 3 L FD 
SLF3_L_FD_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 3 L FC 
SLF3_L_FC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 3 L FDC 
SLF3_L_FDC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 1 R FD 
SLF1_R_FD_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 1 R FC 
SLF1_R_FC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 1 R FDC 
SLF1_R_FDC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 2 R FD 
SLF2_R_FD_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 2 R FC 
SLF2_R_FC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 2 R FDC 
SLF2_R_FDC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 3 R FD 
SLF3_R_FD_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 3 R FC 
SLF3_R_FC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_InhibitionRaw_2covar_mod) #not sig
#Inhibition vs. SLF 3 R FDC 
SLF3_R_FDC_InhibitionRaw_2covar_mod <- lm(Inhibition.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_InhibitionRaw_2covar_mod) #not sig

#3. Trail Making Test (TMT)
#TrailsA vs. SLF FD whole
SLF_FD_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF FC whole
SLF_FC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF FDC whole
SLF_FDC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF L FD 
SLF_L_FD_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF L FC 
SLF_L_FC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF L FDC 
SLF_L_FDC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF R FD 
SLF_R_FD_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF R FC 
SLF_R_FC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF R FDC 
SLF_R_FDC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 1 L FD 
SLF1_L_FD_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 1 L FC 
SLF1_L_FC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 1 L FDC 
SLF1_L_FDC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 2 L FD 
SLF2_L_FD_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 2 L FC
SLF2_L_FC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 2 L FDC 
SLF2_L_FDC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 3 L FD 
SLF3_L_FD_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 3 L FC 
SLF3_L_FC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 3 L FDC 
SLF3_L_FDC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 1 R FD 
SLF1_R_FD_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 1 R FC 
SLF1_R_FC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 1 R FDC 
SLF1_R_FDC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 2 R FD 
SLF2_R_FD_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_TrailsARaw_2covar_mod) #sig
#TrailsA vs. SLF 2 R FC 
SLF2_R_FC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 2 R FDC 
SLF2_R_FDC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 3 R FD 
SLF3_R_FD_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 3 R FC 
SLF3_R_FC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_TrailsARaw_2covar_mod) #not sig
#TrailsA vs. SLF 3 R FDC 
SLF3_R_FDC_TrailsARaw_2covar_mod <- lm(TrailsA.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_TrailsARaw_2covar_mod) #not sig

#TrailsB vs. SLF FD whole
SLF_FD_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF FC whole
SLF_FC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF FDC whole
SLF_FDC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF L FD 
SLF_L_FD_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF L FC 
SLF_L_FC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF L FDC 
SLF_L_FDC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF R FD 
SLF_R_FD_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF R FC 
SLF_R_FC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF R FDC 
SLF_R_FDC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 1 L FD 
SLF1_L_FD_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 1 L FC 
SLF1_L_FC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 1 L FDC 
SLF1_L_FDC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 2 L FD 
SLF2_L_FD_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 2 L FC
SLF2_L_FC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 2 L FDC 
SLF2_L_FDC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 3 L FD 
SLF3_L_FD_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 3 L FC 
SLF3_L_FC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 3 L FDC 
SLF3_L_FDC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 1 R FD 
SLF1_R_FD_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 1 R FC 
SLF1_R_FC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 1 R FDC 
SLF1_R_FDC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 2 R FD 
SLF2_R_FD_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 2 R FC 
SLF2_R_FC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 2 R FDC 
SLF2_R_FDC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 3 R FD 
SLF3_R_FD_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 3 R FC 
SLF3_R_FC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_TrailsBRaw_2covar_mod) #not sig
#TrailsB vs. SLF 3 R FDC 
SLF3_R_FDC_TrailsBRaw_2covar_mod <- lm(TrailsB.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_TrailsBRaw_2covar_mod) #not sig

#LetFluency vs. SLF FD whole
SLF_FD_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF FC whole
SLF_FC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF FDC whole
SLF_FDC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF L FD 
SLF_L_FD_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF L FC 
SLF_L_FC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF L FDC 
SLF_L_FDC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF R FD 
SLF_R_FD_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF R FC 
SLF_R_FC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF R FDC 
SLF_R_FDC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 1 L FD 
SLF1_L_FD_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 1 L FC 
SLF1_L_FC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 1 L FDC 
SLF1_L_FDC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 2 L FD 
SLF2_L_FD_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 2 L FC
SLF2_L_FC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 2 L FDC 
SLF2_L_FDC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 3 L FD 
SLF3_L_FD_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 3 L FC 
SLF3_L_FC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 3 L FDC 
SLF3_L_FDC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 1 R FD 
SLF1_R_FD_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 1 R FC 
SLF1_R_FC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 1 R FDC 
SLF1_R_FDC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 2 R FD 
SLF2_R_FD_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 2 R FC 
SLF2_R_FC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 2 R FDC 
SLF2_R_FDC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 3 R FD 
SLF3_R_FD_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 3 R FC 
SLF3_R_FC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_LetFluencyRaw_2covar_mod) #not sig
#LetFluency vs. SLF 3 R FDC 
SLF3_R_FDC_LetFluencyRaw_2covar_mod <- lm(LetFluency.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_LetFluencyRaw_2covar_mod) #not sig

#CatFluency vs. SLF FD whole
SLF_FD_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF FC whole
SLF_FC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF FDC whole
SLF_FDC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF L FD 
SLF_L_FD_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF L FC 
SLF_L_FC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF L FDC 
SLF_L_FDC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF R FD 
SLF_R_FD_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF R FC 
SLF_R_FC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF R FDC 
SLF_R_FDC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 1 L FD 
SLF1_L_FD_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 1 L FC 
SLF1_L_FC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 1 L FDC 
SLF1_L_FDC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 2 L FD 
SLF2_L_FD_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 2 L FC
SLF2_L_FC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 2 L FDC 
SLF2_L_FDC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 3 L FD 
SLF3_L_FD_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 3 L FC 
SLF3_L_FC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 3 L FDC 
SLF3_L_FDC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 1 R FD 
SLF1_R_FD_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 1 R FC 
SLF1_R_FC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 1 R FDC 
SLF1_R_FDC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 2 R FD 
SLF2_R_FD_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 2 R FC 
SLF2_R_FC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 2 R FDC 
SLF2_R_FDC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 3 R FD 
SLF3_R_FD_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 3 R FC 
SLF3_R_FC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_CatFluencyRaw_2covar_mod) #not sig
#CatFluency vs. SLF 3 R FDC 
SLF3_R_FDC_CatFluencyRaw_2covar_mod <- lm(CatFluency.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_CatFluencyRaw_2covar_mod) #not sig

#Switching vs. SLF FD whole
SLF_FD_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FD_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FD_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF FC whole
SLF_FC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF FDC whole
SLF_FDC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FDC_SLF+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_FDC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF L FD 
SLF_L_FD_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FD_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FD_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF L FC 
SLF_L_FC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF L FDC 
SLF_L_FDC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FDC_SLF_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_L_FDC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF R FD 
SLF_R_FD_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FD_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FD_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF R FC 
SLF_R_FC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF R FDC 
SLF_R_FDC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FDC_SLF_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF_R_FDC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 1 L FD 
SLF1_L_FD_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FD_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FD_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 1 L FC 
SLF1_L_FC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 1 L FDC 
SLF1_L_FDC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FDC_SLF1_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_L_FDC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 2 L FD 
SLF2_L_FD_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FD_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FD_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 2 L FC
SLF2_L_FC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 2 L FDC 
SLF2_L_FDC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FDC_SLF2_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_L_FDC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 3 L FD 
SLF3_L_FD_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FD_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FD_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 3 L FC 
SLF3_L_FC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 3 L FDC 
SLF3_L_FDC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FDC_SLF3_L+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_L_FDC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 1 R FD 
SLF1_R_FD_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FD_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FD_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 1 R FC 
SLF1_R_FC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 1 R FDC 
SLF1_R_FDC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FDC_SLF1_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF1_R_FDC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 2 R FD 
SLF2_R_FD_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FD_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FD_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 2 R FC 
SLF2_R_FC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 2 R FDC 
SLF2_R_FDC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FDC_SLF2_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF2_R_FDC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 3 R FD 
SLF3_R_FD_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FD_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FD_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 3 R FC 
SLF3_R_FC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FC_SwitchingRaw_2covar_mod) #not sig
#Switching vs. SLF 3 R FDC 
SLF3_R_FDC_SwitchingRaw_2covar_mod <- lm(Switching.Raw ~ mn_FDC_SLF3_R+ Age + Sex, data = SLF_data_subtracted)
summary(SLF3_R_FDC_SwitchingRaw_2covar_mod) #not sig





