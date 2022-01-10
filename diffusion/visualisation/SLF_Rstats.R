#Perform statistics on tract of interest (TOI) from diffusion analysis. This script will perform analysis for group 
#differences in the TOI which was selected from the user in the previous pipelines, done mainly from the software
#programme, MRtrix3. The main outputs are the text files of the fixel-based analysis (FBA) metrics from the TOI.m pipleine. 
#For this script, I have chosen to analyse the superior longitudinal fasciculus (SLF) between the 5 participant groups of 
#interest - Controls (C), subjective cognitive decline (SCD), amnestic mild cognitive impairment (aMCI), multiple-domain
#mild cognitive impairment (mMCI), and Alzheimer's disease (AD). 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 8/12/20


#------------------------------Setting up--------------------------------------#
#install packages/open libraries
pacman::p_load(dplyr, ggplot2, psych, car, multcomp, lsr, tidyr, BayesFactor, tidyverse, ppcor, nlme, effectsize, rstatix)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#first read in the covariates group data file: 
#setwd('/yourpathway/')
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')


#-----------------------Cross-sectional analysis-------------------------------#



DPRC_neuropsych_data <- read.csv("cross-sectional_DPRC_neuropsych_data_lined_up_valid_participants.csv")
#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'

#convert variables
DPRC_neuropsych_data$ParticipantID <- as.factor(DPRC_neuropsych_data$ParticipantID)
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Sex_binary <- as.factor(DPRC_neuropsych_data$Sex_binary)
DPRC_neuropsych_data$Sex<- as.factor(DPRC_neuropsych_data$Sex)

#navigate to the correct pathway which contains the SLF metric text files: 
setwd('V:/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/template/TOI')

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


#----------------------------Descriptives--------------------------------------#
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



#----------------------------ANOVA testing-------------------------------------#

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
confint(post_hoc_SLF3_FD_mod_R)
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

post_hoc_SLF2_R_FD_ancova_mod <- glht(SLF2_FD_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_R_FD_ancova_mod)
confint(post_hoc_SLF2_R_FD_ancova_mod)

post_hoc_SLF3_R_FD_ancova_mod <- glht(SLF3_FD_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_R_FD_ancova_mod)
confint(post_hoc_SLF3_R_FD_ancova_mod)

#plot 95% Confidence Interval
SLF_FD_95CI_data <- data.frame(SLF_group_number = c('1','1','1','2','2','3','3','3','4','4','4','5','5','6','6','7','7','7','7'),
                                      SLF_type = c('Whole_SLF', 'Whole_SLF','Whole_SLF','Left_SLF','Left_SLF', 'Right_SLF', 'Right_SLF', 'Right_SLF','Left_SLF2','Left_SLF2','Left_SLF2','Left_SLF3','Left_SLF3','Right_SLF2','Right_SLF2','Right_SLF3', 'Right_SLF3', 'Right_SLF3', 'Right_SLF3'),
                                      Group_contrast = c('CvaMCI', 'CvmMCI', 'CvAD', 'CvmMCI', 'CvAD','CvaMCI','CvmMCI','CvAD','CvmMCI','CvAD','SCDvAD','CvmMCI','CvAD','CvmMCI','CvAD','CvSCD','CvaMCI','CvmMCI','CvAD'),
                                      estimate_diff = c(0.01459, 0.0195, 0.02462, 0.019266, 0.0253978, 0.0152266,0.0197227,0.0238788, 0.0242823,0.0341201,0.0209306,0.0187693,0.0269245,0.022054,0.020964,0.016077,0.017661,0.022054,0.020964),
                                      lower = c(0.02825, 0.03331, 0.0408, 0.0339826, 0.0426387, 0.0295488,0.0342045,0.0408449,0.0438638,0.0570606, 0.0416857,0.0342775,0.0450929,0.036944,0.038408,0.030567,0.032392,0.03695,0.038415), 
                                      upper = c(0.0009389, 0.005692, 0.008448, 0.0045495, 0.0081568, 0.0009045, 0.0052408,0.0069128,0.0047008,0.0111796,0.0001755,0.0032611,0.008756,0.007164,0.00352,0.001586,0.00293,0.007159,0.003514))  
#plot data
ggplot(SLF_FD_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Whole SLF", "2" = "Left SLF", "3" = "Right SLF", "4" = "Left SLF 2", "5" = "Left SLF 3", "6" = "Right SLF 2", "7" = "Right SLF 3"))+
    theme_classic() 

#plot 95% Confidence Interval for ANCOVA (age & sex)
#for FD
#Create dataset: 
SLF_FD_95CI_2covar_data <- data.frame(SLF_group_number = c('1','1','2','2','3','3','4','5','6','7','7','7'),
                               SLF_type = c('Whole_SLF','Whole_SLF','Left_SLF','Left_SLF', 'Right_SLF', 'Right_SLF','Left_SLF2','Left_SLF3','Right_SLF2','Right_SLF3','Right_SLF3','Right_SLF3'),
                               Group_contrast = c('CvmMCI', 'CvAD', 'CvmMCI', 'CvAD','CvmMCI','CvAD','CvAD','CvAD','CvmMCI','CvSCD','CvaMCI','CvmMCI'),
                               estimate_diff = c(0.016502,0.018702,0.015957,0.018917,0.0170249,0.0184954,0.02323,0.02038,0.020692,0.015,0.01504,0.02016),
                               lower = c(0.0302852,0.0353254,0.030638,0.036623,0.0315608,0.0360265,0.04646,0.03906,0.041194,0.02943,0.03,0.03523), 
                               upper = c(0.0027188,0.0020787,0.001277,0.001211,0.0024891,0.0009642,0.00000159,0.001693,0.000191,0.0005687,0.00008587,0.005099))  
#plot data
ggplot(SLF_FD_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Whole SLF", "2" = "Left SLF", "3" = "Right SLF", "4" = "Left SLF 2", "5" = "Left SLF 3", "6" = "Right SLF 2", "7" = "Right SLF 3"))+
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
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_L$Group)
    cohensD(mn_FD_SLF3_L ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF3_L) #this looks like Hedges' g? 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF3_L$Group)
    cohensD(mn_FD_SLF3_L ~ Group, data = DPRC_neuropsych_data_CvAD_SLF3_L) #this looks like Hedges' g? 
#Right SLF2 FD
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
    
#For ANCOVA - effect size for sig. post hoc tests (Cohen's d)
#Whole SLF FD
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_wholeSLF <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_wholeSLF$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_wholeSLF$Group)
    cohensD(mn_FD_SLF ~ Group + Age + Sex, data = DPRC_neuropsych_data_CvmMCI_wholeSLF) #this looks like Hedges' g? 
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
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_L$Group)
    cohensD(mn_FD_SLF3_L ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF3_L) #this looks like Hedges' g? 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF3_L$Group)
    cohensD(mn_FD_SLF3_L ~ Group, data = DPRC_neuropsych_data_CvAD_SLF3_L) #this looks like Hedges' g? 
#Right SLF2 FD
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
summary(post_hoc_SLF2_L_FDC_ancova_mod)
confint(post_hoc_SLF2_L_FDC_ancova_mod)
post_hoc_SLF3_R_FDC_ancova_mod <- glht(SLF3_FDC_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_R_FDC_ancova_mod)
confint(post_hoc_SLF3_R_FDC_ancova_mod)


#plot 95% Confidence Interval
#for FDC
#Create dataset: 
SLF_FDC_95CI_data <- data.frame(SLF_group_number = c('1','2'),
                                SLF_type = c('Left_SLF2','Right_SLF3'),
                                Group_contrast = c('CvAD', 'CvSCD'),
                                estimate_diff = c(0.045921,0.0340832),
                                lower = c(0.0865456,0.0650291), 
                                upper = c(0.0052964,0.0031373)) 
#plot data
ggplot(SLF_FDC_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Left SLF 2", "2" = "Right SLF 3"))+
    theme_classic() 
#plot 95% Confidence Interval - with covariates (age + sex)
#for FDC
#Create dataset: 
SLF_FDC_95CI_2covar_data <- data.frame(SLF_group_number = c('1','2','2','2'),
                                SLF_type = c('Left_SLF2','Right_SLF3','Right_SLF3','Right_SLF3'),
                                Group_contrast = c('CvAD', 'CvSCD','CvaMCI','CvmMCI'),
                                estimate_diff = c(0.0478704,0.0357782,0.0337897,0.0328381),
                                lower = c(0.0897347,0.0664484,0.0655782,0.0648515), 
                                upper = c(0.0060062,0.0051081,0.0020013,0.0008248)) 
 #plot data
ggplot(SLF_FDC_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Left SLF 2", "2" = "Right SLF 3"))+
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
SLF_data_FD <- dplyr::select(SLF_data, 
                               ParticipantID,
                               Group,
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
    xlab("SLF tract") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("mn_FDC_SLF" = "Whole SLF", "mn_FDC_SLF_L" = "Left SLF", "mn_FDC_SLF_R" = "Right SLF", "mn_FDC_SLF1_L" = "Left SLF1","mn_FDC_SLF2_L" = "Left SLF2","mn_FDC_SLF3_L" = "Left SLF3","mn_FDC_SLF1_R" = "Right SLF1","mn_FDC_SLF2_R" = "Right SLF2","mn_FDC_SLF3_R" = "Right SLF3")) + 
    scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
    theme_classic()



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
cor.test(SLF_data$mn_FD_SLF_R, SLF_data$CatFluency.Raw)#not sig.
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
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FD_SLF1_R, SLF_data$HayBCatA.Raw) #not sig.
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
cor.test(SLF_data$mn_FD_SLF3_R, SLF_data$LetFluency.Raw)#not sig.
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
cor.test(SLF_data$mn_FC_SLF1_R, SLF_data$HayBTime2.Raw)
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
cor.test(SLF_data$mn_FDC_SLF, SLF_data$HayBTime2.Raw)#not sig.
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
cor.test(SLF_data$mn_FDC_SLF_L, SLF_data$Switching.Raw)#not sig.
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
cor.test(SLF_data$mn_FDC_SLF2_L, SLF_data$CatFluency.Raw)#not sig.
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
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$LetFluency.Raw)
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$HayBTime2.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$HayBCatA.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF3_L, SLF_data$HayBCatB.Raw) #not sig.
#Right SLF1-
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$TrailsA.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$TrailsB.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$ColorNaming.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$WordReading.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$Inhibition.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$LetFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$CatFluency.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$Switching.Raw)#not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF1_R, SLF_data$HayBTime2.Raw)
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
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$Switching.Raw)
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$HayBTime1.Raw) #not sig.
cor.test(SLF_data$mn_FDC_SLF2_R, SLF_data$HayBTime2.Raw)#not sig.
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
pcor.test(SLF_data_TrailsB$mn_FD_SLF3_L, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_TrailsB$mn_FD_SLF1_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(SLF_data_TrailsB$mn_FD_SLF2_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_TrailsB$mn_FD_SLF3_R, SLF_data_TrailsB$TrailsB.Raw, SLF_data_TrailsB[,c("Age", "Sex_binary")]) #sig
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
pcor.test(SLF_data_ColorNaming$mn_FD_SLF_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) #sig
pcor.test(SLF_data_ColorNaming$mn_FD_SLF_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_ColorNaming$mn_FD_SLF1_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FD_SLF2_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FD_SLF3_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_ColorNaming$mn_FD_SLF1_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(SLF_data_ColorNaming$mn_FD_SLF2_R, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")]) 
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
pcor.test(SLF_data_ColorNaming$mn_FDC_SLF2_L, SLF_data_ColorNaming$ColorNaming.Raw, SLF_data_ColorNaming[,c("Age", "Sex_binary")])
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
pcor.test(SLF_data_Inhibition$mn_FD_SLF2_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
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
pcor.test(SLF_data_Inhibition$mn_FDC_SLF3_L, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])#sig
pcor.test(SLF_data_Inhibition$mn_FDC_SLF1_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(SLF_data_Inhibition$mn_FDC_SLF2_R, SLF_data_Inhibition$Inhibition.Raw, SLF_data_Inhibition[,c("Age", "Sex_binary")]) #sig
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
pcor.test(SLF_data_LetFluency$mn_FD_SLF, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_LetFluency$mn_FD_SLF_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_LetFluency$mn_FD_SLF_R, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")]) 
pcor.test(SLF_data_LetFluency$mn_FD_SLF1_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_LetFluency$mn_FD_SLF2_L, SLF_data_LetFluency$LetFluency.Raw, SLF_data_LetFluency[,c("Age", "Sex_binary")])
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
pcor.test(SLF_data_CatFluency$mn_FD_SLF2_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(SLF_data_CatFluency$mn_FD_SLF3_L, SLF_data_CatFluency$CatFluency.Raw, SLF_data_CatFluency[,c("Age", "Sex_binary")]) #sig
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





#--------------------Longitudinal analysis (F0 vs. F2)-------------------------#

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
setwd('V:/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/template/TOI')

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
SLF_data$Individual_number <- DPRC_neuropsych_data$Individual_number
#add in Group classification column for each participant - data from another worksheet
SLF_data$Group <- DPRC_neuropsych_data$Group
#add in Timepoint (F0 vs. F2)
SLF_data$Timepoint <- DPRC_neuropsych_data$Timepoint
#add in covariates (clinical site) from the covariates_data dataframe to the SLF dataframe
SLF_data$ClinSite_name <- DPRC_neuropsych_data$Clinical_site 
SLF_data$Age <- DPRC_neuropsych_data$Age 
SLF_data$Sex <- DPRC_neuropsych_data$Sex
SLF_data$Sex_binary <- DPRC_neuropsych_data$Sex_binary
#add in neuropsych variables
neuropsych_test_names <- c("TrailsA.Raw","TrailsA.Z","TrailsB.Raw","TrailsB.Z","ColorNaming.Raw","ColorNaming.Z","WordReading.Raw","WordReading.Z","Inhibition.Raw","Inhibition.Z","LetFluency.Raw","LetFluency.Z","CatFluency.Raw","CatFluency.Z","Switching.Raw","Switching.z","HayBTime1.Raw","HayBTime1.z","HayBTime2.Raw","HayBTime2.z","HayBCatA.Raw","HayBCatA.z","HayBCatB.Raw","HayBCatB.z")
SLF_data[neuropsych_test_names] <- DPRC_neuropsych_data[neuropsych_test_names]


#----------------------------Descriptives--------------------------------------#
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
    
#----------------------------ANOVA testing-------------------------------------#
    
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
get_anova_table(aov_SLF2_L_FD_mod) #significant on main effect Timepoint
aov_SLF2_L_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF2_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_FC_mod) 
aov_SLF2_L_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF2_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_FDC_mod)  #significant on main effect Timepoint
#for left SLF 3
aov_SLF3_L_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF3_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_FD_mod) 
aov_SLF3_L_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF3_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_FC_mod) 
aov_SLF3_L_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF3_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_FDC_mod) 
#for right SLF 1
aov_SLF1_R_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF1_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_R_FD_mod)  #sig. interaction
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

#follow up tests with sig. main effects and interactions
#Left SLF 2 FD
#non-sig. interaction - test by Time Point
SLF_data %>%
    pairwise_t_test(
        mn_FD_SLF2_L ~ Timepoint, paired = TRUE,
        p.adjust.method = "fdr"
    )
#effect size for Time Point
SLF_data%>%cohens_d(mn_FD_SLF2_L~Timepoint, paired=TRUE)
#Left SLF 2 FDC
#non-sig. interaction - test by Time Point
SLF_data %>%
    pairwise_t_test(
        mn_FDC_SLF2_L ~ Timepoint, paired = TRUE,
        p.adjust.method = "fdr"
    )
#effect size for Time Point
SLF_data%>%cohens_d(mn_FDC_SLF2_L~Timepoint, paired=TRUE)

#follow up for sig. interaction 
#Right SLF 1
#interaction post hoc follow-up
#view interaction plot
SLF_data %>%
    group_by(Group,Timepoint) %>%
    summarise(s_mean=mean(mn_FD_SLF1_R)) %>%
    ggplot(aes(y=s_mean,x=Group,colour=Timepoint,group=Timepoint))+
    geom_point()+geom_line()+
    ylab("SLF1 Right Fibre Density (FD)") +
    theme_classic()
#Sig. interaction
#Simple main effect w/ interaction - for Group
posthoc_ME_Group_SLF1_R_FD<- SLF_data %>%
    group_by(Timepoint) %>%
    anova_test(dv=mn_FD_SLF1_R,wid=Individual_number,between=Group) %>%
    adjust_pvalue(method="fdr")
posthoc_ME_Group_SLF1_R_FD
#Pairwise comparison between groups levels
posthoc_pairwise_Group_SLF1_R_FD <- SLF_data %>%
    group_by(Timepoint) %>%
    pairwise_t_test(mn_FD_SLF1_R ~ Group, p.adjust.method = "fdr")
posthoc_pairwise_Group_SLF1_R_FD
#Simple main effect w/ ineraction - for Timepoint
posthoc_ME_Timepoint_SLF1_R_FD <- SLF_data %>%
    group_by(Group) %>%
    anova_test(dv=mn_FD_SLF1_R,wid=Individual_number,within=Timepoint,effect.size = "pes") %>%
    get_anova_table() %>%
    adjust_pvalue(method="fdr")
posthoc_ME_Timepoint_SLF1_R_FD
#Pairwise comparison between groups levels
posthoc_pairwise_Timepoint_SLF1_R_FD <- SLF_data %>%
    group_by(Group) %>%
    pairwise_t_test(
    mn_FD_SLF1_R ~ Timepoint, paired = TRUE, 
    p.adjust.method = "fdr") %>%
    select(-df, -statistic, -p) # Remove details
posthoc_pairwise_Timepoint_SLF1_R_FD
#effect size for interaction (SCD Group and Time Point)
SLF_data_SCD_SLF1_R_FD_long <- subset(SLF_data, DPRC_neuropsych_data$Group == 2)
SLF_data_SCD_SLF1_R_FD_long$Group <- droplevels(SLF_data_SCD_SLF1_R_FD_long$Group)
SLF_data_SCD_SLF1_R_FD_long%>%cohens_d(mn_FD_SLF1_R~Timepoint, paired=TRUE)

#CI for Interaction
#Right SLF 1 FD
post_hoc_aov_SLF1_R_FD_mod <- lme(mn_FD_SLF1_R ~ Group*Timepoint, random = ~1 | ParticipantID/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_R_FD_mod, linfct=mcp(Timepoint="Tukey")))
nlme::intervals(post_hoc_aov_SLF1_R_FD_mod, level=0.95)

#plot 95% Confidence Interval (FD)
SLF_FD_95CI_data <- data.frame(SLF_group_number = c('1','2'),
                               SLF_type = c('Left_SLF2', 'Right_SLF1'),
                               Group_contrast = c('F0vF2', 'F0vF2 in SCD Group'),
                               estimate_diff = c(-0.006490455, 0.022892614),
                               lower = c(-0.02402032, 0.002211254), 
                               upper = c(0.011039408, 0.0435739736))  
#plot data
ggplot(SLF_FD_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Left SLF 2", "2" = "Right SLF 1"))+
    theme_classic() 

#plot 95% Confidence Interval (FDC)
SLF_FDC_95CI_data <- data.frame(SLF_group_number = c('1'),
                                SLF_type = c('Left_SLF2'),
                                Group_contrast = c('F0vF2'),
                                estimate_diff = c(-0.008545409),
                                lower = c(-0.04115788), 
                                upper = c(0.02406706))  
#plot data
ggplot(SLF_FDC_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Left SLF 2"))+
    theme_classic() 

#Run ANCOVA With covariates (age & sex)
#whole SLF
#does not work w/ covariate, Age
#SLF_data$Age <- as.factor(SLF_data$Age)
#aov_SLF_FD_2covar_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF, wid=Individual_number, between=Group, within=Timepoint, covariate=c(Sex,Age), effect.size ="pes")
#get_anova_table(aov_SLF_FD_2covar_mod) 
aov_SLF_FD_2covar_mod<- aov(mn_FD_SLF ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
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
summary(aov_SLF2_L_FD_2covar_mod) #Sig. on Timepoint
aov_SLF2_L_FC_2covar_mod<- aov(mn_FC_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_L_FC_2covar_mod)
aov_SLF2_L_FDC_2covar_mod<- aov(mn_FDC_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_L_FDC_2covar_mod) #Sig. on Timepoint
#Left SLF 3
aov_SLF3_L_FD_2covar_mod<- aov(mn_FD_SLF3_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF3_L_FD_2covar_mod)
aov_SLF3_L_FC_2covar_mod<- aov(mn_FC_SLF3_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF3_L_FC_2covar_mod)
aov_SLF3_L_FDC_2covar_mod<- aov(mn_FDC_SLF3_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF3_L_FDC_2covar_mod)
#Right SLF 1
aov_SLF1_R_FD_2covar_mod<- aov(mn_FD_SLF1_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF1_R_FD_2covar_mod) #Sig. interaction (Group x Timepoint)
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
eta_squared(aov_SLF_FD_2covar_mod)
eta_squared(aov_SLF_FC_2covar_mod)
eta_squared(aov_SLF_FDC_2covar_mod)
eta_squared(aov_SLF_L_FD_2covar_mod)
eta_squared(aov_SLF_L_FC_2covar_mod)
eta_squared(aov_SLF_L_FDC_2covar_mod)
eta_squared(aov_SLF_R_FD_2covar_mod)
eta_squared(aov_SLF_R_FC_2covar_mod)
eta_squared(aov_SLF_R_FDC_2covar_mod)
eta_squared(aov_SLF1_L_FD_2covar_mod)
eta_squared(aov_SLF1_L_FC_2covar_mod)
eta_squared(aov_SLF1_L_FDC_2covar_mod)
eta_squared(aov_SLF2_L_FD_2covar_mod)
eta_squared(aov_SLF2_L_FC_2covar_mod)
eta_squared(aov_SLF2_L_FDC_2covar_mod)
eta_squared(aov_SLF3_L_FD_2covar_mod)
eta_squared(aov_SLF3_L_FC_2covar_mod)
eta_squared(aov_SLF3_L_FDC_2covar_mod)
eta_squared(aov_SLF1_R_FD_2covar_mod)
eta_squared(aov_SLF1_R_FC_2covar_mod)
eta_squared(aov_SLF1_R_FDC_2covar_mod)
eta_squared(aov_SLF2_R_FD_2covar_mod)
eta_squared(aov_SLF2_R_FC_2covar_mod)
eta_squared(aov_SLF2_R_FDC_2covar_mod)
eta_squared(aov_SLF3_R_FD_2covar_mod)
eta_squared(aov_SLF3_R_FC_2covar_mod)
eta_squared(aov_SLF3_R_FDC_2covar_mod)

#follow up tests with sig. main effects and interactions
#Left SLF 2 FD
post_hoc_aov_SLF2_L_FD_2covar_mod <- lme(mn_FD_SLF2_L ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_L_FD_2covar_mod, linfct=mcp(Timepoint="Tukey")))
nlme::intervals(post_hoc_aov_SLF2_L_FD_2covar_mod, level=0.95)
summary(glht(post_hoc_aov_SLF2_L_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#Left SLF 2 FDC
post_hoc_aov_SLF2_L_FDC_2covar_mod <- lme(mn_FDC_SLF2_L ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_L_FDC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
nlme::intervals(post_hoc_aov_SLF2_L_FDC_2covar_mod, level=0.95)
summary(glht(post_hoc_aov_SLF2_L_FDC_2covar_mod, linfct=mcp(Group="Tukey")))
#follow up for sig. interaction 
#Right SLF 1
#interaction post hoc follow-up
#view interaction plot
SLF_data %>%
    group_by(Group,Timepoint) %>%
    summarise(s_mean=mean(mn_FD_SLF1_R)) %>%
    ggplot(aes(y=s_mean,x=Group,colour=Timepoint,group=Timepoint))+
    geom_point()+geom_line()+
    ylab("SLF1 Right Fibre Density (FD)") +
    theme_classic()
#view interaction on the variable for 'Group'
SLF_data %>% filter(Group=="1") %>%
    aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
    summary()
SLF_data %>% filter(Group=="2") %>% #significant interaction for SCD group
    aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
    summary()
SLF_data %>% filter(Group=="3") %>% 
    aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
    summary()
SLF_data %>% filter(Group=="4") %>% 
    aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
    summary()
SLF_data %>% filter(Group=="5") %>% 
    aov(mn_FD_SLF1_R~Timepoint+Sex+Age,data=.)%>%
    summary()
#view interaction on the variable for 'Timepoint'
SLF_data %>% filter(Timepoint=="F0") %>%
    aov(mn_FD_SLF1_R~Group+Sex+Age,data=.)%>%
    summary()
SLF_data %>% filter(Timepoint=="F2") %>%
    aov(mn_FD_SLF1_R~Group+Sex+Age,data=.)%>%
    summary()

#CI for Interaction w/ covariates
#Right SLF 1 FD
post_hoc_aov_SLF1_R_FD_2covar_mod <- lme(mn_FD_SLF1_R ~ Group*Timepoint + Age + Sex, random = ~1 | ParticipantID/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_R_FD_2covar_mod, linfct=mcp(Timepoint="Tukey")))
nlme::intervals(post_hoc_aov_SLF1_R_FD_2covar_mod, level=0.95) #doesn't work w/ covariates in the model.

#plot 95% Confidence Interval w/ covariates (FD)
SLF_FD_95CI_2covar_data <- data.frame(SLF_group_number = c('1'),
                               SLF_type = c('Left_SLF2'),
                               Group_contrast = c('F0vF2'),
                               estimate_diff = c(-0.0058891960),
                               lower = c(-0.0234024688), 
                               upper = c(0.0116240768))  
#plot data
ggplot(SLF_FD_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Left SLF 2"))+
    theme_classic() 

#plot 95% Confidence Interval w/ covariates (FDC)
SLF_FDC_95CI_2covar_data <- data.frame(SLF_group_number = c('1'),
                                SLF_type = c('Left_SLF2'),
                                Group_contrast = c('F0vF2'),
                                estimate_diff = c(-0.0077689181),
                                lower = c(-0.040474001), 
                                upper = c(0.0249361650))  
#plot data
ggplot(SLF_FDC_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Left SLF 2"))+
    theme_classic() 

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
#All tracts SLF FD (raincloud plot)
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


#Run linear regression comparison between SLF metrics and neuropsych test scores
#prep dataframe
SLF_data_F0 <- SLF_data[1:124,]
SLF_data_F2 <- SLF_data[125:248,]
SLF_data_subtracted <- SLF_data_F2 [, c(1:189, 197:220)] - SLF_data_F0[, c(1:189, 197:220)]
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

#Simple linear regression test-------------------------------------------------#
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
summary(SLF_L_FC_Hay1Raw_mod) #sig
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





