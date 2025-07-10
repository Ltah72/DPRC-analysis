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
pacman::p_load(dplyr, ggplot2, psych, car, multcomp, lsr, tidyr, BayesFactor, tidyverse, ppcor, nlme, effectsize, rstatix, sjstats, purrr, corrplot, compute.es, writexl)

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
setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/template/TOI')

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
                        "count_FD_SLF3_L", "mn_FC_SLF3_L", "md_FC_SLF3_L", "std_FC_SLF3_L", "std_rv_FC_SLF3_L", "min_FC_SLF3_L", 
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
SLF_FDC_descrip <- describeBy(SLF_data$mn_FDC_SLF, SLF_data$Group)
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



###------------------ ANOVA & ANCOVA (age) + (age + sex) ---------------------------####

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
SLF1_FD_covar_mod_L <- lm(mn_FD_SLF1_L ~ Group + Age, data = SLF_data)
SLF2_FD_covar_mod_L <- lm(mn_FD_SLF2_L ~ Group + Age, data = SLF_data)
SLF3_FD_covar_mod_L <- lm(mn_FD_SLF3_L ~ Group + Age, data = SLF_data)
SLF1_FD_covar_mod_R <- lm(mn_FD_SLF1_R ~ Group + Age, data = SLF_data)
SLF2_FD_covar_mod_R <- lm(mn_FD_SLF2_R ~ Group + Age, data = SLF_data)
SLF3_FD_covar_mod_R <- lm(mn_FD_SLF3_R ~ Group + Age, data = SLF_data)

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


#test to see if the covariates (sex & age) and the treatment variable are independent from one another (assumption). 
#for FD:
  #whole FD
  # sex_FD_check <- aov(mn_FD_SLF ~ Sex+Group, data = SLF_data)
  # summary(sex_FD_check)
  sex_whole_FD_check <- lm(mn_FD_SLF ~ Sex+Group, data = SLF_data)
  anova(sex_whole_FD_check)
  # t.test(SLF_data$mn_FD_SLF ~ SLF_data$Sex, var.equal = FALSE) #you wouldn't run t-test if you want to include a covariate (e.g., Group) into your model - run ANCOVA
  age_whole_FD_check <- lm(mn_FD_SLF ~ Age+Group, data = SLF_data)
  anova(age_whole_FD_check)
  #left SLF
  #sex
  sex_L_FD_check <- lm(mn_FD_SLF_L ~ Sex+Group, data = SLF_data)
  anova(sex_L_FD_check)
  #age
  age_L_FD_check <- lm(mn_FD_SLF_L ~ Age+Group, data = SLF_data)
  anova(age_L_FD_check)
  #right SLF
  #sex
  sex_R_FD_check <- lm(mn_FD_SLF_R ~ Sex+Group, data = SLF_data)
  anova(sex_R_FD_check)
  #age
  age_R_FD_check <- lm(mn_FD_SLF_R ~ Age+Group, data = SLF_data)
  anova(age_R_FD_check)
  #left SLF 1
  #sex
  sex_L1_FD_check <- lm(mn_FD_SLF1_L ~ Sex+Group, data = SLF_data)
  anova(sex_L1_FD_check) #not sig. - no sex diff
  #age
  age_L1_FD_check <- lm(mn_FD_SLF1_L ~ Age+Group, data = SLF_data)
  anova(age_L1_FD_check) #not sig. - no age diff
  #left SLF 2
  #sex
  sex_L2_FD_check <- lm(mn_FD_SLF2_L ~ Sex+Group, data = SLF_data)
  anova(sex_L2_FD_check) 
  #age
  age_L2_FD_check <- lm(mn_FD_SLF2_L ~ Age+Group, data = SLF_data)
  anova(age_L2_FD_check) 
  #left SLF 3
  #sex
  sex_L3_FD_check <- lm(mn_FD_SLF3_L ~ Sex+Group, data = SLF_data)
  anova(sex_L3_FD_check) 
  #age
  age_L3_FD_check <- lm(mn_FD_SLF3_L ~ Age+Group, data = SLF_data)
  anova(age_L3_FD_check) 
  #right SLF 1
  #sex
  sex_R1_FD_check <- lm(mn_FD_SLF1_R ~ Sex+Group, data = SLF_data)
  anova(sex_R1_FD_check) 
  #age
  age_R1_FD_check <- lm(mn_FD_SLF1_R ~ Age+Group, data = SLF_data)
  anova(age_R1_FD_check)
  #right SLF 2
  #sex
  sex_R2_FD_check <- lm(mn_FD_SLF2_R ~ Sex+Group, data = SLF_data)
  anova(sex_R2_FD_check) 
  #age
  age_R2_FD_check <- lm(mn_FD_SLF2_R ~ Age+Group, data = SLF_data)
  anova(age_R2_FD_check) #not sig. - no age diff
  #right SLF 3
  #sex
  sex_R3_FD_check <- lm(mn_FD_SLF3_R ~ Sex+Group, data = SLF_data)
  anova(sex_R3_FD_check) 
  #age
  age_R3_FD_check <- lm(mn_FD_SLF3_R ~ Age+Group, data = SLF_data)
  anova(age_R3_FD_check) #not sig. - no age diff
  
#for FC:
  #whole FC
  #sex
  sex_whole_FC_check <- lm(mn_FC_SLF ~ Sex+Group, data = SLF_data)
  anova(sex_whole_FC_check)
  #age
  age_whole_FC_check <- lm(mn_FC_SLF ~ Age+Group, data = SLF_data)
  anova(age_whole_FC_check) #not sig. - no age diff
  #left SLF
  #sex
  sex_L_FC_check <- lm(mn_FC_SLF_L ~ Sex+Group, data = SLF_data)
  anova(sex_L_FC_check)
  #age
  age_L_FC_check <- lm(mn_FC_SLF_L ~ Age+Group, data = SLF_data)
  anova(age_L_FC_check) #not sig. - no age diff
  #right SLF
  #sex
  sex_R_FC_check <- lm(mn_FC_SLF_R ~ Sex+Group, data = SLF_data)
  anova(sex_R_FC_check) 
  #age
  age_R_FC_check <- lm(mn_FC_SLF_R ~ Age+Group, data = SLF_data)
  anova(age_R_FC_check) #not sig. - no age diff
  #left SLF 1
  #sex
  sex_L1_FC_check <- lm(mn_FC_SLF1_L ~ Sex+Group, data = SLF_data)
  anova(sex_L1_FC_check)
  #age
  age_L1_FC_check <- lm(mn_FC_SLF1_L ~ Age+Group, data = SLF_data)
  anova(age_L1_FC_check) 
  #left SLF 2
  #sex
  sex_L2_FC_check <- lm(mn_FC_SLF2_L ~ Sex+Group, data = SLF_data)
  anova(sex_L2_FC_check) 
  #age
  age_L2_FC_check <- lm(mn_FC_SLF2_L ~ Age+Group, data = SLF_data)
  anova(age_L2_FC_check) #not sig. - no age diff
  #left SLF 3
  #sex
  sex_L3_FC_check <- lm(mn_FC_SLF3_L ~ Sex+Group, data = SLF_data)
  anova(sex_L3_FC_check) 
  #age
  age_L3_FC_check <- lm(mn_FC_SLF3_L ~ Age+Group, data = SLF_data)
  anova(age_L3_FC_check) #not sig. - no age diff
  #right SLF 1
  #sex
  sex_R1_FC_check <- lm(mn_FC_SLF1_R ~ Sex+Group, data = SLF_data)
  anova(sex_R1_FC_check) 
  #age
  age_R1_FC_check <- lm(mn_FC_SLF1_R ~ Age+Group, data = SLF_data)
  anova(age_R1_FC_check)
  #right SLF 2
  #sex
  sex_R2_FC_check <- lm(mn_FC_SLF2_R ~ Sex+Group, data = SLF_data)
  anova(sex_R2_FC_check) 
  #age
  age_R2_FC_check <- lm(mn_FC_SLF2_R ~ Age+Group, data = SLF_data)
  anova(age_R2_FC_check) #not sig. - no age diff
  #right SLF 3
  #sex
  sex_R3_FC_check <- lm(mn_FC_SLF3_R ~ Sex+Group, data = SLF_data)
  anova(sex_R3_FC_check) 
  #age
  age_R3_FC_check <- lm(mn_FC_SLF3_R ~ Age+Group, data = SLF_data)
  anova(age_R3_FC_check) #not sig. - no age diff

#for FDC:
  #whole FDC
  #sex
  sex_whole_FDC_check <- lm(mn_FDC_SLF ~ Sex+Group, data = SLF_data)
  anova(sex_whole_FDC_check) #not sig. - no sex diff
  #age
  age_whole_FDC_check <- lm(mn_FDC_SLF ~ Age+Group, data = SLF_data)
  anova(age_whole_FDC_check)
  #left SLF
  #sex
  sex_L_FDC_check <- lm(mn_FDC_SLF_L ~ Sex+Group, data = SLF_data)
  anova(sex_L_FDC_check) #not sig. - no sex diff
  #age
  age_L_FDC_check <- lm(mn_FDC_SLF_L ~ Age+Group, data = SLF_data)
  anova(age_L_FDC_check) 
  #right SLF
  #sex
  sex_R_FDC_check <- lm(mn_FDC_SLF_R ~ Sex+Group, data = SLF_data)
  anova(sex_R_FDC_check) #not sig. - no sex diff
  #age
  age_R_FDC_check <- lm(mn_FDC_SLF_R ~ Age+Group, data = SLF_data)
  anova(age_R_FDC_check)
  #left SLF 1
  #sex
  sex_L1_FDC_check <- lm(mn_FDC_SLF1_L ~ Sex+Group, data = SLF_data)
  anova(sex_L1_FDC_check) 
  #age
  age_L1_FDC_check <- lm(mn_FDC_SLF1_L ~ Age+Group, data = SLF_data)
  anova(age_L1_FDC_check) #not sig. - no age diff
  #left SLF 2
  #sex
  sex_L2_FDC_check <- lm(mn_FDC_SLF2_L ~ Sex+Group, data = SLF_data)
  anova(sex_L2_FDC_check) #not sig. - no sex diff
  #age
  age_L2_FDC_check <- lm(mn_FDC_SLF2_L ~ Age+Group, data = SLF_data)
  anova(age_L2_FDC_check)
  #left SLF 3
  #sex
  sex_L3_FDC_check <- lm(mn_FDC_SLF3_L ~ Sex+Group, data = SLF_data)
  anova(sex_L3_FDC_check) #not sig. - no sex diff
  #age
  age_L3_FDC_check <- lm(mn_FDC_SLF3_L ~ Age+Group, data = SLF_data)
  anova(age_L3_FDC_check)
  #right SLF 1
  #sex
  sex_R1_FDC_check <- lm(mn_FDC_SLF1_R ~ Sex+Group, data = SLF_data)
  anova(sex_R1_FDC_check) #not sig. - no sex diff
  #age
  age_R1_FDC_check <- lm(mn_FDC_SLF1_R ~ Age+Group, data = SLF_data)
  anova(age_R1_FDC_check)
  #right SLF 2
  #sex
  sex_R2_FDC_check <- lm(mn_FDC_SLF2_R ~ Sex+Group, data = SLF_data)
  anova(sex_R2_FDC_check) 
  #age
  age_R2_FDC_check <- lm(mn_FDC_SLF2_R ~ Age+Group, data = SLF_data)
  anova(age_R2_FDC_check) 
  #right SLF 3
  #sex
  sex_R3_FDC_check <- lm(mn_FDC_SLF3_R ~ Sex+Group, data = SLF_data)
  anova(sex_R3_FDC_check) #not sig. - no sex diff
  #age
  age_R3_FDC_check <- lm(mn_FDC_SLF3_R ~ Age+Group, data = SLF_data)
  anova(age_R3_FDC_check) 
  
# #means examples
#   SLF3_FD_L_descrip <- describeBy(SLF_data$mn_FD_SLF3_R, SLF_data$Sex)
#   SLF_FD_whole_descrip <- describeBy(SLF_data$mn_FD_SLF, SLF_data$Sex)
#   
# #plotting example
#   #whole SLF FD (raincloud plot)
#   ggplot(SLF_data, aes(x = Group, y = mn_FD_SLF, fill = Sex)) + 
#     geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_point(aes(y = mn_FD_SLF, color = Sex), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Sex)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Sex)) + 
#     xlab("Sex") + 
#     ylab("Fibre Density Cross-section (FD)") +
#     scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
#     theme_classic() +
#     #theme(legend.position = "none") +
#     coord_flip()
#   


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
#run ANCOVA (age)
anova(SLF1_FD_covar_mod_L)
anova(SLF2_FD_covar_mod_L)
anova(SLF3_FD_covar_mod_L)
anova(SLF1_FD_covar_mod_R)
anova(SLF2_FD_covar_mod_R)
anova(SLF3_FD_covar_mod_R)

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

#Linear Trend Analysis in Linear Regression w/ covariate age
SLF1_L_FD_LinTrend_covar_mod <- lm(mn_FD_SLF1_L ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF1_L_FD_LinTrend_covar_mod)
summary(SLF1_L_FD_LinTrend_covar_mod) 
SLF2_L_FD_LinTrend_covar_mod <- lm(mn_FD_SLF2_L ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF2_L_FD_LinTrend_covar_mod)
summary(SLF2_L_FD_LinTrend_covar_mod) 
SLF3_L_FD_LinTrend_covar_mod <- lm(mn_FD_SLF3_L ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF3_L_FD_LinTrend_covar_mod)
summary(SLF3_L_FD_LinTrend_covar_mod)
SLF1_R_FD_LinTrend_covar_mod <- lm(mn_FD_SLF1_R ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF1_R_FD_LinTrend_covar_mod)
summary(SLF1_R_FD_LinTrend_covar_mod) 
SLF2_R_FD_LinTrend_covar_mod <- lm(mn_FD_SLF2_R ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF2_R_FD_LinTrend_covar_mod)
summary(SLF2_R_FD_LinTrend_covar_mod) 
SLF3_R_FD_LinTrend_covar_mod <- lm(mn_FD_SLF3_R ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF3_R_FD_LinTrend_covar_mod)
summary(SLF3_R_FD_LinTrend_covar_mod)

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

#post hoc for covariate (age)
#Left SLF 1, 2, and 3
post_hoc_SLF1_FD_covar_mod_L <- glht(SLF1_FD_covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF1_FD_covar_mod_L)
confint(post_hoc_SLF1_FD_covar_mod_L)
post_hoc_SLF2_FD_covar_mod_L <- glht(SLF2_FD_covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_FD_covar_mod_L)
confint(post_hoc_SLF2_FD_covar_mod_L)
post_hoc_SLF3_FD_covar_mod_L <- glht(SLF3_FD_covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_FD_covar_mod_L)
confint(post_hoc_SLF3_FD_covar_mod_L)
#Right SLF 1, 2 and 3
post_hoc_SLF1_FD_covar_mod_R <- glht(SLF1_FD_covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF1_FD_covar_mod_R)
confint(post_hoc_SLF1_FD_covar_mod_R)
post_hoc_SLF2_FD_covar_mod_R <- glht(SLF2_FD_covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_FD_covar_mod_R)
confint(post_hoc_SLF2_FD_covar_mod_R)
post_hoc_SLF3_FD_covar_mod_R <- glht(SLF3_FD_covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_FD_covar_mod_R)
confint(post_hoc_SLF3_FD_covar_mod_R)

#post hoc for covariate (age & sex)
#Whole, Left, and Right SLF
post_hoc_SLF_FD_ancova_mod <- glht(SLF_FD_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_FD_ancova_mod)
confint(post_hoc_SLF_FD_ancova_mod)
post_hoc_SLF_L_FD_ancova_mod <- glht(SLF_FD_2covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_L_FD_ancova_mod)
confint(post_hoc_SLF_L_FD_ancova_mod)
post_hoc_SLF_R_FD_ancova_mod <- glht(SLF_FD_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF_R_FD_ancova_mod)
confint(post_hoc_SLF_R_FD_ancova_mod)
#Left SLF 1, 2, and 3
post_hoc_SLF2_L_FD_ancova_mod <- glht(SLF2_FD_2covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_L_FD_ancova_mod) #not sig. 
confint(post_hoc_SLF2_L_FD_ancova_mod)
post_hoc_SLF3_L_FD_ancova_mod <- glht(SLF3_FD_2covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_L_FD_ancova_mod)
confint(post_hoc_SLF3_L_FD_ancova_mod)
#Right SLF 1, 2, and 3
post_hoc_SLF1_R_FD_ancova_mod <- glht(SLF1_FD_2covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF1_R_FD_ancova_mod) #not sig.
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

#plot 95% Confidence Interval (separated confidence intervals and separated SLF tracts only)
SLF_FD_95CI_data_no_combined_SLF <- data.frame(SLF_group_number = c('1','1','1','2','2','2','2','3','4','4','4','5','5','5','5'),
                               SLF_type = c('Left_SLF2','Left_SLF2','Left_SLF2','Left_SLF3','Left_SLF3','Left_SLF3','Left_SLF3','Right_SLF1','Right_SLF2','Right_SLF2','Right_SLF2','Right_SLF3', 'Right_SLF3', 'Right_SLF3', 'Right_SLF3'),
                               Group_contrast = c('C > mMCI','C > AD','SCD > AD','C > SCD','C > aMCI','C > mMCI','C > AD','C > AD','C > aMCI','C > mMCI','C > AD','C > SCD','C > aMCI','C > mMCI','C > AD'),
                               estimate_diff = c(0.020646,0.032745,0.022266,0.0149789,0.0155556,0.0198357,0.0251135,0.0282489,0.019994,0.0245799,0.0309866,0.0155239,0.0163904,0.0214478,0.0194152),
                               lower = c(0.040536,0.056048,0.043348,0.0280131,0.0288062,0.0332341,0.0408103,0.0490901,0.0394673,0.0442704,0.0540548,0.0293862,0.0304829,0.0356975,0.0361093), 
                               upper = c(0.000755,0.009442,0.001183,0.0019447,0.0023049,0.0064372,0.0094167,0.0074076,0.0005208,0.0048893,0.0079185,0.0016616,0.0022978,0.0071981,0.0027211))  

#plot data
ggplot(SLF_FD_95CI_data_no_combined_SLF, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
  geom_point(position = position_dodge(width=0.5), size=5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), size=1.5, position = position_dodge(width=0.5))+
  xlab("SLF Tract") + 
  ylab("95% Confidence Interval")+
  scale_x_discrete(labels = c("1" = "Left SLF 2", "2" = "Left SLF 3", "3" = "Right SLF 1", "4" = "Right SLF 2", "5" = "Right SLF 3"))+
  #scale_x_discrete(labels = c("1" = "Whole SLF", "2" = "Left SLF", "3" = "Right SLF", "4" = "Right SLF 1", "5" = "Left SLF 2", "6" = "Left SLF 3", "7" = "Right SLF 2", "8" = "Right SLF 3"))+
  scale_color_manual(values = c("#FF00FF", "#33CC66", "#33CCFF", "#999900", "#CCCCCC"))+
  labs(colour='Group Contrast')+   
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()


    
#plot 95% Confidence Interval for ANCOVA (age & sex)
#for FDC
#Create dataset: 
SLF_FDC_95CI_2covar_data <- data.frame(SLF_group_number = c('1','1','3','3','5','5','5','7','8','8'),
                               SLF_type = c('Whole_SLF','Whole_SLF', 'Right_SLF', 'Right_SLF','Left_SLF3','Left_SLF3','Left_SLF3','Right_SLF2','Right_SLF3','Right_SLF3'),
                               Group_contrast = c('CvmMCI', 'CvAD', 'CvmMCI', 'CvAD','CvSCD','CvmMCI','CvAD','CvmMCI','CvSCD','CvmMCI'),
                               estimate_diff = c(0.0149701,0.0171247,0.016138,0.017475,0.013165,0.015781,0.0174047,0.0201535,0.0142024,0.019025),
                               lower = c(0.0286742,0.0336528,0.030519,0.03482,0.0257797,0.0289482,0.0332852,0.0397022,0.02792,0.0333434), 
                               upper = c(0.0012659,0.0005967,0.001757,0.000131,0.0005503,0.0026138,0.0015243,0.0006048,0.0004848,0.0047066))  
#plot data
ggplot(SLF_FDC_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Whole SLF","3" = "Right SLF","5" = "Left SLF 3", "7" = "Right SLF 2", "8" = "Right SLF 3"))+
    theme_classic() 
    

#plot 95% Confidence Interval (separated confidence intervals and separated SLF tracts only) (sex & age covariates)
SLF_FDC_95CI_data_no_combined_SLF <- data.frame(SLF_group_number = c('1','1','1','2','3','3'),
                                               SLF_type = c('Left_SLF3','Left_SLF3','Left_SLF3','Right_SLF2','Right_SLF3', 'Right_SLF3'),
                                               Group_contrast = c('C > SCD','C > mMCI','C > AD','C > mMCI','C > SCD','C > mMCI'),
                                               estimate_diff = c(0.0130476,0.0156277,0.0173008,0.0198848,0.0141504,0.0189202),
                                               lower = c(0.0256490,0.0287923,0.0331592,0.0394443,0.0278481,0.0332301), 
                                               upper = c(0.0004462,0.0024631,0.0014424,0.0003253,0.0004526,0.0046103))  

#plot data
ggplot(SLF_FDC_95CI_data_no_combined_SLF, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
  geom_point(position = position_dodge(width=0.5), size=5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), size=1.5, position = position_dodge(width=0.5))+
  xlab("FDC") + 
  ylab("95% Confidence Interval")+
  scale_x_discrete(labels = c("1" = "Left SLF 3", "2" = "Right SLF 2", "3" = "Right SLF 3"))+
  scale_color_manual(values = c("#FF00FF","#33CCFF","#999900"))+
  labs(colour='Group Contrast')+   
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()


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

#for ANCOVA (age)
#calculate the effect size (eta-squared)
etaSquared(SLF1_FD_covar_mod_L)
etaSquared(SLF2_FD_covar_mod_L)
etaSquared(SLF3_FD_covar_mod_L)
etaSquared(SLF1_FD_covar_mod_R)
etaSquared(SLF2_FD_covar_mod_R)
etaSquared(SLF3_FD_covar_mod_R)

#for ANCOVA (sex & age)
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
    

#For ANCOVA (age) - effect size for sig. post hoc tests (Cohen's d) (using a.tes function from the compute.es package)
#Left SLF2 FD
t_value_effect_size <- summary(post_hoc_SLF2_FD_covar_mod_L) 
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_SLF2_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
  DPRC_neuropsych_data_CvAD_SLF2_L$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF2_L$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_SLF2_L, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF2_L ~ Age, data = DPRC_neuropsych_data_CvAD_SLF2_L)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age   
#Left SLF3 FD
t_value_effect_size <- summary(post_hoc_SLF3_FD_covar_mod_L) 
    #for C vs. SCD
    DPRC_neuropsych_data_CvSCD_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
    DPRC_neuropsych_data_CvSCD_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_L$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_SLF3_L, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age, data = DPRC_neuropsych_data_CvSCD_SLF3_L)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
    #for C vs. aMCI 
    DPRC_neuropsych_data_CvaMCI_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 3)
    DPRC_neuropsych_data_CvaMCI_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_SLF3_L$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_SLF3_L, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age, data = DPRC_neuropsych_data_CvaMCI_SLF3_L)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_L$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF3_L, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age, data = DPRC_neuropsych_data_CvmMCI_SLF3_L)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF3_L <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF3_L$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_SLF3_L, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age, data = DPRC_neuropsych_data_CvAD_SLF3_L)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age   
#Right SLF1 FD
t_value_effect_size <- summary(post_hoc_SLF1_FD_covar_mod_R) 
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF1_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF1_R$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF1_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_SLF1_R, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF1_R ~ Age, data = DPRC_neuropsych_data_CvAD_SLF1_R)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age   
#Right SLF2 FD
t_value_effect_size <- summary(post_hoc_SLF2_FD_covar_mod_R) 
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF2_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF2_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF2_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF2_R, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF2_R ~ Age, data = DPRC_neuropsych_data_CvmMCI_SLF2_R)) #create multiple regression between age, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF2_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF2_R$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF2_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_SLF2_R, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF2_R ~ Age, data = DPRC_neuropsych_data_CvAD_SLF2_R)) #create multiple regression between age, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
#Right SLF3 FD
t_value_effect_size <- summary(post_hoc_SLF3_FD_covar_mod_R) 
    #for C vs. SCD
    DPRC_neuropsych_data_CvSCD_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
    DPRC_neuropsych_data_CvSCD_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_SLF3_R, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_R ~ Age, data = DPRC_neuropsych_data_CvSCD_SLF3_R)) #create multiple regression between age, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
    #for C vs. aMCI
    DPRC_neuropsych_data_CvaMCI_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 3)
    DPRC_neuropsych_data_CvaMCI_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_SLF3_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_SLF3_R, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_R ~ Age, data = DPRC_neuropsych_data_CvaMCI_SLF3_R)) #create multiple regression between age, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
    #for C vs. mMCI 
    DPRC_neuropsych_data_CvmMCI_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 4)
    DPRC_neuropsych_data_CvmMCI_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF3_R, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_R ~ Age, data = DPRC_neuropsych_data_CvmMCI_SLF3_R)) #create multiple regression between age, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
    #for C vs. AD
    DPRC_neuropsych_data_CvAD_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
    DPRC_neuropsych_data_CvAD_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvAD_SLF3_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_SLF3_R, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_R ~ Age, data = DPRC_neuropsych_data_CvAD_SLF3_R)) #create multiple regression between age, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age

#For ANCOVA (age & sex) - effect size for sig. post hoc tests (Cohen's d) (using a.tes function from the compute.es package)
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
    group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_SLF3_L, Group) #count number of participants per group
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

# #include the covariate of clinical site (run an ANCOVA) in model
# SLF_FC_clinsite_mod <- lm(mn_FC_SLF ~ Group + ClinSite_covar, data = SLF_data)
# SLF_FC_clinsite_mod_L <- lm(mn_FC_SLF_L ~ Group + ClinSite_covar, data = SLF_data)
# SLF_FC_clinsite_mod_R <- lm(mn_FC_SLF_R ~ Group + ClinSite_covar, data = SLF_data)

#include the covariate (age) (run an ANCOVA) in model
SLF1_FC_covar_mod_L <- lm(mn_FC_SLF1_L ~ Group + Age, data = SLF_data)
SLF2_FC_covar_mod_L <- lm(mn_FC_SLF2_L ~ Group + Age, data = SLF_data)
SLF3_FC_covar_mod_L <- lm(mn_FC_SLF3_L ~ Group + Age, data = SLF_data)
SLF1_FC_covar_mod_R <- lm(mn_FC_SLF1_R ~ Group + Age, data = SLF_data)
SLF2_FC_covar_mod_R <- lm(mn_FC_SLF2_R ~ Group + Age, data = SLF_data)
SLF3_FC_covar_mod_R <- lm(mn_FC_SLF3_R ~ Group + Age, data = SLF_data)

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
# #for ANCOVA
# anova(SLF_FC_clinsite_mod)
# anova(SLF_FC_clinsite_mod_L)
# anova(SLF_FC_clinsite_mod_R)
#ANCOVA - for age 
anova(SLF1_FC_covar_mod_L)
anova(SLF2_FC_covar_mod_L)
anova(SLF3_FC_covar_mod_L)
anova(SLF1_FC_covar_mod_R)
anova(SLF2_FC_covar_mod_R)
anova(SLF3_FC_covar_mod_R)
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

#Linear Trend Analysis in Linear Regression w/ covariates age 
SLF1_L_FC_LinTrend_covar_mod <- lm(mn_FC_SLF1_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF1_L_FC_LinTrend_covar_mod)
summary(SLF1_L_FC_LinTrend_covar_mod) 
SLF2_L_FC_LinTrend_covar_mod <- lm(mn_FC_SLF2_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF2_L_FC_LinTrend_covar_mod)
summary(SLF2_L_FC_LinTrend_covar_mod) 
SLF3_L_FC_LinTrend_covar_mod <- lm(mn_FC_SLF3_L ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF3_L_FC_LinTrend_covar_mod)
summary(SLF3_L_FC_LinTrend_covar_mod)
SLF1_R_FC_LinTrend_covar_mod <- lm(mn_FC_SLF1_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF1_R_FC_LinTrend_covar_mod)
summary(SLF1_R_FC_LinTrend_covar_mod) 
SLF2_R_FC_LinTrend_covar_mod <- lm(mn_FC_SLF2_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF2_R_FC_LinTrend_covar_mod)
summary(SLF2_R_FC_LinTrend_covar_mod) 
SLF3_R_FC_LinTrend_covar_mod <- lm(mn_FC_SLF3_R ~ Trend_Group + Group + Age + Sex, data = SLF_data)
anova(SLF3_R_FC_LinTrend_covar_mod)
summary(SLF3_R_FC_LinTrend_covar_mod) 

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
# #include the covariate of clinical site (run an ANCOVA) in model
# SLF_FDC_clinsite_mod <- lm(mn_FDC_SLF ~ Group + ClinSite_covar, data = SLF_data)
# SLF_FDC_clinsite_mod_L <- lm(mn_FDC_SLF_L ~ Group + ClinSite_covar, data = SLF_data)
# SLF_FDC_clinsite_mod_R <- lm(mn_FDC_SLF_R ~ Group + ClinSite_covar, data = SLF_data)
#include the covariate (age) (run an ANCOVA) in model
SLF1_FDC_covar_mod_L <- lm(mn_FDC_SLF1_L ~ Group + Age, data = SLF_data)
SLF2_FDC_covar_mod_L <- lm(mn_FDC_SLF2_L ~ Group + Age, data = SLF_data)
SLF3_FDC_covar_mod_L <- lm(mn_FDC_SLF3_L ~ Group + Age, data = SLF_data)
SLF1_FDC_covar_mod_R <- lm(mn_FDC_SLF1_R ~ Group + Age, data = SLF_data)
SLF2_FDC_covar_mod_R <- lm(mn_FDC_SLF2_R ~ Group + Age, data = SLF_data)
SLF3_FDC_covar_mod_R <- lm(mn_FDC_SLF3_R ~ Group + Age, data = SLF_data)
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
#ANCOVA - for age
anova(SLF1_FDC_covar_mod_L)
anova(SLF2_FDC_covar_mod_L)
anova(SLF3_FDC_covar_mod_L)
anova(SLF1_FDC_covar_mod_R)
anova(SLF2_FDC_covar_mod_R)
anova(SLF3_FDC_covar_mod_R)
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

#Linear Trend Analysis in Linear Regression w/ covariate age 
SLF1_L_FDC_LinTrend_covar_mod <- lm(mn_FDC_SLF1_L ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF1_L_FDC_LinTrend_covar_mod)
summary(SLF1_L_FDC_LinTrend_covar_mod) 
SLF2_L_FDC_LinTrend_covar_mod <- lm(mn_FDC_SLF2_L ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF2_L_FDC_LinTrend_covar_mod)
summary(SLF2_L_FDC_LinTrend_covar_mod) 
SLF3_L_FDC_LinTrend_covar_mod <- lm(mn_FDC_SLF3_L ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF3_L_FDC_LinTrend_covar_mod)
summary(SLF3_L_FDC_LinTrend_covar_mod)
SLF1_R_FDC_LinTrend_covar_mod <- lm(mn_FDC_SLF1_R ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF1_R_FDC_LinTrend_covar_mod)
summary(SLF1_R_FDC_LinTrend_covar_mod) 
SLF2_R_FDC_LinTrend_covar_mod <- lm(mn_FDC_SLF2_R ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF2_R_FDC_LinTrend_covar_mod)
summary(SLF2_R_FDC_LinTrend_covar_mod) 
SLF3_R_FDC_LinTrend_covar_mod <- lm(mn_FDC_SLF3_R ~ Trend_Group + Group + Age, data = SLF_data)
anova(SLF3_R_FDC_LinTrend_covar_mod)
summary(SLF3_R_FDC_LinTrend_covar_mod) 

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

#post hoc for covariate (age)
post_hoc_SLF2_L_FDC_covar_mod <- glht(SLF2_FDC_covar_mod_L, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_L_FDC_covar_mod) # no sig. f/u
confint(post_hoc_SLF2_L_FDC_covar_mod)
post_hoc_SLF3_R_FDC_covar_mod <- glht(SLF3_FDC_covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_R_FDC_covar_mod)
confint(post_hoc_SLF3_R_FDC_covar_mod)

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
# #plot data
# ggplot(SLF_FDC_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
#     geom_point()+
#     geom_errorbar(aes(ymin=lower, ymax=upper))+
#     xlab("SLF Tract") + 
#     ylab("95% Confidence Interval")+
#     scale_x_discrete(labels = c("1" = "Right SLF 3"))+
#     theme_classic() 

#plot data (just with the right SLF3)
ggplot(SLF_FDC_95CI_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point(position = position_dodge(width=0.5), size=5)+
    geom_errorbar(aes(ymin=lower, ymax=upper), size = 1.5, position = position_dodge(width=0.5))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Right SLF 3"))+
    labs(colour="Group Contrast")+ 
    scale_color_manual(values = "#999900")+
    theme_classic()+
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    coord_flip()




#plot 95% Confidence Interval - with covariates (age + sex)
#for FDC
#Create dataset: 
SLF_FDC_95CI_2covar_data <- data.frame(SLF_group_number = c('1','1'),
                                SLF_type = c('Right_SLF3','Right_SLF3'),
                                Group_contrast = c('C > SCD', 'C > aMCI'),
                                estimate_diff = c(0.0345572317,0.0313322949),
                                lower = c(0.06474609,0.06256263), 
                                upper = c(0.0043683764,0.0001019582)) 
#plot data
ggplot(SLF_FDC_95CI_2covar_data, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point(position = position_dodge(width=0.5), size=5)+
    geom_errorbar(aes(ymin=lower, ymax=upper), size = 1.5, position = position_dodge(width=0.5))+
    xlab("SLF Tract") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Right SLF 3"))+
    scale_color_manual(values = c("#33CC66","#999900"))+
    labs(colour='Group Contrast')+      
    theme_classic()+
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    coord_flip()

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
#for ANCOVA (age)
#calculate the effect size (eta-squared)
etaSquared(SLF1_FDC_covar_mod_L)
etaSquared(SLF2_FDC_covar_mod_L)
etaSquared(SLF3_FDC_covar_mod_L)
etaSquared(SLF1_FDC_covar_mod_R)
etaSquared(SLF2_FDC_covar_mod_R)
etaSquared(SLF3_FDC_covar_mod_R)
#for ANCOVA (age & sex)
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

#effect size for sig. post hoc tests (with age as covariate)
#Left SLF2 FDC - no sig. f/us
#Right SLF3 FDC
t_value_effect_size <- summary(post_hoc_SLF3_R_FDC_covar_mod) 
    #for C vs. SCD 
    DPRC_neuropsych_data_CvSCD_SLF3_R <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
    DPRC_neuropsych_data_CvSCD_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_R$Group)
    group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_SLF3_R, Group) #count number of participants per group
    mult.r_value_covar_mod<-summary(lm(mn_FDC_SLF3_R ~ Age, data = DPRC_neuropsych_data_CvSCD_SLF3_R)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age

#effect size for sig. post hoc tests (with sex + age as covariates)
#Left SLF2 FDC - no sig. f/us
#Right SLF3 FDC
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
    cohensD(mn_FDC_SLF3_R ~ Group, data = DPRC_neuropsych_data_CvmMCI_SLF3_R) #this looks like Hedges' g? 
    group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF3_R, Group) #count number of participants per group
    mult.r_value_2covar_mod<-summary(lm(mn_FDC_SLF3_R ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_SLF3_R)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
    r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
    a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
    
    

#plot data
#For FDC - mean 
#whole SLF FDC (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FDC)") +
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
    ylab("Fibre Density (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#left SLF FDC (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF_L)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FDC)") +
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
    ylab("Fibre Density (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#right SLF FDC (violin plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF_R)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)

#right SLF FDC (raincloud plot)
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF_R, fill = Group)) + 
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FDC_SLF_R, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()


#plot the 6 tracts out separately (for FDC)
#Left SLF 1: 
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF1_L, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mn_FDC_SLF1_L, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Left SLF 1 Fibre Density (FDC)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  #theme_classic() +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  coord_flip()
#Left SLF 2: 
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF2_L, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mn_FDC_SLF2_L, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Left SLF 2 Fibre Density (FDC)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme(legend.position = "none") +
  theme(panel.background = element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  coord_flip()
#Left SLF 3: 
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF3_L, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mn_FDC_SLF3_L, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Left SLF 3 Fibre Density (FDC)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme(legend.position = "none") +
  theme(panel.background = element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  coord_flip()
#Right SLF 1: 
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF1_R, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mn_FDC_SLF1_R, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Right SLF 1 Fibre Density (FDC)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme(legend.position = "none") +
  theme(panel.background = element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  coord_flip()
#Right SLF 2: 
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF2_R, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mn_FDC_SLF2_R, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Right SLF 2 Fibre Density (FDC)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme(legend.position = "none") +
  theme(panel.background = element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  coord_flip()
#Right SLF 3: 
ggplot(SLF_data, aes(x = Group, y = mn_FDC_SLF3_R, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mn_FDC_SLF3_R, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Right SLF 3 Fibre Density (FDC)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme(legend.position = "none") +
  theme(panel.background = element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
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


#FD - median
# ggplot(SLF_data, aes(x = Group, y = md_FD)) + 
#     geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
#     xlab("Group_status") + 
#     ylab("FD") +
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

#for FD: just with left and right SLF 1, 2, 3
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


#Separate the groups - All tracts SLF FD (raincloud plot)
ggplot(SLF_data_FD_long, aes(x = SLF_type, y = FD_metric, fill = Group)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  #geom_point(aes(y = FD_metric, colour = Group), size = .5, alpha = 0.8, position = position_dodge(width=0.5)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fibre Density (FD)") +
  scale_x_discrete(labels = c("mn_FD_SLF1_L" = "Left SLF1","mn_FD_SLF2_L" = "Left SLF2","mn_FD_SLF3_L" = "Left SLF3","mn_FD_SLF1_R" = "Right SLF1","mn_FD_SLF2_R" = "Right SLF2","mn_FD_SLF3_R" = "Right SLF3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()

#for FC
SLF_data_FC <- dplyr::select(SLF_data, 
                             ParticipantID,
                             Group,
                             mn_FC_SLF1_L,
                             mn_FC_SLF2_L,
                             mn_FC_SLF3_L,
                             mn_FC_SLF1_R,
                             mn_FC_SLF2_R,
                             mn_FC_SLF3_R)
SLF_data_FC_long <- gather(SLF_data_FC, 
                           "SLF_type",
                           "FC_metric",
                           mn_FC_SLF1_L,
                           mn_FC_SLF2_L,
                           mn_FC_SLF3_L,
                           mn_FC_SLF1_R,
                           mn_FC_SLF2_R,
                           mn_FC_SLF3_R)
#All tracts SLF FC (raincloud plot)
ggplot(SLF_data_FC_long, aes(x = SLF_type, y = FC_metric, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = FC_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fibre Cross-section (FC)") +
  scale_x_discrete(labels = c("mn_FC_SLF1_L" = "Left SLF1","mn_FC_SLF2_L" = "Left SLF2","mn_FC_SLF3_L" = "Left SLF3","mn_FC_SLF1_R" = "Right SLF1","mn_FC_SLF2_R" = "Right SLF2","mn_FC_SLF3_R" = "Right SLF3")) + 
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

#just right SLF3
SLF_data_FDC <- dplyr::select(SLF_data, 
                              ParticipantID,
                              Group,
                              mn_FDC_SLF3_R)

SLF_data_FDC_long <- gather(SLF_data_FDC, 
                            "SLF_type",
                            "FDC_metric",
                            mn_FDC_SLF3_R)
#just right SLF 3 FDC (raincloud plot)
ggplot(SLF_data_FDC_long, aes(x = SLF_type, y = FDC_metric, fill = Group)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fibre Density Cross-section (FDC)") +
  scale_x_discrete(labels = c("mn_FDC_SLF3_R" = "Right SLF3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()

#just for left SLF 2 right SLF 3 FDC (for sex and age covariates)
#just right SLF3
SLF_data_FDC <- dplyr::select(SLF_data, 
                              ParticipantID,
                              Group,
                              mn_FDC_SLF2_L,
                              mn_FDC_SLF3_R)

SLF_data_FDC_long <- gather(SLF_data_FDC, 
                            "SLF_type",
                            "FDC_metric",
                            mn_FDC_SLF2_L,
                            mn_FDC_SLF3_R)
#just for left SLF 2 right SLF 3 FDC (for sex and age covariates)
ggplot(SLF_data_FDC_long, aes(x = SLF_type, y = FDC_metric, fill = Group)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fibre Density Cross-section (FDC)") +
  scale_x_discrete(labels = c("mn_FDC_SLF2_L" = "Left SLF 2", "mn_FDC_SLF3_R" = "Right SLF 3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()


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
#plot your correlation plot!
corrplot(cor_matrix_vars_interest, method = "ellipse", tl.col = "black", is.corr=FALSE)


#correlation plot from PLS analysis (from datamatcorrs_lst from the data structure):
#for all 11 behavioural data:
#for FD:
#barplot
neuropsych_test_names_PLS <- c('TMT-A','Colour Naming','Word Reading','HayTime1','TMT-B','Inhibition','Category Switching','HayTime2','HayTotError','Letter Fluency','Category Fluency')
PLS_FD_bootstrap_corr_values<- c(0.010482057,
                                 0.15219887,
                                 0.070730224,
                                 0.057013381,
                                 0.18211010,
                                 0.17940941,
                                 0.13699709,
                                 0.14642851,
                                 0.019853018,
                                 0.17247146,
                                 0.14759645)*-1 
ulimit_PLS_FD_barplot <- c(0.134872026741505,
                           0.252722710371018,
                           0.177258089184761,
                           0.174276560544968,
                           0.307097449898720,
                           0.290097609162331,
                           0.229116745293140,
                           0.256035745143891,
                           0.125543452799320,
                           0.264487668871880,
                           0.256653323769569)*-1
llimit_PLS_FD_barplot <- c(-0.100185252726078,
                           0.0753221325576305,
                           -0.0126400841400027,
                           -0.0466973539441824,
                           0.0691301897168160,
                           0.0889016464352608,
                           0.0615855157375336,
                           0.0451161600649357,
                           -0.0644387267529964,
                           0.100426219403744,
                           0.0627997778356075)*-1

significance_legend<-c('Does Not Reliably Contribute to Latent Variable','Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable','Does Not Reliably Contribute to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable')
df_PLS_FD_barplot <- data.frame(neuropsych_test_names_PLS,significance_legend,PLS_FD_bootstrap_corr_values,ulimit_PLS_FD_barplot,llimit_PLS_FD_barplot)
#convert to factor variables
df_PLS_FD_barplot$neuropsych_test_names_PLS <- factor(df_PLS_FD_barplot$neuropsych_test_names_PLS, levels=c('TMT-A','Colour Naming','Word Reading','HayTime1','TMT-B','Inhibition','Category Switching','HayTime2','HayTotError','Letter Fluency','Category Fluency'))
df_PLS_FD_barplot$significance_legend <- factor(df_PLS_FD_barplot$significance_legend, levels=c('Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable'))
#barplot(PLS_FD_bootstrap_corr_values, xlab = "Neuropsychological Assessment", ylab = "Bootstrap Correlation Value")
ggplot(data=df_PLS_FD_barplot,(aes(x=neuropsych_test_names_PLS, y=PLS_FD_bootstrap_corr_values,fill=significance_legend)))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(x=neuropsych_test_names_PLS,ymin=llimit_PLS_FD_barplot,ymax=ulimit_PLS_FD_barplot))+
  xlab("Neuropsychological Assessment") + 
  ylab("Correlation Values") +
  scale_fill_manual(values=c("steelblue","light grey"))+
  labs(fill=NULL)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) #adjust x-labels
  
  
#for all 11 behavioural data (from datamatcorrs_lst from the data structure) (both neuropsych & SLF data regressed by age):
#heatmap:
PLS_output_FD_data <- c(0.021629179,	-0.0066344580,	0.058540680,	0.020988390,-0.065217622,	-0.010799452,
                        -0.079838991,	-0.13450210,	-0.13028276,	-0.026898282,	-0.15641662,	-0.13970323,
                        -0.025354663,	-0.074115209,	-0.047858231,	0.0029593357,	-0.072819188,	-0.065413959,
                        -0.063448519,	-0.11186280,	-0.030578431,	-0.049415123,	0.0038374958,	-0.0089985961,
                        -0.051998116,	-0.18598628,	-0.094575949,	-0.12659493,	-0.18829104,	-0.15125458,
                        -0.066663720,	-0.17047901,	-0.17376935,	-0.081262514,	-0.15126234,	-0.17010324,
                        -0.013579339,	-0.18600455,	-0.12257311,	-0.041927282,	-0.11327523,	-0.048642337,
                        -0.049725801,	-0.17303550,	-0.071519323,	-0.021273699,	-0.14951769,	-0.11014353,
                        -0.051047891,	-0.074130818,	-0.040458661,	0.017251788,	0.027025491,	0.060011562,
                        -0.014128522,	-0.17156444,	-0.17162377,	-0.057851773,	-0.16890146,	-0.13620061,
                        -0.0086901551,	-0.17256962,	-0.14729679,	0.019446759,	-0.14279711,	-0.095593415) 
#transpose the values
#PLS_FD_corr_matrix <- t(matrix(PLS_output_FD_data,nrow=6,ncol=12))
#create correlation matrix
PLS_FD_corr_matrix <- matrix(PLS_output_FD_data,nrow=6,ncol=11)
#add in row and column names to the matrix
colnames(PLS_FD_corr_matrix) <-c('TMT-A','Colour Naming','Word Reading','Haytime1','TMT-B','Inhibition','Category Switching','HayTime2','HayTotError','Letter Fluency','Category Fluency')
rownames(PLS_FD_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
#corrplot(PLS_FD_corr_matrix, method = "color", tl.col = "black", col=colorRampPalette(c("purple","orange"))(200), cl.lim=c(-.25,.10), is.corr=FALSE)
#cols <-brewer.pal(11,"PuOr")
#pal<-colorRampPalette(cols)
#corrplot(PLS_FD_corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
corrplot(PLS_FD_corr_matrix, method = "color", tl.col = "black",col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)

#for all 5 inhibition tests data:
#for FDC:
#barplot
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FDC_bootstrap_corr_values<- c(0.16149837,
                                 0.16679218,
                                 0.13629955,
                                 0.13754775,
                                 0.022352515)*-1 
ulimit_PLS_FDC_barplot <- c(0.305668428540230,
                           0.273633107542992,
                           0.230966039001942,
                           0.255266830325127,
                           0.125505715608597)*-1
llimit_PLS_FDC_barplot <- c(0.0500048007816076,
                           0.0757585018873215,
                           0.0616370067000389,
                           0.0363335870206356,
                           -0.0591596905142069)*-1

significance_legend<-c('Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable')
df_PLS_FDC_barplot <- data.frame(neuropsych_test_names_PLS,significance_legend,PLS_FDC_bootstrap_corr_values,ulimit_PLS_FDC_barplot,llimit_PLS_FDC_barplot)
#convert to factor variables
df_PLS_FDC_barplot$neuropsych_test_names_PLS <- factor(df_PLS_FDC_barplot$neuropsych_test_names_PLS, levels=c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError'))
df_PLS_FDC_barplot$significance_legend <- factor(df_PLS_FDC_barplot$significance_legend, levels=c('Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable'))
#barplot(PLS_FDC_bootstrap_corr_values, xlab = "Neuropsychological Assessment", ylab = "Bootstrap Correlation Value")
ggplot(data=df_PLS_FDC_barplot,(aes(x=neuropsych_test_names_PLS, y=PLS_FDC_bootstrap_corr_values,fill=significance_legend)))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(x=neuropsych_test_names_PLS,ymin=llimit_PLS_FDC_barplot,ymax=ulimit_PLS_FDC_barplot))+
  xlab("Neuropsychological Assessment") + 
  ylab("Correlation Values") +
  scale_fill_manual(values=c("steelblue","light grey"))+
  labs(fill=NULL)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) #adjust x-labels


#for inhibition 5 FDC:
PLS_output_FDC_inhib5_data <- c(-0.045088284,	-0.16490521,	-0.12426596,	-0.098743699,	-0.13778134,	-0.15541276,
                               -0.056572367,	-0.15917349,	-0.16711017,	-0.089174040,	-0.13429688,	-0.15822305,
                               -0.012648485,	-0.18558364,	-0.12191153,	-0.041127983,	-0.11259479,	-0.047797780,
                               -0.038812712,	-0.15740234,	-0.068443902,	-0.011714508,	-0.15725270,	-0.11359754,
                               -0.051465657,	-0.074342854,	-0.040595248,	0.017224338,	0.027178751,	0.059770372)
#create correlation matrix
PLS_FDC_inhib5_corr_matrix <- matrix(PLS_output_FDC_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FDC_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FDC_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FDC_inhib5_corr_matrix, method = "color", tl.col = "black",col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)

# #for FDC:
# PLS_output_FDC_data <- c(-0.023494439,	-0.10839374,	-0.017945150,	-0.072824180,	-0.097034164,	-0.024255082,
#                          -0.081275865,	-0.15651976,	-0.078872435,	-0.14735182,	-0.16034853,	-0.10748159,
#                          -0.094925098,	-0.19247155,	-0.15371549,	-0.13742670,	-0.13599111,	-0.14026378,
#                          -0.034547165,	-0.13114275,	-0.096889511,	-0.070857622,	-0.097768739,	-0.063849203,
#                          -0.10721113,	-0.20132618,	-0.15642110,	-0.15865122,	-0.14065360,	-0.15279879,
#                          -0.060784571,	-0.17817909,	-0.13939866,	-0.11871763,	-0.12281568,	-0.11398516,
#                          -0.047618702,	-0.15088943,	-0.10603221,	-0.060467798,	-0.098826982,	-0.071290538,
#                          -0.061056864,	-0.17408887,	-0.10280905,	-0.12083188,	-0.11882681,	-0.045651712,
#                          -0.066071272,	-0.12010130,	-0.079085797,	-0.070842467,	-0.016597578,	-0.050753083,
#                          -0.12013579,	-0.18223056,	-0.093465224,	-0.11608600,	-0.15480573,	-0.081962451,
#                          -0.0082467366,	-0.0079095457,	-0.015302321,	-0.060224246,	0.043391790,	0.023135012,
#                          -0.035908546,	-0.037020098,	-0.017377229,	-0.0095059797,	-0.026341366,	0.041452959)
# #create correlation matrix
# PLS_FDC_corr_matrix <- matrix(PLS_output_FDC_data,nrow=6,ncol=12)
# #add in row and column names to the matrix
# colnames(PLS_FDC_corr_matrix) <-c('TMT-A','TMT-B','Colour Naming','Word Reading','Inhibition','Letter Fluency','Category Fluency','Category Switching','Haytime1','HayTime2','HayCatAError','HayCatBError')
# rownames(PLS_FDC_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# 
# #for all 4 inhibition scores (raw):
# #for FDC (raw):
# PLS_output_FDC_inhib_data <- c(-0.059311923,	-0.23900978,	-0.16617303,	-0.21059912,	-0.21506317,	-0.17036039,
#                               -0.055880506,	-0.22353986,	-0.22384585,	-0.13991968,	-0.16974939,	-0.16966677,
#                               -0.010789142,	-0.24867000,	-0.16986607,	-0.10432350,	-0.14096898,	-0.043302536,
#                               -0.056938894,	-0.24372567,	-0.12331681,	-0.083641596,	-0.17522530,	-0.12001289)
# #create correlation matrix
# PLS_FDC_inhib_corr_matrix <- matrix(PLS_output_FDC_inhib_data,nrow=6,ncol=4)
# #add in row and column names to the matrix
# colnames(PLS_FDC_inhib_corr_matrix) <-c('TMT-B','Inhibition','Category Switching','HayTime2')
# rownames(PLS_FDC_inhib_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_inhib_corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #for FDC (raw):
# PLS_output_FDC_inhib_data <- c(-0.081275888,-0.15651974,	-0.078872435,	-0.14735183,	-0.16034850,	-0.10748158,
#                                -0.10721114,	-0.20132613,	-0.15642102,	-0.15865122,	-0.14065360,	-0.15279879,
#                                -0.061056860,-0.17408888,	-0.10280903,	-0.12083188,	-0.11882683,	-0.045651726,
#                                -0.12013578,	-0.18223058,	-0.093465194,	-0.11608600,	-0.15480575,	-0.081962429)
# #create correlation matrix
# PLS_FDC_inhib_corr_matrix <- matrix(PLS_output_FDC_inhib_data,nrow=6,ncol=4)
# #add in row and column names to the matrix
# colnames(PLS_FDC_inhib_corr_matrix) <-c('TMT-B','Inhibition','Category Switching','HayTime2')
# rownames(PLS_FDC_inhib_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_inhib_corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #for all 4 inhibition zscores:
# #for FDC (zscore):
# PLS_output_FDC_inhib_zdata <- c(-0.0013060169,0.17507555,	0.10720466,	0.13353334,	0.20028476,	0.16079105,
#                                0.041074198,	 0.18999904,	0.17024736,	0.094471522,0.16638769,	0.12157550,
#                                0.035442751,	 0.26104787,	0.20363601,	0.10262322,	0.17455740,	0.085880652,
#                                0.012074957,	 0.10491857,	0.050221156,-0.019955613,0.072782405,0.076411277)
# #create correlation matrix
# PLS_FDC_inhib_zcorr_matrix <- matrix(PLS_output_FDC_inhib_zdata,nrow=6,ncol=4)
# #add in row and column names to the matrix
# colnames(PLS_FDC_inhib_zcorr_matrix) <-c('TMT-B','Inhibition','Category Switching','HayTime2')
# rownames(PLS_FDC_inhib_zcorr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_inhib_zcorr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.15,.30), is.corr=FALSE)
# 
# #for FDC (zscore):
# PLS_output_FDC_inhib_zdata <- c(0.00096392963,	0.094259106,-0.0015320014,0.060975134,	0.12675491,	0.052631572,
#                                 0.054371715,	0.16551995,	0.10425925,	  0.083246164,	0.11949693,	0.089968845,
#                                 0.083747402,	0.20385469,	0.14831589,	  0.13770595,	0.17115739,	0.083748452,
#                                 0.034849510,	0.084639393,0.031796504,  0.018525939,	0.054241844,0.022864724)
# #create correlation matrix
# PLS_FDC_inhib_zcorr_matrix <- matrix(PLS_output_FDC_inhib_zdata,nrow=6,ncol=4)
# #add in row and column names to the matrix
# colnames(PLS_FDC_inhib_zcorr_matrix) <-c('TMT-B','Inhibition','Category Switching','HayTime2')
# rownames(PLS_FDC_inhib_zcorr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_inhib_zcorr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.15,.30), is.corr=FALSE)

# #for 7 inhibition (inhibition + inhib rate) scores:
# #for FDC:
# PLS_output_FDC_inhibRate_data <- c(-0.059311923,	-0.23900978,-0.16617303,-0.21059912,-0.21506317,-0.17036039,
#                               -0.055880506,	-0.22353986,-0.22384585,-0.13991968,-0.16974939,-0.16966677,
#                               -0.010789142,	-0.24867000,-0.16986607,-0.10432350,-0.14096898,-0.043302536,
#                               -0.056938894,	-0.24372567,-0.12331681,-0.083641596,-0.17522530,-0.12001289,
#                               0.012781640,	-0.10718581,-0.12468886,-0.077844985,-0.052113518,-0.062034946,
#                               -0.049341843,	-0.15497588,-0.18023990,-0.11758253,-0.11222965,-0.12720178,
#                               -0.023481816,	-0.13349116,-0.12111595,-0.10928466,-0.12700854,-0.14178012)
# #transpose the values
# PLS_FDC_inhibRate_corr_matrix <- t(matrix(PLS_output_FDC_inhibRate_data,nrow=6,ncol=7))
# #add in row and column names to the matrix
# rownames(PLS_FDC_inhibRate_corr_matrix) <-c('TMT-B','Inhibition','Category Switching','HayTime2','TMT-B/TMT-A','Inhibition/Colour Naming','Inhibition/Word Reading')
# colnames(PLS_FDC_inhibRate_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# corrplot(PLS_FDC_inhibRate_corr_matrix, method = "color", tl.col = "black", cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #for FDC:
# PLS_output_FDC_inhib_data <- c(-0.081275888,-0.15651974,	-0.078872435,	-0.14735183,	-0.16034850,	-0.10748158,
#                                -0.10721114,	-0.20132613,	-0.15642102,	-0.15865122,	-0.14065360,	-0.15279879,
#                                -0.061056860,-0.17408888,	-0.10280903,	-0.12083188,	-0.11882683,	-0.045651726,
#                                -0.12013578,	-0.18223058,	-0.093465194,	-0.11608600,	-0.15480575,	-0.081962429,
#                                -0.040008143,-0.076053008,	-0.071811348,	-0.064911842,	-0.050623927,	-0.058439828,
#                                -0.097090632,-0.13097134,	-0.11562419,	-0.12987293,	-0.084710285,	-0.11885823,
#                                -0.056744035,-0.082488939,	-0.065527640,	-0.079254553,	-0.10559794,	-0.093407415)
# #transpose the values
# PLS_FDC_inhib_corr_matrix <- t(matrix(PLS_output_FDC_inhib_data,nrow=6,ncol=7))
# #add in row and column names to the matrix
# rownames(PLS_FDC_inhib_corr_matrix) <-c('TMT-B','Inhibition','Category Switching','HayTime2','TMT-B/TMT-A','Inhibition/Colour Naming','Inhibition/Word Reading')
# colnames(PLS_FDC_inhib_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# corrplot(PLS_FDC_inhib_corr_matrix, method = "color", tl.col = "black", cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #for 3 inhibition rate scores (raw):
# #for FDC (raw):
# PLS_output_FDC_3inhibRate_data <- c(0.012781640,	-0.10718581,	-0.12468886,	-0.077844985,	-0.052113518,	-0.062034946,
#                                    -0.049341843,-0.15497588,	-0.18023990,	-0.11758253,	-0.11222965,	-0.12720178,
#                                    -0.023481816	,-0.13349116,	-0.12111595,	-0.10928466,	-0.12700854,	-0.14178012)
# #create correlation matrix
# PLS_FDC_3inhibRate_corr_matrix <- matrix(PLS_output_FDC_3inhibRate_data,nrow=6,ncol=3)
# #add in row and column names to the matrix
# colnames(PLS_FDC_3inhibRate_corr_matrix) <-c('TMT-B/TMT-A','Inhibition/Colour Naming','Inhibition/Word Reading')
# rownames(PLS_FDC_3inhibRate_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_3inhibRate_corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #for 3 inhibition rate scores (z-score):
# #for FDC (z-score):
# PLS_output_FDC_3inhibRate_zdata <- c(-0.014584403,	0.10287708,	0.11998189,	0.072596632,0.050804809,0.060715336,
#                                    0.046708338,	    0.14920010,	0.17539965,	0.11107749,	0.11057993,	0.12645018,
#                                    0.023481816,	    0.13349117,	0.12111593,	0.10928468,	0.12700854,	0.14178011)
# #create correlation matrix
# PLS_FDC_3inhibRate_zcorr_matrix <- matrix(PLS_output_FDC_3inhibRate_zdata,nrow=6,ncol=3)
# #add in row and column names to the matrix
# colnames(PLS_FDC_3inhibRate_zcorr_matrix) <-c('TMT-B/TMT-A','Inhibition/Colour Naming','Inhibition/Word Reading')
# rownames(PLS_FDC_3inhibRate_zcorr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_3inhibRate_zcorr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.15,.30), is.corr=FALSE)
# 
# 
# 
# #for 4 processing speed scores (raw):
# #for FDC (raw):
# PLS_output_FDC_ProcSpeed_data <- c(-0.052887432,	-0.19881813,	-0.098845303,	-0.15199682,	-0.18395258,	-0.10515217,
#                                   -0.081713542,	-0.19632770,	-0.18399397,	-0.090451822,	-0.17019776,	-0.15496008,
#                                   -0.021202313,	-0.11567444,	-0.084804416,	-0.039748233,	-0.099011421,	-0.065692201,
#                                   -0.074153230,	-0.16982189,	-0.071893364,	-0.10415857,	-0.012070606,	-0.027570801)
# #create correlation matrix
# PLS_FDC_ProcSpeed_corr_matrix <- matrix(PLS_output_FDC_ProcSpeed_data,nrow=6,ncol=4)
# #add in row and column names to the matrix
# colnames(PLS_FDC_ProcSpeed_corr_matrix) <-c('TMT-A','Colour Naming','Word Reading','Haytime1')
# rownames(PLS_FDC_ProcSpeed_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_ProcSpeed_corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #for FDC (raw):
# PLS_output_FDC_ProcSpeed_data <- c(-0.023494413,	-0.10839376,	-0.017945161,	-0.072824195,	-0.097034149,	-0.024255071,
#                                   -0.094925076,	-0.19247152,	-0.15371546,	-0.13742673,	-0.13599110,	-0.14026378,
#                                   -0.034547158,	-0.13114278,	-0.096889503,	-0.070857659,	-0.097768746,	-0.063849173,
#                                   -0.066071294,	-0.12010131,	-0.079085760,	-0.070842393,	-0.016597575,	-0.050753113)
# #create correlation matrix
# PLS_FDC_ProcSpeed_corr_matrix <- matrix(PLS_output_FDC_ProcSpeed_data,nrow=6,ncol=4)
# #add in row and column names to the matrix
# colnames(PLS_FDC_ProcSpeed_corr_matrix) <-c('TMT-A','Colour Naming','Word Reading','Haytime1')
# rownames(PLS_FDC_ProcSpeed_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_ProcSpeed_corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #for 4 processing speed scores (zscores):
# #for FDC (zscores):
# PLS_output_FDC_ProcSpeed_zdata <- c(0.0023706811,0.11799329,	0.019205261,0.066913955,	0.14844789,	0.082191922,
#                                    0.074899413,	0.18496436,	0.15534791,	0.046362143,	0.17685394,	0.15339230,
#                                    0.0091149900,0.080202527,0.035712484,-0.0028785698,	0.078847423,0.050994072,
#                                    0.094311088,	0.17657560,	0.065742359,0.085854836,	0.050265390,0.050757453)
# #create correlation matrix
# PLS_FDC_ProcSpeed_zcorr_matrix <- matrix(PLS_output_FDC_ProcSpeed_zdata,nrow=6,ncol=4)
# #add in row and column names to the matrix
# colnames(PLS_FDC_ProcSpeed_zcorr_matrix) <-c('TMT-A','Colour Naming','Word Reading','Haytime1')
# rownames(PLS_FDC_ProcSpeed_zcorr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_ProcSpeed_zcorr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.15,.30), is.corr=FALSE)
# 
# #for all 13 behavioural scores (no HayError and have included Inhibition Rate scores):
# #for FDC:
# PLS_output_FDC_13data <- c(-0.052887432,	-0.19881806,	-0.098845303,	-0.15199679,	-0.18395256,    -0.10515221,
#                         -0.059311930,	-0.23900978,	-0.16617306,	-0.21059908,	-0.21506318,	-0.17036040,
#                         -0.081713542,	-0.19632770,	-0.18399400,	-0.090451777,	-0.17019773,	-0.15496007,
#                         -0.021202305,	-0.11567443,	-0.084804401,	-0.039748240,	-0.099011421,	-0.065692201,
#                         -0.055880498,	-0.22353993,	-0.22384588,	-0.13991973,	-0.16974936,	-0.16966677,
#                         -0.012210054,	-0.22325370,	-0.21712132,	-0.10130999,	-0.17594329,	-0.14812309,
#                         -0.0098307850,	-0.22235556,	-0.20072870,	-0.048246328,	-0.15880926,	-0.10259723,
#                         -0.010789134,	-0.24867000,	-0.16986604,	-0.10432348,	-0.14096898,	-0.043302517,
#                         -0.074153244,	-0.16982186,	-0.071893349,	-0.10415857,	-0.012070614,	-0.027570795,
#                         -0.056938890,	-0.24372563,	-0.12331680,	-0.083641589,	-0.17522530,	-0.12001289,
#                         0.012781648,	-0.10718579,	-0.12468886,	-0.077844992,	-0.052113514,	-0.062034946,
#                         -0.049341865,	-0.15497589,	-0.18023989,	-0.11758254,	-0.11222964,	-0.12720178,
#                         -0.023481827,	-0.13349116,	-0.12111592,	-0.10928465,	-0.12700856,	-0.14178014)
# #create correlation matrix
# PLS_FDC_13corr_matrix <- matrix(PLS_output_FDC_13data,nrow=6,ncol=13)
# #add in row and column names to the matrix
# colnames(PLS_FDC_13corr_matrix) <-c('TMT-A','TMT-B','Colour Naming','Word Reading','Inhibition','Letter Fluency','Category Fluency','Category Switching','Haytime1','HayTime2','TMT-B/TMT-A','Inhibition/Colour Naming','Inhibition/Word Reading')
# rownames(PLS_FDC_13corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_13corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #reordered variables for FDC (raw): 
# PLS_output_FDC_reordered_13data <- c(-0.052887432,	-0.19881806,	-0.098845303,	-0.15199679,	-0.18395256,	-0.10515221,
#                                     -0.081713542,	-0.19632770,	-0.18399400,	-0.090451777,	-0.17019773,	-0.15496007,
#                                     -0.021202305,	-0.11567443,	-0.084804401,	-0.039748240,	-0.099011421,	-0.065692201,
#                                     -0.074153244,	-0.16982186,	-0.071893349,	-0.10415857,	-0.012070614,	-0.027570795,
#                                     -0.059311930,	-0.23900978,	-0.16617306,	-0.21059908,	-0.21506318,	-0.17036040,
#                                     -0.055880498,	-0.22353993,	-0.22384588,	-0.13991973,	-0.16974936,	-0.16966677,
#                                     -0.010789134,	-0.24867000,	-0.16986604,	-0.10432348,	-0.14096898,	-0.043302517,
#                                     -0.056938890,	-0.24372563,	-0.12331680,	-0.083641589,	-0.17522530,	-0.12001289,
#                                     0.012781648,	-0.10718579,	-0.12468886,	-0.077844992,	-0.052113514,	-0.062034946,
#                                     -0.049341865,	-0.15497589,	-0.18023989,	-0.11758254,	-0.11222964,	-0.12720178,
#                                     -0.023481827,	-0.13349116,	-0.12111592,	-0.10928465,	-0.12700856,	-0.14178014,
#                                     -0.012210054,	-0.22325370,	-0.21712132,	-0.10130999,	-0.17594329,	-0.14812309,
#                                     -0.0098307850,	-0.22235556,	-0.20072870,	-0.048246328,	-0.15880926,	-0.10259723)
# #create correlation matrix
# PLS_FDC_reordered_13corr_matrix <- matrix(PLS_output_FDC_reordered_13data,nrow=6,ncol=13)
# #add in row and column names to the matrix
# colnames(PLS_FDC_reordered_13corr_matrix) <-c('TMT-A','Colour Naming','Word Reading','Haytime1','TMT-B','Inhibition','Category Switching','HayTime2','TMT-B/TMT-A','Inhibition/ColourNaming','Inhibition/WordReading','Letter Fluency','Category Fluency')
# rownames(PLS_FDC_reordered_13corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_reordered_13corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #for FDC (raw):
# PLS_output_FDC_13data <- c(-0.023494439,	-0.10839374,	-0.017945150,	-0.072824180,	-0.097034164,	-0.024255082,
#                            -0.081275865,	-0.15651976,	-0.078872435,	-0.14735182,	-0.16034853,	-0.10748159,
#                            -0.094925098,	-0.19247155,	-0.15371549,	-0.13742670,	-0.13599111,	-0.14026378,
#                            -0.034547165,	-0.13114275,	-0.096889511,	-0.070857622,	-0.097768739,	-0.063849203,
#                            -0.10721113,	    -0.20132618,	-0.15642110,	-0.15865122,	-0.14065360,	-0.15279879,
#                            -0.060784571,	-0.17817909,	-0.13939866,	-0.11871763,	-0.12281568,	-0.11398516,
#                            -0.047618702,	-0.15088943,	-0.10603221,	-0.060467798,	-0.098826982,	-0.071290538,
#                            -0.061056864,	-0.17408887,	-0.10280905,	-0.12083188,	-0.11882681,	-0.045651712,
#                            -0.066071272,	-0.12010130,	-0.079085797,	-0.070842467,	-0.016597578,	-0.050753083,
#                            -0.12013579,	    -0.18223056,	-0.093465224,	-0.11608600,	-0.15480573,	-0.081962451,
#                            -0.040008113,	-0.076052994,	-0.071811356,	-0.064911835,	-0.050623909,	-0.058439810,
#                            -0.097090632,	-0.13097134,	-0.11562417,	-0.12987293,	-0.084710300,	-0.11885826,
#                            -0.056744039,	-0.082488939,	-0.065527633,	-0.079254553,	-0.10559794,	-0.093407415)
# #create correlation matrix
# PLS_FDC_13corr_matrix <- matrix(PLS_output_FDC_13data,nrow=6,ncol=13)
# #add in row and column names to the matrix
# colnames(PLS_FDC_13corr_matrix) <-c('TMT-A','TMT-B','Colour Naming','Word Reading','Inhibition','Letter Fluency','Category Fluency','Category Switching','Haytime1','HayTime2','TMT-B/TMT-A','Inhibition/Colour Naming','Inhibition/Word Reading')
# rownames(PLS_FDC_13corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_13corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #reordered variables for FDC (raw): 
# PLS_output_FDC_reordered_13data <- c(-0.023494439,	-0.10839374,	-0.017945150,	-0.072824180,	-0.097034164,	-0.024255082,
#                                      -0.094925098,	-0.19247155,	-0.15371549,	-0.13742670,	-0.13599111,	-0.14026378,
#                                      -0.034547165,	-0.13114275,	-0.096889511,	-0.070857622,	-0.097768739,	-0.063849203,
#                                      -0.066071272,	-0.12010130,	-0.079085797,	-0.070842467,	-0.016597578,	-0.050753083,
#                                      -0.081275865,	-0.15651976,	-0.078872435,	-0.14735182,	-0.16034853,	-0.10748159,
#                                      -0.10721113,	-0.20132618,	-0.15642110,	-0.15865122,	-0.14065360,	-0.15279879,
#                                      -0.061056864,	-0.17408887,	-0.10280905,	-0.12083188,	-0.11882681,	-0.045651712,
#                                      -0.12013579,	-0.18223056,	-0.093465224,	-0.11608600,	-0.15480573,	-0.081962451,
#                                      -0.040008113,	-0.076052994,	-0.071811356,	-0.064911835,	-0.050623909,	-0.058439810,
#                                      -0.097090632,	-0.13097134,	-0.11562417,	-0.12987293,	-0.084710300,	-0.11885826,
#                                      -0.056744039,	-0.082488939,	-0.065527633,	-0.079254553,	-0.10559794,	-0.093407415,
#                                      -0.060784571,	-0.17817909,	-0.13939866,	-0.11871763,	-0.12281568,	-0.11398516,
#                                      -0.047618702,	-0.15088943,	-0.10603221,	-0.060467798,	-0.098826982,	-0.071290538)
# #create correlation matrix
# PLS_FDC_reordered_13corr_matrix <- matrix(PLS_output_FDC_reordered_13data,nrow=6,ncol=13)
# #add in row and column names to the matrix
# colnames(PLS_FDC_reordered_13corr_matrix) <-c('TMT-A','Colour Naming','Word Reading','Haytime1','TMT-B','Inhibition','Category Switching','HayTime2','TMT-B/TMT-A','Inhibition/ColourNaming','Inhibition/WordReading','Letter Fluency','Category Fluency')
# rownames(PLS_FDC_reordered_13corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_reordered_13corr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.25,.10), is.corr=FALSE)
# 
# #reordered variables for FDC (z-scores): 
# PLS_output_FDC_reordered_13zdata <- c(0.0023706625,	0.11799330,	0.019205257,0.066913947,	0.14844790,	0.082191937,
#                                      0.074899398,	0.18496436,	0.15534791,	0.046362136,	0.17685390,	0.15339229,
#                                      0.0091150058,	0.080202542,0.035712466,-0.0028785770,	0.078847438,0.050994072,
#                                      0.094311066,	0.17657560,	0.065742359,0.085854851,	0.050265402,0.050757475,
#                                      -0.0013060420,	0.17507555,	0.10720468,	0.13353336,	    0.20028478,	0.16079105,
#                                      0.041074220,	0.18999907,	0.17024738,	0.094471514,	0.16638766,	0.12157548,
#                                      0.035442743,	0.26104799,	0.20363599,	0.10262320,	    0.17455740,	0.085880637,
#                                      0.012074946,	0.10491858,	0.050221138,-0.019955594,	0.072782449,0.076411277,
#                                      -0.014584406,	0.10287707,	0.11998191,	0.072596632,	0.050804805,0.060715329,
#                                      0.046708338,	0.14920017,	0.17539963,	0.11107749,	    0.11057995,	0.12645015,
#                                      0.023481827,	0.13349116,	0.12111595,	0.10928469,	    0.12700854,	0.14178012,
#                                      0.0091575170,	0.22464813,	0.20742762,	0.088974096,	0.16642311,	0.13465299,
#                                      0.024852181,	0.21469341,	0.15587012,	0.0050964309,	0.16240194,	0.091314487)
# #create correlation matrix
# PLS_FDC_reordered_13zcorr_matrix <- matrix(PLS_output_FDC_reordered_13zdata,nrow=6,ncol=13)
# #add in row and column names to the matrix
# colnames(PLS_FDC_reordered_13zcorr_matrix) <-c('TMT-A','Colour Naming','Word Reading','Haytime1','TMT-B','Inhibition','Category Switching','HayTime2','TMT-B/TMT-A','Inhibition/ColourNaming','Inhibition/WordReading','Letter Fluency','Category Fluency')
# rownames(PLS_FDC_reordered_13zcorr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_reordered_13zcorr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.15,.30), is.corr=FALSE)
# 
# 
# #reordered variables for FDC (z-scores): 
# PLS_output_FDC_reordered_13zdata <- c(-0.070726469,	0.027615108,-0.079384156,-0.030534575,	0.043896325,-0.045688041,
#                                       0.059925966,	0.17556132,	0.11494567,	  0.082713202,	0.12531365,	0.11277279,
#                                       0.00071619853,0.098389834,0.048477259,  0.025740331,	0.077027030,0.034304988,
#                                       0.065731406,	0.12045781,	0.063242704,  0.045515552,	0.034297023,0.051863890,
#                                       0.00096395350,0.094259106,-0.0015320028,0.060975142,	0.12675494,	0.052631572,
#                                       0.054371707,	0.16551991,	0.10425927,	  0.083246179,	0.11949694,	0.089968838,
#                                       0.083747402,	0.20385467,	0.14831594,	  0.13770594,	0.17115733,	0.083748452,
#                                       0.034849506,	0.084639423,0.031796504,  0.018525930,	0.054241829,0.022864733,
#                                       0.040074956,	0.075706184,0.071926512,  0.063679852,	0.052315231,0.059761856,
#                                       0.096990839,	0.13007848,	0.11667474,	  0.12811339,	0.086751521,0.12117917,
#                                       0.056744024,	0.082488939,0.065527640,  0.079254553,	0.10559792,	0.093407407,
#                                       0.022846149,	0.14405321,	0.10015102,	  0.074285634,	0.083914496,0.074618369,
#                                       0.045237929,	0.14372040,	0.058784656,  0.031497598,	0.098209262,0.039346427)
# #create correlation matrix
# PLS_FDC_reordered_13zcorr_matrix <- matrix(PLS_output_FDC_reordered_13zdata,nrow=6,ncol=13)
# #add in row and column names to the matrix
# colnames(PLS_FDC_reordered_13zcorr_matrix) <-c('TMT-A','Colour Naming','Word Reading','Haytime1','TMT-B','Inhibition','Category Switching','HayTime2','TMT-B/TMT-A','Inhibition/ColourNaming','Inhibition/WordReading','Letter Fluency','Category Fluency')
# rownames(PLS_FDC_reordered_13zcorr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
# #plot correlation heatmap
# cols <-brewer.pal(11,"PuOr")
# pal<-colorRampPalette(cols)
# corrplot(PLS_FDC_reordered_13zcorr_matrix, method = "color", tl.col = "black", col=pal(20), cl.lim=c(-.15,.30), is.corr=FALSE)



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
cor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
#Left SLF vs processing speed
cor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
cor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
#Right SLF vs processing speed
cor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore)
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
cor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore)
#Left SLF vs inhibition
cor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore)
cor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore)
#Right SLF vs inhibition
cor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore)
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
cor.test(generation_zscores_data_noNAS$mn_FDC_SLF, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FC_SLF, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FDC_SLF, generation_zscores_data_noNAS$generation_average_zscore)
#Left SLF vs generation
cor.test(generation_zscores_data_noNAS$mn_FDC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FDC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore)
#Right SLF vs generation
cor.test(generation_zscores_data_noNAS$mn_FDC_SLF_R, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FC_SLF_R, generation_zscores_data_noNAS$generation_average_zscore)
cor.test(generation_zscores_data_noNAS$mn_FDC_SLF_R, generation_zscores_data_noNAS$generation_average_zscore)

#examine partial correlation (age and sex as covariates) for this:
#Processing speed
#Whole SLF vs processing speed -
pcor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Left SLF vs processing speed -
pcor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_L, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Right SLF vs processing speed -
pcor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FC_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(proc_speed_zscores_data_noNAS$mn_FDC_SLF_R, proc_speed_zscores_data_noNAS$proc_speed_average_zscore, proc_speed_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Inhibition
#Whole SLF vs inhibition -
pcor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Left SLF vs inhibition -
pcor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_L, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Right SLF vs inhibition -
pcor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FC_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(inhibition_zscores_data_noNAS$mn_FDC_SLF_R, inhibition_zscores_data_noNAS$inhibition_average_zscore, inhibition_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Generation
#Whole SLF vs generation -
pcor.test(generation_zscores_data_noNAS$mn_FDC_SLF, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FC_SLF, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FDC_SLF, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Left SLF vs generation -
pcor.test(generation_zscores_data_noNAS$mn_FDC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
pcor.test(generation_zscores_data_noNAS$mn_FDC_SLF_L, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
#Right SLF vs generation -
pcor.test(generation_zscores_data_noNAS$mn_FDC_SLF_R, generation_zscores_data_noNAS$generation_average_zscore, generation_zscores_data_noNAS[,c("Age", "Sex_binary")])
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
    map(~corr.test(x = .x %>% dplyr::select(mn_FDC_SLF, TrailsA.Raw),
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
# SLF_FDC_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF, SLF_data$Three_Groups)
# SLF_FC_3groups_descrip <- describeBy(SLF_data$mn_FC_SLF, SLF_data$Three_Groups)
# SLF_FDC_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF, SLF_data$Three_Groups)
# 
# #look at descriptive stats of the Left SLF FBA metrics between groups
# SLF_FDC_L_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF_L, SLF_data$Three_Groups)
# SLF_FC_L_3groups_descrip <- describeBy(SLF_data$mn_FC_SLF_L, SLF_data$Three_Groups)
# SLF_FDC_L_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF_L, SLF_data$Three_Groups)
# 
# 
# #look at descriptive stats of the Right SLF FBA metrics between groups
# SLF_FDC_R_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF_R, SLF_data$Three_Groups)
# SLF_FC_R_3groups_descrip <- describeBy(SLF_data$mn_FC_SLF_R, SLF_data$Three_Groups)
# SLF_FDC_R_3groups_descrip <- describeBy(SLF_data$mn_FDC_SLF_R, SLF_data$Three_Groups)
# 
# 
# #run ANOVA to see if there are significant differences between groups
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
# #calculate the effect size (eta-squared)
# etaSquared(SLF_FDC_3groups_mod)
# etaSquared(SLF_FDC_3groups_mod_L)
# etaSquared(SLF_FDC_3groups_mod_R)
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
# post_hoc_SLF_FDC_3groups_mod <- glht(SLF_FDC_3groups_mod, linfct = mcp(Three_Groups = "Tukey"))
# summary(post_hoc_SLF_FDC_3groups_mod)
# confint(post_hoc_SLF_FDC_3groups_mod)
# 
# post_hoc_SLF_FDC_3groups_mod_L <- glht(SLF_FDC_3groups_mod_L, linfct = mcp(Three_Groups = "Tukey"))
# summary(post_hoc_SLF_FDC_3groups_mod_L)
# confint(post_hoc_SLF_FDC_3groups_mod_L)
# 
# post_hoc_SLF_FDC_3groups_mod_R <- glht(SLF_FDC_3groups_mod_R, linfct = mcp(Three_Groups = "Tukey"))
# summary(post_hoc_SLF_FDC_3groups_mod_R)
# confint(post_hoc_SLF_FDC_3groups_mod_R)
# 
# #conduct power analysis for the whole FDC
# SLF_FDC_3groups_means <- c(SLF_FDC_3groups_descrip$`1`$mean, SLF_FDC_3groups_descrip$`2`$mean, SLF_FDC_3groups_descrip$`3`$mean, SLF_FDC_3groups_descrip$`4`$mean, SLF_FDC_3groups_descrip$`5`$mean)
# power_SLF_FDC_3groups_n <- power.anova.test(groups = length(SLF_FDC_3groups_means), between.var = anova(SLF_FDC_3groups_mod)$`Sum Sq`[1], within.var = anova(SLF_FDC_3groups_mod)$`Sum Sq`[2], power = .8, sig.level = 0.05)
# power_SLF_FDC_3groups_power<- power.anova.test(groups = length(SLF_FDC_3groups_means), between.var = anova(SLF_FDC_3groups_mod)$`Sum Sq`[1], within.var = anova(SLF_FDC_3groups_mod)$`Sum Sq`[2], n = 41, sig.level = 0.05)
# #conduct power analysis for the left FDC
# SLF_FDC_L_3groups_means <- c(SLF_FDC_L_3groups_descrip$`1`$mean, SLF_FDC_L_3groups_descrip$`2`$mean, SLF_FDC_L_3groups_descrip$`3`$mean, SLF_FDC_L_3groups_descrip$`4`$mean, SLF_FDC_L_3groups_descrip$`5`$mean)
# power_SLF_FDC_L_3groups_n <- power.anova.test(groups = length(SLF_FDC_L_3groups_means), between.var = anova(SLF_FDC_3groups_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FDC_3groups_mod_L)$`Sum Sq`[2], power = .8, sig.level = 0.05)
# power_SLF_FDC_L_3groups_power<- power.anova.test(groups = length(SLF_FDC_L_3groups_means), between.var = anova(SLF_FDC_3groups_mod_L)$`Sum Sq`[1], within.var = anova(SLF_FDC_3groups_mod_L)$`Sum Sq`[2], n = 41, sig.level = 0.05)
# #conduct power analysis for the right FDC
# SLF_FDC_R_3groups_means <- c(SLF_FDC_R_3groups_descrip$`1`$mean, SLF_FDC_R_3groups_descrip$`2`$mean, SLF_FDC_R_3groups_descrip$`3`$mean, SLF_FDC_R_3groups_descrip$`4`$mean, SLF_FDC_R_3groups_descrip$`5`$mean)
# power_SLF_FDC_R_3groups_n <- power.anova.test(groups = length(SLF_FDC_R_3groups_means), between.var = anova(SLF_FDC_3groups_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FDC_3groups_mod_R)$`Sum Sq`[2], power = .8, sig.level = 0.05)
# power_SLF_FDC_R_3groups_power<- power.anova.test(groups = length(SLF_FDC_R_3groups_means), between.var = anova(SLF_FDC_3groups_mod_R)$`Sum Sq`[1], within.var = anova(SLF_FDC_3groups_mod_R)$`Sum Sq`[2], n = 41, sig.level = 0.05)
# 
# 
# 
# #plot data
# #for FDC - mean
# #whole SLF FDC (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FDC)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #whole SLF FDC (raincloud plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF, fill = Three_Groups)) + 
#     geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_point(aes(y = mn_FDC_SLF, color = Three_Groups), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FDC)") +
#     scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     coord_flip()
# 
# #left SLF FDC (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF_L)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FDC)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #left SLF FDC (raincloud plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF_L, fill = Three_Groups)) + 
#     geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_point(aes(y = mn_FDC_SLF_L, color = Three_Groups), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FDC)") +
#     scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     coord_flip()
# 
# #right SLF FDC (violin plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF_R)) + 
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FDC)") +
#     scale_x_discrete(labels = c("1" = "Non-impaired", "3" = "MCI", "5" = "AD")) + 
#     theme_classic() +
#     theme(legend.position = "none") +
#     geom_violin(trim = FALSE, alpha = .5, aes(fill = Three_Groups, colour = Three_Groups), size = 1)
# 
# #right SLF FDC (raincloud plot)
# ggplot(SLF_data, aes(x = Three_Groups, y = mn_FDC_SLF_R, fill = Three_Groups)) + 
#     geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#     geom_point(aes(y = mn_FDC_SLF_R, color = Three_Groups), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Three_Groups)) + 
#     stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Three_Groups)) + 
#     xlab("Group") + 
#     ylab("Fibre Density (FDC)") +
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








#removing AD group (group 5) from analysis (covariate age & sex)
SLF_data_noAD<-SLF_data %>%
  group_by(Group) %>%
  filter(!any(Group == "5"))

#descriptive data for total sample (n = 202):
#look at descriptive stats of the Left SLF1 FBA metrics between groups
SLF1_FD_L_descrip <- describeBy(SLF_data_noAD$mn_FD_SLF1_L, SLF_data_noAD$Group)
SLF1_FC_L_descrip <- describeBy(SLF_data_noAD$mn_FC_SLF1_L, SLF_data_noAD$Group)
SLF1_FDC_L_descrip <- describeBy(SLF_data_noAD$mn_FDC_SLF1_L, SLF_data_noAD$Group)
mean(SLF_data_noAD$mn_FD_SLF1_L)
sd(SLF_data_noAD$mn_FD_SLF1_L)
mean(SLF_data_noAD$mn_FC_SLF1_L)
sd(SLF_data_noAD$mn_FC_SLF1_L)
mean(SLF_data_noAD$mn_FDC_SLF1_L)
sd(SLF_data_noAD$mn_FDC_SLF1_L)
#look at descriptive stats of the Left SLF2 FBA metrics between groups
SLF2_FD_L_descrip <- describeBy(SLF_data_noAD$mn_FD_SLF2_L, SLF_data_noAD$Group)
SLF2_FC_L_descrip <- describeBy(SLF_data_noAD$mn_FC_SLF2_L, SLF_data_noAD$Group)
SLF2_FDC_L_descrip <- describeBy(SLF_data_noAD$mn_FDC_SLF2_L, SLF_data_noAD$Group)
mean(SLF_data_noAD$mn_FD_SLF2_L)
sd(SLF_data_noAD$mn_FD_SLF2_L)
mean(SLF_data_noAD$mn_FC_SLF2_L)
sd(SLF_data_noAD$mn_FC_SLF2_L)
mean(SLF_data_noAD$mn_FDC_SLF2_L)
sd(SLF_data_noAD$mn_FDC_SLF2_L)
#look at descriptive stats of the Left SLF2 FBA metrics between groups
SLF3_FD_L_descrip <- describeBy(SLF_data_noAD$mn_FD_SLF3_L, SLF_data_noAD$Group)
SLF3_FC_L_descrip <- describeBy(SLF_data_noAD$mn_FC_SLF3_L, SLF_data_noAD$Group)
SLF3_FDC_L_descrip <- describeBy(SLF_data_noAD$mn_FDC_SLF3_L, SLF_data_noAD$Group)
mean(SLF_data_noAD$mn_FD_SLF3_L)
sd(SLF_data_noAD$mn_FD_SLF3_L)
mean(SLF_data_noAD$mn_FC_SLF3_L)
sd(SLF_data_noAD$mn_FC_SLF3_L)
mean(SLF_data_noAD$mn_FDC_SLF3_L)
sd(SLF_data_noAD$mn_FDC_SLF3_L)
#look at descriptive stats of the Right SLF1 FBA metrics between groups
SLF1_FD_R_descrip <- describeBy(SLF_data_noAD$mn_FD_SLF1_R, SLF_data_noAD$Group)
SLF1_FC_R_descrip <- describeBy(SLF_data_noAD$mn_FC_SLF1_R, SLF_data_noAD$Group)
SLF1_FDC_R_descrip <- describeBy(SLF_data_noAD$mn_FDC_SLF1_R, SLF_data_noAD$Group)
mean(SLF_data_noAD$mn_FD_SLF1_R)
sd(SLF_data_noAD$mn_FD_SLF1_R)
mean(SLF_data_noAD$mn_FC_SLF1_R)
sd(SLF_data_noAD$mn_FC_SLF1_R)
mean(SLF_data_noAD$mn_FDC_SLF1_R)
sd(SLF_data_noAD$mn_FDC_SLF1_R)
#look at descriptive stats of the Right SLF2 FBA metrics between groups
SLF2_FD_R_descrip <- describeBy(SLF_data_noAD$mn_FD_SLF2_R, SLF_data_noAD$Group)
SLF2_FC_R_descrip <- describeBy(SLF_data_noAD$mn_FC_SLF2_R, SLF_data_noAD$Group)
SLF2_FDC_R_descrip <- describeBy(SLF_data_noAD$mn_FDC_SLF2_R, SLF_data_noAD$Group)
mean(SLF_data_noAD$mn_FD_SLF2_R)
sd(SLF_data_noAD$mn_FD_SLF2_R)
mean(SLF_data_noAD$mn_FC_SLF2_R)
sd(SLF_data_noAD$mn_FC_SLF2_R)
mean(SLF_data_noAD$mn_FDC_SLF2_R)
sd(SLF_data_noAD$mn_FDC_SLF2_R)
#look at descriptive stats of the Right SLF2 FBA metrics between groups
SLF3_FD_R_descrip <- describeBy(SLF_data_noAD$mn_FD_SLF3_R, SLF_data_noAD$Group)
SLF3_FC_R_descrip <- describeBy(SLF_data_noAD$mn_FC_SLF3_R, SLF_data_noAD$Group)
SLF3_FDC_R_descrip <- describeBy(SLF_data_noAD$mn_FDC_SLF3_R, SLF_data_noAD$Group)
mean(SLF_data_noAD$mn_FD_SLF3_R)
sd(SLF_data_noAD$mn_FD_SLF3_R)
mean(SLF_data_noAD$mn_FC_SLF3_R)
sd(SLF_data_noAD$mn_FC_SLF3_R)
mean(SLF_data_noAD$mn_FDC_SLF3_R)
sd(SLF_data_noAD$mn_FDC_SLF3_R)



#ANOVA groups differences between 4 groups (C, SCD, aMCI, and mMCI)
#FD
SLF1_FD_2covar_mod_L_noAD <- lm(mn_FD_SLF1_L ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_FD_2covar_mod_L_noAD)
SLF2_FD_2covar_mod_L_noAD <- lm(mn_FD_SLF2_L ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_FD_2covar_mod_L_noAD)
SLF3_FD_2covar_mod_L_noAD <- lm(mn_FD_SLF3_L ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_FD_2covar_mod_L_noAD)
SLF1_FD_2covar_mod_R_noAD <- lm(mn_FD_SLF1_R ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_FD_2covar_mod_R_noAD)
SLF2_FD_2covar_mod_R_noAD <- lm(mn_FD_SLF2_R ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_FD_2covar_mod_R_noAD)
SLF3_FD_2covar_mod_R_noAD <- lm(mn_FD_SLF3_R ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_FD_2covar_mod_R_noAD)
#FC
SLF1_FC_2covar_mod_L_noAD <- lm(mn_FC_SLF1_L ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_FC_2covar_mod_L_noAD)
SLF2_FC_2covar_mod_L_noAD <- lm(mn_FC_SLF2_L ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_FC_2covar_mod_L_noAD)
SLF3_FC_2covar_mod_L_noAD <- lm(mn_FC_SLF3_L ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_FC_2covar_mod_L_noAD)
SLF1_FC_2covar_mod_R_noAD <- lm(mn_FC_SLF1_R ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_FC_2covar_mod_R_noAD)
SLF2_FC_2covar_mod_R_noAD <- lm(mn_FC_SLF2_R ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_FC_2covar_mod_R_noAD)
SLF3_FC_2covar_mod_R_noAD <- lm(mn_FC_SLF3_R ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_FC_2covar_mod_R_noAD)
#FDC
SLF1_FDC_2covar_mod_L_noAD <- lm(mn_FDC_SLF1_L ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_FDC_2covar_mod_L_noAD)
SLF2_FDC_2covar_mod_L_noAD <- lm(mn_FDC_SLF2_L ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_FDC_2covar_mod_L_noAD)
SLF3_FDC_2covar_mod_L_noAD <- lm(mn_FDC_SLF3_L ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_FDC_2covar_mod_L_noAD)
SLF1_FDC_2covar_mod_R_noAD <- lm(mn_FDC_SLF1_R ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_FDC_2covar_mod_R_noAD)
SLF2_FDC_2covar_mod_R_noAD <- lm(mn_FDC_SLF2_R ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_FDC_2covar_mod_R_noAD)
SLF3_FDC_2covar_mod_R_noAD <- lm(mn_FDC_SLF3_R ~ Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_FDC_2covar_mod_R_noAD)

#calculate the effect size (eta-squared) - for ANCOVA (sex & age)
#FD
etaSquared(SLF1_FD_2covar_mod_L_noAD)
etaSquared(SLF2_FD_2covar_mod_L_noAD)
etaSquared(SLF3_FD_2covar_mod_L_noAD)
etaSquared(SLF1_FD_2covar_mod_R_noAD)
etaSquared(SLF2_FD_2covar_mod_R_noAD)
etaSquared(SLF3_FD_2covar_mod_R_noAD)
#FC
etaSquared(SLF1_FC_2covar_mod_L_noAD)
etaSquared(SLF2_FC_2covar_mod_L_noAD)
etaSquared(SLF3_FC_2covar_mod_L_noAD)
etaSquared(SLF1_FC_2covar_mod_R_noAD)
etaSquared(SLF2_FC_2covar_mod_R_noAD)
etaSquared(SLF3_FC_2covar_mod_R_noAD)
#FDC
etaSquared(SLF1_FDC_2covar_mod_L_noAD)
etaSquared(SLF2_FDC_2covar_mod_L_noAD)
etaSquared(SLF3_FDC_2covar_mod_L_noAD)
etaSquared(SLF1_FDC_2covar_mod_R_noAD)
etaSquared(SLF2_FDC_2covar_mod_R_noAD)
etaSquared(SLF3_FDC_2covar_mod_R_noAD)


#linear trend
#add in trend Group variable
Trend_Group <- as.numeric(SLF_data_noAD$Group)
#FD
SLF1_L_FD_LinTrend_2covar_mod_noAD <- lm(mn_FD_SLF1_L ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_L_FD_LinTrend_2covar_mod_noAD)
summary(SLF1_L_FD_LinTrend_2covar_mod_noAD) 
SLF2_L_FD_LinTrend_2covar_mod_noAD <- lm(mn_FD_SLF2_L ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_L_FD_LinTrend_2covar_mod_noAD)
summary(SLF2_L_FD_LinTrend_2covar_mod_noAD) 
SLF3_L_FD_LinTrend_2covar_mod_noAD <- lm(mn_FD_SLF3_L ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_L_FD_LinTrend_2covar_mod_noAD)
summary(SLF3_L_FD_LinTrend_2covar_mod_noAD)
SLF1_R_FD_LinTrend_2covar_mod_noAD <- lm(mn_FD_SLF1_R ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_R_FD_LinTrend_2covar_mod_noAD)
summary(SLF1_R_FD_LinTrend_2covar_mod_noAD) 
SLF2_R_FD_LinTrend_2covar_mod_noAD <- lm(mn_FD_SLF2_R ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_R_FD_LinTrend_2covar_mod_noAD)
summary(SLF2_R_FD_LinTrend_2covar_mod_noAD) 
SLF3_R_FD_LinTrend_2covar_mod_noAD <- lm(mn_FD_SLF3_R ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_R_FD_LinTrend_2covar_mod_noAD)
summary(SLF3_R_FD_LinTrend_2covar_mod_noAD) 
#FC
SLF1_L_FC_LinTrend_2covar_mod_noAD <- lm(mn_FC_SLF1_L ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_L_FC_LinTrend_2covar_mod_noAD)
summary(SLF1_L_FC_LinTrend_2covar_mod_noAD) 
SLF2_L_FC_LinTrend_2covar_mod_noAD <- lm(mn_FC_SLF2_L ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_L_FC_LinTrend_2covar_mod_noAD)
summary(SLF2_L_FC_LinTrend_2covar_mod_noAD) 
SLF3_L_FC_LinTrend_2covar_mod_noAD <- lm(mn_FC_SLF3_L ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_L_FC_LinTrend_2covar_mod_noAD)
summary(SLF3_L_FC_LinTrend_2covar_mod_noAD)
SLF1_R_FC_LinTrend_2covar_mod_noAD <- lm(mn_FC_SLF1_R ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_R_FC_LinTrend_2covar_mod_noAD)
summary(SLF1_R_FC_LinTrend_2covar_mod_noAD) 
SLF2_R_FC_LinTrend_2covar_mod_noAD <- lm(mn_FC_SLF2_R ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_R_FC_LinTrend_2covar_mod_noAD)
summary(SLF2_R_FC_LinTrend_2covar_mod_noAD) 
SLF3_R_FC_LinTrend_2covar_mod_noAD <- lm(mn_FC_SLF3_R ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_R_FC_LinTrend_2covar_mod_noAD)
summary(SLF3_R_FC_LinTrend_2covar_mod_noAD) 
#FDC
SLF1_L_FDC_LinTrend_2covar_mod_noAD <- lm(mn_FDC_SLF1_L ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_L_FDC_LinTrend_2covar_mod_noAD)
summary(SLF1_L_FDC_LinTrend_2covar_mod_noAD) 
SLF2_L_FDC_LinTrend_2covar_mod_noAD <- lm(mn_FDC_SLF2_L ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_L_FDC_LinTrend_2covar_mod_noAD)
summary(SLF2_L_FDC_LinTrend_2covar_mod_noAD) 
SLF3_L_FDC_LinTrend_2covar_mod_noAD <- lm(mn_FDC_SLF3_L ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_L_FDC_LinTrend_2covar_mod_noAD)
summary(SLF3_L_FDC_LinTrend_2covar_mod_noAD)
SLF1_R_FDC_LinTrend_2covar_mod_noAD <- lm(mn_FDC_SLF1_R ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF1_R_FDC_LinTrend_2covar_mod_noAD)
summary(SLF1_R_FDC_LinTrend_2covar_mod_noAD) 
SLF2_R_FDC_LinTrend_2covar_mod_noAD <- lm(mn_FDC_SLF2_R ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF2_R_FDC_LinTrend_2covar_mod_noAD)
summary(SLF2_R_FDC_LinTrend_2covar_mod_noAD) 
SLF3_R_FDC_LinTrend_2covar_mod_noAD <- lm(mn_FDC_SLF3_R ~ Trend_Group + Group + Age + Sex, data = SLF_data_noAD)
anova(SLF3_R_FDC_LinTrend_2covar_mod_noAD)
summary(SLF3_R_FDC_LinTrend_2covar_mod_noAD) 

#post hoc for anova, covariates (age & sex)
#FD
#Left SLF 3
post_hoc_SLF3_L_FD_ancova_mod_noAD <- glht(SLF3_FD_2covar_mod_L_noAD, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_L_FD_ancova_mod_noAD)
confint(post_hoc_SLF3_L_FD_ancova_mod_noAD)
#effect size for sig. post hoc tests (with age as covariate)
#Left SLF3 FD
t_value_effect_size <- summary(post_hoc_SLF3_L_FD_ancova_mod_noAD) 
  #for C vs. SCD 
  DPRC_neuropsych_data_CvSCD_SLF3_L <- subset(SLF_data_noAD, SLF_data_noAD$Group == 1 | SLF_data_noAD$Group == 2)
  DPRC_neuropsych_data_CvSCD_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_L$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_SLF3_L, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age, data = DPRC_neuropsych_data_CvSCD_SLF3_L)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=as.numeric(group_number['1','n']),n.2=as.numeric(group_number['2','n']),R=r_value,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_SLF3_L <- subset(SLF_data_noAD, SLF_data_noAD$Group == 1 | SLF_data_noAD$Group == 4)
  DPRC_neuropsych_data_CvmMCI_SLF3_L$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_L$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF3_L, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age, data = DPRC_neuropsych_data_CvmMCI_SLF3_L)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=as.numeric(group_number['1','n']),n.2=as.numeric(group_number['2','n']),R=r_value,q=1) #calculate Cohen's D with the covariate of age

#Right SLF 2 and 3
post_hoc_SLF2_R_FD_ancova_mod_noAD <- glht(SLF2_FD_2covar_mod_R_noAD, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF2_R_FD_ancova_mod_noAD)
confint(post_hoc_SLF2_R_FD_ancova_mod_noAD)
#effect size for sig. post hoc tests (with age as covariate)
#Right SLF2 FD
t_value_effect_size <- summary(post_hoc_SLF2_R_FD_ancova_mod_noAD) 
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_SLF2_R <- subset(SLF_data_noAD, SLF_data_noAD$Group == 1 | SLF_data_noAD$Group == 2)
  DPRC_neuropsych_data_CvmMCI_SLF2_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF2_R$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF2_R, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF2_R ~ Age, data = DPRC_neuropsych_data_CvmMCI_SLF2_R)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=as.numeric(group_number['1','n']),n.2=as.numeric(group_number['2','n']),R=r_value,q=1) #calculate Cohen's D with the covariate of age

post_hoc_SLF3_R_FD_ancova_mod_noAD <- glht(SLF3_FD_2covar_mod_R_noAD, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_R_FD_ancova_mod_noAD)
confint(post_hoc_SLF3_R_FD_ancova_mod_noAD)
#effect size for sig. post hoc tests (with age as covariate)
#Left SLF3 FD
t_value_effect_size <- summary(post_hoc_SLF3_R_FD_ancova_mod_noAD) 
  #for C vs. SCD 
  DPRC_neuropsych_data_CvSCD_SLF3_R <- subset(SLF_data_noAD, SLF_data_noAD$Group == 1 | SLF_data_noAD$Group == 2)
  DPRC_neuropsych_data_CvSCD_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvSCD_SLF3_R$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_SLF3_R, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_R ~ Age, data = DPRC_neuropsych_data_CvSCD_SLF3_R)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=as.numeric(group_number['1','n']),n.2=as.numeric(group_number['2','n']),R=r_value,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_SLF3_R <- subset(SLF_data_noAD, SLF_data_noAD$Group == 1 | SLF_data_noAD$Group == 4)
  DPRC_neuropsych_data_CvmMCI_SLF3_R$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_SLF3_R$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_SLF3_R, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_L ~ Age, data = DPRC_neuropsych_data_CvmMCI_SLF3_R)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=as.numeric(group_number['1','n']),n.2=as.numeric(group_number['2','n']),R=r_value,q=1) #calculate Cohen's D with the covariate of age


#Plotting graphs via ggplot2
#for FD: just with left and right SLF 1, 2, 3
SLF_data_FD <- dplyr::select(SLF_data_noAD, 
                              ParticipantID,
                              Group,
                              mn_FD_SLF1_L,
                              mn_FD_SLF2_L,
                              mn_FD_SLF3_L,
                              mn_FD_SLF1_R,
                              mn_FD_SLF2_R,
                              mn_FD_SLF3_R)
SLF_data_FD_noAD_long <- gather(SLF_data_noAD, 
                            "SLF_type",
                            "FD_metric",
                            mn_FD_SLF1_L,
                            mn_FD_SLF2_L,
                            mn_FD_SLF3_L,
                            mn_FD_SLF1_R,
                            mn_FD_SLF2_R,
                            mn_FD_SLF3_R)
#Separate the groups - All tracts SLF FD (raincloud plot)
#R default colours for 5 groups (https://www.statology.org/ggplot-default-colors/):
#scale_color_manual(values = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")

group.colours <- c("1" = "#F8766D", "2" = "#A3A500", "3" ="#00BF7D", "4" = "#00B0F6")
#group.colours <- c("C" = "#F8766D", "SCD" = "#A3A500", "aMCI" ="#00BF7D", "mMCI" = "#00B0F6")
#version 1
ggplot(SLF_data_FD_noAD_long, aes(x = SLF_type, y = FD_metric, fill=Group)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  #geom_point(aes(y = FD_metric, colour = Group), size = .5, alpha = 0.8, position = position_dodge(width=0.5)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
  #scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI"))+
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fibre Density (FD)") +
  scale_x_discrete(labels = c("mn_FD_SLF1_L" = "Left SLF1","mn_FD_SLF2_L" = "Left SLF2","mn_FD_SLF3_L" = "Left SLF3","mn_FD_SLF1_R" = "Right SLF1","mn_FD_SLF2_R" = "Right SLF2","mn_FD_SLF3_R" = "Right SLF3")) + 
  scale_fill_manual(values=group.colours)+
  scale_color_manual(values = group.colours)+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()

#version 2
ggplot(SLF_data_FD_noAD_long, aes(x = SLF_type, y = FD_metric)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fibre Density (FD)") +
  scale_x_discrete(labels = c("mn_FD_SLF1_L" = "Left SLF1","mn_FD_SLF2_L" = "Left SLF2","mn_FD_SLF3_L" = "Left SLF3","mn_FD_SLF1_R" = "Right SLF1","mn_FD_SLF2_R" = "Right SLF2","mn_FD_SLF3_R" = "Right SLF3")) + 
  scale_color_manual(values = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()


#plot 95% Confidence Interval (separated confidence intervals and separated SLF tracts only)
SLF_FD_95CI_data_no_combined_SLF <- data.frame(SLF_group_number = c('1','1','2','3','3'),
                                                SLF_type = c('Left_SLF3','Left_SLF3','Right_SLF2','Right_SLF3', 'Right_SLF3'),
                                                Group_contrast = c('C > SCD','C > mMCI','C > mMCI','C > SCD','C > mMCI'),
                                                estimate_diff = c(0.0128785,0.0152952,0.019938,0.0140354,0.0185203),
                                                lower = c(0.0248848,0.0278578,0.038762,0.0271993,0.0322941), 
                                                upper = c(0.0008721,0.0027326,0.001114,0.0008715,0.0047465))  
#plot data
ggplot(SLF_FD_95CI_data_no_combined_SLF, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
  geom_point(position = position_dodge(width=0.5), size=5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), size=1.5, position = position_dodge(width=0.5))+
  xlab("FD") + 
  ylab("95% Confidence Interval")+
  scale_x_discrete(labels = c("1" = "Left SLF 3", "2" = "Right SLF 2", "3" = "Right SLF 3"))+
  scale_color_manual(values = c("#00B0F6", "#A3A500"))+
  labs(colour='Group Contrast')+   
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()


#plot group comparisons separately for sig. SLF tracts
#Left SLF 3 FD (raincloud plot)
ggplot(SLF_data_noAD, aes(x = Group, y = mn_FD_SLF3_L, fill = Group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mn_FD_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) +
  xlab("Group") +
  ylab("Left SLF 3 FD") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) +
  scale_fill_manual(values=group.colours)+
  scale_color_manual(values = group.colours)+
  theme_classic() +
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  theme(legend.position = "none") +
  coord_flip()
#Right SLF 3 FD (raincloud plot)
ggplot(SLF_data_noAD, aes(x = Group, y = mn_FD_SLF3_R, fill = Group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mn_FD_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) +
  xlab("Group") +
  ylab("Right SLF 3 FD") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) +
  scale_fill_manual(values=group.colours)+
  scale_color_manual(values = group.colours)+
  theme_classic() +
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  theme(legend.position = "none") +
  coord_flip()
#Right SLF 2 FD (raincloud plot)
ggplot(SLF_data_noAD, aes(x = Group, y = mn_FD_SLF2_R, fill = Group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mn_FD_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) +
  xlab("Group") +
  ylab("Right SLF 2 FD") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) +
  scale_fill_manual(values=group.colours)+
  scale_color_manual(values = group.colours)+
  theme_classic() +
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  theme(legend.position = "none") +
  coord_flip()







#Run contrast analysis (i.e., one-way ANOVA, b/c covariates) between C vs. AD group
#set up dataframe to only have C and AD group
#removing SCD, aMCI, and mMCI groups (group 2, 3, 4) from analysis
SLF_data_C_AD<- SLF_data %>%
  group_by(Group) %>%
  filter(!any(Group == "2")) %>%
  filter(!any(Group == "3")) %>%
  filter(!any(Group == "4"))
#change factor levels to '2' for calc. effect size
SLF_data_C_AD <- SLF_data_C_AD %>% mutate(Group = factor(Group, levels = c(1, 5)))


#descriptive stats for C and AD groups
#descriptive data for total sample (n = 62):
  #look at descriptive stats of the Left SLF1 FBA metrics between groups
  mean(SLF_data_C_AD$mn_FD_SLF1_L)
  sd(SLF_data_C_AD$mn_FD_SLF1_L)
  mean(SLF_data_C_AD$mn_FC_SLF1_L)
  sd(SLF_data_C_AD$mn_FC_SLF1_L)
  mean(SLF_data_C_AD$mn_FDC_SLF1_L)
  sd(SLF_data_C_AD$mn_FDC_SLF1_L)
  #look at descriptive stats of the Left SLF2 FBA metrics between groups
  mean(SLF_data_C_AD$mn_FD_SLF2_L)
  sd(SLF_data_C_AD$mn_FD_SLF2_L)
  mean(SLF_data_C_AD$mn_FC_SLF2_L)
  sd(SLF_data_C_AD$mn_FC_SLF2_L)
  mean(SLF_data_C_AD$mn_FDC_SLF2_L)
  sd(SLF_data_C_AD$mn_FDC_SLF2_L)
  #look at descriptive stats of the Left SLF3 FBA metrics between groups
  mean(SLF_data_C_AD$mn_FD_SLF3_L)
  sd(SLF_data_C_AD$mn_FD_SLF3_L)
  mean(SLF_data_C_AD$mn_FC_SLF3_L)
  sd(SLF_data_C_AD$mn_FC_SLF3_L)
  mean(SLF_data_C_AD$mn_FDC_SLF3_L)
  sd(SLF_data_C_AD$mn_FDC_SLF3_L)
  #look at descriptive stats of the Right SLF1 FBA metrics between groups
  mean(SLF_data_C_AD$mn_FD_SLF1_R)
  sd(SLF_data_C_AD$mn_FD_SLF1_R)
  mean(SLF_data_C_AD$mn_FC_SLF1_R)
  sd(SLF_data_C_AD$mn_FC_SLF1_R)
  mean(SLF_data_C_AD$mn_FDC_SLF1_R)
  sd(SLF_data_C_AD$mn_FDC_SLF1_R)
  #look at descriptive stats of the Right SLF2 FBA metrics between groups
  mean(SLF_data_C_AD$mn_FD_SLF2_R)
  sd(SLF_data_C_AD$mn_FD_SLF2_R)
  mean(SLF_data_C_AD$mn_FC_SLF2_R)
  sd(SLF_data_C_AD$mn_FC_SLF2_R)
  mean(SLF_data_C_AD$mn_FDC_SLF2_R)
  sd(SLF_data_C_AD$mn_FDC_SLF2_R)
  #look at descriptive stats of the Right SLF3 FBA metrics between groups
  mean(SLF_data_C_AD$mn_FD_SLF3_R)
  sd(SLF_data_C_AD$mn_FD_SLF3_R)
  mean(SLF_data_C_AD$mn_FC_SLF3_R)
  sd(SLF_data_C_AD$mn_FC_SLF3_R)
  mean(SLF_data_C_AD$mn_FDC_SLF3_R)
  sd(SLF_data_C_AD$mn_FDC_SLF3_R)

#SLF differences between C vs. AD group w/ age and sex as covariates ----------#
#ANOVA groups differences between 2 groups (C and AD)
  #FD
  SLF1_FD_2covar_mod_L_C_AD <- lm(mn_FD_SLF1_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF1_FD_2covar_mod_L_C_AD) #not sig
  SLF2_FD_2covar_mod_L_C_AD <- lm(mn_FD_SLF2_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF2_FD_2covar_mod_L_C_AD) #sig
  SLF3_FD_2covar_mod_L_C_AD <- lm(mn_FD_SLF3_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF3_FD_2covar_mod_L_C_AD) #sig
  SLF1_FD_2covar_mod_R_C_AD <- lm(mn_FD_SLF1_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF1_FD_2covar_mod_R_C_AD) #sig
  SLF2_FD_2covar_mod_R_C_AD <- lm(mn_FD_SLF2_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF2_FD_2covar_mod_R_C_AD) #sig
  SLF3_FD_2covar_mod_R_C_AD <- lm(mn_FD_SLF3_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF3_FD_2covar_mod_R_C_AD) #sig
  #FC
  SLF1_FC_2covar_mod_L_C_AD <- lm(mn_FC_SLF1_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF1_FC_2covar_mod_L_C_AD) #not sig
  SLF2_FC_2covar_mod_L_C_AD <- lm(mn_FC_SLF2_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF2_FC_2covar_mod_L_C_AD) #not sig
  SLF3_FC_2covar_mod_L_C_AD <- lm(mn_FC_SLF3_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF3_FC_2covar_mod_L_C_AD) #not sig
  SLF1_FC_2covar_mod_R_C_AD <- lm(mn_FC_SLF1_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF1_FC_2covar_mod_R_C_AD) #not sig
  SLF2_FC_2covar_mod_R_C_AD <- lm(mn_FC_SLF2_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF2_FC_2covar_mod_R_C_AD) #not sig
  SLF3_FC_2covar_mod_R_C_AD <- lm(mn_FC_SLF3_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF3_FC_2covar_mod_R_C_AD) #not sig
  #FDC
  SLF1_FDC_2covar_mod_L_C_AD <- lm(mn_FDC_SLF1_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF1_FDC_2covar_mod_L_C_AD) #not sig
  SLF2_FDC_2covar_mod_L_C_AD <- lm(mn_FDC_SLF2_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF2_FDC_2covar_mod_L_C_AD) #sig
  SLF3_FDC_2covar_mod_L_C_AD <- lm(mn_FDC_SLF3_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF3_FDC_2covar_mod_L_C_AD) #not sig
  SLF1_FDC_2covar_mod_R_C_AD <- lm(mn_FDC_SLF1_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF1_FDC_2covar_mod_R_C_AD) #sig
  SLF2_FDC_2covar_mod_R_C_AD <- lm(mn_FDC_SLF2_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF2_FDC_2covar_mod_R_C_AD) #sig
  SLF3_FDC_2covar_mod_R_C_AD <- lm(mn_FDC_SLF3_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF3_FDC_2covar_mod_R_C_AD) #sig

#calculate the effect size (eta-squared) - for ANCOVA (sex & age)
  #FD
  etaSquared(SLF1_FD_2covar_mod_L_C_AD)
  etaSquared(SLF2_FD_2covar_mod_L_C_AD)
  etaSquared(SLF3_FD_2covar_mod_L_C_AD)
  etaSquared(SLF1_FD_2covar_mod_R_C_AD)
  etaSquared(SLF2_FD_2covar_mod_R_C_AD)
  etaSquared(SLF3_FD_2covar_mod_R_C_AD)
  #FC
  etaSquared(SLF1_FC_2covar_mod_L_C_AD)
  etaSquared(SLF2_FC_2covar_mod_L_C_AD)
  etaSquared(SLF3_FC_2covar_mod_L_C_AD)
  etaSquared(SLF1_FC_2covar_mod_R_C_AD)
  etaSquared(SLF2_FC_2covar_mod_R_C_AD)
  etaSquared(SLF3_FC_2covar_mod_R_C_AD)
  #FDC
  etaSquared(SLF1_FDC_2covar_mod_L_C_AD)
  etaSquared(SLF2_FDC_2covar_mod_L_C_AD)
  etaSquared(SLF3_FDC_2covar_mod_L_C_AD)
  etaSquared(SLF1_FDC_2covar_mod_R_C_AD)
  etaSquared(SLF2_FDC_2covar_mod_R_C_AD)
  etaSquared(SLF3_FDC_2covar_mod_R_C_AD)

#calculate the 95CI for C vs. AD in FD
  #for Left SLF 2
  SLF2_FD_2covar_mod_L_C_AD <- lm(mn_FD_SLF2_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF2_FD_2covar_mod_L_C_AD) #sig with the variable 'Group' has a significant effect overall,
  summary(SLF2_FD_2covar_mod_L_C_AD)  #not sig with the specific contrast between Group5 (AD) and the reference (C)
  confint(SLF2_FD_2covar_mod_L_C_AD) #95%CI passes over zero
  
  #for Left SLF 3
  SLF3_FD_2covar_mod_L_C_AD <- lm(mn_FD_SLF3_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF3_FD_2covar_mod_L_C_AD) #sig
  summary(SLF3_FD_2covar_mod_L_C_AD)  #sig 
  confint(SLF3_FD_2covar_mod_L_C_AD) 
  
  #for Right SLF 1
  SLF1_FD_2covar_mod_R_C_AD <- lm(mn_FD_SLF1_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF1_FD_2covar_mod_R_C_AD) #sig
  summary(SLF1_FD_2covar_mod_R_C_AD) #sig
  confint(SLF1_FD_2covar_mod_R_C_AD)
  
  #for Right SLF 2
  SLF2_FD_2covar_mod_R_C_AD <- lm(mn_FD_SLF2_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF2_FD_2covar_mod_R_C_AD) #sig
  summary(SLF2_FD_2covar_mod_R_C_AD) #sig
  confint(SLF2_FD_2covar_mod_R_C_AD)
  
  #for Right SLF 3
  SLF3_FD_2covar_mod_R_C_AD <- lm(mn_FD_SLF3_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF3_FD_2covar_mod_R_C_AD) #sig
  summary(SLF3_FD_2covar_mod_R_C_AD) #sig
  confint(SLF3_FD_2covar_mod_R_C_AD)

#generate figures for C vs. AD comparison
  #Plotting graphs via ggplot2
  #for FD: just with left and right SLF 1, 2, 3
  SLF_data_FD_C_AD <- dplyr::select(SLF_data_C_AD, 
                               ParticipantID,
                               Group,
                               mn_FD_SLF1_L,
                               mn_FD_SLF2_L,
                               mn_FD_SLF3_L,
                               mn_FD_SLF1_R,
                               mn_FD_SLF2_R,
                               mn_FD_SLF3_R)
  SLF_data_FD_C_AD_long <- gather(SLF_data_C_AD, 
                                  "SLF_type",
                                  "FD_metric",
                                  mn_FD_SLF1_L,
                                  mn_FD_SLF2_L,
                                  mn_FD_SLF3_L,
                                  mn_FD_SLF1_R,
                                  mn_FD_SLF2_R,
                                  mn_FD_SLF3_R)
  #Separate the groups - All tracts SLF FD (raincloud plot)
  #R default colours for 5 groups (https://www.statology.org/ggplot-default-colors/):
  #scale_color_manual(values = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")
  
  group.colours <- c("1" = "#F8766D", "5" = "#E76BF3")
  #group.colours <- c("C" = "#F8766D", "SCD" = "#A3A500", "aMCI" ="#00BF7D", "mMCI" = "#00B0F6","AD" = "#E76BF3")
  #version 1
  ggplot(SLF_data_FD_C_AD_long, aes(x = SLF_type, y = FD_metric, fill=Group)) + 
    geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
    #geom_point(aes(y = FD_metric, colour = Group), size = .5, alpha = 0.8, position = position_dodge(width=0.5)) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
    #scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI","5"="AD"))+
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
    guides(colour=FALSE)+
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF Tract") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("mn_FD_SLF1_L" = "Left SLF1","mn_FD_SLF2_L" = "Left SLF2","mn_FD_SLF3_L" = "Left SLF3","mn_FD_SLF1_R" = "Right SLF1","mn_FD_SLF2_R" = "Right SLF2","mn_FD_SLF3_R" = "Right SLF3")) + 
    scale_fill_manual(values=group.colours, labels = c("1" = "C", "5" = "AD"))+
    scale_color_manual(values = group.colours, labels = c("1" = "C", "5" = "AD"))+    #scale_fill_manual(values=c("1"="#F8766D","5"="#E76BF3"))+
    theme_classic()+
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18)+
    #theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18), legend.position="none")+ #if you want to remove the legend
    coord_flip()
  
  #version 2
  ggplot(SLF_data_FD_C_AD_long, aes(x = SLF_type, y = FD_metric)) + 
    geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
    guides(colour=FALSE)+
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF Tract") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("mn_FD_SLF1_L" = "Left SLF1","mn_FD_SLF2_L" = "Left SLF2","mn_FD_SLF3_L" = "Left SLF3","mn_FD_SLF1_R" = "Right SLF1","mn_FD_SLF2_R" = "Right SLF2","mn_FD_SLF3_R" = "Right SLF3")) + 
    scale_color_manual(values = c("#F8766D", "#E76BF3"))+
    theme_classic()+
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    coord_flip()

  #plot 95% Confidence Interval (separated confidence intervals and separated SLF tracts only)
  SLF_FD_95CI_data_no_combined_SLF_CvAD <- data.frame(SLF_group_number = c('1','2','3','4','5'),
                                                 SLF_type = c('Left_SLF2','Left_SLF3','Right_SLF1','Right_SLF2', 'Right_SLF3'),
                                                 Group_contrast = c('C > AD','C > AD','C > AD','C > AD','C > AD'),
                                                 estimate_diff = c(0.0162383,0.0166872,0.0229439,2.094e-02,0.0148583),
                                                 lower = c(0.034974319,0.029468734,0.038086372,0.0372980340,0.0280043945), 
                                                 upper = c(-2.497785e-03,0.003905601,0.0078013929,0.0045791193,0.0017121763))  
  #plot data
  ggplot(SLF_FD_95CI_data_no_combined_SLF_CvAD, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point(position = position_dodge(width=0.5), size=5)+
    geom_errorbar(aes(ymin=lower, ymax=upper), size=1.5, position = position_dodge(width=0.5))+
    xlab("FD") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Left SLF 2", "2" = "Left SLF 3","3" = "Right SLF 1", "4" = "Right SLF 2", "5" = "Right SLF 3"))+
    scale_color_manual(values = "#E76BF3")+
    labs(colour='Group Contrast')+   
    theme_classic()+
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18), plot.margin = margin(r = 30), legend.position="none")+ #if you want to remove the legend)
    #theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    coord_flip()
  

#plot group comparisons separately for sig. SLF tracts
  #Left SLF 2 FD (raincloud plot)
  ggplot(SLF_data_C_AD, aes(x = Group, y = mn_FD_SLF2_L, fill = Group)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) +
    xlab("Group") +
    ylab("Left SLF 2 FD") +
    scale_x_discrete(labels = c("1" = "Control", "5" = "AD")) +
    scale_fill_manual(values = c("#F8766D", "#E76BF3"))+
    scale_color_manual(values = c("#F8766D", "#E76BF3"))+
    theme_classic() +
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    theme(legend.position = "none") +
    coord_flip()

  #Left SLF 3 FD (raincloud plot)
  ggplot(SLF_data_C_AD, aes(x = Group, y = mn_FD_SLF3_L, fill = Group)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) +
    xlab("Group") +
    ylab("Left SLF 3 FD") +
    scale_x_discrete(labels = c("1" = "Control", "5" = "AD")) +
    scale_fill_manual(values = c("#F8766D", "#E76BF3"))+
    scale_color_manual(values = c("#F8766D", "#E76BF3"))+
    theme_classic() +
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    theme(legend.position = "none") +
    coord_flip()
  
  #Right SLF 1 FD (raincloud plot)
  ggplot(SLF_data_C_AD, aes(x = Group, y = mn_FD_SLF1_R, fill = Group)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) +
    xlab("Group") +
    ylab("Right SLF 1 FD") +
    scale_x_discrete(labels = c("1" = "Control", "5" = "AD")) +
    scale_fill_manual(values = c("#F8766D", "#E76BF3"))+
    scale_color_manual(values = c("#F8766D", "#E76BF3"))+
    theme_classic() +
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    theme(legend.position = "none") +
    coord_flip()
  
  #Right SLF 2 FD (raincloud plot)
  ggplot(SLF_data_C_AD, aes(x = Group, y = mn_FD_SLF2_R, fill = Group)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) +
    xlab("Group") +
    ylab("Right SLF 2 FD") +
    scale_x_discrete(labels = c("1" = "Control", "5" = "AD")) +
    scale_fill_manual(values = c("#F8766D", "#E76BF3"))+
    scale_color_manual(values = c("#F8766D", "#E76BF3"))+
    theme_classic() +
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    theme(legend.position = "none") +
    coord_flip()
  
  #Right SLF 3 FD (raincloud plot)
  ggplot(SLF_data_C_AD, aes(x = Group, y = mn_FD_SLF3_R, fill = Group)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mn_FD_SLF, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) +
    xlab("Group") +
    ylab("Right SLF 3 FD") +
    scale_x_discrete(labels = c("1" = "Control", "5" = "AD")) +
    scale_fill_manual(values = c("#F8766D", "#E76BF3"))+
    scale_color_manual(values = c("#F8766D", "#E76BF3"))+
    theme_classic() +
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    theme(legend.position = "none") +
    coord_flip()
  
  
  
  #calculate the 95CI for C vs. AD in FDC
  #for Left SLF 2
  SLF2_FDC_2covar_mod_L_C_AD <- lm(mn_FDC_SLF2_L ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF2_FDC_2covar_mod_L_C_AD) #sig with the variable 'Group' has a significant effect overall,
  summary(SLF2_FDC_2covar_mod_L_C_AD)  #sig with the specific contrast between Group5 (AD) and the reference (C)
  confint(SLF2_FDC_2covar_mod_L_C_AD) #95%CI passes over zero

  #for Right SLF 1
  SLF1_FDC_2covar_mod_R_C_AD <- lm(mn_FDC_SLF1_R ~ Group + Age + Sex, data = SLF_data_C_AD)
  anova(SLF1_FDC_2covar_mod_R_C_AD) #sig
  summary(SLF1_FDC_2covar_mod_R_C_AD) #sig
  confint(SLF1_FDC_2covar_mod_R_C_AD)
  
  #generate figures for C vs. AD comparison
  #Plotting graphs via ggplot2
  #for FDC: just with left SLF 2 and right SLF 1
  SLF_data_FDC_C_AD <- dplyr::select(SLF_data_C_AD, 
                               ParticipantID,
                               Group,
                               mn_FDC_SLF1_L,
                               mn_FDC_SLF2_L,
                               mn_FDC_SLF3_L,
                               mn_FDC_SLF1_R,
                               mn_FDC_SLF2_R,
                               mn_FDC_SLF3_R)
  SLF_data_FDC_C_AD_long <- gather(SLF_data_C_AD, 
                                  "SLF_type",
                                  "FDC_metric",
                                  mn_FDC_SLF1_L,
                                  mn_FDC_SLF2_L,
                                  mn_FDC_SLF3_L,
                                  mn_FDC_SLF1_R,
                                  mn_FDC_SLF2_R,
                                  mn_FDC_SLF3_R)
  #Separate the groups - All tracts SLF FD (raincloud plot)
  #R default colours for 5 groups (https://www.statology.org/ggplot-default-colors/):
  #scale_color_manual(values = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")
  
  group.colours <- c("1" = "#F8766D", "5" = "#E76BF3")
  #group.colours <- c("C" = "#F8766D", "SCD" = "#A3A500", "aMCI" ="#00BF7D", "mMCI" = "#00B0F6","AD" = "#E76BF3")
  #version 1
  ggplot(SLF_data_FDC_C_AD_long, aes(x = SLF_type, y = FDC_metric, fill=Group)) + 
    geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
    guides(colour=FALSE)+
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF Tract") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("mn_FDC_SLF1_L" = "Left SLF1","mn_FDC_SLF2_L" = "Left SLF2","mn_FDC_SLF3_L" = "Left SLF3","mn_FDC_SLF1_R" = "Right SLF1","mn_FDC_SLF2_R" = "Right SLF2","mn_FDC_SLF3_R" = "Right SLF3")) + 
    scale_fill_manual(values=group.colours, labels = c("1" = "C", "5" = "AD"))+
    scale_color_manual(values = group.colours, labels = c("1" = "C", "5" = "AD"))+    #scale_fill_manual(values=c("1"="#F8766D","5"="#E76BF3"))+
    theme_classic()+
    #theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18)+
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18), legend.position="none")+ #if you want to remove the legend
    coord_flip()
  
  #version 2
  ggplot(SLF_data_FDC_C_AD_long, aes(x = SLF_type, y = FDC_metric)) + 
    geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
    guides(colour=FALSE)+
    guides(fill = guide_legend(override.aes = list(shape = NA)))+
    xlab("SLF Tract") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("mn_FDC_SLF1_L" = "Left SLF1","mn_FDC_SLF2_L" = "Left SLF2","mn_FDC_SLF3_L" = "Left SLF3","mn_FDC_SLF1_R" = "Right SLF1","mn_FDC_SLF2_R" = "Right SLF2","mn_FDC_SLF3_R" = "Right SLF3")) + 
    scale_color_manual(values = c("#F8766D", "#E76BF3"))+
    theme_classic()+
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    coord_flip()
  
  #plot 95% Confidence Interval (separated confidence intervals and separated SLF tracts only)
  SLF_FDC_95CI_data_no_combined_SLF_CvAD <- data.frame(SLF_group_number = c('1','2'),
                                                      SLF_type = c('Left_SLF2','Right_SLF1'),
                                                      Group_contrast = c('C > AD','C > AD'),
                                                      estimate_diff = c(0.0447525,0.048869),
                                                      lower = c(0.076904842,0.074221328), 
                                                      upper = c(0.0126001753,2.351607e-02))  
  #plot data
  ggplot(SLF_FDC_95CI_data_no_combined_SLF_CvAD, aes(x=SLF_group_number, y=estimate_diff, colour=Group_contrast))+
    geom_point(position = position_dodge(width=0.5), size=5)+
    geom_errorbar(aes(ymin=lower, ymax=upper), size=1.5, position = position_dodge(width=0.5))+
    xlab("FD") + 
    ylab("95% Confidence Interval")+
    scale_x_discrete(labels = c("1" = "Left SLF 2", "2" = "Left SLF 3","3" = "Right SLF 1", "4" = "Right SLF 2", "5" = "Right SLF 3"))+
    scale_color_manual(values = "#E76BF3")+
    labs(colour='Group Contrast')+   
    theme_classic()+
    theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18), plot.margin = margin(r = 30), legend.position="none")+ #if you want to remove the legend)
    #theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
    coord_flip()
   

