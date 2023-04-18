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
pacman::p_load(dplyr, ggplot2, psych, car, multcomp, lsr, tidyr, BayesFactor, tidyverse, ppcor, nlme, effectsize, rstatix, sjstats, purrr, corrplot, compute.es)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph
source("H:/ltah262/PhD/Diffusion/script/dprc/diffusion/visualisation/confidence_interval.R")

#first read in the covariates group data file: 
#setwd('/yourpathway/')
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')


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
setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/template/TOI')

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
#95% CI
diff_SLF1_L_FC_data <-  SLF_data[SLF_data[, "Timepoint"] == "F2",]$mn_FC_SLF1_L -  SLF_data[SLF_data[, "Timepoint"] == "F0",]$mn_FC_SLF1_L 
mean(diff_SLF1_L_FC_data)
confidence_interval(diff_SLF1_L_FC_data,0.95)

aov_SLF1_L_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF1_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_L_FDC_mod)
#95% CI
diff_SLF1_L_FDC_data <-  SLF_data[SLF_data[, "Timepoint"] == "F2",]$mn_FDC_SLF1_L -  SLF_data[SLF_data[, "Timepoint"] == "F0",]$mn_FDC_SLF1_L 
mean(diff_SLF1_L_FDC_data)
confidence_interval(diff_SLF1_L_FDC_data,0.95)

#for left SLF 2
aov_SLF2_L_FD_mod <- anova_test(data=SLF_data, dv=mn_FD_SLF2_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_FD_mod) 

aov_SLF2_L_FC_mod <- anova_test(data=SLF_data, dv=mn_FC_SLF2_L, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_FC_mod) 
#95% CI
diff_SLF2_L_FC_data <-  SLF_data[SLF_data[, "Timepoint"] == "F2",]$mn_FC_SLF2_L -  SLF_data[SLF_data[, "Timepoint"] == "F0",]$mn_FC_SLF2_L 
mean(diff_SLF2_L_FC_data)
confidence_interval(diff_SLF2_L_FC_data,0.95)

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
#95% CI
diff_SLF1_R_FC_data <-  SLF_data[SLF_data[, "Timepoint"] == "F2",]$mn_FC_SLF1_R -  SLF_data[SLF_data[, "Timepoint"] == "F0",]$mn_FC_SLF1_R 
mean(diff_SLF1_R_FC_data)
confidence_interval(diff_SLF1_R_FC_data,0.95)

aov_SLF1_R_FDC_mod <- anova_test(data=SLF_data, dv=mn_FDC_SLF1_R, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_R_FDC_mod)
#95% CI
diff_SLF1_R_FDC_data <-  SLF_data[SLF_data[, "Timepoint"] == "F2",]$mn_FDC_SLF1_R -  SLF_data[SLF_data[, "Timepoint"] == "F0",]$mn_FDC_SLF1_R 
mean(diff_SLF1_R_FDC_data)
confidence_interval(diff_SLF1_R_FDC_data,0.95)

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



#FC 95% CI
#plot 95% Confidence Interval (separated confidence intervals and separated SLF tracts only)
SLF_FC_95CI_data_no_combined_SLF <- data.frame(SLF_group_number = c('1','2','3'),
                                               SLF_type = c('Left_SLF1','Left_SLF2','Right_SLF1'),
                                               Timepoint_contrast = c('F0 > F2','F0 > F2','F0 > F2'),
                                               estimate_diff = c(-0.007329525,-0.00368148,-0.007308974)*-1,
                                               lower = c(-0.010635067,-0.0070563139,-0.010512332)*-1, 
                                               upper = c(-0.004023983,-0.0003066471,-0.004105617)*-1)  


# #FC 95% CI (with sex and age as covariates) #does not plot correctly.....
# #plot 95% Confidence Interval (separated confidence intervals and separated SLF tracts only)
# SLF_FC_95CI_data_no_combined_SLF <- data.frame(SLF_group_number = c('1','2','3'),
#                                                SLF_type = c('Left_SLF1','Left_SLF2','Right_SLF1'),
#                                                Timepoint_contrast = c('F0 > F2','F0 > F2','F0 > F2'),
#                                                estimate_diff = c(-0.007329525,-0.00368148,-0.007308974)*-1,#diff. between the means
#                                                lower = c(-0.019460800,-0.0169412925,-0.0204636577)*-1,
#                                                upper = c(-0.0007652705,0.001953156,-0.002335503)*-1)  #why is the second one positive??? not sig. then



#plot data
ggplot(SLF_FC_95CI_data_no_combined_SLF, aes(x=SLF_group_number, y=estimate_diff, colour=Timepoint_contrast))+
  geom_point(position = position_dodge(width=0.5), size=5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), size=1.5, position = position_dodge(width=0.5))+
  xlab("SLF Tract") + 
  ylab("95% Confidence Interval")+
  scale_x_discrete(labels = c("1" = "Left SLF 1", "2" = "Left SLF 2", "3" = "Right SLF 1"))+
  scale_color_manual(values = ("#999999"))+
  labs(colour='Time Point Contrast')+   
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()


#FDC 95% CI
#plot 95% Confidence Interval (separated confidence intervals and separated SLF tracts only)
SLF_FDC_95CI_data_no_combined_SLF <- data.frame(SLF_group_number = c('1','2'),
                                               SLF_type = c('Left_SLF1','Right_SLF1'),
                                               Timepoint_contrast = c('F0 > F2','F0 > F2'),
                                               estimate_diff = c(-0.005480806,-0.004674242)*-1,
                                               lower = c(-0.009500077,-0.0085991466)*-1, 
                                               upper = c(-0.001461536,-0.0007493373)*-1)  
#plot data
ggplot(SLF_FDC_95CI_data_no_combined_SLF, aes(x=SLF_group_number, y=estimate_diff, colour=Timepoint_contrast))+
  geom_point(position = position_dodge(width=0.5), size=5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), size = 1.5, position = position_dodge(width=0.5))+
  xlab("SLF Tract") + 
  ylab("95% Confidence Interval")+
  scale_x_discrete(labels = c("1" = "Left SLF 1", "2" = "Right SLF 1"))+
  scale_color_manual(values = ("#999999"))+
  labs(colour='Time Point Contrast')+   
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()

#Run ANCOVA With covariates (age)
#Left SLF 1
aov_SLF1_L_FD_covar_mod<- aov(mn_FD_SLF1_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF1_L_FD_covar_mod)
aov_SLF1_L_FC_covar_mod<- aov(mn_FC_SLF1_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF1_L_FC_covar_mod)
aov_SLF1_L_FDC_covar_mod<- aov(mn_FDC_SLF1_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF1_L_FDC_covar_mod)
#Left SLF 2
aov_SLF2_L_FD_covar_mod<- aov(mn_FD_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF2_L_FD_covar_mod) 
aov_SLF2_L_FC_covar_mod<- aov(mn_FC_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF2_L_FC_covar_mod)
aov_SLF2_L_FDC_covar_mod<- aov(mn_FDC_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF2_L_FDC_covar_mod) 
#Left SLF 3
aov_SLF3_L_FD_covar_mod<- aov(mn_FD_SLF3_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF3_L_FD_covar_mod)
aov_SLF3_L_FC_covar_mod<- aov(mn_FC_SLF3_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF3_L_FC_covar_mod)
aov_SLF3_L_FDC_covar_mod<- aov(mn_FDC_SLF3_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF3_L_FDC_covar_mod)
#Right SLF 1
aov_SLF1_R_FD_covar_mod<- aov(mn_FD_SLF1_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF1_R_FD_covar_mod)
aov_SLF1_R_FC_covar_mod<- aov(mn_FC_SLF1_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF1_R_FC_covar_mod)
aov_SLF1_R_FDC_covar_mod<- aov(mn_FDC_SLF1_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF1_R_FDC_covar_mod)
#Right SLF 2
aov_SLF2_R_FD_covar_mod<- aov(mn_FD_SLF2_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF2_R_FD_covar_mod) 
aov_SLF2_R_FC_covar_mod<- aov(mn_FC_SLF2_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF2_R_FC_covar_mod)
aov_SLF2_R_FDC_covar_mod<- aov(mn_FDC_SLF2_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF2_R_FDC_covar_mod) 
#Right SLF 3
aov_SLF3_R_FD_covar_mod<- aov(mn_FD_SLF3_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF3_R_FD_covar_mod)
aov_SLF3_R_FC_covar_mod<- aov(mn_FC_SLF3_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF3_R_FC_covar_mod)
aov_SLF3_R_FDC_covar_mod<- aov(mn_FDC_SLF3_R ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_data)
summary(aov_SLF3_R_FDC_covar_mod)
#effect size w/ covariates
sjstats::eta_sq(aov_SLF1_L_FD_covar_mod)
sjstats::eta_sq(aov_SLF1_L_FC_covar_mod)
sjstats::eta_sq(aov_SLF1_L_FDC_covar_mod)
sjstats::eta_sq(aov_SLF2_L_FD_covar_mod)
sjstats::eta_sq(aov_SLF2_L_FC_covar_mod)
sjstats::eta_sq(aov_SLF2_L_FDC_covar_mod)
sjstats::eta_sq(aov_SLF3_L_FD_covar_mod)
sjstats::eta_sq(aov_SLF3_L_FC_covar_mod)
sjstats::eta_sq(aov_SLF3_L_FDC_covar_mod)
sjstats::eta_sq(aov_SLF1_R_FD_covar_mod)
sjstats::eta_sq(aov_SLF1_R_FC_covar_mod)
sjstats::eta_sq(aov_SLF1_R_FDC_covar_mod)
sjstats::eta_sq(aov_SLF2_R_FD_covar_mod)
sjstats::eta_sq(aov_SLF2_R_FC_covar_mod)
sjstats::eta_sq(aov_SLF2_R_FDC_covar_mod)
sjstats::eta_sq(aov_SLF3_R_FD_covar_mod)
sjstats::eta_sq(aov_SLF3_R_FC_covar_mod)
sjstats::eta_sq(aov_SLF3_R_FDC_covar_mod)

#follow up tests with sig. main effects and interactions
#Left SLF1 FC
# post_hoc_aov_SLF1_L_FC_covar_mod <- lme(mn_FC_SLF1_L ~ Group*Timepoint+Age, random = ~1 | Individual_number/Timepoint, data=SLF_data)
# summary(glht(post_hoc_aov_SLF1_L_FC_covar_mod, linfct=mcp(Group="Tukey")))
# summary(glht(post_hoc_aov_SLF1_L_FC_covar_mod, linfct=mcp(Timepoint="Tukey")))
# #Left SLF1 FDC
# post_hoc_aov_SLF1_L_FDC_covar_mod <- lme(mn_FDC_SLF1_L ~ Group*Timepoint+Age, random = ~1 | Individual_number/Timepoint, data=SLF_data)
# summary(glht(post_hoc_aov_SLF1_L_FDC_covar_mod, linfct=mcp(Timepoint="Tukey")))
#Left SLF 2 FD
post_hoc_aov_SLF2_L_FD_covar_mod <- lme(mn_FD_SLF2_L ~ Group*Timepoint + Age, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_L_FD_covar_mod, linfct=mcp(Group="Tukey")))
# #Left SLF 2 FC
# post_hoc_aov_SLF2_L_FC_covar_mod <- lme(mn_FC_SLF2_L ~ Group*Timepoint + Age, random = ~1 | Individual_number/Timepoint, data=SLF_data)
# summary(glht(post_hoc_aov_SLF2_L_FC_covar_mod, linfct=mcp(Timepoint="Tukey")))
#Right SLF1 FD
post_hoc_aov_SLF1_R_FD_covar_mod <- lme(mn_FD_SLF1_R ~ Group*Timepoint+Age, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_R_FD_covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_SLF1_R_FD_covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. SCD
  SLF_data_CvSCD_SLF1_R_FD_covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
  SLF_data_CvSCD_SLF1_R_FD_covar$Group <- droplevels(SLF_data_CvSCD_SLF1_R_FD_covar$Group)
  group_number <-dplyr::count(SLF_data_CvSCD_SLF1_R_FD_covar, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF1_R ~ Age + Sex, data = SLF_data_CvSCD_SLF1_R_FD_covar)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. AD
  SLF_data_CvAD_SLF1_R_FD_covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 5)
  SLF_data_CvAD_SLF1_R_FD_covar$Group <- droplevels(SLF_data_CvAD_SLF1_R_FD_covar$Group)
  group_number <-dplyr::count(SLF_data_CvAD_SLF1_R_FD_covar, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF1_R ~ Age + Sex, data = SLF_data_CvAD_SLF1_R_FD_covar)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
#Right SLF2 FD
post_hoc_aov_SLF2_R_FD_covar_mod <- lme(mn_FD_SLF2_R ~ Group*Timepoint+Age, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_R_FD_covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_SLF2_R_FD_covar_mod, linfct=mcp(Group="Tukey")))
  #calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. SCD
  SLF_data_CvSCD_SLF2_R_FD_covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
  SLF_data_CvSCD_SLF2_R_FD_covar$Group <- droplevels(SLF_data_CvSCD_SLF2_R_FD_covar$Group)
  group_number <-dplyr::count(SLF_data_CvSCD_SLF2_R_FD_covar, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF2_R ~ Age, data = SLF_data_CvSCD_SLF2_R_FD_covar)) #create multiple regression between age, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age
#Right SLF3 FD
post_hoc_aov_SLF3_R_FD_covar_mod <- lme(mn_FD_SLF3_R ~ Group*Timepoint+Age, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF3_R_FD_covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_SLF3_R_FD_covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. SCD
  SLF_data_CvSCD_SLF3_R_FD_covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
  SLF_data_CvSCD_SLF3_R_FD_covar$Group <- droplevels(SLF_data_CvSCD_SLF3_R_FD_covar$Group)
  group_number <-dplyr::count(SLF_data_CvSCD_SLF3_R_FD_covar, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FD_SLF3_R ~ Age, data = SLF_data_CvSCD_SLF3_R_FD_covar)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age 
#Right SLF3 FDC
post_hoc_aov_SLF3_R_FDC_covar_mod <- lme(mn_FDC_SLF3_R ~ Group*Timepoint+Age, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF3_R_FDC_covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_SLF3_R_FDC_covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. SCD
  SLF_data_CvSCD_SLF3_R_FDC_covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
  SLF_data_CvSCD_SLF3_R_FDC_covar$Group <- droplevels(SLF_data_CvSCD_SLF3_R_FDC_covar$Group)
  group_number <-dplyr::count(SLF_data_CvSCD_SLF3_R_FDC_covar, Group) #count number of participants per group
  mult.r_value_covar_mod<-summary(lm(mn_FDC_SLF3_R ~ Age, data = SLF_data_CvSCD_SLF3_R_FDC_covar)) #create multiple regression between age and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=1) #calculate Cohen's D with the covariate of age 
  
  

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
#95%CI
confint(aov_SLF1_L_FC_2covar_mod$`Individual_number:Timepoint`)
aov_SLF1_L_FDC_2covar_mod<- aov(mn_FDC_SLF1_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF1_L_FDC_2covar_mod)
#Left SLF 2
aov_SLF2_L_FD_2covar_mod<- aov(mn_FD_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_L_FD_2covar_mod) 
aov_SLF2_L_FC_2covar_mod<- aov(mn_FC_SLF2_L ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data=SLF_data)
summary(aov_SLF2_L_FC_2covar_mod)
#95%CI
confint(aov_SLF2_L_FC_2covar_mod$`Individual_number:Timepoint`)
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
#95%CI
confint(aov_SLF1_R_FC_2covar_mod$`Individual_number:Timepoint`)
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
t_value_effect_size <- summary(glht(post_hoc_aov_SLF_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. SCD
  SLF_data_CvSCD_SLF_R_FD_2covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
  SLF_data_CvSCD_SLF_R_FD_2covar$Group <- droplevels(SLF_data_CvSCD_SLF_R_FD_2covar$Group)
  group_number <-dplyr::count(SLF_data_CvSCD_SLF_R_FD_2covar, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF_R ~ Age + Sex, data = SLF_data_CvSCD_SLF_R_FD_2covar)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex

#Right SLF FC
post_hoc_aov_SLF_R_FC_2covar_mod <- lme(mn_FC_SLF_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF_R_FC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Left SLF1 FC
post_hoc_aov_SLF1_L_FC_2covar_mod <- lme(mn_FC_SLF1_L ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_L_FC_2covar_mod, linfct=mcp(Group="Tukey"))) #not sig.
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
t_value_effect_size <- summary(glht(post_hoc_aov_SLF1_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
#for C vs. SCD
SLF_data_CvSCD_SLF1_R_FD_2covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
SLF_data_CvSCD_SLF1_R_FD_2covar$Group <- droplevels(SLF_data_CvSCD_SLF1_R_FD_2covar$Group)
group_number <-dplyr::count(SLF_data_CvSCD_SLF1_R_FD_2covar, Group) #count number of participants per group
mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF1_R ~ Age + Sex, data = SLF_data_CvSCD_SLF1_R_FD_2covar)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex

#Right SLF1 FC
post_hoc_aov_SLF1_R_FC_2covar_mod <- lme(mn_FC_SLF1_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_R_FC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Right SLF1 FDC
post_hoc_aov_SLF1_R_FDC_2covar_mod <- lme(mn_FDC_SLF1_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF1_R_FDC_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Right SLF2 FD
post_hoc_aov_SLF2_R_FD_2covar_mod <- lme(mn_FD_SLF2_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_SLF2_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
#for C vs. SCD
SLF_data_CvSCD_SLF2_R_FD_2covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
SLF_data_CvSCD_SLF2_R_FD_2covar$Group <- droplevels(SLF_data_CvSCD_SLF2_R_FD_2covar$Group)
group_number <-dplyr::count(SLF_data_CvSCD_SLF2_R_FD_2covar, Group) #count number of participants per group
mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF2_R ~ Age + Sex, data = SLF_data_CvSCD_SLF2_R_FD_2covar)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex

#Right SLF2 FC
post_hoc_aov_SLF2_R_FC_2covar_mod <- lme(mn_FC_SLF2_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF2_R_FC_2covar_mod, linfct=mcp(Group="Tukey"))) #not sig.
#Right SLF3 FD
post_hoc_aov_SLF3_R_FD_2covar_mod <- lme(mn_FD_SLF3_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF3_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_SLF3_R_FD_2covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
#for C vs. SCD
SLF_data_CvSCD_SLF3_R_FD_2covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
SLF_data_CvSCD_SLF3_R_FD_2covar$Group <- droplevels(SLF_data_CvSCD_SLF3_R_FD_2covar$Group)
group_number <-dplyr::count(SLF_data_CvSCD_SLF3_R_FD_2covar, Group) #count number of participants per group
mult.r_value_2covar_mod<-summary(lm(mn_FD_SLF3_R ~ Age + Sex, data = SLF_data_CvSCD_SLF3_R_FD_2covar)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex

#Right SLF3 FDC
post_hoc_aov_SLF3_R_FDC_2covar_mod <- lme(mn_FDC_SLF3_R ~ Group*Timepoint+Age+Sex, random = ~1 | Individual_number/Timepoint, data=SLF_data)
summary(glht(post_hoc_aov_SLF3_R_FDC_2covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_SLF3_R_FDC_2covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
#for C vs. SCD
SLF_data_CvSCD_SLF3_R_FDC_2covar <- subset(SLF_data, SLF_data$Group == 1 | SLF_data$Group == 2)
SLF_data_CvSCD_SLF3_R_FDC_2covar$Group <- droplevels(SLF_data_CvSCD_SLF3_R_FDC_2covar$Group)
group_number <-dplyr::count(SLF_data_CvSCD_SLF3_R_FDC_2covar, Group) #count number of participants per group
mult.r_value_2covar_mod<-summary(lm(mn_FDC_SLF3_R ~ Age + Sex, data = SLF_data_CvSCD_SLF3_R_FDC_2covar)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex


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

# #Combine all sig. plots on one plot for FC
# #put data into long format
# SLF_data_FC <- dplyr::select(SLF_data, 
#                              ParticipantID,
#                              Group,
#                              Timepoint,
#                              mn_FC_SLF,
#                              mn_FC_SLF_L,
#                              mn_FC_SLF_R,
#                              mn_FC_SLF1_L,
#                              mn_FC_SLF2_L,
#                              mn_FC_SLF3_L,
#                              mn_FC_SLF1_R,
#                              mn_FC_SLF2_R,
#                              mn_FC_SLF3_R)
# SLF_data_FC_long <- gather(SLF_data_FC, 
#                            "SLF_type",
#                            "FC_metric",
#                            mn_FC_SLF,
#                            mn_FC_SLF_L,
#                            mn_FC_SLF_R,
#                            mn_FC_SLF1_L,
#                            mn_FC_SLF2_L,
#                            mn_FC_SLF3_L,
#                            mn_FC_SLF1_R,
#                            mn_FC_SLF2_R,
#                            mn_FC_SLF3_R)
# #All tracts SLF FC (raincloud plot)
# ggplot(SLF_data_FC_long, aes(x = SLF_type, y = FC_metric, fill = Timepoint)) + 
#   geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#   geom_point(aes(y = FC_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
#   geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
#   stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
#   guides(fill = guide_legend(override.aes = list(shape = NA)))+
#   xlab("SLF tract") + 
#   ylab("Fibre Cross-section (FC)") +
#   scale_x_discrete(labels = c("mn_FC_SLF" = "Whole SLF", "mn_FC_SLF_L" = "Left SLF", "mn_FC_SLF_R" = "Right SLF", "mn_FC_SLF1_L" = "Left SLF1","mn_FC_SLF2_L" = "Left SLF2","mn_FC_SLF3_L" = "Left SLF3","mn_FC_SLF1_R" = "Right SLF1","mn_FC_SLF2_R" = "Right SLF2","mn_FC_SLF3_R" = "Right SLF3")) + 
#   scale_colour_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
#   theme_classic()


#FC just with left and right SLF 1, 2, 3
#put data into long format
SLF_data_FC <- dplyr::select(SLF_data, 
                             ParticipantID,
                             Group,
                             Timepoint,
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
ggplot(SLF_data_FC_long, aes(x = SLF_type, y = FC_metric, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = FC_metric, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF tract") + 
  ylab("Fibre Cross-section (FC)") +
  scale_x_discrete(labels = c("mn_FC_SLF1_L" = "Left SLF1","mn_FC_SLF2_L" = "Left SLF2","mn_FC_SLF3_L" = "Left SLF3","mn_FC_SLF1_R" = "Right SLF1","mn_FC_SLF2_R" = "Right SLF2","mn_FC_SLF3_R" = "Right SLF3")) + 
  scale_colour_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()

#Separate the groups - All tracts SLF FC (raincloud plot)
ggplot(SLF_data_FC_long, aes(x = SLF_type, y = FC_metric, fill = Timepoint)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fibre Cross-section (FC)") +
  scale_x_discrete(labels = c("mn_FC_SLF1_L" = "Left SLF1","mn_FC_SLF2_L" = "Left SLF2","mn_FC_SLF3_L" = "Left SLF3","mn_FC_SLF1_R" = "Right SLF1","mn_FC_SLF2_R" = "Right SLF2","mn_FC_SLF3_R" = "Right SLF3")) + 
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()



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


#Just display left SLF 1 and right SLF1 FDC (* sig) put data into long format
SLF_data_FDC <- dplyr::select(SLF_data, 
                              ParticipantID,
                              Group,
                              Timepoint,
                              mn_FDC_SLF1_L,
                              mn_FDC_SLF1_R)

SLF_data_FDC_long <- gather(SLF_data_FDC, 
                            "SLF_type",
                            "FDC_metric",
                            mn_FDC_SLF1_L,
                            mn_FDC_SLF1_R)
#Separate the groups - Just left SLF1 FDC and right SLF1 FDC (raincloud plot)
ggplot(SLF_data_FDC_long, aes(x = SLF_type, y = FDC_metric, fill = Timepoint)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fibre Density Cross-section (FDC)") +
  scale_x_discrete(labels = c("mn_FDC_SLF1_L" = "Left SLF1","mn_FDC_SLF1_R" = "Right SLF1")) + 
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()



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





