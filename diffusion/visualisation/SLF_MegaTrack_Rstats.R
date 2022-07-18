
#Script to analyse MegaTrack SLF data (from Flavio). 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 20/06/22



#------------------------------Setting up--------------------------------------#
#install packages/open libraries
pacman::p_load(dplyr, ggplot2, psych, car, multcomp, lsr, tidyr, BayesFactor, tidyverse, ppcor, nlme, effectsize, rstatix, sjstats, purrr, corrplot, compute.es, readxl)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#Read in neuropsych data file
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')

DPRC_neuropsych_data <- read.csv("cross-sectional_DPRC_neuropsych_data_lined_up_valid_participants.csv")
#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'

#convert variables
DPRC_neuropsych_data$ParticipantID <- as.factor(DPRC_neuropsych_data$ParticipantID)
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Sex_binary <- as.factor(DPRC_neuropsych_data$Sex_binary)
DPRC_neuropsych_data$Sex<- as.factor(DPRC_neuropsych_data$Sex)

#Read in the SLF data file: 
#setwd('/yourpathway/')
setwd('H:/ltah262/PhD/Diffusion/data/MegaTrack_SLF_results(Flavio)/')

MegaTrack_Left_SLF <- read_excel("Left_SLF123_measures.xlsx")
MegaTrack_Right_SLF <- read_excel("Right_SLF123_measures.xlsx")

#re-organise & re-name column names
colnames(MegaTrack_Left_SLF) <- c("count","ParticipantID","VoxelCount_LeftSLF1","VolumeML_LeftSLF1","HMOA1Mean_LeftSLF1","HMOA1std_LeftSLF1","AP6np1Mean_LeftSLF1","AP6np1std_LeftSLF1",
                                  "VoxelCount_LeftSLF2","VolumeML_LeftSLF2","HMOA1Mean_LeftSLF2","HMOA1std_LeftSLF2","AP6np1Mean_LeftSLF2","AP6np1std_LeftSLF2",
                                  "VoxelCount_LeftSLF3","VolumeML_LeftSLF3","HMOA1Mean_LeftSLF3","HMOA1std_LeftSLF3","AP6np1Mean_LeftSLF3","AP6np1std_LeftSLF3")
colnames(MegaTrack_Right_SLF) <- c("count","ParticipantID","VoxelCount_RightSLF1","VolumeML_RightSLF1","HMOA1Mean_RightSLF1","HMOA1std_RightSLF1","AP6np1Mean_RightSLF1","AP6np1std_RightSLF1",
                                  "VoxelCount_RightSLF2","VolumeML_RightSLF2","HMOA1Mean_RightSLF2","HMOA1std_RightSLF2","AP6np1Mean_RightSLF2","AP6np1std_RightSLF2",
                                  "VoxelCount_RightSLF3","VolumeML_RightSLF3","HMOA1Mean_RightSLF3","HMOA1std_RightSLF3","AP6np1Mean_RightSLF3","AP6np1std_RightSLF3")
#combine left and right SLF data files
SLF_MegaTrack_data<-left_join(MegaTrack_Left_SLF,MegaTrack_Right_SLF)

#add in Group classification column for each participant - data from another worksheet
SLF_MegaTrack_data$Group <- DPRC_neuropsych_data$Group
#add in covariates (clinical site) from the covariates_data dataframe to the SLF dataframe
SLF_MegaTrack_data$ClinSite_name <- DPRC_neuropsych_data$Clinical_site 
SLF_MegaTrack_data$Age <- DPRC_neuropsych_data$Age 
SLF_MegaTrack_data$Sex <- DPRC_neuropsych_data$Sex
SLF_MegaTrack_data$Sex_binary <- DPRC_neuropsych_data$Sex_binary
#add in neuropsych variables
neuropsych_test_names <- c("TrailsA.Raw","TrailsA.Z","TrailsB.Raw","TrailsB.Z","ColorNaming.Raw","ColorNaming.Z","WordReading.Raw","WordReading.Z","Inhibition.Raw","Inhibition.Z","LetFluency.Raw","LetFluency.Z","CatFluency.Raw","CatFluency.Z","Switching.Raw","Switching.z","HayBTime1.Raw","HayBTime1.z","HayBTime2.Raw","HayBTime2.z","HayBCatA.Raw","HayBCatA.z","HayBCatB.Raw","HayBCatB.z")
SLF_MegaTrack_data[neuropsych_test_names] <- DPRC_neuropsych_data[neuropsych_test_names]
#add in trend Group variable
Trend_Group <- as.numeric(DPRC_neuropsych_data$Group)

####--------------------------Descriptives----------------------------------####
#Track volume can be voxel count OR volume ml (equivalent to MRtrix's fibre cross-secton (FC)); 
#HMOA is equivalent to MRtrix's fibre density (FD):
#look at descriptive stats of the Left SLF1 between groups
SLF1_VoxelCount_L_descrip <- describeBy(SLF_MegaTrack_data$VoxelCount_LeftSLF1,SLF_MegaTrack_data$Group)
SLF1_VolumeML_L_descrip <- describeBy(SLF_MegaTrack_data$VolumeML_LeftSLF1, SLF_MegaTrack_data$Group)
SLF1_HMOA_L_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_LeftSLF1, SLF_MegaTrack_data$Group)
#SLF1_HMOA_L_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_LeftSLF1, SLF_MegaTrack_data$Group,mat=TRUE,digits=5)
mean(SLF_MegaTrack_data$VoxelCount_LeftSLF1)
sd(SLF_MegaTrack_data$VoxelCount_LeftSLF1)
mean(SLF_MegaTrack_data$VolumeML_LeftSLF1)
sd(SLF_MegaTrack_data$VolumeML_LeftSLF1)
mean(SLF_MegaTrack_data$HMOA1Mean_LeftSLF1)
sd(SLF_MegaTrack_data$HMOA1Mean_LeftSLF1)
#look at descriptive stats of the Left SLF2 between groups
SLF2_VoxelCount_L_descrip <- describeBy(SLF_MegaTrack_data$VoxelCount_LeftSLF2,SLF_MegaTrack_data$Group)
SLF2_VolumeML_L_descrip <- describeBy(SLF_MegaTrack_data$VolumeML_LeftSLF2, SLF_MegaTrack_data$Group)
SLF2_HMOA_L_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_LeftSLF2, SLF_MegaTrack_data$Group)
#SLF2_HMOA_L_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_LeftSLF2, SLF_MegaTrack_data$Group,mat=TRUE,digits=5)
mean(SLF_MegaTrack_data$VoxelCount_LeftSLF2)
sd(SLF_MegaTrack_data$VoxelCount_LeftSLF2)
mean(SLF_MegaTrack_data$VolumeML_LeftSLF2)
sd(SLF_MegaTrack_data$VolumeML_LeftSLF2)
mean(SLF_MegaTrack_data$HMOA1Mean_LeftSLF2)
sd(SLF_MegaTrack_data$HMOA1Mean_LeftSLF2)
#look at descriptive stats of the Left SLF3 between groups
SLF3_VoxelCount_L_descrip <- describeBy(SLF_MegaTrack_data$VoxelCount_LeftSLF3,SLF_MegaTrack_data$Group)
SLF3_VolumeML_L_descrip <- describeBy(SLF_MegaTrack_data$VolumeML_LeftSLF3, SLF_MegaTrack_data$Group)
SLF3_HMOA_L_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_LeftSLF3, SLF_MegaTrack_data$Group)
#SLF3_HMOA_L_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_LeftSLF3, SLF_MegaTrack_data$Group,mat=TRUE,digits=5)
mean(SLF_MegaTrack_data$VoxelCount_LeftSLF3)
sd(SLF_MegaTrack_data$VoxelCount_LeftSLF3)
mean(SLF_MegaTrack_data$VolumeML_LeftSLF3)
sd(SLF_MegaTrack_data$VolumeML_LeftSLF3)
mean(SLF_MegaTrack_data$HMOA1Mean_LeftSLF3)
sd(SLF_MegaTrack_data$HMOA1Mean_LeftSLF3)
#look at descriptive stats of the Right SLF1 between groups
SLF1_VoxelCount_R_descrip <- describeBy(SLF_MegaTrack_data$VoxelCount_RightSLF1,SLF_MegaTrack_data$Group)
SLF1_VolumeML_R_descrip <- describeBy(SLF_MegaTrack_data$VolumeML_RightSLF1, SLF_MegaTrack_data$Group)
SLF1_HMOA_R_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_RightSLF1, SLF_MegaTrack_data$Group)
#SLF1_HMOA_R_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_RightSLF1, SLF_MegaTrack_data$Group,mat=TRUE,digits=5)
mean(SLF_MegaTrack_data$VoxelCount_RightSLF1)
sd(SLF_MegaTrack_data$VoxelCount_RightSLF1)
mean(SLF_MegaTrack_data$VolumeML_RightSLF1)
sd(SLF_MegaTrack_data$VolumeML_RightSLF1)
mean(SLF_MegaTrack_data$HMOA1Mean_RightSLF1)
sd(SLF_MegaTrack_data$HMOA1Mean_RightSLF1)
#look at descriptive stats of the Right SLF2 between groups
SLF2_VoxelCount_R_descrip <- describeBy(SLF_MegaTrack_data$VoxelCount_RightSLF2,SLF_MegaTrack_data$Group)
SLF2_VolumeML_R_descrip <- describeBy(SLF_MegaTrack_data$VolumeML_RightSLF2, SLF_MegaTrack_data$Group)
SLF2_HMOA_R_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_RightSLF2, SLF_MegaTrack_data$Group)
#SLF2_HMOA_R_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_RightSLF2, SLF_MegaTrack_data$Group,mat=TRUE,digits=5)
mean(SLF_MegaTrack_data$VoxelCount_RightSLF2)
sd(SLF_MegaTrack_data$VoxelCount_RightSLF2)
mean(SLF_MegaTrack_data$VolumeML_RightSLF2)
sd(SLF_MegaTrack_data$VolumeML_RightSLF2)
mean(SLF_MegaTrack_data$HMOA1Mean_RightSLF2)
sd(SLF_MegaTrack_data$HMOA1Mean_RightSLF2)
#look at descriptive stats of the Right SLF3 between groups
SLF3_VoxelCount_R_descrip <- describeBy(SLF_MegaTrack_data$VoxelCount_RightSLF3,SLF_MegaTrack_data$Group)
SLF3_VolumeML_R_descrip <- describeBy(SLF_MegaTrack_data$VolumeML_RightSLF3, SLF_MegaTrack_data$Group)
SLF3_HMOA_R_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_RightSLF3, SLF_MegaTrack_data$Group)
#SLF3_HMOA_R_descrip <- describeBy(SLF_MegaTrack_data$HMOA1Mean_RightSLF3, SLF_MegaTrack_data$Group,mat=TRUE,digits=5)
mean(SLF_MegaTrack_data$VoxelCount_RightSLF3)
sd(SLF_MegaTrack_data$VoxelCount_RightSLF3)
mean(SLF_MegaTrack_data$VolumeML_RightSLF3)
sd(SLF_MegaTrack_data$VolumeML_RightSLF3)
mean(SLF_MegaTrack_data$HMOA1Mean_RightSLF3)
sd(SLF_MegaTrack_data$HMOA1Mean_RightSLF3)

###------------------ ANOVA & ANCOVA (age + sex) ---------------------------####

#run ANOVA to see if there are significant differences between groups for each track
#for Left SLF1
SLF1_VoxelCount_mod_L <- lm(VoxelCount_LeftSLF1 ~ Group, data = SLF_MegaTrack_data)
SLF1_VolumeML_mod_L <- lm(VolumeML_LeftSLF1 ~ Group, data = SLF_MegaTrack_data)
SLF1_HMOA1Mean_mod_L <- lm(HMOA1Mean_LeftSLF1 ~ Group, data = SLF_MegaTrack_data)
#for Left SLF2
SLF2_VoxelCount_mod_L <- lm(VoxelCount_LeftSLF2 ~ Group, data = SLF_MegaTrack_data)
SLF2_VolumeML_mod_L <- lm(VolumeML_LeftSLF2 ~ Group, data = SLF_MegaTrack_data)
SLF2_HMOA1Mean_mod_L <- lm(HMOA1Mean_LeftSLF2 ~ Group, data = SLF_MegaTrack_data)
#for Left SLF3
SLF3_VoxelCount_mod_L <- lm(VoxelCount_LeftSLF3 ~ Group, data = SLF_MegaTrack_data)
SLF3_VolumeML_mod_L <- lm(VolumeML_LeftSLF3 ~ Group, data = SLF_MegaTrack_data)
SLF3_HMOA1Mean_mod_L <- lm(HMOA1Mean_LeftSLF3 ~ Group, data = SLF_MegaTrack_data)
#for Right SLF1
SLF1_VoxelCount_mod_R <- lm(VoxelCount_RightSLF1 ~ Group, data = SLF_MegaTrack_data)
SLF1_VolumeML_mod_R <- lm(VolumeML_RightSLF1 ~ Group, data = SLF_MegaTrack_data)
SLF1_HMOA1Mean_mod_R <- lm(HMOA1Mean_RightSLF1 ~ Group, data = SLF_MegaTrack_data)
#for Right SLF2
SLF2_VoxelCount_mod_R <- lm(VoxelCount_RightSLF2 ~ Group, data = SLF_MegaTrack_data)
SLF2_VolumeML_mod_R <- lm(VolumeML_RightSLF2 ~ Group, data = SLF_MegaTrack_data)
SLF2_HMOA1Mean_mod_R <- lm(HMOA1Mean_RightSLF2 ~ Group, data = SLF_MegaTrack_data)
#for Right SLF3
SLF3_VoxelCount_mod_R <- lm(VoxelCount_RightSLF3 ~ Group, data = SLF_MegaTrack_data)
SLF3_VolumeML_mod_R <- lm(VolumeML_RightSLF3 ~ Group, data = SLF_MegaTrack_data)
SLF3_HMOA1Mean_mod_R <- lm(HMOA1Mean_RightSLF3 ~ Group, data = SLF_MegaTrack_data)

#run ANOVAs on the model
#for Left SLF1
anova(SLF1_VoxelCount_mod_L)
anova(SLF1_VolumeML_mod_L)
anova(SLF1_HMOA1Mean_mod_L)
#for Left SLF2
anova(SLF2_VoxelCount_mod_L)
anova(SLF2_VolumeML_mod_L)
anova(SLF2_HMOA1Mean_mod_L)
#for Left SLF3
anova(SLF3_VoxelCount_mod_L)
anova(SLF3_VolumeML_mod_L)
anova(SLF3_HMOA1Mean_mod_L)
#for Right SLF1
anova(SLF1_VoxelCount_mod_R)
anova(SLF1_VolumeML_mod_R)
anova(SLF1_HMOA1Mean_mod_R)
#for Right SLF2
anova(SLF2_VoxelCount_mod_R)
anova(SLF2_VolumeML_mod_R)
anova(SLF2_HMOA1Mean_mod_R)
#for Right SLF3
anova(SLF3_VoxelCount_mod_R)
anova(SLF3_VolumeML_mod_R)
anova(SLF3_HMOA1Mean_mod_R)

#Run ANOVAs with covariates (age & sex)
#for Left SLF1
SLF1_VoxelCount_2covar_mod_L <- lm(VoxelCount_LeftSLF1 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF1_VolumeML_2covar_mod_L <- lm(VolumeML_LeftSLF1 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF1_HMOA1Mean_2covar_mod_L <- lm(HMOA1Mean_LeftSLF1 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
#for Left SLF2
SLF2_VoxelCount_2covar_mod_L <- lm(VoxelCount_LeftSLF2 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF2_VolumeML_2covar_mod_L <- lm(VolumeML_LeftSLF2 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF2_HMOA1Mean_2covar_mod_L <- lm(HMOA1Mean_LeftSLF2 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
#for Left SLF3
SLF3_VoxelCount_2covar_mod_L <- lm(VoxelCount_LeftSLF3 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF3_VolumeML_2covar_mod_L <- lm(VolumeML_LeftSLF3 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF3_HMOA1Mean_2covar_mod_L <- lm(HMOA1Mean_LeftSLF3 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
#for Right SLF1
SLF1_VoxelCount_2covar_mod_R <- lm(VoxelCount_RightSLF1 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF1_VolumeML_2covar_mod_R <- lm(VolumeML_RightSLF1 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF1_HMOA1Mean_2covar_mod_R <- lm(HMOA1Mean_RightSLF1 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
#for Right SLF2
SLF2_VoxelCount_2covar_mod_R <- lm(VoxelCount_RightSLF2 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF2_VolumeML_2covar_mod_R <- lm(VolumeML_RightSLF2 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF2_HMOA1Mean_2covar_mod_R <- lm(HMOA1Mean_RightSLF2 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
#for Right SLF3
SLF3_VoxelCount_2covar_mod_R <- lm(VoxelCount_RightSLF3 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF3_VolumeML_2covar_mod_R <- lm(VolumeML_RightSLF3 ~ Group + Age + Sex, data = SLF_MegaTrack_data)
SLF3_HMOA1Mean_2covar_mod_R <- lm(HMOA1Mean_RightSLF3 ~ Group + Age + Sex, data = SLF_MegaTrack_data)

#run ANCOVAs (w/ age & sex as covariates)
#for Left SLF1
anova(SLF1_VoxelCount_2covar_mod_L)
anova(SLF1_VolumeML_2covar_mod_L)
anova(SLF1_HMOA1Mean_2covar_mod_L)
#for Left SLF2
anova(SLF2_VoxelCount_2covar_mod_L)
anova(SLF2_VolumeML_2covar_mod_L)
anova(SLF2_HMOA1Mean_2covar_mod_L)
#for Left SLF3
anova(SLF3_VoxelCount_2covar_mod_L)
anova(SLF3_VolumeML_2covar_mod_L)
anova(SLF3_HMOA1Mean_2covar_mod_L)
#for Right SLF1
anova(SLF1_VoxelCount_2covar_mod_R)
anova(SLF1_VolumeML_2covar_mod_R)
anova(SLF1_HMOA1Mean_2covar_mod_R)
#for Right SLF2
anova(SLF2_VoxelCount_2covar_mod_R)
anova(SLF2_VolumeML_2covar_mod_R)
anova(SLF2_HMOA1Mean_2covar_mod_R)
#for Right SLF3
anova(SLF3_VoxelCount_2covar_mod_R)
anova(SLF3_VolumeML_2covar_mod_R)
anova(SLF3_HMOA1Mean_2covar_mod_R)


####------------------------- Linear Trend ---------------------------------####
#for Left SLF1
LSLF1_VoxelCount_LinTrend_mod <- lm(VoxelCount_LeftSLF1 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
LSLF1_VolumeML_LinTrend_mod <- lm(VolumeML_LeftSLF1 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
LSLF1_HMOA1Mean_LinTrend_mod <- lm(VoxelCount_LeftSLF1 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
#for Left SLF2
LSLF2_VoxelCount_LinTrend_mod <- lm(VoxelCount_LeftSLF2 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
LSLF2_VolumeML_LinTrend_mod <- lm(VolumeML_LeftSLF2 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
LSLF2_HMOA1Mean_LinTrend_mod <- lm(VoxelCount_LeftSLF2 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
#for Left SLF3
LSLF3_VoxelCount_LinTrend_mod <- lm(VoxelCount_LeftSLF3 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
LSLF3_VolumeML_LinTrend_mod <- lm(VolumeML_LeftSLF3 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
LSLF3_HMOA1Mean_LinTrend_mod <- lm(VoxelCount_LeftSLF3 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
#for Right SLF1
RSLF1_VoxelCount_LinTrend_mod <- lm(VoxelCount_RightSLF1 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
RSLF1_VolumeML_LinTrend_mod <- lm(VolumeML_RightSLF1 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
RSLF1_HMOA1Mean_LinTrend_mod <- lm(VoxelCount_RightSLF1 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
#for Right SLF2
RSLF2_VoxelCount_LinTrend_mod <- lm(VoxelCount_RightSLF2 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
RSLF2_VolumeML_LinTrend_mod <- lm(VolumeML_RightSLF2 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
RSLF2_HMOA1Mean_LinTrend_mod <- lm(VoxelCount_RightSLF2 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
#for Right SLF3
RSLF3_VoxelCount_LinTrend_mod <- lm(VoxelCount_RightSLF3 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
RSLF3_VolumeML_LinTrend_mod <- lm(VolumeML_RightSLF3 ~ Trend_Group + Group, data = SLF_MegaTrack_data)
RSLF3_HMOA1Mean_LinTrend_mod <- lm(VoxelCount_RightSLF3 ~ Trend_Group + Group, data = SLF_MegaTrack_data)

#Run Linear Trend ANOVAs
#for Left SLF1
anova(LSLF1_VoxelCount_LinTrend_mod)
anova(LSLF1_VolumeML_LinTrend_mod)
anova(LSLF1_HMOA1Mean_LinTrend_mod)
#for Left SLF2
anova(LSLF2_VoxelCount_LinTrend_mod)
anova(LSLF2_VolumeML_LinTrend_mod)
anova(LSLF2_HMOA1Mean_LinTrend_mod)
#for Left SLF3
anova(LSLF3_VoxelCount_LinTrend_mod)
anova(LSLF3_VolumeML_LinTrend_mod)
anova(LSLF3_HMOA1Mean_LinTrend_mod)
#for Right SLF1
anova(RSLF1_VoxelCount_LinTrend_mod)
anova(RSLF1_VolumeML_LinTrend_mod)
anova(RSLF1_HMOA1Mean_LinTrend_mod)
#for Right SLF2
anova(RSLF2_VoxelCount_LinTrend_mod)
anova(RSLF2_VolumeML_LinTrend_mod)
anova(RSLF2_HMOA1Mean_LinTrend_mod)
#for Right SLF3
anova(RSLF3_VoxelCount_LinTrend_mod)
anova(RSLF3_VolumeML_LinTrend_mod)
anova(RSLF3_HMOA1Mean_LinTrend_mod)

#Run Linear Trend with covariates (age & sex)
#for Left SLF1
LSLF1_VoxelCount_2covar_LinTrend_mod <- lm(VoxelCount_LeftSLF1 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
LSLF1_VolumeML_2covar_LinTrend_mod <- lm(VolumeML_LeftSLF1 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
LSLF1_HMOA1Mean_2covar_LinTrend_mod <- lm(HMOA1Mean_LeftSLF1 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
#for Left SLF2
LSLF2_VoxelCount_2covar_LinTrend_mod <- lm(VoxelCount_LeftSLF2 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
LSLF2_VolumeML_2covar_LinTrend_mod <- lm(VolumeML_LeftSLF2 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
LSLF2_HMOA1Mean_2covar_LinTrend_mod <- lm(HMOA1Mean_LeftSLF2 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
#for Left SLF3
LSLF3_VoxelCount_2covar_LinTrend_mod <- lm(VoxelCount_LeftSLF3 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
LSLF3_VolumeML_2covar_LinTrend_mod <- lm(VolumeML_LeftSLF3 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
LSLF3_HMOA1Mean_2covar_LinTrend_mod <- lm(HMOA1Mean_LeftSLF3 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
#for Right SLF1
RSLF1_VoxelCount_2covar_LinTrend_mod <- lm(VoxelCount_RightSLF1 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
RSLF1_VolumeML_2covar_LinTrend_mod <- lm(VolumeML_RightSLF1 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
RSLF1_HMOA1Mean_2covar_LinTrend_mod <- lm(HMOA1Mean_RightSLF1 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
#for Right SLF2
RSLF2_VoxelCount_2covar_LinTrend_mod <- lm(VoxelCount_RightSLF2 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
RSLF2_VolumeML_2covar_LinTrend_mod <- lm(VolumeML_RightSLF2 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
RSLF2_HMOA1Mean_2covar_LinTrend_mod <- lm(HMOA1Mean_RightSLF2 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
#for Right SLF3
RSLF3_VoxelCount_2covar_LinTrend_mod <- lm(VoxelCount_RightSLF3 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
RSLF3_VolumeML_2covar_LinTrend_mod <- lm(VolumeML_RightSLF3 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)
RSLF3_HMOA1Mean_2covar_LinTrend_mod <- lm(HMOA1Mean_RightSLF3 ~ Trend_Group + Group + Age + Sex, data = SLF_MegaTrack_data)

#Linear Trend ANOVAs with covariates (sex & age)
#for Left SLF1
anova(LSLF1_VoxelCount_2covar_LinTrend_mod)
anova(LSLF1_VolumeML_2covar_LinTrend_mod)
anova(LSLF1_HMOA1Mean_2covar_LinTrend_mod)
#for Left SLF2
anova(LSLF2_VoxelCount_2covar_LinTrend_mod)
anova(LSLF2_VolumeML_2covar_LinTrend_mod)
anova(LSLF2_HMOA1Mean_2covar_LinTrend_mod)
#for Left SLF3
anova(LSLF3_VoxelCount_2covar_LinTrend_mod)
anova(LSLF3_VolumeML_2covar_LinTrend_mod)
anova(LSLF3_HMOA1Mean_2covar_LinTrend_mod) #sig
#for Right SLF1
anova(RSLF1_VoxelCount_2covar_LinTrend_mod)
anova(RSLF1_VolumeML_2covar_LinTrend_mod)
anova(RSLF1_HMOA1Mean_2covar_LinTrend_mod)
#for Right SLF2
anova(RSLF2_VoxelCount_2covar_LinTrend_mod)
anova(RSLF2_VolumeML_2covar_LinTrend_mod)
anova(RSLF2_HMOA1Mean_2covar_LinTrend_mod)
#for Right SLF3
anova(RSLF3_VoxelCount_2covar_LinTrend_mod)
anova(RSLF3_VolumeML_2covar_LinTrend_mod)
anova(RSLF3_HMOA1Mean_2covar_LinTrend_mod)

####------------------------- plot data ------------------------------------####

#Plot data with left 1,2,3 and right 1,2,3 SLF tracks
#for voxel count
SLF_data_VoxelCount <- dplyr::select(SLF_MegaTrack_data, 
                             ParticipantID,
                             Group,
                             VoxelCount_LeftSLF1,
                             VoxelCount_LeftSLF2,
                             VoxelCount_LeftSLF3,
                             VoxelCount_RightSLF1,
                             VoxelCount_RightSLF2,
                             VoxelCount_RightSLF3)
SLF_data_VoxelCount_long <- gather(SLF_data_VoxelCount, 
                           "SLF_type",
                           "VoxelCount",
                           VoxelCount_LeftSLF1,
                           VoxelCount_LeftSLF2,
                           VoxelCount_LeftSLF3,
                           VoxelCount_RightSLF1,
                           VoxelCount_RightSLF2,
                           VoxelCount_RightSLF3)
#All tracts SLF FD (raincloud plot)
ggplot(SLF_data_VoxelCount_long, aes(x = SLF_type, y = VoxelCount, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = VoxelCount, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Voxel Count") +
  scale_x_discrete(labels = c("VoxelCount_LeftSLF1" = "Left SLF1","VoxelCount_LeftSLF2" = "Left SLF2","VoxelCount_LeftSLF3" = "Left SLF3","VoxelCount_RightSLF1" = "Right SLF1","VoxelCount_RightSLF2" = "Right SLF2","VoxelCount_RightSLF3" = "Right SLF3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()

#for volume ml
SLF_data_VolumeML <- dplyr::select(SLF_MegaTrack_data, 
                                     ParticipantID,
                                     Group,
                                     VolumeML_LeftSLF1,
                                     VolumeML_LeftSLF2,
                                     VolumeML_LeftSLF3,
                                     VolumeML_RightSLF1,
                                     VolumeML_RightSLF2,
                                     VolumeML_RightSLF3)
SLF_data_VolumeML_long <- gather(SLF_data_VolumeML, 
                                   "SLF_type",
                                   "VolumeML",
                                   VolumeML_LeftSLF1,
                                   VolumeML_LeftSLF2,
                                   VolumeML_LeftSLF3,
                                   VolumeML_RightSLF1,
                                   VolumeML_RightSLF2,
                                   VolumeML_RightSLF3)
#if you want to put the tracts in order (from left right 1, left right 2, left right 3): 
SLF_data_VolumeML_long$SLF_type <- factor(SLF_data_VolumeML_long$SLF_type, levels=c("VolumeML_LeftSLF1","VolumeML_RightSLF1","VolumeML_LeftSLF2","VolumeML_RightSLF2","VolumeML_LeftSLF3","VolumeML_RightSLF3"))

#All tracts SLF Volume ML (raincloud plot)
ggplot(SLF_data_VolumeML_long, aes(x = SLF_type, y = VolumeML, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = VolumeML, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Volume ML") +
  scale_x_discrete(labels = c("VolumeML_LeftSLF1" = "Left SLF1","VolumeML_LeftSLF2" = "Left SLF2","VolumeML_LeftSLF3" = "Left SLF3","VolumeML_RightSLF1" = "Right SLF1","VolumeML_RightSLF2" = "Right SLF2","VolumeML_RightSLF3" = "Right SLF3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()

#for HMOA
SLF_data_HMOA1Mean <- dplyr::select(SLF_MegaTrack_data, 
                                     ParticipantID,
                                     Group,
                                     HMOA1Mean_LeftSLF1,
                                     HMOA1Mean_LeftSLF2,
                                     HMOA1Mean_LeftSLF3,
                                     HMOA1Mean_RightSLF1,
                                     HMOA1Mean_RightSLF2,
                                     HMOA1Mean_RightSLF3)
SLF_data_HMOA1Mean_long <- gather(SLF_data_HMOA1Mean, 
                                   "SLF_type",
                                   "HMOA1Mean",
                                   HMOA1Mean_LeftSLF1,
                                   HMOA1Mean_LeftSLF2,
                                   HMOA1Mean_LeftSLF3,
                                   HMOA1Mean_RightSLF1,
                                   HMOA1Mean_RightSLF2,
                                   HMOA1Mean_RightSLF3)
#if you want to put the tracts in order (from left right 1, left right 2, left right 3): 
SLF_data_HMOA1Mean_long$SLF_type <- factor(SLF_data_HMOA1Mean_long$SLF_type, levels=c("HMOA1Mean_LeftSLF1","HMOA1Mean_RightSLF1","HMOA1Mean_LeftSLF2","HMOA1Mean_RightSLF2","HMOA1Mean_LeftSLF3","HMOA1Mean_RightSLF3"))

#All tracts SLF HMOA (raincloud plot)
ggplot(SLF_data_HMOA1Mean_long, aes(x = SLF_type, y = HMOA1Mean, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HMOA1Mean, colour = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Hindrance Modulated Orientation Anisotropy (HMOA)") +
  scale_x_discrete(labels = c("HMOA1Mean_LeftSLF1" = "Left SLF1","HMOA1Mean_LeftSLF2" = "Left SLF2","HMOA1Mean_LeftSLF3" = "Left SLF3","HMOA1Mean_RightSLF1" = "Right SLF1","HMOA1Mean_RightSLF2" = "Right SLF2","HMOA1Mean_RightSLF3" = "Right SLF3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()



####----- Correlation - do MegaTrack and MRtrix3's data SLF metrics correlate w/ one another? -----####


#navigate to the correct pathway which contains the SLF metric text files: 
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/template/TOI')

#read in text file of data
SLF_MRtrix_data <- cbind.data.frame(read.table("FD_SLF_whole_TOI.txt", header = T), 
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
colnames(SLF_MRtrix_data) <- c("mn_FD_SLF", "md_FD_SLF", "std_FD_SLF", "std_rv_FD_SLF", "min_FD_SLF", "max_FD_SLF", "count_FD_SLF", 
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

#check correlations between Megatrack and MRtrix3. We will check these correrlations: 
  #1. MRtrix3's fibre cross-section (FC) vs. MegaTrack's Track Volume (Volume ML) 
  #2. MRtrix3's fibre density (FD) vs. MegaTrack's Hindrance Modulated Orientation Anisotropy (HMOA)

#combine SLF datasets
SLF_combined_data<-cbind(SLF_MRtrix_data,SLF_MegaTrack_data)

#Run correlation test and plot data
#for left SLF1
  cor.test(SLF_combined_data$mn_FC_SLF1_L, SLF_combined_data$VolumeML_LeftSLF1)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FC_SLF1_L, y = VolumeML_LeftSLF1)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Cross-section (FC)", y = "Track Volume (volume ml)")+
    theme_classic()
  cor.test(SLF_combined_data$mn_FD_SLF1_L, SLF_combined_data$HMOA1Mean_LeftSLF1)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FD_SLF1_L, y = HMOA1Mean_LeftSLF1)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Density (FD)", y = "Hindrance Modulated Orientation Anisotropy (HMOA)")+
    theme_classic()

#for left SLF2
  cor.test(SLF_combined_data$mn_FC_SLF2_L, SLF_combined_data$VolumeML_LeftSLF2)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FC_SLF2_L, y = VolumeML_LeftSLF2)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Cross-section (FC)", y = "Track Volume (volume ml)")+
    theme_classic()
  cor.test(SLF_combined_data$mn_FD_SLF2_L, SLF_combined_data$HMOA1Mean_LeftSLF2)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FD_SLF2_L, y = HMOA1Mean_LeftSLF2)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Density (FD)", y = "Hindrance Modulated Orientation Anisotropy (HMOA)")+
    theme_classic()
  
#for left SLF3
  cor.test(SLF_combined_data$mn_FC_SLF3_L, SLF_combined_data$VolumeML_LeftSLF3)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FC_SLF3_L, y = VolumeML_LeftSLF3)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Cross-section (FC)", y = "Track Volume (volume ml)")+
    theme_classic()
  cor.test(SLF_combined_data$mn_FD_SLF3_L, SLF_combined_data$HMOA1Mean_LeftSLF3)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FD_SLF3_L, y = HMOA1Mean_LeftSLF3)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Density (FD)", y = "Hindrance Modulated Orientation Anisotropy (HMOA)")+
    theme_classic()
  
#for Right SLF1
  cor.test(SLF_combined_data$mn_FC_SLF1_R, SLF_combined_data$VolumeML_RightSLF1)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FC_SLF1_R, y = VolumeML_RightSLF1)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Cross-section (FC)", y = "Track Volume (volume ml)")+
    theme_classic()
  cor.test(SLF_combined_data$mn_FD_SLF1_R, SLF_combined_data$HMOA1Mean_RightSLF1)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FD_SLF1_R, y = HMOA1Mean_RightSLF1)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Density (FD)", y = "Hindrance Modulated Orientation Anisotropy (HMOA)")+
    theme_classic()
  
#for Right SLF2
  cor.test(SLF_combined_data$mn_FC_SLF2_R, SLF_combined_data$VolumeML_RightSLF2)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FC_SLF2_R, y = VolumeML_RightSLF2)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Cross-section (FC)", y = "Track Volume (volume ml)")+
    theme_classic()
  cor.test(SLF_combined_data$mn_FD_SLF2_R, SLF_combined_data$HMOA1Mean_RightSLF2)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FD_SLF2_R, y = HMOA1Mean_RightSLF2)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Density (FD)", y = "Hindrance Modulated Orientation Anisotropy (HMOA)")+
    theme_classic()
  
#for Right SLF3
  cor.test(SLF_combined_data$mn_FC_SLF3_R, SLF_combined_data$VolumeML_RightSLF3)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FC_SLF3_R, y = VolumeML_RightSLF3)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Cross-section (FC)", y = "Track Volume (volume ml)")+
    theme_classic()
  cor.test(SLF_combined_data$mn_FD_SLF3_R, SLF_combined_data$HMOA1Mean_RightSLF3)
  #scatterplot
  ggplot(SLF_combined_data, aes(x = mn_FD_SLF3_R, y = HMOA1Mean_RightSLF3)) +
    geom_point()+
    geom_smooth(method = "lm", colour = "purple", se = FALSE)+
    labs(x = "Fibre Density (FD)", y = "Hindrance Modulated Orientation Anisotropy (HMOA)")+
    theme_classic()







