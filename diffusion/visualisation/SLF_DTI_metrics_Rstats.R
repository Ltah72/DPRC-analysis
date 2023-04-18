#Perform statistics on the diffusion tensor imaging (DTI) metrics. This includes fractional anisotropy (MD), mean diffusivity 
#(MD - can also be referred to apparent diffusion coefficient (ADC)), axial diffusivity (AD), and radial diffusivity (RD).

#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 24/07/22


#------------------------------Setting up--------------------------------------#
#install packages/open libraries
pacman::p_load(dplyr, ggplot2, psych, car, purrr, vroom, multcomp, rstatix, tidyr, lsr)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#first read in the demographics group data file: 
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


#navigate to the correct pathway which contains the SLF DTI metric text files: 
#setwd('/yourpathway/')
setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/DTI_metrics/')

#for FA:
setwd('fractional_anisotropy_FA/')
#read in all of the left SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_LSLF1_values
#read in all of the left SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_LSLF2_values
#read in all of the left SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_LSLF3_values
#read in all of the right SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_RSLF1_values
#read in all of the right SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_RSLF2_values
#read in all of the right SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_RSLF3_values

#for MD:
setwd('..')#go up one directory
setwd('mean_diffusivity_MD/')
#read in all of the left SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_LSLF1_values
#read in all of the left SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_LSLF2_values
#read in all of the left SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_LSLF3_values
#read in all of the right SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_RSLF1_values
#read in all of the right SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_RSLF2_values
#read in all of the right SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_RSLF3_values

#for AD:
setwd('..')#go up one directory
setwd('axial_diffusivity_AD/')
#read in all of the left SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_LSLF1_values
#read in all of the left SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_LSLF2_values
#read in all of the left SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_LSLF3_values
#read in all of the right SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_RSLF1_values
#read in all of the right SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_RSLF2_values
#read in all of the right SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_RSLF3_values

#for RD:
setwd('..')#go up one directory
setwd('radial_diffusivity_RD/')
#read in all of the left SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_LSLF1_values
#read in all of the left SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_LSLF2_values
#read in all of the left SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_LSLF3_values
#read in all of the right SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_RSLF1_values
#read in all of the right SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_RSLF2_values
#read in all of the right SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_RSLF3_values


#put all mean DTI metric values into a dataframe
SLF_DTI_metrics_data <- cbind.data.frame(mean_FA_LSLF1_values,mean_FA_LSLF2_values,mean_FA_LSLF3_values,
                                         mean_FA_RSLF1_values,mean_FA_RSLF2_values,mean_FA_RSLF3_values,
                                         mean_MD_LSLF1_values,mean_MD_LSLF2_values,mean_MD_LSLF3_values,
                                         mean_MD_RSLF1_values,mean_MD_RSLF2_values,mean_MD_RSLF3_values,
                                         mean_AD_LSLF1_values,mean_AD_LSLF2_values,mean_AD_LSLF3_values,
                                         mean_AD_RSLF1_values,mean_AD_RSLF2_values,mean_AD_RSLF3_values,
                                         mean_RD_LSLF1_values,mean_RD_LSLF2_values,mean_RD_LSLF3_values,
                                         mean_RD_RSLF1_values,mean_RD_RSLF2_values,mean_RD_RSLF3_values)

#add in Participant ID  - data from another worksheet
SLF_DTI_metrics_data$ParticipantID <- DPRC_neuropsych_data$ParticipantID
#add in Group classification column for each participant - data from another worksheet
SLF_DTI_metrics_data$Group <- DPRC_neuropsych_data$Group
#add in covariates (clinical site) from the covariates_data dataframe to the SLF dataframe
SLF_DTI_metrics_data$ClinSite_name <- DPRC_neuropsych_data$Clinical_site 
SLF_DTI_metrics_data$Age <- DPRC_neuropsych_data$Age 
SLF_DTI_metrics_data$Sex <- DPRC_neuropsych_data$Sex
SLF_DTI_metrics_data$Sex_binary <- DPRC_neuropsych_data$Sex_binary
#add in trend Group variable
Trend_Group <- as.numeric(SLF_DTI_metrics_data$Group)



####--------------------------Descriptives----------------------------------####
#for FA:
SLF1_FA_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_LSLF1_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_FA_LSLF1_values)
sd(SLF_DTI_metrics_data$mean_FA_LSLF1_values)
SLF2_FA_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_LSLF2_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_FA_LSLF2_values)
sd(SLF_DTI_metrics_data$mean_FA_LSLF2_values)
SLF3_FA_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_LSLF3_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_FA_LSLF3_values)
sd(SLF_DTI_metrics_data$mean_FA_LSLF3_values)
SLF1_FA_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_RSLF1_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_FA_RSLF1_values)
sd(SLF_DTI_metrics_data$mean_FA_RSLF1_values)
SLF2_FA_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_RSLF2_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_FA_RSLF2_values)
sd(SLF_DTI_metrics_data$mean_FA_RSLF2_values)
SLF3_FA_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_RSLF3_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_FA_RSLF3_values)
sd(SLF_DTI_metrics_data$mean_FA_RSLF3_values)

#for MD:
SLF1_MD_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_LSLF1_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_MD_LSLF1_values)
sd(SLF_DTI_metrics_data$mean_MD_LSLF1_values)
SLF2_MD_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_LSLF2_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_MD_LSLF2_values)
sd(SLF_DTI_metrics_data$mean_MD_LSLF2_values)
SLF3_MD_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_LSLF3_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_MD_LSLF3_values)
sd(SLF_DTI_metrics_data$mean_MD_LSLF3_values)
SLF1_MD_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_RSLF1_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_MD_RSLF1_values)
sd(SLF_DTI_metrics_data$mean_MD_RSLF1_values)
SLF2_MD_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_RSLF2_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_MD_RSLF2_values)
sd(SLF_DTI_metrics_data$mean_MD_RSLF2_values)
SLF3_MD_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_RSLF3_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_MD_RSLF3_values)
sd(SLF_DTI_metrics_data$mean_MD_RSLF3_values)

#for AD:
SLF1_AD_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_LSLF1_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_AD_LSLF1_values)
sd(SLF_DTI_metrics_data$mean_AD_LSLF1_values)
SLF2_AD_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_LSLF2_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_AD_LSLF2_values)
sd(SLF_DTI_metrics_data$mean_AD_LSLF2_values)
SLF3_AD_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_LSLF3_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_AD_LSLF3_values)
sd(SLF_DTI_metrics_data$mean_AD_LSLF3_values)
SLF1_AD_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_RSLF1_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_AD_RSLF1_values)
sd(SLF_DTI_metrics_data$mean_AD_RSLF1_values)
SLF2_AD_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_RSLF2_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_AD_RSLF2_values)
sd(SLF_DTI_metrics_data$mean_AD_RSLF2_values)
SLF3_AD_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_RSLF3_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_AD_RSLF3_values)
sd(SLF_DTI_metrics_data$mean_AD_RSLF3_values)

#for RD:
SLF1_RD_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_LSLF1_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_RD_LSLF1_values)
sd(SLF_DTI_metrics_data$mean_RD_LSLF1_values)
SLF2_RD_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_LSLF2_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_RD_LSLF2_values)
sd(SLF_DTI_metrics_data$mean_RD_LSLF2_values)
SLF3_RD_L_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_LSLF3_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_RD_LSLF3_values)
sd(SLF_DTI_metrics_data$mean_RD_LSLF3_values)
SLF1_RD_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_RSLF1_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_RD_RSLF1_values)
sd(SLF_DTI_metrics_data$mean_RD_RSLF1_values)
SLF2_RD_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_RSLF2_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_RD_RSLF2_values)
sd(SLF_DTI_metrics_data$mean_RD_RSLF2_values)
SLF3_RD_R_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_RSLF3_values, SLF_DTI_metrics_data$Group)
#find mean & SD from total sample:
mean(SLF_DTI_metrics_data$mean_RD_RSLF3_values)
sd(SLF_DTI_metrics_data$mean_RD_RSLF3_values)

###------------------ ANOVA & ANCOVA (age) + (age + sex) ---------------------------####
#for FA:
SLF1_FA_mod_L <- lm(mean_FA_LSLF1_values ~ Group, data = SLF_DTI_metrics_data)
SLF2_FA_mod_L <- lm(mean_FA_LSLF2_values ~ Group, data = SLF_DTI_metrics_data)
SLF3_FA_mod_L <- lm(mean_FA_LSLF3_values ~ Group, data = SLF_DTI_metrics_data)
SLF1_FA_mod_R <- lm(mean_FA_RSLF1_values ~ Group, data = SLF_DTI_metrics_data)
SLF2_FA_mod_R <- lm(mean_FA_RSLF2_values ~ Group, data = SLF_DTI_metrics_data)
SLF3_FA_mod_R <- lm(mean_FA_RSLF3_values ~ Group, data = SLF_DTI_metrics_data)
anova(SLF1_FA_mod_L)
anova(SLF2_FA_mod_L)
anova(SLF3_FA_mod_L)
anova(SLF1_FA_mod_R)
anova(SLF2_FA_mod_R)
anova(SLF3_FA_mod_R)

#ANCOVA (w/ age)
SLF1_FA_covar_mod_L <- lm(mean_FA_LSLF1_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF2_FA_covar_mod_L <- lm(mean_FA_LSLF2_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF3_FA_covar_mod_L <- lm(mean_FA_LSLF3_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF1_FA_covar_mod_R <- lm(mean_FA_RSLF1_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF2_FA_covar_mod_R <- lm(mean_FA_RSLF2_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF3_FA_covar_mod_R <- lm(mean_FA_RSLF3_values ~ Group+Age, data = SLF_DTI_metrics_data)
anova(SLF1_FA_covar_mod_L)
anova(SLF2_FA_covar_mod_L)
anova(SLF3_FA_covar_mod_L)
anova(SLF1_FA_covar_mod_R)
anova(SLF2_FA_covar_mod_R)
anova(SLF3_FA_covar_mod_R)

#for MD:
SLF1_MD_mod_L <- lm(mean_MD_LSLF1_values ~ Group, data = SLF_DTI_metrics_data)
SLF2_MD_mod_L <- lm(mean_MD_LSLF2_values ~ Group, data = SLF_DTI_metrics_data)
SLF3_MD_mod_L <- lm(mean_MD_LSLF3_values ~ Group, data = SLF_DTI_metrics_data)
SLF1_MD_mod_R <- lm(mean_MD_RSLF1_values ~ Group, data = SLF_DTI_metrics_data)
SLF2_MD_mod_R <- lm(mean_MD_RSLF2_values ~ Group, data = SLF_DTI_metrics_data)
SLF3_MD_mod_R <- lm(mean_MD_RSLF3_values ~ Group, data = SLF_DTI_metrics_data)
anova(SLF1_MD_mod_L)
anova(SLF2_MD_mod_L)
anova(SLF3_MD_mod_L)
anova(SLF1_MD_mod_R)
anova(SLF2_MD_mod_R)
anova(SLF3_MD_mod_R)
#ANCOVA (w/ age)
SLF1_MD_covar_mod_L <- lm(mean_MD_LSLF1_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF2_MD_covar_mod_L <- lm(mean_MD_LSLF2_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF3_MD_covar_mod_L <- lm(mean_MD_LSLF3_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF1_MD_covar_mod_R <- lm(mean_MD_RSLF1_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF2_MD_covar_mod_R <- lm(mean_MD_RSLF2_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF3_MD_covar_mod_R <- lm(mean_MD_RSLF3_values ~ Group+Age, data = SLF_DTI_metrics_data)
anova(SLF1_MD_covar_mod_L)
anova(SLF2_MD_covar_mod_L)
anova(SLF3_MD_covar_mod_L)
anova(SLF1_MD_covar_mod_R)
anova(SLF2_MD_covar_mod_R)
anova(SLF3_MD_covar_mod_R)

#for AD:
SLF1_AD_mod_L <- lm(mean_AD_LSLF1_values ~ Group, data = SLF_DTI_metrics_data)
SLF2_AD_mod_L <- lm(mean_AD_LSLF2_values ~ Group, data = SLF_DTI_metrics_data)
SLF3_AD_mod_L <- lm(mean_AD_LSLF3_values ~ Group, data = SLF_DTI_metrics_data)
SLF1_AD_mod_R <- lm(mean_AD_RSLF1_values ~ Group, data = SLF_DTI_metrics_data)
SLF2_AD_mod_R <- lm(mean_AD_RSLF2_values ~ Group, data = SLF_DTI_metrics_data)
SLF3_AD_mod_R <- lm(mean_AD_RSLF3_values ~ Group, data = SLF_DTI_metrics_data)
anova(SLF1_AD_mod_L)
anova(SLF2_AD_mod_L)
anova(SLF3_AD_mod_L)
anova(SLF1_AD_mod_R)
anova(SLF2_AD_mod_R)
anova(SLF3_AD_mod_R) #sig
#effect size omnibus ANOVA
etaSquared(SLF3_AD_mod_R)

#ANCOVA (w/ age)
SLF1_AD_covar_mod_L <- lm(mean_AD_LSLF1_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF2_AD_covar_mod_L <- lm(mean_AD_LSLF2_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF3_AD_covar_mod_L <- lm(mean_AD_LSLF3_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF1_AD_covar_mod_R <- lm(mean_AD_RSLF1_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF2_AD_covar_mod_R <- lm(mean_AD_RSLF2_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF3_AD_covar_mod_R <- lm(mean_AD_RSLF3_values ~ Group+Age, data = SLF_DTI_metrics_data)
anova(SLF1_AD_covar_mod_L)
anova(SLF2_AD_covar_mod_L)
anova(SLF3_AD_covar_mod_L)
anova(SLF1_AD_covar_mod_R)
anova(SLF2_AD_covar_mod_R)
anova(SLF3_AD_covar_mod_R) #sig

#for RD:
SLF1_RD_mod_L <- lm(mean_RD_LSLF1_values ~ Group, data = SLF_DTI_metrics_data)
SLF2_RD_mod_L <- lm(mean_RD_LSLF2_values ~ Group, data = SLF_DTI_metrics_data)
SLF3_RD_mod_L <- lm(mean_RD_LSLF3_values ~ Group, data = SLF_DTI_metrics_data)
SLF1_RD_mod_R <- lm(mean_RD_RSLF1_values ~ Group, data = SLF_DTI_metrics_data)
SLF2_RD_mod_R <- lm(mean_RD_RSLF2_values ~ Group, data = SLF_DTI_metrics_data)
SLF3_RD_mod_R <- lm(mean_RD_RSLF3_values ~ Group, data = SLF_DTI_metrics_data)
anova(SLF1_RD_mod_L)
anova(SLF2_RD_mod_L)
anova(SLF3_RD_mod_L)
anova(SLF1_RD_mod_R)
anova(SLF2_RD_mod_R)
anova(SLF3_RD_mod_R)
#ANCOVA (w/ age)
SLF1_RD_covar_mod_L <- lm(mean_RD_LSLF1_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF2_RD_covar_mod_L <- lm(mean_RD_LSLF2_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF3_RD_covar_mod_L <- lm(mean_RD_LSLF3_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF1_RD_covar_mod_R <- lm(mean_RD_RSLF1_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF2_RD_covar_mod_R <- lm(mean_RD_RSLF2_values ~ Group+Age, data = SLF_DTI_metrics_data)
SLF3_RD_covar_mod_R <- lm(mean_RD_RSLF3_values ~ Group+Age, data = SLF_DTI_metrics_data)
anova(SLF1_RD_covar_mod_L)
anova(SLF2_RD_covar_mod_L)
anova(SLF3_RD_covar_mod_L)
anova(SLF1_RD_covar_mod_R)
anova(SLF2_RD_covar_mod_R)
anova(SLF3_RD_covar_mod_R)

#post hoc test:
post_hoc_SLF3_AD_mod_R <- glht(SLF3_AD_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_AD_mod_R) #not sig
confint(post_hoc_SLF3_AD_mod_R)

#post hoc test (with age as covariate):
post_hoc_SLF3_AD_covar_mod_R <- glht(SLF3_AD_covar_mod_R, linfct = mcp(Group = "Tukey"))
summary(post_hoc_SLF3_AD_covar_mod_R) #not sig
confint(post_hoc_SLF3_AD_covar_mod_R)

#linear trend, for FA:
SLF1_L_FA_LinTrend_mod <- lm(mean_FA_LSLF1_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF1_L_FA_LinTrend_mod)
SLF2_L_FA_LinTrend_mod <- lm(mean_FA_LSLF2_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF2_L_FA_LinTrend_mod) 
SLF3_L_FA_LinTrend_mod <- lm(mean_FA_LSLF3_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF3_L_FA_LinTrend_mod)
SLF1_R_FA_LinTrend_mod <- lm(mean_FA_RSLF1_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF1_R_FA_LinTrend_mod)
SLF2_R_FA_LinTrend_mod <- lm(mean_FA_RSLF2_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF2_R_FA_LinTrend_mod)
SLF3_R_FA_LinTrend_mod <- lm(mean_FA_RSLF3_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF3_R_FA_LinTrend_mod)
#linear trend, for MD:
SLF1_L_MD_LinTrend_mod <- lm(mean_MD_LSLF1_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF1_L_MD_LinTrend_mod)
SLF2_L_MD_LinTrend_mod <- lm(mean_MD_LSLF2_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF2_L_MD_LinTrend_mod)
SLF3_L_MD_LinTrend_mod <- lm(mean_MD_LSLF3_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF3_L_MD_LinTrend_mod)
SLF1_R_MD_LinTrend_mod <- lm(mean_MD_RSLF1_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF1_R_MD_LinTrend_mod)
SLF2_R_MD_LinTrend_mod <- lm(mean_MD_RSLF2_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF2_R_MD_LinTrend_mod)
SLF3_R_MD_LinTrend_mod <- lm(mean_MD_RSLF3_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF3_R_MD_LinTrend_mod) #sig
#linear trend, for AD:
SLF1_L_AD_LinTrend_mod <- lm(mean_AD_LSLF1_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF1_L_AD_LinTrend_mod)
SLF2_L_AD_LinTrend_mod <- lm(mean_AD_LSLF2_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF2_L_AD_LinTrend_mod)
SLF3_L_AD_LinTrend_mod <- lm(mean_AD_LSLF3_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF3_L_AD_LinTrend_mod) #sig
SLF1_R_AD_LinTrend_mod <- lm(mean_AD_RSLF1_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF1_R_AD_LinTrend_mod)
SLF2_R_AD_LinTrend_mod <- lm(mean_AD_RSLF2_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF2_R_AD_LinTrend_mod)
SLF3_R_AD_LinTrend_mod <- lm(mean_AD_RSLF3_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF3_R_AD_LinTrend_mod) #sig
#linear trend, for RD:
SLF1_L_RD_LinTrend_mod <- lm(mean_RD_LSLF1_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF1_L_RD_LinTrend_mod)
SLF2_L_RD_LinTrend_mod <- lm(mean_RD_LSLF2_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF2_L_RD_LinTrend_mod)
SLF3_L_RD_LinTrend_mod <- lm(mean_RD_LSLF3_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF3_L_RD_LinTrend_mod)
SLF1_R_RD_LinTrend_mod <- lm(mean_RD_RSLF1_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF1_R_RD_LinTrend_mod)
SLF2_R_RD_LinTrend_mod <- lm(mean_RD_RSLF2_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF2_R_RD_LinTrend_mod)
SLF3_R_RD_LinTrend_mod <- lm(mean_RD_RSLF3_values ~ Trend_Group + Group, data = SLF_DTI_metrics_data)
anova(SLF3_R_RD_LinTrend_mod) #sig



#### -------------------------- Plotting -----------------------------------####
#for FA: just with left and right SLF 1, 2, 3
SLF_data_FA <- dplyr::select(SLF_DTI_metrics_data, 
                             ParticipantID,
                             Group,
                             mean_FA_LSLF1_values,
                             mean_FA_LSLF2_values,
                             mean_FA_LSLF3_values,
                             mean_FA_RSLF1_values,
                             mean_FA_RSLF2_values,
                             mean_FA_RSLF3_values)
SLF_data_FA_long <- gather(SLF_data_FA, 
                           "SLF_type",
                           "FA_metric",
                           mean_FA_LSLF1_values,
                           mean_FA_LSLF2_values,
                           mean_FA_LSLF3_values,
                           mean_FA_RSLF1_values,
                           mean_FA_RSLF2_values,
                           mean_FA_RSLF3_values)
#Separate the groups - All tracts SLF FA (raincloud plot)
ggplot(SLF_data_FA_long, aes(x = SLF_type, y = FA_metric, fill = Group)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fractional Anisotropy (FA)") +
  scale_x_discrete(labels = c("mean_FA_LSLF1_values" = "Left SLF1","mean_FA_LSLF2_values" = "Left SLF2","mean_FA_LSLF3_values" = "Left SLF3","mean_FA_RSLF1_values" = "Right SLF1","mean_FA_RSLF2_values" = "Right SLF2","mean_FA_RSLF3_values" = "Right SLF3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()

#for MD: just with left and right SLF 1, 2, 3
SLF_data_MD <- dplyr::select(SLF_DTI_metrics_data, 
                             ParticipantID,
                             Group,
                             mean_MD_LSLF1_values,
                             mean_MD_LSLF2_values,
                             mean_MD_LSLF3_values,
                             mean_MD_RSLF1_values,
                             mean_MD_RSLF2_values,
                             mean_MD_RSLF3_values)
SLF_data_MD_long <- gather(SLF_data_MD, 
                           "SLF_type",
                           "MD_metric",
                           mean_MD_LSLF1_values,
                           mean_MD_LSLF2_values,
                           mean_MD_LSLF3_values,
                           mean_MD_RSLF1_values,
                           mean_MD_RSLF2_values,
                           mean_MD_RSLF3_values)
#Separate the groups - All tracts SLF MD (raincloud plot)
ggplot(SLF_data_MD_long, aes(x = SLF_type, y = MD_metric, fill = Group)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Mean Diffusivity (MD)") +
  scale_x_discrete(labels = c("mean_MD_LSLF1_values" = "Left SLF1","mean_MD_LSLF2_values" = "Left SLF2","mean_MD_LSLF3_values" = "Left SLF3","mean_MD_RSLF1_values" = "Right SLF1","mean_MD_RSLF2_values" = "Right SLF2","mean_MD_RSLF3_values" = "Right SLF3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 18))+
  coord_flip()

#for AD: just with left and right SLF 1, 2, 3
SLF_data_AD <- dplyr::select(SLF_DTI_metrics_data, 
                             ParticipantID,
                             Group,
                             mean_AD_LSLF1_values,
                             mean_AD_LSLF2_values,
                             mean_AD_LSLF3_values,
                             mean_AD_RSLF1_values,
                             mean_AD_RSLF2_values,
                             mean_AD_RSLF3_values)
SLF_data_AD_long <- gather(SLF_data_AD, 
                           "SLF_type",
                           "AD_metric",
                           mean_AD_LSLF1_values,
                           mean_AD_LSLF2_values,
                           mean_AD_LSLF3_values,
                           mean_AD_RSLF1_values,
                           mean_AD_RSLF2_values,
                           mean_AD_RSLF3_values)
#Separate the groups - All tracts SLF AD (raincloud plot)
ggplot(SLF_data_AD_long, aes(x = SLF_type, y = AD_metric, fill = Group)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Axial Diffusivity (AD)") +
  scale_x_discrete(labels = c("mean_AD_LSLF1_values" = "Left SLF1","mean_AD_LSLF2_values" = "Left SLF2","mean_AD_LSLF3_values" = "Left SLF3","mean_AD_RSLF1_values" = "Right SLF1","mean_AD_RSLF2_values" = "Right SLF2","mean_AD_RSLF3_values" = "Right SLF3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 18))+
  coord_flip()

#for RD: just with left and right SLF 1, 2, 3
SLF_data_RD <- dplyr::select(SLF_DTI_metrics_data, 
                             ParticipantID,
                             Group,
                             mean_RD_LSLF1_values,
                             mean_RD_LSLF2_values,
                             mean_RD_LSLF3_values,
                             mean_RD_RSLF1_values,
                             mean_RD_RSLF2_values,
                             mean_RD_RSLF3_values)
SLF_data_RD_long <- gather(SLF_data_RD, 
                           "SLF_type",
                           "RD_metric",
                           mean_RD_LSLF1_values,
                           mean_RD_LSLF2_values,
                           mean_RD_LSLF3_values,
                           mean_RD_RSLF1_values,
                           mean_RD_RSLF2_values,
                           mean_RD_RSLF3_values)
#Separate the groups - All tracts SLF RD (raincloud plot)
ggplot(SLF_data_RD_long, aes(x = SLF_type, y = RD_metric, fill = Group)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Group), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Radial Diffusivity (RD)") +
  scale_x_discrete(labels = c("mean_RD_LSLF1_values" = "Left SLF1","mean_RD_LSLF2_values" = "Left SLF2","mean_RD_LSLF3_values" = "Left SLF3","mean_RD_RSLF1_values" = "Right SLF1","mean_RD_RSLF2_values" = "Right SLF2","mean_RD_RSLF3_values" = "Right SLF3")) + 
  scale_fill_discrete(labels=c("1" = "C", "2" = "SCD", "3"= "aMCI", "4"= "mMCI", "5" = "AD"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 18))+
  coord_flip()












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

#navigate to the correct pathway which contains the SLF DTI metric text files: 
#setwd('/yourpathway/')
setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/DTI_metrics/')

#for FA:
setwd('fractional_anisotropy_FA/')
#read in all of the left SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_LSLF1_values
#read in all of the left SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_LSLF2_values
#read in all of the left SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_LSLF3_values
#read in all of the right SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_RSLF1_values
#read in all of the right SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_RSLF2_values
#read in all of the right SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_FA_RSLF3_values

#for MD:
setwd('..')#go up one directory
setwd('mean_diffusivity_MD/')
#read in all of the left SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_LSLF1_values
#read in all of the left SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_LSLF2_values
#read in all of the left SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_LSLF3_values
#read in all of the right SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_RSLF1_values
#read in all of the right SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_RSLF2_values
#read in all of the right SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_MD_RSLF3_values

#for AD:
setwd('..')#go up one directory
setwd('axial_diffusivity_AD/')
#read in all of the left SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_LSLF1_values
#read in all of the left SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_LSLF2_values
#read in all of the left SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_LSLF3_values
#read in all of the right SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_RSLF1_values
#read in all of the right SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_RSLF2_values
#read in all of the right SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_AD_RSLF3_values

#for RD:
setwd('..')#go up one directory
setwd('radial_diffusivity_RD/')
#read in all of the left SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_LSLF1_values
#read in all of the left SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_LSLF2_values
#read in all of the left SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("LSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_LSLF3_values
#read in all of the right SLF 1 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF1"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_RSLF1_values
#read in all of the right SLF 2 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF2"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_RSLF2_values
#read in all of the right SLF 3 files
Filenames <- list.files(path = getwd(), 
                        pattern = ("RSLF3"),
                        full.names = T) 
#calculate the mean for each participant
map_dbl(Filenames, ~ vroom(.x) %>% 
          select(where(is.numeric)) %>% 
          unlist %>% mean(na.rm = TRUE)) -> mean_RD_RSLF3_values


#put all mean DTI metric values into a dataframe
SLF_DTI_metrics_data <- cbind.data.frame(mean_FA_LSLF1_values,mean_FA_LSLF2_values,mean_FA_LSLF3_values,
                                         mean_FA_RSLF1_values,mean_FA_RSLF2_values,mean_FA_RSLF3_values,
                                         mean_MD_LSLF1_values,mean_MD_LSLF2_values,mean_MD_LSLF3_values,
                                         mean_MD_RSLF1_values,mean_MD_RSLF2_values,mean_MD_RSLF3_values,
                                         mean_AD_LSLF1_values,mean_AD_LSLF2_values,mean_AD_LSLF3_values,
                                         mean_AD_RSLF1_values,mean_AD_RSLF2_values,mean_AD_RSLF3_values,
                                         mean_RD_LSLF1_values,mean_RD_LSLF2_values,mean_RD_LSLF3_values,
                                         mean_RD_RSLF1_values,mean_RD_RSLF2_values,mean_RD_RSLF3_values)
#add in Timepoint (alternating pattern for DTI metric analysis)
Timepoint <- rep(c('F0','F2'), 124)
SLF_DTI_metrics_data$Timepoint<-Timepoint
#add in Participant ID in alternating order for DTI metric analysis
ParticipantID_SLF <- sort(as.character(DPRC_neuropsych_data$ParticipantID))
SLF_DTI_metrics_data$ParticipantID_SLF<-ParticipantID_SLF
#sort dataframe by Time Point + PT ID
SLF_DTI_metrics_data <- SLF_DTI_metrics_data[order(SLF_DTI_metrics_data$Timepoint),]
#rename Participant ID 
names(SLF_DTI_metrics_data)[names(SLF_DTI_metrics_data)=='ParticipantID_SLF'] <- 'ParticipantID'
#convert to factor variables
SLF_DTI_metrics_data$ParticipantID <- as.factor(SLF_DTI_metrics_data$ParticipantID)
SLF_DTI_metrics_data$Timepoint <- as.factor(SLF_DTI_metrics_data$Timepoint)
#add in Individual number
SLF_DTI_metrics_data$Individual_number <- DPRC_neuropsych_data$Individual_number
#add in Group classification column for each participant - data from another worksheet
SLF_DTI_metrics_data$Group <- DPRC_neuropsych_data$Group
#add in covariates (clinical site) from the covariates_data dataframe to the SLF dataframe
SLF_DTI_metrics_data$ClinSite_name <- DPRC_neuropsych_data$Clinical_site 
SLF_DTI_metrics_data$Age <- DPRC_neuropsych_data$Age 
SLF_DTI_metrics_data$Sex <- DPRC_neuropsych_data$Sex
SLF_DTI_metrics_data$Sex_binary <- DPRC_neuropsych_data$Sex_binary


###----------------------------Descriptives----------------------------------####
#look at descriptive stats of the DTI metrics between groups and over time for F0 and F2 time points
#FA
#Left SLF 1
SLF1_L_FA_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_LSLF1_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF1_L_FA_descrip <-as.data.frame.matrix(SLF1_L_FA_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF1_L_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF1_L_FA <- na.omit(F0_SLF1_L_FA$mean_FA_LSLF1_values)
mean(noNAsF0_SLF1_L_FA)
sd(noNAsF0_SLF1_L_FA)
#F2
F2_SLF1_L_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF1_L_FA <- na.omit(F2_SLF1_L_FA$mean_FA_LSLF1_values)
mean(noNAsF2_SLF1_L_FA)
sd(noNAsF2_SLF1_L_FA)
#Left SLF 2
SLF2_L_FA_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_LSLF2_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF2_L_FA_descrip <-as.data.frame.matrix(SLF2_L_FA_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF2_L_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF2_L_FA <- na.omit(F0_SLF2_L_FA$mean_FA_LSLF2_values)
mean(noNAsF0_SLF2_L_FA)
sd(noNAsF0_SLF2_L_FA)
#F2
F2_SLF2_L_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF2_L_FA <- na.omit(F2_SLF2_L_FA$mean_FA_LSLF2_values)
mean(noNAsF2_SLF2_L_FA)
sd(noNAsF2_SLF2_L_FA)
#Left SLF 3
SLF3_L_FA_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_LSLF3_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF3_L_FA_descrip <-as.data.frame.matrix(SLF3_L_FA_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF3_L_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF3_L_FA <- na.omit(F0_SLF3_L_FA$mean_FA_LSLF3_values)
mean(noNAsF0_SLF3_L_FA)
sd(noNAsF0_SLF3_L_FA)
#F2
F2_SLF3_L_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF3_L_FA <- na.omit(F2_SLF3_L_FA$mean_FA_LSLF3_values)
mean(noNAsF2_SLF3_L_FA)
sd(noNAsF2_SLF3_L_FA)
#Right SLF 1
SLF1_R_FA_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_RSLF1_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF1_R_FA_descrip <-as.data.frame.matrix(SLF1_R_FA_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF1_R_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF1_R_FA <- na.omit(F0_SLF1_R_FA$mean_FA_RSLF1_values)
mean(noNAsF0_SLF1_R_FA)
sd(noNAsF0_SLF1_R_FA)
#F2
F2_SLF1_R_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF1_R_FA <- na.omit(F2_SLF1_R_FA$mean_FA_RSLF1_values)
mean(noNAsF2_SLF1_R_FA)
sd(noNAsF2_SLF1_R_FA)
#Right SLF 2
SLF2_R_FA_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_RSLF2_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF2_R_FA_descrip <-as.data.frame.matrix(SLF2_R_FA_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF2_R_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF2_R_FA <- na.omit(F0_SLF2_R_FA$mean_FA_RSLF2_values)
mean(noNAsF0_SLF2_R_FA)
sd(noNAsF0_SLF2_R_FA)
#F2
F2_SLF2_R_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF2_R_FA <- na.omit(F2_SLF2_R_FA$mean_FA_RSLF2_values)
mean(noNAsF2_SLF2_R_FA)
sd(noNAsF2_SLF2_R_FA)
#Right SLF 3
SLF3_R_FA_descrip <- describeBy(SLF_DTI_metrics_data$mean_FA_RSLF3_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF3_R_FA_descrip <-as.data.frame.matrix(SLF3_R_FA_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF3_R_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF3_R_FA <- na.omit(F0_SLF3_R_FA$mean_FA_RSLF3_values)
mean(noNAsF0_SLF3_R_FA)
sd(noNAsF0_SLF3_R_FA)
#F2
F2_SLF3_R_FA <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF3_R_FA <- na.omit(F2_SLF3_R_FA$mean_FA_RSLF3_values)
mean(noNAsF2_SLF3_R_FA)
sd(noNAsF2_SLF3_R_FA)

#MD
#Left SLF 1
SLF1_L_MD_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_LSLF1_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF1_L_MD_descrip <-as.data.frame.matrix(SLF1_L_MD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF1_L_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF1_L_MD <- na.omit(F0_SLF1_L_MD$mean_MD_LSLF1_values)
mean(noNAsF0_SLF1_L_MD)
sd(noNAsF0_SLF1_L_MD)
#F2
F2_SLF1_L_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF1_L_MD <- na.omit(F2_SLF1_L_MD$mean_MD_LSLF1_values)
mean(noNAsF2_SLF1_L_MD)
sd(noNAsF2_SLF1_L_MD)
#Left SLF 2
SLF2_L_MD_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_LSLF2_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF2_L_MD_descrip <-as.data.frame.matrix(SLF2_L_MD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF2_L_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF2_L_MD <- na.omit(F0_SLF2_L_MD$mean_MD_LSLF2_values)
mean(noNAsF0_SLF2_L_MD)
sd(noNAsF0_SLF2_L_MD)
#F2
F2_SLF2_L_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF2_L_MD <- na.omit(F2_SLF2_L_MD$mean_MD_LSLF2_values)
mean(noNAsF2_SLF2_L_MD)
sd(noNAsF2_SLF2_L_MD)
#Left SLF 3
SLF3_L_MD_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_LSLF3_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF3_L_MD_descrip <-as.data.frame.matrix(SLF3_L_MD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF3_L_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF3_L_MD <- na.omit(F0_SLF3_L_MD$mean_MD_LSLF3_values)
mean(noNAsF0_SLF3_L_MD)
sd(noNAsF0_SLF3_L_MD)
#F2
F2_SLF3_L_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF3_L_MD <- na.omit(F2_SLF3_L_MD$mean_MD_LSLF3_values)
mean(noNAsF2_SLF3_L_MD)
sd(noNAsF2_SLF3_L_MD)
#Right SLF 1
SLF1_R_MD_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_RSLF1_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF1_R_MD_descrip <-as.data.frame.matrix(SLF1_R_MD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF1_R_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF1_R_MD <- na.omit(F0_SLF1_R_MD$mean_MD_RSLF1_values)
mean(noNAsF0_SLF1_R_MD)
sd(noNAsF0_SLF1_R_MD)
#F2
F2_SLF1_R_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF1_R_MD <- na.omit(F2_SLF1_R_MD$mean_MD_RSLF1_values)
mean(noNAsF2_SLF1_R_MD)
sd(noNAsF2_SLF1_R_MD)
#Right SLF 2
SLF2_R_MD_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_RSLF2_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF2_R_MD_descrip <-as.data.frame.matrix(SLF2_R_MD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF2_R_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF2_R_MD <- na.omit(F0_SLF2_R_MD$mean_MD_RSLF2_values)
mean(noNAsF0_SLF2_R_MD)
sd(noNAsF0_SLF2_R_MD)
#F2
F2_SLF2_R_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF2_R_MD <- na.omit(F2_SLF2_R_MD$mean_MD_RSLF2_values)
mean(noNAsF2_SLF2_R_MD)
sd(noNAsF2_SLF2_R_MD)
#Right SLF 3
SLF3_R_MD_descrip <- describeBy(SLF_DTI_metrics_data$mean_MD_RSLF3_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF3_R_MD_descrip <-as.data.frame.matrix(SLF3_R_MD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF3_R_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF3_R_MD <- na.omit(F0_SLF3_R_MD$mean_MD_RSLF3_values)
mean(noNAsF0_SLF3_R_MD)
sd(noNAsF0_SLF3_R_MD)
#F2
F2_SLF3_R_MD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF3_R_MD <- na.omit(F2_SLF3_R_MD$mean_MD_RSLF3_values)
mean(noNAsF2_SLF3_R_MD)
sd(noNAsF2_SLF3_R_MD)

#AD
#Left SLF 1
SLF1_L_AD_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_LSLF1_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF1_L_AD_descrip <-as.data.frame.matrix(SLF1_L_AD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF1_L_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF1_L_AD <- na.omit(F0_SLF1_L_AD$mean_AD_LSLF1_values)
mean(noNAsF0_SLF1_L_AD)
sd(noNAsF0_SLF1_L_AD)
#F2
F2_SLF1_L_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF1_L_AD <- na.omit(F2_SLF1_L_AD$mean_AD_LSLF1_values)
mean(noNAsF2_SLF1_L_AD)
sd(noNAsF2_SLF1_L_AD)
#Left SLF 2
SLF2_L_AD_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_LSLF2_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF2_L_AD_descrip <-as.data.frame.matrix(SLF2_L_AD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF2_L_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF2_L_AD <- na.omit(F0_SLF2_L_AD$mean_AD_LSLF2_values)
mean(noNAsF0_SLF2_L_AD)
sd(noNAsF0_SLF2_L_AD)
#F2
F2_SLF2_L_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF2_L_AD <- na.omit(F2_SLF2_L_AD$mean_AD_LSLF2_values)
mean(noNAsF2_SLF2_L_AD)
sd(noNAsF2_SLF2_L_AD)
#Left SLF 3
SLF3_L_AD_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_LSLF3_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF3_L_AD_descrip <-as.data.frame.matrix(SLF3_L_AD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF3_L_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF3_L_AD <- na.omit(F0_SLF3_L_AD$mean_AD_LSLF3_values)
mean(noNAsF0_SLF3_L_AD)
sd(noNAsF0_SLF3_L_AD)
#F2
F2_SLF3_L_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF3_L_AD <- na.omit(F2_SLF3_L_AD$mean_AD_LSLF3_values)
mean(noNAsF2_SLF3_L_AD)
sd(noNAsF2_SLF3_L_AD)
#Right SLF 1
SLF1_R_AD_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_RSLF1_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF1_R_AD_descrip <-as.data.frame.matrix(SLF1_R_AD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF1_R_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF1_R_AD <- na.omit(F0_SLF1_R_AD$mean_AD_RSLF1_values)
mean(noNAsF0_SLF1_R_AD)
sd(noNAsF0_SLF1_R_AD)
#F2
F2_SLF1_R_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF1_R_AD <- na.omit(F2_SLF1_R_AD$mean_AD_RSLF1_values)
mean(noNAsF2_SLF1_R_AD)
sd(noNAsF2_SLF1_R_AD)
#Right SLF 2
SLF2_R_AD_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_RSLF2_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF2_R_AD_descrip <-as.data.frame.matrix(SLF2_R_AD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF2_R_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF2_R_AD <- na.omit(F0_SLF2_R_AD$mean_AD_RSLF2_values)
mean(noNAsF0_SLF2_R_AD)
sd(noNAsF0_SLF2_R_AD)
#F2
F2_SLF2_R_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF2_R_AD <- na.omit(F2_SLF2_R_AD$mean_AD_RSLF2_values)
mean(noNAsF2_SLF2_R_AD)
sd(noNAsF2_SLF2_R_AD)
#Right SLF 3
SLF3_R_AD_descrip <- describeBy(SLF_DTI_metrics_data$mean_AD_RSLF3_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF3_R_AD_descrip <-as.data.frame.matrix(SLF3_R_AD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF3_R_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF3_R_AD <- na.omit(F0_SLF3_R_AD$mean_AD_RSLF3_values)
mean(noNAsF0_SLF3_R_AD)
sd(noNAsF0_SLF3_R_AD)
#F2
F2_SLF3_R_AD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF3_R_AD <- na.omit(F2_SLF3_R_AD$mean_AD_RSLF3_values)
mean(noNAsF2_SLF3_R_AD)
sd(noNAsF2_SLF3_R_AD)

#RD
#Left SLF 1
SLF1_L_RD_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_LSLF1_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF1_L_RD_descrip <-as.data.frame.matrix(SLF1_L_RD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF1_L_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF1_L_RD <- na.omit(F0_SLF1_L_RD$mean_RD_LSLF1_values)
mean(noNAsF0_SLF1_L_RD)
sd(noNAsF0_SLF1_L_RD)
#F2
F2_SLF1_L_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF1_L_RD <- na.omit(F2_SLF1_L_RD$mean_RD_LSLF1_values)
mean(noNAsF2_SLF1_L_RD)
sd(noNAsF2_SLF1_L_RD)
#Left SLF 2
SLF2_L_RD_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_LSLF2_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF2_L_RD_descrip <-as.data.frame.matrix(SLF2_L_RD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF2_L_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF2_L_RD <- na.omit(F0_SLF2_L_RD$mean_RD_LSLF2_values)
mean(noNAsF0_SLF2_L_RD)
sd(noNAsF0_SLF2_L_RD)
#F2
F2_SLF2_L_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF2_L_RD <- na.omit(F2_SLF2_L_RD$mean_RD_LSLF2_values)
mean(noNAsF2_SLF2_L_RD)
sd(noNAsF2_SLF2_L_RD)
#Left SLF 3
SLF3_L_RD_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_LSLF3_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF3_L_RD_descrip <-as.data.frame.matrix(SLF3_L_RD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF3_L_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF3_L_RD <- na.omit(F0_SLF3_L_RD$mean_RD_LSLF3_values)
mean(noNAsF0_SLF3_L_RD)
sd(noNAsF0_SLF3_L_RD)
#F2
F2_SLF3_L_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF3_L_RD <- na.omit(F2_SLF3_L_RD$mean_RD_LSLF3_values)
mean(noNAsF2_SLF3_L_RD)
sd(noNAsF2_SLF3_L_RD)
#Right SLF 1
SLF1_R_RD_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_RSLF1_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF1_R_RD_descrip <-as.data.frame.matrix(SLF1_R_RD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF1_R_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF1_R_RD <- na.omit(F0_SLF1_R_RD$mean_RD_RSLF1_values)
mean(noNAsF0_SLF1_R_RD)
sd(noNAsF0_SLF1_R_RD)
#F2
F2_SLF1_R_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF1_R_RD <- na.omit(F2_SLF1_R_RD$mean_RD_RSLF1_values)
mean(noNAsF2_SLF1_R_RD)
sd(noNAsF2_SLF1_R_RD)
#Right SLF 2
SLF2_R_RD_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_RSLF2_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF2_R_RD_descrip <-as.data.frame.matrix(SLF2_R_RD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF2_R_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF2_R_RD <- na.omit(F0_SLF2_R_RD$mean_RD_RSLF2_values)
mean(noNAsF0_SLF2_R_RD)
sd(noNAsF0_SLF2_R_RD)
#F2
F2_SLF2_R_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF2_R_RD <- na.omit(F2_SLF2_R_RD$mean_RD_RSLF2_values)
mean(noNAsF2_SLF2_R_RD)
sd(noNAsF2_SLF2_R_RD)
#Right SLF 3
SLF3_R_RD_descrip <- describeBy(SLF_DTI_metrics_data$mean_RD_RSLF3_values, list(SLF_DTI_metrics_data$Group, SLF_DTI_metrics_data$Timepoint))
SLF3_R_RD_descrip <-as.data.frame.matrix(SLF3_R_RD_descrip)
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SLF3_R_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F0",]
noNAsF0_SLF3_R_RD <- na.omit(F0_SLF3_R_RD$mean_RD_RSLF3_values)
mean(noNAsF0_SLF3_R_RD)
sd(noNAsF0_SLF3_R_RD)
#F2
F2_SLF3_R_RD <- SLF_DTI_metrics_data[SLF_DTI_metrics_data[, "Timepoint"] == "F2",]
noNAsF2_SLF3_R_RD <- na.omit(F2_SLF3_R_RD$mean_RD_RSLF3_values)
mean(noNAsF2_SLF3_R_RD)
sd(noNAsF2_SLF3_R_RD)



###------------------ ANOVA & ANCOVA (age) + (age + sex) ---------------------------####
#run ANOVA to see if there are significant differences between groups: run mixed design, 2 x 5 ANOVA

#for FA:
aov_SLF1_L_FA_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_FA_LSLF1_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_L_FA_mod)
aov_SLF2_L_FA_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_FA_LSLF2_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_FA_mod)
aov_SLF3_L_FA_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_FA_LSLF3_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_FA_mod)
aov_SLF1_R_FA_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_FA_RSLF1_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_R_FA_mod)
aov_SLF2_R_FA_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_FA_RSLF2_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_R_FA_mod)
aov_SLF3_R_FA_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_FA_RSLF3_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_R_FA_mod)
#ANCOVA (age)
aov_SLF1_L_FA_covar_mod<- aov(mean_FA_LSLF1_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF1_L_FA_covar_mod)
aov_SLF2_L_FA_covar_mod<- aov(mean_FA_LSLF2_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF2_L_FA_covar_mod)
aov_SLF3_L_FA_covar_mod<- aov(mean_FA_LSLF3_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF3_L_FA_covar_mod)
aov_SLF1_R_FA_covar_mod<- aov(mean_FA_RSLF1_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF1_R_FA_covar_mod)
aov_SLF2_R_FA_covar_mod<- aov(mean_FA_RSLF2_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF2_R_FA_covar_mod)
aov_SLF3_R_FA_covar_mod<- aov(mean_FA_RSLF3_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF3_R_FA_covar_mod)

#for MD:
aov_SLF1_L_MD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_MD_LSLF1_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_L_MD_mod)
aov_SLF2_L_MD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_MD_LSLF2_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_MD_mod)
aov_SLF3_L_MD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_MD_LSLF3_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_MD_mod)
aov_SLF1_R_MD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_MD_RSLF1_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_R_MD_mod)
aov_SLF2_R_MD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_MD_RSLF2_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_R_MD_mod) #sig. Timepoint
aov_SLF3_R_MD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_MD_RSLF3_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_R_MD_mod)
#ANCOVA (age)
aov_SLF1_L_MD_covar_mod<- aov(mean_MD_LSLF1_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF1_L_MD_covar_mod)
aov_SLF2_L_MD_covar_mod<- aov(mean_MD_LSLF2_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF2_L_MD_covar_mod)
aov_SLF3_L_MD_covar_mod<- aov(mean_MD_LSLF3_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF3_L_MD_covar_mod)
aov_SLF1_R_MD_covar_mod<- aov(mean_MD_RSLF1_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF1_R_MD_covar_mod)
aov_SLF2_R_MD_covar_mod<- aov(mean_MD_RSLF2_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF2_R_MD_covar_mod) #sig. Timepoint
aov_SLF3_R_MD_covar_mod<- aov(mean_MD_RSLF3_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF3_R_MD_covar_mod)

#for AD:
aov_SLF1_L_AD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_AD_LSLF1_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_L_AD_mod)
aov_SLF2_L_AD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_AD_LSLF2_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_AD_mod) #sig. Timepoint
aov_SLF3_L_AD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_AD_LSLF3_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_AD_mod)
aov_SLF1_R_AD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_AD_RSLF1_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_R_AD_mod)
aov_SLF2_R_AD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_AD_RSLF2_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_R_AD_mod) #sig. Timepoint
aov_SLF3_R_AD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_AD_RSLF3_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_R_AD_mod)
#ANCOVA (age)
aov_SLF1_L_AD_covar_mod<- aov(mean_AD_LSLF1_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF1_L_AD_covar_mod)
aov_SLF2_L_AD_covar_mod<- aov(mean_AD_LSLF2_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF2_L_AD_covar_mod)
aov_SLF3_L_AD_covar_mod<- aov(mean_AD_LSLF3_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF3_L_AD_covar_mod)
aov_SLF1_R_AD_covar_mod<- aov(mean_AD_RSLF1_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF1_R_AD_covar_mod)
aov_SLF2_R_AD_covar_mod<- aov(mean_AD_RSLF2_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF2_R_AD_covar_mod) #sig. Timepoint
aov_SLF3_R_AD_covar_mod<- aov(mean_AD_RSLF3_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF3_R_AD_covar_mod)

#for RD:
aov_SLF1_L_RD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_RD_LSLF1_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_L_RD_mod)
aov_SLF2_L_RD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_RD_LSLF2_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_L_RD_mod)
aov_SLF3_L_RD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_RD_LSLF3_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_L_RD_mod)
aov_SLF1_R_RD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_RD_RSLF1_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF1_R_RD_mod)
aov_SLF2_R_RD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_RD_RSLF2_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF2_R_RD_mod) #sig. Timepoint
aov_SLF3_R_RD_mod <- anova_test(data=SLF_DTI_metrics_data, dv=mean_RD_RSLF3_values, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SLF3_R_RD_mod)
#ANCOVA (age)
aov_SLF1_L_RD_covar_mod<- aov(mean_RD_LSLF1_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF1_L_RD_covar_mod)
aov_SLF2_L_RD_covar_mod<- aov(mean_RD_LSLF2_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF2_L_RD_covar_mod)
aov_SLF3_L_RD_covar_mod<- aov(mean_RD_LSLF3_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF3_L_RD_covar_mod)
aov_SLF1_R_RD_covar_mod<- aov(mean_RD_RSLF1_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF1_R_RD_covar_mod)
aov_SLF2_R_RD_covar_mod<- aov(mean_RD_RSLF2_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF2_R_RD_covar_mod) #sig. Timepoint
aov_SLF3_R_RD_covar_mod<- aov(mean_RD_RSLF3_values ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age, data=SLF_DTI_metrics_data)
summary(aov_SLF3_R_RD_covar_mod)

#post-hoc follow up tests with sig. main effects and interactions: only Timepoint sig. different, so no need for 
#f/u tests



#### ---------------------------- Plotting -------------------------------- ####

#Combine all sig. plots on one plot 
#put data into long format for FA:
SLF_DTI_metrics_data_FA <- dplyr::select(SLF_DTI_metrics_data, 
                             ParticipantID,
                             Group,
                             Timepoint,
                             mean_FA_LSLF1_values,
                             mean_FA_LSLF2_values,
                             mean_FA_LSLF3_values,
                             mean_FA_RSLF1_values,
                             mean_FA_RSLF2_values,
                             mean_FA_RSLF3_values)
SLF_DTI_metrics_data_FA_long <- gather(SLF_DTI_metrics_data_FA, 
                           "SLF_type",
                           "FA_metric",
                           mean_FA_LSLF1_values,
                           mean_FA_LSLF2_values,
                           mean_FA_LSLF3_values,
                           mean_FA_RSLF1_values,
                           mean_FA_RSLF2_values,
                           mean_FA_RSLF3_values)
#Separate the groups - All tracts SLF FC (raincloud plot)
ggplot(SLF_DTI_metrics_data_FA_long, aes(x = SLF_type, y = FA_metric, fill = Timepoint)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Fractional Anisotropy (FA)") +
  scale_x_discrete(labels = c("mean_FA_LSLF1_values" = "Left SLF1","mean_FA_LSLF2_values" = "Left SLF2","mean_FA_LSLF3_values" = "Left SLF3","mean_FA_RSLF1_values" = "Right SLF1","mean_FA_RSLF2_values" = "Right SLF2","mean_FA_RSLF3_values" = "Right SLF3")) + 
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))+
  coord_flip()

#put data into long format for MD:
SLF_DTI_metrics_data_MD <- dplyr::select(SLF_DTI_metrics_data, 
                                         ParticipantID,
                                         Group,
                                         Timepoint,
                                         mean_MD_LSLF1_values,
                                         mean_MD_LSLF2_values,
                                         mean_MD_LSLF3_values,
                                         mean_MD_RSLF1_values,
                                         mean_MD_RSLF2_values,
                                         mean_MD_RSLF3_values)
SLF_DTI_metrics_data_MD_long <- gather(SLF_DTI_metrics_data_MD, 
                                       "SLF_type",
                                       "MD_metric",
                                       mean_MD_LSLF1_values,
                                       mean_MD_LSLF2_values,
                                       mean_MD_LSLF3_values,
                                       mean_MD_RSLF1_values,
                                       mean_MD_RSLF2_values,
                                       mean_MD_RSLF3_values)
#Separate the groups - All tracts SLF FC (raincloud plot)
ggplot(SLF_DTI_metrics_data_MD_long, aes(x = SLF_type, y = MD_metric, fill = Timepoint)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Mean Diffusivity (MD)") +
  scale_x_discrete(labels = c("mean_MD_LSLF1_values" = "Left SLF1","mean_MD_LSLF2_values" = "Left SLF2","mean_MD_LSLF3_values" = "Left SLF3","mean_MD_RSLF1_values" = "Right SLF1","mean_MD_RSLF2_values" = "Right SLF2","mean_MD_RSLF3_values" = "Right SLF3")) + 
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 18))+
  coord_flip()

#put data into long format for AD:
SLF_DTI_metrics_data_AD <- dplyr::select(SLF_DTI_metrics_data, 
                                         ParticipantID,
                                         Group,
                                         Timepoint,
                                         mean_AD_LSLF1_values,
                                         mean_AD_LSLF2_values,
                                         mean_AD_LSLF3_values,
                                         mean_AD_RSLF1_values,
                                         mean_AD_RSLF2_values,
                                         mean_AD_RSLF3_values)
SLF_DTI_metrics_data_AD_long <- gather(SLF_DTI_metrics_data_AD, 
                                       "SLF_type",
                                       "AD_metric",
                                       mean_AD_LSLF1_values,
                                       mean_AD_LSLF2_values,
                                       mean_AD_LSLF3_values,
                                       mean_AD_RSLF1_values,
                                       mean_AD_RSLF2_values,
                                       mean_AD_RSLF3_values)
#Separate the groups - All tracts SLF FC (raincloud plot)
ggplot(SLF_DTI_metrics_data_AD_long, aes(x = SLF_type, y = AD_metric, fill = Timepoint)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Mean Diffusivity (AD)") +
  scale_x_discrete(labels = c("mean_AD_LSLF1_values" = "Left SLF1","mean_AD_LSLF2_values" = "Left SLF2","mean_AD_LSLF3_values" = "Left SLF3","mean_AD_RSLF1_values" = "Right SLF1","mean_AD_RSLF2_values" = "Right SLF2","mean_AD_RSLF3_values" = "Right SLF3")) + 
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 18))+
  coord_flip()

#put data into long format for RD:
SLF_DTI_metrics_data_RD <- dplyr::select(SLF_DTI_metrics_data, 
                                         ParticipantID,
                                         Group,
                                         Timepoint,
                                         mean_RD_LSLF1_values,
                                         mean_RD_LSLF2_values,
                                         mean_RD_LSLF3_values,
                                         mean_RD_RSLF1_values,
                                         mean_RD_RSLF2_values,
                                         mean_RD_RSLF3_values)
SLF_DTI_metrics_data_RD_long <- gather(SLF_DTI_metrics_data_RD, 
                                       "SLF_type",
                                       "RD_metric",
                                       mean_RD_LSLF1_values,
                                       mean_RD_LSLF2_values,
                                       mean_RD_LSLF3_values,
                                       mean_RD_RSLF1_values,
                                       mean_RD_RSLF2_values,
                                       mean_RD_RSLF3_values)
#Separate the groups - All tracts SLF FC (raincloud plot)
ggplot(SLF_DTI_metrics_data_RD_long, aes(x = SLF_type, y = RD_metric, fill = Timepoint)) + 
  geom_flat_violin(scale="width", alpha = .8, position = position_dodge(width=0.8)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint), position = position_dodge(width=0.8)) + 
  guides(colour=FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  xlab("SLF Tract") + 
  ylab("Mean Diffusivity (RD)") +
  scale_x_discrete(labels = c("mean_RD_LSLF1_values" = "Left SLF1","mean_RD_LSLF2_values" = "Left SLF2","mean_RD_LSLF3_values" = "Left SLF3","mean_RD_RSLF1_values" = "Right SLF1","mean_RD_RSLF2_values" = "Right SLF2","mean_RD_RSLF3_values" = "Right SLF3")) + 
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 18))+
  coord_flip()







