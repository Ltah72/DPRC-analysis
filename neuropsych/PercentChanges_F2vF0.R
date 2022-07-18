#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 20/06/21

#load libraries via pacman
pacman::p_load(dplyr, ggplot2, psych, car, multcomp, lsr, BayesFactor, tidyr, GPArotation, corrplot, lsmeans, TukeyC, lme4, lmerTest, emmeans)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#set up pathway
#setwd('/yourpathway/')
#file.choose function
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')


#######----------------- Longitudinal (F0 vs. F2) analysis -----------------########

#read in csv files (participant file)
DPRC_neuropsych_data <- read.csv("longitudinal_DPRC_neuropsych_data_lined_up_valid_participants.csv")

#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'

#convert variables
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Sex_binary <- as.factor(DPRC_neuropsych_data$Sex_binary)
DPRC_neuropsych_data$Timepoint <- as.factor(DPRC_neuropsych_data$Timepoint)


#Percentage changes ((F2 - F0) / F0)
#put into long format
percent_change_data <- DPRC_neuropsych_data %>% 
  filter(Timepoint == "F2") %>%
  select(ParticipantID, 
         Age,
         Classification,
         Group, 
         Sex,
         Sex_binary,
         Timepoint,
         PercentChange_TMTA,
         PercentChange_TMTB,
         PercentChange_ColorNaming,
         PercentChange_WordReading,
         PercentChange_Inhibition,
         PercentChange_LetFluency,
         PercentChange_CatFluency,
         PercentChange_Switching,
         PercentChange_HayBTime1z,
         PercentChange_HayBTime2z,
         PercentChange_HayBCatAz,
         PercentChange_HayBCatBz)

#rename variables 
percent_change_data <- rename(percent_change_data, TMTA = PercentChange_TMTA)
percent_change_data <- rename(percent_change_data, TMTB = PercentChange_TMTB)
percent_change_data <- rename(percent_change_data, ColorNaming = PercentChange_ColorNaming)
percent_change_data <- rename(percent_change_data, WordReading = PercentChange_WordReading)
percent_change_data <- rename(percent_change_data, Inhibition = PercentChange_Inhibition)
percent_change_data <- rename(percent_change_data, LetFluency = PercentChange_LetFluency)
percent_change_data <- rename(percent_change_data, CatFluency = PercentChange_CatFluency)
percent_change_data <- rename(percent_change_data, Switching = PercentChange_Switching)
percent_change_data <- rename(percent_change_data, HayBTime1z = PercentChange_HayBTime1z)
percent_change_data <- rename(percent_change_data, HayBTime2z = PercentChange_HayBTime2z)
percent_change_data <- rename(percent_change_data, HayBCatAz = PercentChange_HayBCatAz)
percent_change_data <- rename(percent_change_data, HayBCatBz = PercentChange_HayBCatBz)

#put into long format
percent_change_data_long <- gather(percent_change_data, 
                                   "Test",
                                   "Percent_change", 
                                   TMTA,
                                   TMTB,
                                   ColorNaming,
                                   WordReading,
                                   Inhibition,
                                   LetFluency,
                                   CatFluency,
                                   Switching,
                                   HayBTime1z,
                                   HayBTime2z,
                                   HayBCatAz,
                                   HayBCatBz)
#plot percent changes (whole)
ggplot(percent_change_data_long, aes(x = Test, y = Percent_change)) + 
  geom_point(aes(y = Percent_change, color = Percent_change), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  #geom_boxplot(width = 0.1, fill = "white", outlier.size = .15, aes(colour = Percent_change)) + 
  xlab("Test") + 
  ylab("Percent Change") +
  theme_classic()+
  theme(legend.position = "none")


#put into long format (no Hayling)
percent_change_data_noHay_long <- gather(percent_change_data, 
                                   "Test",
                                   "Percent_change", 
                                   TMTA,
                                   TMTB,
                                   ColorNaming,
                                   WordReading,
                                   Inhibition,
                                   LetFluency,
                                   CatFluency,
                                   Switching)
#plot percent changes (no Hayling)
ggplot(percent_change_data_noHay_long, aes(x = Test, y = Percent_change)) + 
  geom_point(aes(y = Percent_change, color = Percent_change), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  #geom_boxplot(width = 0.1, fill = "white", outlier.size = .15, aes(colour = Percent_change)) + 
  xlab("Tests") + 
  ylab("Percent Change") +
  theme_classic() +
  theme(legend.position = "none")
                               
#put into long format (only Hayling)
percent_change_Hay_long <- gather(percent_change_data, 
                                   "Test",
                                   "Percent_change", 
                                   HayBTime1z,
                                   HayBTime2z,
                                   HayBCatAz,
                                   HayBCatBz)

#plot percent changes (only Hayling)
ggplot(percent_change_Hay_long, aes(x = Test, y = Percent_change)) + 
  geom_point(aes(y = Percent_change, color = Percent_change), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  #geom_boxplot(width = 0.1, fill = "white", outlier.size = .15, aes(colour = Percent_change)) + 
  xlab("Hayling's Sentence Completion Test") + 
  ylab("Percent Change") +
  theme_classic() +
  theme(legend.position = "none")


















#plot percent changes (no Hayling, only aMCI)
ggplot(subset(percent_change_data_noHay_long, Group %in% ("3")))+ 
  geom_point(aes(x = Test, y = Percent_change, color = Percent_change), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  #geom_text() +
  #geom_boxplot(width = 0.1, fill = "white", outlier.size = .15, aes(colour = Percent_change)) + 
  xlab("Tests") + 
  ylab("Percent Change") +
  theme_classic() +
  theme(legend.position = "none")




#plot by participant percent changes (no Hayling)
ggplot(percent_change_data_noHay_long, aes(x = ParticipantID, y = Percent_change)) + 
  geom_point(aes(y = Percent_change, color = Percent_change), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  #geom_boxplot(width = 0.1, fill = "white", outlier.size = .15, aes(colour = Percent_change)) + 
  xlab("Participant") + 
  ylab("Percent Change") +
  theme_classic() +
  theme(legend.position = "none")


