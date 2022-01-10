#Analysis for the Travelling Heads (TH) for the Dementia Prevention Research 
#Clinic (DPRC) data. 

#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 29/11/21


#load libraries via pacman
pacman::p_load(psych, ggplot2, rstatix)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

# choose & set to your directory. This is where the TH files are.
#setwd('/yourpathway/')
setwd('V:/NECTAR_data/LENORE/derivatives/groups/TH/diff_data/cross-sectional/')  

#Create TH dataframe (since it is a small dataset).
ParticipantID <- c('HQC_ADPRC0001F4', 'HQC_ADPRC0002F5', 'HQC_ADPRC0004F3', 
                  'HQC_ADPRC0005F4', 'HQC_CDPRC0001F1', 'HQC_CDPRC0002F1', 
                  'HQC_CDPRC0004F1', 'HQC_CDPRC0005F1', 'HQC_DDPRC0001F2',
                  'HQC_DDPRC0002F2','HQC_DDPRC0004F2','HQC_DDPRC0005F1')

Age <- c(37.58333333, 40.41666667, 37.91666667, 55.5, 37.5, 40.41666667,38, 55.5, 37.5, 40.41666667, 38, 55.5)

Sex <- c('M', 'F', 'M', 'M', 'M', 'F', 'M', 'M', 'M', 'F', 'M', 'M')

Clinical_Site <- c('Auckland', 'Auckland', 'Auckland', 'Auckland', 
                   'Christchurch', 'Christchurch', 'Christchurch', 'Christchurch',
                   'Dunedin', 'Dunedin', 'Dunedin', 'Dunedin')

Individual_number <- c('1','2','3','4','1','2','3','4','1','2','3','4')

#Whole-brain FBA metrics
#FD
FD_Whole_brain <- read.delim('FD.txt', sep=' ')
FD_Whole_brain <- FD_Whole_brain[1:12, 1:7]
names(FD_Whole_brain)[1] <- 'mn'
names(FD_Whole_brain)[2] <- 'md'
names(FD_Whole_brain)[3] <- 'std'
names(FD_Whole_brain)[4] <- 'std_rv'
names(FD_Whole_brain)[5] <- 'min'
names(FD_Whole_brain)[6] <- 'max'
names(FD_Whole_brain)[7] <- 'count'
#FC
FC_Whole_brain <- read.delim('FC_log.txt', sep=' ')
FC_Whole_brain <- FC_Whole_brain[1:12, 1:7]
names(FC_Whole_brain)[1] <- 'mn'
names(FC_Whole_brain)[2] <- 'md'
names(FC_Whole_brain)[3] <- 'std'
names(FC_Whole_brain)[4] <- 'std_rv'
names(FC_Whole_brain)[5] <- 'min'
names(FC_Whole_brain)[6] <- 'max'
names(FC_Whole_brain)[7] <- 'count'
#FDC
FDC_Whole_brain <- read.delim('FDC.txt', sep=' ')
FDC_Whole_brain <- FDC_Whole_brain[1:12, 1:7]
names(FDC_Whole_brain)[1] <- 'mn'
names(FDC_Whole_brain)[2] <- 'md'
names(FDC_Whole_brain)[3] <- 'std'
names(FDC_Whole_brain)[4] <- 'std_rv'
names(FDC_Whole_brain)[5] <- 'min'
names(FDC_Whole_brain)[6] <- 'max'
names(FDC_Whole_brain)[7] <- 'count'

#SLF FBA metrics
#change directory to SLF results
setwd('V:/NECTAR_data/LENORE/derivatives/groups/TH/diff_data/cross-sectional/template/TOI/')  
#FD
FD_SLF <- read.delim('FD_SLF_whole_TOI.txt', sep=' ')
FD_SLF <- FD_SLF[1:12, 1:7]
names(FD_SLF)[1] <- 'mn'
names(FD_SLF)[2] <- 'md'
names(FD_SLF)[3] <- 'std'
names(FD_SLF)[4] <- 'std_rv'
names(FD_SLF)[5] <- 'min'
names(FD_SLF)[6] <- 'max'
names(FD_SLF)[7] <- 'count'
#FC
FC_SLF <- read.delim('FC_log_SLF_whole_TOI.txt', sep=' ')
FC_SLF <- FC_SLF[1:12, 1:7]
names(FC_SLF)[1] <- 'mn'
names(FC_SLF)[2] <- 'md'
names(FC_SLF)[3] <- 'std'
names(FC_SLF)[4] <- 'std_rv'
names(FC_SLF)[5] <- 'min'
names(FC_SLF)[6] <- 'max'
names(FC_SLF)[7] <- 'count'
#FDC
FDC_SLF <- read.delim('FDC_SLF_whole_TOI.txt', sep=' ')
FDC_SLF <- FDC_SLF[1:12, 1:7]
names(FDC_SLF)[1] <- 'mn'
names(FDC_SLF)[2] <- 'md'
names(FDC_SLF)[3] <- 'std'
names(FDC_SLF)[4] <- 'std_rv'
names(FDC_SLF)[5] <- 'min'
names(FDC_SLF)[6] <- 'max'
names(FDC_SLF)[7] <- 'count'

#Put all data together into a dataframe
TH_data <- data.frame(ParticipantID, Age, Sex, Clinical_Site, Individual_number, 
                      FD_Whole_brain$mn, FD_Whole_brain$std, FC_Whole_brain$mn, 
                      FC_Whole_brain$std, FDC_Whole_brain$mn, FDC_Whole_brain$std, 
                      FD_SLF$mn,FD_SLF$std, FC_SLF$mn, FC_SLF$std, FDC_SLF$mn, FDC_SLF$std)

#convert variables to factor
TH_data$ParticipantID <- as.factor(TH_data$ParticipantID)
TH_data$Clinical_Site <- as.factor(TH_data$Clinical_Site)
TH_data$Sex <- as.factor(TH_data$Sex)
TH_data$Individual_number <- as.factor(TH_data$Individual_number)

#Look at descriptive statistics
mean(TH_data$Age)
sd(TH_data$Age)
age_descrip <- describeBy(TH_data$Age, TH_data$Clinical_Site)
#gender_descrip <- by(TH_data$Clinical_Site, TH_data$Clinical_Site, summary)
#gender_descrip_detail <- describeBy(TH_data ~ Sex + Clinical_Site, skew=FALSE, ranges=FALSE)
#Whole-brain FBA descriptives
mean(TH_data$FD_Whole_brain.mn)
sd(TH_data$FD_Whole_brain.mn)
mean(TH_data$FC_Whole_brain.mn)
sd(TH_data$FC_Whole_brain.mn)
mean(TH_data$FDC_Whole_brain.mn)
sd(TH_data$FDC_Whole_brain.mn)
FD_whole_descrip <- describeBy(TH_data$FD_Whole_brain.mn, TH_data$Clinical_Site)
FC_whole_descrip <- describeBy(TH_data$FC_Whole_brain.mn, TH_data$Clinical_Site)
FDC_whole_descrip <- describeBy(TH_data$FDC_Whole_brain.mn, TH_data$Clinical_Site)
#SLF FBA descriptives
mean(TH_data$FD_SLF.mn)
sd(TH_data$FD_SLF.mn)
mean(TH_data$FC_SLF.mn)
sd(TH_data$FC_SLF.mn)
mean(TH_data$FDC_SLF.mn)
sd(TH_data$FDC_SLF.mn)
FD_SLF_descrip <- describeBy(TH_data$FD_SLF.mn, TH_data$Clinical_Site)
FC_SLF_descrip <- describeBy(TH_data$FC_SLF.mn, TH_data$Clinical_Site)
FDC_SLF_descrip <- describeBy(TH_data$FDC_SLF.mn, TH_data$Clinical_Site)

#Run statistic tests- repeated measures ANOVA
#check for significant difference in age between groups 
rep.aov_age_mod <- anova_test(data=TH_data, dv=Age, wid=Individual_number, within=Clinical_Site)
get_anova_table(rep.aov_age_mod)
#age_mod <- lm(Age ~ Clinical_Site, data = TH_data)
#age_mod <- lm(Age ~ 0 + Clinical_Site, data = TH_data) #test against y-intercept
#anova(age_mod)

#check for significant difference in gender between groups 
#reformat data for chi-square test
gender_data_chisq <- rbind(c(1,1,1), c(3,3,3))
colnames(gender_data_chisq) <- c("AKL", "ChCh", "Ddn")
rownames(gender_data_chisq) <- c("F", "M")
#run chi-square test
gender_chi_test <- chisq.test(gender_data_chisq)

# #check for significant difference in FBA metrics between groups
#whole-brain FBA
rep.aov_FD_whole_mod <- anova_test(data=TH_data, dv=FD_Whole_brain.mn, wid=Individual_number, within=Clinical_Site)
get_anova_table(rep.aov_FD_whole_mod)
rep.aov_FC_whole_mod <- anova_test(data=TH_data, dv=FC_Whole_brain.mn, wid=Individual_number, within=Clinical_Site)
get_anova_table(rep.aov_FC_whole_mod) #sig
rep.aov_FDC_whole_mod <- anova_test(data=TH_data, dv=FDC_Whole_brain.mn, wid=Individual_number, within=Clinical_Site)
get_anova_table(rep.aov_FDC_whole_mod)
#SLF FBA
rep.aov_FD_SLF_mod <- anova_test(data=TH_data, dv=FD_SLF.mn, wid=Individual_number, within=Clinical_Site)
get_anova_table(rep.aov_FD_SLF_mod)
rep.aov_FC_SLF_mod <- anova_test(data=TH_data, dv=FC_SLF.mn, wid=Individual_number, within=Clinical_Site)
get_anova_table(rep.aov_FC_SLF_mod)
rep.aov_FDC_SLF_mod <- anova_test(data=TH_data, dv=FDC_SLF.mn, wid=Individual_number, within=Clinical_Site)
get_anova_table(rep.aov_FDC_SLF_mod)

#post hoc 
#test by Time Point
TH_data %>%
  pairwise_t_test(
    FC_Whole_brain.mn ~ Clinical_Site, paired = TRUE,
    p.adjust.method = "fdr"
  )
#effect size for Time Point
TH_data%>%cohens_d(FC_Whole_brain.mn~Clinical_Site, paired=TRUE)

#plot data
#Whole brain FD (raincloud plot)
ggplot(TH_data, aes(x = Clinical_Site, y = FD_Whole_brain.mn, fill = Clinical_Site)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2) + 
  geom_line(aes(group = Individual_number, colour = Individual_number)) +
  geom_point(aes(y = FD_Whole_brain.mn, colour = Individual_number), position = position_jitter(width = .15), size = 5, alpha = 0.8) +
  xlab("Clinical Site") + 
  ylab("Whole-brain Fibre Density (FD)") +
  theme_classic()+
  theme(legend.position = "none") 
#Whole brain FC (raincloud plot)
ggplot(TH_data, aes(x = Clinical_Site, y = FC_Whole_brain.mn, fill = Clinical_Site)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2) + 
  geom_line(aes(group = Individual_number, colour = Individual_number)) +
  geom_point(aes(y = FC_Whole_brain.mn, colour = Individual_number), position = position_jitter(width = .15), size = 5, alpha = 0.8) +
  xlab("Clinical Site") + 
  ylab("Whole-brain Fibre Cross-section (FC)") +
  theme_classic()+
  theme(legend.position = "none") 
#Whole brain FDC (raincloud plot)
ggplot(TH_data, aes(x = Clinical_Site, y = FDC_Whole_brain.mn, fill = Clinical_Site)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2) + 
  geom_line(aes(group = Individual_number, colour = Individual_number)) +
  geom_point(aes(y = FDC_Whole_brain.mn, colour = Individual_number), position = position_jitter(width = .15), size = 5, alpha = 0.8) +
  xlab("Clinical Site") + 
  ylab("Whole-brain Fibre Density Cross-section (FDC)") +
  theme_classic()+
  theme(legend.position = "none") 
#SLF FD (raincloud plot)
ggplot(TH_data, aes(x = Clinical_Site, y = FD_SLF.mn, fill = Clinical_Site)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2) + 
  geom_line(aes(group = Individual_number, colour = Individual_number)) +
  geom_point(aes(y = FD_SLF.mn, colour = Individual_number), position = position_jitter(width = .15), size = 5, alpha = 0.8) +
  xlab("Clinical Site") + 
  ylab("Superior Longitudinal Fasciculus (SLF) Fibre Density (FD)") +
  theme_classic()+
  theme(legend.position = "none") 
#SLF FC (raincloud plot)
ggplot(TH_data, aes(x = Clinical_Site, y = FC_SLF.mn, fill = Clinical_Site)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2) + 
  geom_line(aes(group = Individual_number, colour = Individual_number)) +
  geom_point(aes(y = FC_SLF.mn, colour = Individual_number), position = position_jitter(width = .15), size = 5, alpha = 0.8) +
  xlab("Clinical Site") + 
  ylab("Superior Longitudinal Fasciculus (SLF) Fibre Cross-section (FC)") +
  theme_classic()+
  theme(legend.position = "none") 
#SLF FD (raincloud plot)
ggplot(TH_data, aes(x = Clinical_Site, y = FDC_SLF.mn, fill = Clinical_Site)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2) + 
  geom_line(aes(group = Individual_number, colour = Individual_number)) +
  geom_point(aes(y = FDC_SLF.mn, colour = Individual_number), position = position_jitter(width = .15), size = 5, alpha = 0.8) +
  xlab("Clinical Site") + 
  ylab("Superior Longitudinal Fasciculus (SLF) Fibre Density Cross-section (FDC)") +
  theme_classic()+
  theme(legend.position = "none") 



