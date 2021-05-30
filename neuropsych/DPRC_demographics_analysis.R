#Analysing demographics information from the DPRC cohort. Here, we will be looking at the age and sex of the participants 
#and the clinical sites (Auckland, Christchurch and Dunedin) of where the participant MRI scans took place. We have 5 
#participant groups of interest, which are: Controls (C), subjective cognitive decline (SCD), amnestic mild cognitive 
#impairment (aMCI), multiple-domain mild cognitive impairment (mMCI), and Alzheimer's disease (AD).   


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 1/12/20


#load libraries via pacman
pacman::p_load(dplyr, ggplot2, psych, car, BayesFactor)

#load in pathway for functions
#source('path/tofile/here.R')
source('H:/ltah262/PhD/Diffusion/script/dprc/neuropsych/insertRow.R')

#set up pathway
#setwd('/yourpathway/')
setwd('V:/PartInfo')

#read in excel files (participant file)
DPRC_demographics <- read.csv("DPRC_subject_list.csv")
covariates_data <- read.csv("covariates-participants-lined-up_update.csv")

#rename some columns
names(DPRC_demographics)[1] <- "Subject_ID"
names(DPRC_demographics)[5] <- "Age"
names(DPRC_demographics)[6] <- "Sex"
names(DPRC_demographics)[7] <- "Sex_binary"
names(DPRC_demographics)[9] <- "Classification"
names(DPRC_demographics)[10] <- "ACE_score"


#make edits to include / exclude certain participants (MRI data exclusion) - this excel sheet updates, so you will need to update this code too. 
ADPRC0004F2_info <- as.data.frame(cbind('ADRPC_0004_F2', 1, '25/05/2018', '1/05/1939', 79.00, 'M', 0, 'AD(mod)', 5, 52, 'Note this was first scan. F1 ACE=73; F0 ACE=87', '', '', '', '')) #ADPRC0004F2
#use function below to insert new edited row into your existing dataframe
#insertRow(existingDF, newrow, r)
DPRC_demographics <- insertRow(DPRC_demographics, ADPRC0004F2_info, 4) #ADPRC0004F2 will be used as baseline (i.e. ADPRC0004F0)
#exclude participants because of warped mask from subject space to template space error
DPRC_demographics[198,'Classification'] <- 0 #CDPRC0002F0
DPRC_demographics[201,'Classification'] <- 0 #CDPRC0005F0
DPRC_demographics[211,'Classification'] <- 0 #CDPRC0019F0
DPRC_demographics[224,'Classification'] <- 0 #CDPRC0035F0
DPRC_demographics[225,'Classification'] <- 0 #CDPRC0036F0
DPRC_demographics[233,'Classification'] <- 0 #CDPRC0045F0
DPRC_demographics[239,'Classification'] <- 0 #CDPRC0051F0
DPRC_demographics[268,'Classification'] <- 0 #DDPRC0036F0
DPRC_demographics[269,'Classification'] <- 0 #DDPRC0037F0


#remove 0s (excluded) and -1s (not yet classified) from the data (classifications of participants).
DPRC_demographics_classification <- replace(DPRC_demographics$Classification, DPRC_demographics$Classification<=0, NA)
DPRC_demographics$Classification <- DPRC_demographics_classification


#convert categorical variables to a factor
DPRC_demographics$Classification <- as.factor(DPRC_demographics$Classification)
covariates_data$Group <- as.factor(covariates_data$Group)
DPRC_demographics$Sex_binary <- as.factor(DPRC_demographics$Sex_binary)

#convert variables to numeric 
DPRC_demographics$Age<- as.numeric(DPRC_demographics$Age)
DPRC_demographics$ACE_score<- as.numeric(DPRC_demographics$ACE_score)


#remove NAs from variables + put make a new dataset for these
age_NAs_omitted<- na.omit(DPRC_demographics$Age)
ACE_NAs_omitted<- na.omit(DPRC_demographics$ACE_score)
classification_NAs_omitted<- na.omit(DPRC_demographics$Classification)

#DPRC_demographic_NAs_omitted <- cbind(classification_NAs_omitted, age_NAs_omitted)

#DPRC_demographic_NAs_omitted <- cbind(classification_NAs_omitted, age_NAs_omitted, ACE_NAs_omitted)

#look at descriptive statistics
age_descrip <- describeBy(DPRC_demographics$Age, DPRC_demographics$Classification)
ACE_descrip <- describeBy(DPRC_demographics$ACE_score, DPRC_demographics$Classification)
gender_descrip <- describeBy(DPRC_demographics ~ Sex_binary + Classification, skew=FALSE, ranges=FALSE)
#clinsite_descrip <- describeBy(DPRC_demographics ~ Sex_binary + Classification, skew=FALSE, ranges=FALSE)


#plot the data to visualise



#plot age
ggplot(subset(DPRC_demographics, Classification %in% c("1", "2", "3", "4", "5")), aes(x = Classification, y = Age)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Classification)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Classification)) + 
    ylim(50, 95) +
    xlab("Group") + 
    ylab("Age") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Classification, colour = Classification), size = 1)

#plot gender
gender_data <- data.frame(table(covariates_data$Classification, covariates_data$Sex))
names(gender_data) <- c("Group", "Sex", "Count")
Group_order <- c("C", "SCD", "aMCI", "mMCI", "AD")
gender_data <- gender_data %>% arrange(factor(Group, levels=Group_order))
gender_data$Group <- factor(gender_data$Group, levels=c("C", "SCD", "aMCI", "mMCI", "AD"))

ggplot(data=gender_data, aes(x=Group, y=Count, fill=Sex)) +
    geom_bar(stat="identity") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#plot clinical site
location_data <- data.frame(table(covariates_data$Classification, covariates_data$Clinical_site))
names(location_data) <- c("Group", "Clinical_site", "Count")
location_data <- location_data %>% arrange(factor(Group, levels=Group_order))
location_data$Group <- factor(location_data$Group, levels=c("C", "SCD", "aMCI", "mMCI", "AD"))

ggplot(data=location_data, aes(x=Group, y=Count, fill=Clinical_site)) +
    geom_bar(stat="identity") +
    scale_fill_discrete(name = "Clinical Site") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



#plot ACE
ggplot(subset(DPRC_demographics, Classification %in% c("1", "2", "3", "4", "5")), aes(x = Classification, y = ACE_score)) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Classification)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Classification)) + 
    xlab("Group") + 
    ylab("ACE Score") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Classification, colour = Classification), size = 1)



#plot FBA metrics--
#for FD:
ggplot(subset(covariates_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Mean_FD)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density (FD)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")
#for FC:
ggplot(subset(covariates_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Mean_FC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Cross-section (FC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")
#for FDC:
ggplot(subset(covariates_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Mean_FDC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Group)) + 
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
    xlab("Group") + 
    ylab("Fibre Density Cross-section (FDC)") +
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
    theme_classic() +
    theme(legend.position = "none")



#check for significant difference in age between groups 
#need to exclude the naMCI group in the model
DPRC_demographics_classification <- replace(DPRC_demographics$Classification, DPRC_demographics$Classification==6, NA)
DPRC_demographics$Classification <- DPRC_demographics_classification

#convert continuous variables to a numeric
#DPRC_demographics$Age <- as.numeric(as.character(DPRC_demographics$Age))
age_mod <- lm(Age ~ Classification, data = DPRC_demographics)
anova(age_mod)
#anovaBF(Age ~ Classification, data = DPRC_demographics) #Bayesian (put into a new dataset b/c can't have any NA values)



#check for significant difference in gender between groups 
#reformat data for chi-square test
gender_data_chisq <- rbind(c(23,35,23,22,6), c(7,24,27,24,15))
colnames(gender_data_chisq) <- c("C", "SCD", "aMCI", "mMCI", "AD")
rownames(gender_data_chisq) <- c("F", "M")
#run chi-square test
gender_chi_test <- chisq.test(gender_data_chisq)


#check for significant difference in clinical site between groups 
#reformat data for chi-square test
location_data_chisq <- rbind(c(29,43,32,35,15), c(1,6,12,9,2), c(0,10,6,2,4))
colnames(location_data_chisq) <- c("C", "SCD", "aMCI", "mMCI", "AD")
rownames(location_data_chisq) <- c("Auckland", "Christchurch", "Dunedin")
#run chi-square test
location_chi_test <- chisq.test(location_data_chisq)
#if expected values in cells are too small, simulate more p-values:
#location_chi_test <- chisq.test(location_data_chisq, simulate.p.value = TRUE)


#check for significant difference in ACE between groups 
ACE_mod <- lm(ACE_score ~ Classification, data = DPRC_demographics)
anova(ACE_mod)




#check for significant difference in FBA metrics between groups 
FD_mod <- lm(Mean_FD ~ Group, data = covariates_data)
anova(FD_mod)
FC_mod <- lm(Mean_FC ~ Group, data = covariates_data)
anova(FC_mod)
FDC_mod <- lm(Mean_FDC ~ Group, data = covariates_data)
anova(FDC_mod)



#quick view of the data
plot(age_mod)


#Perform Levene's Test for homogenity of variances
leveneTest(Age ~ Classification, data = DPRC_demographics)


leveneTest(Mean_FD ~ Group, data = covariates_data)
leveneTest(Mean_FD ~ Sex, data = covariates_data)





#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(age_mod$residuals)


#check to see if variables have a linear relationship (e.g. age vs. white matter integrity metrics)
#for Mean FD
ggplot(subset(covariates_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Age, y = Mean_FD, shape = Group, color = Group)) +
    geom_point() +
    geom_smooth(method="lm", formula= y~x, se=FALSE) +
    ggtitle("Age vs. Fibre Density (FD)") +
    theme(plot.title = element_text(hjust=0.5))
    #scale_x_discrete(labels = c("1"="Control", "2"="SCD", "3"="aMCI", "4"="mMCI", "5"="AD")) 
    

#for Mean FC
ggplot(subset(covariates_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Age, y = Mean_FC, shape = Group, color = Group)) +
    geom_point() +
    geom_smooth(method="lm", formula= y~x, se=FALSE) +
    ggtitle("Age vs. Fibre Cross-section (FC)") +
    theme(plot.title = element_text(hjust=0.5))


#for Mean FDC
ggplot(subset(covariates_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Age, y = Mean_FDC, shape = Group, color = Group)) +
    geom_point() +
    geom_smooth(method="lm", formula= y~x, se=FALSE) +
    ggtitle("Age vs. Fibre Density Cross-section (FDC)") +
    theme(plot.title = element_text(hjust=0.5))




#all 3 metrics (FD, FC, FDC)
ggplot(subset(covariates_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Age)) +
    geom_point(aes(y=Mean_FD, color="Mean_FD")) + 
    geom_smooth(aes(y=Mean_FD, color="Mean_FD"), method="lm", formula= y~x, se=FALSE) +
    geom_point(aes(y=Mean_FC, color="Mean_FC")) + 
    geom_smooth(aes(y=Mean_FC, color="Mean_FC"), method="lm", formula= y~x, se=FALSE) +
    geom_point(aes(y=Mean_FDC, color="Mean_FDC")) + 
    geom_smooth(aes(y=Mean_FDC, color="Mean_FDC"), method="lm", formula= y~x, se=FALSE) +
    ylab("Fixel-Based Analysis (FBA) metrics") +
    ggtitle("Age vs. all 3 metrics (FD, FC, & FDC)") +
    theme(plot.title = element_text(hjust=0.5))



plot(covariates_data$Age, covariates_data$Mean_FD)



#check for significance with age and between age x group and white matter integrity (FBA metrics)
#run correlation + simple linear regression between age and white matter metrics
#for FD
Age_FBA_FD_cor <- cor.test(covariates_data$Age, covariates_data$Mean_FD)
Age_FBA_FD_mod <- lm(Mean_FD ~ Age, data = covariates_data)
summary(Age_FBA_FD_mod)
anova(Age_FBA_FD_mod)
#for FC
Age_FBA_FC_cor <- cor.test(covariates_data$Age, covariates_data$Mean_FC)
Age_FBA_FC_mod <- lm(Mean_FC ~ Age, data = covariates_data)
summary(Age_FBA_FC_mod)
#for FDC
Age_FBA_FDC_cor <- cor.test(covariates_data$Age, covariates_data$Mean_FDC)
Age_FBA_FDC_mod <- lm(Mean_FDC ~ Age, data = covariates_data)
summary(Age_FBA_FDC_mod)

#run multiple linear regression on age x group and white matter integrity (FBA metrics)
#for FD
AgexGroup_FD_mult_mod <- lm(Mean_FD ~ Age + Group + Age:Group, data = covariates_data)
summary(AgexGroup_FD_mult_mod)
anova(AgexGroup_FD_mult_mod)
#for FC
AgexGroup_FC_mult_mod <- lm(Mean_FC ~ Age + Group + Age:Group, data = covariates_data)
summary(AgexGroup_FC_mult_mod)
anova(AgexGroup_FC_mult_mod)
#for FDC
AgexGroup_FDC_mult_mod <- lm(Mean_FDC ~ Age + Group + Age:Group, data = covariates_data)
summary(AgexGroup_FDC_mult_mod)
anova(AgexGroup_FDC_mult_mod)



Sex_Group_FD_mod <- lm(Mean_FD ~ Group + Sex + Group:Sex, data = covariates_data)



#plot gender against white matter integrity metrics (FBA)
#for FD
ggplot(covariates_data, aes(x = Sex, y = Mean_FD)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Sex)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Sex)) +
    xlab("Sex") + 
    ylab("Mean Fibre Density (FD)")
#for FC
ggplot(covariates_data, aes(x = Sex, y = Mean_FC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Sex)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Sex)) +
    xlab("Sex") + 
    ylab("Mean Fibre Cross-section (FC)")
#for FDC
ggplot(covariates_data, aes(x = Sex, y = Mean_FDC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Sex)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Sex)) +
    xlab("Sex") + 
    ylab("Mean Fibre Density Cross-section (FDC)")
     
    
#check if there's a difference between the white matter integrity metrics (FBA) between gender
#for FD
t.test(covariates_data$Mean_FD ~ covariates_data$Sex, var.equal = TRUE)
#for FC
t.test(covariates_data$Mean_FC ~ covariates_data$Sex, var.equal = TRUE)
#for FDC
t.test(covariates_data$Mean_FDC ~ covariates_data$Sex, var.equal = TRUE)



#plot group status and gender against white matter integrity metrics (FBA) - 2x5 factorial ANOVA
#for FD
ggplot(covariates_data, aes(x = Group, y = Mean_FD, color = Sex)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Sex)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(.5), shape = 19, size = 2, aes(colour = Sex)) +
    xlab("Group") + 
    ylab("Mean Fibre Density (FD)")+
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) 
#for FC
ggplot(covariates_data, aes(x = Group, y = Mean_FC, color = Sex)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Sex)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(.5), shape = 19, size = 2, aes(colour = Sex)) +
    xlab("Group") + 
    ylab("Mean Fibre Cross-section (FC)")+
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) 
#for FDC
ggplot(covariates_data, aes(x = Group, y = Mean_FDC, color = Sex)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Sex)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(.5), shape = 19, size = 2, aes(colour = Sex)) +
    xlab("Group") + 
    ylab("Mean Fibre Density Cross-section (FDC)")+
    scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) 



#check for significant difference in white matter integrity between groups, gender, and groups x gender
# 2x5 factorial ANOVA
#for FD
Sex_Group_FD_mod <- lm(Mean_FD ~ Group + Sex + Group:Sex, data = covariates_data)
anova(Sex_Group_FD_mod)
#for FC
Sex_Group_FC_mod <- lm(Mean_FC ~ Group + Sex + Group:Sex, data = covariates_data)
anova(Sex_Group_FC_mod)
#for FDC
Sex_Group_FDC_mod <- lm(Mean_FDC ~ Group + Sex + Group:Sex, data = covariates_data)
anova(Sex_Group_FDC_mod)





#plot clinical site against white matter integrity metrics (FBA)
#for FD
ggplot(covariates_data, aes(x = Clinical_site, y = Mean_FD)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Clinical_site)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Clinical_site)) +
    xlab("Clinical Site") + 
    ylab("Mean Fibre Density (FD)")
#for FC
ggplot(covariates_data, aes(x = Clinical_site, y = Mean_FC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Clinical_site)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Clinical_site)) +
    xlab("Clinical Site") + 
    ylab("Mean Fibre Cross-section")
#for FDC
ggplot(covariates_data, aes(x = Clinical_site, y = Mean_FDC)) + 
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 1, aes(colour = Clinical_site)) +
    stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Clinical_site)) +
    xlab("Clinical Site") + 
    ylab("Mean Fibre Density Cross-section (FDC)")


#check for significant difference in white matter integrity between clinical sites
#for FD
Clinical_site_FD_mod <- lm(Mean_FD ~ Clinical_site, data = covariates_data)
anova(Clinical_site_FD_mod)
#for FC
Clinical_site_FC_mod <- lm(Mean_FC ~ Clinical_site, data = covariates_data)
anova(Clinical_site_FC_mod)
#for FDC
Clinical_site_FDC_mod <- lm(Mean_FDC ~ Clinical_site, data = covariates_data)
anova(Clinical_site_FDC_mod)














