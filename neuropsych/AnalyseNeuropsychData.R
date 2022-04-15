#This script will analyse the DPRC neuropsychological assessment data. Will be 
#looking at executive function / cognitive control data. Statistical tests that 
#have been run are an exploratory factor analysis (EFA), multivariate analysis 
#of variance (MANOVAs), and analysis of variance (ANOVAs).

#To deal with unbalanced ANOVA designs, you can use the anova function from the car package. (https://rpubs.com/mhanauer/300976)

#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 20/06/21

#load libraries via pacman
pacman::p_load(dplyr, ggplot2, psych, car, multcomp, lsr, BayesFactor, tidyr, GPArotation, corrplot, lsmeans, TukeyC, lme4, lmerTest, emmeans, effectsize, nlme, rstatix, sjstats, EMAtools, phia, compute.es)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#set up pathway
#setwd('/yourpathway/')
#file.choose function
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')



#######----------------- Cross-sectional (F0) analysis -----------------########

#read in csv files (participant file)
DPRC_neuropsych_data <- read.csv("cross-sectional_DPRC_neuropsych_data_lined_up_valid_participants.csv")

#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'

#convert variables
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Sex <- as.factor(DPRC_neuropsych_data$Sex)
DPRC_neuropsych_data$Sex_binary <- as.factor(DPRC_neuropsych_data$Sex_binary)

#add in trend Group variable
Trend_Group <- as.numeric(DPRC_neuropsych_data$Group)
#Trend group with contrasts that sum to zero
Trend_Group_equate_zero_contrast<- vector(mode='numeric',length=length(Trend_Group))
for (i in seq(Trend_Group)) {
  if ((Trend_Group[i] >= 1) && (Trend_Group[i]<= 1)) {
    Trend_Group_equate_zero_contrast[i] <- 2
  }  else if   ((Trend_Group[i] >= 2) && (Trend_Group[i]<= 2)) {
    Trend_Group_equate_zero_contrast[i] <- 1
  }  else if   ((Trend_Group[i] >= 3) && (Trend_Group[i]<= 3)) {
    Trend_Group_equate_zero_contrast[i] <- 0
  }  else if   ((Trend_Group[i] >= 4) && (Trend_Group[i]<= 4)) {
    Trend_Group_equate_zero_contrast[i] <- -1
  }  else if   ((Trend_Group[i] >= 5) && (Trend_Group[i]<= 5)) {
    Trend_Group_equate_zero_contrast[i] <- -2
  }
}

#----------------plot data to visualise & run stat tests ----------------------#


####--------------------------------EFA-------------------------------------####
#view a correlation matrix of the variables
#for raw data
#need to remove the columns without continuous values - should only include the 
#neuropsych test values for evaluation. 
DPRC_neuropsych_data_raw_vars <- subset(DPRC_neuropsych_data, select = c(TrailsA.Raw,
                                                                         TrailsB.Raw,
                                                                         ColorNaming.Raw, 
                                                                         WordReading.Raw, 
                                                                         Inhibition.Raw, 
                                                                         LetFluency.Raw, 
                                                                         CatFluency.Raw, 
                                                                         Switching.Raw,
                                                                         HayBTime1.Raw,
                                                                         HayBTime2.Raw,
                                                                         HayBCatA.Raw,
                                                                         HayBCatB.Raw))
#need to omit NA values from the dataframe
DPRC_neuropsych_data_raw_vars_noNAs<- na.omit(DPRC_neuropsych_data_raw_vars)

exec_raw_cor_matrix <- cor(DPRC_neuropsych_data_raw_vars_noNAs)
exec_raw_cor_matrix

#View the correlation matrix as a 'heatmap': 
corrplot(exec_raw_cor_matrix, type = "lower", method = "color", tl.col = "black")

#Run the **Kaiser Meyer Olkin (KMO) test**  to see if our data seems suitable 
#for factor analysis.The statistic shows the proportion of variance among the 
#variables that is "common variance" (might be explained by factors).
KMO(exec_raw_cor_matrix)

#Now, let's view a **scree plot** to see how much variance in the data is 
#explained by each component. This is the **eigen value** of each component.

# Calculate the eigen values
eigen_exec_raw <- eigen(exec_raw_cor_matrix)

# Create the scree plot
plot(eigen_exec_raw$values, xlab = "Component", ylab = "Eigen Values", type = "b")

# Check out variances as percentages
exec_raw_var_prop <- eigen_exec_raw$values/sum(eigen_exec_raw$values)
exec_raw_var_prop

#Run an EFA using the fa() function, extracting just 1 factor.
efa1_raw <- fa(exec_raw_cor_matrix, nfactors = 1, rotate = "oblimin")
efa1_raw
#path diagram
fa.diagram(efa1_raw, digits = 2)

#Run an EFA using the fa() function, extracting with 2 factors.
efa2_raw <- fa(exec_raw_cor_matrix, nfactors = 2, rotate = "oblimin")
efa2_raw
#path diagram
fa.diagram(efa2_raw, digits = 2)

#Run an EFA using the fa() function, extracting with 3 factors.
efa3_raw <- fa(exec_raw_cor_matrix, nfactors = 3, rotate = "oblimin")
efa3_raw
#path diagram
fa.diagram(efa3_raw, digits = 2)

#Simplify the output
print(efa3_raw$loadings, cutoff = 0.3)

#calculate factor scores from efa to run further analysis on
fscores_exec_raw <- factor.scores(DPRC_neuropsych_data_raw_vars_noNAs, efa3_raw)

#run ANOVA for the factor scores
exec_raw_data <- subset(DPRC_neuropsych_data, select = c(ParticipantID,
                                                         Age,
                                                         Classification,
                                                         Group,
                                                         Sex,
                                                         Sex_binary,
                                                         Clinical_site,
                                                         TrailsA.Raw,
                                                         TrailsB.Raw,
                                                         ColorNaming.Raw, 
                                                         WordReading.Raw, 
                                                         Inhibition.Raw, 
                                                         LetFluency.Raw, 
                                                         CatFluency.Raw, 
                                                         Switching.Raw,
                                                         HayBTime1.Raw,
                                                         HayBTime2.Raw,
                                                         HayBCatA.Raw,
                                                         HayBCatB.Raw))

#need to omit NA values from the dataframe
exec_raw_data_noNAs<- na.omit(exec_raw_data)
exec_raw_fs_data_noNAs <- cbind(exec_raw_data_noNAs, fscores_exec_raw$scores)
#run ANOVA with the factor 1
exec_raw_fs1_mod <- lm(MR1 ~ Group, data = exec_raw_fs_data_noNAs)
anova(exec_raw_fs1_mod)
#run ANOVA with the factor 2
exec_raw_fs2_mod <- lm(MR2 ~ Group, data = exec_raw_fs_data_noNAs)
anova(exec_raw_fs2_mod)
#run ANOVA with the factor 3
exec_raw_fs3_mod <- lm(MR3 ~ Group, data = exec_raw_fs_data_noNAs)
anova(exec_raw_fs3_mod)

#plot exec raw factor score 1
ggplot(subset(exec_raw_fs_data_noNAs, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = MR1, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = MR1, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("MR1 (factor score 1)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#plot exec raw factor score 2
ggplot(subset(exec_raw_fs_data_noNAs, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = MR2, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = MR2, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("MR2 (factor score 2)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#plot exec raw factor score 3
ggplot(subset(exec_raw_fs_data_noNAs, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = MR3, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = MR3, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("MR3 (factor score 3)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()


#for z-scores data
#need to remove the columns without continuous values - should only include the 
#neuropsych test values for evaluation. 
DPRC_neuropsych_data_zscores_vars <- subset(DPRC_neuropsych_data, select = c(TrailsA.Z,
                                                                         TrailsB.Z,
                                                                         ColorNaming.Z, 
                                                                         WordReading.Z, 
                                                                         Inhibition.Z, 
                                                                         LetFluency.Z, 
                                                                         CatFluency.Z, 
                                                                         Switching.z,
                                                                         HayBTime1.z,
                                                                         HayBTime2.z,
                                                                         HayBCatA.z,
                                                                         HayBCatB.z))

#need to omit NA values from the dataframe
DPRC_neuropsych_data_zscores_vars_noNAs<- na.omit(DPRC_neuropsych_data_zscores_vars)

exec_zscores_cor_matrix <- cor(DPRC_neuropsych_data_zscores_vars_noNAs)
exec_zscores_cor_matrix

#View the correlation matrix as a 'heatmap': 
corrplot(exec_zscores_cor_matrix, type = "lower", method = "color", tl.col = "black")

#Run the **Kaiser Meyer Olkin (KMO) test**  to see if our data seems suitable 
#for factor analysis.The statistic shows the proportion of variance among the 
#variables that is "common variance" (might be explained by factors).
KMO(exec_zscores_cor_matrix)

#Now, let's view a **scree plot** to see how much variance in the data is 
#explained by each component. This is the **eigen value** of each component.

# Calculate the eigen values
eigen_exec_zscores<- eigen(exec_zscores_cor_matrix)

# Create the scree plot
plot(eigen_exec_zscores$values, xlab = "Component", ylab = "Eigen Values", type = "b")

# Check out variances as percentages
exec_zscores_var_prop <- eigen_exec_zscores$values/sum(eigen_exec_zscores$values)
exec_zscores_var_prop

#Run an EFA using the fa() function, extracting just 1 factor.
efa1_z <- fa(exec_zscores_cor_matrix, nfactors = 1, rotate = "oblimin")
efa1_z
#path diagram
fa.diagram(efa1, digits = 2)

#Run an EFA using the fa() function, extracting with 2 factors.
efa2_z <- fa(exec_zscores_cor_matrix, nfactors = 2, rotate = "oblimin")
efa2_z
#path diagram
fa.diagram(efa2_z, digits = 2)

#Run an EFA using the fa() function, extracting with 3 factors.
efa3_z <- fa(exec_zscores_cor_matrix, nfactors = 3, rotate = "oblimin")
efa3_z
#path diagram
fa.diagram(efa3_z, digits = 2)

#Run an EFA using the fa() function, extracting with 4 factors.
efa4_z <- fa(exec_zscores_cor_matrix, nfactors = 4, rotate = "oblimin")
efa4_z
#path diagram
fa.diagram(efa4_z, digits = 2)

#Run an EFA using the fa() function, extracting with 5 factors.
efa5_z <- fa(exec_zscores_cor_matrix, nfactors = 5, rotate = "oblimin")
efa5_z
#path diagram
fa.diagram(efa5_z, digits = 2)

#Simplify the output
print(efa5_z$loadings, cutoff = 0.3)

#calculate factor scores from efa to run further analysis on
fscores_exec_zscores <- factor.scores(DPRC_neuropsych_data_zscores_vars_noNAs, efa5_z)

#run ANOVA for the factor scores
exec_zscores_data <- subset(DPRC_neuropsych_data, select = c(ParticipantID,
                                                         Age,
                                                         Classification,
                                                         Group,
                                                         Sex,
                                                         Sex_binary,
                                                         Clinical_site,
                                                         TrailsA.Z,
                                                         TrailsB.Z,
                                                         ColorNaming.Z, 
                                                         WordReading.Z, 
                                                         Inhibition.Z, 
                                                         LetFluency.Z, 
                                                         CatFluency.Z, 
                                                         Switching.z,
                                                         HayBTime1.z,
                                                         HayBTime2.z,
                                                         HayBCatA.z,
                                                         HayBCatB.z))

#need to omit NA values from the dataframe
exec_zscores_data_noNAs<- na.omit(exec_zscores_data)
exec_zscores_fs_data_noNAs <- cbind(exec_zscores_data_noNAs, fscores_exec_zscores$scores)
#run ANOVA with the factor 1
exec_zscores_fs1_mod <- lm(MR1 ~ Group, data = exec_zscores_fs_data_noNAs)
anova(exec_zscores_fs1_mod)
#run ANOVA with the factor 2
exec_zscores_fs2_mod <- lm(MR2 ~ Group, data = exec_zscores_fs_data_noNAs)
anova(exec_zscores_fs2_mod)
#run ANOVA with the factor 3
exec_zscores_fs3_mod <- lm(MR3 ~ Group, data = exec_zscores_fs_data_noNAs)
anova(exec_zscores_fs3_mod)
#run ANOVA with the factor 4
exec_zscores_fs4_mod <- lm(MR4 ~ Group, data = exec_zscores_fs_data_noNAs)
anova(exec_zscores_fs4_mod)
#run ANOVA with the factor 5
exec_zscores_fs5_mod <- lm(MR5 ~ Group, data = exec_zscores_fs_data_noNAs)
anova(exec_zscores_fs5_mod)

#plot exec zscores factor score 1
ggplot(subset(exec_zscores_fs_data_noNAs, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = MR1, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = MR1, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("MR1 (factor score 1)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#plot exec zscores factor score 2
ggplot(subset(exec_zscores_fs_data_noNAs, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = MR2, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = MR2, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("MR2 (factor score 2)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#plot exec zscores factor score 3
ggplot(subset(exec_zscores_fs_data_noNAs, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = MR3, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = MR3, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("MR3 (factor score 3)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#plot exec zscores factor score 4
ggplot(subset(exec_zscores_fs_data_noNAs, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = MR4, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = MR4, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("MR4 (factor score 4)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#plot exec zscores factor score 5
ggplot(subset(exec_zscores_fs_data_noNAs, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = MR5, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = MR5, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("MR5 (factor score 5)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()





####-----------------------------MANOVAs------------------------------------####
#for raw data
manova_exec_func_raw_mod <- manova(cbind(TrailsA.Raw, 
                                     TrailsB.Raw, 
                                     ColorNaming.Raw, 
                                     WordReading.Raw, 
                                     Inhibition.Raw, 
                                     LetFluency.Raw, 
                                     CatFluency.Raw, 
                                     Switching.Raw, 
                                     HayBTime1.Raw, 
                                     HayBTime2.Raw, 
                                     HayBCatA.Raw, 
                                     HayBCatB.Raw) ~ Group, data = DPRC_neuropsych_data) 
#overall model
summary(manova_exec_func_raw_mod)
#anova outputs
summary.aov(manova_exec_func_raw_mod)

#for z-scores
manova_exec_func_z_mod <- manova(cbind(TrailsA.Z, 
                                         TrailsB.Z, 
                                         ColorNaming.Z, 
                                         WordReading.Z, 
                                         Inhibition.Z, 
                                         LetFluency.Z, 
                                         CatFluency.Z, 
                                         Switching.z, 
                                         HayBTime1.z, 
                                         HayBTime2.z, 
                                         HayBCatA.z, 
                                         HayBCatB.z) ~ Group, data = DPRC_neuropsych_data) 
#overall model
summary(manova_exec_func_z_mod)
#anova outputs
summary.aov(manova_exec_func_z_mod)

#plot the manova data using the z-scores of the variables, so that it is all on the same scale. 
#put exec func variable z-scores onto a new dataset, as long format
exec_func_zscores_data <- dplyr::select(DPRC_neuropsych_data, 
                                        ParticipantID,
                                        Group,
                                        TrailsA.Z, 
                                        TrailsB.Z, 
                                        ColorNaming.Z, 
                                        WordReading.Z, 
                                        Inhibition.Z, 
                                        LetFluency.Z, 
                                        CatFluency.Z,
                                        Switching.z,
                                        HayBTime1.z,
                                        HayBTime2.z,
                                        HayBCatA.z,
                                        HayBCatB.z)
#put into long format
exec_func_zscores_data_long <- gather(exec_func_zscores_data, 
                                      "Test",
                                      "Z_scores", 
                                      TrailsA.Z, 
                                      TrailsB.Z, 
                                      ColorNaming.Z, 
                                      WordReading.Z, 
                                      Inhibition.Z, 
                                      LetFluency.Z, 
                                      CatFluency.Z, 
                                      Switching.z, 
                                      HayBTime1.z, 
                                      HayBTime2.z, 
                                      HayBCatA.z, 
                                      HayBCatB.z)

#plot data - z-score for overall 
ggplot(exec_func_zscores_data_long, aes(x = Group, y = Z_scores, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Z_scores, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Executive Functioning (Z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()


#for raw data - processing speeds variables
manova_proc_speed_raw_mod <- manova(cbind(TrailsA.Raw, 
                                         ColorNaming.Raw, 
                                         WordReading.Raw,
                                         HayBTime1.Raw) ~ Group, data = DPRC_neuropsych_data) 

#overall model
summary(manova_proc_speed_raw_mod)
#anova outputs
summary.aov(manova_proc_speed_raw_mod)

#for z-scores
manova_proc_speed_z_mod <- manova(cbind(TrailsA.Z, 
                                       ColorNaming.Z, 
                                       WordReading.Z, 
                                       HayBTime1.z) ~ Group, data = DPRC_neuropsych_data) 


#overall model
summary(manova_proc_speed_z_mod)
#anova outputs
summary.aov(manova_proc_speed_z_mod)




#plot the manova data using the z-scores of the variables, so that it is all on the same scale. 
#put proc speed variable z-scores onto a new dataset, as long format
proc_speed_zscores_data <- dplyr::select(DPRC_neuropsych_data, 
                                        ParticipantID,
                                        Group,
                                        TrailsA.Z, 
                                        ColorNaming.Z, 
                                        WordReading.Z, 
                                        HayBTime1.z)
#put into long format
proc_speed_zscores_data_long <- gather(proc_speed_zscores_data, 
                                      "Processing_Speeds",
                                      "Z_scores", 
                                      TrailsA.Z, 
                                      ColorNaming.Z, 
                                      WordReading.Z, 
                                      HayBTime1.z)


#plot data - z-score for processing speeds 
ggplot(proc_speed_zscores_data_long, aes(x = Group, y = Z_scores, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Z_scores, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Processing Speeds (Z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

#get descriptives of the collated z-scores - e.g. 'processing speeds score' mean beta values
proc_speed_z_descrip <- describeBy((proc_speed_zscores_data$TrailsA.Z + 
                                      proc_speed_zscores_data$ColorNaming.Z +
                                      proc_speed_zscores_data$WordReading.Z +
                                      proc_speed_zscores_data$HayBTime1.z) / (length(proc_speed_zscores_data) - 2), proc_speed_zscores_data$Group)

#find mean & SD from total sample: 
all_proc_speed_zscores_data <- proc_speed_zscores_data[,3:6]
noNAsproc_speed_zscores_data <- na.omit(all_proc_speed_zscores_data)
mean(as.matrix(noNAsproc_speed_zscores_data))
sd(as.matrix(noNAsproc_speed_zscores_data))

#for inhibition variables
manova_inhibition_raw_mod <- manova(cbind(TrailsB.Raw, 
                                         Inhibition.Raw, 
                                         Switching.Raw, 
                                         HayBTime2.Raw, 
                                         HayBCatA.Raw, 
                                         HayBCatB.Raw) ~ Group, data = DPRC_neuropsych_data) 
#overall model
summary(manova_inhibition_raw_mod)
#anova outputs
summary.aov(manova_inhibition_raw_mod)

#for z-scores
manova_inhibition_z_mod <- manova(cbind(TrailsB.Z, 
                                       Inhibition.Z, 
                                       Switching.z, 
                                       HayBTime2.z, 
                                       HayBCatA.z, 
                                       HayBCatB.z) ~ Group, data = DPRC_neuropsych_data) 
#overall model
summary(manova_inhibition_z_mod)
#anova outputs
summary.aov(manova_inhibition_z_mod)

#plot the manova data using the z-scores of the variables, so that it is all on the same scale. 
#put exec func variable z-scores onto a new dataset, as long format
inhibition_zscores_data <- dplyr::select(DPRC_neuropsych_data, 
                                        ParticipantID,
                                        Group,
                                        TrailsB.Z, 
                                        Inhibition.Z, 
                                        Switching.z,
                                        HayBTime2.z,
                                        HayBCatA.z,
                                        HayBCatB.z)
#put into long format
inhibition_zscores_data_long <- gather(inhibition_zscores_data, 
                                      "Inhibition",
                                      "Z_scores", 
                                      TrailsB.Z, 
                                      Inhibition.Z, 
                                      Switching.z, 
                                      HayBTime2.z, 
                                      HayBCatA.z, 
                                      HayBCatB.z)

#plot data - z-score for overall 
ggplot(inhibition_zscores_data_long, aes(x = Group, y = Z_scores, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Z_scores, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Inhibition (Z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

#get descriptives of the collated z-scores - e.g. 'executive function score' mean beta values
inhibition_z_descrip <- describeBy((inhibition_zscores_data$TrailsB.Z + 
                                    inhibition_zscores_data$Inhibition.Z +
                                    inhibition_zscores_data$Switching.z +
                                    inhibition_zscores_data$HayBTime2.z +
                                    inhibition_zscores_data$HayBCatA.z +
                                    inhibition_zscores_data$HayBCatB.z) / (length(inhibition_zscores_data) - 2), inhibition_zscores_data$Group)
#find mean & SD from total sample: 
all_inhibition_zscores_data <- inhibition_zscores_data[,3:8]
noNAsinhibition_zscores_data <- na.omit(all_inhibition_zscores_data)
mean(as.matrix(noNAsinhibition_zscores_data))
sd(as.matrix(noNAsinhibition_zscores_data))


####---------------------------------ANOVAs---------------------------------####
#run ANOVAs on your data:
#1.Test of premorbid functioning (TOPF) ---------------------------------------#
#plot TOPF.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TOPF.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TOPF.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("TOPF Score (Raw)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for TOPF.Raw
TOPF.Raw_mod <- lm(TOPF.Raw ~ Group, data = DPRC_neuropsych_data)
anova(TOPF.Raw_mod)
#check descriptive statistics per each group
TOPF.Raw_descrip <- describeBy(DPRC_neuropsych_data$TOPF.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_TOPFRaw <- DPRC_neuropsych_data$TOPF.Raw
noNAsTOPFRaw <- na.omit(all_TOPFRaw)
mean(noNAsTOPFRaw)
sd(noNAsTOPFRaw)

#plot TOPF.Z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TOPF.Z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TOPF.Z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("TOPF Score (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for TOPF.Z
TOPF.Z_mod <- lm(TOPF.Z ~ Group, data = DPRC_neuropsych_data)
anova(TOPF.Z_mod)

#2.Hayling Sentence completion test -------------------------------------------#
#plot HayBTime1.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBTime1.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBTime1.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 1 Time Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for HayBTime1.Raw
HayBTime1.Raw_mod <- lm(HayBTime1.Raw ~ Group, data = DPRC_neuropsych_data)
anova(HayBTime1.Raw_mod)
#check descriptive statistics per each group
HayBTime1.Raw_descrip <- describeBy(DPRC_neuropsych_data$HayBTime1.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_HayBTime1Raw <- DPRC_neuropsych_data$HayBTime1.Raw
noNAsHayBTime1Raw <- na.omit(all_HayBTime1Raw)
mean(noNAsHayBTime1Raw)
sd(noNAsHayBTime1Raw)
#add in Linear Trend Analysis in Linear Regression
HayBTime1Raw_LinTrend_mod <- lm(HayBTime1.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(HayBTime1Raw_LinTrend_mod)
summary(HayBTime1Raw_LinTrend_mod) 

#plot HayBTime1.z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBTime1.z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBTime1.z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 1 Time (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
HayBTime1.z_descrip <- describeBy(DPRC_neuropsych_data$HayBTime1.z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_HayBTime1z <- DPRC_neuropsych_data$HayBTime1.z
noNAsHayBTime1z <- na.omit(all_HayBTime1z)
mean(noNAsHayBTime1z)
sd(noNAsHayBTime1z)
#run ANOVA for HayBTime1.z
HayBTime1.Z_mod <- lm(HayBTime1.z ~ Group, data = DPRC_neuropsych_data)
anova(HayBTime1.Z_mod)

#plot HayBTime2.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBTime2.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBTime2.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 2 Time Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for HayBTime2.Raw
HayBTime2.Raw_mod <- lm(HayBTime2.Raw ~ Group, data = DPRC_neuropsych_data)
anova(HayBTime2.Raw_mod)
#check descriptive statistics per each group
HayBTime2.Raw_descrip <- describeBy(DPRC_neuropsych_data$HayBTime2.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_HayBTime2Raw <- DPRC_neuropsych_data$HayBTime2.Raw
noNAsHayBTime2Raw <- na.omit(all_HayBTime2Raw)
mean(noNAsHayBTime2Raw)
sd(noNAsHayBTime2Raw)
#add in Linear Trend Analysis in Linear Regression
HayBTime2Raw_LinTrend_mod <- lm(HayBTime2.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(HayBTime2Raw_LinTrend_mod)
summary(HayBTime2Raw_LinTrend_mod) 

#plot HayBTime2.z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBTime2.z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBTime2.z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 2 Time (z-score)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
HayBTime2.z_descrip <- describeBy(DPRC_neuropsych_data$HayBTime2.z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_HayBTime2z <- DPRC_neuropsych_data$HayBTime2.z
noNAsHayBTime2z <- na.omit(all_HayBTime2z)
mean(noNAsHayBTime2z)
sd(noNAsHayBTime2z)
#run ANOVA for HayBTime2.Z
HayBTime2.Z_mod <- lm(HayBTime2.z ~ Group, data = DPRC_neuropsych_data)
anova(HayBTime2.Z_mod)

#plot HayBCatA.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBCatA.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBCatA.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 2 Error Rate - Category A (total number of incorrect 'plausible' responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for HayBCatA.Raw
HayBCatA.Raw_mod <- lm(HayBCatA.Raw ~ Group, data = DPRC_neuropsych_data)
anova(HayBCatA.Raw_mod)
#effect size omnibus ANOVA
etaSquared(HayBCatA.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBCatA.Raw_mod <- glht(HayBCatA.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBCatA.Raw_mod)
confint(post_hoc_HayBCatA.Raw_mod)
#check descriptive statistics per each group
HayBCatA.Raw_descrip <- describeBy(DPRC_neuropsych_data$HayBCatA.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_HayBCatARaw <- DPRC_neuropsych_data$HayBCatA.Raw
noNAsHayBCatARaw <- na.omit(all_HayBCatARaw)
mean(noNAsHayBCatARaw)
sd(noNAsHayBCatARaw)
#effect size for sig. post hoc tests
    #for C vs. AD 
    DPRC_neuropsych_data_CvAD_CatARaw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
    DPRC_neuropsych_data_CvAD_CatARaw$Group <- droplevels(DPRC_neuropsych_data_CvAD_CatARaw$Group)
    cohensD(HayBCatA.Raw ~ Group, data = DPRC_neuropsych_data_CvAD_CatARaw) #this looks like Hedges' g? 
    #for SCD vs. AD 
    DPRC_neuropsych_data_SCDvAD_CatARaw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
    DPRC_neuropsych_data_SCDvAD_CatARaw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatARaw$Group)
    cohensD(HayBCatA.Raw ~ Group, data = DPRC_neuropsych_data_SCDvAD_CatARaw) #this looks like Hedges' g? 
    #for aMCI vs. AD 
    DPRC_neuropsych_data_aMCIvAD_CatARaw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
    DPRC_neuropsych_data_aMCIvAD_CatARaw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_CatARaw$Group)
    cohensD(HayBCatA.Raw ~ Group, data = DPRC_neuropsych_data_aMCIvAD_CatARaw) #this looks like Hedges' g? 
    #for mMCI vs. AD 
    DPRC_neuropsych_data_mMCIvAD_CatARaw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 4 | DPRC_neuropsych_data$Group == 5)
    DPRC_neuropsych_data_mMCIvAD_CatARaw$Group <- droplevels(DPRC_neuropsych_data_mMCIvAD_CatARaw$Group)
    cohensD(HayBCatA.Raw ~ Group, data = DPRC_neuropsych_data_mMCIvAD_CatARaw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_HayBCatA.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, HayBCatA.Raw)
For_Bay_data_noNas_neuropsych_HayBCatA.Raw <- na.omit(For_Bay_data_noNas_neuropsych_HayBCatA.Raw)
anovaBF(HayBCatA.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_HayBCatA.Raw) 
#add in Linear Trend Analysis in Linear Regression
HayBCatARaw_LinTrend_mod <- lm(HayBCatA.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(HayBCatARaw_LinTrend_mod)
summary(HayBCatARaw_LinTrend_mod) 

#plot HayBCatA.z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = HayBCatA.z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBCatA.z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 2 Error Rate - Category A (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
HayBCatA.z_descrip <- describeBy(DPRC_neuropsych_data$HayBCatA.z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_HayBCatAz <- DPRC_neuropsych_data$HayBCatA.z
noNAsHayBCatAz <- na.omit(all_HayBCatAz)
mean(noNAsHayBCatAz)
sd(noNAsHayBCatAz)
#run ANOVA for HayBCatA.z
HayBCatA.z_mod <- lm(HayBCatA.z ~ Group, data = DPRC_neuropsych_data)
anova(HayBCatA.z_mod)
#effect size omnibus ANOVA
etaSquared(HayBCatA.z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBCatA.z_mod <- glht(HayBCatA.z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBCatA.z_mod)
confint(post_hoc_HayBCatA.z_mod)
#check descriptive statistics per each group
HayBCatA.z_descrip <- describeBy(DPRC_neuropsych_data$HayBCatA.z, DPRC_neuropsych_data$Group)
#effect size for sig. post hoc tests
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_CatAZ <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_CatAZ $Group <- droplevels(DPRC_neuropsych_data_CvAD_CatAZ $Group)
  cohensD(HayBCatA.z ~ Group, data = DPRC_neuropsych_data_CvAD_CatAZ ) #this looks like Hedges' g? 
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_CatAZ  <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatAZ$Group <- droplevels(DPRC_neuropsych_data_SCDvAD.Z$Group)
  cohensD(HayBCatA.z ~ Group, data = DPRC_neuropsych_data_SCDvAD_CatAZ) #this looks like Hedges' g? 
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_CatAZ <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_CatAZ$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_CatAZ$Group)
  cohensD(HayBCatA.z ~ Group, data = DPRC_neuropsych_data_aMCIvAD_CatAZ) #this looks like Hedges' g? 
  #for mMCI vs. AD 
  DPRC_neuropsych_data_mMCIvAD_CatAZ<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 4 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_mMCIvAD_CatAZ$Group <- droplevels(DPRC_neuropsych_data_mMCIvAD_CatAZ$Group)
  cohensD(HayBCatA.z ~ Group, data = DPRC_neuropsych_data_mMCIvAD_CatAZ) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_HayBCatA.Z <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, HayBCatA.z)
For_Bay_data_noNas_neuropsych_HayBCatA.Z <- na.omit(For_Bay_data_noNas_neuropsych_HayBCatA.Z)
anovaBF(HayBCatA.z ~ Group, data = For_Bay_data_noNas_neuropsych_HayBCatA.Z) 

#plot HayBCatB.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBCatB.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBCatB.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 2 Error Rate - Category B (total number of incorrect 'somewhat plausible' responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for HayBCatB.Raw
HayBCatB.Raw_mod <- lm(HayBCatB.Raw ~ Group, data = DPRC_neuropsych_data)
anova(HayBCatB.Raw_mod)
#effect size omnibus ANOVA
etaSquared(HayBCatB.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBCatB.Raw_mod <- glht(HayBCatB.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBCatB.Raw_mod)
confint(post_hoc_HayBCatB.Raw_mod)
#check descriptive statistics per each group
HayBCatB.Raw_descrip <- describeBy(DPRC_neuropsych_data$HayBCatB.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_HayBCatBRaw <- DPRC_neuropsych_data$HayBCatB.Raw
noNAsHayBCatBRaw <- na.omit(all_HayBCatBRaw)
mean(noNAsHayBCatBRaw)
sd(noNAsHayBCatBRaw)
#effect size for sig. post hoc tests
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_CatBRaw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_CatBRaw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_CatBRaw$Group)
  cohensD(HayBCatB.Raw ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_CatBRaw) #this looks like Hedges' g? 
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_CatBRaw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatBRaw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatBRaw$Group)
  cohensD(HayBCatB.Raw ~ Group, data = DPRC_neuropsych_data_SCDvAD_CatBRaw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_HayBCatB.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, HayBCatB.Raw)
For_Bay_data_noNas_neuropsych_HayBCatB.Raw <- na.omit(For_Bay_data_noNas_neuropsych_HayBCatB.Raw)
anovaBF(HayBCatB.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_HayBCatB.Raw) 
#add in Linear Trend Analysis in Linear Regression
HayBCatBRaw_LinTrend_mod <- lm(HayBCatB.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(HayBCatBRaw_LinTrend_mod)
summary(HayBCatBRaw_LinTrend_mod) 

#plot HayBCatB.z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBCatB.z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBCatB.z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 2 Error Rate - Category B (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
HayBCatB.z_descrip <- describeBy(DPRC_neuropsych_data$HayBCatB.z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_HayBCatBz<- DPRC_neuropsych_data$HayBCatB.z
noNAsHayBCatBz <- na.omit(all_HayBCatBz)
mean(noNAsHayBCatBz)
sd(noNAsHayBCatBz)
#run ANOVA for HayBCatB.Z
HayBCatB.Z_mod <- lm(HayBCatB.z ~ Group, data = DPRC_neuropsych_data)
anova(HayBCatB.Z_mod)
#effect size omnibus ANOVA
etaSquared(HayBCatB.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBCatB.Z_mod <- glht(HayBCatB.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBCatB.Z_mod)
confint(post_hoc_HayBCatB.Z_mod)
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_HayBCatB.Z <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, HayBCatB.z)
For_Bay_data_noNas_neuropsych_HayBCatB.Z <- na.omit(For_Bay_data_noNas_neuropsych_HayBCatB.Z)
anovaBF(HayBCatB.z ~ Group, data = For_Bay_data_noNas_neuropsych_HayBCatB.Z) 

#3.D-KEFS Stroop Task ---------------------------------------------------------#
#plot ColorNaming.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = ColorNaming.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = ColorNaming.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Color Naming Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for ColorNaming.Raw
ColorNaming.Raw_mod <- lm(ColorNaming.Raw ~ Group, data = DPRC_neuropsych_data)
anova(ColorNaming.Raw_mod)
#effect size omnibus ANOVA
etaSquared(ColorNaming.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_ColorNaming.Raw_mod <- glht(ColorNaming.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_ColorNaming.Raw_mod)
confint(post_hoc_ColorNaming.Raw_mod)
#check descriptive statistics per each group
ColorNaming.Raw_descrip <- describeBy(DPRC_neuropsych_data$ColorNaming.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_ColorNamingRaw <- DPRC_neuropsych_data$ColorNaming.Raw
noNAsColorNamingRaw <- na.omit(all_ColorNamingRaw)
mean(noNAsColorNamingRaw)
sd(noNAsColorNamingRaw)
#effect size for sig. post hoc tests
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_ColorNaming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_ColorNaming.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_ColorNaming.Raw$Group)
  cohensD(ColorNaming.Raw ~ Group, data = DPRC_neuropsych_data_CvmMCI_ColorNaming.Raw) #this looks like Hedges' g? 
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_ColorNaming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_ColorNaming.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_ColorNaming.Raw$Group)
  cohensD(ColorNaming.Raw ~ Group, data = DPRC_neuropsych_data_CvAD_ColorNaming.Raw) #this looks like Hedges' g? 
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_ColorNaming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_ColorNaming.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_ColorNaming.Raw$Group)
  cohensD(ColorNaming.Raw ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_ColorNaming.Raw) #this looks like Hedges' g? 
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw$Group)
  cohensD(ColorNaming.Raw ~ Group, data = DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_ColorNaming.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, ColorNaming.Raw)
For_Bay_data_noNas_neuropsych_ColorNaming.Raw <- na.omit(For_Bay_data_noNas_neuropsych_ColorNaming.Raw)
anovaBF(ColorNaming.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_ColorNaming.Raw) 
#add in Linear Trend Analysis in Linear Regression
ColorNamingRaw_LinTrend_mod <- lm(ColorNaming.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(ColorNamingRaw_LinTrend_mod)
summary(ColorNamingRaw_LinTrend_mod) 

#plot ColorNaming.Z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = ColorNaming.Z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = ColorNaming.Z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Color Naming (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
ColorNaming.Z_descrip <- describeBy(DPRC_neuropsych_data$ColorNaming.Z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_ColorNamingZ<- DPRC_neuropsych_data$ColorNaming.Z
noNAsColorNamingZ <- na.omit(all_ColorNamingZ)
mean(noNAsColorNamingZ)
sd(noNAsColorNamingZ)
#run ANOVA for ColorNaming.Z
ColorNaming.Z_mod <- lm(ColorNaming.Z ~ Group, data = DPRC_neuropsych_data)
anova(ColorNaming.Z_mod)
#effect size omnibus ANOVA
etaSquared(ColorNaming.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_ColorNaming.Z_mod <- glht(ColorNaming.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_ColorNaming.Z_mod)
confint(post_hoc_ColorNaming.Z_mod)
#effect size for sig. post hoc tests
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_ColorNaming.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_ColorNaming.Z$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_ColorNaming.Z$Group)
  cohensD(ColorNaming.Z ~ Group, data = DPRC_neuropsych_data_CvmMCI_ColorNaming.Z) #this looks like Hedges' g? 
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_ColorNaming.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_ColorNaming.Z$Group <- droplevels(DPRC_neuropsych_data_CvAD_ColorNaming.Z$Group)
  cohensD(ColorNaming.Z ~ Group, data = DPRC_neuropsych_data_CvAD_ColorNaming.Z) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_ColorNaming.Z <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, ColorNaming.Z)
For_Bay_data_noNas_neuropsych_ColorNaming.Z <- na.omit(For_Bay_data_noNas_neuropsych_ColorNaming.Z)
anovaBF(ColorNaming.Z ~ Group, data = For_Bay_data_noNas_neuropsych_ColorNaming.Z) 

#plot WordReading.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = WordReading.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = WordReading.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Word Reading Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for WordReading.Raw
WordReading.Raw_mod <- lm(WordReading.Raw ~ Group, data = DPRC_neuropsych_data)
anova(WordReading.Raw_mod)
#effect size omnibus ANOVA
etaSquared(WordReading.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_WordReading.Raw_mod <- glht(WordReading.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_WordReading.Raw_mod)
confint(post_hoc_WordReading.Raw_mod)
#check descriptive statistics per each group
WordReading.Raw_descrip <- describeBy(DPRC_neuropsych_data$WordReading.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_WordReadingRaw <- DPRC_neuropsych_data$WordReading.Raw
noNAsWordReadingRaw <- na.omit(all_WordReadingRaw)
mean(noNAsWordReadingRaw)
sd(noNAsWordReadingRaw)
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_WordReading.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, WordReading.Raw)
For_Bay_data_noNas_neuropsych_WordReading.Raw <- na.omit(For_Bay_data_noNas_neuropsych_WordReading.Raw)
anovaBF(WordReading.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_WordReading.Raw) 
#add in Linear Trend Analysis in Linear Regression
WordReadingRaw_LinTrend_mod <- lm(WordReading.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(WordReadingRaw_LinTrend_mod)
summary(WordReadingRaw_LinTrend_mod) 

#plot WordReading.Z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = WordReading.Z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = WordReading.Z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Word Reading (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
WordReading.Z_descrip <- describeBy(DPRC_neuropsych_data$WordReading.Z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_WordReadingZ<- DPRC_neuropsych_data$WordReading.Z
noNAsWordReadingZ <- na.omit(all_WordReadingZ)
mean(noNAsWordReadingZ)
sd(noNAsWordReadingZ)
#run ANOVA for WordReading.Z
WordReading.Z_mod <- lm(WordReading.Z ~ Group, data = DPRC_neuropsych_data)
anova(WordReading.Z_mod)

#plot Inhibition.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Inhibition.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Inhibition.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Inhibition Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for Inhibition.Raw
Inhibition.Raw_mod <- lm(Inhibition.Raw ~ Group, data = DPRC_neuropsych_data)
anova(Inhibition.Raw_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Raw_mod <- glht(Inhibition.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Raw_mod)
confint(post_hoc_Inhibition.Raw_mod)
#check descriptive statistics per each group
Inhibition.Raw_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_InhibitionRaw <- DPRC_neuropsych_data$Inhibition.Raw
noNAsInhibitionRaw <- na.omit(all_InhibitionRaw)
mean(noNAsInhibitionRaw)
sd(noNAsInhibitionRaw)
#effect size for sig. post hoc tests
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Raw$Group)
  cohensD(Inhibition.Raw ~ Group, data = DPRC_neuropsych_data_CvAD_Inhibition.Raw) #this looks like Hedges' g? 
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_Inhibition.Raw  <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Inhibition.Raw$Group)
  cohensD(Inhibition.Raw ~ Group, data = DPRC_neuropsych_data_SCDvAD_Inhibition.Raw) #this looks like Hedges' g? 
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Inhibition.Raw$Group)
  cohensD(Inhibition.Raw ~ Group, data = DPRC_neuropsych_data_CvmMCI_Inhibition.Raw) #this looks like Hedges' g? 
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw$Group)
  cohensD(Inhibition.Raw ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw) #this looks like Hedges' g?
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw$Group)
  cohensD(Inhibition.Raw ~ Group, data = DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_Inhibition.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, Inhibition.Raw)
For_Bay_data_noNas_neuropsych_Inhibition.Raw<- na.omit(For_Bay_data_noNas_neuropsych_Inhibition.Raw)
anovaBF(Inhibition.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_Inhibition.Raw) 
#add in Linear Trend Analysis in Linear Regression
InhibitionRaw_LinTrend_mod <- lm(Inhibition.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(InhibitionRaw_LinTrend_mod)
summary(InhibitionRaw_LinTrend_mod) 

#plot Inhibition.Z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Inhibition.Z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Inhibition.Z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Inhibition (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
Inhibition.Z_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_InhibitionZ<- DPRC_neuropsych_data$Inhibition.Z
noNAsInhibitionZ <- na.omit(all_InhibitionZ)
mean(noNAsInhibitionZ)
sd(noNAsInhibitionZ)
#run ANOVA for Inhibition.Z
Inhibition.Z_mod <- lm(Inhibition.Z ~ Group, data = DPRC_neuropsych_data)
anova(Inhibition.Z_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Z_mod <- glht(Inhibition.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Z_mod)
confint(post_hoc_Inhibition.Z_mod)
#check descriptive statistics per each group
Inhibition.Z_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Z, DPRC_neuropsych_data$Group)
#effect size for sig. post hoc tests
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_Inhibition.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Z$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Z$Group)
  cohensD(Inhibition.Z ~ Group, data = DPRC_neuropsych_data_CvAD_Inhibition.Z) #this looks like Hedges' g? 
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Inhibition.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Inhibition.Z$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Inhibition.Z$Group)
  cohensD(Inhibition.Z ~ Group, data = DPRC_neuropsych_data_CvmMCI_Inhibition.Z) #this looks like Hedges' g? 
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Z<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Inhibition.Z$Group)
  cohensD(Inhibition.Z ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_Inhibition.Z) #this looks like Hedges' g?
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_Inhibition.Z <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, Inhibition.Z)
For_Bay_data_noNas_neuropsych_Inhibition.Z<- na.omit(For_Bay_data_noNas_neuropsych_Inhibition.Z)
anovaBF(Inhibition.Z ~ Group, data = For_Bay_data_noNas_neuropsych_Inhibition.Z) 

#plot Inhibition.Colour.Naming (Inhibition / Colour Naming) (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Inhibition.Colour.Naming, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Inhibition.Colour.Naming, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Inhibition / Colour Naming") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
Inhibition.Colour.Naming_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Colour.Naming, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_InhibitionColourNaming <- DPRC_neuropsych_data$Inhibition.Colour.Naming
noNAsInhibitionColourNaming <- na.omit(all_InhibitionColourNaming)
mean(noNAsInhibitionColourNaming)
sd(noNAsInhibitionColourNaming)
#run ANOVA for Inhibition.Colour.Naming (Inhibition / Colour Naming)
Inhibition.Colour.Naming_mod <- lm(Inhibition.Colour.Naming ~ Group, data = DPRC_neuropsych_data)
anova(Inhibition.Colour.Naming_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Colour.Naming_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Colour.Naming_mod <- glht(Inhibition.Colour.Naming_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Colour.Naming_mod)
confint(post_hoc_Inhibition.Colour.Naming_mod)
#effect size for sig. post hoc tests
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming$Group)
  cohensD(Inhibition.Colour.Naming ~ Group, data = DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming) #this looks like Hedges' g? 
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming  <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming$Group)
  cohensD(Inhibition.Colour.Naming ~ Group, data = DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming) #this looks like Hedges' g? 
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming$Group)
  cohensD(Inhibition.Colour.Naming ~ Group, data = DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_Inhibition.Colour.Naming <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, Inhibition.Colour.Naming)
For_Bay_data_noNas_neuropsych_Inhibition.Colour.Naming<- na.omit(For_Bay_data_noNas_neuropsych_Inhibition.Colour.Naming)
anovaBF(Inhibition.Colour.Naming ~ Group, data = For_Bay_data_noNas_neuropsych_Inhibition.Colour.Naming) 
#add in Linear Trend Analysis in Linear Regression
Inhibition.Colour.Naming_LinTrend_mod <- lm(Inhibition.Colour.Naming ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(Inhibition.Colour.Naming_LinTrend_mod)
summary(Inhibition.Colour.Naming_LinTrend_mod) 

#plot Inhibition.Colour.Naming z-scores (Inhibition / Colour Naming.zscores.calc) (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Inhibition.Colour.Naming.zscores.calc., fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Inhibition.Colour.Naming.zscores.calc., color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Inhibition / Colour Naming Z-scores") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
Inhibition.Colour.NamingZ_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Colour.Naming.zscores.calc., DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_InhibitionColourNamingZ <- DPRC_neuropsych_data$Inhibition.Colour.Naming.zscores.calc.
noNAsInhibitionColourNamingZ <- na.omit(all_InhibitionColourNamingZ)
mean(noNAsInhibitionColourNamingZ)
sd(noNAsInhibitionColourNamingZ)
#run ANOVA for Inhibition.Colour.Naming (Inhibition / Colour Naming)
Inhibition.Colour.NamingZ_mod <- lm(Inhibition.Colour.Naming.zscores.calc. ~ Group, data = DPRC_neuropsych_data)
anova(Inhibition.Colour.NamingZ_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Colour.NamingZ_mod)

#plot Inhibition.Word.Reading (Inhibition / Word Reading) (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Inhibition.Word.Reading, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Inhibition.Word.Reading, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Inhibition / Word Reading") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for Inhibition.Word.Reading (Inhibition / Word Reading)
Inhibition.Word.Reading_mod <- lm(Inhibition.Word.Reading ~ Group, data = DPRC_neuropsych_data)
anova(Inhibition.Word.Reading_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Word.Reading_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Word.Reading_mod <- glht(Inhibition.Word.Reading_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Word.Reading_mod)
confint(post_hoc_Inhibition.Word.Reading_mod)
#check descriptive statistics per each group
Inhibition.Word.Reading_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Word.Reading, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_InhibitionWordReading <- DPRC_neuropsych_data$Inhibition.Word.Reading
noNAsInhibitionWordReading <- na.omit(all_InhibitionWordReading)
mean(noNAsInhibitionWordReading)
sd(noNAsInhibitionWordReading)
#effect size for sig. post hoc tests
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading$Group)
  cohensD(Inhibition.Word.Reading ~ Group, data = DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading) #this looks like Hedges' g? 
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading$Group)
  cohensD(Inhibition.Word.Reading ~ Group, data = DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading) #this looks like Hedges' g? 
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading  <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading$Group)
  cohensD(Inhibition.Word.Reading ~ Group, data = DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading) #this looks like Hedges' g? 
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Word.Reading<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Inhibition.Word.Reading$Group)
  cohensD(Inhibition.Word.Reading ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_Inhibition.Word.Reading) #this looks like Hedges' g?
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading$Group)
  cohensD(Inhibition.Word.Reading ~ Group, data = DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_Inhibition.Word.Reading <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, Inhibition.Word.Reading)
For_Bay_data_noNas_neuropsych_Inhibition.Word.Reading <- na.omit(For_Bay_data_noNas_neuropsych_Inhibition.Word.Reading)
anovaBF(Inhibition.Word.Reading ~ Group, data = For_Bay_data_noNas_neuropsych_Inhibition.Word.Reading) 
#add in Linear Trend Analysis in Linear Regression
Inhibition.Word.Reading_LinTrend_mod <- lm(Inhibition.Word.Reading ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(Inhibition.Word.Reading_LinTrend_mod)
summary(Inhibition.Word.Reading_LinTrend_mod) 

#plot Inhibition.Word.Reading z-scores (Inhibition / Colour Naming.zscores.calc) (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Inhibition.Word.Reading.zscores.calc., fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Inhibition.Word.Reading.zscores.calc., color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Stroop - Inhibition / Word Reading Z-scores") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
Inhibition.Word.ReadingZ_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Word.Reading.zscores.calc., DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_InhibitionWordReadingZ <- DPRC_neuropsych_data$Inhibition.Word.Reading.zscores.calc.
noNAsInhibitionWordReadingZ <- na.omit(all_InhibitionWordReadingZ)
mean(noNAsInhibitionWordReadingZ)
sd(noNAsInhibitionWordReadingZ)
#run ANOVA for Inhibition.Word.Reading (Inhibition / Colour Naming)
Inhibition.Word.ReadingZ_mod <- lm(Inhibition.Word.Reading.zscores.calc. ~ Group, data = DPRC_neuropsych_data)
anova(Inhibition.Word.ReadingZ_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Word.ReadingZ_mod)


#4.D-KEFS Verbal Fluency + Category Fluency Task ------------------------------#
#plot LetFluency.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = LetFluency.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = LetFluency.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Letter Fluency (total correct responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for LetFluency.Raw
LetFluency.Raw_mod <- lm(LetFluency.Raw ~ Group, data = DPRC_neuropsych_data)
anova(LetFluency.Raw_mod)
#effect size omnibus ANOVA
etaSquared(LetFluency.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_LetFluency.Raw_mod <- glht(LetFluency.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_LetFluency.Raw_mod)
confint(post_hoc_LetFluency.Raw_mod)
#check descriptive statistics per each group
LetFluency.Raw_descrip <- describeBy(DPRC_neuropsych_data$LetFluency.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_LetFluencyRaw <- DPRC_neuropsych_data$LetFluency.Raw
noNAsLetFluencyRaw <- na.omit(all_LetFluencyRaw)
mean(noNAsLetFluencyRaw)
sd(noNAsLetFluencyRaw)
#effect size for sig. post hoc tests
  #for C vs. SCD 
  DPRC_neuropsych_data_CvSCD_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 2)
  DPRC_neuropsych_data_CvSCD_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvSCD_LetFluency.Raw$Group)
  cohensD(LetFluency.Raw ~ Group, data = DPRC_neuropsych_data_CvSCD_LetFluency.Raw) #this looks like Hedges' g? 
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group)
  cohensD(LetFluency.Raw ~ Group, data = DPRC_neuropsych_data_CvaMCI_LetFluency.Raw) #this looks like Hedges' g? 
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_LetFluency.Raw<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_LetFluency.Raw$Group)
  cohensD(LetFluency.Raw ~ Group, data = DPRC_neuropsych_data_CvmMCI_LetFluency.Raw) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_LetFluency.Raw<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_LetFluency.Raw$Group)
  cohensD(LetFluency.Raw ~ Group, data = DPRC_neuropsych_data_CvAD_LetFluency.Raw) #this looks like Hedges' g?
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_LetFluency.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, LetFluency.Raw)
For_Bay_data_noNas_neuropsych_LetFluency.Raw<- na.omit(For_Bay_data_noNas_neuropsych_LetFluency.Raw)
anovaBF(LetFluency.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_LetFluency.Raw) 
#add in Linear Trend Analysis in Linear Regression
LetFluencyRaw_LinTrend_mod <- lm(LetFluency.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(LetFluencyRaw_LinTrend_mod)
summary(LetFluencyRaw_LinTrend_mod) 

#plot LetFluency.Z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = LetFluency.Z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = LetFluency.Z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Letter Fluency (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
LetFluency.Z_descrip <- describeBy(DPRC_neuropsych_data$LetFluency.Z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_LetFluencyZ <- DPRC_neuropsych_data$LetFluency.Z
noNAsLetFluencyZ <- na.omit(all_LetFluencyZ)
mean(noNAsLetFluencyZ)
sd(noNAsLetFluencyZ)
#run ANOVA for LetFluency.Raw
LetFluency.Z_mod <- lm(LetFluency.Z ~ Group, data = DPRC_neuropsych_data)
anova(LetFluency.Z_mod)
#effect size omnibus ANOVA
etaSquared(LetFluency.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_LetFluency.Z_mod <- glht(LetFluency.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_LetFluency.Z_mod)
confint(post_hoc_LetFluency.Z_mod)
#effect size for sig. post hoc tests
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_LetFluency.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_LetFluency.Z$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_LetFluency.Z$Group)
  cohensD(LetFluency.Z ~ Group, data = DPRC_neuropsych_data_CvaMCI_LetFluency.Z) #this looks like Hedges' g? 
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_LetFluency.Z<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_LetFluency.Z$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_LetFluency.Z$Group)
  cohensD(LetFluency.Z ~ Group, data = DPRC_neuropsych_data_CvmMCI_LetFluency.Z) #this looks like Hedges' g?
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_LetFluency.Z <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, LetFluency.Z)
For_Bay_data_noNas_neuropsych_LetFluency.Z<- na.omit(For_Bay_data_noNas_neuropsych_LetFluency.Z)
anovaBF(LetFluency.Z ~ Group, data = For_Bay_data_noNas_neuropsych_LetFluency.Z) 

#plot CatFluency.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = CatFluency.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = CatFluency.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Category Fluency (total correct responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for CatFluency.Raw
CatFluency.Raw_mod <- lm(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data)
anova(CatFluency.Raw_mod)
#effect size omnibus ANOVA
etaSquared(CatFluency.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_CatFluency.Raw_mod <- glht(CatFluency.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_CatFluency.Raw_mod)
confint(post_hoc_CatFluency.Raw_mod)
#check descriptive statistics per each group
CatFluency.Raw_descrip <- describeBy(DPRC_neuropsych_data$CatFluency.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_CatFluencyRaw <- DPRC_neuropsych_data$CatFluency.Raw
noNAsCatFluencyRaw <- na.omit(all_CatFluencyRaw)
mean(noNAsCatFluencyRaw)
sd(noNAsCatFluencyRaw)
#effect size for sig. post hoc tests
  #for C vs. aMCI  
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group)
  cohensD(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data_CvaMCI_CatFluency.Raw) #this looks like Hedges' g? 
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group)
  cohensD(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data_CvmMCI_CatFluency.Raw) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group)
  cohensD(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data_CvAD_CatFluency.Raw) #this looks like Hedges' g?
  #for SCD vs. aMCI 
  DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw$Group)
  cohensD(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw) #this looks like Hedges' g?
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group)
  cohensD(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group)
  cohensD(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data_SCDvAD_CatFluency.Raw) #this looks like Hedges' g? 
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_CatFluency.Raw$Group)
  cohensD(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data_aMCIvAD_CatFluency.Raw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_CatFluency.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, CatFluency.Raw)
For_Bay_data_noNas_neuropsych_CatFluency.Raw <- na.omit(For_Bay_data_noNas_neuropsych_CatFluency.Raw)
anovaBF(CatFluency.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_CatFluency.Raw) 
#add in Linear Trend Analysis in Linear Regression
CatFluencyRaw_LinTrend_mod <- lm(CatFluency.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(CatFluencyRaw_LinTrend_mod)
summary(CatFluencyRaw_LinTrend_mod) 

#plot CatFluency.Z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = CatFluency.Z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = CatFluency.Z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Category Fluency (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
CatFluency.Z_descrip <- describeBy(DPRC_neuropsych_data$CatFluency.Z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_CatFluencyZ <- DPRC_neuropsych_data$CatFluency.Z
noNAsCatFluencyZ <- na.omit(all_CatFluencyZ)
mean(noNAsCatFluencyZ)
sd(noNAsCatFluencyZ)
#run ANOVA for CatFluency.Z
CatFluency.Z_mod <- lm(CatFluency.Z ~ Group, data = DPRC_neuropsych_data)
anova(CatFluency.Z_mod)
#effect size omnibus ANOVA
etaSquared(CatFluency.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_CatFluency.Z_mod <- glht(CatFluency.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_CatFluency.Z_mod)
confint(post_hoc_CatFluency.Z_mod)
#effect size for sig. post hoc tests
  #for C vs. aMCI  
  DPRC_neuropsych_data_CvaMCI_CatFluency.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_CatFluency.Z$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_CatFluency.Z$Group)
  cohensD(CatFluency.Z ~ Group, data = DPRC_neuropsych_data_CvaMCI_CatFluency.Z) #this looks like Hedges' g? 
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_CatFluency.Z<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_CatFluency.Z$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_CatFluency.Z$Group)
  cohensD(CatFluency.Z ~ Group, data = DPRC_neuropsych_data_CvmMCI_CatFluency.Z) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_CatFluency.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_CatFluency.Z$Group <- droplevels(DPRC_neuropsych_data_CvAD_CatFluency.Z$Group)
  cohensD(CatFluency.Z ~ Group, data = DPRC_neuropsych_data_CvAD_CatFluency.Z) #this looks like Hedges' g?
  #for SCD vs. aMCI 
  DPRC_neuropsych_data_SCDvaMCI_CatFluency.Z<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_CatFluency.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_CatFluency.Z$Group)
  cohensD(CatFluency.Z ~ Group, data = DPRC_neuropsych_data_SCDvaMCI_CatFluency.Z) #this looks like Hedges' g?
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Z$Group)
  cohensD(CatFluency.Z ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_CatFluency.Z) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_CatFluency.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatFluency.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatFluency.Z$Group)
  cohensD(CatFluency.Z ~ Group, data = DPRC_neuropsych_data_SCDvAD_CatFluency.Z) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_CatFluency.Z <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, CatFluency.Z)
For_Bay_data_noNas_neuropsych_CatFluency.Z <- na.omit(For_Bay_data_noNas_neuropsych_CatFluency.Z)
anovaBF(CatFluency.Z ~ Group, data = For_Bay_data_noNas_neuropsych_CatFluency.Z) 

#plot Switching.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Switching.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Switching.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Category Switching (total correct responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for Switching.Raw
Switching.Raw_mod <- lm(Switching.Raw ~ Group, data = DPRC_neuropsych_data)
anova(Switching.Raw_mod)
#effect size omnibus ANOVA
etaSquared(Switching.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Switching.Raw_mod <- glht(Switching.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Switching.Raw_mod)
confint(post_hoc_Switching.Raw_mod)
#check descriptive statistics per each group
Switching.Raw_descrip <- describeBy(DPRC_neuropsych_data$Switching.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_SwitchingRaw <- DPRC_neuropsych_data$Switching.Raw
noNAsSwitchingRaw <- na.omit(all_SwitchingRaw)
mean(noNAsSwitchingRaw)
sd(noNAsSwitchingRaw)
#effect size for sig. post hoc tests
  #for C vs. aMCI  
  DPRC_neuropsych_data_CvaMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_Switching.Raw$Group)
  cohensD(Switching.Raw ~ Group, data = DPRC_neuropsych_data_CvaMCI_Switching.Raw) #this looks like Hedges' g? 
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group)
  cohensD(Switching.Raw ~ Group, data = DPRC_neuropsych_data_CvmMCI_Switching.Raw) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Switching.Raw$Group)
  cohensD(Switching.Raw ~ Group, data = DPRC_neuropsych_data_CvAD_Switching.Raw) #this looks like Hedges' g?
  #for SCD vs. aMCI 
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group)
  cohensD(Switching.Raw ~ Group, data = DPRC_neuropsych_data_SCDvaMCI_Switching.Raw) #this looks like Hedges' g?
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group)
  cohensD(Switching.Raw ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_Switching.Raw) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_CatFluency.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatFluency.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatFluency.Z$Group)
  cohensD(CatFluency.Z ~ Group, data = DPRC_neuropsych_data_SCDvAD_CatFluency.Z) #this looks like Hedges' g? 
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Switching.Raw$Group)
  cohensD(Switching.Raw ~ Group, data = DPRC_neuropsych_data_aMCIvAD_Switching.Raw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_Switching.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, Switching.Raw)
For_Bay_data_noNas_neuropsych_Switching.Raw <- na.omit(For_Bay_data_noNas_neuropsych_Switching.Raw)
anovaBF(Switching.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_Switching.Raw) 
#add in Linear Trend Analysis in Linear Regression
SwitchingRaw_LinTrend_mod <- lm(Switching.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(SwitchingRaw_LinTrend_mod)
summary(SwitchingRaw_LinTrend_mod) 

#plot Switching.z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Switching.z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Switching.z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Category Switching (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
Switching.z_descrip <- describeBy(DPRC_neuropsych_data$Switching.z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_Switchingz <- DPRC_neuropsych_data$Switching.z
noNAsSwitchingz <- na.omit(all_Switchingz)
mean(noNAsSwitchingz)
sd(noNAsSwitchingz)
#run ANOVA for Switching.Z
Switching.z_mod <- lm(Switching.z ~ Group, data = DPRC_neuropsych_data)
anova(Switching.z_mod)
#effect size omnibus ANOVA
etaSquared(Switching.z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Switching.z_mod <- glht(Switching.z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Switching.z_mod)
confint(post_hoc_Switching.z_mod)
#effect size for sig. post hoc tests
  #for C vs. aMCI  
  DPRC_neuropsych_data_CvaMCI_Switching.z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_Switching.z$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_Switching.z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_CvaMCI_Switching.z) #this looks like Hedges' g? 
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_Switching.z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Switching.z$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Switching.z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_CvmMCI_Switching.z) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_Switching.z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Switching.z$Group <- droplevels(DPRC_neuropsych_data_CvAD_Switching.z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_CvAD_Switching.z) #this looks like Hedges' g?
  #for SCD vs. aMCI 
  DPRC_neuropsych_data_SCDvaMCI_Switching.z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_Switching.z$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_Switching.z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_SCDvaMCI_Switching.z) #this looks like Hedges' g?
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_Switching.z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Switching.z$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Switching.z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_Switching.z) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_Switching.z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Switching.z$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Switching.z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_SCDvAD_Switching.z) #this looks like Hedges' g? 
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_Switching.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Switching.Z$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Switching.Z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_aMCIvAD_Switching.Z) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_Switching.z <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, Switching.z)
For_Bay_data_noNas_neuropsych_Switching.z <- na.omit(For_Bay_data_noNas_neuropsych_Switching.z)
anovaBF(Switching.z ~ Group, data = For_Bay_data_noNas_neuropsych_Switching.z) 

#5.Trail Making Test (TMT) A & B ----------------------------------------------#
#plot TrailsA.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TrailsA.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TrailsA.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Trail Making Test A Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for TrailsA.Raw
TrailsA.Raw_mod <- lm(TrailsA.Raw~ Group, data = DPRC_neuropsych_data)
anova(TrailsA.Raw_mod)
#effect size omnibus ANOVA
etaSquared(TrailsA.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsA.Raw_mod <- glht(TrailsA.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsA.Raw_mod)
confint(post_hoc_TrailsA.Raw_mod)
#check descriptive statistics per each group
TrailsA.Raw_descrip <- describeBy(DPRC_neuropsych_data$TrailsA.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_TrailsARaw <- DPRC_neuropsych_data$TrailsA.Raw
noNAsTrailsARaw <- na.omit(all_TrailsARaw)
mean(noNAsTrailsARaw)
sd(noNAsTrailsARaw)
#effect size for sig. post hoc tests
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TrailsA.Raw$Group)
  cohensD(TrailsA.Raw ~ Group, data = DPRC_neuropsych_data_CvmMCI_TrailsA.Raw) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_TrailsA.Raw$Group)
  cohensD(TrailsA.Raw ~ Group, data = DPRC_neuropsych_data_CvAD_TrailsA.Raw) #this looks like Hedges' g?
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TrailsA.Raw$Group)
  cohensD(TrailsA.Raw ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_TrailsA.Raw) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TrailsA.Raw$Group)
  cohensD(TrailsA.Raw ~ Group, data = DPRC_neuropsych_data_SCDvAD_TrailsA.Raw) #this looks like Hedges' g? 
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw$Group)
  cohensD(TrailsA.Raw ~ Group, data = DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_TrailsA.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, TrailsA.Raw)
For_Bay_data_noNas_neuropsych_TrailsA.Raw <- na.omit(For_Bay_data_noNas_neuropsych_TrailsA.Raw)
anovaBF(TrailsA.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_TrailsA.Raw) 
#add in Linear Trend Analysis in Linear Regression
TrailsARaw_LinTrend_mod <- lm(TrailsA.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(TrailsARaw_LinTrend_mod)
summary(TrailsARaw_LinTrend_mod) 

#plot TrailsA.Z (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TrailsA.Z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TrailsA.Z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Trail Making Test A (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
TrailsA.Z_descrip <- describeBy(DPRC_neuropsych_data$TrailsA.Z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_TrailsAZ <- DPRC_neuropsych_data$TrailsA.Z
noNAsTrailsAZ <- na.omit(all_TrailsAZ)
mean(noNAsTrailsAZ)
sd(noNAsTrailsAZ)
#run ANOVA for TrailsA.Z
TrailsA.Z_mod <- lm(TrailsA.Z ~ Group, data = DPRC_neuropsych_data)
anova(TrailsA.Z_mod)
#effect size omnibus ANOVA
etaSquared(TrailsA.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsA.Z_mod <- glht(TrailsA.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsA.Z_mod)
confint(post_hoc_TrailsA.Z_mod)
#effect size for sig. post hoc tests
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_TrailsA.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TrailsA.Z$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TrailsA.Z$Group)
  cohensD(TrailsA.Z ~ Group, data = DPRC_neuropsych_data_CvmMCI_TrailsA.Z) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_TrailsA.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TrailsA.Z$Group <- droplevels(DPRC_neuropsych_data_CvAD_TrailsA.Z$Group)
  cohensD(TrailsA.Z ~ Group, data = DPRC_neuropsych_data_CvAD_TrailsA.Z) #this looks like Hedges' g?
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_TrailsA.Z <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, TrailsA.Z)
For_Bay_data_noNas_neuropsych_TrailsA.Z <- na.omit(For_Bay_data_noNas_neuropsych_TrailsA.Z)
anovaBF(TrailsA.Z ~ Group, data = For_Bay_data_noNas_neuropsych_TrailsA.Z) 

#plot TrailsB.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TrailsB.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TrailsB.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Trail Making Test B Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for TrailsB.Raw
TrailsB.Raw_mod <- lm(TrailsB.Raw ~ Group, data = DPRC_neuropsych_data)
anova(TrailsB.Raw_mod)
#effect size omnibus ANOVA
etaSquared(TrailsB.Raw_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsB.Raw_mod <- glht(TrailsB.Raw_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsB.Raw_mod)
confint(post_hoc_TrailsB.Raw_mod)
#check descriptive statistics per each group
TrailsB.Raw_descrip <- describeBy(DPRC_neuropsych_data$TrailsB.Raw, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_TrailsBRaw <- DPRC_neuropsych_data$TrailsB.Raw
noNAsTrailsBRaw <- na.omit(all_TrailsBRaw)
mean(noNAsTrailsBRaw)
sd(noNAsTrailsBRaw)
#effect size for sig. post hoc tests
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TrailsB.Raw$Group)
  cohensD(TrailsB.Raw ~ Group, data = DPRC_neuropsych_data_CvmMCI_TrailsB.Raw) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_TrailsB.Raw$Group)
  cohensD(TrailsB.Raw ~ Group, data = DPRC_neuropsych_data_CvAD_TrailsB.Raw) #this looks like Hedges' g?
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw$Group)
  cohensD(TrailsB.Raw ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TrailsB.Raw$Group)
  cohensD(TrailsB.Raw ~ Group, data = DPRC_neuropsych_data_SCDvAD_TrailsB.Raw) #this looks like Hedges' g? 
  #for aMCI vs. mMCI 
  DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw$Group)
  cohensD(TrailsB.Raw ~ Group, data = DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw) #this looks like Hedges' g? 
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw$Group)
  cohensD(TrailsB.Raw ~ Group, data = DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_TrailsB.Raw <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, TrailsB.Raw)
For_Bay_data_noNas_neuropsych_TrailsB.Raw <- na.omit(For_Bay_data_noNas_neuropsych_TrailsB.Raw)
anovaBF(TrailsB.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_TrailsB.Raw) 
#add in Linear Trend Analysis in Linear Regression
TrailsBRaw_LinTrend_mod <- lm(TrailsB.Raw ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(TrailsBRaw_LinTrend_mod)
summary(TrailsBRaw_LinTrend_mod) 

#plot TrailsB.Z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = TrailsB.Z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TrailsB.Z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Trail Making Test B (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#check descriptive statistics per each group
TrailsB.Z_descrip <- describeBy(DPRC_neuropsych_data$TrailsB.Z, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_TrailsBZ <- DPRC_neuropsych_data$TrailsB.Z
noNAsTrailsBZ <- na.omit(all_TrailsBZ)
mean(noNAsTrailsBZ)
sd(noNAsTrailsBZ)
#run ANOVA for TrailsB.Z
TrailsB.Z_mod <- lm(TrailsB.Z ~ Group, data = DPRC_neuropsych_data)
anova(TrailsB.Z_mod)
#effect size omnibus ANOVA
etaSquared(TrailsB.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsB.Z_mod <- glht(TrailsB.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsB.Z_mod)
confint(post_hoc_TrailsB.Z_mod)
#effect size for sig. post hoc tests
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_TrailsB.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TrailsB.Z$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TrailsB.Z$Group)
  cohensD(TrailsB.Z ~ Group, data = DPRC_neuropsych_data_CvmMCI_TrailsB.Z) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_TrailsB.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TrailsB.Z$Group <- droplevels(DPRC_neuropsych_data_CvAD_TrailsB.Z$Group)
  cohensD(TrailsB.Z ~ Group, data = DPRC_neuropsych_data_CvAD_TrailsB.Z) #this looks like Hedges' g?
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_TrailsB.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TrailsB.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TrailsB.Z$Group)
  cohensD(TrailsB.Z ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_TrailsB.Z) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_TrailsB.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TrailsB.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TrailsB.Z$Group)
  cohensD(TrailsB.Z ~ Group, data = DPRC_neuropsych_data_SCDvAD_TrailsB.Z) #this looks like Hedges' g? 
  #for aMCI vs. mMCI 
  DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Z$Group <- droplevels(DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Z$Group)
  cohensD(TrailsB.Z ~ Group, data = DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Z) #this looks like Hedges' g? 
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_TrailsB.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_TrailsB.Z$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_TrailsB.Z$Group)
  cohensD(TrailsB.Z ~ Group, data = DPRC_neuropsych_data_aMCIvAD_TrailsB.Z) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_TrailsB.Z <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, TrailsB.Z)
For_Bay_data_noNas_neuropsych_TrailsB.Z <- na.omit(For_Bay_data_noNas_neuropsych_TrailsB.Z)
anovaBF(TrailsB.Z ~ Group, data = For_Bay_data_noNas_neuropsych_TrailsB.Z) 

#plot TMT.B.TMT.A (TMT-B / TMT-A) (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TMT.B.TMT.A, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TMT.B.TMT.A, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Trail Making Test B / Trail Making Test A") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for TMT.B.TMT.A
TMT.B.TMT.A_mod <- lm(TMT.B.TMT.A ~ Group, data = DPRC_neuropsych_data)
anova(TMT.B.TMT.A_mod)
#effect size omnibus ANOVA
etaSquared(TMT.B.TMT.A_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TMT.B.TMT.A_mod <- glht(TMT.B.TMT.A_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TMT.B.TMT.A_mod)
confint(post_hoc_TMT.B.TMT.A_mod)
#check descriptive statistics per each group
TMT.B.TMT.A_descrip <- describeBy(DPRC_neuropsych_data$TMT.B.TMT.A, DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_TMTBTMTA <- DPRC_neuropsych_data$TMT.B.TMT.A
noNAsTMTBTMTA <- na.omit(all_TMTBTMTA)
mean(noNAsTMTBTMTA)
sd(noNAsTMTBTMTA)
#effect size for sig. post hoc tests
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A$Group)
  cohensD(TMT.B.TMT.A ~ Group, data = DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_CvAD_TMT.B.TMT.A$Group)
  cohensD(TMT.B.TMT.A ~ Group, data = DPRC_neuropsych_data_CvAD_TMT.B.TMT.A) #this looks like Hedges' g?
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A$Group)
  cohensD(TMT.B.TMT.A ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A$Group)
  cohensD(TMT.B.TMT.A ~ Group, data = DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_TMT.B.TMT.A <- dplyr::select(DPRC_neuropsych_data, ParticipantID, Group, TMT.B.TMT.A)
For_Bay_data_noNas_neuropsych_TMT.B.TMT.A <- na.omit(For_Bay_data_noNas_neuropsych_TMT.B.TMT.A)
anovaBF(TMT.B.TMT.A ~ Group, data = For_Bay_data_noNas_neuropsych_TMT.B.TMT.A) 
#add in Linear Trend Analysis in Linear Regression
TMTBTMTA_LinTrend_mod <- lm(TMT.B.TMT.A ~ Trend_Group + Group, data = DPRC_neuropsych_data)
anova(TMTBTMTA_LinTrend_mod)
summary(TMTBTMTA_LinTrend_mod)

#plot TMT.B.TMT.A z-scores (TMT-B / TMT-A zscore.calc.) (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TMT.B.TMT.A.zscores.calc., fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TMT.B.TMT.A.zscores.calc., color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Trail Making Test B / Trail Making Test A Z-scores") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for TMT.B.TMT.A zscores
TMT.B.TMT.Az_mod <- lm(TMT.B.TMT.A.zscores.calc. ~ Group, data = DPRC_neuropsych_data)
anova(TMT.B.TMT.Az_mod)
#effect size omnibus ANOVA
etaSquared(TMT.B.TMT.Az_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TMT.B.TMT.Az_mod <- glht(TMT.B.TMT.Az_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TMT.B.TMT.Az_mod)
confint(post_hoc_TMT.B.TMT.Az_mod)
#check descriptive statistics per each group
TMT.B.TMT.Az_descrip <- describeBy(DPRC_neuropsych_data$TMT.B.TMT.A.zscores.calc., DPRC_neuropsych_data$Group)
#find mean & SD from total sample:
all_TMTBTMTAz <- DPRC_neuropsych_data$TMT.B.TMT.A.zscores.calc.
noNAsTMTBTMTAz <- na.omit(all_TMTBTMTAz)
mean(noNAsTMTBTMTAz)
sd(noNAsTMTBTMTAz)

#####----------------ANCOVA - test with covariates (age)---------------------########

#2.Hayling Sentence completion test -------------------------------------------#
#run ANCOVA for HayBTime1.Raw
HayBTime1.Raw_covar_mod <- lm(HayBTime1.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(HayBTime1.Raw_covar_mod) #not sig.
#add in Linear Trend Analysis in Linear Regression
HayBTime1Raw_LinTrend_covar_mod <- lm(HayBTime1.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(HayBTime1Raw_LinTrend_covar_mod)
summary(HayBTime1Raw_LinTrend_covar_mod) 

#run ANCOVA for HayBTime2.Raw
HayBTime2.Raw_covar_mod <- lm(HayBTime2.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(HayBTime2.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(HayBTime2.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBTime2.Raw_covar_mod <- glht(HayBTime2.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBTime2.Raw_covar_mod) #no sig. contrasts
confint(post_hoc_HayBTime2.Raw_covar_mod)
#add in Linear Trend Analysis in Linear Regression
HayBTime2Raw_LinTrend_mod <- lm(HayBTime2.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(HayBTime2Raw_LinTrend_covar_mod)
summary(HayBTime2Raw_LinTrend_covar_mod) 

#run ANCOVA for HayBCatA.Raw
HayBCatA.Raw_covar_mod <- lm(HayBCatA.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(HayBCatA.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(HayBCatA.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBCatA.Raw_covar_mod <- glht(HayBCatA.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBCatA.Raw_covar_mod) 
confint(post_hoc_HayBCatA.Raw_covar_mod)
t_value_effect_size <- summary(post_hoc_HayBCatA.Raw_covar_mod) 
#effect size for sig. post hoc tests
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_HayBCatA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_HayBCatA.Raw$HayBCatA.Raw,DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$HayBCatA.Raw,DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$HayBCatA.Raw,DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for mMCI vs. AD
  DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 4 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$HayBCatA.Raw,DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 4'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
HayBCatARaw_LinTrend_covar_mod <- lm(HayBCatA.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(HayBCatARaw_LinTrend_covar_mod)
summary(HayBCatARaw_LinTrend_covar_mod) 
  
#run ANCOVA for HayBCatB.Raw
HayBCatB.Raw_covar_mod <- lm(HayBCatB.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(HayBCatB.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(HayBCatB.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBCatB.Raw_covar_mod <- glht(HayBCatB.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBCatB.Raw_covar_mod) 
confint(post_hoc_HayBCatB.Raw_covar_mod)
t_value_effect_size <- summary(post_hoc_HayBCatB.Raw_covar_mod) 
#effect size with covariate
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$HayBCatB.Raw,DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
HayBCatBRaw_LinTrend_covar_mod <- lm(HayBCatB.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(HayBCatBRaw_LinTrend_covar_mod)
summary(HayBCatBRaw_LinTrend_covar_mod) 

#3.Stroop Test ----------------------------------------------------------------#
#run ANCOVA for ColourNaming.Raw
ColourNaming.Raw_covar_mod <- lm(ColorNaming.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(ColourNaming.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(ColourNaming.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_ColourNaming.Raw_covar_mod <- glht(ColourNaming.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_ColourNaming.Raw_covar_mod) 
confint(post_hoc_ColourNaming.Raw_covar_mod)
t_value_effect_size <- summary(post_hoc_ColourNaming.Raw_covar_mod) 
#effect size with covariate
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw$ColorNaming.Raw,DPRC_neuropsych_data_SCDvAD_ColorNaming.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
ColorNamingRaw_LinTrend_covar_mod <- lm(ColorNaming.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(ColorNamingRaw_LinTrend_covar_mod)
summary(ColorNamingRaw_LinTrend_covar_mod) 

#run ANCOVA for WordReading.Raw
WordReading.Raw_covar_mod <- lm(WordReading.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(WordReading.Raw_covar_mod) #not sig.
#add in Linear Trend Analysis in Linear Regression
WordReadingRaw_LinTrend_covar_mod <- lm(WordReading.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(WordReadingRaw_LinTrend_covar_mod)
summary(WordReadingRaw_LinTrend_covar_mod) 

#run ANCOVA for Inhibition.Raw
Inhibition.Raw_covar_mod <- lm(Inhibition.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(Inhibition.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Raw_covar_mod <- glht(Inhibition.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Raw_covar_mod) 
confint(post_hoc_Inhibition.Raw_covar_mod)
t_value_effect_size <- summary(post_hoc_Inhibition.Raw_covar_mod) 
#effect size with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_Inhibition.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_Inhibition.Raw$Inhibition.Raw,DPRC_neuropsych_data_CvmMCI_Inhibition.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Inhibition.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_Inhibition.Raw$Inhibition.Raw,DPRC_neuropsych_data_CvAD_Inhibition.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw$Inhibition.Raw,DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Inhibition.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_Inhibition.Raw$Inhibition.Raw,DPRC_neuropsych_data_SCDvAD_Inhibition.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvAD_Inhibition$Inhibition.Raw,DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
InhibitionRaw_LinTrend_covar_mod <- lm(Inhibition.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(InhibitionRaw_LinTrend_covar_mod)
summary(InhibitionRaw_LinTrend_covar_mod) 

#run ANCOVA for Inhibition.Colour.Naming (Inhibition / ColourNaming)
Inhibition.Colour.Naming_covar_mod <- lm(Inhibition.Colour.Naming ~ Group+Age, data = DPRC_neuropsych_data)
anova(Inhibition.Colour.Naming_covar_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Colour.Naming_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Colour.Naming_covar_mod <- glht(Inhibition.Colour.Naming_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Colour.Naming_covar_mod) 
confint(post_hoc_Inhibition.Colour.Naming_covar_mod)
t_value_effect_size <- summary(post_hoc_Inhibition.Colour.Naming_covar_mod) 
#effect size with covariate
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw$Inhibition.Colour.Naming,DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Inhibition.Colour.Naming,DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Inhibition.Colour.Naming,DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
Inhibition.Colour.Naming_LinTrend_covar_mod <- lm(Inhibition.Colour.Naming ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(Inhibition.Colour.Naming_LinTrend_covar_mod)
summary(Inhibition.Colour.Naming_LinTrend_covar_mod) 

#run ANCOVA for Inhibition.Word.Reading (Inhibition / WordReading)
Inhibition.Word.Reading_covar_mod <- lm(Inhibition.Word.Reading ~ Group+Age, data = DPRC_neuropsych_data)
anova(Inhibition.Word.Reading_covar_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Word.Reading_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Word.Reading_covar_mod <- glht(Inhibition.Word.Reading_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Word.Reading_covar_mod) 
confint(post_hoc_Inhibition.Word.Reading_covar_mod)
t_value_effect_size <- summary(post_hoc_Inhibition.Word.Reading_covar_mod) 
#effect size with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading$Inhibition.Word.Reading,DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading$Inhibition.Word.Reading,DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading$Inhibition.Word.Reading,DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvAD_Inhibition$Inhibition.Word.Reading,DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
Inhibition.Word.Reading_LinTrend_covar_mod <- lm(Inhibition.Word.Reading ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(Inhibition.Word.Reading_LinTrend_covar_mod)
summary(Inhibition.Word.Reading_LinTrend_covar_mod) 


#4.TMT ------------------------------------------------------------------------#
#run ANCOVA for TrailsA.Raw
TrailsA.Raw_covar_mod <- lm(TrailsA.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(TrailsA.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(TrailsA.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsA.Raw_covar_mod <- glht(TrailsA.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsA.Raw_covar_mod) 
confint(post_hoc_TrailsA.Raw_covar_mod)
t_value_effect_size <- summary(post_hoc_TrailsA.Raw_covar_mod) 
#effect size with covariate
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_TrailsA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_TrailsA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_TrailsA.Raw$TrailsA.Raw,DPRC_neuropsych_data_CvAD_TrailsA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TrailsA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_TrailsA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_TrailsA.Raw$TrailsA.Raw,DPRC_neuropsych_data_SCDvAD_TrailsA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvAD_TrailsA$TrailsA.Raw,DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
TrailsA.Raw_LinTrend_covar_mod <- lm(TrailsA.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(TrailsA.Raw_LinTrend_covar_mod)
summary(TrailsA.Raw_LinTrend_covar_mod) 

#run ANCOVA for TrailsB.Raw
TrailsB.Raw_covar_mod <- lm(TrailsB.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(TrailsB.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(TrailsB.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsB.Raw_covar_mod <- glht(TrailsB.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsB.Raw_covar_mod) 
confint(post_hoc_TrailsB.Raw_covar_mod)
t_value_effect_size <- summary(post_hoc_TrailsB.Raw_covar_mod) 
#effect size with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_TrailsB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_TrailsB.Raw$TrailsB.Raw,DPRC_neuropsych_data_CvmMCI_TrailsB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_TrailsB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_TrailsB.Raw$TrailsB.Raw,DPRC_neuropsych_data_CvAD_TrailsB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw$TrailsB.Raw,DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_TrailsB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_TrailsB.Raw$TrailsB.Raw,DPRC_neuropsych_data_SCDvAD_TrailsB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. mMCI
  DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvmMCI_TrailsB$TrailsB.Raw,DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvAD_TrailsB$TrailsB.Raw,DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
TrailsB.Raw_LinTrend_covar_mod <- lm(TrailsB.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(TrailsB.Raw_LinTrend_covar_mod)
summary(TrailsB.Raw_LinTrend_covar_mod) 

#run ANCOVA for TMT.B.TMT.A
TMT.B.TMT.A_covar_mod <- lm(TMT.B.TMT.A ~ Group+Age, data = DPRC_neuropsych_data)
anova(TMT.B.TMT.A_covar_mod)
#effect size omnibus ANOVA
etaSquared(TMT.B.TMT.A_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TMT.B.TMT.A_covar_mod <- glht(TMT.B.TMT.A_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TMT.B.TMT.A_covar_mod) 
confint(post_hoc_TMT.B.TMT.A_covar_mod)
t_value_effect_size <- summary(post_hoc_TMT.B.TMT.A_covar_mod) 
#effect size with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A$TMT.B.TMT.A,DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_CvAD_TMT.B.TMT.A$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_TMT.B.TMT.A, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_TMT.B.TMT.A$TMT.B.TMT.A,DPRC_neuropsych_data_CvAD_TMT.B.TMT.A$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A$TMT.B.TMT.A,DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A$TMT.B.TMT.A,DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
TMT.B.TMT.A_LinTrend_covar_mod <- lm(TMT.B.TMT.A ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(TMT.B.TMT.A_LinTrend_covar_mod)
summary(TMT.B.TMT.A_LinTrend_covar_mod) 

#5. Letter and Verbal Fluency -------------------------------------------------#
#run ANCOVA for LetFluency.Raw
LetFluency.Raw_covar_mod <- lm(LetFluency.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(LetFluency.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(LetFluency.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_LetFluency.Raw_covar_mod <- glht(LetFluency.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_LetFluency.Raw_covar_mod) 
confint(post_hoc_LetFluency.Raw_covar_mod)
t_value_effect_size<-summary(post_hoc_LetFluency.Raw_covar_mod) 
#effect size with covariate
  #for C vs. SCD
  DPRC_neuropsych_data_CvSCD_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 2)
  DPRC_neuropsych_data_CvSCD_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvSCD_LetFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_LetFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvSCD_LetFluency.Raw$LetFluency.Raw,DPRC_neuropsych_data_CvSCD_LetFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$LetFluency.Raw,DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_LetFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_LetFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_LetFluency.Raw$LetFluency.Raw,DPRC_neuropsych_data_CvmMCI_LetFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
LetFluency.Raw_LinTrend_covar_mod <- lm(LetFluency.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(LetFluency.Raw_LinTrend_covar_mod)
summary(LetFluency.Raw_LinTrend_covar_mod) 

#run ANCOVA for CatFluency.Raw
CatFluency.Raw_covar_mod <- lm(CatFluency.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(CatFluency.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(CatFluency.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_CatFluency.Raw_covar_mod <- glht(CatFluency.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_CatFluency.Raw_covar_mod) 
confint(post_hoc_CatFluency.Raw_covar_mod)
t_value_effect_size<-summary(post_hoc_CatFluency.Raw_covar_mod) 
#effect size with covariate
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_CvAD_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
CatFluency.Raw_LinTrend_covar_mod <- lm(CatFluency.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(CatFluency.Raw_LinTrend_covar_mod)
summary(CatFluency.Raw_LinTrend_covar_mod) 

#run ANCOVA for Switching.Raw
Switching.Raw_covar_mod <- lm(Switching.Raw ~ Group+Age, data = DPRC_neuropsych_data)
anova(Switching.Raw_covar_mod)
#effect size omnibus ANOVA
etaSquared(Switching.Raw_covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Switching.Raw_covar_mod <- glht(Switching.Raw_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Switching.Raw_covar_mod) 
confint(post_hoc_Switching.Raw_covar_mod)
t_value_effect_size <- summary(post_hoc_Switching.Raw_covar_mod) 
#effect size with covariate
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvaMCI_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_CvaMCI_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_CvmMCI_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_CvAD_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. aMCI
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_SCDvAD_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvAD_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_aMCIvAD_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#add in Linear Trend Analysis in Linear Regression
Switching.Raw_LinTrend_covar_mod <- lm(Switching.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(Switching.Raw_LinTrend_covar_mod)
summary(Switching.Raw_LinTrend_covar_mod) 

#----------ANCOVA - test with covariates (age & sex)---------------------------------#####

#2.Hayling Sentence completion test -------------------------------------------#
#run ANCOVA for HayBTime1.Raw
HayBTime1.Raw_2covar_mod <- lm(HayBTime1.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(HayBTime1.Raw_2covar_mod) #not sig.
#add in Linear Trend Analysis in Linear Regression
HayBTime1Raw_LinTrend_2covar_mod <- lm(HayBTime1.Raw ~ Trend_Group + Group + Age + Sex, data = DPRC_neuropsych_data)
anova(HayBTime1Raw_LinTrend_2covar_mod)
summary(HayBTime1Raw_LinTrend_2covar_mod) 

#run ANCOVA for HayBTime2.Raw
HayBTime2.Raw_2covar_mod <- lm(HayBTime2.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(HayBTime2.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(HayBTime2.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBTime2.Raw_2covar_mod <- glht(HayBTime2.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBTime2.Raw_2covar_mod) #no sig. contrasts
confint(post_hoc_HayBTime2.Raw_2covar_mod)
#add in Linear Trend Analysis in Linear Regression
HayBTime2Raw_LinTrend_2covar_mod <- lm(HayBTime2.Raw ~ Trend_Group + Group + Age + Sex, data = DPRC_neuropsych_data)
anova(HayBTime2Raw_LinTrend_2covar_mod)
summary(HayBTime2Raw_LinTrend_2covar_mod) 

#run ANCOVA for HayBCatA.Raw
HayBCatA.Raw_2covar_mod <- lm(HayBCatA.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(HayBCatA.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(HayBCatA.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBCatA.Raw_2covar_mod <- glht(HayBCatA.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBCatA.Raw_2covar_mod) 
confint(post_hoc_HayBCatA.Raw_2covar_mod)
t_value_effect_size <- summary(post_hoc_HayBCatA.Raw_2covar_mod) 
#effect size for sig. post hoc tests
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_HayBCatA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_HayBCatA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for mMCI vs. AD
  DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 4 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 4'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
HayBCatARaw_LinTrend_2covar_mod <- lm(HayBCatA.Raw ~ Trend_Group + Group + Age + Sex, data = DPRC_neuropsych_data)
anova(HayBCatARaw_LinTrend_2covar_mod)
summary(HayBCatARaw_LinTrend_2covar_mod) 

#run ANCOVA for HayBCatB.Raw
HayBCatB.Raw_2covar_mod <- lm(HayBCatB.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(HayBCatB.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(HayBCatB.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_HayBCatB.Raw_2covar_mod <- glht(HayBCatB.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_HayBCatB.Raw_2covar_mod) 
confint(post_hoc_HayBCatB.Raw_2covar_mod)
t_value_effect_size <- summary(post_hoc_HayBCatB.Raw_2covar_mod) 
#effect size with covariate
#for SCD vs. AD
DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Group)
group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw, Group) #count number of participants per group
mult.r_value_2covar_mod<-summary(lm(HayBCatB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
HayBCatBRaw_LinTrend_2covar_mod <- lm(HayBCatB.Raw ~ Trend_Group + Group + Age, data = DPRC_neuropsych_data)
anova(HayBCatBRaw_LinTrend_2covar_mod)
summary(HayBCatBRaw_LinTrend_2covar_mod) 

#3.Stroop Test ----------------------------------------------------------------#
#run ANCOVA for ColourNaming.Raw
ColourNaming.Raw_2covar_mod <- lm(ColorNaming.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(ColourNaming.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(ColourNaming.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_ColourNaming.Raw_2covar_mod <- glht(ColourNaming.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_ColourNaming.Raw_2covar_mod) 
confint(post_hoc_ColourNaming.Raw_2covar_mod)
t_value_effect_size <- summary(post_hoc_ColourNaming.Raw_2covar_mod) 
#add in Linear Trend Analysis in Linear Regression
ColorNamingRaw_LinTrend_2covar_mod <- lm(ColorNaming.Raw ~ Trend_Group + Group+Age+Sex, data = DPRC_neuropsych_data)
anova(ColorNamingRaw_LinTrend_2covar_mod)
summary(ColorNamingRaw_LinTrend_2covar_mod) 

#run ANCOVA for WordReading.Raw
WordReading.Raw_2covar_mod <- lm(WordReading.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(WordReading.Raw_2covar_mod) #not sig.
#add in Linear Trend Analysis in Linear Regression
WordReadingRaw_LinTrend_2covar_mod <- lm(WordReading.Raw ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(WordReadingRaw_LinTrend_2covar_mod)
summary(WordReadingRaw_LinTrend_2covar_mod) 

#run ANCOVA for Inhibition.Raw
Inhibition.Raw_2covar_mod <- lm(Inhibition.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(Inhibition.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Raw_2covar_mod <- glht(Inhibition.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Raw_2covar_mod) 
confint(post_hoc_Inhibition.Raw_2covar_mod)
t_value_effect_size <- summary(post_hoc_Inhibition.Raw_2covar_mod) 
#effect size with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_Inhibition.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_Inhibition.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Inhibition.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_Inhibition.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvmMCI_Inhibition.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Inhibition.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_Inhibition.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Raw ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvAD_Inhibition.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
InhibitionRaw_LinTrend_2covar_mod <- lm(Inhibition.Raw ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(InhibitionRaw_LinTrend_2covar_mod)
summary(InhibitionRaw_LinTrend_2covar_mod) 

#run ANCOVA for Inhibition.Colour.Naming (Inhibition / ColourNaming)
Inhibition.Colour.Naming_2covar_mod <- lm(Inhibition.Colour.Naming ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(Inhibition.Colour.Naming_2covar_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Colour.Naming_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Colour.Naming_2covar_mod <- glht(Inhibition.Colour.Naming_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Colour.Naming_2covar_mod) 
confint(post_hoc_Inhibition.Colour.Naming_2covar_mod)
t_value_effect_size <- summary(post_hoc_Inhibition.Colour.Naming_2covar_mod) 
#effect size with covariate
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Colour.Naming ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_Inhibition.Colour.Naming.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Colour.Naming ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_Inhibition.Colour.Naming.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Colour.Naming ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvAD_Inhibition.Colour.Naming.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
Inhibition.Colour.Naming_LinTrend_2covar_mod <- lm(Inhibition.Colour.Naming ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(Inhibition.Colour.Naming_LinTrend_2covar_mod)
summary(Inhibition.Colour.Naming_LinTrend_2covar_mod) 

#run ANCOVA for Inhibition.Word.Reading (Inhibition / WordReading)
Inhibition.Word.Reading_2covar_mod <- lm(Inhibition.Word.Reading ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(Inhibition.Word.Reading_2covar_mod)
#effect size omnibus ANOVA
etaSquared(Inhibition.Word.Reading_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Inhibition.Word.Reading_2covar_mod <- glht(Inhibition.Word.Reading_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Inhibition.Word.Reading_2covar_mod) 
confint(post_hoc_Inhibition.Word.Reading_2covar_mod)
t_value_effect_size <- summary(post_hoc_Inhibition.Word.Reading_2covar_mod) 
#effect size with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Word.Reading ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_Inhibition.Word.Reading)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Word.Reading ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_Inhibition.Word.Reading)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Word.Reading ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_Inhibition.Word.Reading)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Inhibition.Word.Reading ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvAD_Inhibition.Word.Reading)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
Inhibition.Word.Reading_LinTrend_2covar_mod <- lm(Inhibition.Word.Reading ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(Inhibition.Word.Reading_LinTrend_2covar_mod)
summary(Inhibition.Word.Reading_LinTrend_2covar_mod) 

#4.TMT ------------------------------------------------------------------------#
#run ANCOVA for TrailsA.Raw
TrailsA.Raw_2covar_mod <- lm(TrailsA.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(TrailsA.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(TrailsA.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsA.Raw_2covar_mod <- glht(TrailsA.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsA.Raw_2covar_mod) 
confint(post_hoc_TrailsA.Raw_2covar_mod)
t_value_effect_size <- summary(post_hoc_TrailsA.Raw_2covar_mod) 
#effect size with covariate
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_TrailsA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_TrailsA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_TrailsA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TrailsA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_TrailsA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_TrailsA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvAD_TrailsA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
TrailsA.Raw_LinTrend_2covar_mod <- lm(TrailsA.Raw ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(TrailsA.Raw_LinTrend_2covar_mod)
summary(TrailsA.Raw_LinTrend_2covar_mod) 

#run ANCOVA for TrailsB.Raw
TrailsB.Raw_2covar_mod <- lm(TrailsB.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(TrailsB.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(TrailsB.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsB.Raw_2covar_mod <- glht(TrailsB.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsB.Raw_2covar_mod) 
confint(post_hoc_TrailsB.Raw_2covar_mod)
t_value_effect_size <- summary(post_hoc_TrailsB.Raw_2covar_mod) 
#effect size with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_TrailsB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_TrailsB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_TrailsB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_TrailsB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvmMCI_TrailsB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_TrailsB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_TrailsB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. mMCI
  DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvmMCI_TrailsB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvAD_TrailsB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
TrailsB.Raw_LinTrend_2covar_mod <- lm(TrailsB.Raw ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(TrailsB.Raw_LinTrend_2covar_mod)
summary(TrailsB.Raw_LinTrend_2covar_mod) 

#run ANCOVA for TMT.B.TMT.A
TMT.B.TMT.A_2covar_mod <- lm(TMT.B.TMT.A ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(TMT.B.TMT.A_2covar_mod)
#effect size omnibus ANOVA
etaSquared(TMT.B.TMT.A_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TMT.B.TMT.A_2covar_mod <- glht(TMT.B.TMT.A_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TMT.B.TMT.A_2covar_mod) 
confint(post_hoc_TMT.B.TMT.A_2covar_mod)
t_value_effect_size <- summary(post_hoc_TMT.B.TMT.A_2covar_mod) 
#effect size with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TMT.B.TMT.A ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_TMT.B.TMT.A)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_CvAD_TMT.B.TMT.A$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_TMT.B.TMT.A, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TMT.B.TMT.A ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_TMT.B.TMT.A)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TMT.B.TMT.A ~ Age + Sex, data = DPRC_neuropsych_data_SCDvmMCI_TMT.B.TMT.A)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TMT.B.TMT.A ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_TMT.B.TMT.A)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
TMT.B.TMT.A_LinTrend_2covar_mod <- lm(TMT.B.TMT.A ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(TMT.B.TMT.A_LinTrend_2covar_mod)
summary(TMT.B.TMT.A_LinTrend_2covar_mod) 

#5. Letter and Verbal Fluency -------------------------------------------------#
#run ANCOVA for LetFluency.Raw
LetFluency.Raw_2covar_mod <- lm(LetFluency.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(LetFluency.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(LetFluency.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_LetFluency.Raw_2covar_mod <- glht(LetFluency.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_LetFluency.Raw_2covar_mod) 
confint(post_hoc_LetFluency.Raw_2covar_mod)
t_value_effect_size<-summary(post_hoc_LetFluency.Raw_2covar_mod) 
  #effect size with covariate
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(LetFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvaMCI_LetFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_LetFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_LetFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(LetFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_LetFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
LetFluency.Raw_LinTrend_2covar_mod <- lm(LetFluency.Raw ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(LetFluency.Raw_LinTrend_2covar_mod)
summary(LetFluency.Raw_LinTrend_2covar_mod) 

#run ANCOVA for CatFluency.Raw
CatFluency.Raw_2covar_mod <- lm(CatFluency.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(CatFluency.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(CatFluency.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_CatFluency.Raw_2covar_mod <- glht(CatFluency.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_CatFluency.Raw_2covar_mod) 
confint(post_hoc_CatFluency.Raw_2covar_mod)
t_value_effect_size<-summary(post_hoc_CatFluency.Raw_2covar_mod) 
#effect size with covariate
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvaMCI_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
CatFluency.Raw_LinTrend_2covar_mod <- lm(CatFluency.Raw ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(CatFluency.Raw_LinTrend_2covar_mod)
summary(CatFluency.Raw_LinTrend_2covar_mod) 

#run ANCOVA for Switching.Raw
Switching.Raw_2covar_mod <- lm(Switching.Raw ~ Group+Age+Sex, data = DPRC_neuropsych_data)
anova(Switching.Raw_2covar_mod)
#effect size omnibus ANOVA
etaSquared(Switching.Raw_2covar_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Switching.Raw_2covar_mod <- glht(Switching.Raw_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Switching.Raw_2covar_mod) 
confint(post_hoc_Switching.Raw_2covar_mod)
t_value_effect_size <- summary(post_hoc_Switching.Raw_2covar_mod) 
#effect size with covariate
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvaMCI_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. aMCI
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvaMCI_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvmMCI_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvAD_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#add in Linear Trend Analysis in Linear Regression
Switching.Raw_LinTrend_2covar_mod <- lm(Switching.Raw ~ Trend_Group+Group+Age+Sex, data = DPRC_neuropsych_data)
anova(Switching.Raw_LinTrend_2covar_mod)
summary(Switching.Raw_LinTrend_2covar_mod) 





















#------------------------------------------------------------------------------#





#######----------------- Longitudinal (F0 vs. F2) analysis -----------------########

#read in csv files (participant file)
DPRC_neuropsych_data <- read.csv("longitudinal_DPRC_neuropsych_data_lined_up_valid_participants.csv")

#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'

#Add in longitudinal values for participants
Individual_number <- c(1:124, 1:124)
DPRC_neuropsych_data$Individual_number <- Individual_number

#convert variables
DPRC_neuropsych_data$ParticipantID <- as.factor(DPRC_neuropsych_data$ParticipantID)
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Sex <- as.factor(DPRC_neuropsych_data$Sex)
DPRC_neuropsych_data$Sex_binary <- as.factor(DPRC_neuropsych_data$Sex_binary)
DPRC_neuropsych_data$Timepoint <- as.factor(DPRC_neuropsych_data$Timepoint)
DPRC_neuropsych_data$Individual_number <- as.factor(DPRC_neuropsych_data$Individual_number)

#----------------plot data to visualise & run stat tests ----------------------#

#longitudinal example for age (F0 vs. F2)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Age, group=interaction(Group, Timepoint))) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group), position = position_dodge(.9)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group), position = position_dodge(.9)) + 
  ylim(50, 95) +
  xlab("Group") + 
  ylab("Age") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  geom_violin(trim = FALSE, alpha = .5, aes(fill = Group, colour = Group), size = 1)

#take a subset of the DPRC data - only baseline data (F0) for certain graphs or analysis
baseline_DPRC_neuropsych <- DPRC_neuropsych_data[ which(DPRC_neuropsych_data$Timepoint=='F0'), ]


####-----------------------------MANOVAs------------------------------------####
#for raw data
manova_exec_func_raw_mod <- manova(cbind(TrailsA.Raw, 
                                         TrailsB.Raw, 
                                         ColorNaming.Raw, 
                                         WordReading.Raw, 
                                         Inhibition.Raw, 
                                         LetFluency.Raw, 
                                         CatFluency.Raw, 
                                         Switching.Raw, 
                                         HayBTime1.Raw, 
                                         HayBTime2.Raw, 
                                         HayBCatA.Raw, 
                                         HayBCatB.Raw) ~ Group, data = baseline_DPRC_neuropsych) 
#overall model
summary(manova_exec_func_raw_mod)
#anova outputs
summary.aov(manova_exec_func_raw_mod)

#for z-scores
manova_exec_func_z_mod <- manova(cbind(TrailsA.Z, 
                                       TrailsB.Z, 
                                       ColorNaming.Z, 
                                       WordReading.Z, 
                                       Inhibition.Z, 
                                       LetFluency.Z, 
                                       CatFluency.Z, 
                                       Switching.z, 
                                       HayBTime1.z, 
                                       HayBTime2.z, 
                                       HayBCatA.z, 
                                       HayBCatB.z) ~ Group, data = baseline_DPRC_neuropsych) 
#overall model
summary(manova_exec_func_z_mod)
#anova outputs
summary.aov(manova_exec_func_z_mod)

#plot the manova data using the z-scores of the variables, so that it is all on the same scale. 
#put exec func variable z-scores onto a new dataset, as long format
exec_func_zscores_data <- dplyr::select(baseline_DPRC_neuropsych, 
                                        ParticipantID,
                                        Group,
                                        TrailsA.Z, 
                                        TrailsB.Z, 
                                        ColorNaming.Z, 
                                        WordReading.Z, 
                                        Inhibition.Z, 
                                        LetFluency.Z, 
                                        CatFluency.Z,
                                        Switching.z,
                                        HayBTime1.z,
                                        HayBTime2.z,
                                        HayBCatA.z,
                                        HayBCatB.z)
#put into long format
exec_func_zscores_data_long <- gather(exec_func_zscores_data, 
                                      "Test",
                                      "Z_scores", 
                                      TrailsA.Z, 
                                      TrailsB.Z, 
                                      ColorNaming.Z, 
                                      WordReading.Z, 
                                      Inhibition.Z, 
                                      LetFluency.Z, 
                                      CatFluency.Z, 
                                      Switching.z, 
                                      HayBTime1.z, 
                                      HayBTime2.z, 
                                      HayBCatA.z, 
                                      HayBCatB.z)

#plot data - z-score for overall 
ggplot(exec_func_zscores_data_long, aes(x = Group, y = Z_scores, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Z_scores, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Executive Functioning (Z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()


#for raw data - processing speeds variables
manova_proc_speed_raw_mod <- manova(cbind(TrailsA.Raw, 
                                          ColorNaming.Raw, 
                                          WordReading.Raw, 
                                          HayBTime1.Raw) ~ Group, data = baseline_DPRC_neuropsych) 

#overall model
summary(manova_proc_speed_raw_mod)
#anova outputs
summary.aov(manova_proc_speed_raw_mod)

#for z-scores
manova_proc_speed_z_mod <- manova(cbind(TrailsA.Z, 
                                        ColorNaming.Z, 
                                        WordReading.Z, 
                                        HayBTime1.z) ~ Group, data = baseline_DPRC_neuropsych) 


#overall model
summary(manova_proc_speed_z_mod)
#anova outputs
summary.aov(manova_proc_speed_z_mod)




#plot the manova data using the z-scores of the variables, so that it is all on the same scale. 
#put proc speed variable z-scores onto a new dataset, as long format
proc_speed_zscores_data <- dplyr::select(baseline_DPRC_neuropsych, 
                                         ParticipantID,
                                         Group,
                                         TrailsA.Z, 
                                         ColorNaming.Z, 
                                         WordReading.Z, 
                                         HayBTime1.z)
#put into long format
proc_speed_zscores_data_long <- gather(exec_func_zscores_data, 
                                       "Processing_Speeds",
                                       "Z_scores", 
                                       TrailsA.Z, 
                                       ColorNaming.Z, 
                                       WordReading.Z, 
                                       HayBTime1.z)


#plot data - z-score for processing speeds 
ggplot(proc_speed_zscores_data_long, aes(x = Group, y = Z_scores, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Z_scores, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Processing Speeds (Z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

#get descriptives of the collated z-scores - e.g. 'processing speeds score' mean beta values
proc_speed_z_descrip <- describeBy((proc_speed_zscores_data$TrailsA.Z + 
                                      proc_speed_zscores_data$ColorNaming.Z +
                                      proc_speed_zscores_data$WordReading.Z +
                                      proc_speed_zscores_data$HayBTime1.z) / (length(proc_speed_zscores_data) - 2), proc_speed_zscores_data$Group)
#find mean & SD from total sample: 
all_proc_speed_zscores_data <- proc_speed_zscores_data[,3:6]
noNAsproc_speed_zscores_data <- na.omit(all_proc_speed_zscores_data)
mean(as.matrix(noNAsproc_speed_zscores_data))
sd(as.matrix(noNAsproc_speed_zscores_data))

#for inhibition variables
manova_inhibition_raw_mod <- manova(cbind(TrailsB.Raw, 
                                          Inhibition.Raw, 
                                          Switching.Raw, 
                                          HayBTime2.Raw, 
                                          HayBCatA.Raw, 
                                          HayBCatB.Raw) ~ Group, data = baseline_DPRC_neuropsych) 
#overall model
summary(manova_inhibition_raw_mod)
#anova outputs
summary.aov(manova_inhibition_raw_mod)

#for z-scores
manova_inhibition_z_mod <- manova(cbind(TrailsB.Z, 
                                        Inhibition.Z, 
                                        Switching.z, 
                                        HayBTime2.z, 
                                        HayBCatA.z, 
                                        HayBCatB.z) ~ Group, data = baseline_DPRC_neuropsych) 
#overall model
summary(manova_inhibition_z_mod)
#anova outputs
summary.aov(manova_inhibition_z_mod)

#plot the manova data using the z-scores of the variables, so that it is all on the same scale. 
#put exec func variable z-scores onto a new dataset, as long format
inhibition_zscores_data <- dplyr::select(baseline_DPRC_neuropsych, 
                                         ParticipantID,
                                         Group,
                                         TrailsB.Z, 
                                         Inhibition.Z, 
                                         Switching.z,
                                         HayBTime2.z,
                                         HayBCatA.z,
                                         HayBCatB.z)
#put into long format
inhibition_zscores_data_long <- gather(inhibition_zscores_data, 
                                       "Inhibition",
                                       "Z_scores", 
                                       TrailsB.Z, 
                                       Inhibition.Z, 
                                       Switching.z, 
                                       HayBTime2.z, 
                                       HayBCatA.z, 
                                       HayBCatB.z)

#plot data - z-score for overall 
ggplot(inhibition_zscores_data_long, aes(x = Group, y = Z_scores, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Z_scores, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Inhibition (Z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

#get descriptives of the collated z-scores - e.g. 'executive function score' mean beta values
inhibition_z_descrip <- describeBy((inhibition_zscores_data$TrailsB.Z + 
                                      inhibition_zscores_data$Inhibition.Z +
                                      inhibition_zscores_data$Switching.z +
                                      inhibition_zscores_data$HayBTime2.z +
                                      inhibition_zscores_data$HayBCatA.z +
                                      inhibition_zscores_data$HayBCatB.z) / (length(inhibition_zscores_data) - 2), inhibition_zscores_data$Group)
#find mean & SD from total sample: 
all_inhibition_zscores_data <- inhibition_zscores_data[,3:8]
noNAsinhibition_zscores_data <- na.omit(all_inhibition_zscores_data)
mean(as.matrix(noNAsinhibition_zscores_data))
sd(as.matrix(noNAsinhibition_zscores_data))


####---------------------------------ANOVAs---------------------------------####
#run ANOVAs on your data:
#1.Test of premorbid functioning (TOPF) ---------------------------------------#
#plot TOPF.Raw (raincloud plot) - only for baseline
ggplot(baseline_DPRC_neuropsych, aes(x = Group, y = TOPF.Raw, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TOPF.Raw, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("TOPF Score (Raw)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for TOPF.Raw
TOPF.Raw_mod <- lm(TOPF.Raw ~ Group, data = baseline_DPRC_neuropsych)
anova(TOPF.Raw_mod)
#check descriptive statistics per each group
TOPF.Raw_descrip <- describeBy(baseline_DPRC_neuropsych$TOPF.Raw, baseline_DPRC_neuropsych$Group)
#whole sample descriptives
noNAs_baseline_DPRC_neuropsych <- na.omit(baseline_DPRC_neuropsych)
mean(noNAs_baseline_DPRC_neuropsych$TOPF.Raw)
sd(noNAs_baseline_DPRC_neuropsych$TOPF.Raw)
#plot TOPF.Z (raincloud plot)
ggplot(baseline_DPRC_neuropsych, aes(x = Group, y = TOPF.Z, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TOPF.Z, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("TOPF Score (z-scores)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
#run ANOVA for TOPF.Z
TOPF.Z_mod <- lm(TOPF.Z ~ Group, data = baseline_DPRC_neuropsych)
anova(TOPF.Z_mod)

#2.Hayling Sentence completion test -------------------------------------------#
#plot HayBTime1.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBTime1.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBTime1.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 1 Time Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()

#run mixed design, 2 x 5 ANOVA for HayBTime1.Raw
#aov_HayBTime1Raw <- aov(HayBTime1.Raw ~ Group*Timepoint + Error(ParticipantID/Timepoint), data = DPRC_neuropsych_data)
#summary(aov_HayBTime1Raw)
#aov_HayBTime1Raw <- Anova(lm(HayBTime1.Raw ~ Group*Timepoint, data=DPRC_neuropsych_data), type = "III") #use this anova test to account for unbalanced designs/sample sizes
#aov_HayBTime1Raw
aov_HayBTime1Raw <- anova_test(data=DPRC_neuropsych_data, dv=HayBTime1.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_HayBTime1Raw)
#effect size
#eta_squared(aov_HayBTime1Raw)
#check descriptive statistics per each group, per each timepoint
HayBTime1.Raw_descrip <- describeBy(DPRC_neuropsych_data$HayBTime1.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_HayBTime1Raw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_HayBTime1Raw <- na.omit(F0_HayBTime1Raw$HayBTime1.Raw)
mean(noNAsF0_HayBTime1Raw)
sd(noNAsF0_HayBTime1Raw)
#F2
F2_HayBTime1Raw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_HayBTime1Raw <- na.omit(F2_HayBTime1Raw$HayBTime1.Raw)
mean(noNAsF2_HayBTime1Raw)
sd(noNAsF2_HayBTime1Raw)
#post hoc test
#compare main effects - Group
#GroupME_HayBTime1 <- emmeans(aov_HayBTime1Raw, ~ Group)
#pairs(GroupME_HayBTime1ME)
#compare main effects - Timepoint
#TimepointME_HayBTime1 <- emmeans(aov_HayBTime1Raw, ~ Timepoint)
#pairs(TimepointME_HayBTime1ME)
#check mean values
#lsmeans(aov_HayBTime1Raw, pairwise ~ Group | Timepoint)
#lsmeans(aov_HayBTime1Raw, pairwise ~ Timepoint | Group)
#lsmeans(aov_HayBTime1Raw, c("Group", "Timepoint"))
#run linear mixed modelling version of ANOVA & posthoc
#lmeHayBTime1mod <- lmer(HayBTime1.Raw ~ Group + (1|Timepoint), data = DPRC_neuropsych_data, REML=TRUE)
#anova(lmeHayBTime1mod)
#Group lme posthoc
#posthoc_Group_lmeHayBTime1mod <- glht(lmeHayBTime1mod, linfct=mcp(Group ="Tukey"))
#summary(posthoc_Group_lmeHayBTime1mod)
#remove NAs from dataset for given variable
#noNAs_HayBTime1Raw <- DPRC_neuropsych_data[,c("ParticipantID","Group","Timepoint","Age","Sex","HayBTime1.Raw")]
#noNAs_HayBTime1Raw <- noNAs_HayBTime1Raw[complete.cases(noNAs_HayBTime1Raw), ]
#run post hoc w/ Tukey correction
#post_hoc_aov_HayBTime1Raw_mod <- lme(HayBTime1.Raw ~ Group*Timepoint, random = ~1 | ParticipantID/Timepoint, data=noNAs_HayBTime1Raw)
#summary(glht(post_hoc_aov_HayBTime1Raw_mod, linfct=mcp(Timepoint="Tukey")))
#nlme::intervals(post_hoc_aov_HayBTime1Raw_mod, level=0.95) #not working
#summary(glht(post_hoc_aov_HayBTime1Raw_mod, linfct=mcp(Group="Tukey")))
#non-sig. interaction - test by Group (FDR)
# DPRC_neuropsych_data %>%
#   pairwise_t_test(
#     HayBTime1.Raw ~ Group, 
#     p.adjust.method = "fdr"
#   )
#non-sig. interaction & Timepoint, Test by Group (Tukey test)
aov(HayBTime1.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(HayBTime1.Raw~Group)
#plot HayBTime2.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBTime2.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBTime2.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 2 Time Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for HayBTime2.Raw
aov_HayBTime2Raw <- anova_test(data=DPRC_neuropsych_data, dv=HayBTime2.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_HayBTime2Raw)
#check descriptive statistics per each group, per each timepoint
HayBTime2.Raw_descrip <- describeBy(DPRC_neuropsych_data$HayBTime2.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_HayBTime2Raw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_HayBTime2Raw <- na.omit(F0_HayBTime2Raw$HayBTime2.Raw)
mean(noNAsF0_HayBTime2Raw)
sd(noNAsF0_HayBTime2Raw)
#F2
F2_HayBTime2Raw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_HayBTime2Raw <- na.omit(F2_HayBTime2Raw$HayBTime2.Raw)
mean(noNAsF2_HayBTime2Raw)
sd(noNAsF2_HayBTime2Raw)
#Post hoc tests
#non-sig. interaction & Timepoint, Test by Group (Tukey test)
aov(HayBTime2.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(HayBTime2.Raw~Group)
#plot HayBCatA.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBCatA.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBCatA.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 2 Error Rate - Category A (total number of incorrect 'plausible' responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for HayBCatA.Raw
aov_HayBCatARaw <- anova_test(data=DPRC_neuropsych_data, dv=HayBCatA.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_HayBCatARaw)
#check descriptive statistics per each group, per each timepoint
HayBCatA.Raw_descrip <- describeBy(DPRC_neuropsych_data$HayBCatA.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_HayBCatARaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_HayBCatARaw <- na.omit(F0_HayBCatARaw$HayBCatA.Raw)
mean(noNAsF0_HayBCatARaw)
sd(noNAsF0_HayBCatARaw)
#F2
F2_HayBCatARaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_HayBCatARaw <- na.omit(F2_HayBCatARaw$HayBCatA.Raw)
mean(noNAsF2_HayBCatARaw)
sd(noNAsF2_HayBCatARaw)
#Post hoc tests
#non-sig. interaction & Timepoint, Test by Group (Tukey test)
aov(HayBCatA.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(HayBCatA.Raw~Group)

#plot HayBCatB.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = HayBCatB.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = HayBCatB.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Hayling Sentence Set 2 Error Rate - Category B (total number of incorrect 'somewhat plausible' responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for HayBTime1.Raw
aov_HayBCatBRaw <- anova_test(data=DPRC_neuropsych_data, dv=HayBCatB.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_HayBCatBRaw)
#check descriptive statistics per each group, per each timepoint
HayBCatB.Raw_descrip <- describeBy(DPRC_neuropsych_data$HayBCatB.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_HayBCatBRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_HayBCatBRaw <- na.omit(F0_HayBCatBRaw$HayBCatB.Raw)
mean(noNAsF0_HayBCatBRaw)
sd(noNAsF0_HayBCatBRaw)
#F2
F2_HayBCatBRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_HayBCatBRaw <- na.omit(F2_HayBCatBRaw$HayBCatB.Raw)
mean(noNAsF2_HayBCatBRaw)
sd(noNAsF2_HayBCatBRaw)
#Post hoc tests
#non-sig. interaction & Timepoint, Test by Group (Tukey test)
aov(HayBCatB.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(HayBCatB.Raw~Group)

#3.D-KEFS Stroop Task ---------------------------------------------------------#
#plot ColorNaming.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = ColorNaming.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = ColorNaming.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Stroop - Color Naming Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for ColorNaming.Raw
aov_ColorNamingRaw <- anova_test(data=DPRC_neuropsych_data, dv=ColorNaming.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_ColorNamingRaw)
#check descriptive statistics per each group, per each timepoint
ColorNaming.Raw_descrip <- describeBy(DPRC_neuropsych_data$ColorNaming.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_ColorNamingRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_ColorNamingRaw <- na.omit(F0_ColorNamingRaw$ColorNaming.Raw)
mean(noNAsF0_ColorNamingRaw)
sd(noNAsF0_ColorNamingRaw)
#F2
F2_ColorNamingRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_ColorNamingRaw <- na.omit(F2_ColorNamingRaw$ColorNaming.Raw)
mean(noNAsF2_ColorNamingRaw)
sd(noNAsF2_ColorNamingRaw)
#Post hoc tests
#Main Effect Test by Group (Tukey test)
aov(ColorNaming.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(ColorNaming.Raw~Group)
#non-sig. interaction - test by Time Point
# DPRC_neuropsych_data %>%
#   pairwise_t_test(
#     ColorNaming.Raw ~ Timepoint, paired = TRUE,
#     p.adjust.method = "fdr"
#   )
#Main Effect Test by Timepoint (Tukey test)
aov(ColorNaming.Raw ~ Timepoint, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Time Point
DPRC_neuropsych_data%>%cohens_d(ColorNaming.Raw~Timepoint, paired=TRUE)
#Sig. Interaction testing
#Simple main effect (aka simple effects tests) w/ interaction - for Group
# posthoc_ME_Group_ColorNaming <- DPRC_neuropsych_data %>%
#   group_by(Timepoint) %>%
#   anova_test(dv=ColorNaming.Raw,wid=Individual_number,between=Group) %>%
#   tukey_hsd() #doesn't work
# posthoc_ME_Group_ColorNaming ##within each timepoint (F0 & F2), there are Group differences##
posthoc_ME_Group_ColorNaming <- DPRC_neuropsych_data %>%
  group_by(Timepoint) %>%
  anova_test(dv=ColorNaming.Raw,wid=Individual_number,between=Group) %>%
  adjust_pvalue(method="fdr")
posthoc_ME_Group_ColorNaming ##within each timepoint (F0 & F2), there are Group differences##
#Run Pairwise comparison between groups levels if simple main effects (above) is sig.
posthoc_pairwise_Group_ColorNaming <- DPRC_neuropsych_data %>%
  group_by(Timepoint) %>%
  tukey_hsd(ColorNaming.Raw ~ Group)
posthoc_pairwise_Group_ColorNaming ##examines all Group contrasts within each timepoint (F0 & F2) ##
# posthoc_pairwise_Group_ColorNaming <- DPRC_neuropsych_data %>%
#   group_by(Timepoint) %>%
#   pairwise_t_test(ColorNaming.Raw ~ Group, p.adjust.method = "fdr")
# posthoc_pairwise_Group_ColorNaming ##examines all Group contrasts within each timepoint (F0 & F2) ##
#Simple main effect (aka simple effects tests) w/ interaction - for Timepoint
# posthoc_ME_Timepoint_ColorNaming <- DPRC_neuropsych_data %>%
#   group_by(Group) %>%
#   anova_test(dv=ColorNaming.Raw,wid=Individual_number,within=Timepoint,effect.size = "pes") %>%
#   get_anova_table() %>%
#   tukey_hsd() #doesn't work 
# posthoc_ME_Timepoint_ColorNaming ## Across each Group, where were there sig. differences between timepoints (F0 vs. F2)? ##
posthoc_ME_Timepoint_ColorNaming <- DPRC_neuropsych_data %>%
  group_by(Group) %>%
  anova_test(dv=ColorNaming.Raw,wid=Individual_number,within=Timepoint,effect.size = "pes") %>%
  get_anova_table() %>%
  adjust_pvalue(method="fdr")
posthoc_ME_Timepoint_ColorNaming ## Across each Group, where were there sig. differences between timepoints (F0 vs. F2)? ##
#Pairwise comparison between groups levels if simple main effects (above) is sig.
posthoc_pairwise_Timepoint_ColorNaming <- DPRC_neuropsych_data %>%
  group_by(Group) %>%
  tukey_hsd(ColorNaming.Raw ~ Timepoint, paired = TRUE)
posthoc_pairwise_Timepoint_ColorNaming ##examines all Timepoint contrasts between each Group (1,2,3,4,5) ##
# posthoc_pairwise_Timepoint_ColorNaming <- DPRC_neuropsych_data %>%
#   group_by(Group) %>%
#   pairwise_t_test(
#     ColorNaming.Raw ~ Timepoint, paired = TRUE,
#     p.adjust.methods = "fdr") %>%
#   select(-df, -statistic, -p) # Remove details
# posthoc_pairwise_Timepoint_ColorNaming  ##examines all Timepoint contrasts between each Group (1,2,3,4,5) ##


#all pairwise comparisons, Tukey test - from: https://online.stat.psu.edu/stat485/lesson/12/12.7
summary(lm(ColorNaming.Raw~Group+Timepoint+Group:Timepoint,data=DPRC_neuropsych_data)) 
TukeyHSD(aov(ColorNaming.Raw~Group+Timepoint+Group:Timepoint,data=DPRC_neuropsych_data)) 
#effect size for interaction (aMCI Group and Time Point)
DPRC_neuropsych_data_aMCI_ColorNaming_long <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3)
DPRC_neuropsych_data_aMCI_ColorNaming_long$Group <- droplevels(DPRC_neuropsych_data_aMCI_ColorNaming_long$Group)
DPRC_neuropsych_data_aMCI_ColorNaming_long%>%cohens_d(ColorNaming.Raw~Timepoint, paired=TRUE)

#plot interaction 
#remove NAs for ColorNaming variable
noNAs_ColorNamingRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","ColorNaming.Raw")]
noNAs_ColorNamingRaw <- noNAs_ColorNamingRaw[complete.cases(noNAs_ColorNamingRaw), ]
#view interaction plot - by Group
noNAs_ColorNamingRaw%>%
  group_by(Group,Timepoint) %>%
  summarise(s_mean=mean(ColorNaming.Raw)) %>%
  ggplot(aes(y=s_mean,x=Group,colour=Timepoint,group=Timepoint))+
  geom_point()+geom_line()+
  scale_x_discrete(labels = c("1" = "C", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  ylab("Stroop - Colour Naming Raw (seconds)") +
  theme_classic()
#view interaction plot - by Timepoint
noNAs_ColorNamingRaw %>%
  group_by(Group,Timepoint) %>%
  summarise(s_mean=mean(ColorNaming.Raw)) %>%
  ggplot(aes(y=s_mean,x=Timepoint,colour=Group,group=Group))+
  geom_point()+geom_line()+
  scale_color_discrete(labels = c("1" = "C", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  ylab("Stroop - Colour Naming Raw (seconds)") +
  theme_classic()

#test if there is a difference among differences (e.g., is the aMCI group difference for color naming between the time points 
#(F0 - F2) sig. greater compared to the other groups?)
#subtract differences in F2 - F0 for color naming
ColorNaming_diff_data<-DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","ColorNaming.Raw")]
F0_ColorNaming_data<-ColorNaming_diff_data[ColorNaming_diff_data$Timepoint=='F0', ]
F2_ColorNaming_data<-ColorNaming_diff_data[ColorNaming_diff_data$Timepoint=='F2', ]
ColorNamingDiff<-F2_ColorNaming_data$ColorNaming.Raw - F0_ColorNaming_data$ColorNaming.Raw
F0_ColorNaming_data$ColorNamingDiff<-ColorNamingDiff
#run one-way ANOVA on the F2-F0 differences in color naming between groups
ColorNamingDiff_mod <- lm(ColorNamingDiff ~ Group, data = F0_ColorNaming_data)
anova(ColorNamingDiff_mod) #sig. difference
#effect size omnibus ANOVA
etaSquared(ColorNamingDiff_mod)
#post hoc f/u test
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_ColorNamingDiff_mod <- glht(ColorNamingDiff_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_ColorNamingDiff_mod)
confint(post_hoc_ColorNamingDiff_mod)
#effect size for sig. post hoc tests
  #for SCD vs. aMCI 
  SCDvaMCI_ColorNamingDiff <- subset(F0_ColorNaming_data, F0_ColorNaming_data$Group == 2 | F0_ColorNaming_data$Group == 3)
  SCDvaMCI_ColorNamingDiff$Group <- droplevels(SCDvaMCI_ColorNamingDiff$Group)
  cohensD(ColorNamingDiff ~ Group, data = SCDvaMCI_ColorNamingDiff) #this looks like Hedges' g? 
#plot Color Naming differences (raincloud plot)
ggplot(F0_ColorNaming_data, aes(x = Group, y = ColorNamingDiff, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = ColorNamingDiff, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Color Naming Difference (F2 - F0)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none")+
  coord_flip()

#plot WordReading.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = WordReading.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = WordReading.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Stroop - Word Reading Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for WordReading.Raw
aov_WordReadingRaw <- anova_test(data=DPRC_neuropsych_data, dv=WordReading.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_WordReadingRaw)
#check descriptive statistics per each group, per each timepoint
WordReading.Raw_descrip <- describeBy(DPRC_neuropsych_data$WordReading.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_WordReadingRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_WordReadingRaw <- na.omit(F0_WordReadingRaw$WordReading.Raw)
mean(noNAsF0_WordReadingRaw)
sd(noNAsF0_WordReadingRaw)
#F2
F2_WordReadingRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_WordReadingRaw <- na.omit(F2_WordReadingRaw$WordReading.Raw)
mean(noNAsF2_WordReadingRaw)
sd(noNAsF2_WordReadingRaw)
#Post hoc tests
#non-sig. interaction - test by Group
aov(WordReading.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(WordReading.Raw~Group)
#non-sig. interaction - test by Time Point
aov(WordReading.Raw ~ Timepoint, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Time Point
DPRC_neuropsych_data%>%cohens_d(WordReading.Raw~Timepoint, paired=TRUE)

#plot Inhibition.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Inhibition.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Inhibition.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Stroop - Inhibition Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for Inhibition.Raw
aov_InhibitionRaw <- anova_test(data=DPRC_neuropsych_data, dv=Inhibition.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_InhibitionRaw)
#check descriptive statistics per each group, per each timepoint
Inhibition.Raw_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_InhibitionRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_InhibitionRaw <- na.omit(F0_InhibitionRaw$Inhibition.Raw)
mean(noNAsF0_InhibitionRaw)
sd(noNAsF0_InhibitionRaw)
#F2
F2_InhibitionRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_InhibitionRaw <- na.omit(F2_InhibitionRaw$Inhibition.Raw)
mean(noNAsF2_InhibitionRaw)
sd(noNAsF2_InhibitionRaw)
#Post hoc tests
#non-sig. interaction - test by Group
aov(Inhibition.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(Inhibition.Raw~Group)
#non-sig. interaction - test by Time Point
aov(Inhibition.Raw ~ Timepoint, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Time Point
DPRC_neuropsych_data%>%cohens_d(Inhibition.Raw~Timepoint, paired=TRUE)

#plot Inhibition.Colour.Naming (Inhibition / Colour Naming) (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Inhibition.Colour.Naming, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Inhibition.Colour.Naming, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Stroop - Inhibition / Colour Naming") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for Inhibition.Colour.Naming
aov_InhibitionColourNaming <- anova_test(data=DPRC_neuropsych_data, dv=Inhibition.Colour.Naming, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_InhibitionColourNaming)
#check descriptive statistics per each group, per each timepoint
Inhibition.Colour.Naming_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Colour.Naming, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_InhibitionColourNaming <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_InhibitionColourNaming <- na.omit(F0_InhibitionColourNaming$Inhibition.Colour.Naming)
mean(noNAsF0_InhibitionColourNaming)
sd(noNAsF0_InhibitionColourNaming)
#F2
F2_InhibitionColourNaming <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_InhibitionColourNaming <- na.omit(F2_InhibitionColourNaming$Inhibition.Colour.Naming)
mean(noNAsF2_InhibitionColourNaming)
sd(noNAsF2_InhibitionColourNaming)

#plot Inhibition.Word.Reading (Inhibition / Word Reading) (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Inhibition.Word.Reading, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Inhibition.Word.Reading, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Stroop - Inhibition / Word Reading") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for Inhibition.Word.Reading
aov_InhibitionWordReading <- anova_test(data=DPRC_neuropsych_data, dv=Inhibition.Word.Reading, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_InhibitionWordReading)
#check descriptive statistics per each group, per each timepoint
Inhibition.Word.Reading_descrip <- describeBy(DPRC_neuropsych_data$Inhibition.Word.Reading, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_InhibitionWordReading <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_InhibitionWordReading <- na.omit(F0_InhibitionWordReading$Inhibition.Word.Reading)
mean(noNAsF0_InhibitionWordReading)
sd(noNAsF0_InhibitionWordReading)
#F2
F2_InhibitionWordReading <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_InhibitionWordReading <- na.omit(F2_InhibitionWordReading$Inhibition.Word.Reading)
mean(noNAsF2_InhibitionWordReading)
sd(noNAsF2_InhibitionWordReading)

#4.D-KEFS Verbal Fluency + Category Fluency Task ------------------------------#
#plot LetFluency.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = LetFluency.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = LetFluency.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Letter Fluency (total correct responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for LetFluency.Raw
aov_LetFluencyRaw <- anova_test(data=DPRC_neuropsych_data, dv=LetFluency.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_LetFluencyRaw)
#check descriptive statistics per each group, per each timepoint
LetFluency.Raw_descrip <- describeBy(DPRC_neuropsych_data$LetFluency.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_LetFluencyRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_LetFluencyRaw <- na.omit(F0_LetFluencyRaw$LetFluency.Raw)
mean(noNAsF0_LetFluencyRaw)
sd(noNAsF0_LetFluencyRaw)
#F2
F2_LetFluencyRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_LetFluencyRaw <- na.omit(F2_LetFluencyRaw$LetFluency.Raw)
mean(noNAsF2_LetFluencyRaw)
sd(noNAsF2_LetFluencyRaw)
#Post hoc tests
#non-sig. interaction - test by Group
aov(LetFluency.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(LetFluency.Raw~Group)
#plot CatFluency.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = CatFluency.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = CatFluency.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Category Fluency (total correct responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for CatFluency.Raw
aov_CatFluencyRaw <- anova_test(data=DPRC_neuropsych_data, dv=CatFluency.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_CatFluencyRaw)
#check descriptive statistics per each group, per each timepoint
CatFluency.Raw_descrip <- describeBy(DPRC_neuropsych_data$CatFluency.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_CatFluencyRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_CatFluencyRaw <- na.omit(F0_CatFluencyRaw$CatFluency.Raw)
mean(noNAsF0_CatFluencyRaw)
sd(noNAsF0_CatFluencyRaw)
#F2
F2_CatFluencyRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_CatFluencyRaw <- na.omit(F2_CatFluencyRaw$CatFluency.Raw)
mean(noNAsF2_CatFluencyRaw)
sd(noNAsF2_CatFluencyRaw)
#Post hoc tests
#non-sig. interaction - test by Group
aov(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(CatFluency.Raw~Group)

#plot Switching.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = Switching.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Switching.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Category Switching (total correct responses)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for CatFluency.Raw
aov_SwitchingRaw <- anova_test(data=DPRC_neuropsych_data, dv=Switching.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_SwitchingRaw)
#check descriptive statistics per each group, per each timepoint
Switching.Raw_descrip <- describeBy(DPRC_neuropsych_data$Switching.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_SwitchingRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_SwitchingRaw <- na.omit(F0_SwitchingRaw$Switching.Raw)
mean(noNAsF0_SwitchingRaw)
sd(noNAsF0_SwitchingRaw)
#F2
F2_SwitchingRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_SwitchingRaw <- na.omit(F2_SwitchingRaw$Switching.Raw)
mean(noNAsF2_SwitchingRaw)
sd(noNAsF2_SwitchingRaw)
#Post hoc tests
#non-sig. interaction - test by Group
aov(Switching.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(Switching.Raw~Group)
#non-sig. interaction - test by Time Point
aov(Switching.Raw ~ Timepoint, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Time Point
DPRC_neuropsych_data%>%cohens_d(Switching.Raw~Timepoint, paired=TRUE)

#5.Trail Making Test (TMT) A & B ----------------------------------------------#
#plot TrailsA.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TrailsA.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TrailsA.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Trail Making Test A Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for TrailsA.Raw
aov_TrailsARaw <- anova_test(data=DPRC_neuropsych_data, dv=TrailsA.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_TrailsARaw)
#check descriptive statistics per each group, per each timepoint
TrailsA.Raw_descrip <- describeBy(DPRC_neuropsych_data$TrailsA.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_TrailsARaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_TrailsARaw <- na.omit(F0_TrailsARaw$TrailsA.Raw)
mean(noNAsF0_TrailsARaw)
sd(noNAsF0_TrailsARaw)
#F2
F2_TrailsARaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_TrailsARaw <- na.omit(F2_TrailsARaw$TrailsA.Raw)
mean(noNAsF2_TrailsARaw)
sd(noNAsF2_TrailsARaw)
#Post hoc tests
#non-sig. interaction - test by Group
aov(TrailsA.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(TrailsA.Raw~Group)

#plot TrailsB.Raw (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TrailsB.Raw, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TrailsB.Raw, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("Trail Making Test B Raw (seconds)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for TrailsB.Raw
aov_TrailsBRaw <- anova_test(data=DPRC_neuropsych_data, dv=TrailsB.Raw, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_TrailsBRaw)
#check descriptive statistics per each group, per each timepoint
TrailsB.Raw_descrip <- describeBy(DPRC_neuropsych_data$TrailsB.Raw, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_TrailsBRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_TrailsBRaw <- na.omit(F0_TrailsBRaw$TrailsB.Raw)
mean(noNAsF0_TrailsBRaw)
sd(noNAsF0_TrailsBRaw)
#F2
F2_TrailsBRaw <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_TrailsBRaw <- na.omit(F2_TrailsBRaw$TrailsB.Raw)
mean(noNAsF2_TrailsBRaw)
sd(noNAsF2_TrailsBRaw)
#Post hoc tests
#Test by Group
aov(TrailsB.Raw ~ Group, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Groups
DPRC_neuropsych_data%>%cohens_d(TrailsB.Raw~Group)
#Test by Time Point
aov(TrailsB.Raw ~ Timepoint, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Time Point
DPRC_neuropsych_data%>%cohens_d(TrailsB.Raw~Timepoint, paired=TRUE)
#Sig. interaction
#Simple main effect w/ interaction - for Group
posthoc_ME_Group_TrailsB <- DPRC_neuropsych_data %>%
  group_by(Timepoint) %>%
  anova_test(dv=TrailsB.Raw,wid=Individual_number,between=Group) %>%
  adjust_pvalue(method="fdr")
posthoc_ME_Group_TrailsB ##within each timepoint (F0 & F2), there are Group differences##
#Run Pairwise comparison between groups levels if simple main effects (above) is sig.
posthoc_pairwise_Group_TrailsB <- DPRC_neuropsych_data %>%
  group_by(Timepoint) %>%
  tukey_hsd(TrailsB.Raw ~ Group)
posthoc_pairwise_Group_TrailsB ##examines all Group contrasts within each timepoint (F0 & F2) ##
#Simple main effect w/ ineraction - for Timepoint
posthoc_ME_Timepoint_TrailsB <- DPRC_neuropsych_data %>%
  group_by(Group) %>%
  anova_test(dv=TrailsB.Raw,wid=Individual_number,within=Timepoint,effect.size = "pes") %>%
  get_anova_table() %>%
  adjust_pvalue(method="fdr")
posthoc_ME_Timepoint_TrailsB ## Across each Group, where were there sig. differences between timepoints (F0 vs. F2)? ##
#Pairwise comparison between groups levels if simple main effects (above) is sig.
posthoc_pairwise_Timepoint_TrailsB <- DPRC_neuropsych_data %>%
  group_by(Group) %>%
  tukey_hsd(TrailsB.Raw ~ Timepoint, paired = TRUE)
posthoc_pairwise_Timepoint_TrailsB ##examines all Timepoint contrasts between each Group (1,2,3,4,5) ##
#All pairwise comparisons with Tukey test
summary(lm(TrailsB.Raw~Group+Timepoint+Group:Timepoint,data=DPRC_neuropsych_data))
TukeyHSD(aov(TrailsB.Raw~Group+Timepoint+Group:Timepoint,data=DPRC_neuropsych_data))
#effect size for interaction (mMCI Group and Time Point)
DPRC_neuropsych_data_mMCI_TrailsB_long <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 4)
DPRC_neuropsych_data_mMCI_TrailsB_long$Group <- droplevels(DPRC_neuropsych_data_mMCI_TrailsB_long$Group)
DPRC_neuropsych_data_mMCI_TrailsB_long%>%cohens_d(TrailsB.Raw~Timepoint, paired=TRUE)
#plot interaction 
#remove NAs for Trails B variable
noNAs_TrailsBRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","TrailsB.Raw")]
noNAs_TrailsBRaw<- noNAs_TrailsBRaw[complete.cases(noNAs_TrailsBRaw), ]
#view interaction plot - by Group
noNAs_TrailsBRaw%>%
  group_by(Group,Timepoint) %>%
  summarise(s_mean=mean(TrailsB.Raw)) %>%
  ggplot(aes(y=s_mean,x=Group,colour=Timepoint,group=Timepoint))+
  geom_point()+geom_line()+
  scale_x_discrete(labels = c("1" = "C", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  ylab("Trail Making Test B Raw (seconds)") +
  theme_classic()
#view interaction plot - by Timepoint
noNAs_TrailsBRaw %>%
  group_by(Group,Timepoint) %>%
  summarise(s_mean=mean(TrailsB.Raw)) %>%
  ggplot(aes(y=s_mean,x=Timepoint,colour=Group,group=Group))+
  geom_point()+geom_line()+
  scale_color_discrete(labels = c("1" = "C", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  ylab("Trail Making Test B Raw (seconds)") +
  theme_classic()

#test if there is a difference among differences (e.g., is the aMCI group difference for TMT-B between the time points 
#(F0 - F2) sig. greater compared to the other groups?)
#subtract differences in F2 - F0 for TMT-B
TrailsB_diff_data<-DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","TrailsB.Raw")]
F0_TrailsB_data<-TrailsB_diff_data[TrailsB_diff_data$Timepoint=='F0', ]
F2_TrailsB_data<-TrailsB_diff_data[TrailsB_diff_data$Timepoint=='F2', ]
TrailsBDiff<-F2_TrailsB_data$TrailsB.Raw - F0_TrailsB_data$TrailsB.Raw
F0_TrailsB_data$TrailsBDiff<-TrailsBDiff
#run one-way ANOVA on the F2-F0 differences in color naming between groups
TrailsBDiff_mod <- lm(TrailsBDiff ~ Group, data = F0_TrailsB_data)
anova(TrailsBDiff_mod) #sig. difference
#effect size omnibus ANOVA
etaSquared(TrailsBDiff_mod)
#post hoc f/u test
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsBDiff_mod <- glht(TrailsBDiff_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsBDiff_mod)
confint(post_hoc_TrailsBDiff_mod)
#effect size for sig. post hoc tests
#for SCD vs. aMCI 
SCDvaMCI_TrailsBDiff <- subset(F0_TrailsB_data, F0_TrailsB_data$Group == 2 | F0_TrailsB_data$Group == 3)
SCDvaMCI_TrailsBDiff$Group <- droplevels(SCDvaMCI_TrailsBDiff$Group)
cohensD(TrailsBDiff ~ Group, data = SCDvaMCI_TrailsBDiff) #this looks like Hedges' g? 
#plot Color Naming differences (raincloud plot)
ggplot(F0_TrailsB_data, aes(x = Group, y = TrailsBDiff, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TrailsBDiff, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Trails B Difference (F2 - F0)") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none")+
  coord_flip()

#plot TMT.B.TMT.A (TMT-B / TMT-A) (raincloud plot)
ggplot(DPRC_neuropsych_data, aes(x = Group, y = TMT.B.TMT.A, fill = Timepoint)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = TMT.B.TMT.A, color = Timepoint), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Timepoint)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Timepoint)) + 
  xlab("Group") + 
  ylab("TMT-B / TMT-A") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  coord_flip()
#run mixed design, 2 x 5 ANOVA for TrailsB.Raw
aov_TMTBTMTA <- anova_test(data=DPRC_neuropsych_data, dv=TMT.B.TMT.A, wid=Individual_number, between=Group, within=Timepoint, effect.size = "pes")
get_anova_table(aov_TMTBTMTA)
#check descriptive statistics per each group, per each timepoint
TMT.B.TMT.A_descrip <- describeBy(DPRC_neuropsych_data$TMT.B.TMT.A, list(DPRC_neuropsych_data$Group, DPRC_neuropsych_data$Timepoint))
#find mean & SD from total sample in F0 and F2 timepoints:
F0_TMTBTMTA <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F0",]
noNAsF0_TMTBTMTA <- na.omit(F0_TMTBTMTA$TMT.B.TMT.A)
mean(noNAsF0_TMTBTMTA)
sd(noNAsF0_TMTBTMTA)
#F2
F2_TMTBTMTA <- DPRC_neuropsych_data[DPRC_neuropsych_data[, "Timepoint"] == "F2",]
noNAsF2_TMTBTMTA <- na.omit(F2_TMTBTMTA$TMT.B.TMT.A)
mean(noNAsF2_TMTBTMTA)
sd(noNAsF2_TMTBTMTA)
#Post hoc
#non-sig. interaction - test by Time Point
aov(TMT.B.TMT.A ~ Timepoint, data = DPRC_neuropsych_data) %>% tukey_hsd()
#effect size for Time Point
DPRC_neuropsych_data%>%cohens_d(TMT.B.TMT.A~Timepoint, paired=TRUE)

############-----------ANCOVA Testing (w/ age) -----------------------------------------####### 
#run mixed design, 2 x 5 ANCOVA for HayBTime1.Raw
aov_HayBTime1Raw <- aov(HayBTime1.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_HayBTime1Raw)
# aov_HayBTime1Raw <- anova_test(data=DPRC_neuropsych_data, dv=HayBTime1.Raw, wid=Individual_number, between=Group, within=Timepoint, covariate=Age, effect.size = "pes") #will not accept 'Age' as a covariate
# get_anova_table(aov_HayBTime1Raw)
#aov_HayBTime1Raw <- Anova(lm(HayBTime1.Raw ~ Group*Timepoint + Age  data=DPRC_neuropsych_data), type = "III") #use this anova test to account for unbalanced designs/sample sizes
#aov_HayBTime1Raw
#effect size w/ covariates
#eta_squared(aov_HayBTime1Raw) #doesn't work - doesn't like the class of the model
sjstats::eta_sq(aov_HayBTime1Raw)
#run mixed design, 2 x 5 ANCOVA for HayBTime2.Raw
aov_HayBTime2Raw <- aov(HayBTime2.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_HayBTime2Raw)
sjstats::eta_sq(aov_HayBTime2Raw)
#run mixed design, 2 x 5 ANCOVA for HayBCatA.Raw
aov_HayBCatARaw <- aov(HayBCatA.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_HayBCatARaw)
sjstats::eta_sq(aov_HayBCatARaw)
#run mixed design, 2 x 5 ANCOVA for HayBCatB.Raw
aov_HayBCatBRaw <- aov(HayBCatB.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_HayBCatBRaw)
sjstats::eta_sq(aov_HayBCatBRaw)
#run mixed design, 2 x 5 ANCOVA for ColorNaming.Raw
aov_ColorNamingRaw <- aov(ColorNaming.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_ColorNamingRaw)
sjstats::eta_sq(aov_ColorNamingRaw)
#run mixed design, 2 x 5 ANCOVA for WordReading.Raw
aov_WordReadingRaw <- aov(WordReading.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_WordReadingRaw)
sjstats::eta_sq(aov_WordReadingRaw)
#run mixed design, 2 x 5 ANCOVA for WordReading.Raw
aov_InhibitionRaw <- aov(Inhibition.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_InhibitionRaw)
sjstats::eta_sq(aov_InhibitionRaw)
#run mixed design, 2 x 5 ANCOVA for Inhibition.Colour.Naming
aov_InhibitionColourNaming <- aov(Inhibition.Colour.Naming~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_InhibitionColourNaming)
#run mixed design, 2 x 5 ANCOVA for Inhibition.Word.Reading
aov_InhibitionWordReading <- aov(Inhibition.Word.Reading~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_InhibitionWordReading)
sjstats::eta_sq(aov_InhibitionWordReading)
#run mixed design, 2 x 5 ANCOVA for TrailsA.Raw
aov_TrailsARaw <- aov(TrailsA.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_TrailsARaw)
sjstats::eta_sq(aov_TrailsARaw)
#run mixed design, 2 x 5 ANCOVA for TrailsB.Raw
aov_TrailsBRaw <- aov(TrailsB.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_TrailsBRaw)
sjstats::eta_sq(aov_TrailsBRaw)
#run mixed design, 2 x 5 ANCOVA for TMT.B.TMT.A
aov_TMTBTMTA <- aov(TMT.B.TMT.A ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_TMTBTMTA)
sjstats::eta_sq(aov_TMTBTMTA)
#run mixed design, 2 x 5 ANCOVA for LetFluency.Raw
aov_LetFluencyRaw <- aov(LetFluency.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_LetFluencyRaw)
sjstats::eta_sq(aov_LetFluencyRaw)
#run mixed design, 2 x 5 ANCOVA for CatFluency.Raw
aov_CatFluencyRaw <- aov(CatFluency.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_CatFluencyRaw)
sjstats::eta_sq(aov_CatFluencyRaw)
#run mixed design, 2 x 5 ANCOVA for Switching.Raw
aov_SwitchingRaw <- aov(Switching.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age,  data = DPRC_neuropsych_data)
summary(aov_SwitchingRaw)
sjstats::eta_sq(aov_SwitchingRaw)

#f/u post hoc tests
#HayBTime1
#remove NAs from dataset for given variable
noNAs_HayBTime1Raw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","HayBTime1.Raw")]
noNAs_HayBTime1Raw <- noNAs_HayBTime1Raw[complete.cases(noNAs_HayBTime1Raw), ]
post_hoc_aov_HayBTime1_covar_mod <- lme(HayBTime1.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_HayBTime1Raw)
summary(glht(post_hoc_aov_HayBTime1_covar_mod, linfct=mcp(Group="Tukey"))) #no sig.
#HayBTime2
#remove NAs from dataset for given variable
noNAs_HayBTime2Raw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","HayBTime2.Raw")]
noNAs_HayBTime2Raw <- noNAs_HayBTime2Raw[complete.cases(noNAs_HayBTime2Raw), ]
post_hoc_aov_HayBTime2_covar_mod <- lme(HayBTime2.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_HayBTime2Raw)
summary(glht(post_hoc_aov_HayBTime2_covar_mod, linfct=mcp(Group="Tukey"))) #no sig.
#HayBCatA
#remove NAs from dataset for given variable
noNAs_HayBCatARaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","HayBCatA.Raw")]
noNAs_HayBCatARaw <- noNAs_HayBCatARaw[complete.cases(noNAs_HayBCatARaw), ]
post_hoc_aov_HayBCatA_covar_mod <- lme(HayBCatA.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_HayBCatARaw)
t_value_effect_size <- summary(glht(post_hoc_aov_HayBCatA_covar_mod, linfct=mcp(Group="Tukey")))
summary(glht(post_hoc_aov_HayBCatA_covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate  
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_HayBCatA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_HayBCatA.Raw$HayBCatA.Raw,DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$HayBCatA.Raw,DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$HayBCatA.Raw,DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for mMCI vs. AD
  DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 4 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$HayBCatA.Raw,DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 4'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#HayBCatB
#remove NAs from dataset for given variable
noNAs_HayBCatBRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","HayBCatB.Raw")]
noNAs_HayBCatBRaw <- noNAs_HayBCatBRaw[complete.cases(noNAs_HayBCatBRaw), ]
post_hoc_aov_HayBCatB_covar_mod <- lme(HayBCatB.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_HayBCatBRaw)
t_value_effect_size <- summary(glht(post_hoc_aov_HayBCatB_covar_mod, linfct=mcp(Group="Tukey")))
summary(glht(post_hoc_aov_HayBCatB_covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate  
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_HayBCatB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_HayBCatB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_HayBCatB.Raw$HayBCatB.Raw,DPRC_neuropsych_data_CvAD_HayBCatB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$HayBCatB.Raw,DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw$HayBCatB.Raw,DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for mMCI vs. AD
  DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 4 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw$HayBCatB.Raw,DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 4'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#ColorNaming
#remove NAs from dataset for given variable
noNAs_ColorNamingRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","ColorNaming.Raw")]
noNAs_ColorNamingRaw <- noNAs_ColorNamingRaw[complete.cases(noNAs_ColorNamingRaw), ]
post_hoc_aov_ColorNaming_covar_mod <- lme(ColorNaming.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_ColorNamingRaw)
summary(post_hoc_aov_ColorNaming_covar_mod)
summary(glht(post_hoc_aov_ColorNaming_covar_mod, linfct=mcp(Group="Tukey")))
summary(glht(post_hoc_aov_ColorNaming_covar_mod, linfct=mcp(Timepoint="Tukey")))
#interaction - ColorNaming
#Tukey test
TukeyHSD(aov(ColorNaming.Raw~Group*Timepoint+as.factor(Age), data = DPRC_neuropsych_data)) 
# #Simple main effect w/ interaction - for Group
# posthoc_ME_Group_ColorNaming_Age <- DPRC_neuropsych_data %>%
#   group_by(Timepoint) %>%
#   anova_test(dv=TrailsB.Raw,wid=Individual_number,between=Group,covariate=Age) %>%
#   adjust_pvalue(method="fdr")
# posthoc_ME_Group_ColorNaming_Age  ##within each timepoint (F0 & F2), there are Group differences##
# #Run Pairwise comparison between groups levels if simple main effects (above) is sig.
# posthoc_pairwise_Group_ColorNaming_Age <- DPRC_neuropsych_data %>%
#   group_by(Timepoint) %>%
#   tukey_hsd(ColorNaming.Raw ~ Group+as.factor(Age))
# posthoc_pairwise_Group_ColorNaming_Age  ##examines all Group contrasts within each timepoint (F0 & F2) ##
# #Simple main effect w/ ineraction - for Timepoint
# posthoc_ME_Timepoint_ColorNaming_Age <- DPRC_neuropsych_data %>%
#   group_by(Group) %>%
#   anova_test(dv=ColorNaming.Raw,wid=Individual_number,within=Timepoint,covariate=Age,effect.size = "pes") %>%
#   get_anova_table() %>%
#   adjust_pvalue(method="fdr")
# posthoc_ME_Timepoint_ColorNaming_Age ## Across each Group, where were there sig. differences between timepoints (F0 vs. F2)? ##
# #Pairwise comparison between groups levels if simple main effects (above) is sig.
# posthoc_pairwise_Timepoint_ColorNaming_Age <- DPRC_neuropsych_data %>%
#   group_by(Group) %>%
#   tukey_hsd(ColorNaming.Raw ~ Timepoint+as.factor(Age), paired = TRUE)
# posthoc_pairwise_Timepoint_ColorNaming_Age ##examines all Timepoint contrasts between each Group (1,2,3,4,5) ##
#test if there is a difference among differences (e.g., is the aMCI group difference for color naming between the time points 
#(F0 - F2) sig. greater compared to the other groups?)

#subtract differences in F2 - F0 for color naming with covariate age
ColorNaming_diff_data<-DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","ColorNaming.Raw")]
F0_ColorNaming_data<-ColorNaming_diff_data[ColorNaming_diff_data$Timepoint=='F0', ]
F2_ColorNaming_data<-ColorNaming_diff_data[ColorNaming_diff_data$Timepoint=='F2', ]
ColorNamingDiff<-F2_ColorNaming_data$ColorNaming.Raw - F0_ColorNaming_data$ColorNaming.Raw
F0_ColorNaming_data$ColorNamingDiff<-ColorNamingDiff
#run one-way ANOVA on the F2-F0 differences in color naming between groups
ColorNamingDiff_covar_mod <- lm(ColorNamingDiff ~ Group+Age, data = F0_ColorNaming_data)
anova(ColorNamingDiff_covar_mod) #sig. difference
#effect size omnibus ANOVA
etaSquared(ColorNamingDiff_covar_mod)
#post hoc f/u test
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_ColorNamingDiff_covar_mod <- glht(ColorNamingDiff_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_ColorNamingDiff_covar_mod)
confint(post_hoc_ColorNamingDiff_covar_mod)
t_value_effect_size <- summary(glht(ColorNamingDiff_covar_mod, linfct = mcp(Group = "Tukey")))
#effect size for sig. post hoc tests with covariate
  #for SCD vs. aMCI
  DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff$ColorNaming.Raw,DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#WordReading
#remove NAs from dataset for given variable
noNAs_WordReadingRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","WordReading.Raw")]
noNAs_WordReadingRaw <- noNAs_WordReadingRaw[complete.cases(noNAs_WordReadingRaw), ]
post_hoc_aov_WordReading_covar_mod <- lme(WordReading.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_WordReadingRaw)
summary(glht(post_hoc_aov_WordReading_covar_mod, linfct=mcp(Group="Tukey"))) #no sig 
summary(glht(post_hoc_aov_WordReading_covar_mod, linfct=mcp(Timepoint="Tukey")))
#Inhibition
#remove NAs from dataset for given variable
noNAs_InhibitionRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","Inhibition.Raw")]
noNAs_InhibitionRaw <- noNAs_InhibitionRaw[complete.cases(noNAs_InhibitionRaw), ]
post_hoc_aov_Inhibition_covar_mod <- lme(Inhibition.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_InhibitionRaw)
summary(glht(post_hoc_aov_Inhibition_covar_mod, linfct=mcp(Group="Tukey")))
summary(glht(post_hoc_aov_Inhibition_covar_mod, linfct=mcp(Timepoint="Tukey")))
#Inhibition.Word.Reading
#remove NAs from dataset for given variable
noNAs_InhibitionWordReading <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","Inhibition.Word.Reading")]
noNAs_InhibitionWordReading <- noNAs_InhibitionWordReading[complete.cases(noNAs_InhibitionWordReading), ]
post_hoc_aov_InhibitionWordReading_covar_mod <- lme(Inhibition.Word.Reading ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_InhibitionWordReading)
summary(glht(post_hoc_aov_InhibitionWordReading_covar_mod, linfct=mcp(Group="Tukey")))
#TrailsA
#remove NAs from dataset for given variable
noNAs_TrailsARaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","TrailsA.Raw")]
noNAs_TrailsARaw <- noNAs_TrailsARaw[complete.cases(noNAs_TrailsARaw), ]
post_hoc_aov_TrailsA_covar_mod <- lme(TrailsA.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_TrailsARaw)
summary(glht(post_hoc_aov_TrailsA_covar_mod, linfct=mcp(Group="Tukey")))
#TrailsB
#remove NAs from dataset for given variable
noNAs_TrailsBRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","TrailsB.Raw")]
noNAs_TrailsBRaw <- noNAs_TrailsBRaw[complete.cases(noNAs_TrailsBRaw), ]
post_hoc_aov_TrailsB_covar_mod <- lme(TrailsB.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_TrailsBRaw)
summary(glht(post_hoc_aov_TrailsB_covar_mod, linfct=mcp(Group="Tukey")))
summary(glht(post_hoc_aov_TrailsB_covar_mod, linfct=mcp(Timepoint="Tukey")))
#interaction - TrailsB
#Tukey test
TukeyHSD(aov(TrailsB.Raw~Group*Timepoint+as.factor(Age), data = DPRC_neuropsych_data)) 
#subtract differences in F2 - F0 for color naming with covariate age
TrailsB_diff_data<-DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","TrailsB.Raw")]
F0_TrailsB_data<-TrailsB_diff_data[TrailsB_diff_data$Timepoint=='F0', ]
F2_TrailsB_data<-TrailsB_diff_data[TrailsB_diff_data$Timepoint=='F2', ]
TrailsBDiff<-F2_TrailsB_data$TrailsB.Raw - F0_TrailsB_data$TrailsB.Raw
F0_TrailsB_data$TrailsBDiff<-TrailsBDiff
#run one-way ANOVA on the F2-F0 differences in color naming between groups
TrailsBDiff_covar_mod <- lm(TrailsBDiff ~ Group+Age, data = F0_TrailsB_data)
anova(TrailsBDiff_covar_mod) #sig. difference
#effect size omnibus ANOVA
etaSquared(TrailsBDiff_covar_mod)
#post hoc f/u test
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsBDiff_covar_mod <- glht(TrailsBDiff_covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsBDiff_covar_mod)
confint(post_hoc_TrailsBDiff_covar_mod)
t_value_effect_size <- summary(glht(TrailsBDiff_covar_mod, linfct = mcp(Group = "Tukey")))
#effect size for sig. post hoc tests with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_TrailsBDiff <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TrailsBDiff$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TrailsBDiff$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_TrailsBDiff, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_TrailsBDiff$TrailsB.Raw,DPRC_neuropsych_data_CvmMCI_TrailsBDiff$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff$TrailsB.Raw,DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#TrailsB/TrailsA
noNAs_TMTBTMTA <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","TMT.B.TMT.A")]
noNAs_TMTBTMTA <- noNAs_TMTBTMTA[complete.cases(noNAs_TMTBTMTA), ]
post_hoc_aov_TMTBTMTA_covar_mod <- lme(TMT.B.TMT.A ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_TMTBTMTA)
summary(glht(post_hoc_aov_TMTBTMTA_covar_mod, linfct=mcp(Timepoint="Tukey")))
#LetFluency
#remove NAs from dataset for given variable
noNAs_LetFluencyRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","LetFluency.Raw")]
noNAs_LetFluencyRaw <- noNAs_LetFluencyRaw[complete.cases(noNAs_LetFluencyRaw), ]
post_hoc_aov_LetFluency_covar_mod <- lme(LetFluency.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_LetFluencyRaw)
summary(glht(post_hoc_aov_LetFluency_covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size<-summary(glht(post_hoc_aov_LetFluency_covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect size for the significant Group levels
  #effect size for sig. post hoc tests - need to use a.tes (from the compute.es package) to account for covariates
  #for C vs. SCD 
  DPRC_neuropsych_data_CvSCD_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 2)
  DPRC_neuropsych_data_CvSCD_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvSCD_LetFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvSCD_LetFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvSCD_LetFluency.Raw$LetFluency.Raw,DPRC_neuropsych_data_CvSCD_LetFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['2 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. aMCI 
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$LetFluency.Raw,DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#CatFluency
#remove NAs from dataset for given variable
noNAs_CatFluencyRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","CatFluency.Raw")]
noNAs_CatFluencyRaw <- noNAs_CatFluencyRaw[complete.cases(noNAs_CatFluencyRaw), ]
post_hoc_aov_CatFluency_covar_mod <- lme(CatFluency.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_CatFluencyRaw)
summary(glht(post_hoc_aov_CatFluency_covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size<-summary(glht(post_hoc_aov_CatFluency_covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect size for the significant Group levels
  #effect size for sig. post hoc tests - need to use a.tes (from the compute.es package) to account for covariates
  #for C vs. aMCI 
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_CvAD_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. aMCI 
  DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$CatFluency.Raw,DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
#Switching
#remove NAs from dataset for given variable
noNAs_SwitchingRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","Switching.Raw")]
noNAs_SwitchingRaw <- noNAs_SwitchingRaw[complete.cases(noNAs_SwitchingRaw), ]
post_hoc_aov_Switching_covar_mod <- lme(Switching.Raw ~ Group*Timepoint + Age,  random = ~1 | Individual_number/Timepoint, data=noNAs_SwitchingRaw)
summary(glht(post_hoc_aov_Switching_covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_Switching_covar_mod, linfct=mcp(Group="Tukey"))) #save to variable for effect size testing
summary(glht(post_hoc_aov_Switching_covar_mod, linfct=mcp(Timepoint="Tukey")))
#calculate effect size for the significant Group levels
#effect size for sig. post hoc tests - need to use a.tes (from the compute.es package) to account for covariates
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvmMCI_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_CvmMCI_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_CvAD_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_CvAD_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. aMCI 
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Switching.Raw, Group) #count number of participants per group
  r_value <- cor.test(DPRC_neuropsych_data_SCDvAD_Switching.Raw$Switching.Raw,DPRC_neuropsych_data_SCDvAD_Switching.Raw$Age) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value$estimate,q=1) #calculate Cohen's D with the covariate of age
  
  
#####--------------ANCOVA Testing (w/ sex and age) --------------------------------######### 
#run mixed design, 2 x 5 ANCOVA for HayBTime1.Raw
aov_HayBTime1Raw <- aov(HayBTime1.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_HayBTime1Raw)
#aov_HayBTime1Raw <- Anova(lm(HayBTime1.Raw ~ Group*Timepoint + Age + Sex, data=DPRC_neuropsych_data), type = "III") #use this anova test to account for unbalanced designs/sample sizes
#aov_HayBTime1Raw
#effect size w/ covariates
#eta_squared(aov_HayBTime1Raw) #doesn't work - doesn't like the class of the model
sjstats::eta_sq(aov_HayBTime1Raw)
#run mixed design, 2 x 5 ANCOVA for HayBTime2.Raw
aov_HayBTime2Raw <- aov(HayBTime2.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_HayBTime2Raw)
#effect size w/ covariates
sjstats::eta_sq(aov_HayBTime2Raw)
#run mixed design, 2 x 5 ANCOVA for HayBCatA.Raw
aov_HayBCatARaw <- aov(HayBCatA.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_HayBCatARaw)
#effect size w/ covariates
sjstats::eta_sq(aov_HayBCatARaw)
#run mixed design, 2 x 5 ANCOVA for HayBCatB.Raw
aov_HayBCatBRaw <- aov(HayBCatB.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_HayBCatBRaw)
#effect size w/ covariates
sjstats::eta_sq(aov_HayBCatBRaw)
#run mixed design, 2 x 5 ANCOVA for ColorNaming.Raw
aov_ColorNamingRaw <- aov(ColorNaming.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_ColorNamingRaw)
#effect size w/ covariates
sjstats::eta_sq(aov_ColorNamingRaw)
#run mixed design, 2 x 5 ANCOVA for WordReading.Raw
aov_WordReadingRaw <- aov(WordReading.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_WordReadingRaw)
#effect size w/ covariates
sjstats::eta_sq(aov_WordReadingRaw)
#run mixed design, 2 x 5 ANCOVA for WordReading.Raw
aov_InhibitionRaw <- aov(Inhibition.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_InhibitionRaw)
#effect size w/ covariates
sjstats::eta_sq(aov_InhibitionRaw)
#run mixed design, 2 x 5 ANCOVA for Inhibition.Colour.Naming
aov_InhibitionColourNaming <- aov(Inhibition.Colour.Naming~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_InhibitionColourNaming)
#run mixed design, 2 x 5 ANCOVA for Inhibition.Word.Reading
aov_InhibitionWordReading <- aov(Inhibition.Word.Reading~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_InhibitionWordReading)
#effect size w/ covariates
sjstats::eta_sq(aov_InhibitionWordReading)
#run mixed design, 2 x 5 ANCOVA for TrailsA.Raw
aov_TrailsARaw <- aov(TrailsA.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_TrailsARaw)
#effect size w/ covariates
sjstats::eta_sq(aov_TrailsARaw)
#run mixed design, 2 x 5 ANCOVA for TrailsB.Raw
aov_TrailsBRaw <- aov(TrailsB.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_TrailsBRaw)
#effect size w/ covariates
sjstats::eta_sq(aov_TrailsBRaw)
#run mixed design, 2 x 5 ANCOVA for TMT.B.TMT.A
aov_TMTBTMTA <- aov(TMT.B.TMT.A ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_TMTBTMTA)
#effect size w/ covariates
sjstats::eta_sq(aov_TMTBTMTA)
#run mixed design, 2 x 5 ANCOVA for LetFluency.Raw
aov_LetFluencyRaw <- aov(LetFluency.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_LetFluencyRaw)
#effect size w/ covariates
sjstats::eta_sq(aov_LetFluencyRaw)
#run mixed design, 2 x 5 ANCOVA for CatFluency.Raw
aov_CatFluencyRaw <- aov(CatFluency.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_CatFluencyRaw)
#effect size w/ covariates
sjstats::eta_sq(aov_CatFluencyRaw)
#run mixed design, 2 x 5 ANCOVA for Switching.Raw
aov_SwitchingRaw <- aov(Switching.Raw ~ Group*Timepoint + Error(Individual_number/Timepoint) + Age + Sex, data = DPRC_neuropsych_data)
summary(aov_SwitchingRaw)
#effect size w/ covariates
sjstats::eta_sq(aov_SwitchingRaw)

#f/u post hoc tests
#HayBTime1
#remove NAs from dataset for given variable
noNAs_HayBTime1Raw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","HayBTime1.Raw")]
noNAs_HayBTime1Raw <- noNAs_HayBTime1Raw[complete.cases(noNAs_HayBTime1Raw), ]
post_hoc_aov_HayBTime1_2covar_mod <- lme(HayBTime1.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_HayBTime1Raw)
summary(glht(post_hoc_aov_HayBTime1_2covar_mod, linfct=mcp(Group="Tukey"))) #no sig.
#HayBTime2
#remove NAs from dataset for given variable
noNAs_HayBTime2Raw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","HayBTime2.Raw")]
noNAs_HayBTime2Raw <- noNAs_HayBTime2Raw[complete.cases(noNAs_HayBTime2Raw), ]
post_hoc_aov_HayBTime2_2covar_mod <- lme(HayBTime2.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_HayBTime2Raw)
summary(glht(post_hoc_aov_HayBTime2_2covar_mod, linfct=mcp(Group="Tukey"))) #no sig
#HayBCatA
#remove NAs from dataset for given variable
noNAs_HayBCatARaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","HayBCatA.Raw")]
noNAs_HayBCatARaw <- noNAs_HayBCatARaw[complete.cases(noNAs_HayBCatARaw), ]
post_hoc_aov_HayBCatA_2covar_mod <- lme(HayBCatA.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_HayBCatARaw)
summary(glht(post_hoc_aov_HayBCatA_2covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_HayBCatA_2covar_mod, linfct=mcp(Group="Tukey")))
#effect size for sig. post hoc tests
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_HayBCatA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_HayBCatA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_HayBCatA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvAD_HayBCatA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for mMCI vs. AD
  DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 4 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Group <- droplevels(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatA.Raw ~ Age + Sex, data = DPRC_neuropsych_data_mMCIvAD_HayBCatA.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 4'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#HayBCatB
#remove NAs from dataset for given variable
noNAs_HayBCatBRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","HayBCatB.Raw")]
noNAs_HayBCatBRaw <- noNAs_HayBCatBRaw[complete.cases(noNAs_HayBCatBRaw), ]
post_hoc_aov_HayBCatB_2covar_mod <- lme(HayBCatB.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_HayBCatBRaw)
summary(glht(post_hoc_aov_HayBCatB_2covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_HayBCatB_2covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_HayBCatB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_HayBCatB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_HayBCatB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_HayBCatB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. AD
  DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvAD_HayBCatB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for mMCI vs. AD
  DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 4 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw$Group <- droplevels(DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(HayBCatB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_mMCIvAD_HayBCatB.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 4'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#ColorNaming
#remove NAs from dataset for given variable
noNAs_ColorNamingRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","ColorNaming.Raw")]
noNAs_ColorNamingRaw <- noNAs_ColorNamingRaw[complete.cases(noNAs_ColorNamingRaw), ]
post_hoc_aov_ColorNaming_2covar_mod <- lme(ColorNaming.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_ColorNamingRaw)
summary(glht(post_hoc_aov_ColorNaming_2covar_mod, linfct=mcp(Group="Tukey")))
summary(glht(post_hoc_aov_ColorNaming_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#interaction - ColorNaming
#Tukey test
TukeyHSD(aov(ColorNaming.Raw~Group*Timepoint+as.factor(Age)+Sex, data = DPRC_neuropsych_data))   
#subtract differences in F2 - F0 for color naming with covariate age
ColorNaming_diff_data<-DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","ColorNaming.Raw")]
F0_ColorNaming_data<-ColorNaming_diff_data[ColorNaming_diff_data$Timepoint=='F0', ]
F2_ColorNaming_data<-ColorNaming_diff_data[ColorNaming_diff_data$Timepoint=='F2', ]
ColorNamingDiff<-F2_ColorNaming_data$ColorNaming.Raw - F0_ColorNaming_data$ColorNaming.Raw
F0_ColorNaming_data$ColorNamingDiff<-ColorNamingDiff
#run one-way ANOVA on the F2-F0 differences in color naming between groups
ColorNamingDiff_2covar_mod <- lm(ColorNamingDiff ~ Group+Age+Sex, data = F0_ColorNaming_data)
anova(ColorNamingDiff_2covar_mod) #sig. difference
#effect size omnibus ANOVA
etaSquared(ColorNamingDiff_2covar_mod)
#post hoc f/u test
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_ColorNamingDiff_2covar_mod <- glht(ColorNamingDiff_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_ColorNamingDiff_2covar_mod)
confint(post_hoc_ColorNamingDiff_2covar_mod)
t_value_effect_size <- summary(glht(ColorNamingDiff_2covar_mod, linfct = mcp(Group = "Tukey")))
#effect size for sig. post hoc tests with covariate
  #for SCD vs. aMCI
  DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvaMCI_TrailsBDiff, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(ColorNaming.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvaMCI_ColorNamingDiff)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#WordReading
#remove NAs from dataset for given variable
noNAs_WordReadingRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","WordReading.Raw")]
noNAs_WordReadingRaw <- noNAs_WordReadingRaw[complete.cases(noNAs_WordReadingRaw), ]
post_hoc_aov_WordReading_2covar_mod <- lme(WordReading.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_WordReadingRaw)
summary(glht(post_hoc_aov_WordReading_2covar_mod, linfct=mcp(Group="Tukey")))
summary(glht(post_hoc_aov_WordReading_2covar_mod, linfct=mcp(Timepoint="Tukey")))
#Inhibition
#remove NAs from dataset for given variable
noNAs_InhibitionRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","Inhibition.Raw")]
noNAs_InhibitionRaw <- noNAs_InhibitionRaw[complete.cases(noNAs_InhibitionRaw), ]
post_hoc_aov_Inhibition_2covar_mod <- lme(Inhibition.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_InhibitionRaw)
summary(glht(post_hoc_aov_Inhibition_2covar_mod, linfct=mcp(Group="Tukey")))
#Inhibition.Word.Reading
#remove NAs from dataset for given variable
noNAs_InhibitionWordReading <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","Inhibition.Word.Reading")]
noNAs_InhibitionWordReading <- noNAs_InhibitionWordReading[complete.cases(noNAs_InhibitionWordReading), ]
post_hoc_aov_InhibitionWordReading_2covar_mod <- lme(Inhibition.Word.Reading ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_InhibitionWordReading)
summary(glht(post_hoc_aov_InhibitionWordReading_2covar_mod, linfct=mcp(Group="Tukey")))
#TrailsA
#remove NAs from dataset for given variable
noNAs_TrailsARaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","TrailsA.Raw")]
noNAs_TrailsARaw <- noNAs_TrailsARaw[complete.cases(noNAs_TrailsARaw), ]
post_hoc_aov_TrailsA_2covar_mod <- lme(TrailsA.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_TrailsARaw)
summary(glht(post_hoc_aov_TrailsA_2covar_mod, linfct=mcp(Group="Tukey")))
#TrailsB
#remove NAs from dataset for given variable
noNAs_TrailsBRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","TrailsB.Raw")]
noNAs_TrailsBRaw <- noNAs_TrailsBRaw[complete.cases(noNAs_TrailsBRaw), ]
post_hoc_aov_TrailsB_2covar_mod <- lme(TrailsB.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_TrailsBRaw)
summary(glht(post_hoc_aov_TrailsB_2covar_mod, linfct=mcp(Group="Tukey")))
#interaction - TrailsB
#Tukey test
TukeyHSD(aov(TrailsB.Raw~Group*Timepoint+as.factor(Age)+Sex, data = DPRC_neuropsych_data))   
#subtract differences in F2 - F0 for color naming with covariate age
TrailsB_diff_data<-DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","TrailsB.Raw")]
F0_TrailsB_data<-TrailsB_diff_data[TrailsB_diff_data$Timepoint=='F0', ]
F2_TrailsB_data<-TrailsB_diff_data[TrailsB_diff_data$Timepoint=='F2', ]
TrailsBDiff<-F2_TrailsB_data$TrailsB.Raw - F0_TrailsB_data$TrailsB.Raw
F0_TrailsB_data$TrailsBDiff<-TrailsBDiff
#run one-way ANOVA on the F2-F0 differences in color naming between groups
TrailsBDiff_2covar_mod <- lm(TrailsBDiff ~ Group+Age+Sex, data = F0_TrailsB_data)
anova(TrailsBDiff_2covar_mod) #sig. difference
#effect size omnibus ANOVA
etaSquared(TrailsBDiff_2covar_mod)
#post hoc f/u test
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsBDiff_2covar_mod <- glht(TrailsBDiff_2covar_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsBDiff_2covar_mod)
confint(post_hoc_TrailsBDiff_2covar_mod)
t_value_effect_size <- summary(glht(TrailsBDiff_2covar_mod, linfct = mcp(Group = "Tukey")))
#effect size for sig. post hoc tests with covariate
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_TrailsBDiff <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_TrailsBDiff$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_TrailsBDiff$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_TrailsBDiff, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_TrailsBDiff)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvmMCI_TrailsBDiff)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for aMCI vs. mMCI
  DPRC_neuropsych_data_aMCIvmMCI_TrailsBDiff <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_aMCIvmMCI_TrailsBDiff$Group <- droplevels(DPRC_neuropsych_data_aMCIvmMCI_TrailsBDiff$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_aMCIvmMCI_TrailsBDiff, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(TrailsB.Raw ~ Age + Sex, data = DPRC_neuropsych_data_aMCIvmMCI_TrailsBDiff)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 3'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#LetFluency
#remove NAs from dataset for given variable
noNAs_LetFluencyRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","LetFluency.Raw")]
noNAs_LetFluencyRaw <- noNAs_LetFluencyRaw[complete.cases(noNAs_LetFluencyRaw), ]
post_hoc_aov_LetFluency_2covar_mod <- lme(LetFluency.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_LetFluencyRaw)
summary(glht(post_hoc_aov_LetFluency_2covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_LetFluency_2covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_LetFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(LetFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvaMCI_LetFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#CatFluency
#remove NAs from dataset for given variable
noNAs_CatFluencyRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","CatFluency.Raw")]
noNAs_CatFluencyRaw <- noNAs_CatFluencyRaw[complete.cases(noNAs_CatFluencyRaw), ]
post_hoc_aov_CatFluency_2covar_mod <- lme(CatFluency.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_CatFluencyRaw)
summary(glht(post_hoc_aov_CatFluency_2covar_mod, linfct=mcp(Group="Tukey")))
summary(glht(post_hoc_aov_CatFluency_2covar_mod, linfct=mcp(Timepoint="Tukey"))) #not sig.
t_value_effect_size <- summary(glht(post_hoc_aov_CatFluency_2covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. aMCI
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvaMCI_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvaMCI_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. aMCI
  DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvaMCI_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(CatFluency.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_CatFluency.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
#Switching
#remove NAs from dataset for given variable
noNAs_SwitchingRaw <- DPRC_neuropsych_data[,c("ParticipantID","Individual_number","Group","Timepoint","Age","Sex","Switching.Raw")]
noNAs_SwitchingRaw <- noNAs_SwitchingRaw[complete.cases(noNAs_SwitchingRaw), ]
post_hoc_aov_Switching_2covar_mod <- lme(Switching.Raw ~ Group*Timepoint + Age + Sex, random = ~1 | Individual_number/Timepoint, data=noNAs_SwitchingRaw)
summary(glht(post_hoc_aov_Switching_2covar_mod, linfct=mcp(Group="Tukey")))
t_value_effect_size <- summary(glht(post_hoc_aov_Switching_2covar_mod, linfct=mcp(Group="Tukey")))
#calculate effect sizes (Cohen's D) with covariate (using a.tes function from the compute.es package)
  #for C vs. mMCI
  DPRC_neuropsych_data_CvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvmMCI_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvmMCI_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for C vs. AD
  DPRC_neuropsych_data_CvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_CvAD_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_CvAD_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 1'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. aMCI
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvaMCI_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvaMCI_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['3 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. mMCI
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvmMCI_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvmMCI_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['4 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex
  #for SCD vs. AD
  DPRC_neuropsych_data_SCDvAD_Switching.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_Switching.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_Switching.Raw$Group)
  group_number <-dplyr::count(DPRC_neuropsych_data_SCDvAD_Switching.Raw, Group) #count number of participants per group
  mult.r_value_2covar_mod<-summary(lm(Switching.Raw ~ Age + Sex, data = DPRC_neuropsych_data_SCDvAD_Switching.Raw)) #create multiple regression between age, sex, and y-var, and get square root of mult-r squared as the r-value
  r_value <- sqrt(mult.r_value_2covar_mod$r.squared) #find correlation value (r) between dependent variable
  a.tes(t=t_value_effect_size$test$tstat['5 - 2'],n.1=group_number['1','n'],n.2=group_number['2','n'],R=r_value,q=2) #calculate Cohen's D with the covariate of age & sex








  
  
  
  
  
  
  

  
  
  
  
  
  