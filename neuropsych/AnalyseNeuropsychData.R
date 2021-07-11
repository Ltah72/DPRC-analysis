#This script will analyse the DPRC neuropsychological assessment data. Will be 
#looking at executive function / cognitive control data. Statistical tests that 
#have been run are an exploratory factor analysis (EFA), multivariate analysis 
#of variance (MANOVAs), and analysis of variance (ANOVAs). 

#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 20/06/21

#load libraries via pacman
pacman::p_load(dplyr, ggplot2, psych, car, multcomp, lsr, BayesFactor, tidyr, GPArotation, corrplot)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#set up pathway
#setwd('/yourpathway/')
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/June/')

#read in csv files (participant file)
DPRC_neuropsych_data <- read.csv("DPRC_neuropsych_data_lined_up_valid_participants.csv")

#rename some columns
names(DPRC_neuropsych_data)[1] <- "Subject_ID"

#convert variables
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Binary_sex <- as.factor(DPRC_neuropsych_data$Binary_sex)


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
exec_raw_data <- subset(DPRC_neuropsych_data, select = c(Subject_ID,
                                                         Age,
                                                         Classification,
                                                         Group,
                                                         Sex,
                                                         Binary_sex,
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
exec_zscores_data <- subset(DPRC_neuropsych_data, select = c(Subject_ID,
                                                         Age,
                                                         Classification,
                                                         Group,
                                                         Sex,
                                                         Binary_sex,
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
                                        Subject_ID,
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
ggplot(subset(exec_func_zscores_data_long, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Z_scores, fill = Group)) + 
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
                                         LetFluency.Raw, 
                                         CatFluency.Raw, 
                                         HayBTime1.Raw) ~ Group, data = DPRC_neuropsych_data) 

#overall model
summary(manova_proc_speed_raw_mod)
#anova outputs
summary.aov(manova_proc_speed_raw_mod)

#for z-scores
manova_proc_speed_z_mod <- manova(cbind(TrailsA.Z, 
                                       ColorNaming.Z, 
                                       WordReading.Z, 
                                       LetFluency.Z, 
                                       CatFluency.Z, 
                                       HayBTime1.z) ~ Group, data = DPRC_neuropsych_data) 


#overall model
summary(manova_proc_speed_z_mod)
#anova outputs
summary.aov(manova_proc_speed_z_mod)


#plot the manova data using the z-scores of the variables, so that it is all on the same scale. 
#put proc speed variable z-scores onto a new dataset, as long format
proc_speed_zscores_data <- dplyr::select(DPRC_neuropsych_data, 
                                        Subject_ID,
                                        Group,
                                        TrailsA.Z, 
                                        ColorNaming.Z, 
                                        WordReading.Z, 
                                        LetFluency.Z, 
                                        CatFluency.Z,
                                        HayBTime1.z)
#put into long format
proc_speed_zscores_data_long <- gather(exec_func_zscores_data, 
                                      "Processing_Speeds",
                                      "Z_scores", 
                                      TrailsA.Z, 
                                      ColorNaming.Z, 
                                      WordReading.Z, 
                                      LetFluency.Z, 
                                      CatFluency.Z, 
                                      HayBTime1.z)


#plot data - z-score for processing speeds 
ggplot(subset(proc_speed_zscores_data_long, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Z_scores, fill = Group)) + 
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
                                        Subject_ID,
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
ggplot(subset(inhibition_zscores_data_long, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Z_scores, fill = Group)) + 
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


####---------------------------------ANOVAs---------------------------------####
#run ANOVAs on your data:
#1.Test of premorbid functioning (TOPF) ---------------------------------------#
#plot TOPF.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = TOPF.Raw, fill = Group)) + 
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

#plot TOPF.Z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = TOPF.Z, fill = Group)) + 
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
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = HayBTime1.Raw, fill = Group)) + 
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

#plot HayBTime1.z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = HayBTime1.z, fill = Group)) + 
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
#run ANOVA for HayBTime1.z
HayBTime1.Z_mod <- lm(HayBTime1.z ~ Group, data = DPRC_neuropsych_data)
anova(HayBTime1.Z_mod)

#plot HayBTime2.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = HayBTime2.Raw, fill = Group)) + 
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

#plot HayBTime2.z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = HayBTime2.z, fill = Group)) + 
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
#run ANOVA for HayBTime2.Z
HayBTime2.Z_mod <- lm(HayBTime2.z ~ Group, data = DPRC_neuropsych_data)
anova(HayBTime2.Z_mod)

#plot HayBCatA.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = HayBCatA.Raw, fill = Group)) + 
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
For_Bay_data_noNas_neuropsych_HayBCatA.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, HayBCatA.Raw)
For_Bay_data_noNas_neuropsych_HayBCatA.Raw <- na.omit(For_Bay_data_noNas_neuropsych_HayBCatA.Raw)
anovaBF(HayBCatA.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_HayBCatA.Raw) 

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
For_Bay_data_noNas_neuropsych_HayBCatA.Z <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, HayBCatA.z)
For_Bay_data_noNas_neuropsych_HayBCatA.Z <- na.omit(For_Bay_data_noNas_neuropsych_HayBCatA.Z)
anovaBF(HayBCatA.z ~ Group, data = For_Bay_data_noNas_neuropsych_HayBCatA.Z) 

#plot HayBCatB.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = HayBCatB.Raw, fill = Group)) + 
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
For_Bay_data_noNas_neuropsych_HayBCatB.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, HayBCatB.Raw)
For_Bay_data_noNas_neuropsych_HayBCatB.Raw <- na.omit(For_Bay_data_noNas_neuropsych_HayBCatB.Raw)
anovaBF(HayBCatB.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_HayBCatB.Raw) 

#plot HayBCatB.z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = HayBCatB.z, fill = Group)) + 
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
For_Bay_data_noNas_neuropsych_HayBCatB.Z <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, HayBCatB.z)
For_Bay_data_noNas_neuropsych_HayBCatB.Z <- na.omit(For_Bay_data_noNas_neuropsych_HayBCatB.Z)
anovaBF(HayBCatB.z ~ Group, data = For_Bay_data_noNas_neuropsych_HayBCatB.Z) 


#3.D-KEFS Stroop Task ---------------------------------------------------------#
#plot ColorNaming.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = ColorNaming.Raw, fill = Group)) + 
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
#effect size for sig. post hoc tests
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_ColorNaming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_ColorNaming.Raw$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_ColorNaming.Raw$Group)
  cohensD(ColorNaming.Raw ~ Group, data = DPRC_neuropsych_data_CvmMCI_ColorNaming.Raw) #this looks like Hedges' g? 
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_ColorNaming.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_ColorNaming.Raw$Group <- droplevels(DPRC_neuropsych_data_CvAD_ColorNaming.Raw$Group)
  cohensD(ColorNaming.Raw ~ Group, data = DPRC_neuropsych_data_CvAD_ColorNaming.Raw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_ColorNaming.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, ColorNaming.Raw)
For_Bay_data_noNas_neuropsych_ColorNaming.Raw <- na.omit(For_Bay_data_noNas_neuropsych_ColorNaming.Raw)
anovaBF(ColorNaming.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_ColorNaming.Raw) 

#plot ColorNaming.Z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = ColorNaming.Z, fill = Group)) + 
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
#run ANOVA for ColorNaming.Z
ColorNaming.Z_mod <- lm(ColorNaming.Z ~ Group, data = DPRC_neuropsych_data)
anova(ColorNaming.Z_mod)
#effect size omnibus ANOVA
etaSquared(ColorNaming.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_ColorNaming.Z_mod <- glht(ColorNaming.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_ColorNaming.Z_mod)
confint(post_hoc_ColorNaming.Z_mod)
#check descriptive statistics per each group
ColorNaming.Z_descrip <- describeBy(DPRC_neuropsych_data$ColorNaming.Z, DPRC_neuropsych_data$Group)
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
For_Bay_data_noNas_neuropsych_ColorNaming.Z <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, ColorNaming.Z)
For_Bay_data_noNas_neuropsych_ColorNaming.Z <- na.omit(For_Bay_data_noNas_neuropsych_ColorNaming.Z)
anovaBF(ColorNaming.Z ~ Group, data = For_Bay_data_noNas_neuropsych_ColorNaming.Z) 

#plot WordReading.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = WordReading.Raw, fill = Group)) + 
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
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_WordReading.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, WordReading.Raw)
For_Bay_data_noNas_neuropsych_WordReading.Raw <- na.omit(For_Bay_data_noNas_neuropsych_WordReading.Raw)
anovaBF(WordReading.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_WordReading.Raw) 

#plot WordReading.Z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = WordReading.Z, fill = Group)) + 
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
#run ANOVA for WordReading.Z
WordReading.Z_mod <- lm(WordReading.Z ~ Group, data = DPRC_neuropsych_data)
anova(WordReading.Z_mod)

#plot Inhibition.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Inhibition.Raw, fill = Group)) + 
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
For_Bay_data_noNas_neuropsych_Inhibition.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, Inhibition.Raw)
For_Bay_data_noNas_neuropsych_Inhibition.Raw<- na.omit(For_Bay_data_noNas_neuropsych_Inhibition.Raw)
anovaBF(Inhibition.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_Inhibition.Raw) 

#plot Inhibition.Z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Inhibition.Z, fill = Group)) + 
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
For_Bay_data_noNas_neuropsych_Inhibition.Z <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, Inhibition.Z)
For_Bay_data_noNas_neuropsych_Inhibition.Z<- na.omit(For_Bay_data_noNas_neuropsych_Inhibition.Z)
anovaBF(Inhibition.Z ~ Group, data = For_Bay_data_noNas_neuropsych_Inhibition.Z) 


#4.D-KEFS Verbal Fluency + Category Fluency Task ------------------------------#
#plot LetFluency.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = LetFluency.Raw, fill = Group)) + 
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
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_LetFluency.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, LetFluency.Raw)
For_Bay_data_noNas_neuropsych_LetFluency.Raw<- na.omit(For_Bay_data_noNas_neuropsych_LetFluency.Raw)
anovaBF(LetFluency.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_LetFluency.Raw) 

#plot LetFluency.Z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = LetFluency.Z, fill = Group)) + 
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
#run ANOVA for LetFluency.Raw
LetFluency.Z_mod <- lm(LetFluency.Z ~ Group, data = DPRC_neuropsych_data)
anova(LetFluency.Z_mod)
#effect size omnibus ANOVA
etaSquared(LetFluency.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_LetFluency.Z_mod <- glht(LetFluency.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_LetFluency.Z_mod)
confint(post_hoc_LetFluency.Z_mod)
#check descriptive statistics per each group
LetFluency.Z_descrip <- describeBy(DPRC_neuropsych_data$LetFluency.Z, DPRC_neuropsych_data$Group)
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
For_Bay_data_noNas_neuropsych_LetFluency.Z <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, LetFluency.Z)
For_Bay_data_noNas_neuropsych_LetFluency.Z<- na.omit(For_Bay_data_noNas_neuropsych_LetFluency.Z)
anovaBF(LetFluency.Z ~ Group, data = For_Bay_data_noNas_neuropsych_LetFluency.Z) 

#plot CatFluency.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = CatFluency.Raw, fill = Group)) + 
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
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw<- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw$Group)
  cohensD(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_CatFluency.Raw) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatFluency.Raw$Group)
  cohensD(CatFluency.Raw ~ Group, data = DPRC_neuropsych_data_SCDvAD_CatFluency.Raw) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_CatFluency.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, CatFluency.Raw)
For_Bay_data_noNas_neuropsych_CatFluency.Raw <- na.omit(For_Bay_data_noNas_neuropsych_CatFluency.Raw)
anovaBF(CatFluency.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_CatFluency.Raw) 

#plot CatFluency.Z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = CatFluency.Z, fill = Group)) + 
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
#run ANOVA for CatFluency.Z
CatFluency.Z_mod <- lm(CatFluency.Z ~ Group, data = DPRC_neuropsych_data)
anova(CatFluency.Z_mod)
#effect size omnibus ANOVA
etaSquared(CatFluency.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_CatFluency.Z_mod <- glht(CatFluency.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_CatFluency.Z_mod)
confint(post_hoc_CatFluency.Z_mod)
#check descriptive statistics per each group
CatFluency.Z_descrip <- describeBy(DPRC_neuropsych_data$CatFluency.Z, DPRC_neuropsych_data$Group)
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
For_Bay_data_noNas_neuropsych_CatFluency.Z <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, CatFluency.Z)
For_Bay_data_noNas_neuropsych_CatFluency.Z <- na.omit(For_Bay_data_noNas_neuropsych_CatFluency.Z)
anovaBF(CatFluency.Z ~ Group, data = For_Bay_data_noNas_neuropsych_CatFluency.Z) 

#plot Switching.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Switching.Raw, fill = Group)) + 
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
For_Bay_data_noNas_neuropsych_Switching.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, Switching.Raw)
For_Bay_data_noNas_neuropsych_Switching.Raw <- na.omit(For_Bay_data_noNas_neuropsych_Switching.Raw)
anovaBF(Switching.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_Switching.Raw) 

#plot Switching.z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = Switching.z, fill = Group)) + 
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
#run ANOVA for Switching.Z
Switching.Z_mod <- lm(Switching.z ~ Group, data = DPRC_neuropsych_data)
anova(Switching.Z_mod)
#effect size omnibus ANOVA
etaSquared(Switching.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_Switching.Z_mod <- glht(Switching.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_Switching.Z_mod)
confint(post_hoc_Switching.Z_mod)
#check descriptive statistics per each group
Switching.Z_descrip <- describeBy(DPRC_neuropsych_data$Switching.Z, DPRC_neuropsych_data$Group)
#effect size for sig. post hoc tests
  #for C vs. aMCI  
  DPRC_neuropsych_data_CvaMCI_Switching.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_CvaMCI_Switching.Z$Group <- droplevels(DPRC_neuropsych_data_CvaMCI_Switching.Z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_CvaMCI_Switching.Z) #this looks like Hedges' g? 
  #for C vs. mMCI 
  DPRC_neuropsych_data_CvmMCI_Switching.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_CvmMCI_Switching.Z$Group <- droplevels(DPRC_neuropsych_data_CvmMCI_Switching.Z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_CvmMCI_Switching.Z) #this looks like Hedges' g?
  #for C vs. AD 
  DPRC_neuropsych_data_CvAD_Switching.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 1 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_CvAD_Switching.Z$Group <- droplevels(DPRC_neuropsych_data_CvAD_Switching.Z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_CvAD_Switching.Z) #this looks like Hedges' g?
  #for SCD vs. aMCI 
  DPRC_neuropsych_data_SCDvaMCI_Switching.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 3)
  DPRC_neuropsych_data_SCDvaMCI_Switching.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvaMCI_Switching.Z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_SCDvaMCI_Switching.Z) #this looks like Hedges' g?
  #for SCD vs. mMCI 
  DPRC_neuropsych_data_SCDvmMCI_Switching.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 4)
  DPRC_neuropsych_data_SCDvmMCI_Switching.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvmMCI_Switching.Z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_SCDvmMCI_Switching.Z) #this looks like Hedges' g?
  #for SCD vs. AD 
  DPRC_neuropsych_data_SCDvAD_CatFluency.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 2 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_SCDvAD_CatFluency.Z$Group <- droplevels(DPRC_neuropsych_data_SCDvAD_CatFluency.Z$Group)
  cohensD(CatFluency.Z ~ Group, data = DPRC_neuropsych_data_SCDvAD_CatFluency.Z) #this looks like Hedges' g? 
  #for aMCI vs. AD 
  DPRC_neuropsych_data_aMCIvAD_Switching.Z <- subset(DPRC_neuropsych_data, DPRC_neuropsych_data$Group == 3 | DPRC_neuropsych_data$Group == 5)
  DPRC_neuropsych_data_aMCIvAD_Switching.Z$Group <- droplevels(DPRC_neuropsych_data_aMCIvAD_Switching.Z$Group)
  cohensD(Switching.z ~ Group, data = DPRC_neuropsych_data_aMCIvAD_Switching.Z) #this looks like Hedges' g? 
#run Bayesian ANOVA
For_Bay_data_noNas_neuropsych_Switching.Z <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, Switching.z)
For_Bay_data_noNas_neuropsych_Switching.Z <- na.omit(For_Bay_data_noNas_neuropsych_Switching.Z)
anovaBF(Switching.z ~ Group, data = For_Bay_data_noNas_neuropsych_Switching.Z) 

#5.Trail Making Test (TMT) A & B ----------------------------------------------#
#plot TrailsA.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = TrailsA.Raw, fill = Group)) + 
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
For_Bay_data_noNas_neuropsych_TrailsA.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, TrailsA.Raw)
For_Bay_data_noNas_neuropsych_TrailsA.Raw <- na.omit(For_Bay_data_noNas_neuropsych_TrailsA.Raw)
anovaBF(TrailsA.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_TrailsA.Raw) 

#plot TrailsA.Z (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = TrailsA.Z, fill = Group)) + 
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
#run ANOVA for TrailsA.Z
TrailsA.Z_mod <- lm(TrailsA.Z ~ Group, data = DPRC_neuropsych_data)
anova(TrailsA.Z_mod)
#effect size omnibus ANOVA
etaSquared(TrailsA.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsA.Z_mod <- glht(TrailsA.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsA.Z_mod)
confint(post_hoc_TrailsA.Z_mod)
#check descriptive statistics per each group
TrailsA.Z_descrip <- describeBy(DPRC_neuropsych_data$TrailsA.Z, DPRC_neuropsych_data$Group)
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
For_Bay_data_noNas_neuropsych_TrailsA.Z <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, TrailsA.Z)
For_Bay_data_noNas_neuropsych_TrailsA.Z <- na.omit(For_Bay_data_noNas_neuropsych_TrailsA.Z)
anovaBF(TrailsA.Z ~ Group, data = For_Bay_data_noNas_neuropsych_TrailsA.Z) 

#plot TrailsB.Raw (raincloud plot)
ggplot(subset(DPRC_neuropsych_data, Group %in% c("1", "2", "3", "4", "5")), aes(x = Group, y = TrailsB.Raw, fill = Group)) + 
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
For_Bay_data_noNas_neuropsych_TrailsB.Raw <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, TrailsB.Raw)
For_Bay_data_noNas_neuropsych_TrailsB.Raw <- na.omit(For_Bay_data_noNas_neuropsych_TrailsB.Raw)
anovaBF(TrailsB.Raw ~ Group, data = For_Bay_data_noNas_neuropsych_TrailsB.Raw) 

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
#run ANOVA for TrailsB.Z
TrailsB.Z_mod <- lm(TrailsB.Z ~ Group, data = DPRC_neuropsych_data)
anova(TrailsB.Z_mod)
#effect size omnibus ANOVA
etaSquared(TrailsB.Z_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_TrailsB.Z_mod <- glht(TrailsB.Z_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_TrailsB.Z_mod)
confint(post_hoc_TrailsB.Z_mod)
#check descriptive statistics per each group
TrailsB.Z_descrip <- describeBy(DPRC_neuropsych_data$TrailsB.Z, DPRC_neuropsych_data$Group)
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
For_Bay_data_noNas_neuropsych_TrailsB.Z <- dplyr::select(DPRC_neuropsych_data, Subject_ID, Group, TrailsB.Z)
For_Bay_data_noNas_neuropsych_TrailsB.Z <- na.omit(For_Bay_data_noNas_neuropsych_TrailsB.Z)
anovaBF(TrailsB.Z ~ Group, data = For_Bay_data_noNas_neuropsych_TrailsB.Z) 









