#Use imputation to fill in the missing data from the DPRC neuropsych data. Will 
#be using the MICE (Multivariate Imputation via Chained Equations) multiple 
#imputation method for this (van Buuren & Groothuis-Oudshoorn, 2010). This is 
#based upon Rubin's Multiple Imputation (Rubin, 2004)

#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 28/06/21

#load libraries via pacman
pacman::p_load(mice, VIM)

#set up pathway
#setwd('/yourpathway/')
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/June/')

#read in csv files (participant file)
DPRC_neuropsych_data <- read.csv("DPRC_neuropsych_data_lined_up_valid_participants_imputation.csv")

#rename some columns
names(DPRC_neuropsych_data)[1] <- "Subject_ID"

#convert variables
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Binary_sex <- as.factor(DPRC_neuropsych_data$Binary_sex)


#-------------------run MICE imputation on the dataset-------------------------#

#check the missing values present in each variable in the dataset
md.pattern(DPRC_neuropsych_data)

#visualise this missing data
mice_plot <- aggr(DPRC_neuropsych_data, col=c('navyblue', 'yellow'), numbers=TRUE, sortVars=TRUE, labels=names(DPRC_neuropsych_data), cex.axis=.7, gap=3, ylab=c("Missing data", "Pattern"))

#MICE details: 
#m = 5, for 5 multiple imputations datasets being done. 
#maxit = 50, for 50 iterations being taken to impute the missing values.The higher value you have, then the better (more accurate) predictions you will have. 
#In regards to the iterations, use 'Bodner's rule of thumb' (White et al., 2011) as a guide --> e.g. 73 imputations w/ 10 iterations
#method = method used in imputation. 'PMM' for Predictive Mean Matching, for numeric variables. LogReg is for categorical data. 
#set seed for 500
imputed_Data <- mice(DPRC_neuropsych_data, m=5, maxit = 50, method = 'pmm', seed = 500)

#check imputted values
summary(imputed_Data)
#e.g., check the TOPF score
imputed_Data$imp$TOPF.Raw

#5 datasets are available, and so you can select any using the complete() function. 
#You want to choose the dataset which represents the mean most closely. 
#e.g. select the 2nd dataset
#completeData <- mice::complete(imputed_Data, 2)

#Build your model using all 5 datasets, and then pool together your results from 
#your models 

#build your predictive models
fit0 <- with(data = imputed_Data, exp = lm(Inhibition.Raw ~ ACE + Group + TOPF.Raw + TrailsA.Raw + TrailsB.Raw + ColorNaming.Raw + WordReading.Raw + LetFluency.Raw + CatFluency.Raw + Switching.Raw + HayBTime1.Raw + HayBTime2.Raw + HayBCatA.Raw + HayBCatB.Raw)) 
fit1 <- with(data = imputed_Data, exp = lm(Inhibition.Raw ~ TrailsB.Raw + Switching.Raw + HayBCatA.Raw)) 


#compare your models to see which variables would contribute best to its predictive power.
stat <- pool.compare(fit1, fit0, method="Wald")

#combine results of all 5 models into 1 model (pooled model)
combine <- pool(fit)
summary(combine)

combine2 <- pool(fit2)
summary(combine2)


#Finally, use this model to fill in the missing values from your dataset. 







