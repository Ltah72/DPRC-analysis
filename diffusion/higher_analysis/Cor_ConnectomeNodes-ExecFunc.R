#This script will analyse the correlation between the connectome node values and
#and the neuropsychological assesssment data. 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 03/01/22

#------------------------------Setting up--------------------------------------#
#install packages/open libraries
pacman::p_load(ggplot2, ppcor, dplyr, tidyr)

#first read in the neuropsych data file: 
#setwd('/yourpathway/')
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')

#----- Neuropsych data --------------------------------------------------------#

#read in csv files (participant file)
DPRC_neuropsych_data <- read.csv("cross-sectional_DPRC_neuropsych_data_lined_up_valid_participants.csv")

#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'

#convert variables
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
DPRC_neuropsych_data$Sex <- as.factor(DPRC_neuropsych_data$Sex)
DPRC_neuropsych_data$Sex_binary <- as.factor(DPRC_neuropsych_data$Sex_binary)

#for the connectome data, omit 2 participants from the dataset (DDPRC0025F0 (row 221) & DDPRC0029F0 (row 223))
omit_participants <- c("221","223")
DPRC_neuropsych_data <- DPRC_neuropsych_data[!(row.names(DPRC_neuropsych_data) %in% omit_participants), ]

#create new columns for the connectome data
small_node_data <- setNames(data.frame(matrix(ncol = 10, nrow = 227)), c("L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                                    "R-DLPFC--MidCing","L-DLPFC--L-Par",
                                                    "R-DLPFC--R-Par","L-DLPFC--R-Par",
                                                    "R-DLPFC--L-Par","L-Par--R-Par",
                                                    "L-Par--MidCing","R-Par--MidCing"))

#navigate to connectome data pathway
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/FPN_Files/weighted/small-nodes')

#collate connectome values into dataset
filenames <- list.files(path='V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/FPN_Files/weighted/small-nodes')
r.nms <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
AllData <- lapply(filenames, read.table, header=FALSE, row.names=r.nms, col.names="V1")

for(i in sequence(length(filenames))){
  for (j in 1:10){
    small_node_data[i,j]<- AllData[[i]][[1]][j] 
  }
}

#combine neuropsych and connectome dataset
Connectome_data <- cbind(DPRC_neuropsych_data,small_node_data)


#### ----------------------- PLS visualisation ---------------------------- ####


#for all 5 inhibition tests data (barplot):
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.16374236,
                                 0.19805148,
                                 0.19091377,
                                 0.040685397,
                                 0.059221849)*-1 
ulimit_PLS_FD_barplot <- c(0.303537055850029,
                           0.326308354735375,
                           0.301895633339882,
                           0.174151211977005,
                           0.168247200548649)*-1
llimit_PLS_FD_barplot <- c(0.0683764740824699,
                           0.0977761261165142,
                           0.110962077975273,
                           -0.0429181344807148,
                           -0.00895849987864494)*-1

significance_legend<-c('Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable','Does Not Reliably Contribute to Latent Variable')
df_PLS_FD_barplot <- data.frame(neuropsych_test_names_PLS,significance_legend,PLS_FD_bootstrap_corr_values,ulimit_PLS_FD_barplot,llimit_PLS_FD_barplot)
#convert to factor variables
df_PLS_FD_barplot$neuropsych_test_names_PLS <- factor(df_PLS_FD_barplot$neuropsych_test_names_PLS, levels=c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError'))
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


#for inhibition 5 FD:
PLS_output_FD_inhib5_data <- c(-0.084862299,	0.093119934,	0.14467922,	0.051879663,	0.18415737,	  -0.040196970,	-0.047457095,	-0.068977371,	0.066325560,	0.045272447,
                               -0.074255034,	0.13599660,	  0.17942502,	0.050874077,	0.082799844,	-0.061809991,	-0.017830579,	-0.095650934,	0.13915233,	0.068397991,
                               -0.018515429,	0.11484726,	  0.18373598,	0.10852513,	  0.17443596,	  -0.0013519998, 0.050445873,	-0.055006236,	0.090899222,	0.067352593,
                               -0.084564656,	0.025944147,	0.056567095,-0.090561010,	-0.020252027,	-0.022282815,	-0.056460965,	-0.044417012,	-0.021510171,	-0.090607285,
                               -0.0060276738,	0.030635329,	0.051555585, 0.10613545,	0.091820262,	-0.050619118,	 0.034042846,	 0.086551219,	0.038677871,	-0.015485675)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=10,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c("L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                        "R-DLPFC--MidCing","L-DLPFC--L-Par",
                                        "R-DLPFC--R-Par","L-DLPFC--R-Par",
                                        "R-DLPFC--L-Par","L-Par--R-Par",
                                        "L-Par--MidCing","R-Par--MidCing")
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black",col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)






#run correlation tests between variables
#Left DLPFC -- Right DLPFC
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$HayBTime2.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$ColorNaming.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$Inhibition.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$TrailsA.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$TrailsB.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$CatFluency.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$Switching.Raw) #not sig
#Left DLPFC -- Mid Cingulate
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$Inhibition.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$TrailsA.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$TrailsB.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$CatFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$Switching.Raw) #sig
#Right DLPFC -- Mid Cingulate
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$Inhibition.Raw) #sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$TrailsA.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$TrailsB.Raw) #sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$CatFluency.Raw) #sig
cor.test(Connectome_data$`R-DLPFC--MidCing`, Connectome_data$Switching.Raw) #sig
#Left DLPFC -- Left Parietal
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$HayBCatA.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$Inhibition.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$TrailsA.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$TrailsB.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$CatFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$Switching.Raw) # sig
#Right DLPFC -- Right Parietal
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$Inhibition.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$TrailsA.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$TrailsB.Raw) #sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$CatFluency.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--R-Par`, Connectome_data$Switching.Raw) #not sig
#Left DLPFC -- Right Parietal
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$Inhibition.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$TrailsA.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$TrailsB.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$CatFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-Par`, Connectome_data$Switching.Raw) #not sig
#Right DLPFC -- Left Parietal
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$Inhibition.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$TrailsA.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$TrailsB.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$CatFluency.Raw) #not sig
cor.test(Connectome_data$`R-DLPFC--L-Par`, Connectome_data$Switching.Raw) #not sig
#Left Parietal -- Right Parietal
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$HayBTime1.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$HayBTime2.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$ColorNaming.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$Inhibition.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$TrailsA.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$TrailsB.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$CatFluency.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$Switching.Raw) # sig
#Left Parietal -- Mid Cingulate
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$HayBCatA.Raw) #sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$Inhibition.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$TrailsA.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$TrailsB.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$CatFluency.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$Switching.Raw) #not sig
#Right Parietal -- Mid Cingulate
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$Inhibition.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$TrailsA.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$TrailsB.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$CatFluency.Raw) #not sig
cor.test(Connectome_data$`R-Par--MidCing`, Connectome_data$Switching.Raw) #not sig

#Run partial correlation (controlling for covariates, age and sex) - note that
#missing values are not allowed.
#for TrailsA:
Connectome_data_TrailsA <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "TrailsA.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_TrailsA <- Connectome_data_TrailsA %>% drop_na(TrailsA.Raw)
Connectome_data_TrailsA$Sex_binary <- as.numeric(Connectome_data_TrailsA$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_TrailsA$`L-DLPFC--R-DLPFC`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")]) 
pcor.test(Connectome_data_TrailsA$`L-DLPFC--MidCing`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsA$`R-DLPFC--MidCing`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsA$`L-DLPFC--L-Par`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsA$`R-DLPFC--R-Par`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsA$`L-DLPFC--R-Par`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsA$`R-DLPFC--L-Par`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsA$`L-Par--R-Par`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsA$`L-Par--MidCing`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsA$`R-Par--MidCing`, Connectome_data_TrailsA$TrailsA.Raw, Connectome_data_TrailsA[,c("Age", "Sex_binary")])
#for TrailsB:
Connectome_data_TrailsB <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "TrailsB.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_TrailsB <- Connectome_data_TrailsB %>% drop_na(TrailsB.Raw)
Connectome_data_TrailsB$Sex_binary <- as.numeric(Connectome_data_TrailsB$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_TrailsB$`L-DLPFC--R-DLPFC`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsB$`L-DLPFC--MidCing`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsB$`R-DLPFC--MidCing`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(Connectome_data_TrailsB$`L-DLPFC--L-Par`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsB$`R-DLPFC--R-Par`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")]) #sig
pcor.test(Connectome_data_TrailsB$`L-DLPFC--R-Par`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsB$`R-DLPFC--L-Par`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsB$`L-Par--R-Par`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsB$`L-Par--MidCing`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_TrailsB$`R-Par--MidCing`, Connectome_data_TrailsB$TrailsB.Raw, Connectome_data_TrailsB[,c("Age", "Sex_binary")])
#for HayBTime1:
Connectome_data_HayBTime1 <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "HayBTime1.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_HayBTime1 <- Connectome_data_HayBTime1 %>% drop_na(HayBTime1.Raw)
Connectome_data_HayBTime1$Sex_binary <- as.numeric(Connectome_data_HayBTime1$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_HayBTime1$`L-DLPFC--R-DLPFC`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime1$`L-DLPFC--MidCing`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime1$`R-DLPFC--MidCing`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime1$`L-DLPFC--L-Par`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime1$`R-DLPFC--R-Par`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime1$`L-DLPFC--R-Par`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime1$`R-DLPFC--L-Par`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime1$`L-Par--R-Par`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime1$`L-Par--MidCing`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime1$`R-Par--MidCing`, Connectome_data_HayBTime1$HayBTime1.Raw, Connectome_data_HayBTime1[,c("Age", "Sex_binary")])
#for HayBTime2:
Connectome_data_HayBTime2 <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "HayBTime2.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_HayBTime2 <- Connectome_data_HayBTime2 %>% drop_na(HayBTime2.Raw)
Connectome_data_HayBTime2$Sex_binary <- as.numeric(Connectome_data_HayBTime2$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_HayBTime2$`L-DLPFC--R-DLPFC`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime2$`L-DLPFC--MidCing`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime2$`R-DLPFC--MidCing`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime2$`L-DLPFC--L-Par`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime2$`R-DLPFC--R-Par`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime2$`L-DLPFC--R-Par`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime2$`R-DLPFC--L-Par`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime2$`L-Par--R-Par`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime2$`L-Par--MidCing`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBTime2$`R-Par--MidCing`, Connectome_data_HayBTime2$HayBTime2.Raw, Connectome_data_HayBTime2[,c("Age", "Sex_binary")])
#for HayBCatA:
Connectome_data_HayBCatA <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "HayBCatA.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_HayBCatA <- Connectome_data_HayBCatA %>% drop_na(HayBCatA.Raw)
Connectome_data_HayBCatA$Sex_binary <- as.numeric(Connectome_data_HayBCatA$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_HayBCatA$`L-DLPFC--R-DLPFC`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatA$`L-DLPFC--MidCing`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatA$`R-DLPFC--MidCing`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatA$`L-DLPFC--L-Par`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")]) #sig
pcor.test(Connectome_data_HayBCatA$`R-DLPFC--R-Par`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatA$`L-DLPFC--R-Par`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatA$`R-DLPFC--L-Par`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatA$`L-Par--R-Par`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")]) 
pcor.test(Connectome_data_HayBCatA$`L-Par--MidCing`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")]) #sig
pcor.test(Connectome_data_HayBCatA$`R-Par--MidCing`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")])
#for HayBCatB:
Connectome_data_HayBCatB <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "HayBCatB.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_HayBCatB <- Connectome_data_HayBCatB %>% drop_na(HayBCatB.Raw)
Connectome_data_HayBCatB$Sex_binary <- as.numeric(Connectome_data_HayBCatB$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_HayBCatB$`L-DLPFC--R-DLPFC`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatB$`L-DLPFC--MidCing`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatB$`R-DLPFC--MidCing`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatB$`L-DLPFC--L-Par`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatB$`R-DLPFC--R-Par`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatB$`L-DLPFC--R-Par`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatB$`R-DLPFC--L-Par`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatB$`L-Par--R-Par`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatB$`L-Par--MidCing`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_HayBCatB$`R-Par--MidCing`, Connectome_data_HayBCatB$HayBCatB.Raw, Connectome_data_HayBCatB[,c("Age", "Sex_binary")])
#for ColorNaming:
Connectome_data_ColorNaming <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "ColorNaming.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_ColorNaming <- Connectome_data_ColorNaming %>% drop_na(ColorNaming.Raw)
Connectome_data_ColorNaming$Sex_binary <- as.numeric(Connectome_data_ColorNaming$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_ColorNaming$`L-DLPFC--R-DLPFC`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_ColorNaming$`L-DLPFC--MidCing`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_ColorNaming$`R-DLPFC--MidCing`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_ColorNaming$`L-DLPFC--L-Par`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_ColorNaming$`R-DLPFC--R-Par`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_ColorNaming$`L-DLPFC--R-Par`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_ColorNaming$`R-DLPFC--L-Par`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_ColorNaming$`L-Par--R-Par`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_ColorNaming$`L-Par--MidCing`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_ColorNaming$`R-Par--MidCing`, Connectome_data_ColorNaming$ColorNaming.Raw, Connectome_data_ColorNaming[,c("Age", "Sex_binary")])
#for WordReading:
Connectome_data_WordReading <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "WordReading.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_WordReading <- Connectome_data_WordReading %>% drop_na(WordReading.Raw)
Connectome_data_WordReading$Sex_binary <- as.numeric(Connectome_data_WordReading$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_WordReading$`L-DLPFC--R-DLPFC`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_WordReading$`L-DLPFC--MidCing`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_WordReading$`R-DLPFC--MidCing`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_WordReading$`L-DLPFC--L-Par`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_WordReading$`R-DLPFC--R-Par`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_WordReading$`L-DLPFC--R-Par`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_WordReading$`R-DLPFC--L-Par`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_WordReading$`L-Par--R-Par`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_WordReading$`L-Par--MidCing`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_WordReading$`R-Par--MidCing`, Connectome_data_WordReading$WordReading.Raw, Connectome_data_WordReading[,c("Age", "Sex_binary")])
#for Inhibition:
Connectome_data_Inhibition <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "Inhibition.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_Inhibition <- Connectome_data_Inhibition %>% drop_na(Inhibition.Raw)
Connectome_data_Inhibition$Sex_binary <- as.numeric(Connectome_data_Inhibition$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_Inhibition$`L-DLPFC--R-DLPFC`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`L-DLPFC--MidCing`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`R-DLPFC--MidCing`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(Connectome_data_Inhibition$`L-DLPFC--L-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`R-DLPFC--R-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`L-DLPFC--R-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`R-DLPFC--L-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`L-Par--R-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`L-Par--MidCing`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")]) #sig
pcor.test(Connectome_data_Inhibition$`R-Par--MidCing`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
#for LetFluency:
Connectome_data_LetFluency <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "LetFluency.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_LetFluency <- Connectome_data_LetFluency %>% drop_na(LetFluency.Raw)
Connectome_data_LetFluency$Sex_binary <- as.numeric(Connectome_data_LetFluency$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_LetFluency$`L-DLPFC--R-DLPFC`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_LetFluency$`L-DLPFC--MidCing`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_LetFluency$`R-DLPFC--MidCing`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_LetFluency$`L-DLPFC--L-Par`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_LetFluency$`R-DLPFC--R-Par`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_LetFluency$`L-DLPFC--R-Par`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_LetFluency$`R-DLPFC--L-Par`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_LetFluency$`L-Par--R-Par`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_LetFluency$`L-Par--MidCing`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_LetFluency$`R-Par--MidCing`, Connectome_data_LetFluency$LetFluency.Raw, Connectome_data_LetFluency[,c("Age", "Sex_binary")])
#for CatFluency:
Connectome_data_CatFluency <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "CatFluency.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_CatFluency <- Connectome_data_CatFluency %>% drop_na(CatFluency.Raw)
Connectome_data_CatFluency$Sex_binary <- as.numeric(Connectome_data_CatFluency$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_CatFluency$`L-DLPFC--R-DLPFC`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_CatFluency$`L-DLPFC--MidCing`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_CatFluency$`R-DLPFC--MidCing`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")]) #sig
pcor.test(Connectome_data_CatFluency$`L-DLPFC--L-Par`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_CatFluency$`R-DLPFC--R-Par`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_CatFluency$`L-DLPFC--R-Par`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_CatFluency$`R-DLPFC--L-Par`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_CatFluency$`L-Par--R-Par`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_CatFluency$`L-Par--MidCing`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_CatFluency$`R-Par--MidCing`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
#for Switching:
Connectome_data_Switching <- Connectome_data[c("Group", "Age", "Sex", "Sex_binary", "Switching.Raw", "L-DLPFC--R-DLPFC","L-DLPFC--MidCing",
                                             "R-DLPFC--MidCing","L-DLPFC--L-Par","R-DLPFC--R-Par","L-DLPFC--R-Par","R-DLPFC--L-Par",
                                             "L-Par--R-Par","L-Par--MidCing","R-Par--MidCing")]
Connectome_data_Switching <- Connectome_data_Switching %>% drop_na(Switching.Raw)
Connectome_data_Switching$Sex_binary <- as.numeric(Connectome_data_Switching$Sex_binary)
#for L-DLPFC--R-DLPFC
pcor.test(Connectome_data_Switching$`L-DLPFC--R-DLPFC`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`L-DLPFC--MidCing`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`R-DLPFC--MidCing`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")]) #sig
pcor.test(Connectome_data_Switching$`L-DLPFC--L-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`R-DLPFC--R-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`L-DLPFC--R-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`R-DLPFC--L-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`L-Par--R-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`L-Par--MidCing`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`R-Par--MidCing`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])


#Look at correlations per group
#for L-DLPFC--R-DLPFC
#TMT-A
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime2
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-DLPFC` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#for L-DLPFC--MidCing
#TMT-A
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `L-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `L-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `L-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `L-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `L-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `L-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `L-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime2
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#for R-DLPFC--MidCing
#TMT-A
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `R-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `R-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `R-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `R-DLPFC--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `R-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `R-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `R-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime2
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#for L-DLPFC--L-Par
#TMT-A
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `L-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `L-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #sig
cor.test(formula= ~ `L-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `L-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `L-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `L-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `L-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime2
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#for R-DLPFC--R-Par
#TMT-A
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `R-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `R-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `R-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `R-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `R-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `R-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `R-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `R-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #sig
#HayBTime2
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#for L-DLPFC--R-Par
#TMT-A
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `L-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `L-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `L-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `L-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `L-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `L-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `L-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime2
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-DLPFC--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#for R-DLPFC--L-Par
#TMT-A
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `R-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `R-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `R-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `R-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #sig
cor.test(formula= ~ `R-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `R-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `R-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime2
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-DLPFC--L-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#for L-Par--R-Par
#TMT-A
cor.test(formula= ~ `L-Par--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `L-Par--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `L-Par--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `L-Par--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `L-Par--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `L-Par--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `L-Par--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `L-Par--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `L-Par--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `L-Par--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `L-Par--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime2
cor.test(formula= ~ `L-Par--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `L-Par--R-Par` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `L-Par--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `L-Par--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--R-Par` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#for L-Par--MidCing
#TMT-A
cor.test(formula= ~ `L-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `L-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `L-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `L-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `L-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `L-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `L-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `L-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `L-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `L-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime2
cor.test(formula= ~ `L-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `L-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `L-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `L-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#for R-Par--MidCing
#TMT-A
cor.test(formula= ~ `R-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `R-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + TrailsA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#TMT-B
cor.test(formula= ~ `R-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + TrailsB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - ColourNaming
cor.test(formula= ~ `R-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + ColorNaming.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - WordReading
cor.test(formula= ~ `R-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `R-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + WordReading.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Stroop - Inhibition
cor.test(formula= ~ `R-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `R-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + Inhibition.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#LetFluency
cor.test(formula= ~ `R-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + LetFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#CatFluency
cor.test(formula= ~ `R-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `R-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + CatFluency.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#Switching
cor.test(formula= ~ `R-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("2")) #sig
cor.test(formula= ~ `R-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + Switching.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBTime1
cor.test(formula= ~ `R-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBTime1.Raw, data = Connectome_data, subset = Group == c("5")) #sig
#HayBTime2
cor.test(formula= ~ `R-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("1")) #sig
cor.test(formula= ~ `R-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBTime2.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatA
cor.test(formula= ~ `R-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("4")) #sig
cor.test(formula= ~ `R-Par--MidCing` + HayBCatA.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
#HayBCatB
cor.test(formula= ~ `R-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("1")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("2")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("3")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("4")) #not sig
cor.test(formula= ~ `R-Par--MidCing` + HayBCatB.Raw, data = Connectome_data, subset = Group == c("5")) #not sig
