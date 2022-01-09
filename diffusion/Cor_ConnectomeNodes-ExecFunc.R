#This script will analyse the correlation between the connectome node values and
#and the neuropsychological assesssment data. 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 03/01/22

#------------------------------Setting up--------------------------------------#
#install packages/open libraries
pacman::p_load(ggplot2, ppcor)

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
setwd('V:/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/FPN/small-nodes')

#collate connectome values into dataset
filenames <- list.files(path='V:/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/FPN/small-nodes')
r.nms <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
AllData <- lapply(filenames, read.table, header=FALSE, row.names=r.nms, col.names="V1")

for(i in sequence(length(filenames))){
  for (j in 1:10){
    small_node_data[i,j]<- AllData[[i]][[1]][j] 
  }
}

#combine neuropsych and connectome dataset
Connectome_data <- cbind(DPRC_neuropsych_data,small_node_data)

#run correlation tests between variables
#Left DLPFC -- Right DLPFC
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$ColorNaming.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$Inhibition.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$TrailsA.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$TrailsB.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$CatFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--R-DLPFC`, Connectome_data$Switching.Raw) #not sig
#Left DLPFC -- Mid Cingulate
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$HayBTime1.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$ColorNaming.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$Inhibition.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$TrailsA.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$TrailsB.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$CatFluency.Raw) #sig
cor.test(Connectome_data$`L-DLPFC--MidCing`, Connectome_data$Switching.Raw) #not sig
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
cor.test(Connectome_data$`L-DLPFC--L-Par`, Connectome_data$Switching.Raw) #not sig
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
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$Inhibition.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$TrailsA.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$TrailsB.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$LetFluency.Raw) #not sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$CatFluency.Raw) #sig
cor.test(Connectome_data$`L-Par--R-Par`, Connectome_data$Switching.Raw) #not sig
#Left DLPFC -- Mid Cingulate
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$HayBTime1.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$HayBTime2.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$HayBCatA.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$HayBCatB.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$ColorNaming.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$WordReading.Raw) #not sig
cor.test(Connectome_data$`L-Par--MidCing`, Connectome_data$Inhibition.Raw) #sig
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
pcor.test(Connectome_data_HayBCatA$`L-Par--R-Par`, Connectome_data_HayBCatA$HayBCatA.Raw, Connectome_data_HayBCatA[,c("Age", "Sex_binary")]) #sig
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
pcor.test(Connectome_data_Inhibition$`R-DLPFC--MidCing`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`L-DLPFC--L-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`R-DLPFC--R-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`L-DLPFC--R-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`R-DLPFC--L-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`L-Par--R-Par`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Inhibition$`L-Par--MidCing`, Connectome_data_Inhibition$Inhibition.Raw, Connectome_data_Inhibition[,c("Age", "Sex_binary")])
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
pcor.test(Connectome_data_CatFluency$`R-DLPFC--MidCing`, Connectome_data_CatFluency$CatFluency.Raw, Connectome_data_CatFluency[,c("Age", "Sex_binary")])
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
pcor.test(Connectome_data_Switching$`R-DLPFC--MidCing`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`L-DLPFC--L-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`R-DLPFC--R-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`L-DLPFC--R-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`R-DLPFC--L-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`L-Par--R-Par`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`L-Par--MidCing`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])
pcor.test(Connectome_data_Switching$`R-Par--MidCing`, Connectome_data_Switching$Switching.Raw, Connectome_data_Switching[,c("Age", "Sex_binary")])














