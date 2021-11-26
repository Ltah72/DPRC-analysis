#calculating the average for the group connectome files. This will allow us to
#display the .csv file for each group. 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 19/01/21


#for Control (C) group
# choose & set to your directory
setwd('V:/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/C')  

files_C <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_C))){
  Data_C <- read.csv(files_C[i], header = FALSE)
  Means_C <- apply(Data_C, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}

write.csv(Means_C, 'Control_avg_connectome.csv')

#for Subjective Cognitive Decline (SCD) group
# choose & set to your directory
setwd('V:/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/SCD')  

files_SCD <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_SCD))){
  Data_SCD <- read.csv(files_SCD[i], header = FALSE)
  Means_SCD <- apply(Data_SCD, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}

write.csv(Means_SCD, 'SCD_avg_connectome.csv')


#for Amnestic Mild Cognitive Impairment (aMCI) group
# choose & set to your directory
setwd('V:/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/aMCI')  

files_aMCI <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_aMCI))){
  Data_aMCI <- read.csv(files_aMCI[i], header = FALSE)
  Means_aMCI <- apply(Data_aMCI, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}

write.csv(Means_aMCI, 'aMCI_avg_connectome.csv')


#for Multiple-Domain Mild Cognitive Impairment group
# choose & set to your directory
setwd('V:/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/mMCI')  

files_mMCI <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_mMCI))){
  Data_mMCI <- read.csv(files_mMCI[i], header = FALSE)
  Means_mMCI <- apply(Data_mMCI, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}

write.csv(Means_mMCI, 'mMCI_avg_connectome.csv')


#for Alzheimer's Disease (AD) group
# choose & set to your directory
setwd('V:/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/AD')  

files_AD <- list.files()  #load the file names into the work space

for(i in sequence(length(files_AD))){
  Data_AD <- read.csv(files_AD[i], header = FALSE)
  Means_AD <- apply(Data_AD, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}

write.csv(Means_AD, 'AD_avg_connectome.csv')

