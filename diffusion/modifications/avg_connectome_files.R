#Calculating the whole-brain average connectome and specific frontoparietal 
#network (FPN) for the group connectome files. This will allow us to display the 
#.csv file for each group. These files can then be visualised through Matlab as 
#a connectome matrix and/or MNE Python as a circular plot. 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 19/01/21



#define FPN nodes - and put in order
nodes <- sort(c(73,253,67,247,97,277,98,278,26,206,70,250,71,251,87,267,68,248,83,263,85,265,84,264,86,266,40,220,41,221,55,235,44,224,43,223,36,216,39,219,37,217,48,228,95,275,49,229,117,297,50,230,47,227,42,222,45,225,46,226,29,209,143,323,151,331,150,330,149,329,148,328,116,296,147,327,146,326,145,325,144,324))


#for Control (C) group
# choose & set to your directory
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/C')  

files_C <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_C))){
  Data_C <- read.csv(files_C[i], header = FALSE)
  Means_C <- apply(Data_C, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}
write.table(Means_C, 'Control_avg_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#for circle plot (360 nodes), remove rows and columns from 361-379.
Means_C_360nodes <- Means_C[1:360, 1:360]
write.table(Means_C_360nodes, 'Control_avg_connectome_360nodes.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#delete rows/columns from FPN nodes for comparing FPN averages
Means_FPN_C <-Means_C[(nodes), (nodes)]
write.table(Means_FPN_C, 'Control_avg_FPN_connectome.csv', sep = ",", row.names = FALSE, col.names = FALSE)


#for Subjective Cognitive Decline (SCD) group
# choose & set to your directory
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/SCD')  

files_SCD <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_SCD))){
  Data_SCD <- read.csv(files_SCD[i], header = FALSE)
  Means_SCD <- apply(Data_SCD, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}
write.table(Means_SCD, 'SCD_avg_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#for circle plot (360 nodes), remove rows and columns from 361-379.
Means_SCD_360nodes <- Means_SCD[1:360, 1:360]
write.table(Means_SCD_360nodes, 'SCD_avg_connectome_360nodes.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#delete rows/columns from FPN nodes for comparing FPN averages
Means_FPN_SCD <-Means_SCD[(nodes), (nodes)]
write.table(Means_FPN_SCD, 'SCD_avg_FPN_connectome.csv', sep = ",", row.names = FALSE, col.names = FALSE)

#for Amnestic Mild Cognitive Impairment (aMCI) group
# choose & set to your directory
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/aMCI')  

files_aMCI <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_aMCI))){
  Data_aMCI <- read.csv(files_aMCI[i], header = FALSE)
  Means_aMCI <- apply(Data_aMCI, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}
write.table(Means_aMCI, 'aMCI_avg_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#for circle plot (360 nodes), remove rows and columns from 361-379.
Means_aMCI_360nodes <- Means_aMCI[1:360, 1:360]
write.table(Means_aMCI_360nodes, 'aMCI_avg_connectome_360nodes.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#delete rows/columns from FPN nodes for comparing FPN averages
Means_FPN_aMCI <-Means_aMCI[(nodes), (nodes)]
write.table(Means_FPN_aMCI, 'aMCI_avg_FPN_connectome.csv', sep = ",", row.names = FALSE, col.names = FALSE)

#for Multiple-Domain Mild Cognitive Impairment group
# choose & set to your directory
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/mMCI')  

files_mMCI <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_mMCI))){
  Data_mMCI <- read.csv(files_mMCI[i], header = FALSE)
  Means_mMCI <- apply(Data_mMCI, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}
write.table(Means_mMCI, 'mMCI_avg_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#for circle plot (360 nodes), remove rows and columns from 361-379.
Means_mMCI_360nodes <- Means_mMCI[1:360, 1:360]
write.table(Means_mMCI_360nodes, 'mMCI_avg_connectome_360nodes.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#delete rows/columns from FPN nodes for comparing FPN averages
Means_FPN_mMCI <-Means_mMCI[(nodes), (nodes)]
write.table(Means_FPN_mMCI, 'mMCI_avg_FPN_connectome.csv', sep = ",", row.names = FALSE, col.names = FALSE)

#for Alzheimer's Disease (AD) group
# choose & set to your directory
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/AD')  

files_AD <- list.files()  #load the file names into the work space

for(i in sequence(length(files_AD))){
  Data_AD <- read.csv(files_AD[i], header = FALSE)
  Means_AD <- apply(Data_AD, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}
write.table(Means_AD, 'AD_avg_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#for circle plot (360 nodes), remove rows and columns from 361-379.
Means_AD_360nodes <- Means_AD[1:360, 1:360]
write.table(Means_AD_360nodes, 'AD_avg_connectome_360nodes.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#delete rows/columns from FPN nodes for comparing FPN averages
Means_FPN_AD <-Means_AD[(nodes), (nodes)]
write.table(Means_FPN_AD, 'AD_avg_FPN_connectome.csv', sep = ",", row.names = FALSE, col.names = FALSE)


#for Group average (FBA population template cohort)
# choose & set to your directory
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/group_template')  

files_G <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_G))){
  Data_G <- read.csv(files_G[i], header = FALSE)
  Means_G <- apply(Data_G, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}
write.table(Means_G, 'Group_avg_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#for circle plot (360 nodes), remove rows and columns from 361-379.
Means_G_360nodes <- Means_G[1:360, 1:360]
write.table(Means_G_360nodes, 'Group_avg_connectome_360nodes.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#delete rows/columns from FPN nodes for comparing FPN averages
Means_FPN_G <-Means_G[(nodes), (nodes)]
write.table(Means_FPN_G, 'Group_avg_FPN_connectome.csv', sep = ",", row.names = FALSE, col.names = FALSE)

