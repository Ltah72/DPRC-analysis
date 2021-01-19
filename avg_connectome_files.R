#calculating the average for the group connectome files. This will allow us to
#display the .csv file for each group. 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 19/01/21


#for control group

# choose & set to your directory
setwd("/Connectome_test/pt_connectome_csv_files/C")  

files_C <- list.files()  #load the file names into the workspace

for(i in sequence(length(files_C))){
  Data_C <- read.csv(files_C[i], header = FALSE)
  Means_C <- apply(Data_C, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}

write.csv(Means_C, 'Control_avg_connectome.csv')


#for AD group
# choose & set to your directory
setwd("/Connectome_test/pt_connectome_csv_files/AD")  

files_AD <- list.files()  #load the file names into the work space

for(i in sequence(length(files_AD))){
  Data_AD <- read.csv(files_AD[i], header = FALSE)
  Means_AD <- apply(Data_AD, c(1,2), mean)
  #save your means in some meaningful way from each csv. 
}

write.csv(Means_AD, 'AD_avg_connectome.csv')

