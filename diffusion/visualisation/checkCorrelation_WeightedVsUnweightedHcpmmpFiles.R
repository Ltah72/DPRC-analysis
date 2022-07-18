#Check correlation values between unweighted (x) vs. weighted files from hcpmmp 
#atlas



#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 04/05/22


# choose & set to your directory for unweighted first. This is where all of your participant 
#whole-brain connectome files should be. 
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/all_files/')  


files_all <- list.files()  #load the file names into the work space

for(i in sequence(length(files_all))){
  
  #extract participant name
  PAR_NAME = substr(files_all[i], 9, 23) 
  
  #fullname <- paste('FPN_BigNodeData_', PAR_NAME, '.csv', sep='')
  
  #read .csv files for every participant
  pt_data <- read.csv(files_all[i], header = FALSE)
  
  #put data into a vector
  unweighted_values_vec <- 