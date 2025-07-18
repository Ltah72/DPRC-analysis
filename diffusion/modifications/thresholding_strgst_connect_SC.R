#Threshold structural connectivity (SC) connectome matrix by the strongest connections (e.g., 75% of the highest connections).

#set up directory where connectome matrix files are stored: 
#dk
input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered/all_files/matrix_files/'
#set working directory
setwd(input_directory)
#set up the output directory (where the thresholded connectome matrices will be stored)
output_directory <-'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered/two_threshold_steps/strng_connect(step1)&consensus(step2)/threshd_by_strng_connect(step1)/'
#place all csvs into one list (e.g., should be 227 files for cross-sectional cohort; 228 if you also want to include the sample avg file)
all_files_csvs <- lapply(list.files(input_directory), read.csv, header=FALSE)

#keep values that are the 75% strongest values (~) or more - discard the lowest 25% (~) or less. Turn all other values which are smaller than this set threshold as a zero (0).
for (i in 1:length(all_files_csvs)) {
  #start with empty new matrix filled with zeros
  new_mat <- matrix(rep(0,42),nrow=84,ncol=84) #dk
  
  #read in the current participant file
  current_csv_file <-all_files_csvs[[i]] #select one csv file
  PAR_NAME <- substr(list.files()[i],4,18)#dk
  
  #calculate the threshold value of participant dk connectome matrix file (e.g., lower 25% quantile)
  lower_values <- current_csv_file[lower.tri(current_csv_file, diag = FALSE)] #just take the values from one side (e.g,. lower half triangle) 
  numeric_values_no_zeros <- lower_values[lower_values != 0] #do not count the zeros
  threshold_25<-quantile(numeric_values_no_zeros , 0.25) #calculate the lower 25% quantile
  
  #apply changes to the new matrix
  for (j in 1:length(current_csv_file)) {
    for (k in 1:length(current_csv_file)) {
      
      if (current_csv_file[j,k] > threshold_25) {
        new_mat[j,k] <- current_csv_file[j,k]
      }
    }
  }
  
  #write new csv output  file with with thresholded matrix 
  setwd(output_directory)
  write.table(new_mat, paste('dk_', PAR_NAME, '_strng_connect_threshd.csv',sep=""), sep = ",", row.names= FALSE, col.names = FALSE) #dk
  setwd(input_directory)
}





