#Threshold structural connectivity (SC) connectome matrix by the number of 
#streamlines.



#set up directory where connectome matrix files are stored: 
#dk
#cross-sectional
input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered/all_files/'
#longitudinal
#input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/all_files/'
setwd(input_directory)

output_directory <-'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered/two_threshold_steps/test/thresholded_by_streamlines(step1)/'
#output_directory <-'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/two_threshold_steps/test/thresholded_by_streamlines(step1)/'

#place all csvs into one list (e.g., should be 227 files for cross-sectional cohort)
all_files_csvs <- lapply(list.files(input_directory), read.csv, header=FALSE)

#keep values that are .0001 (~100 streamlines) or more. Turn all other values which are smaller than this as a zero (0).
for (i in 1:length(all_files_csvs)) {
  #start with empty new matrix
  new_mat <- matrix(rep(0,42),nrow=84,ncol=84) #dk
  
  current_csv_file <-all_files_csvs[[i]] #select one csv file
  PAR_NAME <- substr(list.files()[i],4,18)#dk
  
  for (j in 1:length(current_csv_file)) {
    for (k in 1:length(current_csv_file)) {
      
      if (current_csv_file[j,k] >= .0001) {
        #print("hi")
        new_mat[j,k] <- current_csv_file[j,k]
      }
    }
  }
  
  #write new csv output  file with with thresholded matrix 
  setwd(output_directory)
  write.table(new_mat, paste('dk_', PAR_NAME, '_streamline_thresholded.csv',sep=""), sep = ",", row.names= FALSE, col.names = FALSE) #dk
  setwd(input_directory)
  
}



