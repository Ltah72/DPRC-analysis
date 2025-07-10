#subtract connectome matrices to create interaction files (F2 - F0).


#specify directory
input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/two_threshold_steps/thresholded_by_consensus(step2)/'
#input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/FPN_Files/weighted/big-nodes/'
setwd(input_directory)


#place all csvs into one list (should be 95 files with C and SCD group)
all_csvs_files <- lapply(list.files(input_directory), read.csv, header=FALSE)

#directory where the interaction files (e.g., F2 - F0) will be outputted to
output_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/two_threshold_steps_interaction/'
#output_directory<- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/FPN_Files/weighted/interaction_study(F2-F0)'

#set up indices for the for loop below
j = 1
k = 2
#for loop here or lapply
for (i in 1:length(all_csvs_files)/2) {
  #for (j in 2:length(all_csvs_files)/2) {
  setwd(input_directory)
  current_csv_file_F0 <-all_csvs_files[[j]] #select one csv file
  current_csv_file_F2 <-all_csvs_files[[k]] #select the other csv file
  current_diff_csv_file <- current_csv_file_F2 - current_csv_file_F0 #subtract them!
  PAR_NAME <- substr(list.files()[j],4,16) #for whole connectome
  #PAR_NAME <- substr(list.files()[j],17,29) #for FPN connectome
  
  #write out the diff csv file
  setwd(output_directory)
  write.table(current_diff_csv_file, paste('dk_', PAR_NAME, '_F2-F0diff.csv',sep=""), sep = ",", row.names= FALSE, col.names = FALSE)
  #write.table(current_diff_csv_file, paste('FPN_BigNodeData_', PAR_NAME, '_F2-F0diff.csv',sep=""), sep = ",", row.names= FALSE, col.names = FALSE)
  
  j = j + 2 #add correct indexing to j and k to maintain subtraction pairing
  k = k + 2 
}

#this will result in an error at the end of the loop, but oh well, you should get all of your needed files. 





