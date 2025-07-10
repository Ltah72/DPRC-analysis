#Count number of zeros in participants (e.g., C and SCD groups) structural 
#connectome matrices


#for cross-sectional study, C + SCD groups, n = 95; or just C, n = 35
#for longitudinal study, C + SCD groups, n = 61 (F0); or just C, n = 21 (F0)

#set up directory where C and SCD files are stored: 
#dk
#directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered/C&SCD_groups/'
#directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/C&SCD_F0/'
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered_thr/C&SCD_groups/'
#directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/C&SCD_F0/'

#hcp
#directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/C&SCD_groups/'
#directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/hcpmmpFiles/weighted/C&SCD_F0/'
setwd(directory)

#place all csvs into one list (e.g., should be 95 files with C and SCD group)
threshold_group_csvs <- lapply(list.files(directory), read.csv, header=FALSE)
#empty matrix
num_mat <- matrix(rep(0,42),nrow=84,ncol=84) #dk
#num_mat <- matrix(rep(0,190),nrow=380,ncol=380) #hcp

#count number of zeros with a for loop here
for (i in 1:length(threshold_group_csvs)) {
  current_csv_file <-threshold_group_csvs[[i]] #select one csv file
  
  for (j in 1:length(current_csv_file)) {
    for (k in 1:length(current_csv_file)) {
      
      if (current_csv_file[j,k] == 0) {
        
        num_mat[j,k] <- num_mat[j,k]+1
      }
    }
  }
}

#write csv file with the total number of zeros in each element of the matrix
write.table(num_mat, 'number_of_zeros_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)
#write csv file as a proportion of the number of zeros in each element of the matrix as well      
prop_mat <- (num_mat/95)*100 #C&SCD group cross-sectional
#prop_mat <- (num_mat/35)*100 #C group cross-sectional
#prop_mat <- (num_mat/61)*100 #C&SCD group longitudinal
#prop_mat <- (num_mat/21)*100 # C group longitudinal

write.table(prop_mat, 'prop_of_zeros_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)



#------------------------------------------------------------------------------#
#now, we want to create the thresholded matrices of each participant based on our findings
#dk
#input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered/all_files/'
#input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/all_files/'
input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered_thr/all_files/'
#input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/two_threshold_steps/thresholded_by_streamlines(step1)/'

#hcp
#input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/all_files/'
#input_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/hcpmmpFiles/weighted/all_files/'
setwd(input_directory)

#where to store the thresholded participant files
#dk
#output_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered/thresholded_connectomes/'
#output_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/thresholded_connectomes/'
output_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered_thr/consensus_thresholded(step2)/'
#output_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/two_threshold_steps/thresholded_by_consensus(step2)/'

#hcp
#output_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/thresholded_connectomes/'
#output_directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/hcpmmpFiles/weighted/thresholded_connectomes/'

#place all csvs into one list (e.g., should be all 227 files for cross-sectional study with all groups)
all_csvs <- lapply(list.files(input_directory), read.csv, header=FALSE)


#for loop here - create newly thresholded matrices for each participant
for (i in 1:length(all_csvs)) {
  setwd(input_directory)
  current_csv_file <- all_csvs[[i]] #select one csv file
  
  PAR_NAME <- substr(list.files()[i],8,22)#dk wgt_str_thr
  #PAR_NAME <- substr(list.files()[i],4,18)#dk
  #PAR_NAME <- substr(list.files()[i],10,24)#hcp
   
  for (j in 1:length(current_csv_file)) {
    for (k in 1:length(current_csv_file)) {
      
      #50% threshold
      if (prop_mat[j,k] >= 50) {
        current_csv_file[j,k] <- 0
      }
      
    }
  }

  setwd(output_directory)
  #write.table(current_csv_file, paste('dk_', PAR_NAME, '_thresholded.csv',sep=""), sep = ",", row.names= FALSE, col.names = FALSE) #dk
  #write.table(current_csv_file, paste('hcpmmp_', PAR_NAME, '_thresholded.csv',sep=""), sep = ",", row.names= FALSE, col.names = FALSE) #hcp
  write.table(current_csv_file, paste('dk_', PAR_NAME, '_wgt_strmln&cns_thr.csv',sep=""), sep = ",", row.names= FALSE, col.names = FALSE) #dk
  
}

  





