#Take the average of the connectome files. This is being done for the 
#interaction files (post-pre), but should applicable for other files too. 



#for C
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/hcpmmpFiles/weighted/interaction_study(F2-F0)/groups/C'
setwd(directory)

csvs <- lapply(list.files(directory), read.csv, header=FALSE)
C_avg_csvs<-Reduce("+", csvs) / length(csvs)
write.table(C_avg_csvs, 'C_avg_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)

#for SCD

#for aMCI

#for mMCI

#for AD


#for all (collapse across groups)
#directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/all_files/'
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered_thr/all_files/'
directory <- "V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered_thr/consensus_thresholded(step2)"

setwd(directory)

csvs <- lapply(list.files(directory), read.csv, header=FALSE)
all_avg_csvs<-Reduce("+", csvs) / length(csvs)
write.table(all_avg_csvs, 'all_avg_dk_connectome.csv', sep = ",", row.names= FALSE, col.names = FALSE)



