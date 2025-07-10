#After running a script to extract the interested values from the connectome 
#network, (create_FPN_connectome_files.R) and performing statistical tests to 
#compare groups (via MRtrix's vectorstats function), this script will now take 
#each of the significantly different values, and place them back into the 
#connectome file in order to visualise them. In other words, we will visualise 
#the interested connectome networks between our interested groups.  


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 21/01/21


#choose & set to your directory. This is where each of your participant's 
#values for the interested connectome network files should be. 
#setwd("/Connectome_test/FPN_connectome_csv_files/done/output")  

setwd("V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/stats_results/weighted/thresholded/C_SCD_threshold/FPN_stats/bigNode_group_diff_cov-age_sex")
#setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/stats_results/weighted/thresholded/C_SCD_threshold/FPN_stats/bigNode_linear_trend_cov-age_sex')

#read in the family-wise error (fwe) stat files
FPN_stats_connectome_fwe0 <- read.csv('bigNode_group_diff_cov-age_sex_fwe_1mpvalue_t3.csv', header = FALSE)
#FPN_stats_connectome_fwe0 <- read.csv('big-node_linear_trend_FPN_stats-Zstat_t1.csv', header = TRUE) 
#FPN_stats_connectome_fwe0 <- read.csv('output_FPN_2tail_stats-fwe_1mpvalue.csv', header = FALSE) #two-tailed (C != AD)
#FPN_stats_connectome_fwe1 <- read.csv('output_FPN_stats-fwe_1mpvalue_t1.csv', header = FALSE) #one-tailed (C > AD)
#FPN_stats_connectome_fwe2 <- read.csv('output_FPN_stats-fwe_1mpvalue_t2.csv', header = FALSE) #one-tailed (AD > C)


#use the same node pair as in the previous script (create_FPN_big_node_connectome_files.R) 
# - my chosen nodes are below:
#nodes <- c(" ")
nodes <- c(73,253,67,247,97,277,98,278,26,206,70,250,71,251,87,267,68,248,83,263,85,265,84,264,86,266,40,220,41,221,55,235,44,224,43,223,36,216,39,219,37,217,48,228,95,275,49,229,117,297,50,230,47,227,42,222,45,225,46,226,29,209,143,323,151,331,150,330,149,329,148,328,116,296,147,327,146,326,145,325,144,324)
#create a matrix of unique node pairs. The formula is n(n-1)/2, where n is the 
#number of chosen nodes. 
node_combos <- combn(nodes,2)

#Preallocate zeros to .csv files, per each participant. 
#In my case,I have a 379 x 379 matrix. 
connectome_matrix0 <- matrix(data=0, nrow=379, ncol=379)
#connectome_matrix1 <- matrix(data=0, nrow=379, ncol=379)
#connectome_matrix2 <- matrix(data=0, nrow=379, ncol=379)

#put connectome values from your stats output file into the connectome matrix, 
#using the coordinates from the node combination vector. 

for (i in 1:(length(node_combos)/2)) {
  
  current_pair1 <- node_combos[1,i]
  current_pair2 <- node_combos[2,i]
  
  node_value0 <- FPN_stats_connectome_fwe0[2,i]
  #node_value1 <- FPN_stats_connectome_fwe1[2,i]
  #node_value2 <- FPN_stats_connectome_fwe2[2,i]
  
  #fill both combinations of the pairs with the node value from the stats sheet. 
  connectome_matrix0[current_pair1, current_pair2] <- node_value0
  connectome_matrix0[current_pair2, current_pair1] <- node_value0
  
  #connectome_matrix1[current_pair1, current_pair2] <- node_value1
  #connectome_matrix1[current_pair2, current_pair1] <- node_value1
  
  #connectome_matrix2[current_pair1, current_pair2] <- node_value2
  #connectome_matrix2[current_pair2, current_pair1] <- node_value2
 
  
}

#write matrices as a .csv file. 
write.table(format(connectome_matrix0, digits=20), 'transferred_bigNode_group_diff_cov-age_sex_fwe_1mpvalue_t3.csv', row.names=FALSE, col.names=FALSE, quote=FALSE)
#write.table(connectome_matrix0, 'transferred_values_FPN_connectome_2tailed_fwe_1mpvalue.csv', row.names=FALSE, col.names=FALSE, quote=FALSE) #two-tailed (C != AD)
#write.table(connectome_matrix1, 'transferred_values_FPN_connectome_stats_fwe_1mpvalue_t1.csv', row.names=FALSE, col.names=FALSE, quote=FALSE) #one-tailed (C > AD)
#write.table(connectome_matrix2, 'transferred_values_FPN_connectome_stats_fwe_1mpvalue_t2.csv', row.names=FALSE, col.names=FALSE, quote=FALSE) #one-tailed (AD > C)






