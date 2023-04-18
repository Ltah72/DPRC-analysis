#compare participant hcp connectome files with sig. p-values


#set up director where p-value files are stored: 
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/stats_results/weighted/ConnectomeWhole_linear_trend_stats/'
setwd(directory)
#read in pvalue file
sig_pvalues_matrix <- read.csv("outputWhole_connectome_linear_trend_fwe_1mpvalue_t1.csv", skip=1, header=FALSE)


#set up directory where hcp files are stored:
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/all_files/'
setwd(directory)

#create a dataframe to record the edges and participants
df_columns <- c("edges","number_of_participants")
df_edges_participants_bad <- data.frame(matrix(nrow=0,ncol=length(df_columns)))
colnames(df_edges_participants_bad) <- df_columns

#create a matrix of unique node pairs. The formula is n(n-1)/2, where n is the 
#number of chosen nodes. 
hcp_nodes <- seq(from=1, to=379, by=1)
hcp_node_combos <- combn(hcp_nodes,2)

# #create recording vector outputs
# #edges
# edge_vec <- c();
# #number of participants
# par_number_vec <- c();

#compare each hcp participant file to the pvalue file, and make note of which values do not match. 
#Specifically, if the sig. pvalue is >= .95, while there is 0 value in the participant file in 
#the same location. 

files_all <- list.files(pattern="hcp") 


for (i in 1:length(files_all) {}

par_matrix_file <- read.csv(files_all[i], header=FALSE)

for (j in 1:length(par_matrix_file)) {
  for (k in 1:length(par_matrix_file)) {
   
     if(sig_pvalues_matrix[j,k] >= 0.95 & par_matrix_file[j,k] == 0) {
      
      #record the bad edge's index location and the number of participants who have this
      print("BAD!")
      
      df_edges_participants_bad[k, "edges"] <- as.character(paste(j,',',k))
      #assign a value (zero) to the number of participants at that new index
      df_edges_participants_bad[k,"number_of_participants"] <- 0
      df_edges_participants_bad[k,"number_of_participants"] <- df_edges_participants_bad[k,"number_of_participants"]+1
      
    }
    
  }
}


#record number of participants who have this 'bad edge'
if df_edges_participants_bad[1,"edges"]

df_edges_participants_bad[1,"number_of_participants"] <- +1


#testing
sig_pvalues_matrix[1,1] <- .95
i = 1
j = 1
k = 1

sig_pvalues_matrix[1,2] <- .95
par_matrix_file[1,2] <- 0
i = 1
j = 1
k = 2

k = 3

sig_pvalues_matrix[1,4] <- .95
par_matrix_file[1,4] <- 0
k = 4






