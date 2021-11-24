#This script will choose certain interested nodes from a whole-brain connectome 
#file, and will then take the average of those nodes as a certain cluster area. 
#For this analysis, I am interested in the Frontoparietal network (FPN), and I 
#have identified 5 regions of interest (Left Frontal, Right Frontal, 
#Mid-cingulate, Left Parietal, and Right Parietal).


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 24/11/21


# choose & set to your directory. This is where all of your participant 
#whole-brain connectome files should be. 
#setwd('/yourpathway/')
setwd('V:/NECTAR_data/LENORE/test/connectome_test/')  

#define the 5 regions' nodes of interest
nodes_LFront <- c(73,67,97,98,26,70,71,87,68,83,85,84,86)
nodes_RFront <- c(253,246,277,278,206,250,251,267,248,263,265,264,266)
nodes_MidCing <- c(40,41,55,44,43,36,39,37,220,221,235,224,223,216,219,217)
nodes_LParietal <- c(48,95,49,117,50,47,42,45,46,29,143,151,150,149,148,116,147,146,145,144)
nodes_RParietal <- c(228,275,229,297,230,227,222,225,226,209,323,331,330,329,328,296,327,326,325,324)

#create unique combination between the regions of interest
Front_L_Front_R <- expand.grid(nodes_LFront, nodes_RFront)
Front_L_MidCing <- expand.grid(nodes_LFront, nodes_MidCing)
Front_R_MidCing <- expand.grid(nodes_RFront, nodes_MidCing)
Front_L_Par_L <- expand.grid(nodes_LFront, nodes_LParietal)
Front_R_Par_R <- expand.grid(nodes_RFront, nodes_RParietal)
Front_L_Par_R <- expand.grid(nodes_LFront, nodes_RParietal)
Front_R_Par_L <- expand.grid(nodes_RFront, nodes_LParietal)
Par_L_Par_R <- expand.grid(nodes_LParietal, nodes_RParietal)
Par_L_MidCing <- expand.grid(nodes_LParietal, nodes_MidCing)
Par_R_MidCing <- expand.grid(nodes_RParietal, nodes_MidCing)

#load the file names into the work space
files_all <- list.files() 

for(i in sequence(length(files_all))){
  
  #extract participant name
  PAR_NAME = substr(files_all[i], 9, 23) 
  fullname <- paste('FPN_SmallNodeData_', PAR_NAME, '.csv', sep='')
  file_data <- read.csv(paste('hcpmmp1_', PAR_NAME, '.csv', sep=''),  header=FALSE)
  
  #create empty vectors and matrices to put values in
  Front_L_Front_R_vec <- vector()
  Front_L_MidCing_vec <- vector()
  Front_R_MidCing_vec <- vector()
  Front_L_Par_L_vec <- vector()
  Front_R_Par_R_vec <- vector()
  Front_L_Par_R_vec <- vector()
  Front_R_Par_L_vec <- vector()
  Par_L_Par_R_vec <- vector()
  Par_L_MidCing_vec <- vector()
  Par_R_MidCing_vec <- vector()
  
  #create 10 x 1 empty matrix
  avg_values_mat <- matrix(, nrow=10, ncol=1)

  #extract the values from each node pair for each region of interest
  for (j in 1:(nrow(Front_L_Front_R))){
    
    current_node1 <- Front_L_Front_R[j,1]
    current_node2 <- Front_L_Front_R[j,2]
    
    #add values to vector
    Front_L_Front_R_vec[j] <- file_data[current_node1,current_node2]
  }
  for (j in 1:(nrow(Front_L_MidCing))){
    
    current_node1 <- Front_L_MidCing[j,1]
    current_node2 <- Front_L_MidCing[j,2]
    
    #add values to vector
    Front_L_MidCing_vec[j] <- file_data[current_node1,current_node2]
  }
  for (j in 1:(nrow(Front_R_MidCing))){
    
    current_node1 <- Front_R_MidCing[j,1]
    current_node2 <- Front_R_MidCing[j,2]
    
    #add values to vector
    Front_R_MidCing_vec[j] <- file_data[current_node1,current_node2]
  }
  for (j in 1:(nrow(Front_L_Par_L))){
    
    current_node1 <- Front_L_Par_L[j,1]
    current_node2 <- Front_L_Par_L[j,2]
    
    #add values to vector
    Front_L_Par_L_vec[j] <- file_data[current_node1,current_node2]
  }
  for (j in 1:(nrow(Front_R_Par_R))){
    
    current_node1 <- Front_R_Par_R[j,1]
    current_node2 <- Front_R_Par_R[j,2]
    
    #add values to vector
    Front_R_Par_R_vec[j] <- file_data[current_node1,current_node2]
  }
  for (j in 1:(nrow(Front_L_Par_R))){
    
    current_node1 <- Front_L_Par_R[j,1]
    current_node2 <- Front_L_Par_R[j,2]
    
    #add values to vector
    Front_L_Par_R_vec[j] <- file_data[current_node1,current_node2]
  }
  for (j in 1:(nrow(Front_R_Par_L))){
    
    current_node1 <- Front_R_Par_L[j,1]
    current_node2 <- Front_R_Par_L[j,2]
    
    #add values to vector
    Front_R_Par_L_vec[j] <- file_data[current_node1,current_node2]
  }
  for (j in 1:(nrow(Par_L_Par_R))){
    
    current_node1 <- Par_L_Par_R[j,1]
    current_node2 <- Par_L_Par_R[j,2]
    
    #add values to vector
    Par_L_Par_R_vec[j] <- file_data[current_node1,current_node2]
  }
  for (j in 1:(nrow(Par_L_MidCing))){
    
    current_node1 <- Par_L_MidCing[j,1]
    current_node2 <- Par_L_MidCing[j,2]
    
    #add values to vector
    Par_L_MidCing_vec[j] <- file_data[current_node1,current_node2]
  }
  for (j in 1:(nrow(Par_R_MidCing))){
    
    current_node1 <- Par_R_MidCing[j,1]
    current_node2 <- Par_R_MidCing[j,2]
    
    #add values to vector
    Par_R_MidCing_vec[j] <- file_data[current_node1,current_node2]
  }
  
  #calculate averages of all of the vectors and put into one vector
  avg_values_mat[1,1] <- mean(Front_L_Front_R_vec)
  avg_values_mat[2,1] <- mean(Front_L_MidCing_vec) 
  avg_values_mat[3,1] <- mean(Front_R_MidCing_vec) 
  avg_values_mat[4,1] <- mean(Front_L_Par_L_vec) 
  avg_values_mat[5,1] <- mean(Front_R_Par_R_vec) 
  avg_values_mat[6,1] <- mean(Front_L_Par_R_vec) 
  avg_values_mat[7,1] <- mean(Front_R_Par_L_vec) 
  avg_values_mat[8,1] <- mean(Par_L_Par_R_vec) 
  avg_values_mat[9,1] <- mean(Par_L_MidCing_vec) 
  avg_values_mat[10,1] <- mean(Par_R_MidCing_vec) 
  
  #write these average values into a .csv file for every participant
  write.table(avg_values_mat, fullname, append=TRUE, row.names=FALSE, col.names=FALSE)
  
}
