#This script takes the interested values (or nodes) based upon the interested 
#network. These values are extracted from each of the participant's connectome 
#file and placed into another vector file to conduct statistical tests on them. 
#And so, after running this script, you will then need to run MRtrix's 
#'vectorstats' on the newly generated excel files. 


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 21/01/21


# choose & set to your directory. This is where all of your participant 
#whole-brain connectome files should be. 
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/all_files/')  

#calculate your node pairs - my chosen nodes are below:
#nodes <- c(" ")
nodes <- c(73,253,67,247,97,277,98,278,26,206,70,250,71,251,87,267,68,248,83,263,85,265,84,264,86,266,40,220,41,221,55,235,44,224,43,223,36,216,39,219,37,217,48,228,95,275,49,229,117,297,50,230,47,227,42,222,45,225,46,226,29,209,143,323,151,331,150,330,149,329,148,328,116,296,147,327,146,326,145,325,144,324)
#create a matrix of unique node pairs. The formula is n(n-1)/2, where n is the 
#number of chosen nodes. 
node_combos <- combn(nodes,2)


files_all <- list.files()  #load the file names into the work space

for(i in sequence(length(files_all))){
  
  #extract participant name
  PAR_NAME = substr(files_all[i], 10, 24) 
  
  fullname <- paste('FPN_BigNodeData_', PAR_NAME, '.csv', sep='')
  
  pt_data <- read.csv(files_all[i], header = FALSE)
  
  #extract the cells for each unique node pair
  for (j in 1:(length(node_combos)/2)) {
    
    current_pair1 <- node_combos[1,j]
    current_pair2 <- node_combos[2,j]
  
    node_data <- pt_data[current_pair1, current_pair2]
    
    write.table(node_data, fullname, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  
  }
}