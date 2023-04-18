#Modify .mat fMRI ROI connectivity files to have 379 nodes (surface-space) or 
#426 nodes (volumetric-space) for hcpmmp1 parcellation atlas. After organising 
#the files, you can run the 'create_FPN_big_node_connectome_file.R to extract 
#just the FPN node values (82 nodes, 3,321 edges).


pacman::p_load(R.matlab)

#go into working directory with all participant ROI.mat files from the CONN analysis 
#setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/fMRI_data/connectome/data/fMRI_ROI_mat_files/surface-space/C3_post/')  
setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/fMRI_data/connectome/data/fMRI_ROI_mat_files/volumetric-space/')  

#add any necessary sources: 
source('H:/ltah262/PhD/Diffusion/script/dprc/neuropsych/insertRow.R')
source('H:/ltah262/PhD/Diffusion/script/dprc/neuropsych/insertColumn.R')



files_all <- list.files()  #load the file names into the work space

#create a vector of NaN values (to be inserted into the matrix, later, 
#for the missing node)
#surface-space
#vector_NaNs1 <- rep(NaN,378) #missing 1 ROI for some reason...
#vector_NaNs2 <- rep(NaN,379)
#volumetric-space
vector_NaNs1 <- rep(NaN,426)
vector_NaNs2 <- rep(NaN,426)


for(i in sequence(length(files_all))){
  
  PAR_NAME = substr(files_all[i], 12, 21) 
  
  #this will be read in as a list
  data <- readMat(paste('resultsROI_', PAR_NAME, '_Condition001.mat',sep=''))
  #modify the 'Z' values from the list  
  Z_data <- data[["Z"]]
  #extract rows and columns 149 - 526 (378 x 378 matrix)
  #hcp_data <- Z_data[149:526,149:526]#surface
  hcp_data <- Z_data[559:984,559:984]#volume
  
  #convert to dataframe
  hcp_data <- as.data.frame(hcp_data)
  hcp_complete_data <- as.data.frame(hcp_data)
  #insert row and column 300 with NaN values (for the missing node 300 - hcpmmp1.R_H) - needed for the surface-space analysis
  #hcp_complete_data<-insertRow(hcp_data,vector_NaNs1,300) #replace row
  #hcp_complete_data<-insertColumn(hcp_complete_data,vector_NaNs2,300) #replace column
  
  #set working directory output folder
  #setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/fMRI_data/connectome/data/hcpmmpFiles/surface-space/') #surface  
  setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/fMRI_data/connectome/data/hcpmmpFiles/volumetric-space/') #volume 
  
  #write complete hcp file to .csv
  fullname <- paste('hcpmmp1w_', PAR_NAME, 'conC1.csv', sep='') #to fit 13 characters, so consistent with structural connectivity analysis
  #write.table(hcp_complete_data, fullname, append=TRUE, row.names=FALSE, col.names=FALSE) #surface
  write.table(hcp_complete_data, fullname, append=TRUE, row.names=FALSE, col.names=FALSE,sep=',') #volume
  
  #reset working directory to input folder
  #setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/fMRI_data/connectome/data/fMRI_ROI_mat_files/surface-space/C2_pre/') #surface  
  setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/fMRI_data/connectome/data/fMRI_ROI_mat_files/volumetric-space/') #volume 
}

  
  
  
  