#This script will perform two main analysis. It will analyse the relationship 
#between the fMRI timeseries from track end points (the anterior and the 
#posterior end points) of the superior longitudinal fasciculus (SLF) - SLF1, 
#SLF2, and SLF3. 

#This script will also analyse the superior longitudinal fasciculus (SLF) white 
#matter tracts and the fMRI timeseries from the tract end points.


#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 


#------------------------------Setting up--------------------------------------#
#install packages/open libraries
pacman::p_load(dplyr, ggplot2, ppcor)

#first read in the covariates group data file: 
#setwd('/yourpathway/')
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')




#-----------------------Cross-sectional analysis-------------------------------#


DPRC_neuropsych_data <- read.csv("connectome_DPRC_neuropsych_data_lined_up_valid_participants.csv")



#navigate to the correct pathway which contains the fMRI timeseries text files: 
setwd('V:/NECTAR_data/LENORE/derivatives/fMRI_denoised/unzipped_files')

#read in text file of data
fMRI_SLF_data <- cbind.data.frame(read.table("Ptpt_corr_pval.txt", header = T)) 



                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
