#Organise DPRC neuropsychological data into one spreadsheet, containing all 
#participant values. 

#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 24/01/21


#install package manager to handle all other package installations and dependencies. 
install.packages("pacman")
pacman::p_load(xlsx, stringr, data.table)

#load in libraries
library("xlsx")
library("stringr")
library("data.table") 

# choose & set to your directory. This is where each of your participant's 
#files should be. 
setwd("N:/DPRC/Neuropsych Summary Files/Participant Files/")  
#setwd("/Participant Files/")  

#load the file names from directory into the work space
files_all <- list.files() 

#set up the neuropsych matrix headers to carry all participant data
neuropsych_matrixF0 <- matrix(nrow = length(files_all), ncol=134)
neuropsych_matrixF1 <- matrix(nrow = length(files_all), ncol=134)
neuropsych_matrixF2 <- matrix(nrow = length(files_all), ncol=134)
neuropsych_matrixF3 <- matrix(nrow = length(files_all), ncol=134)
neuropsych_matrixF4 <- matrix(nrow = length(files_all), ncol=134)

#set up column names (134 variables total) - need to add clin premorbid scores? 
colnames(neuropsych_matrixF0) <- c('ParticipantID', 'Age', 'EduLevel', 'Group', 'Ethnicity', 'Neuropsych_TestDate', 'TOPF-Raw', 'TOPF-Z', 'TOPF-Scaled', 'ClinPremorbid-z', 'ClinPremorbid-Scaled', 'ClinPremorbid-Cutoff', 'DSF-Raw', 'DSF-Z', 'DSF-Scaled', 'DSB-Raw', 'DSB-Z', 'DSB-Scaled', 'TrailsA-Raw', 'TrailsA-Z', 'TrailsA-Scaled', 'TrailsB-Raw', 'TrailsB-Z', 'TrialsB-Scaled', 'Coding-Raw', 'Coding-Z', 'Coding-Scaled', 'ColorNaming-Raw', 'ColorNaming-Z', 'ColorNaming-Scaled', 'WordReading-Raw', 'WordReading-Z', 'WordReading-Scaled', 'Inhibition-Raw', 'Inhibition-Z', 'Inhibition-Scaled', 'SYDBATNaming-Raw', 'SYDBATNaming-Z', 'SYDBATNaming-Scaled', 'SYDBATComp-Raw', 'SYDBATComp-Z', 'SYDBATComp-Scaled', 'SYDBATSemAss-Raw', 'SYDBATSemAss-Z', 'SYDBATSemAss-Scaled', 'BNTManual-Raw', 'BNTManual-Z', 'BNTManual-Scaled', 'BNTIvnik-Raw', 'BNTIvnik-Z', 'BNTIvnik-Scaled', 'Similarities-Raw', 'Similarities-Z', 'Similarities-Scaled', 'LineO-Raw', 'LineO-Z', 'LineO-Scaled', 'BD-Raw', 'BD-Z', 'BD-Scaled', 'Matrix-Raw', 'Matrix-Z', 'Matrix-Scaled', 'LM_I-Raw', 'LM_I-Z', 'LM_I-Scaled', 'LM_II-Raw', 'LM_II-Z', 'LM_II-Scaled', 'StoryAI-Raw', 'StoryAI-Z', 'StoryAI-Scaled', 'StoryAII-Raw', 'StoryAII-Z', 'StoryAII-Scaled',  'CVLT-II_Total-Raw', 'CVLT-II_Total-Z', 'CVLT-II_Total-Scaled', 'CVLT-II_Short-Raw', 'CVLT-II_Short-Z', 'CVLT-II_Short-Scaled', 'CVLT-II_Long-Raw', 'CVLT-II_Long-Z', 'CVLT-II_Long-Scaled', 'CVLT-II_RecogD-Raw', 'CVLT-II_RecogD-Z', 'CVLT-II_RecogD-Scaled', 'RCFT_lmm-Raw', 'RCFT_lmm-Z', 'RCFT_lmm-Scaled', 'RCFT_Del-Raw', 'RCFT_Del-Z', 'RCFT_Del-Scaled', 'RCFT_Recog-Raw', 'RCFT_Recog-Z', 'RCFT_Recog-Scaled', 'BVMT-R_Total-Raw', 'BVMT-R_Total-Z', 'BVMT-R_Total-Scaled', 'BVMT-R_Total-zGale', 'BVMT-R_Del-Raw', 'BVMT-R_Del-Z', 'BVMT-R_Del-Scaled', 'BVMT-R_Del-zGale', 'BVMT-R_RecogDis-Raw', 'BVMT-R_RecogDis-Z', 'BVMT-R_RecogDis-Scaled', 'BVMT-R_RecogDis-zGale', 'RCFT_Copy-Raw', 'RCFT_Copy-Z', 'RCFT-Copy-Scaled', 'LetFluency-Raw', 'LetFluency-Z', 'LetFluency-Scaled', 'CatFluency-Raw', 'CatFluency-Z', 'CatFluency-Scaled', 'Fluency-Raw', 'Fluency-Z', 'Fluency-Scaled', 'C/W_Inhib-Raw', 'C/W_Inhib-Z', 'C/W_Inhib-Scaled', 'Switching-Raw', 'Switching-z', 'Switching-Scaled', 'HayBTime1-Raw', 'HayBTime1-z', 'HayBTime2-Raw', 'HayBTime2-z', 'HayBCatA-Raw', 'HayBCatA-z', 'HayBCatB-Raw', 'HayBCatB-z')
colnames(neuropsych_matrixF1) <- c('ParticipantID', 'Age', 'EduLevel', 'Group', 'Ethnicity', 'Neuropsych_TestDate', 'TOPF-Raw', 'TOPF-Z', 'TOPF-Scaled', 'ClinPremorbid-z', 'ClinPremorbid-Scaled', 'ClinPremorbid-Cutoff', 'DSF-Raw', 'DSF-Z', 'DSF-Scaled', 'DSB-Raw', 'DSB-Z', 'DSB-Scaled', 'TrailsA-Raw', 'TrailsA-Z', 'TrailsA-Scaled', 'TrailsB-Raw', 'TrailsB-Z', 'TrialsB-Scaled', 'Coding-Raw', 'Coding-Z', 'Coding-Scaled', 'ColorNaming-Raw', 'ColorNaming-Z', 'ColorNaming-Scaled', 'WordReading-Raw', 'WordReading-Z', 'WordReading-Scaled', 'Inhibition-Raw', 'Inhibition-Z', 'Inhibition-Scaled', 'SYDBATNaming-Raw', 'SYDBATNaming-Z', 'SYDBATNaming-Scaled', 'SYDBATComp-Raw', 'SYDBATComp-Z', 'SYDBATComp-Scaled', 'SYDBATSemAss-Raw', 'SYDBATSemAss-Z', 'SYDBATSemAss-Scaled', 'BNTManual-Raw', 'BNTManual-Z', 'BNTManual-Scaled', 'BNTIvnik-Raw', 'BNTIvnik-Z', 'BNTIvnik-Scaled', 'Similarities-Raw', 'Similarities-Z', 'Similarities-Scaled', 'LineO-Raw', 'LineO-Z', 'LineO-Scaled', 'BD-Raw', 'BD-Z', 'BD-Scaled', 'Matrix-Raw', 'Matrix-Z', 'Matrix-Scaled', 'LM_I-Raw', 'LM_I-Z', 'LM_I-Scaled', 'LM_II-Raw', 'LM_II-Z', 'LM_II-Scaled', 'StoryAI-Raw', 'StoryAI-Z', 'StoryAI-Scaled', 'StoryAII-Raw', 'StoryAII-Z', 'StoryAII-Scaled',  'CVLT-II_Total-Raw', 'CVLT-II_Total-Z', 'CVLT-II_Total-Scaled', 'CVLT-II_Short-Raw', 'CVLT-II_Short-Z', 'CVLT-II_Short-Scaled', 'CVLT-II_Long-Raw', 'CVLT-II_Long-Z', 'CVLT-II_Long-Scaled', 'CVLT-II_RecogD-Raw', 'CVLT-II_RecogD-Z', 'CVLT-II_RecogD-Scaled', 'RCFT_lmm-Raw', 'RCFT_lmm-Z', 'RCFT_lmm-Scaled', 'RCFT_Del-Raw', 'RCFT_Del-Z', 'RCFT_Del-Scaled', 'RCFT_Recog-Raw', 'RCFT_Recog-Z', 'RCFT_Recog-Scaled', 'BVMT-R_Total-Raw', 'BVMT-R_Total-Z', 'BVMT-R_Total-Scaled', 'BVMT-R_Total-zGale', 'BVMT-R_Del-Raw', 'BVMT-R_Del-Z', 'BVMT-R_Del-Scaled', 'BVMT-R_Del-zGale', 'BVMT-R_RecogDis-Raw', 'BVMT-R_RecogDis-Z', 'BVMT-R_RecogDis-Scaled', 'BVMT-R_RecogDis-zGale', 'RCFT_Copy-Raw', 'RCFT_Copy-Z', 'RCFT-Copy-Scaled', 'LetFluency-Raw', 'LetFluency-Z', 'LetFluency-Scaled', 'CatFluency-Raw', 'CatFluency-Z', 'CatFluency-Scaled', 'Fluency-Raw', 'Fluency-Z', 'Fluency-Scaled', 'C/W_Inhib-Raw', 'C/W_Inhib-Z', 'C/W_Inhib-Scaled', 'Switching-Raw', 'Switching-z', 'Switching-Scaled', 'HayBTime1-Raw', 'HayBTime1-z', 'HayBTime2-Raw', 'HayBTime2-z', 'HayBCatA-Raw', 'HayBCatA-z', 'HayBCatB-Raw', 'HayBCatB-z')
colnames(neuropsych_matrixF2) <- c('ParticipantID', 'Age', 'EduLevel', 'Group', 'Ethnicity', 'Neuropsych_TestDate', 'TOPF-Raw', 'TOPF-Z', 'TOPF-Scaled', 'ClinPremorbid-z', 'ClinPremorbid-Scaled', 'ClinPremorbid-Cutoff', 'DSF-Raw', 'DSF-Z', 'DSF-Scaled', 'DSB-Raw', 'DSB-Z', 'DSB-Scaled', 'TrailsA-Raw', 'TrailsA-Z', 'TrailsA-Scaled', 'TrailsB-Raw', 'TrailsB-Z', 'TrialsB-Scaled', 'Coding-Raw', 'Coding-Z', 'Coding-Scaled', 'ColorNaming-Raw', 'ColorNaming-Z', 'ColorNaming-Scaled', 'WordReading-Raw', 'WordReading-Z', 'WordReading-Scaled', 'Inhibition-Raw', 'Inhibition-Z', 'Inhibition-Scaled', 'SYDBATNaming-Raw', 'SYDBATNaming-Z', 'SYDBATNaming-Scaled', 'SYDBATComp-Raw', 'SYDBATComp-Z', 'SYDBATComp-Scaled', 'SYDBATSemAss-Raw', 'SYDBATSemAss-Z', 'SYDBATSemAss-Scaled', 'BNTManual-Raw', 'BNTManual-Z', 'BNTManual-Scaled', 'BNTIvnik-Raw', 'BNTIvnik-Z', 'BNTIvnik-Scaled', 'Similarities-Raw', 'Similarities-Z', 'Similarities-Scaled', 'LineO-Raw', 'LineO-Z', 'LineO-Scaled', 'BD-Raw', 'BD-Z', 'BD-Scaled', 'Matrix-Raw', 'Matrix-Z', 'Matrix-Scaled', 'LM_I-Raw', 'LM_I-Z', 'LM_I-Scaled', 'LM_II-Raw', 'LM_II-Z', 'LM_II-Scaled', 'StoryAI-Raw', 'StoryAI-Z', 'StoryAI-Scaled', 'StoryAII-Raw', 'StoryAII-Z', 'StoryAII-Scaled',  'CVLT-II_Total-Raw', 'CVLT-II_Total-Z', 'CVLT-II_Total-Scaled', 'CVLT-II_Short-Raw', 'CVLT-II_Short-Z', 'CVLT-II_Short-Scaled', 'CVLT-II_Long-Raw', 'CVLT-II_Long-Z', 'CVLT-II_Long-Scaled', 'CVLT-II_RecogD-Raw', 'CVLT-II_RecogD-Z', 'CVLT-II_RecogD-Scaled', 'RCFT_lmm-Raw', 'RCFT_lmm-Z', 'RCFT_lmm-Scaled', 'RCFT_Del-Raw', 'RCFT_Del-Z', 'RCFT_Del-Scaled', 'RCFT_Recog-Raw', 'RCFT_Recog-Z', 'RCFT_Recog-Scaled', 'BVMT-R_Total-Raw', 'BVMT-R_Total-Z', 'BVMT-R_Total-Scaled', 'BVMT-R_Total-zGale', 'BVMT-R_Del-Raw', 'BVMT-R_Del-Z', 'BVMT-R_Del-Scaled', 'BVMT-R_Del-zGale', 'BVMT-R_RecogDis-Raw', 'BVMT-R_RecogDis-Z', 'BVMT-R_RecogDis-Scaled', 'BVMT-R_RecogDis-zGale', 'RCFT_Copy-Raw', 'RCFT_Copy-Z', 'RCFT-Copy-Scaled', 'LetFluency-Raw', 'LetFluency-Z', 'LetFluency-Scaled', 'CatFluency-Raw', 'CatFluency-Z', 'CatFluency-Scaled', 'Fluency-Raw', 'Fluency-Z', 'Fluency-Scaled', 'C/W_Inhib-Raw', 'C/W_Inhib-Z', 'C/W_Inhib-Scaled', 'Switching-Raw', 'Switching-z', 'Switching-Scaled', 'HayBTime1-Raw', 'HayBTime1-z', 'HayBTime2-Raw', 'HayBTime2-z', 'HayBCatA-Raw', 'HayBCatA-z', 'HayBCatB-Raw', 'HayBCatB-z')
colnames(neuropsych_matrixF3) <- c('ParticipantID', 'Age', 'EduLevel', 'Group', 'Ethnicity', 'Neuropsych_TestDate', 'TOPF-Raw', 'TOPF-Z', 'TOPF-Scaled', 'ClinPremorbid-z', 'ClinPremorbid-Scaled', 'ClinPremorbid-Cutoff', 'DSF-Raw', 'DSF-Z', 'DSF-Scaled', 'DSB-Raw', 'DSB-Z', 'DSB-Scaled', 'TrailsA-Raw', 'TrailsA-Z', 'TrailsA-Scaled', 'TrailsB-Raw', 'TrailsB-Z', 'TrialsB-Scaled', 'Coding-Raw', 'Coding-Z', 'Coding-Scaled', 'ColorNaming-Raw', 'ColorNaming-Z', 'ColorNaming-Scaled', 'WordReading-Raw', 'WordReading-Z', 'WordReading-Scaled', 'Inhibition-Raw', 'Inhibition-Z', 'Inhibition-Scaled', 'SYDBATNaming-Raw', 'SYDBATNaming-Z', 'SYDBATNaming-Scaled', 'SYDBATComp-Raw', 'SYDBATComp-Z', 'SYDBATComp-Scaled', 'SYDBATSemAss-Raw', 'SYDBATSemAss-Z', 'SYDBATSemAss-Scaled', 'BNTManual-Raw', 'BNTManual-Z', 'BNTManual-Scaled', 'BNTIvnik-Raw', 'BNTIvnik-Z', 'BNTIvnik-Scaled', 'Similarities-Raw', 'Similarities-Z', 'Similarities-Scaled', 'LineO-Raw', 'LineO-Z', 'LineO-Scaled', 'BD-Raw', 'BD-Z', 'BD-Scaled', 'Matrix-Raw', 'Matrix-Z', 'Matrix-Scaled', 'LM_I-Raw', 'LM_I-Z', 'LM_I-Scaled', 'LM_II-Raw', 'LM_II-Z', 'LM_II-Scaled', 'StoryAI-Raw', 'StoryAI-Z', 'StoryAI-Scaled', 'StoryAII-Raw', 'StoryAII-Z', 'StoryAII-Scaled',  'CVLT-II_Total-Raw', 'CVLT-II_Total-Z', 'CVLT-II_Total-Scaled', 'CVLT-II_Short-Raw', 'CVLT-II_Short-Z', 'CVLT-II_Short-Scaled', 'CVLT-II_Long-Raw', 'CVLT-II_Long-Z', 'CVLT-II_Long-Scaled', 'CVLT-II_RecogD-Raw', 'CVLT-II_RecogD-Z', 'CVLT-II_RecogD-Scaled', 'RCFT_lmm-Raw', 'RCFT_lmm-Z', 'RCFT_lmm-Scaled', 'RCFT_Del-Raw', 'RCFT_Del-Z', 'RCFT_Del-Scaled', 'RCFT_Recog-Raw', 'RCFT_Recog-Z', 'RCFT_Recog-Scaled', 'BVMT-R_Total-Raw', 'BVMT-R_Total-Z', 'BVMT-R_Total-Scaled', 'BVMT-R_Total-zGale', 'BVMT-R_Del-Raw', 'BVMT-R_Del-Z', 'BVMT-R_Del-Scaled', 'BVMT-R_Del-zGale', 'BVMT-R_RecogDis-Raw', 'BVMT-R_RecogDis-Z', 'BVMT-R_RecogDis-Scaled', 'BVMT-R_RecogDis-zGale', 'RCFT_Copy-Raw', 'RCFT_Copy-Z', 'RCFT-Copy-Scaled', 'LetFluency-Raw', 'LetFluency-Z', 'LetFluency-Scaled', 'CatFluency-Raw', 'CatFluency-Z', 'CatFluency-Scaled', 'Fluency-Raw', 'Fluency-Z', 'Fluency-Scaled', 'C/W_Inhib-Raw', 'C/W_Inhib-Z', 'C/W_Inhib-Scaled', 'Switching-Raw', 'Switching-z', 'Switching-Scaled', 'HayBTime1-Raw', 'HayBTime1-z', 'HayBTime2-Raw', 'HayBTime2-z', 'HayBCatA-Raw', 'HayBCatA-z', 'HayBCatB-Raw', 'HayBCatB-z')
colnames(neuropsych_matrixF4) <- c('ParticipantID', 'Age', 'EduLevel', 'Group', 'Ethnicity', 'Neuropsych_TestDate', 'TOPF-Raw', 'TOPF-Z', 'TOPF-Scaled', 'ClinPremorbid-z', 'ClinPremorbid-Scaled', 'ClinPremorbid-Cutoff', 'DSF-Raw', 'DSF-Z', 'DSF-Scaled', 'DSB-Raw', 'DSB-Z', 'DSB-Scaled', 'TrailsA-Raw', 'TrailsA-Z', 'TrailsA-Scaled', 'TrailsB-Raw', 'TrailsB-Z', 'TrialsB-Scaled', 'Coding-Raw', 'Coding-Z', 'Coding-Scaled', 'ColorNaming-Raw', 'ColorNaming-Z', 'ColorNaming-Scaled', 'WordReading-Raw', 'WordReading-Z', 'WordReading-Scaled', 'Inhibition-Raw', 'Inhibition-Z', 'Inhibition-Scaled', 'SYDBATNaming-Raw', 'SYDBATNaming-Z', 'SYDBATNaming-Scaled', 'SYDBATComp-Raw', 'SYDBATComp-Z', 'SYDBATComp-Scaled', 'SYDBATSemAss-Raw', 'SYDBATSemAss-Z', 'SYDBATSemAss-Scaled', 'BNTManual-Raw', 'BNTManual-Z', 'BNTManual-Scaled', 'BNTIvnik-Raw', 'BNTIvnik-Z', 'BNTIvnik-Scaled', 'Similarities-Raw', 'Similarities-Z', 'Similarities-Scaled', 'LineO-Raw', 'LineO-Z', 'LineO-Scaled', 'BD-Raw', 'BD-Z', 'BD-Scaled', 'Matrix-Raw', 'Matrix-Z', 'Matrix-Scaled', 'LM_I-Raw', 'LM_I-Z', 'LM_I-Scaled', 'LM_II-Raw', 'LM_II-Z', 'LM_II-Scaled', 'StoryAI-Raw', 'StoryAI-Z', 'StoryAI-Scaled', 'StoryAII-Raw', 'StoryAII-Z', 'StoryAII-Scaled',  'CVLT-II_Total-Raw', 'CVLT-II_Total-Z', 'CVLT-II_Total-Scaled', 'CVLT-II_Short-Raw', 'CVLT-II_Short-Z', 'CVLT-II_Short-Scaled', 'CVLT-II_Long-Raw', 'CVLT-II_Long-Z', 'CVLT-II_Long-Scaled', 'CVLT-II_RecogD-Raw', 'CVLT-II_RecogD-Z', 'CVLT-II_RecogD-Scaled', 'RCFT_lmm-Raw', 'RCFT_lmm-Z', 'RCFT_lmm-Scaled', 'RCFT_Del-Raw', 'RCFT_Del-Z', 'RCFT_Del-Scaled', 'RCFT_Recog-Raw', 'RCFT_Recog-Z', 'RCFT_Recog-Scaled', 'BVMT-R_Total-Raw', 'BVMT-R_Total-Z', 'BVMT-R_Total-Scaled', 'BVMT-R_Total-zGale', 'BVMT-R_Del-Raw', 'BVMT-R_Del-Z', 'BVMT-R_Del-Scaled', 'BVMT-R_Del-zGale', 'BVMT-R_RecogDis-Raw', 'BVMT-R_RecogDis-Z', 'BVMT-R_RecogDis-Scaled', 'BVMT-R_RecogDis-zGale', 'RCFT_Copy-Raw', 'RCFT_Copy-Z', 'RCFT-Copy-Scaled', 'LetFluency-Raw', 'LetFluency-Z', 'LetFluency-Scaled', 'CatFluency-Raw', 'CatFluency-Z', 'CatFluency-Scaled', 'Fluency-Raw', 'Fluency-Z', 'Fluency-Scaled', 'C/W_Inhib-Raw', 'C/W_Inhib-Z', 'C/W_Inhib-Scaled', 'Switching-Raw', 'Switching-z', 'Switching-Scaled', 'HayBTime1-Raw', 'HayBTime1-z', 'HayBTime2-Raw', 'HayBTime2-z', 'HayBCatA-Raw', 'HayBCatA-z', 'HayBCatB-Raw', 'HayBCatB-z')

for(i in sequence(length(files_all))){
  #Reset the task sheet triggers
  taskF0 <- 'RESET'
  taskF1 <- 'RESET'
  taskF2 <- 'RESET'
  taskF3 <- 'RESET'
  taskF4 <- 'RESET'
  
  #read in each of the participant files 
  setwd(files_all[i])
  pt_directory <- list.files() 
  xlsx_filename <- pt_directory[grepl('.xlsx', pt_directory)]
  #check the number of sheets (i.e. timepoints) in each participant excel sheet
  WorkBookLoad <- loadWorkbook(xlsx_filename)
  NumSheets <- names(getSheets(WorkBookLoad))
  
  #check if the sheets exist on the excel file. 
  for (j in 1:length(NumSheets)){
    if (NumSheets[j] == 'F0') {
      taskF0 <- 'RUN'
    }
    if (NumSheets[j] == 'F1') {
      taskF1 <- 'RUN'
    }
    if (NumSheets[j] == 'F2') {
      taskF2 <- 'RUN'
    }
    if (NumSheets[j] == 'F3') {
      taskF3 <- 'RUN' 
    }
    if (NumSheets[j] == 'F4') {
      taskF4 <- 'RUN'
    }
  }
  
  #if the sheet of the specific timepoint exists, then read in the data from the file.
  #for F0: 
  if (taskF0 == 'RUN') {
    pt_dataF0 <- read.xlsx(xlsx_filename, sheetName = 'F0')
    
    #find participant details (e.g ID, age, group) and add to matrix
    #ID
    pt_IDF0 <- pt_dataF0$Raw.1[grep("F0", pt_dataF0$Raw.1) ]
    pt_IDF0 <- paste('ADPRC_', pt_IDF0, sep='') #add in 'ADPRC' ID tag at the beginning. 
    neuropsych_matrixF0[i,1] <- pt_IDF0 
    #Age
    pt_Age_locF0 <- which(pt_dataF0 == 'Age', arr.ind=TRUE) #find participant age location on the sheet
    pt_AgeF0 <- pt_dataF0[pt_Age_locF0[1], pt_Age_locF0[2]+1]
    neuropsych_matrixF0[,2] <- pt_AgeF0
    #Education Level
    pt_EduLvl_locF0 <- which(pt_dataF0 == "Ed'n", arr.ind=TRUE) #find participant education level location on the sheet
    pt_EduLvlF0 <- pt_dataF0[pt_EduLvl_locF0[1], pt_EduLvl_locF0[2]+1]
    neuropsych_matrixF0[i,3] <- pt_EduLvlF0
    #Group status
    pt_Group_locF0 <- which(pt_dataF0 == 'Dx', arr.ind=TRUE) #find participant group status location on the sheet
    pt_GroupF0 <- pt_dataF0[pt_Group_locF0[1], pt_Group_locF0[2]+1]
    neuropsych_matrixF0[i,4] <- pt_GroupF0
    #Ethnicity
    pt_ethnicity_locF0 <- which(pt_dataF0 == 'Ethnicity', arr.ind=TRUE) #find participant ethnicity location on the sheet
    pt_ethnicityF0 <- pt_dataF0[pt_ethnicity_locF0[1], pt_ethnicity_locF0[2]+1]
    neuropsych_matrixF0[i,5] <- pt_ethnicityF0
    #Neuropsych test date
    pt_testDateF0 <- as.Date(as.integer(substr(colnames(pt_dataF0[1]), 2, 6)), origin = "1899-12-30")
    neuropsych_matrixF0[i,6] <- as.character(pt_testDateF0)
    #extract the neuropsych scores (114 total) and add into a matrix
    neuropsych_matrixF0[i,7] <- pt_dataF0[1,'Raw'] #TOPF raw
    neuropsych_matrixF0[i,8] <- pt_dataF0[1,'Z'] #TOPF Z
    neuropsych_matrixF0[i,9] <- pt_dataF0[1,'Scaled'] #TOPF Scaled
    neuropsych_matrixF0[i,10] <- pt_dataF0[2,'Raw'] #DSF raw
    neuropsych_matrixF0[i,11] <- pt_dataF0[2,'Z'] #DSF Z
    neuropsych_matrixF0[i,12] <- pt_dataF0[2,'Scaled'] #DSF Scaled
    neuropsych_matrixF0[i,13] <- pt_dataF0[3,'Raw'] #DSB raw
    neuropsych_matrixF0[i,14] <- pt_dataF0[3,'Z'] #DSB Z
    neuropsych_matrixF0[i,15] <- pt_dataF0[3,'Scaled'] #DSB Scaled
    neuropsych_matrixF0[i,16] <- pt_dataF0[4,'Raw'] #Trails A raw
    neuropsych_matrixF0[i,17] <- pt_dataF0[4,'Z'] #Trails A Z
    neuropsych_matrixF0[i,18] <- pt_dataF0[4,'Scaled'] #Trails A Scaled
    neuropsych_matrixF0[i,19] <- pt_dataF0[5,'Raw'] #Trails B raw
    neuropsych_matrixF0[i,20] <- pt_dataF0[5,'Z'] #Trails B Z
    neuropsych_matrixF0[i,21] <- pt_dataF0[5,'Scaled'] #Trails B Scaled
    neuropsych_matrixF0[i,22] <- pt_dataF0[2,'Raw.1'] #Coding raw
    neuropsych_matrixF0[i,23] <- pt_dataF0[2,'Z.1'] #Coding Z
    neuropsych_matrixF0[i,24] <- pt_dataF0[2,'Scaled.1'] #Coding Scaled
    neuropsych_matrixF0[i,25] <- pt_dataF0[3,'Raw.1'] #Colour Naming raw
    neuropsych_matrixF0[i,26] <- pt_dataF0[3,'Z.1'] #Colour Naming Z
    neuropsych_matrixF0[i,27] <- pt_dataF0[3,'Scaled.1'] #Colour Naming Scaled
    neuropsych_matrixF0[i,28] <- pt_dataF0[4,'Raw.1'] #Word Reading raw
    neuropsych_matrixF0[i,29] <- pt_dataF0[4,'Z.1'] #Word Reading Z
    neuropsych_matrixF0[i,30] <- pt_dataF0[4,'Scaled.1'] #Word Reading Scaled
    neuropsych_matrixF0[i,31] <- pt_dataF0[5,'Raw.1'] #Inhibition raw
    neuropsych_matrixF0[i,32] <- pt_dataF0[5,'Z.1'] #Inhibition Z
    neuropsych_matrixF0[i,33] <- pt_dataF0[5,'Scaled.1'] #Inhibition Scaled
    neuropsych_matrixF0[i,34] <- pt_dataF0[6,'Raw'] #SYDBAT Naming raw
    neuropsych_matrixF0[i,35] <- pt_dataF0[6,'Z'] #SYDBAT Naming Z
    neuropsych_matrixF0[i,36] <- pt_dataF0[6,'Scaled'] #SYDBAT Naming Scaled
    neuropsych_matrixF0[i,37] <- pt_dataF0[7,'Raw'] #SYDBAT Comp raw
    neuropsych_matrixF0[i,38] <- pt_dataF0[7,'Z'] #SYDBAT Comp Z
    neuropsych_matrixF0[i,39] <- pt_dataF0[7,'Scaled'] #SYDBAT Comp Scaled
    neuropsych_matrixF0[i,40] <- pt_dataF0[8,'Raw'] #SYDBAT SemAss raw
    neuropsych_matrixF0[i,41] <- pt_dataF0[8,'Z'] #SYDBAT SemAss  Z
    neuropsych_matrixF0[i,42] <- pt_dataF0[8,'Scaled'] #SYDBAT SemAss Scaled
    neuropsych_matrixF0[i,43] <- pt_dataF0[6,'Raw.1'] #BNT Man raw
    neuropsych_matrixF0[i,44] <- pt_dataF0[6,'Z.1'] #BNT Man  Z
    neuropsych_matrixF0[i,45] <- pt_dataF0[6,'Scaled.1'] #BNT Man Scaled
    neuropsych_matrixF0[i,46] <- pt_dataF0[7,'Raw.1'] #BNT Ivnik raw
    neuropsych_matrixF0[i,47] <- pt_dataF0[7,'Z.1'] #BNT Ivnik  Z
    neuropsych_matrixF0[i,48] <- pt_dataF0[7,'Scaled.1'] #BNT Ivnik Scaled
    neuropsych_matrixF0[i,49] <- pt_dataF0[8,'Raw.1'] #Similarities raw
    neuropsych_matrixF0[i,50] <- pt_dataF0[8,'Z.1'] #Similarities Z
    neuropsych_matrixF0[i,51] <- pt_dataF0[8,'Scaled.1'] #Similarities Scaled
    neuropsych_matrixF0[i,52] <- pt_dataF0[9,'Raw'] #Line O raw
    neuropsych_matrixF0[i,53] <- pt_dataF0[9,'Z'] #Line O  Z
    neuropsych_matrixF0[i,54] <- pt_dataF0[9,'Scaled'] #Line O Scaled
    neuropsych_matrixF0[i,55] <- pt_dataF0[9,'Raw.1'] #BD raw
    neuropsych_matrixF0[i,56] <- pt_dataF0[9,'Z.1'] #BD Z
    neuropsych_matrixF0[i,57] <- pt_dataF0[9,'Scaled.1'] #BD Scaled
    neuropsych_matrixF0[i,58] <- pt_dataF0[10,'Raw.1'] #Matrix raw
    neuropsych_matrixF0[i,59] <- pt_dataF0[10,'Z.1'] #Matrix Z
    neuropsych_matrixF0[i,60] <- pt_dataF0[10,'Scaled.1'] #Matrix Scaled
    neuropsych_matrixF0[i,61] <- pt_dataF0[11,'Raw'] #LM I raw
    neuropsych_matrixF0[i,62] <- pt_dataF0[11,'Z'] #LM I Z
    neuropsych_matrixF0[i,63] <- pt_dataF0[11,'Scaled'] #LM I Scaled
    neuropsych_matrixF0[i,64] <- pt_dataF0[12,'Raw'] #LM II raw
    neuropsych_matrixF0[i,65] <- pt_dataF0[12,'Z'] #LM II Z
    neuropsych_matrixF0[i,66] <- pt_dataF0[12,'Scaled'] #LM II Scaled
    neuropsych_matrixF0[i,67] <- pt_dataF0[13,'Raw'] #Story AI raw
    neuropsych_matrixF0[i,68] <- pt_dataF0[13,'Z'] #Story AI Z
    neuropsych_matrixF0[i,69] <- pt_dataF0[13,'Scaled'] #Story AI Scaled
    neuropsych_matrixF0[i,70] <- pt_dataF0[14,'Raw'] #Story AII raw
    neuropsych_matrixF0[i,71] <- pt_dataF0[14,'Z'] #Story AII Z
    neuropsych_matrixF0[i,72] <- pt_dataF0[14,'Scaled'] #Story AII Scaled
    neuropsych_matrixF0[i,73] <- pt_dataF0[11,'Raw.1'] #CVLT-II Total raw
    neuropsych_matrixF0[i,74] <- pt_dataF0[11,'Z.1'] #CVLT-II Total Z
    neuropsych_matrixF0[i,75] <- pt_dataF0[11,'Scaled.1'] #CVLT-II Total Scaled
    neuropsych_matrixF0[i,76] <- pt_dataF0[12,'Raw.1'] #CVLT-II Short raw
    neuropsych_matrixF0[i,77] <- pt_dataF0[12,'Z.1'] #CVLT-II Short Z
    neuropsych_matrixF0[i,78] <- pt_dataF0[12,'Scaled.1'] #CVLT-II Short Scaled
    neuropsych_matrixF0[i,79] <- pt_dataF0[13,'Raw.1'] #CVLT-II Long raw
    neuropsych_matrixF0[i,80] <- pt_dataF0[13,'Z.1'] #CVLT-II Long Z
    neuropsych_matrixF0[i,81] <- pt_dataF0[13,'Scaled.1'] #CVLT-II Long Scaled
    neuropsych_matrixF0[i,82] <- pt_dataF0[14,'Raw.1'] #CVLT-II Recog D raw
    neuropsych_matrixF0[i,83] <- pt_dataF0[14,'Z.1'] #CVLT-II Recog D Z
    neuropsych_matrixF0[i,84] <- pt_dataF0[14,'Scaled.1'] #CVLT-II Recog D Scaled
    neuropsych_matrixF0[i,85] <- pt_dataF0[15,'Raw'] #RCFT lmm raw
    neuropsych_matrixF0[i,86] <- pt_dataF0[15,'Z'] #RCFT lmm  Z
    neuropsych_matrixF0[i,87] <- pt_dataF0[15,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF0[i,88] <- pt_dataF0[16,'Raw'] #RCFT Del raw
    neuropsych_matrixF0[i,89] <- pt_dataF0[16,'Z'] #RCFT Del  Z
    neuropsych_matrixF0[i,90] <- pt_dataF0[16,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF0[i,91] <- pt_dataF0[17,'Raw'] #RCFT Recog raw
    neuropsych_matrixF0[i,92] <- pt_dataF0[17,'Z'] #RCFT Recog  Z
    neuropsych_matrixF0[i,93] <- pt_dataF0[17,'Scaled'] #RCFT Recog Scaled
    neuropsych_matrixF0[i,94] <- pt_dataF0[15,'Raw.1'] #BVMT Total raw
    neuropsych_matrixF0[i,95] <- pt_dataF0[15,'Z.1'] #BVMT Total Z
    neuropsych_matrixF0[i,96] <- pt_dataF0[15,'Scaled.1'] #BVMT Total Scaled
    neuropsych_matrixF0[i,97] <- pt_dataF0[15,'NA.'] #BVMT Total z Gale 
    neuropsych_matrixF0[i,98] <- pt_dataF0[16,'Raw.1'] #BVMT Del raw
    neuropsych_matrixF0[i,99] <- pt_dataF0[16,'Z.1'] #BVMT Del Z
    neuropsych_matrixF0[i,100] <- pt_dataF0[16,'Scaled.1'] #BVMT Del Scaled
    neuropsych_matrixF0[i,101] <- pt_dataF0[16,'NA.'] #BVMT Del z Gale
    neuropsych_matrixF0[i,102] <- pt_dataF0[17,'Raw.1'] #BVMT Recog Dis raw
    neuropsych_matrixF0[i,103] <- pt_dataF0[17,'Z.1'] #BVMT Recog Dis Z
    neuropsych_matrixF0[i,104] <- pt_dataF0[17,'Scaled.1'] #BVMT Recog Dis Scaled
    neuropsych_matrixF0[i,105] <- pt_dataF0[17,'NA.'] #BVMT Recog Dis z Gale
    neuropsych_matrixF0[i,106] <- pt_dataF0[18,'Raw'] #RCFT Copy raw
    neuropsych_matrixF0[i,107] <- pt_dataF0[18,'Z'] #RCFT Copy Z
    neuropsych_matrixF0[i,108] <- pt_dataF0[18,'Scaled'] #RCFT Copy Scaled
    neuropsych_matrixF0[i,109] <- pt_dataF0[19,'Raw'] #Letter Fluency raw
    neuropsych_matrixF0[i,110] <- pt_dataF0[19,'Z'] #Letter Fluency Z
    neuropsych_matrixF0[i,111] <- pt_dataF0[19,'Scaled'] #Letter Fluency Scaled
    neuropsych_matrixF0[i,112] <- pt_dataF0[20,'Raw'] #Category Fluency raw
    neuropsych_matrixF0[i,113] <- pt_dataF0[20,'Z'] #Category Fluency Z
    neuropsych_matrixF0[i,114] <- pt_dataF0[20,'Scaled'] #Category Fluency Scaled
    #neuropsych_matrixF0[i,115] <- pt_dataF0[18,'Raw'] #Fluency raw
    #neuropsych_matrixF0[i,116] <- pt_dataF0[18,'Z'] #Fluency Z
    #neuropsych_matrixF0[i,117] <- pt_dataF0[18,'Scaled'] #Fluency Scaled
    neuropsych_matrixF0[i,118] <- pt_dataF0[21,'Raw'] #C/W Inhibition raw
    neuropsych_matrixF0[i,119] <- pt_dataF0[21,'Z'] #C/W Inhibition Z
    neuropsych_matrixF0[i,120] <- pt_dataF0[21,'Scaled'] #C/W Inhibition Scaled
    neuropsych_matrixF0[i,121] <- pt_dataF0[19, 'Raw.1'] #Switching raw
    neuropsych_matrixF0[i,122] <- pt_dataF0[19, 'Z.1'] #Switching Z
    neuropsych_matrixF0[i,123] <- pt_dataF0[19, 'Scaled.1'] #Switching Scaled
    Hay_locF0 <- which(pt_dataF0 == 'Hayling B', arr.ind=TRUE) #find participant Hayling scores location on the sheet
    HayBTime1RawF0 <- pt_dataF0[Hay_locF0[1]+1, Hay_locF0[2]+1] #Hayling B Time 1 raw
    HayBTime1zF0 <- pt_dataF0[Hay_locF0[1]+1, Hay_locF0[2]+2] #Hayling B Time 1 z
    HayBTime2RawF0 <- pt_dataF0[Hay_locF0[1]+2, Hay_locF0[2]+1] #Hayling B Time 2 raw
    HayBTime2zF0 <- pt_dataF0[Hay_locF0[1]+2, Hay_locF0[2]+2] #Hayling B Time 2 z
    HayBCatARawF0 <- pt_dataF0[Hay_locF0[1]+3, Hay_locF0[2]+1] #Hayling B Cat A raw
    HayBCatAzF0 <- pt_dataF0[Hay_locF0[1]+3, Hay_locF0[2]+2] #Hayling B Cat A z
    HayBCatBRawF0 <- pt_dataF0[Hay_locF0[1]+4, Hay_locF0[2]+1] #Hayling B Cat B raw
    HayBCatBzF0 <- pt_dataF0[Hay_locF0[1]+4, Hay_locF0[2]+2] #Hayling B Cat B z
    if (class(HayBTime1RawF0) == 'NULL') {
      neuropsych_matrixF0[i,124] <- 'N/A'
    } else {
      neuropsych_matrixF0[i,124] <- HayBTime1RawF0
    }
    if (class(HayBTime1zF0) == 'NULL') {
      neuropsych_matrixF0[i,125] <- 'N/A'
    } else {
      neuropsych_matrixF0[i,125] <- HayBTime1zF0
    }
    if (class(HayBTime2RawF0) == 'NULL') {
      neuropsych_matrixF0[i,126] <- 'N/A'
    } else {
      neuropsych_matrixF0[i,126] <- HayBTime2RawF0
    }
    if (class(HayBTime2zF0) == 'NULL') {
      neuropsych_matrixF0[i,127] <- 'N/A'
    } else {
      neuropsych_matrixF0[i,127] <- HayBTime2zF0
    }
    if (class(HayBCatARawF0) == 'NULL') {
      neuropsych_matrixF0[i,128] <- 'N/A'
    } else {
      neuropsych_matrixF0[i,128] <- HayBCatARawF0
    }
    if (class(HayBCatAzF0) == 'NULL') {
      neuropsych_matrixF0[i,129] <- 'N/A'
    } else {
      neuropsych_matrixF0[i,129] <- HayBCatAzF0
    }
    if (class(HayBCatBRawF0) == 'NULL') {
      neuropsych_matrixF0[i,130] <- 'N/A'
    } else {
      neuropsych_matrixF0[i,130] <- HayBCatBRawF0
    }
    if (class(HayBCatBzF0) == 'NULL') {
      neuropsych_matrixF0[i,131] <- 'N/A'
    } else {
      neuropsych_matrixF0[i,131] <- HayBCatBzF0
    }
  }
  
  #For F1: 
  if (taskF1 == 'RUN') {
    pt_dataF1 <- read.xlsx(xlsx_filename, sheetName = 'F1')
    
    #find participant details (e.g ID, age, group) and add to matrix
    #ID
    pt_IDF1 <- pt_dataF1$Raw.1[grep("F1", pt_dataF1$Raw.1) ]
    pt_IDF1 <- paste('ADPRC_', pt_IDF1, sep='') #add in 'ADPRC' ID tag at the beginning. 
    neuropsych_matrixF1[i,1] <- pt_IDF1
    #Age
    pt_Age_locF1 <- which(pt_dataF1 == 'Age', arr.ind=TRUE) #find participant age location on the sheet
    pt_AgeF1 <- pt_dataF1[pt_Age_locF1[1], pt_Age_locF1[2]+1]
    neuropsych_matrixF1[i,2] <- pt_AgeF1
    #Education Level
    pt_EduLvl_locF1 <- which(pt_dataF1 == "Ed'n", arr.ind=TRUE) #find participant education level location on the sheet
    pt_EduLvlF1 <- pt_dataF1[pt_EduLvl_locF1[1], pt_EduLvl_locF1[2]+1]
    neuropsych_matrixF1[i,3] <- pt_EduLvlF1
    #Group status
    pt_Group_locF1 <- which(pt_dataF1 == 'Dx', arr.ind=TRUE) #find participant group status location on the sheet
    pt_GroupF1 <- pt_dataF1[pt_Group_locF1[1], pt_Group_locF1[2]+1]
    neuropsych_matrixF1[i,4] <- pt_GroupF1
    #Ethnicity
    pt_ethnicity_locF1 <- which(pt_dataF1 == 'Ethnicity', arr.ind=TRUE) #find participant ethnicity location on the sheet
    pt_ethnicityF1 <- pt_dataF0[pt_ethnicity_locF1[1], pt_ethnicity_locF1[2]+1]
    neuropsych_matrixF1[i,5] <- pt_ethnicityF1
    #Neuropsych test date
    pt_testDateF1 <- as.Date(as.integer(substr(colnames(pt_dataF1[1]), 2, 6)), origin = "1899-12-30")
    neuropsych_matrixF1[i,6] <- as.character(pt_testDateF1)
    #extract the neuropsych scores (114 total) and add into a matrix
    neuropsych_matrixF1[i,7] <- pt_dataF1[1,'Raw'] #TOPF raw
    neuropsych_matrixF1[i,8] <- pt_dataF1[1,'Z'] #TOPF Z
    neuropsych_matrixF1[i,9] <- pt_dataF1[1,'Scaled'] #TOPF Scaled
    neuropsych_matrixF1[i,10] <- pt_dataF1[2,'Raw'] #DSF raw
    neuropsych_matrixF1[i,11] <- pt_dataF1[2,'Z'] #DSF Z
    neuropsych_matrixF1[i,12] <- pt_dataF1[2,'Scaled'] #DSF Scaled
    neuropsych_matrixF1[i,13] <- pt_dataF1[3,'Raw'] #DSB raw
    neuropsych_matrixF1[i,14] <- pt_dataF1[3,'Z'] #DSB Z
    neuropsych_matrixF1[i,15] <- pt_dataF1[3,'Scaled'] #DSB Scaled
    neuropsych_matrixF1[i,16] <- pt_dataF1[4,'Raw'] #Trails A raw
    neuropsych_matrixF1[i,17] <- pt_dataF1[4,'Z'] #Trails A Z
    neuropsych_matrixF1[i,18] <- pt_dataF1[4,'Scaled'] #Trails A Scaled
    neuropsych_matrixF1[i,19] <- pt_dataF1[5,'Raw'] #Trails B raw
    neuropsych_matrixF1[i,20] <- pt_dataF1[5,'Z'] #Trails B Z
    neuropsych_matrixF1[i,21] <- pt_dataF1[5,'Scaled'] #Trails B Scaled
    neuropsych_matrixF1[i,22] <- pt_dataF1[2,'Raw.1'] #Coding raw
    neuropsych_matrixF1[i,23] <- pt_dataF1[2,'Z.1'] #Coding Z
    neuropsych_matrixF1[i,24] <- pt_dataF1[2,'Scaled.1'] #Coding Scaled
    neuropsych_matrixF1[i,25] <- pt_dataF1[3,'Raw.1'] #Colour Naming raw
    neuropsych_matrixF1[i,26] <- pt_dataF1[3,'Z.1'] #Colour Naming Z
    neuropsych_matrixF1[i,27] <- pt_dataF1[3,'Scaled.1'] #Colour Naming Scaled
    neuropsych_matrixF1[i,28] <- pt_dataF1[4,'Raw.1'] #Word Reading raw
    neuropsych_matrixF1[i,29] <- pt_dataF1[4,'Z.1'] #Word Reading Z
    neuropsych_matrixF1[i,30] <- pt_dataF1[4,'Scaled.1'] #Word Reading Scaled
    neuropsych_matrixF1[i,31] <- pt_dataF1[5,'Raw.1'] #Inhibition raw
    neuropsych_matrixF1[i,32] <- pt_dataF1[5,'Z.1'] #Inhibition Z
    neuropsych_matrixF1[i,33] <- pt_dataF1[5,'Scaled.1'] #Inhibition Scaled
    neuropsych_matrixF1[i,34] <- pt_dataF1[6,'Raw'] #SYDBAT Naming raw
    neuropsych_matrixF1[i,35] <- pt_dataF1[6,'Z'] #SYDBAT Naming Z
    neuropsych_matrixF1[i,36] <- pt_dataF1[6,'Scaled'] #SYDBAT Naming Scaled
    neuropsych_matrixF1[i,37] <- pt_dataF1[7,'Raw'] #SYDBAT Comp raw
    neuropsych_matrixF1[i,38] <- pt_dataF1[7,'Z'] #SYDBAT Comp Z
    neuropsych_matrixF1[i,39] <- pt_dataF1[7,'Scaled'] #SYDBAT Comp Scaled
    neuropsych_matrixF1[i,40] <- pt_dataF1[8,'Raw'] #SYDBAT SemAss raw
    neuropsych_matrixF1[i,41] <- pt_dataF1[8,'Z'] #SYDBAT SemAss  Z
    neuropsych_matrixF1[i,42] <- pt_dataF1[8,'Scaled'] #SYDBAT SemAss Scaled
    neuropsych_matrixF1[i,43] <- pt_dataF1[6,'Raw.1'] #BNT Man raw
    neuropsych_matrixF1[i,44] <- pt_dataF1[6,'Z.1'] #BNT Man  Z
    neuropsych_matrixF1[i,45] <- pt_dataF1[6,'Scaled.1'] #BNT Man Scaled
    neuropsych_matrixF1[i,46] <- pt_dataF1[7,'Raw.1'] #BNT Ivnik raw
    neuropsych_matrixF1[i,47] <- pt_dataF1[7,'Z.1'] #BNT Ivnik  Z
    neuropsych_matrixF1[i,48] <- pt_dataF1[7,'Scaled.1'] #BNT Ivnik Scaled
    neuropsych_matrixF1[i,49] <- pt_dataF1[8,'Raw.1'] #Similarities raw
    neuropsych_matrixF1[i,50] <- pt_dataF1[8,'Z.1'] #Similarities Z
    neuropsych_matrixF1[i,51] <- pt_dataF1[8,'Scaled.1'] #Similarities Scaled
    neuropsych_matrixF1[i,52] <- pt_dataF1[9,'Raw'] #Line O raw
    neuropsych_matrixF1[i,53] <- pt_dataF1[9,'Z'] #Line O  Z
    neuropsych_matrixF1[i,54] <- pt_dataF1[9,'Scaled'] #Line O Scaled
    neuropsych_matrixF1[i,55] <- pt_dataF1[9,'Raw.1'] #BD raw
    neuropsych_matrixF1[i,56] <- pt_dataF1[9,'Z.1'] #BD Z
    neuropsych_matrixF1[i,57] <- pt_dataF1[9,'Scaled.1'] #BD Scaled
    neuropsych_matrixF1[i,58] <- pt_dataF1[10,'Raw.1'] #Matrix raw
    neuropsych_matrixF1[i,59] <- pt_dataF1[10,'Z.1'] #Matrix Z
    neuropsych_matrixF1[i,60] <- pt_dataF1[10,'Scaled.1'] #Matrix Scaled
    neuropsych_matrixF1[i,61] <- pt_dataF1[11,'Raw'] #LM I raw
    neuropsych_matrixF1[i,62] <- pt_dataF1[11,'Z'] #LM I Z
    neuropsych_matrixF1[i,63] <- pt_dataF1[11,'Scaled'] #LM I Scaled
    neuropsych_matrixF1[i,64] <- pt_dataF1[12,'Raw'] #LM II raw
    neuropsych_matrixF1[i,65] <- pt_dataF1[12,'Z'] #LM II Z
    neuropsych_matrixF1[i,66] <- pt_dataF1[12,'Scaled'] #LM II Scaled
    neuropsych_matrixF1[i,67] <- pt_dataF1[13,'Raw'] #Story AI raw
    neuropsych_matrixF1[i,68] <- pt_dataF1[13,'Z'] #Story AI Z
    neuropsych_matrixF1[i,69] <- pt_dataF1[13,'Scaled'] #Story AI Scaled
    neuropsych_matrixF1[i,70] <- pt_dataF1[14,'Raw'] #Story AII raw
    neuropsych_matrixF1[i,71] <- pt_dataF1[14,'Z'] #Story AII Z
    neuropsych_matrixF1[i,72] <- pt_dataF1[14,'Scaled'] #Story AII Scaled
    neuropsych_matrixF1[i,73] <- pt_dataF1[11,'Raw.1'] #CVLT-II Total raw
    neuropsych_matrixF1[i,74] <- pt_dataF1[11,'Z.1'] #CVLT-II Total Z
    neuropsych_matrixF1[i,75] <- pt_dataF1[11,'Scaled.1'] #CVLT-II Total Scaled
    neuropsych_matrixF1[i,76] <- pt_dataF1[12,'Raw.1'] #CVLT-II Short raw
    neuropsych_matrixF1[i,77] <- pt_dataF1[12,'Z.1'] #CVLT-II Short Z
    neuropsych_matrixF1[i,78] <- pt_dataF1[12,'Scaled.1'] #CVLT-II Short Scaled
    neuropsych_matrixF1[i,79] <- pt_dataF1[13,'Raw.1'] #CVLT-II Long raw
    neuropsych_matrixF1[i,80] <- pt_dataF1[13,'Z.1'] #CVLT-II Long Z
    neuropsych_matrixF1[i,81] <- pt_dataF1[13,'Scaled.1'] #CVLT-II Long Scaled
    neuropsych_matrixF1[i,82] <- pt_dataF1[14,'Raw.1'] #CVLT-II Recog D raw
    neuropsych_matrixF1[i,83] <- pt_dataF1[14,'Z.1'] #CVLT-II Recog D Z
    neuropsych_matrixF1[i,84] <- pt_dataF1[14,'Scaled.1'] #CVLT-II Recog D Scaled
    neuropsych_matrixF1[i,85] <- pt_dataF1[15,'Raw'] #RCFT lmm raw
    neuropsych_matrixF1[i,86] <- pt_dataF1[15,'Z'] #RCFT lmm  Z
    neuropsych_matrixF1[i,87] <- pt_dataF1[15,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF1[i,88] <- pt_dataF1[16,'Raw'] #RCFT Del raw
    neuropsych_matrixF1[i,89] <- pt_dataF1[16,'Z'] #RCFT Del  Z
    neuropsych_matrixF1[i,90] <- pt_dataF1[16,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF1[i,91] <- pt_dataF1[17,'Raw'] #RCFT Recog raw
    neuropsych_matrixF1[i,92] <- pt_dataF1[17,'Z'] #RCFT Recog  Z
    neuropsych_matrixF1[i,93] <- pt_dataF1[17,'Scaled'] #RCFT Recog Scaled
    neuropsych_matrixF1[i,94] <- pt_dataF1[15,'Raw.1'] #BVMT Total raw
    neuropsych_matrixF1[i,95] <- pt_dataF1[15,'Z.1'] #BVMT Total Z
    neuropsych_matrixF1[i,96] <- pt_dataF1[15,'Scaled.1'] #BVMT Total Scaled
    neuropsych_matrixF1[i,97] <- pt_dataF1[15,'NA.'] #BVMT Total z Gale 
    neuropsych_matrixF1[i,98] <- pt_dataF1[16,'Raw.1'] #BVMT Del raw
    neuropsych_matrixF1[i,99] <- pt_dataF1[16,'Z.1'] #BVMT Del Z
    neuropsych_matrixF1[i,100] <- pt_dataF1[16,'Scaled.1'] #BVMT Del Scaled
    neuropsych_matrixF1[i,101] <- pt_dataF1[16,'NA.'] #BVMT Del z Gale
    neuropsych_matrixF1[i,102] <- pt_dataF1[17,'Raw.1'] #BVMT Recog Dis raw
    neuropsych_matrixF1[i,103] <- pt_dataF1[17,'Z.1'] #BVMT Recog Dis Z
    neuropsych_matrixF1[i,104] <- pt_dataF1[17,'Scaled.1'] #BVMT Recog Dis Scaled
    neuropsych_matrixF1[i,105] <- pt_dataF1[17,'NA.'] #BVMT Recog Dis z Gale
    neuropsych_matrixF1[i,106] <- pt_dataF1[18,'Raw'] #RCFT Copy raw
    neuropsych_matrixF1[i,107] <- pt_dataF1[18,'Z'] #RCFT Copy Z
    neuropsych_matrixF1[i,108] <- pt_dataF1[18,'Scaled'] #RCFT Copy Scaled
    neuropsych_matrixF1[i,109] <- pt_dataF1[19,'Raw'] #Letter Fluency raw
    neuropsych_matrixF1[i,110] <- pt_dataF1[19,'Z'] #Letter Fluency Z
    neuropsych_matrixF1[i,111] <- pt_dataF1[19,'Scaled'] #Letter Fluency Scaled
    neuropsych_matrixF1[i,112] <- pt_dataF1[20,'Raw'] #Category Fluency raw
    neuropsych_matrixF1[i,113] <- pt_dataF1[20,'Z'] #Category Fluency Z
    neuropsych_matrixF1[i,114] <- pt_dataF1[20,'Scaled'] #Category Fluency Scaled
    #neuropsych_matrixF1[i,115] <- pt_dataF1[18,'Raw'] #Fluency raw
    #neuropsych_matrixF1[i,116] <- pt_dataF1[18,'Z'] #Fluency Z
    #neuropsych_matrixF1[i,117] <- pt_dataF1[18,'Scaled'] #Fluency Scaled
    neuropsych_matrixF1[i,118] <- pt_dataF1[21,'Raw'] #C/W Inhibition raw
    neuropsych_matrixF1[i,119] <- pt_dataF1[21,'Z'] #C/W Inhibition Z
    neuropsych_matrixF1[i,120] <- pt_dataF1[21,'Scaled'] #C/W Inhibition Scaled
    neuropsych_matrixF1[i,121] <- pt_dataF1[19, 'Raw.1'] #Switching raw
    neuropsych_matrixF1[i,122] <- pt_dataF1[19, 'Z.1'] #Switching Z
    neuropsych_matrixF1[i,123] <- pt_dataF1[19, 'Scaled.1'] #Switching Scaled
    Hay_locF1 <- which(pt_dataF1 == 'Hayling B', arr.ind=TRUE) #find participant Hayling scores location on the sheet
    HayBTime1RawF1 <- pt_dataF1[Hay_locF1[1]+1, Hay_locF1[2]+1] #Hayling B Time 1 raw
    HayBTime1zF1 <- pt_dataF1[Hay_locF1[1]+1, Hay_locF1[2]+2] #Hayling B Time 1 z
    HayBTime2RawF1 <- pt_dataF1[Hay_locF1[1]+2, Hay_locF1[2]+1] #Hayling B Time 2 raw
    HayBTime2zF1 <- pt_dataF1[Hay_locF1[1]+2, Hay_locF1[2]+2] #Hayling B Time 2 z
    HayBCatARawF1 <- pt_dataF1[Hay_locF1[1]+3, Hay_locF1[2]+1] #Hayling B Cat A raw
    HayBCatAzF1 <- pt_dataF1[Hay_locF1[1]+3, Hay_locF1[2]+2] #Hayling B Cat A z
    HayBCatBRawF1 <- pt_dataF1[Hay_locF1[1]+4, Hay_locF1[2]+1] #Hayling B Cat B raw
    HayBCatBzF1 <- pt_dataF1[Hay_locF1[1]+4, Hay_locF1[2]+2] #Hayling B Cat B z
    if (class(HayBTime1RawF1) == 'NULL') {
      neuropsych_matrixF1[i,124] <- 'N/A'
    } else {
      neuropsych_matrixF1[i,124] <- HayBTime1RawF1
    }
    if (class(HayBTime1zF1) == 'NULL') {
      neuropsych_matrixF1[i,125] <- 'N/A'
    } else {
      neuropsych_matrixF1[i,125] <- HayBTime1zF1
    }
    if (class(HayBTime2RawF1) == 'NULL') {
      neuropsych_matrixF1[i,126] <- 'N/A'
    } else {
      neuropsych_matrixF1[i,126] <- HayBTime2RawF1
    }
    if (class(HayBTime2zF1) == 'NULL') {
      neuropsych_matrixF1[i,127] <- 'N/A'
    } else {
      neuropsych_matrixF1[i,127] <- HayBTime2zF1
    }
    if (class(HayBCatARawF1) == 'NULL') {
      neuropsych_matrixF1[i,128] <- 'N/A'
    } else {
      neuropsych_matrixF1[i,128] <- HayBCatARawF1
    }
    if (class(HayBCatAzF1) == 'NULL') {
      neuropsych_matrixF1[i,129] <- 'N/A'
    } else {
      neuropsych_matrixF1[i,129] <- HayBCatAzF1
    }
    if (class(HayBCatBRawF1) == 'NULL') {
      neuropsych_matrixF1[i,130] <- 'N/A'
    } else {
      neuropsych_matrixF1[i,130] <- HayBCatBRawF1
    }
    if (class(HayBCatBzF1) == 'NULL') {
      neuropsych_matrixF1[i,131] <- 'N/A'
    } else {
      neuropsych_matrixF1[i,131] <- HayBCatBzF1
    }
  }
  
  #For F2: 
  if (taskF2 == 'RUN') {
    pt_dataF2 <- read.xlsx(xlsx_filename, sheetName = 'F2')
    
    #find participant details (e.g ID, age, group) and add to matrix
    #ID
    pt_IDF2 <- pt_dataF2$Raw.1[grep("F2", pt_dataF2$Raw.1) ]
    pt_IDF2 <- paste('ADPRC_', pt_IDF2, sep='') #add in 'ADPRC' ID tag at the beginning. 
    neuropsych_matrixF2[i,1] <- pt_IDF2
    #Age
    pt_Age_locF2 <- which(pt_dataF2 == 'Age', arr.ind=TRUE) #find participant age location on the sheet
    pt_AgeF2 <- pt_dataF2[pt_Age_locF2[1], pt_Age_locF2[2]+1]
    neuropsych_matrixF2[i,2] <- pt_AgeF2
    #Education Level
    pt_EduLvl_locF2 <- which(pt_dataF2 == "Ed'n", arr.ind=TRUE) #find participant education level location on the sheet
    pt_EduLvlF2 <- pt_dataF2[pt_EduLvl_locF2[1], pt_EduLvl_locF2[2]+1]
    neuropsych_matrixF2[i,3] <- pt_EduLvlF2
    #Group status
    pt_Group_locF2 <- which(pt_dataF2 == 'Dx', arr.ind=TRUE) #find participant group status location on the sheet
    pt_GroupF2 <- pt_dataF2[pt_Group_locF2[1], pt_Group_locF2[2]+1]
    neuropsych_matrixF2[i,4] <- pt_GroupF2
    #Ethnicity
    pt_ethnicity_locF2 <- which(pt_dataF2 == 'Ethnicity', arr.ind=TRUE) #find participant ethnicity location on the sheet
    pt_ethnicityF2 <- pt_dataF2[pt_ethnicity_locF2[1], pt_ethnicity_locF2[2]+1]
    neuropsych_matrixF2[i,5] <- pt_ethnicityF2
    #Neuropsych test date
    pt_testDateF2 <- as.Date(as.integer(substr(colnames(pt_dataF2[1]), 2, 6)), origin = "1899-12-30")
    neuropsych_matrixF2[i,6] <- as.character(pt_testDateF2)
    #extract the neuropsych scores (114 total) and add into a matrix
    neuropsych_matrixF2[i,7] <- pt_dataF2[1,'Raw'] #TOPF raw
    neuropsych_matrixF2[i,8] <- pt_dataF2[1,'Z'] #TOPF Z
    neuropsych_matrixF2[i,9] <- pt_dataF2[1,'Scaled'] #TOPF Scaled
    neuropsych_matrixF2[i,10] <- pt_dataF2[2,'Raw'] #DSF raw
    neuropsych_matrixF2[i,11] <- pt_dataF2[2,'Z'] #DSF Z
    neuropsych_matrixF2[i,12] <- pt_dataF2[2,'Scaled'] #DSF Scaled
    neuropsych_matrixF2[i,13] <- pt_dataF2[3,'Raw'] #DSB raw
    neuropsych_matrixF2[i,14] <- pt_dataF2[3,'Z'] #DSB Z
    neuropsych_matrixF2[i,15] <- pt_dataF2[3,'Scaled'] #DSB Scaled
    neuropsych_matrixF2[i,16] <- pt_dataF2[4,'Raw'] #Trails A raw
    neuropsych_matrixF2[i,17] <- pt_dataF2[4,'Z'] #Trails A Z
    neuropsych_matrixF2[i,18] <- pt_dataF2[4,'Scaled'] #Trails A Scaled
    neuropsych_matrixF2[i,19] <- pt_dataF2[5,'Raw'] #Trails B raw
    neuropsych_matrixF2[i,20] <- pt_dataF2[5,'Z'] #Trails B Z
    neuropsych_matrixF2[i,21] <- pt_dataF2[5,'Scaled'] #Trails B Scaled
    neuropsych_matrixF2[i,22] <- pt_dataF2[2,'Raw.1'] #Coding raw
    neuropsych_matrixF2[i,23] <- pt_dataF2[2,'Z.1'] #Coding Z
    neuropsych_matrixF2[i,24] <- pt_dataF2[2,'Scaled.1'] #Coding Scaled
    neuropsych_matrixF2[i,25] <- pt_dataF2[3,'Raw.1'] #Colour Naming raw
    neuropsych_matrixF2[i,26] <- pt_dataF2[3,'Z.1'] #Colour Naming Z
    neuropsych_matrixF2[i,27] <- pt_dataF2[3,'Scaled.1'] #Colour Naming Scaled
    neuropsych_matrixF2[i,28] <- pt_dataF2[4,'Raw.1'] #Word Reading raw
    neuropsych_matrixF2[i,29] <- pt_dataF2[4,'Z.1'] #Word Reading Z
    neuropsych_matrixF2[i,30] <- pt_dataF2[4,'Scaled.1'] #Word Reading Scaled
    neuropsych_matrixF2[i,31] <- pt_dataF2[5,'Raw.1'] #Inhibition raw
    neuropsych_matrixF2[i,32] <- pt_dataF2[5,'Z.1'] #Inhibition Z
    neuropsych_matrixF2[i,33] <- pt_dataF2[5,'Scaled.1'] #Inhibition Scaled
    neuropsych_matrixF2[i,34] <- pt_dataF2[6,'Raw'] #SYDBAT Naming raw
    neuropsych_matrixF2[i,35] <- pt_dataF2[6,'Z'] #SYDBAT Naming Z
    neuropsych_matrixF2[i,36] <- pt_dataF2[6,'Scaled'] #SYDBAT Naming Scaled
    neuropsych_matrixF2[i,37] <- pt_dataF2[7,'Raw'] #SYDBAT Comp raw
    neuropsych_matrixF2[i,38] <- pt_dataF2[7,'Z'] #SYDBAT Comp Z
    neuropsych_matrixF2[i,39] <- pt_dataF2[7,'Scaled'] #SYDBAT Comp Scaled
    neuropsych_matrixF2[i,40] <- pt_dataF2[8,'Raw'] #SYDBAT SemAss raw
    neuropsych_matrixF2[i,41] <- pt_dataF2[8,'Z'] #SYDBAT SemAss  Z
    neuropsych_matrixF2[i,42] <- pt_dataF2[8,'Scaled'] #SYDBAT SemAss Scaled
    neuropsych_matrixF2[i,43] <- pt_dataF2[6,'Raw.1'] #BNT Man raw
    neuropsych_matrixF2[i,44] <- pt_dataF2[6,'Z.1'] #BNT Man  Z
    neuropsych_matrixF2[i,45] <- pt_dataF2[6,'Scaled.1'] #BNT Man Scaled
    neuropsych_matrixF2[i,46] <- pt_dataF2[7,'Raw.1'] #BNT Ivnik raw
    neuropsych_matrixF2[i,47] <- pt_dataF2[7,'Z.1'] #BNT Ivnik  Z
    neuropsych_matrixF2[i,48] <- pt_dataF2[7,'Scaled.1'] #BNT Ivnik Scaled
    neuropsych_matrixF2[i,49] <- pt_dataF2[8,'Raw.1'] #Similarities raw
    neuropsych_matrixF2[i,50] <- pt_dataF2[8,'Z.1'] #Similarities Z
    neuropsych_matrixF2[i,51] <- pt_dataF2[8,'Scaled.1'] #Similarities Scaled
    neuropsych_matrixF2[i,52] <- pt_dataF2[9,'Raw'] #Line O raw
    neuropsych_matrixF2[i,53] <- pt_dataF2[9,'Z'] #Line O  Z
    neuropsych_matrixF2[i,54] <- pt_dataF2[9,'Scaled'] #Line O Scaled
    neuropsych_matrixF2[i,55] <- pt_dataF2[9,'Raw.1'] #BD raw
    neuropsych_matrixF2[i,56] <- pt_dataF2[9,'Z.1'] #BD Z
    neuropsych_matrixF2[i,57] <- pt_dataF2[9,'Scaled.1'] #BD Scaled
    neuropsych_matrixF2[i,58] <- pt_dataF2[10,'Raw.1'] #Matrix raw
    neuropsych_matrixF2[i,59] <- pt_dataF2[10,'Z.1'] #Matrix Z
    neuropsych_matrixF2[i,60] <- pt_dataF2[10,'Scaled.1'] #Matrix Scaled
    neuropsych_matrixF2[i,61] <- pt_dataF2[11,'Raw'] #LM I raw
    neuropsych_matrixF2[i,62] <- pt_dataF2[11,'Z'] #LM I Z
    neuropsych_matrixF2[i,63] <- pt_dataF2[11,'Scaled'] #LM I Scaled
    neuropsych_matrixF2[i,64] <- pt_dataF2[12,'Raw'] #LM II raw
    neuropsych_matrixF2[i,65] <- pt_dataF2[12,'Z'] #LM II Z
    neuropsych_matrixF2[i,66] <- pt_dataF2[12,'Scaled'] #LM II Scaled
    neuropsych_matrixF2[i,67] <- pt_dataF2[13,'Raw'] #Story AI raw
    neuropsych_matrixF2[i,68] <- pt_dataF2[13,'Z'] #Story AI Z
    neuropsych_matrixF2[i,69] <- pt_dataF2[13,'Scaled'] #Story AI Scaled
    neuropsych_matrixF2[i,70] <- pt_dataF2[14,'Raw'] #Story AII raw
    neuropsych_matrixF2[i,71] <- pt_dataF2[14,'Z'] #Story AII Z
    neuropsych_matrixF2[i,72] <- pt_dataF2[14,'Scaled'] #Story AII Scaled
    neuropsych_matrixF2[i,73] <- pt_dataF2[11,'Raw.1'] #CVLT-II Total raw
    neuropsych_matrixF2[i,74] <- pt_dataF2[11,'Z.1'] #CVLT-II Total Z
    neuropsych_matrixF2[i,75] <- pt_dataF2[11,'Scaled.1'] #CVLT-II Total Scaled
    neuropsych_matrixF2[i,76] <- pt_dataF2[12,'Raw.1'] #CVLT-II Short raw
    neuropsych_matrixF2[i,77] <- pt_dataF2[12,'Z.1'] #CVLT-II Short Z
    neuropsych_matrixF2[i,78] <- pt_dataF2[12,'Scaled.1'] #CVLT-II Short Scaled
    neuropsych_matrixF2[i,79] <- pt_dataF2[13,'Raw.1'] #CVLT-II Long raw
    neuropsych_matrixF2[i,80] <- pt_dataF2[13,'Z.1'] #CVLT-II Long Z
    neuropsych_matrixF2[i,81] <- pt_dataF2[13,'Scaled.1'] #CVLT-II Long Scaled
    neuropsych_matrixF2[i,82] <- pt_dataF2[14,'Raw.1'] #CVLT-II Recog D raw
    neuropsych_matrixF2[i,83] <- pt_dataF2[14,'Z.1'] #CVLT-II Recog D Z
    neuropsych_matrixF2[i,84] <- pt_dataF2[14,'Scaled.1'] #CVLT-II Recog D Scaled
    neuropsych_matrixF2[i,85] <- pt_dataF2[15,'Raw'] #RCFT lmm raw
    neuropsych_matrixF2[i,86] <- pt_dataF2[15,'Z'] #RCFT lmm  Z
    neuropsych_matrixF2[i,87] <- pt_dataF2[15,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF2[i,88] <- pt_dataF2[16,'Raw'] #RCFT Del raw
    neuropsych_matrixF2[i,89] <- pt_dataF2[16,'Z'] #RCFT Del  Z
    neuropsych_matrixF2[i,90] <- pt_dataF2[16,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF2[i,91] <- pt_dataF2[17,'Raw'] #RCFT Recog raw
    neuropsych_matrixF2[i,92] <- pt_dataF2[17,'Z'] #RCFT Recog  Z
    neuropsych_matrixF2[i,93] <- pt_dataF2[17,'Scaled'] #RCFT Recog Scaled
    neuropsych_matrixF2[i,94] <- pt_dataF2[15,'Raw.1'] #BVMT Total raw
    neuropsych_matrixF2[i,95] <- pt_dataF2[15,'Z.1'] #BVMT Total Z
    neuropsych_matrixF2[i,96] <- pt_dataF2[15,'Scaled.1'] #BVMT Total Scaled
    neuropsych_matrixF2[i,97] <- pt_dataF2[15,'NA.'] #BVMT Total z Gale 
    neuropsych_matrixF2[i,98] <- pt_dataF2[16,'Raw.1'] #BVMT Del raw
    neuropsych_matrixF2[i,99] <- pt_dataF2[16,'Z.1'] #BVMT Del Z
    neuropsych_matrixF2[i,100] <- pt_dataF2[16,'Scaled.1'] #BVMT Del Scaled
    neuropsych_matrixF2[i,101] <- pt_dataF2[16,'NA.'] #BVMT Del z Gale
    neuropsych_matrixF2[i,102] <- pt_dataF2[17,'Raw.1'] #BVMT Recog Dis raw
    neuropsych_matrixF2[i,103] <- pt_dataF2[17,'Z.1'] #BVMT Recog Dis Z
    neuropsych_matrixF2[i,104] <- pt_dataF2[17,'Scaled.1'] #BVMT Recog Dis Scaled
    neuropsych_matrixF2[i,105] <- pt_dataF2[17,'NA.'] #BVMT Recog Dis z Gale
    neuropsych_matrixF2[i,106] <- pt_dataF2[18,'Raw'] #RCFT Copy raw
    neuropsych_matrixF2[i,107] <- pt_dataF2[18,'Z'] #RCFT Copy Z
    neuropsych_matrixF2[i,108] <- pt_dataF2[18,'Scaled'] #RCFT Copy Scaled
    neuropsych_matrixF2[i,109] <- pt_dataF2[19,'Raw'] #Letter Fluency raw
    neuropsych_matrixF2[i,110] <- pt_dataF2[19,'Z'] #Letter Fluency Z
    neuropsych_matrixF2[i,111] <- pt_dataF2[19,'Scaled'] #Letter Fluency Scaled
    neuropsych_matrixF2[i,112] <- pt_dataF2[20,'Raw'] #Category Fluency raw
    neuropsych_matrixF2[i,113] <- pt_dataF2[20,'Z'] #Category Fluency Z
    neuropsych_matrixF2[i,114] <- pt_dataF2[20,'Scaled'] #Category Fluency Scaled
    #neuropsych_matrixF2[i,115] <- pt_dataF2[18,'Raw'] #Fluency raw
    #neuropsych_matrixF2[i,116] <- pt_dataF2[18,'Z'] #Fluency Z
    #neuropsych_matrixF2[i,117] <- pt_dataF2[18,'Scaled'] #Fluency Scaled
    neuropsych_matrixF2[i,118] <- pt_dataF2[21,'Raw'] #C/W Inhibition raw
    neuropsych_matrixF2[i,119] <- pt_dataF2[21,'Z'] #C/W Inhibition Z
    neuropsych_matrixF2[i,120] <- pt_dataF2[21,'Scaled'] #C/W Inhibition Scaled
    neuropsych_matrixF2[i,121] <- pt_dataF2[19, 'Raw.1'] #Switching raw
    neuropsych_matrixF2[i,122] <- pt_dataF2[19, 'Z.1'] #Switching Z
    neuropsych_matrixF2[i,123] <- pt_dataF2[19, 'Scaled.1'] #Switching Scaled
    Hay_locF2 <- which(pt_dataF2 == 'Hayling B', arr.ind=TRUE) #find participant Hayling scores location on the sheet
    HayBTime1RawF2 <- pt_dataF2[Hay_locF2[1]+1, Hay_locF2[2]+1] #Hayling B Time 1 raw
    HayBTime1zF2 <- pt_dataF2[Hay_locF2[1]+1, Hay_locF2[2]+2] #Hayling B Time 1 z
    HayBTime2RawF2 <- pt_dataF2[Hay_locF2[1]+2, Hay_locF2[2]+1] #Hayling B Time 2 raw
    HayBTime2zF2 <- pt_dataF2[Hay_locF2[1]+2, Hay_locF2[2]+2] #Hayling B Time 2 z
    HayBCatARawF2 <- pt_dataF2[Hay_locF2[1]+3, Hay_locF2[2]+1] #Hayling B Cat A raw
    HayBCatAzF2 <- pt_dataF2[Hay_locF2[1]+3, Hay_locF2[2]+2] #Hayling B Cat A z
    HayBCatBRawF2 <- pt_dataF2[Hay_locF2[1]+4, Hay_locF2[2]+1] #Hayling B Cat B raw
    HayBCatBzF2 <- pt_dataF2[Hay_locF2[1]+4, Hay_locF2[2]+2] #Hayling B Cat B z
    if (class(HayBTime1RawF2) == 'NULL') {
      neuropsych_matrixF2[i,124] <- 'N/A'
    } else {
      neuropsych_matrixF2[i,124] <- HayBTime1RawF2
    }
    if (class(HayBTime1zF2) == 'NULL') {
      neuropsych_matrixF2[i,125] <- 'N/A'
    } else {
      neuropsych_matrixF2[i,125] <- HayBTime1zF2
    }
    if (class(HayBTime2RawF2) == 'NULL') {
      neuropsych_matrixF2[i,126] <- 'N/A'
    } else {
      neuropsych_matrixF2[i,126] <- HayBTime2RawF2
    }
    if (class(HayBTime2zF2) == 'NULL') {
      neuropsych_matrixF2[i,127] <- 'N/A'
    } else {
      neuropsych_matrixF2[i,127] <- HayBTime2zF2
    }
    if (class(HayBCatARawF2) == 'NULL') {
      neuropsych_matrixF2[i,128] <- 'N/A'
    } else {
      neuropsych_matrixF2[i,128] <- HayBCatARawF2
    }
    if (class(HayBCatAzF2) == 'NULL') {
      neuropsych_matrixF2[i,129] <- 'N/A'
    } else {
      neuropsych_matrixF2[i,129] <- HayBCatAzF2
    }
    if (class(HayBCatBRawF2) == 'NULL') {
      neuropsych_matrixF2[i,130] <- 'N/A'
    } else {
      neuropsych_matrixF2[i,130] <- HayBCatBRawF2
    }
    if (class(HayBCatBzF2) == 'NULL') {
      neuropsych_matrixF2[i,131] <- 'N/A'
    } else {
      neuropsych_matrixF2[i,131] <- HayBCatBzF2
    }
  }
  
  #For F3: 
  if (taskF3 == 'RUN') {
    pt_dataF3 <- read.xlsx(xlsx_filename, sheetName = 'F3')
    
    #find participant details (e.g ID, age, group) and add to matrix
    #ID
    pt_IDF3 <- pt_dataF3$Raw.1[grep("F3", pt_dataF3$Raw.1) ]
    pt_IDF3 <- paste('ADPRC_', pt_IDF3, sep='') #add in 'ADPRC' ID tag at the beginning. 
    neuropsych_matrixF3[i,1] <- pt_IDF3
    #Age
    pt_Age_locF3 <- which(pt_dataF3 == 'Age', arr.ind=TRUE) #find participant age location on the sheet
    pt_AgeF3 <- pt_dataF3[pt_Age_locF3[1], pt_Age_locF3[2]+1]
    neuropsych_matrixF3[i,2] <- pt_AgeF3
    #Education Level
    pt_EduLvl_locF3 <- which(pt_dataF3 == "Ed'n", arr.ind=TRUE) #find participant education level location on the sheet
    pt_EduLvlF3 <- pt_dataF3[pt_EduLvl_locF3[1], pt_EduLvl_locF3[2]+1]
    neuropsych_matrixF3[i,3] <- pt_EduLvlF3
    #Group status
    pt_Group_locF3 <- which(pt_dataF3 == 'Dx', arr.ind=TRUE) #find participant group status location on the sheet
    pt_GroupF3 <- pt_dataF3[pt_Group_locF3[1], pt_Group_locF3[2]+1]
    neuropsych_matrixF3[i,4] <- pt_GroupF3
    #Ethnicity
    pt_ethnicity_locF3 <- which(pt_dataF3 == 'Ethnicity', arr.ind=TRUE) #find participant ethnicity location on the sheet
    pt_ethnicityF3 <- pt_dataF3[pt_ethnicity_locF3[1], pt_ethnicity_locF3[2]+1]
    neuropsych_matrixF3[i,5] <- pt_ethnicityF3
    #Neuropsych test date
    pt_testDateF3 <- as.Date(as.integer(substr(colnames(pt_dataF3[1]), 2, 6)), origin = "1899-12-30")
    neuropsych_matrixF3[i,6] <- as.character(pt_testDateF3)
    #extract the neuropsych scores (114 total) and add into a matrix
    neuropsych_matrixF3[i,7] <- pt_dataF3[1,'Raw'] #TOPF raw
    neuropsych_matrixF3[i,8] <- pt_dataF3[1,'Z'] #TOPF Z
    neuropsych_matrixF3[i,9] <- pt_dataF3[1,'Scaled'] #TOPF Scaled
    neuropsych_matrixF3[i,10] <- pt_dataF3[2,'Raw'] #DSF raw
    neuropsych_matrixF3[i,11] <- pt_dataF3[2,'Z'] #DSF Z
    neuropsych_matrixF3[i,12] <- pt_dataF3[2,'Scaled'] #DSF Scaled
    neuropsych_matrixF3[i,13] <- pt_dataF3[3,'Raw'] #DSB raw
    neuropsych_matrixF3[i,14] <- pt_dataF3[3,'Z'] #DSB Z
    neuropsych_matrixF3[i,15] <- pt_dataF3[3,'Scaled'] #DSB Scaled
    neuropsych_matrixF3[i,16] <- pt_dataF3[4,'Raw'] #Trails A raw
    neuropsych_matrixF3[i,17] <- pt_dataF3[4,'Z'] #Trails A Z
    neuropsych_matrixF3[i,18] <- pt_dataF3[4,'Scaled'] #Trails A Scaled
    neuropsych_matrixF3[i,19] <- pt_dataF3[5,'Raw'] #Trails B raw
    neuropsych_matrixF3[i,20] <- pt_dataF3[5,'Z'] #Trails B Z
    neuropsych_matrixF3[i,21] <- pt_dataF3[5,'Scaled'] #Trails B Scaled
    neuropsych_matrixF3[i,22] <- pt_dataF3[2,'Raw.1'] #Coding raw
    neuropsych_matrixF3[i,23] <- pt_dataF3[2,'Z.1'] #Coding Z
    neuropsych_matrixF3[i,24] <- pt_dataF3[2,'Scaled.1'] #Coding Scaled
    neuropsych_matrixF3[i,25] <- pt_dataF3[3,'Raw.1'] #Colour Naming raw
    neuropsych_matrixF3[i,26] <- pt_dataF3[3,'Z.1'] #Colour Naming Z
    neuropsych_matrixF3[i,27] <- pt_dataF3[3,'Scaled.1'] #Colour Naming Scaled
    neuropsych_matrixF3[i,28] <- pt_dataF3[4,'Raw.1'] #Word Reading raw
    neuropsych_matrixF3[i,29] <- pt_dataF3[4,'Z.1'] #Word Reading Z
    neuropsych_matrixF3[i,30] <- pt_dataF3[4,'Scaled.1'] #Word Reading Scaled
    neuropsych_matrixF3[i,31] <- pt_dataF3[5,'Raw.1'] #Inhibition raw
    neuropsych_matrixF3[i,32] <- pt_dataF3[5,'Z.1'] #Inhibition Z
    neuropsych_matrixF3[i,33] <- pt_dataF3[5,'Scaled.1'] #Inhibition Scaled
    neuropsych_matrixF3[i,34] <- pt_dataF3[6,'Raw'] #SYDBAT Naming raw
    neuropsych_matrixF3[i,35] <- pt_dataF3[6,'Z'] #SYDBAT Naming Z
    neuropsych_matrixF3[i,36] <- pt_dataF3[6,'Scaled'] #SYDBAT Naming Scaled
    neuropsych_matrixF3[i,37] <- pt_dataF3[7,'Raw'] #SYDBAT Comp raw
    neuropsych_matrixF3[i,38] <- pt_dataF3[7,'Z'] #SYDBAT Comp Z
    neuropsych_matrixF3[i,39] <- pt_dataF3[7,'Scaled'] #SYDBAT Comp Scaled
    neuropsych_matrixF3[i,40] <- pt_dataF3[8,'Raw'] #SYDBAT SemAss raw
    neuropsych_matrixF3[i,41] <- pt_dataF3[8,'Z'] #SYDBAT SemAss  Z
    neuropsych_matrixF3[i,42] <- pt_dataF3[8,'Scaled'] #SYDBAT SemAss Scaled
    neuropsych_matrixF3[i,43] <- pt_dataF3[6,'Raw.1'] #BNT Man raw
    neuropsych_matrixF3[i,44] <- pt_dataF3[6,'Z.1'] #BNT Man  Z
    neuropsych_matrixF3[i,45] <- pt_dataF3[6,'Scaled.1'] #BNT Man Scaled
    neuropsych_matrixF3[i,46] <- pt_dataF3[7,'Raw.1'] #BNT Ivnik raw
    neuropsych_matrixF3[i,47] <- pt_dataF3[7,'Z.1'] #BNT Ivnik  Z
    neuropsych_matrixF3[i,48] <- pt_dataF3[7,'Scaled.1'] #BNT Ivnik Scaled
    neuropsych_matrixF3[i,49] <- pt_dataF3[8,'Raw.1'] #Similarities raw
    neuropsych_matrixF3[i,50] <- pt_dataF3[8,'Z.1'] #Similarities Z
    neuropsych_matrixF3[i,51] <- pt_dataF3[8,'Scaled.1'] #Similarities Scaled
    neuropsych_matrixF3[i,52] <- pt_dataF3[9,'Raw'] #Line O raw
    neuropsych_matrixF3[i,53] <- pt_dataF3[9,'Z'] #Line O  Z
    neuropsych_matrixF3[i,54] <- pt_dataF3[9,'Scaled'] #Line O Scaled
    neuropsych_matrixF3[i,55] <- pt_dataF3[9,'Raw.1'] #BD raw
    neuropsych_matrixF3[i,56] <- pt_dataF3[9,'Z.1'] #BD Z
    neuropsych_matrixF3[i,57] <- pt_dataF3[9,'Scaled.1'] #BD Scaled
    neuropsych_matrixF3[i,58] <- pt_dataF3[10,'Raw.1'] #Matrix raw
    neuropsych_matrixF3[i,59] <- pt_dataF3[10,'Z.1'] #Matrix Z
    neuropsych_matrixF3[i,60] <- pt_dataF3[10,'Scaled.1'] #Matrix Scaled
    neuropsych_matrixF3[i,61] <- pt_dataF3[11,'Raw'] #LM I raw
    neuropsych_matrixF3[i,62] <- pt_dataF3[11,'Z'] #LM I Z
    neuropsych_matrixF3[i,63] <- pt_dataF3[11,'Scaled'] #LM I Scaled
    neuropsych_matrixF3[i,64] <- pt_dataF3[12,'Raw'] #LM II raw
    neuropsych_matrixF3[i,65] <- pt_dataF3[12,'Z'] #LM II Z
    neuropsych_matrixF3[i,66] <- pt_dataF3[12,'Scaled'] #LM II Scaled
    neuropsych_matrixF3[i,67] <- pt_dataF3[13,'Raw'] #Story AI raw
    neuropsych_matrixF3[i,68] <- pt_dataF3[13,'Z'] #Story AI Z
    neuropsych_matrixF3[i,69] <- pt_dataF3[13,'Scaled'] #Story AI Scaled
    neuropsych_matrixF3[i,70] <- pt_dataF3[14,'Raw'] #Story AII raw
    neuropsych_matrixF3[i,71] <- pt_dataF3[14,'Z'] #Story AII Z
    neuropsych_matrixF3[i,72] <- pt_dataF3[14,'Scaled'] #Story AII Scaled
    neuropsych_matrixF3[i,73] <- pt_dataF3[11,'Raw.1'] #CVLT-II Total raw
    neuropsych_matrixF3[i,74] <- pt_dataF3[11,'Z.1'] #CVLT-II Total Z
    neuropsych_matrixF3[i,75] <- pt_dataF3[11,'Scaled.1'] #CVLT-II Total Scaled
    neuropsych_matrixF3[i,76] <- pt_dataF3[12,'Raw.1'] #CVLT-II Short raw
    neuropsych_matrixF3[i,77] <- pt_dataF3[12,'Z.1'] #CVLT-II Short Z
    neuropsych_matrixF3[i,78] <- pt_dataF3[12,'Scaled.1'] #CVLT-II Short Scaled
    neuropsych_matrixF3[i,79] <- pt_dataF3[13,'Raw.1'] #CVLT-II Long raw
    neuropsych_matrixF3[i,80] <- pt_dataF3[13,'Z.1'] #CVLT-II Long Z
    neuropsych_matrixF3[i,81] <- pt_dataF3[13,'Scaled.1'] #CVLT-II Long Scaled
    neuropsych_matrixF3[i,82] <- pt_dataF3[14,'Raw.1'] #CVLT-II Recog D raw
    neuropsych_matrixF3[i,83] <- pt_dataF3[14,'Z.1'] #CVLT-II Recog D Z
    neuropsych_matrixF3[i,84] <- pt_dataF3[14,'Scaled.1'] #CVLT-II Recog D Scaled
    neuropsych_matrixF3[i,85] <- pt_dataF3[15,'Raw'] #RCFT lmm raw
    neuropsych_matrixF3[i,86] <- pt_dataF3[15,'Z'] #RCFT lmm  Z
    neuropsych_matrixF3[i,87] <- pt_dataF3[15,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF3[i,88] <- pt_dataF3[16,'Raw'] #RCFT Del raw
    neuropsych_matrixF3[i,89] <- pt_dataF3[16,'Z'] #RCFT Del  Z
    neuropsych_matrixF3[i,90] <- pt_dataF3[16,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF3[i,91] <- pt_dataF3[17,'Raw'] #RCFT Recog raw
    neuropsych_matrixF3[i,92] <- pt_dataF3[17,'Z'] #RCFT Recog  Z
    neuropsych_matrixF3[i,93] <- pt_dataF3[17,'Scaled'] #RCFT Recog Scaled
    neuropsych_matrixF3[i,94] <- pt_dataF3[15,'Raw.1'] #BVMT Total raw
    neuropsych_matrixF3[i,95] <- pt_dataF3[15,'Z.1'] #BVMT Total Z
    neuropsych_matrixF3[i,96] <- pt_dataF3[15,'Scaled.1'] #BVMT Total Scaled
    neuropsych_matrixF3[i,97] <- pt_dataF3[15,'NA.'] #BVMT Total z Gale 
    neuropsych_matrixF3[i,98] <- pt_dataF3[16,'Raw.1'] #BVMT Del raw
    neuropsych_matrixF3[i,99] <- pt_dataF3[16,'Z.1'] #BVMT Del Z
    neuropsych_matrixF3[i,100] <- pt_dataF3[16,'Scaled.1'] #BVMT Del Scaled
    neuropsych_matrixF3[i,101] <- pt_dataF3[16,'NA.'] #BVMT Del z Gale
    neuropsych_matrixF3[i,102] <- pt_dataF3[17,'Raw.1'] #BVMT Recog Dis raw
    neuropsych_matrixF3[i,103] <- pt_dataF3[17,'Z.1'] #BVMT Recog Dis Z
    neuropsych_matrixF3[i,104] <- pt_dataF3[17,'Scaled.1'] #BVMT Recog Dis Scaled
    neuropsych_matrixF3[i,105] <- pt_dataF3[17,'NA.'] #BVMT Recog Dis z Gale
    neuropsych_matrixF3[i,106] <- pt_dataF3[18,'Raw'] #RCFT Copy raw
    neuropsych_matrixF3[i,107] <- pt_dataF3[18,'Z'] #RCFT Copy Z
    neuropsych_matrixF3[i,108] <- pt_dataF3[18,'Scaled'] #RCFT Copy Scaled
    neuropsych_matrixF3[i,109] <- pt_dataF3[19,'Raw'] #Letter Fluency raw
    neuropsych_matrixF3[i,110] <- pt_dataF3[19,'Z'] #Letter Fluency Z
    neuropsych_matrixF3[i,111] <- pt_dataF3[19,'Scaled'] #Letter Fluency Scaled
    neuropsych_matrixF3[i,112] <- pt_dataF3[20,'Raw'] #Category Fluency raw
    neuropsych_matrixF3[i,113] <- pt_dataF3[20,'Z'] #Category Fluency Z
    neuropsych_matrixF3[i,114] <- pt_dataF3[20,'Scaled'] #Category Fluency Scaled
    #neuropsych_matrixF3[i,115] <- pt_dataF3[18,'Raw'] #Fluency raw
    #neuropsych_matrixF3[i,116] <- pt_dataF3[18,'Z'] #Fluency Z
    #neuropsych_matrixF3[i,117] <- pt_dataF3[18,'Scaled'] #Fluency Scaled
    neuropsych_matrixF3[i,118] <- pt_dataF3[21,'Raw'] #C/W Inhibition raw
    neuropsych_matrixF3[i,119] <- pt_dataF3[21,'Z'] #C/W Inhibition Z
    neuropsych_matrixF3[i,120] <- pt_dataF3[21,'Scaled'] #C/W Inhibition Scaled
    neuropsych_matrixF3[i,121] <- pt_dataF3[19, 'Raw.1'] #Switching raw
    neuropsych_matrixF3[i,122] <- pt_dataF3[19, 'Z.1'] #Switching Z
    neuropsych_matrixF3[i,123] <- pt_dataF3[19, 'Scaled.1'] #Switching Scaled
    Hay_locF3 <- which(pt_dataF3 == 'Hayling B', arr.ind=TRUE) #find participant Hayling scores location on the sheet
    HayBTime1RawF3 <- pt_dataF3[Hay_locF3[1]+1, Hay_locF3[2]+1] #Hayling B Time 1 raw
    HayBTime1zF3 <- pt_dataF3[Hay_locF3[1]+1, Hay_locF3[2]+2] #Hayling B Time 1 z
    HayBTime2RawF3 <- pt_dataF3[Hay_locF3[1]+2, Hay_locF3[2]+1] #Hayling B Time 2 raw
    HayBTime2zF3 <- pt_dataF3[Hay_locF3[1]+2, Hay_locF3[2]+2] #Hayling B Time 2 z
    HayBCatARawF3 <- pt_dataF3[Hay_locF3[1]+3, Hay_locF3[2]+1] #Hayling B Cat A raw
    HayBCatAzF3 <- pt_dataF3[Hay_locF3[1]+3, Hay_locF3[2]+2] #Hayling B Cat A z
    HayBCatBRawF3 <- pt_dataF3[Hay_locF3[1]+4, Hay_locF3[2]+1] #Hayling B Cat B raw
    HayBCatBzF3 <- pt_dataF3[Hay_locF3[1]+4, Hay_locF3[2]+2] #Hayling B Cat B z
    if (class(HayBTime1RawF3) == 'NULL') {
      neuropsych_matrixF3[i,124] <- 'N/A'
    } else {
      neuropsych_matrixF3[i,124] <- HayBTime1RawF3
    }
    if (class(HayBTime1zF3) == 'NULL') {
      neuropsych_matrixF3[i,125] <- 'N/A'
    } else {
      neuropsych_matrixF3[i,125] <- HayBTime1zF3
    }
    if (class(HayBTime2RawF3) == 'NULL') {
      neuropsych_matrixF3[i,126] <- 'N/A'
    } else {
      neuropsych_matrixF3[i,126] <- HayBTime2RawF3
    }
    if (class(HayBTime2zF3) == 'NULL') {
      neuropsych_matrixF3[i,127] <- 'N/A'
    } else {
      neuropsych_matrixF3[i,127] <- HayBTime2zF3
    }
    if (class(HayBCatARawF3) == 'NULL') {
      neuropsych_matrixF3[i,128] <- 'N/A'
    } else {
      neuropsych_matrixF3[i,128] <- HayBCatARawF3
    }
    if (class(HayBCatAzF3) == 'NULL') {
      neuropsych_matrixF3[i,129] <- 'N/A'
    } else {
      neuropsych_matrixF3[i,129] <- HayBCatAzF3
    }
    if (class(HayBCatBRawF3) == 'NULL') {
      neuropsych_matrixF3[i,130] <- 'N/A'
    } else {
      neuropsych_matrixF3[i,130] <- HayBCatBRawF3
    }
    if (class(HayBCatBzF3) == 'NULL') {
      neuropsych_matrixF3[i,131] <- 'N/A'
    } else {
      neuropsych_matrixF3[i,131] <- HayBCatBzF3
    }
  } 
  
  #For F4
  if (taskF4 == 'RUN') {
    pt_dataF4 <- read.xlsx(xlsx_filename, sheetName = 'F4')
    
    #find participant details (e.g ID, age, group) and add to matrix
    #ID
    pt_IDF4 <- pt_dataF4$Raw.1[grep("F4", pt_dataF4$Raw.1) ]
    pt_IDF4 <- paste('ADPRC_', pt_IDF4, sep='') #add in 'ADPRC' ID tag at the beginning. 
    neuropsych_matrixF4[i,1] <- pt_IDF4
    #Age
    pt_Age_locF4 <- which(pt_dataF4 == 'Age', arr.ind=TRUE) #find participant age location on the sheet
    pt_AgeF4 <- pt_dataF4[pt_Age_locF4[1], pt_Age_locF4[2]+1]
    neuropsych_matrixF4[i,2] <- pt_AgeF4
    #Education Level
    pt_EduLvl_locF4 <- which(pt_dataF4 == "Ed'n", arr.ind=TRUE) #find participant education level location on the sheet
    pt_EduLvlF4 <- pt_dataF4[pt_EduLvl_locF4[1], pt_EduLvl_locF4[2]+1]
    neuropsych_matrixF4[i,3] <- pt_EduLvlF4
    #Group status
    pt_Group_locF4 <- which(pt_dataF4 == 'Dx', arr.ind=TRUE) #find participant group status location on the sheet
    pt_GroupF4 <- pt_dataF4[pt_Group_locF4[1], pt_Group_locF4[2]+1]
    neuropsych_matrixF4[i,4] <- pt_GroupF4
    #Ethnicity
    pt_ethnicity_locF4 <- which(pt_dataF4 == 'Ethnicity', arr.ind=TRUE) #find participant ethnicity location on the sheet
    pt_ethnicityF4 <- pt_dataF0[pt_ethnicity_locF4[1], pt_ethnicity_locF4[2]+1]
    neuropsych_matrixF0[i,5] <- pt_ethnicityF4
    #Neuropsych test date
    pt_testDateF4 <- as.Date(as.integer(substr(colnames(pt_dataF4[1]), 2, 6)), origin = "1899-12-30")
    neuropsych_matrixF4[i,6] <- as.character(pt_testDateF4)
    #extract the neuropsych scores (114 total) and add into a matrix
    neuropsych_matrixF4[i,7] <- pt_dataF4[1,'Raw'] #TOPF raw
    neuropsych_matrixF4[i,8] <- pt_dataF4[1,'Z'] #TOPF Z
    neuropsych_matrixF4[i,9] <- pt_dataF4[1,'Scaled'] #TOPF Scaled
    neuropsych_matrixF4[i,10] <- pt_dataF4[2,'Raw'] #DSF raw
    neuropsych_matrixF4[i,11] <- pt_dataF4[2,'Z'] #DSF Z
    neuropsych_matrixF4[i,12] <- pt_dataF4[2,'Scaled'] #DSF Scaled
    neuropsych_matrixF4[i,13] <- pt_dataF4[3,'Raw'] #DSB raw
    neuropsych_matrixF4[i,14] <- pt_dataF4[3,'Z'] #DSB Z
    neuropsych_matrixF4[i,15] <- pt_dataF4[3,'Scaled'] #DSB Scaled
    neuropsych_matrixF4[i,16] <- pt_dataF4[4,'Raw'] #Trails A raw
    neuropsych_matrixF4[i,17] <- pt_dataF4[4,'Z'] #Trails A Z
    neuropsych_matrixF4[i,18] <- pt_dataF4[4,'Scaled'] #Trails A Scaled
    neuropsych_matrixF4[i,19] <- pt_dataF4[5,'Raw'] #Trails B raw
    neuropsych_matrixF4[i,20] <- pt_dataF4[5,'Z'] #Trails B Z
    neuropsych_matrixF4[i,21] <- pt_dataF4[5,'Scaled'] #Trails B Scaled
    neuropsych_matrixF4[i,22] <- pt_dataF4[2,'Raw.1'] #Coding raw
    neuropsych_matrixF4[i,23] <- pt_dataF4[2,'Z.1'] #Coding Z
    neuropsych_matrixF4[i,24] <- pt_dataF4[2,'Scaled.1'] #Coding Scaled
    neuropsych_matrixF4[i,25] <- pt_dataF4[3,'Raw.1'] #Colour Naming raw
    neuropsych_matrixF4[i,26] <- pt_dataF4[3,'Z.1'] #Colour Naming Z
    neuropsych_matrixF4[i,27] <- pt_dataF4[3,'Scaled.1'] #Colour Naming Scaled
    neuropsych_matrixF4[i,28] <- pt_dataF4[4,'Raw.1'] #Word Reading raw
    neuropsych_matrixF4[i,29] <- pt_dataF4[4,'Z.1'] #Word Reading Z
    neuropsych_matrixF4[i,30] <- pt_dataF4[4,'Scaled.1'] #Word Reading Scaled
    neuropsych_matrixF4[i,31] <- pt_dataF4[5,'Raw.1'] #Inhibition raw
    neuropsych_matrixF4[i,32] <- pt_dataF4[5,'Z.1'] #Inhibition Z
    neuropsych_matrixF4[i,33] <- pt_dataF4[5,'Scaled.1'] #Inhibition Scaled
    neuropsych_matrixF4[i,34] <- pt_dataF4[6,'Raw'] #SYDBAT Naming raw
    neuropsych_matrixF4[i,35] <- pt_dataF4[6,'Z'] #SYDBAT Naming Z
    neuropsych_matrixF4[i,36] <- pt_dataF4[6,'Scaled'] #SYDBAT Naming Scaled
    neuropsych_matrixF4[i,37] <- pt_dataF4[7,'Raw'] #SYDBAT Comp raw
    neuropsych_matrixF4[i,38] <- pt_dataF4[7,'Z'] #SYDBAT Comp Z
    neuropsych_matrixF4[i,39] <- pt_dataF4[7,'Scaled'] #SYDBAT Comp Scaled
    neuropsych_matrixF4[i,40] <- pt_dataF4[8,'Raw'] #SYDBAT SemAss raw
    neuropsych_matrixF4[i,41] <- pt_dataF4[8,'Z'] #SYDBAT SemAss  Z
    neuropsych_matrixF4[i,42] <- pt_dataF4[8,'Scaled'] #SYDBAT SemAss Scaled
    neuropsych_matrixF4[i,43] <- pt_dataF4[6,'Raw.1'] #BNT Man raw
    neuropsych_matrixF4[i,44] <- pt_dataF4[6,'Z.1'] #BNT Man  Z
    neuropsych_matrixF4[i,45] <- pt_dataF4[6,'Scaled.1'] #BNT Man Scaled
    neuropsych_matrixF4[i,46] <- pt_dataF4[7,'Raw.1'] #BNT Ivnik raw
    neuropsych_matrixF4[i,47] <- pt_dataF4[7,'Z.1'] #BNT Ivnik  Z
    neuropsych_matrixF4[i,48] <- pt_dataF4[7,'Scaled.1'] #BNT Ivnik Scaled
    neuropsych_matrixF4[i,49] <- pt_dataF4[8,'Raw.1'] #Similarities raw
    neuropsych_matrixF4[i,50] <- pt_dataF4[8,'Z.1'] #Similarities Z
    neuropsych_matrixF4[i,51] <- pt_dataF4[8,'Scaled.1'] #Similarities Scaled
    neuropsych_matrixF4[i,52] <- pt_dataF4[9,'Raw'] #Line O raw
    neuropsych_matrixF4[i,53] <- pt_dataF4[9,'Z'] #Line O  Z
    neuropsych_matrixF4[i,54] <- pt_dataF4[9,'Scaled'] #Line O Scaled
    neuropsych_matrixF4[i,55] <- pt_dataF4[9,'Raw.1'] #BD raw
    neuropsych_matrixF4[i,56] <- pt_dataF4[9,'Z.1'] #BD Z
    neuropsych_matrixF4[i,57] <- pt_dataF4[9,'Scaled.1'] #BD Scaled
    neuropsych_matrixF4[i,58] <- pt_dataF4[10,'Raw.1'] #Matrix raw
    neuropsych_matrixF4[i,59] <- pt_dataF4[10,'Z.1'] #Matrix Z
    neuropsych_matrixF4[i,60] <- pt_dataF4[10,'Scaled.1'] #Matrix Scaled
    neuropsych_matrixF4[i,61] <- pt_dataF4[11,'Raw'] #LM I raw
    neuropsych_matrixF4[i,62] <- pt_dataF4[11,'Z'] #LM I Z
    neuropsych_matrixF4[i,63] <- pt_dataF4[11,'Scaled'] #LM I Scaled
    neuropsych_matrixF4[i,64] <- pt_dataF4[12,'Raw'] #LM II raw
    neuropsych_matrixF4[i,65] <- pt_dataF4[12,'Z'] #LM II Z
    neuropsych_matrixF4[i,66] <- pt_dataF4[12,'Scaled'] #LM II Scaled
    neuropsych_matrixF4[i,67] <- pt_dataF4[13,'Raw'] #Story AI raw
    neuropsych_matrixF4[i,68] <- pt_dataF4[13,'Z'] #Story AI Z
    neuropsych_matrixF4[i,69] <- pt_dataF4[13,'Scaled'] #Story AI Scaled
    neuropsych_matrixF4[i,70] <- pt_dataF4[14,'Raw'] #Story AII raw
    neuropsych_matrixF4[i,71] <- pt_dataF4[14,'Z'] #Story AII Z
    neuropsych_matrixF4[i,72] <- pt_dataF4[14,'Scaled'] #Story AII Scaled
    neuropsych_matrixF4[i,73] <- pt_dataF4[11,'Raw.1'] #CVLT-II Total raw
    neuropsych_matrixF4[i,74] <- pt_dataF4[11,'Z.1'] #CVLT-II Total Z
    neuropsych_matrixF4[i,75] <- pt_dataF4[11,'Scaled.1'] #CVLT-II Total Scaled
    neuropsych_matrixF4[i,76] <- pt_dataF4[12,'Raw.1'] #CVLT-II Short raw
    neuropsych_matrixF4[i,77] <- pt_dataF4[12,'Z.1'] #CVLT-II Short Z
    neuropsych_matrixF4[i,78] <- pt_dataF4[12,'Scaled.1'] #CVLT-II Short Scaled
    neuropsych_matrixF4[i,79] <- pt_dataF4[13,'Raw.1'] #CVLT-II Long raw
    neuropsych_matrixF4[i,80] <- pt_dataF4[13,'Z.1'] #CVLT-II Long Z
    neuropsych_matrixF4[i,81] <- pt_dataF4[13,'Scaled.1'] #CVLT-II Long Scaled
    neuropsych_matrixF4[i,82] <- pt_dataF4[14,'Raw.1'] #CVLT-II Recog D raw
    neuropsych_matrixF4[i,83] <- pt_dataF4[14,'Z.1'] #CVLT-II Recog D Z
    neuropsych_matrixF4[i,84] <- pt_dataF4[14,'Scaled.1'] #CVLT-II Recog D Scaled
    neuropsych_matrixF4[i,85] <- pt_dataF4[15,'Raw'] #RCFT lmm raw
    neuropsych_matrixF4[i,86] <- pt_dataF4[15,'Z'] #RCFT lmm  Z
    neuropsych_matrixF4[i,87] <- pt_dataF4[15,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF4[i,88] <- pt_dataF4[16,'Raw'] #RCFT Del raw
    neuropsych_matrixF4[i,89] <- pt_dataF4[16,'Z'] #RCFT Del  Z
    neuropsych_matrixF4[i,90] <- pt_dataF4[16,'Scaled'] #RCFT Del Scaled
    neuropsych_matrixF4[i,91] <- pt_dataF4[17,'Raw'] #RCFT Recog raw
    neuropsych_matrixF4[i,92] <- pt_dataF4[17,'Z'] #RCFT Recog  Z
    neuropsych_matrixF4[i,93] <- pt_dataF4[17,'Scaled'] #RCFT Recog Scaled
    neuropsych_matrixF4[i,94] <- pt_dataF4[15,'Raw.1'] #BVMT Total raw
    neuropsych_matrixF4[i,95] <- pt_dataF4[15,'Z.1'] #BVMT Total Z
    neuropsych_matrixF4[i,96] <- pt_dataF4[15,'Scaled.1'] #BVMT Total Scaled
    neuropsych_matrixF4[i,97] <- pt_dataF4[15,'NA.'] #BVMT Total z Gale 
    neuropsych_matrixF4[i,98] <- pt_dataF4[16,'Raw.1'] #BVMT Del raw
    neuropsych_matrixF4[i,99] <- pt_dataF4[16,'Z.1'] #BVMT Del Z
    neuropsych_matrixF4[i,100] <- pt_dataF4[16,'Scaled.1'] #BVMT Del Scaled
    neuropsych_matrixF4[i,101] <- pt_dataF4[16,'NA.'] #BVMT Del z Gale
    neuropsych_matrixF4[i,102] <- pt_dataF4[17,'Raw.1'] #BVMT Recog Dis raw
    neuropsych_matrixF4[i,103] <- pt_dataF4[17,'Z.1'] #BVMT Recog Dis Z
    neuropsych_matrixF4[i,104] <- pt_dataF4[17,'Scaled.1'] #BVMT Recog Dis Scaled
    neuropsych_matrixF4[i,105] <- pt_dataF4[17,'NA.'] #BVMT Recog Dis z Gale
    neuropsych_matrixF4[i,106] <- pt_dataF4[18,'Raw'] #RCFT Copy raw
    neuropsych_matrixF4[i,107] <- pt_dataF4[18,'Z'] #RCFT Copy Z
    neuropsych_matrixF4[i,108] <- pt_dataF4[18,'Scaled'] #RCFT Copy Scaled
    neuropsych_matrixF4[i,109] <- pt_dataF4[19,'Raw'] #Letter Fluency raw
    neuropsych_matrixF4[i,110] <- pt_dataF4[19,'Z'] #Letter Fluency Z
    neuropsych_matrixF4[i,111] <- pt_dataF4[19,'Scaled'] #Letter Fluency Scaled
    neuropsych_matrixF4[i,112] <- pt_dataF4[20,'Raw'] #Category Fluency raw
    neuropsych_matrixF4[i,113] <- pt_dataF4[20,'Z'] #Category Fluency Z
    neuropsych_matrixF4[i,114] <- pt_dataF4[20,'Scaled'] #Category Fluency Scaled
    #neuropsych_matrixF4[i,115] <- pt_dataF4[18,'Raw'] #Fluency raw
    #neuropsych_matrixF4[i,116] <- pt_dataF4[18,'Z'] #Fluency Z
    #neuropsych_matrixF4[i,117] <- pt_dataF4[18,'Scaled'] #Fluency Scaled
    neuropsych_matrixF4[i,118] <- pt_dataF4[21,'Raw'] #C/W Inhibition raw
    neuropsych_matrixF4[i,119] <- pt_dataF4[21,'Z'] #C/W Inhibition Z
    neuropsych_matrixF4[i,120] <- pt_dataF4[21,'Scaled'] #C/W Inhibition Scaled
    neuropsych_matrixF4[i,121] <- pt_dataF4[19, 'Raw.1'] #Switching raw
    neuropsych_matrixF4[i,122] <- pt_dataF4[19, 'Z.1'] #Switching Z
    neuropsych_matrixF4[i,123] <- pt_dataF4[19, 'Scaled.1'] #Switching Scaled
    Hay_locF4 <- which(pt_dataF4 == 'Hayling B', arr.ind=TRUE) #find participant Hayling scores location on the sheet
    HayBTime1RawF4 <- pt_dataF4[Hay_locF4[1]+1, Hay_locF4[2]+1] #Hayling B Time 1 raw
    HayBTime1zF4 <- pt_dataF4[Hay_locF4[1]+1, Hay_locF4[2]+2] #Hayling B Time 1 z
    HayBTime2RawF4 <- pt_dataF4[Hay_locF4[1]+2, Hay_locF4[2]+1] #Hayling B Time 2 raw
    HayBTime2zF4 <- pt_dataF4[Hay_locF4[1]+2, Hay_locF4[2]+2] #Hayling B Time 2 z
    HayBCatARawF4 <- pt_dataF4[Hay_locF4[1]+3, Hay_locF4[2]+1] #Hayling B Cat A raw
    HayBCatAzF4 <- pt_dataF4[Hay_locF4[1]+3, Hay_locF4[2]+2] #Hayling B Cat A z
    HayBCatBRawF4 <- pt_dataF4[Hay_locF4[1]+4, Hay_locF4[2]+1] #Hayling B Cat B raw
    HayBCatBzF4 <- pt_dataF4[Hay_locF4[1]+4, Hay_locF4[2]+2] #Hayling B Cat B z
    if (class(HayBTime1RawF4) == 'NULL') {
      neuropsych_matrixF4[i,124] <- 'N/A'
    } else {
      neuropsych_matrixF4[i,124] <- HayBTime1RawF4
    }
    if (class(HayBTime1zF4) == 'NULL') {
      neuropsych_matrixF4[i,125] <- 'N/A'
    } else {
      neuropsych_matrixF4[i,125] <- HayBTime1zF4
    }
    if (class(HayBTime2RawF4) == 'NULL') {
      neuropsych_matrixF4[i,126] <- 'N/A'
    } else {
      neuropsych_matrixF4[i,126] <- HayBTime2RawF4
    }
    if (class(HayBTime2zF4) == 'NULL') {
      neuropsych_matrixF4[i,127] <- 'N/A'
    } else {
      neuropsych_matrixF4[i,127] <- HayBTime2zF4
    }
    if (class(HayBCatARawF4) == 'NULL') {
      neuropsych_matrixF4[i,128] <- 'N/A'
    } else {
      neuropsych_matrixF4[i,128] <- HayBCatARawF4
    }
    if (class(HayBCatAzF4) == 'NULL') {
      neuropsych_matrixF4[i,129] <- 'N/A'
    } else {
      neuropsych_matrixF4[i,129] <- HayBCatAzF4
    }
    if (class(HayBCatBRawF4) == 'NULL') {
      neuropsych_matrixF4[i,130] <- 'N/A'
    } else {
      neuropsych_matrixF4[i,130] <- HayBCatBRawF4
    }
    if (class(HayBCatBzF4) == 'NULL') {
      neuropsych_matrixF4[i,131] <- 'N/A'
    } else {
      neuropsych_matrixF4[i,131] <- HayBCatBzF4
    }
  }
  
  #Upon finishing with the participant file, move up one folder, to the main participant list folder
  setwd('..')
  
}

#combine all data into one matrix
neuropsych_matrixAll <- rbind(neuropsych_matrixF0, neuropsych_matrixF1, neuropsych_matrixF2, neuropsych_matrixF3, neuropsych_matrixF4)

#write data to the 'master excel sheet' containing all data
write.xlsx(neuropsych_matrixF0, 'DPRC_neuropsych_data.xlsx', sheetName='F0', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF1, 'DPRC_neuropsych_data.xlsx', sheetName='F1', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF2, 'DPRC_neuropsych_data.xlsx', sheetName='F2', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF3, 'DPRC_neuropsych_data.xlsx', sheetName='F3', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF4, 'DPRC_neuropsych_data.xlsx', sheetName='F4', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixAll, 'DPRC_neuropsych_data.xlsx', sheetName='All', row.names=FALSE, append=TRUE)