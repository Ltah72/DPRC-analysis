#Organise DPRC neuropsychological data into one spreadsheet, containing all 
#participant values. This is the NEW version of the participant files, 
#including the REDcap summary sheet, in which the data are extracted from. 

#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 13/01/22


#install package manager to handle all other package installations and dependencies. 
#install.packages("pacman")
pacman::p_load(xlsx, stringr, data.table)

#load in libraries
library("xlsx")
library("stringr")
library("data.table") 

# choose & set to your directory. This is where each of your participant's 
#files should be. 
setwd("C:/Users/ltah262/Desktop/DPRC_data_organisation/Participant Files/")  
#setwd("Z:/DPRC Neuropsychologists/Participant Files/")  

#load the file names from directory into the work space
files_all_whole <- list.files() 
files_all <- files_all_whole
#files_all <- files_all_whole[4:198]
#files_all <- files_all[-c(74, 75, 131)]#exclude certain participants (e.g. don't have excel file)

#vector for time points: 
timepoints <- c('F0','F1','F2','F3','F4','F5','F6','F8','F10') #The 9 timepoints

#put in data for all available time points
for (j in 1:length(timepoints)){                               
  
  current_matrix <- matrix(nrow = length(files_all), ncol=133)
  colnames(current_matrix) <- c('ParticipantID','Date of Assessment','TOPF Raw','TOPF z-score','TOPF St Score','Agreed premorbid z-score','Education','DF max span','DF max span z-score','DB max span','DB max span z-score','DF total raw','DF total z-score','DF total scaled score','DB total raw','DB total z-score','DB total scaled score','Trails A raw','Trails A z-score','Trails A errors','Trails B raw','Trails B z-score','Trails B errors','Coding raw','Coding z-score','Coding scaled score','Stroop Color Naming raw','Stroop Color Naming z-score','Stroop Color Naming scaled score','Stroop Word raw','Stroop Word z-score','Stroop Word scaled score','Stroop Inhibition raw','Stroop Inhibition z-score','Stroop Inhibition scaled score','Stroop Inhibition errors','Stroop Inhibition errors z-score','Stroop Inhibition errors scaled score','SYDBAT Naming raw','SYDBAT Naming z-score','SYDBAT Naming cut-off','SYDBAT Comprehension raw','SYDBAT Comprehension z-score','SYDBAT Comprehension cut-off','SYDBAT Semantic raw','SYDBAT Semantic z-score','SYDBAT Semantic cut-off','Boston 30 raw','Boston 60 raw','Boston Manual z-score','Boston Ivnik z-score','Boston Ivnik scaled score','JOL raw','JOL z-score','JOL scaled score','RCFT Copy raw','RCFT Copy z-score','RCFT Copy Time','RCFT Copy Time z-score','RCFT Immediate raw','RCFT Immediate z-score','RCFT Delayed raw','RCFT Delayed z-score','RCFT Recognition raw','RCFT Recognition z-score','BVMT Trail 1','BVMT Trail 2','BVMT Trail 3','BVMT Total raw','BVMT Total manual z-score','BVMT Delayed raw','BVMT Delayed manual z-score','BVMT Recognition Discrimination raw','BVMT Recognition Discrimination manual z-score','BVMT Total gale z-score','BVMT Delayed gale z-score','LM I raw','LM I z-score','LM I scaled score','LM II raw','LM II z-score','LM II scaled score','Story A I','Story A II','CVLT Trial One','CVLT Trial One z-score','CVLT Trial Two','CVLT Trial Three','CVLT Trial Four','CVLT Trial Five','CVLT Total raw','CVLT Total z-score','CVLT Short Delay raw','CVLT Short Delay raw z-score','CVLT Long Delay raw','CVLT Long Delay raw z-score','CVLT Recog Discrimination raw','CVLT Recog Discrimination raw z-score','Letter Fluency F','Letter Fluency A','Letter Fluency S','Letter Fluency Total raw','Letter Fluency z-score','Letter Fluency scaled score','Category Fluency Animals','Category Fluency Boys names','Category Fluency Total raw','Category Fluency Total z-score','Category Fluency Total scaled score','Category Switching Accuracy raw','Category Switching Accuracy z-score','Category Switching Accuracy scaled score','Hayling Manual Overall Scaled Score','Hayling Bielak Time 1 raw','Hayling Bielak Time 1 z-score','Hayling Bielak Time 2 raw','Hayling Bielak Time 2 z-score','Hayling Category A Errors raw','Hayling Category A Errors z-score','Hayling Category B Errors raw','Hayling Category A Errors z-score','BD raw','BD z-score','BD scaled score','BD No Time Bonus raw','BD No Time Bonus z-score','BD No Time Bonus scaled score','Matrix raw','Matrix z-score','Matrix scaled score','Similaritites raw','Similaritites z-score','Similaritites scaled score')
  
  #loop through every participant
  for(i in 1:length(files_all)){
    
    #read in each of the participant files 
    setwd(files_all[i])
    pt_directory <- list.files(pattern = '^[^~]') #pattern will ignore the temp files 
    xlsx_filename <- pt_directory[grepl('.xlsx', pt_directory) & !grepl('.xlsx.sb', pt_directory) & grepl('NEW', pt_directory)] #choose the NEW excel file and ignore the temp folder named with .xlsx
    
    #if the sheet of the specific timepoint exists, then read in the data from the file.
    pt_data <- read.xlsx(xlsx_filename, sheetName = 'REDcap')
    
    #for F0
    if (grepl('ADPRC', colnames(pt_data))[j+1]==TRUE) {
      
      #find participant details (e.g ID, test date) and add to matrix
      #ID
      pt_ID <- colnames(pt_data[j+1])
      current_matrix[i,'ParticipantID'] <- pt_ID 
      #Neuropsych test date
      pt_testDate <- as.Date(as.integer(substr(pt_data[2,pt_ID], 1, 6)), origin = "1899-12-30")
      current_matrix[i,'Date of Assessment'] <- as.character(pt_testDate)
      
      #extract the neuropsych scores (130 total) and add into a matrix
      current_matrix[i,'TOPF Raw'] <- pt_data[3,pt_ID] #TOPF raw
      current_matrix[i,'TOPF z-score'] <- pt_data[4,pt_ID] #TOPF Z
      current_matrix[i,'TOPF St Score'] <- pt_data[5,pt_ID] #TOPF Scaled
      current_matrix[i,'Agreed premorbid z-score'] <- pt_data[6,pt_ID] #Agreed premorbid z-score
      current_matrix[i,'Education'] <- pt_data[7,pt_ID] #Education
      current_matrix[i,'DF max span'] <- pt_data[10,pt_ID] 
      current_matrix[i,'DF max span z-score'] <- pt_data[11,pt_ID]
      current_matrix[i,'DB max span'] <- pt_data[12,pt_ID] 
      current_matrix[i,'DB max span z-score'] <- pt_data[13,pt_ID]
      current_matrix[i,'DF total raw'] <- pt_data[14,pt_ID] 
      current_matrix[i,'DF total z-score'] <- pt_data[15,pt_ID] 
      current_matrix[i,'DF total scaled score'] <- pt_data[16,pt_ID] 
      current_matrix[i,'DB total raw'] <- pt_data[17,pt_ID] 
      current_matrix[i,'DB total z-score'] <- pt_data[18,pt_ID] 
      current_matrix[i,'DB total scaled score'] <- pt_data[19,pt_ID] 
      current_matrix[i,'Trails A raw'] <- pt_data[20,pt_ID] 
      current_matrix[i,'Trails A z-score'] <- pt_data[21,pt_ID] 
      current_matrix[i,'Trails A errors'] <- pt_data[22,pt_ID] 
      current_matrix[i,'Trails B raw'] <- pt_data[23,pt_ID] 
      current_matrix[i,'Trails B z-score'] <- pt_data[24,pt_ID] 
      current_matrix[i,'Trails B errors'] <- pt_data[25,pt_ID] 
      current_matrix[i,'Coding raw'] <- pt_data[26,pt_ID] 
      current_matrix[i,'Coding z-score'] <- pt_data[27,pt_ID] 
      current_matrix[i,'Coding scaled score'] <- pt_data[28,pt_ID] 
      current_matrix[i,'Stroop Color Naming raw'] <- pt_data[29,pt_ID] 
      current_matrix[i,'Stroop Color Naming z-score'] <- pt_data[30,pt_ID] 
      current_matrix[i,'Stroop Color Naming scaled score'] <- pt_data[31,pt_ID] 
      current_matrix[i,'Stroop Word raw'] <- pt_data[32,pt_ID] 
      current_matrix[i,'Stroop Word z-score'] <- pt_data[33,pt_ID] 
      current_matrix[i,'Stroop Word scaled score'] <- pt_data[34,pt_ID] 
      current_matrix[i,'Stroop Inhibition raw'] <- pt_data[35,pt_ID] 
      current_matrix[i,'Stroop Inhibition z-score'] <- pt_data[36,pt_ID] 
      current_matrix[i,'Stroop Inhibition scaled score'] <- pt_data[37,pt_ID] 
      current_matrix[i,'Stroop Inhibition errors'] <- pt_data[38,pt_ID] 
      current_matrix[i,'Stroop Inhibition errors z-score'] <- pt_data[39,pt_ID] 
      current_matrix[i,'Stroop Inhibition errors scaled score'] <- pt_data[40,pt_ID] 
      current_matrix[i,'SYDBAT Naming raw'] <- pt_data[43,pt_ID] 
      current_matrix[i,'SYDBAT Naming z-score'] <- pt_data[44,pt_ID] 
      current_matrix[i,'SYDBAT Naming cut-off'] <- pt_data[45,pt_ID] 
      current_matrix[i,'SYDBAT Comprehension raw'] <- pt_data[46,pt_ID]
      current_matrix[i,'SYDBAT Comprehension z-score'] <- pt_data[47,pt_ID]
      current_matrix[i,'SYDBAT Comprehension cut-off'] <- pt_data[48,pt_ID] 
      current_matrix[i,'SYDBAT Semantic raw'] <- pt_data[49,pt_ID] 
      current_matrix[i,'SYDBAT Semantic z-score'] <- pt_data[50,pt_ID] 
      current_matrix[i,'SYDBAT Semantic cut-off'] <- pt_data[51,pt_ID] 
      current_matrix[i,'Boston 30 raw'] <- pt_data[52,pt_ID] 
      current_matrix[i,'Boston 60 raw'] <- pt_data[53,pt_ID] 
      current_matrix[i,'Boston Manual z-score'] <- pt_data[54,pt_ID] 
      current_matrix[i,'Boston Ivnik z-score'] <- pt_data[55,pt_ID] 
      current_matrix[i,'Boston Ivnik scaled score'] <- pt_data[56,pt_ID] 
      current_matrix[i,'JOL raw'] <- pt_data[59,pt_ID] 
      current_matrix[i,'JOL z-score'] <- pt_data[60,pt_ID] 
      current_matrix[i,'JOL scaled score'] <- pt_data[61,pt_ID] 
      current_matrix[i,'RCFT Copy raw'] <- pt_data[62,pt_ID] 
      current_matrix[i,'RCFT Copy z-score'] <- pt_data[63,pt_ID]
      current_matrix[i,'RCFT Copy Time'] <- pt_data[64,pt_ID] 
      current_matrix[i,'RCFT Copy Time z-score'] <- pt_data[65,pt_ID] 
      current_matrix[i,'RCFT Immediate raw'] <- pt_data[66,pt_ID] 
      current_matrix[i,'RCFT Immediate z-score'] <- pt_data[67,pt_ID] 
      current_matrix[i,'RCFT Delayed raw'] <- pt_data[68,pt_ID] 
      current_matrix[i,'RCFT Delayed z-score'] <- pt_data[69,pt_ID]
      current_matrix[i,'RCFT Recognition raw'] <- pt_data[70,pt_ID] 
      current_matrix[i,'RCFT Recognition z-score'] <- pt_data[71,pt_ID]
      current_matrix[i,'BVMT Trail 1'] <- pt_data[72,pt_ID] 
      current_matrix[i,'BVMT Trail 2'] <- pt_data[73,pt_ID] 
      current_matrix[i,'BVMT Trail 3'] <- pt_data[74,pt_ID] 
      current_matrix[i,'BVMT Total raw'] <- pt_data[75,pt_ID] 
      current_matrix[i,'BVMT Total manual z-score'] <- pt_data[76,pt_ID] 
      current_matrix[i,'BVMT Delayed raw'] <- pt_data[77,pt_ID] 
      current_matrix[i,'BVMT Delayed manual z-score'] <- pt_data[78,pt_ID] 
      current_matrix[i,'BVMT Recognition Discrimination raw'] <- pt_data[79,pt_ID] 
      current_matrix[i,'BVMT Recognition Discrimination manual z-score'] <- pt_data[80,pt_ID] 
      current_matrix[i,'BVMT Total gale z-score'] <- pt_data[81,pt_ID] 
      current_matrix[i,'BVMT Delayed gale z-score'] <- pt_data[82,pt_ID] 
      current_matrix[i,'LM I raw'] <- pt_data[85,pt_ID] 
      current_matrix[i,'LM I z-score'] <- pt_data[86,pt_ID] 
      current_matrix[i,'LM I scaled score'] <- pt_data[87,pt_ID] 
      current_matrix[i,'LM II raw'] <- pt_data[88,pt_ID] 
      current_matrix[i,'LM II z-score'] <- pt_data[89,pt_ID] 
      current_matrix[i,'LM II scaled score'] <- pt_data[90,pt_ID] 
      current_matrix[i,'Story A I'] <- pt_data[91,pt_ID] 
      current_matrix[i,'Story A II'] <- pt_data[92,pt_ID] 
      current_matrix[i,'CVLT Trial One'] <- pt_data[93,pt_ID] 
      current_matrix[i,'CVLT Trial One z-score'] <- pt_data[94,pt_ID] 
      current_matrix[i,'CVLT Trial Two'] <- pt_data[95,pt_ID] 
      current_matrix[i,'CVLT Trial Three'] <- pt_data[96,pt_ID] 
      current_matrix[i,'CVLT Trial Four'] <- pt_data[97,pt_ID] 
      current_matrix[i,'CVLT Trial Five'] <- pt_data[98,pt_ID] 
      current_matrix[i,'CVLT Total raw'] <- pt_data[99,pt_ID] 
      current_matrix[i,'CVLT Total z-score'] <- pt_data[100,pt_ID] 
      current_matrix[i,'CVLT Short Delay raw'] <- pt_data[101,pt_ID] 
      current_matrix[i,'CVLT Short Delay raw z-score'] <- pt_data[102,pt_ID] 
      current_matrix[i,'CVLT Long Delay raw'] <- pt_data[103,pt_ID] 
      current_matrix[i,'CVLT Long Delay raw z-score'] <- pt_data[104,pt_ID] 
      current_matrix[i,'CVLT Recog Discrimination raw'] <- pt_data[105,pt_ID] 
      current_matrix[i,'CVLT Recog Discrimination raw z-score'] <- pt_data[106,pt_ID] 
      current_matrix[i,'Letter Fluency F'] <- pt_data[109,pt_ID] 
      current_matrix[i,'Letter Fluency A'] <- pt_data[110,pt_ID] 
      current_matrix[i,'Letter Fluency S'] <- pt_data[111,pt_ID] 
      current_matrix[i,'Letter Fluency Total raw'] <- pt_data[112,pt_ID] 
      current_matrix[i,'Letter Fluency z-score'] <- pt_data[113,pt_ID] 
      current_matrix[i,'Letter Fluency scaled score'] <- pt_data[114,pt_ID] 
      current_matrix[i,'Category Fluency Animals'] <- pt_data[115,pt_ID] 
      current_matrix[i,'Category Fluency Boys names'] <- pt_data[116,pt_ID] 
      current_matrix[i,'Category Fluency Total raw'] <- pt_data[117,pt_ID] 
      current_matrix[i,'Category Fluency Total z-score'] <- pt_data[118,pt_ID] 
      current_matrix[i,'Category Fluency Total scaled score'] <- pt_data[119,pt_ID] 
      current_matrix[i,'Category Switching Accuracy raw'] <- pt_data[120,pt_ID] 
      current_matrix[i,'Category Switching Accuracy z-score'] <- pt_data[121,pt_ID] 
      current_matrix[i,'Category Switching Accuracy scaled score'] <- pt_data[122,pt_ID] 
      current_matrix[i,'Hayling Manual Overall Scaled Score'] <- pt_data[123,pt_ID] 
      current_matrix[i,'Hayling Bielak Time 1 raw'] <- pt_data[124,pt_ID] 
      current_matrix[i,'Hayling Bielak Time 1 z-score'] <- pt_data[125,pt_ID] 
      current_matrix[i,'Hayling Bielak Time 2 raw'] <- pt_data[126,pt_ID] 
      current_matrix[i,'Hayling Bielak Time 2 z-score'] <- pt_data[127,pt_ID] 
      current_matrix[i,'Hayling Category A Errors raw'] <- pt_data[128,pt_ID] 
      #current_matrix[i,'Hayling Category A Errors z-score'] <- pt_data[129,pt_ID] 
      current_matrix[i,'Hayling Category B Errors raw'] <- pt_data[130,pt_ID] 
      #current_matrix[i,'Hayling Category B Errors z-score'] <- pt_data[131,pt_ID] 
      current_matrix[i,'BD raw'] <- pt_data[134,pt_ID] 
      current_matrix[i,'BD z-score'] <- pt_data[135,pt_ID] 
      current_matrix[i,'BD scaled score'] <- pt_data[136,pt_ID] 
      current_matrix[i,'BD No Time Bonus raw'] <- pt_data[137,pt_ID] 
      current_matrix[i,'BD No Time Bonus z-score'] <- pt_data[138,pt_ID] 
      current_matrix[i,'BD No Time Bonus scaled score'] <- pt_data[139,pt_ID] 
      current_matrix[i,'Matrix raw'] <- pt_data[140,pt_ID] 
      current_matrix[i,'Matrix z-score'] <- pt_data[141,pt_ID] 
      current_matrix[i,'Matrix scaled score'] <- pt_data[142,pt_ID] 
      current_matrix[i,'Similaritites raw'] <- pt_data[143,pt_ID] 
      current_matrix[i,'Similaritites z-score'] <- pt_data[144,pt_ID] 
      current_matrix[i,'Similaritites scaled score'] <- pt_data[145,pt_ID]
    }
    setwd('..')
  }
  #create new matrix for new timepoint
  current_matrix_name <- paste('neuropsych_matrix', timepoints[j], sep='')
  assign(current_matrix_name, current_matrix)
  
}
setwd('..')


#combine all data into one matrix
neuropsych_matrixAll <- rbind(neuropsych_matrixF0, neuropsych_matrixF1, neuropsych_matrixF2, neuropsych_matrixF3, neuropsych_matrixF4, neuropsych_matrixF5, neuropsych_matrixF6, neuropsych_matrixF8, neuropsych_matrixF10)

#write data to the 'master excel sheet' containing all data
write.xlsx(neuropsych_matrixF0, 'DPRC_neuropsych_data.xlsx', sheetName='F0', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF1, 'DPRC_neuropsych_data.xlsx', sheetName='F1', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF2, 'DPRC_neuropsych_data.xlsx', sheetName='F2', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF3, 'DPRC_neuropsych_data.xlsx', sheetName='F3', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF4, 'DPRC_neuropsych_data.xlsx', sheetName='F4', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF5, 'DPRC_neuropsych_data.xlsx', sheetName='F5', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF6, 'DPRC_neuropsych_data.xlsx', sheetName='F6', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF8, 'DPRC_neuropsych_data.xlsx', sheetName='F8', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixF10, 'DPRC_neuropsych_data.xlsx', sheetName='F10', row.names=FALSE, append=TRUE)
write.xlsx(neuropsych_matrixAll, 'DPRC_neuropsych_data.xlsx', sheetName='All', row.names=FALSE, append=TRUE)