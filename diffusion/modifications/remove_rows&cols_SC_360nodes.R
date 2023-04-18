#Remove 19 rows and columns from each participant's structural connectome. Transform the connectome from 379x379
#to 360x360.  

library(tidyverse)

#set up pathway (where you hcpmmp files are):
setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/all_files_360nodes')


files <- list.files(".", pattern = ".csv")

for (i in seq(files)) {
  PAR_NAME <- substr(files[i],10,24)
  
  dat <- files[i] %>% 
    map_dfr(
      ~ read_csv(.x, col_names=FALSE) %>% 
        slice(1:360) #only include the first 360 rows
    ) %>% 
    select(1:360) #only include the first 360 columns
  
  #write to a new .csv file
  write.table(dat, paste('hcpmmp1w_', PAR_NAME, '_360nodes.csv', sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE,sep=",")
}

