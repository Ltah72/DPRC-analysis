#Show Histogram of FBC values



#for all
#set up directory
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/all_files/'
setwd(directory)
#read in data file
all_data <- read.csv("all_avg_connectome.csv", header=FALSE)
#plot histogram
#par(mfrow=c(2,1))
#histout=apply(all_data,2,hist)
vec_all_data <- unlist(all_data,use.names=FALSE)
hist(vec_all_data)
hist(vec_all_data, ylim=c(0,1000))


#for C
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/C/'
setwd(directory)
#read in data file
C_data <- read.csv("Control_avg_connectome.csv", header=FALSE)
#plot histogram
vec_C_data <- unlist(C_data,use.names=FALSE)
hist(vec_C_data)
hist(vec_C_data, ylim=c(0,1000))


#for SCD
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/SCD/'
setwd(directory)
#read in data file
SCD_data <- read.csv("SCD_avg_connectome.csv", header=FALSE)
#plot histogram
vec_SCD_data <- unlist(SCD_data,use.names=FALSE)
hist(vec_SCD_data)
hist(vec_SCD_data, ylim=c(0,1000))


#for aMCI
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/aMCI/'
setwd(directory)
#read in data file
aMCI_data <- read.csv("aMCI_avg_connectome.csv", header=FALSE)
#plot histogram
vec_aMCI_data <- unlist(aMCI_data,use.names=FALSE)
hist(vec_aMCI_data)
hist(vec_aMCI_data, ylim=c(0,1000))


#for mMCI
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/mMCI/'
setwd(directory)
#read in data file
mMCI_data <- read.csv("mMCI_avg_connectome.csv", header=FALSE)
#plot histogram
vec_mMCI_data <- unlist(mMCI_data,use.names=FALSE)
hist(vec_mMCI_data)
hist(vec_mMCI_data, ylim=c(0,1000))


#for AD
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/hcpmmpFiles/weighted/AD/'
setwd(directory)
#read in data file
AD_data <- read.csv("AD_avg_connectome.csv", header=FALSE)
#plot histogram
vec_AD_data <- unlist(AD_data,use.names=FALSE)
hist(vec_AD_data)
hist(vec_AD_data, ylim=c(0,1000))


