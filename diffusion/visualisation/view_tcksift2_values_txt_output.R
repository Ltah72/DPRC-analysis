#Plot a histogram to view the distribution of the streamline weight values from
#the output text files of MRtrix3's command tcksift2. This is from a single 
#participant (ADPRC0001F0). Each txt file contains 10 millions values. 



#for one participant - e.g., ADPRC0001F0
#define and set up file path and variables
file_path <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/sift_sub-ADPRC0004F2.txt'
sift_text_file <- scan(file_path, what = 'numeric')

#Remove MRtrix's command line from file - make sure to check the numbers for this, as it can vary.
sift_text_file <- sift_text_file[15:10000014]
sift_text_file <- sift_text_file[14:10000013]

#put into dataframe
dataframe_sift <- read.csv(textConnection(sift_text_file),header=FALSE)

#check values--
#range (min - 0.1147106, max- 8.5923341)
range(dataframe_sift[['V1']])

#mean (mean - 1.126185)
mean(dataframe_sift[['V1']])

#create histogram with 91 bins (from 0 to 9; 10 values per each integer [e.g., 0.1, 0.2, 0.3...] and n+1)
hist(dataframe_sift[['V1']], breaks=seq(min(dataframe_sift[['V1']]), max(dataframe_sift[['V1']]), length.out=91), labels=TRUE, las = 0, xlab="Streamline weight value", main="Histogram of P1s streamline weight values")
axis(side=1, at=seq(0,9,0.1), labels=seq(0,9,0.1))
text(x = 1:length(dataframe_sift[['V1']]), srt=35)






#for the average
#define and set up directory to take average of the text files
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/tcksift2_text_files'
setwd(directory)

txt_files <- lapply(list.files(file_path), read.table, header=FALSE)
all_avg_txts<-Reduce("+", txt_files) / length(txt_files)
write.table(all_avg_txts, 'all_avg_sift.txt', sep = ",", row.names= FALSE, col.names = FALSE)

#set up variables
file_path <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/all_avg_sift.txt'
sift_text_file <- scan(file_path, what = 'numeric')

#Remove MRtrix's command line from file
sift_text_file <- sift_text_file[15:10000014]

#put into dataframe
dataframe_sift <- read.csv(textConnection(sift_text_file),header=FALSE)

#check values--
#range (min - 0.1147106, max- 8.5923341)
range(dataframe_sift[['V1']])
#mean (mean - 1.126185)
mean(dataframe_sift[['V1']])

#create histogram with 91 bins (from 0 to 9; 10 values per each integer [e.g., 0.1, 0.2, 0.3...] and n+1)
hist(dataframe_sift[['V1']], breaks=seq(min(dataframe_sift[['V1']]), max(dataframe_sift[['V1']]), length.out=91), labels=TRUE, las = 0, xlab="Streamline weight value", main="Histogram of P1s streamline weight values")
axis(side=1, at=seq(0,9,0.1), labels=seq(0,9,0.1))
text(x = 1:length(dataframe_sift[['V1']]), srt=35)







