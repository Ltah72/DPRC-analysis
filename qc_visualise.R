#Visualise eddy quality control (qc) values from DPRC study

#We want to check for any outliers or abnormal patterns among our study participants with the specific 
#qc measures from the eddy qc output. These specific measurements are: movement across all volumes 
#(eddyqc_movement_all_vols), average movement (eddyqc_movement_average), percentage of outliers 
#(eddyqc_outliers.txt) and signal-to-noise (SNR) and contrast-to-noise (CNR) ratios (eddyqc_SNR&CNR.txt). 

#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 24/08/20


#open the necessary libraries to do the analysis - you may need to install the packages if you don't have them. 
#install.packages("readxl")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("tidyr")

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

#define pathway of where you excel sheet is
xl_data <- "H:/ltah262/PhD/Diffusion/data/eddy/eddy_qc_movement.xlsx"


#Read in the excel file which contains the eddy qc values generated from the preprocessing pipeline. 
movement_all_vols <- read_excel(path = xl_data, sheet= 1)
movement_avg <- read_excel(path = xl_data, sheet= 2)
outliers <- read_excel(path = xl_data, sheet= 3)
SNR_CNR <- read_excel(path = xl_data, sheet= 4)



#graph your data in ggplot:


#1. for movement_all_vols (from vol 1 - absolute motion)
ggplot(data = movement_all_vols, aes(x = Volume, y = mvmt_vol_1, colour = Participant, label = Participant)) +
    geom_line(aes(group = Participant), alpha = 0.5) +
    geom_text(size = 3, vjust = 0, nudge_x = 0.225, nudge_y = 0.225, check_overlap = TRUE) + #display the participants with this line.
    theme(legend.position = "none")
   

#2. for movement_all_vols (from previous vol - relative motion)
ggplot(data = movement_all_vols, aes(x = Volume, y = mvmt_prev_vol, colour = Participant, label = Participant)) +
    geom_line(aes(group = Participant), alpha = 0.5) +
    geom_text(size = 3, vjust = 0, nudge_x = 0.225, nudge_y = 0.225, check_overlap = TRUE) + #display the participants with this line.
    theme(legend.position = "none")


#3. for movement_avg (from vol 1 - absolute motion)
ggplot(data = movement_avg, aes(x = Condition, y = mvmt_vol_1, label = Participant)) +
    geom_point() +
    geom_text(size = 3, vjust = 0, nudge_x = 0.225, nudge_y = 0.225, check_overlap = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = 1)) + # add boxplot
    stat_summary(fun = mean, geom = "point", size = 2, aes(colour = 1))+ # add means
    ggtitle("Abs.motion") +
    theme(plot.title = element_text(hjust = 0.5)) + # centre the plot title
    ylab("mm(avg)") +
    scale_x_discrete(labels = c("0" = "Participant")) +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Condition, colour = Condition), size = 1)

#4. for movement_avg (from previous vol - relative motion)
ggplot(data = movement_avg, aes(x = Condition, y = mvmt_prev_vol, label = Participant)) +
    geom_point() +
    geom_text(size = 3, vjust = 0, nudge_x = 0.225, nudge_y = 0.225, check_overlap = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = 1)) + # add boxplot
    stat_summary(fun = mean, geom = "point", size = 2, aes(colour = 1))+ # add means
    ggtitle("Rel.motion") +
    theme(plot.title = element_text(hjust = 0.5)) + # centre the plot title
    ylab("mm(avg)") +
    scale_x_discrete(labels = c("0" = "Participant")) +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Condition, colour = Condition), size = 1)


#5. for percentage of outliers
ggplot(data = outliers, aes(x = Condition, y = outlier_percentage, label = Participant)) +
    geom_point() +
    geom_text(size = 3, vjust = 0, nudge_x = 0.225, nudge_y = 0.225, check_overlap = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = 1)) + # add boxplot
    stat_summary(fun = mean, geom = "point", size = 2, aes(colour = 1))+ # add means
    ggtitle("Total Outliers") +
    theme(plot.title = element_text(hjust = 0.5)) + # centre the plot title
    ylab("%") +
    scale_x_discrete(labels = c("0" = "Participant")) +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Condition, colour = Condition), size = 1)


#6. for SNR (b0)
ggplot(data = SNR_CNR, aes(x = Condition, y = SNR_b0, label = Participant)) +
    geom_point() +
    geom_text(size = 3, vjust = 0, nudge_x = 0.225, nudge_y = 0.225, check_overlap = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = 1)) + # add boxplot
    stat_summary(fun = mean, geom = "point", size = 2, aes(colour = 1))+ # add means
    ggtitle("SNR (avg)") +
    theme(plot.title = element_text(hjust = 0.5)) + # centre the plot title
    xlab("b-value (0)") +
    ylab("signal") +
    scale_x_discrete(labels = c("0" = "Participant")) +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Condition, colour = Condition), size = 1)


#7. for CNRs (b1000, b2000)
#rearrange data for this plot. Extract interested variables. 
CNR_data_only <- select(SNR_CNR, Participant, CNR_b1000, CNR_b2000)
#change dataset from wide format to long format
CNR_long <- gather(CNR_data_only, "condition", "signal", CNR_b1000, CNR_b2000)
#plot the graph
ggplot(data = CNR_long, aes(x = condition, y = signal, label = Participant)) +
    geom_point() +
    geom_text(size = 3, vjust = 0, nudge_x = 0.225, nudge_y = 0.225, check_overlap = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = condition)) + # add boxplot
    stat_summary(fun = mean, geom = "point", size = 2, aes(colour = condition)) + # add means
    ggtitle("CNR (avg)") +
    theme(plot.title = element_text(hjust = 0.5)) + # centre the plot title
    xlab("bvalue") +
    ylab("signal") +
    scale_x_discrete(labels = c("CNR_b1000" = "b1000", "CNR_b2000" = "b2000")) +
    theme(legend.position = "none") +
    geom_violin(trim = FALSE, alpha = .5, aes(fill = condition, colour = condition), size = 1)

