#Visualise eddy quality control (qc) values from DPRC study

#We want to check for any outliers or abnormal patterns among our study participants with the specific 
#qc measures from the eddy qc output. These specific measurements are: movement across all volumes 
#(eddyqc_movement_all_vols), average movement (eddyqc_movement_average), percentage of outliers 
#(eddyqc_outliers.txt) and signal-to-noise (SNR) and contrast-to-noise (CNR) ratios (eddyqc_SNR&CNR.txt). 

#Author: Lenore Tahara-Eckl
#Email: Ltah262@aucklanduni.ac.nz
#Date: 24/08/20


#use pacman to install & load packages
#install.packages("pacman")


#Load libraries via pacman
pacman::p_load(readxl, ggplot2, dplyr, tidyr, plyr, tibble)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#load in pathway for functions
#source('path/tofile/here.R')
source('H:/ltah262/PhD/Diffusion/script/dprc/diffusion/visualisation/is_outlier.R')

#define pathway of where you excel sheet is
xl_data <- "H:/ltah262/PhD/Diffusion/data/eddy/eddy_qc_movement2.xlsx"


#Read in the excel file which contains the eddy qc values generated from the preprocessing pipeline. 
movement_all_vols <- read_excel(path = xl_data, sheet= 1)
movement_avg <- read_excel(path = xl_data, sheet= 2)
outliers <- read_excel(path = xl_data, sheet= 3)
SNR_CNR <- read_excel(path = xl_data, sheet= 4)



#graph your data in ggplot:

#1. for movement_all_vols (from vol 1 - absolute motion) (line plot)
ggplot(data = movement_all_vols, aes(x = Volume, y = mvmt_vol_1, colour = Participant, label = Participant)) +
    geom_line(aes(group = Participant), alpha = 0.5) +
    geom_text(size = 3, vjust = 0, nudge_x = 0.225, nudge_y = 0.225, check_overlap = TRUE) + #display the participants with this line.
    theme(legend.position = "none") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#2. for movement_all_vols (from previous vol - relative motion) (line plot)
ggplot(data = movement_all_vols, aes(x = Volume, y = mvmt_prev_vol, colour = Participant, label = Participant)) +
    geom_line(aes(group = Participant), alpha = 0.5) +
    geom_text(size = 3, vjust = 0, nudge_x = 0.225, nudge_y = 0.225, check_overlap = TRUE) + #display the participants with this line.
    theme(legend.position = "none") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#3. for movement_avg (from vol 1 - absolute motion) (violin plot)
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
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Condition, colour = Condition), size = 1) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#4. for movement_avg (from previous vol - relative motion) (violin plot)
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
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Condition, colour = Condition), size = 1) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#5. for absolute (move from vol 0) and relative motion (move from prev. vol) (raincloud plot)

#need to reshape the data to long format
names(movement_avg)[1] <- "id"
movement_avg_long <- reshape(movement_avg, varying = c("mvmt_vol_1", "mvmt_prev_vol"), v.names = "movement", timevar = "mvmt_type", times = c("mvmt_vol_1", "mvmt_prev_vol"), new.row.names = 1:1000, direction = "long")
#to plot outliers, use this  
dat <-  movement_avg_long %>% tibble::rownames_to_column(var = "outlier") %>% group_by(mvmt_type) %>% mutate(is_outlier=ifelse(is_outlier(movement), movement, as.numeric(NA)))
dat$outlier[which(is.na(dat$is_outlier))] <- as.numeric(NA)
#plot data
ggplot(data = dat, aes(y = movement, x = mvmt_type, fill = mvmt_type, label = id)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = movement, color = mvmt_type), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    #geom_text(aes(label=outlier), na.rm = TRUE, size = 3, vjust = 0, nudge_x = 0.05, check_overlap = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 1) + # add boxplot
    stat_summary(fun = mean, geom = "point", size = 2) + #add means
    guides(fill = FALSE) +
    guides(color = FALSE) +
    xlab("Movement Type") +
    ylab("Movement (mm)") +
    scale_x_discrete(labels = c("Relative Movement", "Absolute Movement")) +
    coord_flip() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#6. for percentage of outliers (violin plot)
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
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Condition, colour = Condition), size = 1) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#7. for percentage of outliers (raincloud plot)
ggplot(data = outliers, aes(y = outlier_percentage, x = Condition, label = Participant)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, aes(fill = Condition, colour = Condition)) +
    geom_point(aes(y = outlier_percentage, color = Condition), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    #geom_text(aes(label=outlier), na.rm = TRUE, size = 3, vjust = 0, nudge_x = 0.05, check_overlap = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 1, aes(colour = Condition)) + # add boxplot
    stat_summary(fun = mean, geom = "point", size = 2) + #add means
    guides(fill = FALSE) +
    guides(color = FALSE) +
    xlab("Outliers") +
    ylab("Percentage (%)") +
    #scale_x_discrete(labels = c("Relative Movement", "Absolute Movement")) +
    coord_flip() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#8. for SNR (b0) (violin plot)
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
    geom_violin(trim = FALSE, alpha = .5, aes(fill = Condition, colour = Condition), size = 1) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#9. for CNRs (b1000, b2000) (violin plot)
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
    geom_violin(trim = FALSE, alpha = .5, aes(fill = condition, colour = condition), size = 1) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#10. for SNR (b0) & CNRs (b1000, b2000) (raincloud plot)
#put data into long format
names(SNR_CNR)[1] <- "id"
SNR_CNR_long <- reshape(SNR_CNR, varying = c("SNR_b0", "CNR_b1000", "CNR_b2000"), v.names = "signal", timevar = "signal_type", times = c("SNR_b0", "CNR_b1000", "CNR_b2000"), new.row.names = 1:1000, direction = "long")
#plot data
ggplot(data = SNR_CNR_long, aes(y = signal, x = signal_type, fill = signal_type, label = id)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = signal, color = signal_type), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    #geom_text(aes(label=outlier), na.rm = TRUE, size = 3, vjust = 0, nudge_x = 0.05, check_overlap = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 1) + # add boxplot
    stat_summary(fun = mean, geom = "point", size = 2) + #add means
    guides(fill = FALSE) +
    guides(color = FALSE) +
    xlab("Signal Type") +
    ylab("Signal (SNR or CNR)") +
    scale_x_discrete(labels = c("b2000", "b1000", "b0")) +
    coord_flip() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



