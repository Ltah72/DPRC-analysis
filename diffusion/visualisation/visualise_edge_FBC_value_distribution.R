#compare edges for a given node -- e.g, 
#compare the edges between nodes 27 (left superior frontal) and 30 (left 
#supramarginal gyrus)

#load packages
pacman::p_load(ggplot2, dplyr, plyr)

#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph

#set directory
directory <- 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/dkFiles/weighted/fs_default_ordered/thresholded_connectomes/C&SCD'
setwd(directory)

#read in all participant .csv files
all_csvs <- lapply(list.files(directory), read.csv, header=FALSE)

#extract edge from node 27 and 30 for each participant and put into a dataframe
edges_value_vec <- unlist(lapply(all_csvs, `[`, 36, 56))
#This elegant line of code comes from: https://stackoverflow.com/questions/23758858/how-can-i-extract-elements-from-lists-of-lists-in-r

#create other columns of the dataframehttp://127.0.0.1:21613/graphics/plot_zoom_png?width=1051&height=900
participant_vec <- substr(list.files(directory),4,18)
group_numbers_vec <- c(2,	5,	3,	5,	3,	1,	2,	4,	2,	2,	2,	4,	4,	1,	1,	1,	3,	2,	3,	1,	2,	3,	3,	3,	1,	4,	5,	3,	3,	3,	4,	4,	3,	2,	3,	4,	5,	3,	1,	3,	3,	1,	1,	1,	5,	4,	4,	5,	2,	2,	4,	2,	3,	1,	1,	1,	2,	4,	4,	3,	2,	1,	2,	1,	2,	2,	2,	1,	5,	2,	4,	4,	2,	5,	2,	4,	3,	2,	4,	2,	2,	2,	2,	2,	3,	3,	3,	2,	2,	5,	1,	1,	4,	2,	4,	2,	3,	4,	3,	5,	5,	2,	5,	2,	2,	3,	2,	5,	4,	1,	5,	4,	4,	2,	1,	2,	5,	3,	2,	4,	2,	2,	4,	1,	1,	3,	4,	4,	4,	2,	4,	4,	4,	1,	2,	3,	2,	4,	3,	3,	3,	4,	1,	3,	3,	1,	1,	1,	4,	3,	1,	4,	2,	3,	1,	1,	2,	4,	5,	1,	3,	4,	1,	4,	5,	2,	1,	4,	5,	5,	3,	5,	3,	3,	4,	1,	4,	3,	1,	2,	2,	4,	5,	2,	4,	3,	3,	3,	2,	3,	3,	2,	3,	4,	4,	3,	4,	5,	3,	3,	2,	4,	3,	4,	4,	4,	3,	3,	5,	2,	3,	2,	4,	3,	2,	2,	2,	2,	2,	3,	4,	4,	2,	5,	2,	5,	5)

#create and combine into one dataframe
dataframe_edge_compare <- data.frame(matrix(ncol=3,nrow=227))
colnames(dataframe_edge_compare) <- c('ParticipantID','Group','FBC_edge_value')
dataframe_edge_compare$ParticipantID <- participant_vec 
dataframe_edge_compare$Group <- as.factor(group_numbers_vec) 
dataframe_edge_compare$FBC_edge_value <- edges_value_vec


#violin plot
ggplot(dataframe_edge_compare, aes(x = Group, y = FBC_edge_value, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = FBC_edge_value, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("FBC Edge Value") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=18))+
  theme(axis.title=element_text(size=24))#+
  #coord_flip()

#find max FBC edge value by group
result_max <- dataframe_edge_compare %>% 
  group_by(Group) %>%
  filter(FBC_edge_value == max(FBC_edge_value))

# #find min FBC edge value by group
#   result_min <- dataframe_edge_compare %>% 
#   group_by(Group) %>%
#   filter(FBC_edge_value == min(FBC_edge_value))

#find min FBC edge value by group (non-zero)
  result_min <- dataframe_edge_compare %>% 
    group_by(Group) %>%
    filter(FBC_edge_value == min(FBC_edge_value[FBC_edge_value>0]))




