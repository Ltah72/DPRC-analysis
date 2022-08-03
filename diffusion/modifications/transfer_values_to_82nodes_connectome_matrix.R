

#set up packages:
pacman::p_load(ggplot2, corrplot, viridis)


# choose & set to your directory. This is where each of your participant's 
#values for the interested connectome network files should be. 
#setwd("/Connectome_test/FPN_connectome_csv_files/done/output")  

setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/stats_results/weighted/FPN_stats/bigNode_linear_trend')

#read in the family-wise error (fwe) stat files
FPN_stats_connectome_fwe0 <- read.csv('transferred_values_big-node_linear_trend_FPN_stats-Zstat_t1.csv', sep = "", header = FALSE)

#remove the rows and columns if they contain all zeros in them.
#z<-y[rowSums(y[])>0,colSums(y[])>0]

#use the same node pair as in the previous script (create_FPN_big_node_connectome_files.R) 
# - my chosen nodes are below:
#nodes <- c(" ")
nodes <- c(73,253,67,247,97,277,98,278,26,206,70,250,71,251,87,267,68,248,83,263,85,265,84,264,86,266,40,220,41,221,55,235,44,224,43,223,36,216,39,219,37,217,48,228,95,275,49,229,117,297,50,230,47,227,42,222,45,225,46,226,29,209,143,323,151,331,150,330,149,329,148,328,116,296,147,327,146,326,145,325,144,324)
nodes_sorted <- sort(nodes)

#remove all rows and columns except the assigned nodes of interest (82 nodes)
#if you want to view the nodes by network (good to see how networks relate to one another)
FPN_stats_82nodes_connectome_by_network_fwe0 <- FPN_stats_connectome_fwe0[nodes,nodes]
#if you want to view the nodes in order (good to assess centrality & laterality)
FPN_stats_82nodes_connectome_node_order_fwe0 <- FPN_stats_connectome_fwe0[nodes_sorted,nodes_sorted]

#view matrices with correlation matrix plot:
corrplot(data.matrix(FPN_stats_82nodes_connectome_by_network_fwe0),is.corr=FALSE,method = "color",col=viridis(200), tl.col = "black")
#cl.lim=c(.95,1)
#is.corr=FALSE

#if you want to display only sig. values (>.95), delete values which are lower than .95 and replace them with zeros (0) in the matrix.
FPN_stats_sig_82nodes_connectome_by_network<-replace(data.matrix(FPN_stats_82nodes_connectome_by_network_fwe0), data.matrix(FPN_stats_82nodes_connectome_by_network_fwe0)<.95, 0) 
#display sig. nodes only (>.95)
corrplot(FPN_stats_sig_82nodes_connectome_by_network,method = "color", tl.col = "black")





#write matrices as a .csv file. 
write.table(format(FPN_stats_82nodes_connectome_fwe0, digits=20), 'transferred_values_big-node_linear_trend_82_nodes_FPN_stats-fwe_1mpvalue_t1.csv', row.names=FALSE, col.names=FALSE, quote=FALSE)

