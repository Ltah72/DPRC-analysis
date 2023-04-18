#Post hoc pairwise comparison p-value with False Discovery Rate (FDR) correction. 
#Output from the whole-brain and network of interest (frontoparietal network 
#(FPN)) functional connectivity results from CONN analysis. 
#This script will take in the p-values in vector form and calculate the FDR 
#correction. 

pacman::p_load(fuzzySim)
#webiste: https://www.rdocumentation.org/packages/fuzzySim/versions/4.9.2/topics/FDR


#example from website:
#pvalues <- data.frame(var = letters[1:5], pval = c(0.02, 0.004, 0.07, 0.03, 0.05))
#FDR(pvalues=pvalues)

#20 p-values values
# Wholebrain_pvalues <- data.frame(group_comparison = c('C>SCD','C<SCD','C>aMCI',
#                                                       'C<aMCI','C>mMCI','C<mMCI',
#                                                       'C>AD','C<AD','SCD>aMCI',
#                                                       'SCD<aMCI','SCD>mMCI',
#                                                       'SCD<mMCI','SCD>AD','SCD<AD',
#                                                       'aMCI>mMCI','aMCI<mMCI',
#                                                       'aMCI>AD','aMCI<AD','mMCI>AD',
#                                                       'mMCI<AD'), pval = c(0.02, 0.004, 0.07, 0.03, 0.05))


#For the whole-brain analysis (2 clusters) (no covariates)
Wholebrain_pvalues_cluster1 <- data.frame(group_comparison = c('C>SCD','C>aMCI',
                                                      'C>mMCI','C>AD',
                                                      'SCD>aMCI','SCD<mMCI',
                                                      'SCD>AD','aMCI<mMCI',
                                                      'aMCI>AD','mMCI>AD'), pval = c(0.007988,0.000533,0.067818,0.000066,0.837643,0.411705,0.039883,0.107733,0.103977,0.002283))

Wholebrain_pvalues_cluster2 <- data.frame(group_comparison = c('C>SCD','C>aMCI',
                                                               'C>mMCI','C>AD',
                                                               'SCD>aMCI','SCD<mMCI',
                                                               'SCD>AD','aMCI<mMCI',
                                                               'aMCI>AD','mMCI>AD'), pval = c(0.001257,0.000533,0.206694,0.000089,0.837643,0.022614,0.172624,0.025601,0.103977,0.002113))
FDR(pvalues=Wholebrain_pvalues_cluster1)
FDR(pvalues=Wholebrain_pvalues_cluster2)




# #For the frontoparietal network analysis (1 cluster) (sex & age covariates) (corrected twice)
# FPN_pvalues_covars <- data.frame(group_comparison = c('C<SCD','C>aMCI','C<mMCI','C>AD',
#                                               'SCD>aMCI','SCD<mMCI','SCD>AD',
#                                                'aMCI<mMCI','aMCI>AD','mMCI>AD'), pval = c(0.877294,0.369212,0.372981,0.485386,0.276214,0.231096,0.414047,0.030412,0.88582,0.148099))
# FDR(pvalues=FPN_pvalues_covars)


#For the frontoparietal network analysis (1 cluster) (no covariates)
FPN_pvalues <- data.frame(group_comparison = c('C<SCD','C>aMCI','C<mMCI','C>AD',
                                               'SCD>aMCI','SCD<mMCI','SCD>AD',
                                               'aMCI<mMCI','aMCI>AD','mMCI>AD'), pval = c(0.904991,0.096093,0.025783,0.114336,0.053691,0.027302,0.095954,0.000084,0.909784,0.000769))
FDR(pvalues=FPN_pvalues)



#For the frontoparietal network analysis (1 cluster) (sex & age covariates) (corrected once)
FPN_pvalues_covars <- data.frame(group_comparison = c('C<SCD','C>aMCI','C<mMCI','C>AD',
                                                      'SCD>aMCI','SCD<mMCI','SCD>AD',
                                                      'aMCI<mMCI','aMCI>AD','mMCI>AD'), pval = c(0.895513,0.105265,0.023934,0.126579,0.056133,0.024683,0.100664,0.000090,0.922691,0.000822))
FDR(pvalues=FPN_pvalues_covars)





