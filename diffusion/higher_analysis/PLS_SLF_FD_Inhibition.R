#PLS output measures between SLF FD and inhibition measures. SLF FD has been covaried for sex and age, and inhibition 
#measures have been covaried for age. 

#barplot
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.16234683,
                                 0.15593764,
                                 0.095833875,
                                 0.11250146,
                                 0.013961068)*-1 
ulimit_PLS_FD_barplot <- c(0.323669552803040,
                           0.261669516563416,
                           0.196356505155563,
                           0.236582353711128,
                           0.126927271485329)*-1
llimit_PLS_FD_barplot <- c(0.0347465574741364,
                           0.0718666836619377,
                           0.0243388582020998,
                           0.0174392154440284,
                           -0.0636002123355866)*-1

significance_legend<-c('Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable')
df_PLS_FD_barplot <- data.frame(neuropsych_test_names_PLS,significance_legend,PLS_FD_bootstrap_corr_values,ulimit_PLS_FD_barplot,llimit_PLS_FD_barplot)
#convert to factor variables
df_PLS_FD_barplot$neuropsych_test_names_PLS <- factor(df_PLS_FD_barplot$neuropsych_test_names_PLS, levels=c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError'))
df_PLS_FD_barplot$significance_legend <- factor(df_PLS_FD_barplot$significance_legend, levels=c('Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable'))
#barplot(PLS_FD_bootstrap_corr_values, xlab = "Neuropsychological Assessment", ylab = "Bootstrap Correlation Value")
ggplot(data=df_PLS_FD_barplot,(aes(x=neuropsych_test_names_PLS, y=PLS_FD_bootstrap_corr_values,fill=significance_legend)))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(x=neuropsych_test_names_PLS,ymin=llimit_PLS_FD_barplot,ymax=ulimit_PLS_FD_barplot))+
  xlab("Neuropsychological Assessment") + 
  ylab("Correlation Values") +
  scale_fill_manual(values=c("steelblue","light grey"))+
  labs(fill=NULL)+
  theme_classic()+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18,angle = 45, vjust = 1, hjust=1),axis.text.y = element_text(size = 18))
  




#Correlation heatmap plot:
#for inhibition 5 FD:
PLS_output_FD_inhib5_data <- c(-0.044300914,	-0.16606779,	-0.12265084,	-0.096109018,	-0.13687924,	-0.15432121,
                               -0.054412790,	-0.14751273,	-0.15650904,	-0.076666139,	-0.12212638,	-0.14817488,
                               -0.0070511983,	-0.14434420,	-0.083447576,	-0.0034717475,	-0.072038278,	-0.012402611,
                               -0.034953527,	-0.13032354,	-0.041190058,	0.014895153,	-0.13241680,	-0.091274053,
                               -0.050062761,	-0.065170266,	-0.031090828,	0.027725603,	0.039880473,	0.070810914)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)


#average correlation (effect size) for significant SLF tracts:
#for left SLF 2 (**)
leftSLF2_avg_corr <- mean(c(-0.16606779,-0.14751273,-0.14434420,-0.13032354,-0.065170266))
leftSLF2_avg_corr 

#for left SLF 3 (*)
leftSLF3_avg_corr <- mean(c(-0.12265084,-0.15650904,-0.083447576,-0.041190058,-0.031090828))
leftSLF3_avg_corr 




