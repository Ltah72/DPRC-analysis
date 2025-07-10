#PLS output measures between SLF FD and inhibition measures. SLF FD has been covaried for sex and age, and inhibition 
#measures have been covaried for age. 

#barplot (behavioral covaried for age; SLF covaried for age & sex)
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.16229625,
                                 0.15414658,
                                 0.093738705,
                                 0.10999521,
                                 0.011870223)*-1 
ulimit_PLS_FD_barplot <- c(0.322998702526093,
                           0.258648023009300,
                           0.200298555195332,
                           0.236799158155918,
                           0.128134340047836)*-1
llimit_PLS_FD_barplot <- c(0.0373809263110161,
                           0.0702661052346230,
                           0.0267098769545555,
                           0.0179583150893450,
                           -0.0657086633145809)*-1
#old
# PLS_FD_bootstrap_corr_values<- c(0.16234683,
#                                  0.15593764,
#                                  0.095833875,
#                                  0.11250146,
#                                  0.013961068)*-1 
# ulimit_PLS_FD_barplot <- c(0.323669552803040,
#                            0.261669516563416,
#                            0.196356505155563,
#                            0.236582353711128,
#                            0.126927271485329)*-1
# llimit_PLS_FD_barplot <- c(0.0347465574741364,
#                            0.0718666836619377,
#                            0.0243388582020998,
#                            0.0174392154440284,
#                            -0.0636002123355866)*-1

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
PLS_output_FD_inhib5_data <-c(-0.042637907,	-0.16580543,	-0.12276131,	-0.095673196,	  -0.13533525,	-0.15349399,
                              -0.050800979,	-0.14589985,	-0.15577170,	-0.074549116,	  -0.11840498,	-0.14621796,
                              -0.0027877691,-0.14222196,	-0.082403108,	-0.00064342923,	-0.067541681,	-0.010230582,
                              -0.031099768,	-0.12820543,	-0.039936353,	0.017875250,	  -0.12833764,	-0.089232132,
                              -0.046119776,	-0.063668542,	-0.030657532,	0.029734945,	  0.043695252,	0.072532997)

#old
# PLS_output_FD_inhib5_data <- c(-0.044300914,	-0.16606779,	-0.12265084,	-0.096109018,	-0.13687924,	-0.15432121,
#                                -0.054412790,	-0.14751273,	-0.15650904,	-0.076666139,	-0.12212638,	-0.14817488,
#                                -0.0070511983,	-0.14434420,	-0.083447576,	-0.0034717475,	-0.072038278,	-0.012402611,
#                                -0.034953527,	-0.13032354,	-0.041190058,	0.014895153,	-0.13241680,	-0.091274053,
#                                -0.050062761,	-0.065170266,	-0.031090828,	0.027725603,	0.039880473,	0.070810914)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, cl.cex = 1, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)


#average correlation (effect size) for significant SLF tracts:
#for left SLF 2 (**)
leftSLF2_avg_corr <- mean(c(-0.16606779,-0.14751273,-0.14434420,-0.13032354,-0.065170266))
leftSLF2_avg_corr 

#for left SLF 3 (*)
leftSLF3_avg_corr <- mean(c(-0.12265084,-0.15650904,-0.083447576,-0.041190058,-0.031090828))
leftSLF3_avg_corr


#barplot (behavioral z-score (one value); SLF covaried for age & sex) - FD
neuropsych_test_names_PLS <- ('Inhibition')
PLS_FD_bootstrap_corr_values<- (0.17696142)*-1 
ulimit_PLS_FD_barplot <- (0.285595774650574)*-1
llimit_PLS_FD_barplot <- c(0.126650236546993)*-1

significance_legend<-('Reliably Contributes to Latent Variable')
df_PLS_FD_barplot <- data.frame(neuropsych_test_names_PLS,significance_legend,PLS_FD_bootstrap_corr_values,ulimit_PLS_FD_barplot,llimit_PLS_FD_barplot)
#convert to factor variables
df_PLS_FD_barplot$neuropsych_test_names_PLS <- factor(df_PLS_FD_barplot$neuropsych_test_names_PLS)
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
PLS_output_FD_inhib5_data <-c(0.056728754,	0.21190090,	0.14251685,	0.039392568,	0.14188629,	0.10806088)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=1)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <- ('Inhibition')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, cl.cex = 1, col=colorRampPalette(c("white","red"))(200),is.corr=FALSE)


#barplot (behavioral z-score (one value); SLF covaried for age & sex) - FDC
neuropsych_test_names_PLS <- ('Inhibition')
PLS_FD_bootstrap_corr_values<- (0.14231814)*-1 
ulimit_PLS_FD_barplot <- (0.242426343262196)*-1
llimit_PLS_FD_barplot <- c(0.0806632302701473)*-1

significance_legend<-('Reliably Contributes to Latent Variable')
df_PLS_FD_barplot <- data.frame(neuropsych_test_names_PLS,significance_legend,PLS_FD_bootstrap_corr_values,ulimit_PLS_FD_barplot,llimit_PLS_FD_barplot)
#convert to factor variables
df_PLS_FD_barplot$neuropsych_test_names_PLS <- factor(df_PLS_FD_barplot$neuropsych_test_names_PLS)
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
PLS_output_FD_inhib5_data <-c(0.10208745,	0.17762658,	0.096085623,	0.10457813,	0.13942733,	0.082684681)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=1)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <- ('Inhibition')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, cl.cex = 1, col=colorRampPalette(c("white","red"))(200),is.corr=FALSE)



#barplot (behavioral covaried for age; SLF not covaried)
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.15886237,
                                 0.16486149,
                                 0.13376220,
                                 0.13644467,
                                 0.022582190)*-1 
ulimit_PLS_FD_barplot <- c(0.291011184453964,
                           0.278761029243469,
                           0.225330673158169,
                           0.255926370620728,
                           0.124554652720690)*-1
llimit_PLS_FD_barplot <- c(0.0558530669659376,
                           0.0723173096776009,
                           0.0593150034546852,
                           0.0341500919312239,
                           -0.0601383447647095)*-1

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
PLS_output_FD_inhib5_data <- c(-0.045279719,	-0.16150342,	-0.12181339,	-0.096407726,	-0.13711637,	-0.15515956,
                               -0.056973971,	-0.15664959,	-0.16471541,	-0.087941214,	-0.13418181,	-0.15835433,
                               -0.013154528,	-0.18154615,	-0.11966192,	-0.039840333,	-0.11224336,	-0.048079960,
                               -0.039271459,	-0.15517741,	-0.068745255,	-0.012835450,	-0.15710880,	-0.11401216,
                               -0.051900920,	-0.073506840,	-0.041012574,	0.015953202,	0.026231347,	0.058858778)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)





#barplot (behavioural z-scores; SLF not covaried - FDC)
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.10962801,
                                 0.16415893,
                                 0.14063519,
                                 0.15420002,
                                 0.018525185)*-1 
ulimit_PLS_FD_barplot <- c(0.253190353512764,
                           0.256629750132561,
                           0.244233921170235,
                           0.261911794543266,
                           0.122737228870392)*-1
llimit_PLS_FD_barplot <- c(-0.00137866870500147,
                           0.0792596749961376,
                           0.0495626311749220,
                           0.0620858334004879,
                           -0.0714475810527802)*-1

significance_legend<-c('Does Not Reliably Contribute to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable')
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
#for inhibition 5 FDC:
PLS_output_FD_inhib5_data <- c(0.064852960,	0.11341975,	0.080024593,0.092696965,0.11264170,	0.11087885,
                               0.10454190,	0.18770264,	0.14087959,	0.14718716,	0.13023023,	0.14519972,
                               0.059960660,	0.19022231,	0.11004117,	0.10936873,	0.14592695,	0.070020035,
                               0.11747230,	0.17367887,	0.10094100,	0.12497623,	0.16797434,	0.098836198,
                               0.022131911,	0.058294345,0.029687544,0.033476204,-0.020268029,	-0.061509814)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)

#barplot (behavioural z-scores; SLF not covaried - FD)
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.17899950,
                                 0.20345035,
                                 0.20374426,
                                 0.17642093,
                                 0.069172703)*-1 
ulimit_PLS_FD_barplot <- c(0.314796999096871,
                           0.311231881380081,
                           0.297280699014664,
                           0.289171978831291,
                           0.171910740435123)*-1
llimit_PLS_FD_barplot <- c(0.0763193368911743,
                           0.112822651863098,
                           0.120641369372606,
                           0.0764428526163101,
                           -0.02226649690419445)*-1

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
PLS_output_FD_inhib5_data <- c(0.046035424,	0.18345921,	0.14521618,	0.12139831,	0.14940077,	0.16308190,
                               0.057158880,	0.20064151,	0.20741850,	0.13992518,	0.15571930,	0.16800283,
                               0.032022443,	0.25123808,	0.18563139,	0.085934617,0.18105575,	0.098016694,
                               0.041379441,	0.20365563,	0.12091937,	0.069599882,0.17998947,	0.13123187,
                               0.053439595,	0.12372985,	0.092991687,0.040500723,0.0021138890,	-0.036282860)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)



#barplot (behavioural z-scores; SLF covaried for age and sex - FD)
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.16033821,
                                 0.14836954,
                                 0.13596185,
                                 0.10663046,
                                 0.010477081)*-1 
ulimit_PLS_FD_barplot <- c(0.325984403491020,
                           0.255182549357414,
                           0.240801744163036,
                           0.228078797459602,
                           0.132367067039013)*-1
llimit_PLS_FD_barplot <- c(0.0348211675882340,
                           0.0678415782749653,
                           0.0603142268955708,
                           0.0118325101211667,
                           -0.0684439241886139)*-1
#old 
# neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
# PLS_FD_bootstrap_corr_values<- c(0.16046065,
#                                  0.15053931,
#                                  0.13510700,
#                                  0.10968954,
#                                  0.013025515)*-1 
# ulimit_PLS_FD_barplot <- c(0.326225951313973,
#                            0.259664237499237,
#                            0.239632077515125,
#                            0.231231451034546,
#                            0.131857842206955)*-1
# llimit_PLS_FD_barplot <- c(0.0331384390592575,
#                            0.0710192359983921,
#                            0.0610979888588190,
#                            0.0150311659090221,
#                            -0.0670690909028053)*-1

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
PLS_output_FD_inhib5_data <- c(0.042084035,	0.16358517,	0.12158168,	0.094686612,	0.13451275,	0.15239871,
                               0.048555989,	0.14196922,	0.14934002,	0.074866898,	0.11427382,	0.13819958,
                               0.022056503,	0.17871775,	0.11385572,	0.0095995478,	0.11998393,	0.047392085,
                               0.030547176,	0.12445536,	0.038804423,-0.017527144,	0.12337276,	0.087547928,
                               0.045066155,	0.062330406,0.029861860,-0.029197706,	-0.043008354,-0.068908811)
#old
# PLS_output_FD_inhib5_data <- c(0.043890029,	0.16370192,	0.12157211,	0.095043927,	0.13602009,	0.15321481,
#                                0.052335888,	0.14367492,	0.15072666,	0.077326253,  0.11795849, 0.14012694,
#                                0.024066202,	0.17797065,	0.11291276,	0.0091826478,	0.12140537,	0.048193432,
#                                0.034642413,	0.12693506,	0.041143700,-0.013883173,	0.12754828,	0.089740954,
#                                0.049181901,	0.063852258,0.031049276,-0.026813276,	-0.038979039,-0.066828750)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, cl.cex = 1, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)




#barplot (behavioural z-scores; SLF covaried for age and sex - FDC)
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.093294308,
                                 0.12997514,
                                 0.13337146,
                                 0.12381479,
                                 -0.025506983)*-1 
ulimit_PLS_FD_barplot <- c(0.240156590938568,
                           0.224212989211083,
                           0.240339815616608,
                           0.236042931675911,
                           0.0951292105019093)*-1
llimit_PLS_FD_barplot <- c(-0.0232583703473210,
                           0.0459701474756002,
                           0.0482189077883959,
                           0.0330338496714830,
                           -0.105526130646467)*-1

significance_legend<-c('Does Not Reliably Contribute to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable')
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
PLS_output_FD_inhib5_data <- c(0.057139907,	0.094673462,	0.060951993,	0.071233489,0.10173357,	0.099570341,
                               0.090797715,	0.14964050,	0.10287382,	  0.10383445,	 0.10808627,	0.12244032,
                               0.077993840,	0.17248057,	0.094847724,	0.096409470,	0.14169085,	0.066316210,
                               0.11138196,	0.13374364,	0.061539553,	0.081456915,	0.14736915,	0.077049054,
                               -0.0019120640,	0.012367142,-0.016840154,	-0.022175740,-0.048378058,	-0.091995649)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)





#for PLS analysis w/o AD group (C, SCD, aMCI, and mMCI). 
  #-behav covaried for age
  #-SLF data covaried for age and sex
#barplot
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.148880079882814,
                                 0.137124038738739,
                                 0.101550018677087,
                                 0.149035431390626,
                                 0.00648994844356021)*-1 
ulimit_PLS_FD_barplot <- c(0.332177586562676,
                           0.255169763854637,
                           0.229473085008027,
                           0.285512006535744,
                           0.151981813940229)*-1
llimit_PLS_FD_barplot <- c(0.00554867445998682,
                           0.0578572121438813,
                           0.0337071125146369,
                           0.0412993773172031,
                           -0.0691796851668426)*-1
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
PLS_output_FD_inhib5_data <-c(-0.0196642262335174,-0.152088060635176,	-0.128926626571659,	-0.0620909220847840,	-0.119166676439566,	-0.128401628141687,
                              -0.0114568977693691,-0.134983679784785,	-0.140790936557189,	-0.00787941859339306,	-0.0975260318307024,-0.105319000874097,
                              -0.00803788882983792,	-0.152222095722902,	-0.0757395075342480,	0.0454082352257158,	-0.0605951470576168,0.0110023751879161,
                              -0.0369711825104523, -0.146867397542555,	-0.108970858994496,	-0.00314833855235335,	-0.143965023688079,	-0.0980391050243214,
                              -0.0378784214190710,-0.00534288915545500,	-0.0330311766052289,	0.117147907660564,	0.0315519659235159,	0.0290982862760147)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, cl.cex = 2, cl.align.text = "l", col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)


#-z-scores behav data
#-SLF data covaried for age and sex
#barplot
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_FD_bootstrap_corr_values<- c(0.146748041120187,
                                 0.128040063088731,
                                 0.148476519732522,
                                 0.140760810216962,
                                 0.00407689008853109)*-1 
ulimit_PLS_FD_barplot <- c(0.329634264035936,
                           0.245677634519249,
                           0.264859775689149,
                           0.277493061806009,
                           0.152381044180645)*-1
llimit_PLS_FD_barplot <- c(0.00541678749029709,
                           0.0528104097887443,
                           0.0703059942538717,
                           0.0384947899035643,
                           -0.0744643838730628)*-1
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
PLS_output_FD_inhib5_data <-c(0.0195670999776906,	0.151336861797753,	0.128289827525179,	0.0617842403617754,	  0.118578084090981,	0.127767422186398,
                              0.00917705939517311,0.129302450364990,	0.132932444704687,	0.00833198358593093,	0.0937874420497071,	0.0985221681346647,
                              0.0314372831289750,	0.191565333934513,	0.116557906183891,	-0.0233135381994107,	0.118432538153667,	0.0308403736635725,
                              0.0359598711891473,	0.141154838302969,	0.104243955033958,	0.00265411691803735,	0.136427273660256,	0.0953374825707689,
                              0.0376337867959614,	0.00597564797081891,0.0326321455806566,	-0.115794488449617,	 -0.0319331931983686,	-0.0273344958530731)
#create correlation matrix
PLS_FD_inhib5_corr_matrix <- matrix(PLS_output_FD_inhib5_data,nrow=6,ncol=5)
#add in row and column names to the matrix
colnames(PLS_FD_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
rownames(PLS_FD_inhib5_corr_matrix) <-c('Left SLF1','Left SLF2','Left SLF3','Right SLF1','Right SLF2','Right SLF3')
#plot correlation heatmap
corrplot(PLS_FD_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, cl.cex = 2, cl.align.text = "l", col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)











