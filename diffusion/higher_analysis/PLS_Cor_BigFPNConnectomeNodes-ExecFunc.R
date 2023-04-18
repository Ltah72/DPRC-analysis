#PLS output measures between the big FPN nodes (3,321) and executive function (inhibition and proc speed) measures. 


pacman::p_load(ggplot2, corrplot, viridis, magrittr)


#Inhibition: 

#barplot
neuropsych_test_names_PLS <- c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
PLS_bootstrap_corr_values<- c(0.18347713,
                              0.31900072,
                              0.31740135,
                              0.16431601,
                              0.24171278)*-1 
ulimit_PLS_barplot <- c(0.348343864083290,
                        0.474500849843025,
                        0.467083469033241,
                        0.315400630235672,
                        0.394739598035812)*-1
llimit_PLS_barplot <- c(0.141206495463848,
                        0.290422126650810,
                        0.288500636816025,
                        0.123571436852217,
                        0.216989666223526)*-1

significance_legend<-c('Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable')
df_PLS_barplot <- data.frame(neuropsych_test_names_PLS,significance_legend,PLS_bootstrap_corr_values,ulimit_PLS_barplot,llimit_PLS_barplot)
#df_PLS_barplot <- data.frame(neuropsych_test_names_PLS,PLS_bootstrap_corr_values,ulimit_PLS_barplot,llimit_PLS_barplot)

#convert to factor variables
df_PLS_barplot$neuropsych_test_names_PLS <- factor(df_PLS_barplot$neuropsych_test_names_PLS, levels=c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError'))
df_PLS_barplot$significance_legend <- factor(df_PLS_barplot$significance_legend, levels=c('Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable'))
#barplot(PLS_FD_bootstrap_corr_values, xlab = "Neuropsychological Assessment", ylab = "Bootstrap Correlation Value")

ggplot(data=df_PLS_barplot,(aes(x=neuropsych_test_names_PLS, y=PLS_bootstrap_corr_values,fill=significance_legend)))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(x=neuropsych_test_names_PLS,ymin=llimit_PLS_barplot,ymax=ulimit_PLS_barplot))+
  xlab("Neuropsychological Assessment") + 
  ylab("Correlation Values") +
  scale_fill_manual(values=c("steelblue","light grey"))+
  labs(fill=NULL)+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18,angle = 45, vjust = 1, hjust=1),axis.text.y = element_text(size = 18))


#Correlation heatmap plot:
#for inhibition 5 FBC:

#choose directory for data
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/PLS/MDI/cross-sectional/dMRI/group_contrast/correlations')
#read in data file
PLS_output_inhib5_data <- read.csv('PLS_correlation_BigNodeFPN_raw-inhibition.csv')

#create correlation matrix
#PLS_inhib5_corr_matrix <- matrix(PLS_output_inhib5_data,nrow=3321,ncol=5)
PLS_inhib5_corr_matrix <- as.matrix(PLS_output_inhib5_data[4:8])
#add in row and column names to the matrix
colnames(PLS_inhib5_corr_matrix) <-c('TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotError')
#rownames(PLS_inhib5_corr_matrix) <-c()
rownames(PLS_inhib5_corr_matrix) <- PLS_output_inhib5_data$Connection_Types



#plot correlation heatmap
corrplot(PLS_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)
corrplot(PLS_inhib5_corr_matrix, method = "color", tl.col = "black", tl.cex = 0.25, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)


#y<-as.matrix(x)
#corrplot(z)
#corrplot(z, method = "color", tl.col = "black", tl.cex = 1.5, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)

#z<-y[c(1:50),]
#average correlation (effect size) for significant SLF tracts:
#for left SLF 2 (**)
leftSLF2_avg_corr <- mean(c(-0.16606779,-0.14751273,-0.14434420,-0.13032354,-0.065170266))
leftSLF2_avg_corr 

#for left SLF 3 (*)
leftSLF3_avg_corr <- mean(c(-0.12265084,-0.15650904,-0.083447576,-0.041190058,-0.031090828))
leftSLF3_avg_corr 




#Proc Speed:

#barplot
neuropsych_test_names_PLS <- c('TMT-A','Colour Naming','Word Reading')
PLS_bootstrap_corr_values<- c(0.21657050,
                              0.30868170,
                              0.29746792)*-1 
ulimit_PLS_barplot <- c(0.387492910027504,
                        0.463861480355263,
                        0.466919660568237)*-1
llimit_PLS_barplot <- c(0.173085004091263,
                        0.312210470438004,
                        0.275467216968536)*-1

significance_legend<-c('Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable','Reliably Contributes to Latent Variable')
df_PLS_barplot <- data.frame(neuropsych_test_names_PLS,significance_legend,PLS_bootstrap_corr_values,ulimit_PLS_barplot,llimit_PLS_barplot)
#df_PLS_barplot <- data.frame(neuropsych_test_names_PLS,PLS_bootstrap_corr_values,ulimit_PLS_barplot,llimit_PLS_barplot)

#convert to factor variables
df_PLS_barplot$neuropsych_test_names_PLS <- factor(df_PLS_barplot$neuropsych_test_names_PLS, levels=c('TMT-A','Colour Naming','Word Reading'))
df_PLS_barplot$significance_legend <- factor(df_PLS_barplot$significance_legend, levels=c('Reliably Contributes to Latent Variable','Does Not Reliably Contribute to Latent Variable'))
#barplot(PLS_FD_bootstrap_corr_values, xlab = "Neuropsychological Assessment", ylab = "Bootstrap Correlation Value")

ggplot(data=df_PLS_barplot,(aes(x=neuropsych_test_names_PLS, y=PLS_bootstrap_corr_values,fill=significance_legend)))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(x=neuropsych_test_names_PLS,ymin=llimit_PLS_barplot,ymax=ulimit_PLS_barplot))+
  xlab("Neuropsychological Assessment") + 
  ylab("Correlation Values") +
  scale_fill_manual(values=c("steelblue","light grey"))+
  labs(fill=NULL)+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x = element_text(size = 18,angle = 45, vjust = 1, hjust=1),axis.text.y = element_text(size = 18))


#Correlation heatmap plot:
#for ProcSpeed 3 FBC:

#choose directory for data
setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/PLS/MDI/cross-sectional/dMRI/group_contrast/correlations')
#read in data file
PLS_output_ProcSpeed3_data <- read.csv('PLS_correlation_BigNodeFPN_raw-ProcSpeed_BSR3stars_smaller.csv')

#create correlation matrix
#PLS_inhib5_corr_matrix <- matrix(PLS_output_inhib5_data,nrow=3321,ncol=5)
PLS_ProcSpeed3_corr_matrix <- as.matrix(PLS_output_ProcSpeed3_data[6:8])
#add in row and column names to the matrix
colnames(PLS_ProcSpeed3_corr_matrix) <-c('TMT-A','Colour Naming','Word Reading')
#rownames(PLS_inhib5_corr_matrix) <-c()
rownames(PLS_ProcSpeed3_corr_matrix) <- PLS_output_ProcSpeed3_data$Connection_Types



#plot correlation heatmap
corrplot(PLS_ProcSpeed3_corr_matrix, method = "color", tl.col = "black", tl.cex = 1, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)
corrplot(PLS_ProcSpeed3_corr_matrix, method = "color", tl.col = "black", tl.cex = 0.25, tl.srt = 45, col=colorRampPalette(c("blue","white","red"))(200),is.corr=FALSE)






