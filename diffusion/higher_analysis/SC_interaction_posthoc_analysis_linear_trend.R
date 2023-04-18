pacman::p_load(ggplot2, multcomp, lsr, dplyr)
#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph


#for the whole brain structural connectome - 44 significant edges for interaction
#load in data, copied from matlab
C_mean <- c(0.00575161049781501,
            -0.00796977732441966,
            0.000289573750317234,
            0.00205078152243688,
            -0.00123035034667045,
            0.00104554481199603,
            0.000359204775328804,
            -0.00226100166297372,
            -0.00447251942695996,
            0.00523671559010567,
            0.0133165766722218,
            0.0105884697434933,
            0.00492785007337450,
            -0.000380097137886851,
            0.0184972707252998,
            -4.48964961603924e-05,
            -0.000826924562141787,
            0.00501883485793570,
            0.00394591011881223,
            0.000279468167402942,
            0.0178273759347515) #n = 21

SCD_mean <- c(-0.00350955823533134,
              0.0123754181080009,
              0.00471027134096341,
              0.00185568994670126,
              -0.000583671020069227,
              8.31627802125145e-05,
              0.00256068549338065,
              0.00423054816850166,
              0.00113478617358578,
              -0.00424429463918512,
              -0.00390826534401634,
              -0.00129466809000573,
              0.00540542588740180,
              0.00924146447905399,
              -3.12948882568216e-05,
              0.00131147672407757,
              0.00426292820600626,
              0.00523857258396970,
              -0.00180628786289965,
              0.00510743823356043,
              0.00147243834575272,
              0.00228344971396824,
              0.00901428463732863,
              0.00624877069109778,
              0.00309027204285516,
              0.00721825446664238,
              0.0130452809894248,
              4.46335610312916e-05,
              0.000931454294013367,
              -0.00407731078826188,
              0.00251730693617412,
              0.00139675651446285,
              0.0134596327388468,
              -6.61441566618147e-06,
              -0.00452650504912824,
              -0.00161589092369292,
              0.00154227114005151,
              -0.000978179726377265,
              0.000926001479337587,
              0.00700402288610714) #n = 40

aMCI_mean <- c(0.00980132713143291,
               0.00341963749862397,
               0.00372520897933067,
               -0.000145808850535400,
               -0.00137980716419119,
               0.00508241010844436,
               0.0114884625097343,
               -0.00790278493071357,
               0.0105147425437085,
               0.00139791635383797,
               0.00140109192191856,
               0.00629866292214689,
               0.00238924213083292,
               -0.00120513161663469,
               0.0140118360938688,
               -0.00205111722435118,
               0.00656340332470773,
               -0.00373962593212453,
               -0.00123895954997558,
               0.00503760937830165,
               0.00171011216916547,
               0.0154809242835818,
               0.00280247875687003,
               0.0101795455802117,
               -0.00393617718851942,
               0.000730651209756147,
               0.0109203112620216,
               0.00874157473925577,
               -0.0207415686732292) #n = 29

mMCI_mean <- c(0.00606874887283650,
               0.0124365912508551,
               0.00366257992558272,
               0.00109281822758944,
               -0.000845867271652063,
               0.00306980578810869,
               0.000394703907313753,
               0.000658076840555110,
               0.00471414127290306,
               0.000336712097141123,
               -0.00378725073927189,
               -0.00305998161095170,
               0.00986868356005840,
               0.00334978816743654,
               -0.00584912870383866,
               -0.0169643574926343,
               -0.000943513138980518,
               -0.000952097548261628,
               0.00464816196375605) #n = 19

AD_mean <- c(-0.00726364409998651,
             -0.0153610078110010,
             -0.00238496769323297,
             -0.0113031640005995,
             -0.00953914083514844,
             -0.00173669289945752,
             -0.0242630558506964,
             -0.0116481100449644,
             -0.00309543154831533,
             -0.0193902389796867) #n = 10

#combine into one dataframe
all_groups_mean <-data.frame(matrix(ncol=2,nrow=119))
colnames(all_groups_mean) <- c('Group','change_in_FBC(post-pre)')
all_groups_mean$`change_in_FBC(post-pre)`<-(c(C_mean,SCD_mean,aMCI_mean,mMCI_mean,AD_mean))
C_number <- rep('1',21)
SCD_number <- rep('2',40)
aMCI_number <- rep('3',29)
mMCI_number <- rep('4',19)
AD_number <- rep('5',10)
all_groups_mean$Group<-(c(C_number,SCD_number,aMCI_number,mMCI_number,AD_number))
all_groups_mean$Group<-as.factor(all_groups_mean$Group)


#F-test  - is there a difference between the groups?
all_groups_mean_mod <- lm(`change_in_FBC(post-pre)` ~ Group, data = all_groups_mean)
anova(all_groups_mean_mod)
etaSquared(all_groups_mean_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_all_groups_mean_mod <- glht(all_groups_mean_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_all_groups_mean_mod)
confint(post_hoc_all_groups_mean_mod)
#effect size for sig. post hoc tests
#for AD vs. SCD
all_groups_mean_ADvSCD_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 5 | all_groups_mean$Group == 2)
all_groups_mean_ADvSCD_whole_brain$Group <- droplevels(all_groups_mean_ADvSCD_whole_brain$Group)
cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_ADvSCD_whole_brain) #this looks like Hedges' g? 
#for AD vs. aMCI
all_groups_mean_ADvaMCI_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 5 | all_groups_mean$Group == 3)
all_groups_mean_ADvaMCI_whole_brain$Group <- droplevels(all_groups_mean_ADvaMCI_whole_brain$Group)
cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_ADvaMCI_whole_brain) #this looks like Hedges' g? 
#for AD vs. mMCI
all_groups_mean_ADvmMCI_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 5 | all_groups_mean$Group == 4)
all_groups_mean_ADvmMCI_whole_brain$Group <- droplevels(all_groups_mean_ADvmMCI_whole_brain$Group)
cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_ADvmMCI_whole_brain) #this looks like Hedges' g? 
#for C vs. AD
all_groups_mean_CvAD_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 1 | all_groups_mean$Group == 5)
all_groups_mean_CvAD_whole_brain$Group <- droplevels(all_groups_mean_CvAD_whole_brain$Group)
cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_CvAD_whole_brain) #this looks like Hedges' g? 

#get confidence intervals for each group
group_95CI_model <- lm(`change_in_FBC(post-pre)` ~ Group-1, all_groups_mean)
CI_95_group <- confint(group_95CI_model,level=0.95)
CI_95_group<-as.data.frame(CI_95_group)
#https://stats.stackexchange.com/questions/210515/confidence-intervals-for-group-means-r
# #modify to get the correct the group CIs (intercept + beta estimates)
# CI_95_group$`2.5 %`[2]<-CI_95_group$`2.5 %`[2]+CI_95_group$`2.5 %`[1]
# CI_95_group$`2.5 %`[3]<-CI_95_group$`2.5 %`[3]+CI_95_group$`2.5 %`[1]
# CI_95_group$`2.5 %`[4]<-CI_95_group$`2.5 %`[4]+CI_95_group$`2.5 %`[1]
# CI_95_group$`2.5 %`[5]<-CI_95_group$`2.5 %`[5]+CI_95_group$`2.5 %`[1]
# CI_95_group$`97.5 %`[2]<-CI_95_group$`97.5 %`[2]+CI_95_group$`97.5 %`[1]
# CI_95_group$`97.5 %`[3]<-CI_95_group$`97.5 %`[3]+CI_95_group$`97.5 %`[1]
# CI_95_group$`97.5 %`[4]<-CI_95_group$`97.5 %`[4]+CI_95_group$`97.5 %`[1]
# CI_95_group$`97.5 %`[5]<-CI_95_group$`97.5 %`[5]+CI_95_group$`97.5 %`[1]

#plot data in R
mean_group_data <-data.frame(matrix(ncol=4,nrow=5))
colnames(mean_group_data) <- c('Group','average_change_in_FBC(post-pre)','ymin','ymax')
mean_group_data$`average_change_in_FBC(post-pre)`<-tapply(all_groups_mean$`change_in_FBC(post-pre)`,all_groups_mean$Group,mean)
mean_group_data$Group<-c('1','2','3','4','5')
mean_group_data$Group<-as.factor(mean_group_data$Group)
mean_group_data$ymin <-CI_95_group$`2.5 %`
mean_group_data$ymax <-CI_95_group$`97.5 %`

#bar graph
ggplot(mean_group_data, aes(x = Group, y = `average_change_in_FBC(post-pre)`, fill = Group)) + 
  geom_bar(position = "dodge",stat="identity")+
  geom_errorbar( aes(x=Group, ymin=ymin, ymax=ymax),width=0.4,alpha=0.9,size=1.3)+
  ylab("Difference scores in structural connectivity") +
  scale_x_discrete(labels = c("C","SCD","aMCI","mMCI","AD"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(size=18))+
  theme(axis.title=element_text(size=24))

#violin plot
ggplot(all_groups_mean, aes(x = Group, y = `change_in_FBC(post-pre)`, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = `change_in_FBC(post-pre)`, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Difference scores in structural connectivity") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=18))+
  theme(axis.title=element_text(size=24))#+
#coord_flip()


#add participant ID to dataframe - need to run some code from below for this
mod_dataframe <- DPRC_neuropsych_data[seq(1, nrow(DPRC_neuropsych_data), 2), ]
all_groups_mean$ParticipantID <- mod_dataframe$ParticipantID

#find max FBC edge value by group
result_max <- all_groups_mean %>% 
  group_by(Group) %>%
  filter(`change_in_FBC(post-pre)` == max(`change_in_FBC(post-pre)`))

# #find min FBC edge value by group
#   result_min <- dataframe_edge_compare %>% 
#   group_by(Group) %>%
#   filter(FBC_edge_value == min(FBC_edge_value))

#find min FBC edge value by group (non-zero)
result_min <- all_groups_mean %>% 
  group_by(Group) %>%
  filter(`change_in_FBC(post-pre)`== min(`change_in_FBC(post-pre)`[`change_in_FBC(post-pre)`>0]))




#show interaction plot (F0 mean and F2 mean plotted)


setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')
#for functional-diffusion connectome, longitudinal study (n = 119)
DPRC_neuropsych_data <- read.csv("fMRI_dMRI_longitudinal_DPRC_neuropsych_data_lined_up_valid_participants.csv")
#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'
#Add in longitudinal values for participants
Individual_number <- c(1:119, 1:119)
DPRC_neuropsych_data$Individual_number <- Individual_number
#convert variables
DPRC_neuropsych_data$ParticipantID <- as.factor(DPRC_neuropsych_data$ParticipantID)
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
#add variables to a new dataframe
interaction_data <- subset(DPRC_neuropsych_data, select=c('ParticipantID','Group'))

#go into directory with F0 and F2 files
setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/thresholded_connectomes/C&SCD')

#extract values from the sig. edges (e.g., 11) from each F0 and F2 file, and put into a dataframe
#read in all participant .csv files
all_csvs <- lapply(list.files(), read.csv, header=FALSE)

#extract edge from node 27 and 30 for each participant and put into a dataframe
edges_value_vec1 <- unlist(lapply(all_csvs, `[`, 2, 27))
edges_value_vec2 <- unlist(lapply(all_csvs, `[`, 2,	77))
edges_value_vec3 <- unlist(lapply(all_csvs, `[`, 3,	10))
edges_value_vec4 <- unlist(lapply(all_csvs, `[`, 3,	27))
edges_value_vec5 <- unlist(lapply(all_csvs, `[`, 3,	71))
edges_value_vec6 <- unlist(lapply(all_csvs, `[`, 9,	77))
edges_value_vec7 <- unlist(lapply(all_csvs, `[`, 9,	79))
edges_value_vec8 <- unlist(lapply(all_csvs, `[`, 10, 27))
edges_value_vec9 <- unlist(lapply(all_csvs, `[`, 10, 74))
edges_value_vec10 <- unlist(lapply(all_csvs, `[`, 10, 75))
edges_value_vec11 <- unlist(lapply(all_csvs, `[`, 11,	61))
edges_value_vec12 <- unlist(lapply(all_csvs, `[`, 14,	23))
edges_value_vec13 <- unlist(lapply(all_csvs, `[`, 20,	40))
edges_value_vec14 <- unlist(lapply(all_csvs, `[`, 20,	41))
edges_value_vec15 <- unlist(lapply(all_csvs, `[`, 22,	50))
edges_value_vec16 <- unlist(lapply(all_csvs, `[`, 22,	55))
edges_value_vec17 <- unlist(lapply(all_csvs, `[`, 22,	57))
edges_value_vec18 <- unlist(lapply(all_csvs, `[`, 22,	77))
edges_value_vec19 <- unlist(lapply(all_csvs, `[`, 23,	71))
edges_value_vec20 <- unlist(lapply(all_csvs, `[`, 27,	37))
edges_value_vec21 <- unlist(lapply(all_csvs, `[`, 27,	47))
edges_value_vec22 <- unlist(lapply(all_csvs, `[`, 27,	63))
edges_value_vec23 <- unlist(lapply(all_csvs, `[`, 27,	79))
edges_value_vec24 <- unlist(lapply(all_csvs, `[`, 28,	30))
edges_value_vec25 <- unlist(lapply(all_csvs, `[`, 28,	55))
edges_value_vec26 <- unlist(lapply(all_csvs, `[`, 29,	62))
edges_value_vec27 <- unlist(lapply(all_csvs, `[`, 29,	73))
edges_value_vec28 <- unlist(lapply(all_csvs, `[`, 29,	74))
edges_value_vec29 <- unlist(lapply(all_csvs, `[`, 36,	61))
edges_value_vec30 <- unlist(lapply(all_csvs, `[`, 37,	41))
edges_value_vec31 <- unlist(lapply(all_csvs, `[`, 37,	56))
edges_value_vec32 <- unlist(lapply(all_csvs, `[`, 37,	64))
edges_value_vec33 <- unlist(lapply(all_csvs, `[`, 41,	57))
edges_value_vec34 <- unlist(lapply(all_csvs, `[`, 52,	79))
edges_value_vec35 <- unlist(lapply(all_csvs, `[`, 55,	56))
edges_value_vec36 <- unlist(lapply(all_csvs, `[`, 55,	70))
edges_value_vec37 <- unlist(lapply(all_csvs, `[`, 55,	78))
edges_value_vec38 <- unlist(lapply(all_csvs, `[`, 56,	58))
edges_value_vec39 <- unlist(lapply(all_csvs, `[`, 56,	61))
edges_value_vec40 <- unlist(lapply(all_csvs, `[`, 57,	58))
edges_value_vec41 <- unlist(lapply(all_csvs, `[`, 57,	69))
edges_value_vec42 <- unlist(lapply(all_csvs, `[`, 62,	64))
edges_value_vec43<- unlist(lapply(all_csvs, `[`, 73,	82))
edges_value_vec44 <- unlist(lapply(all_csvs, `[`, 74,	82))


timepoint <- as.factor(rep(c('F0','F2'),119))

interaction_data$Timepoint <- timepoint 
interaction_data$edges_value_vec1 <- edges_value_vec1
interaction_data$edges_value_vec2 <- edges_value_vec2
interaction_data$edges_value_vec3 <- edges_value_vec3
interaction_data$edges_value_vec4 <- edges_value_vec4
interaction_data$edges_value_vec5 <- edges_value_vec5
interaction_data$edges_value_vec6 <- edges_value_vec6
interaction_data$edges_value_vec7 <- edges_value_vec7
interaction_data$edges_value_vec8 <- edges_value_vec8
interaction_data$edges_value_vec9 <- edges_value_vec9
interaction_data$edges_value_vec10 <- edges_value_vec10
interaction_data$edges_value_vec11 <- edges_value_vec11
interaction_data$edges_value_vec12 <- edges_value_vec12
interaction_data$edges_value_vec13 <- edges_value_vec13
interaction_data$edges_value_vec14<- edges_value_vec14
interaction_data$edges_value_vec15 <- edges_value_vec15
interaction_data$edges_value_vec16<- edges_value_vec16
interaction_data$edges_value_vec17<- edges_value_vec17
interaction_data$edges_value_vec18<- edges_value_vec18
interaction_data$edges_value_vec19 <- edges_value_vec19
interaction_data$edges_value_vec20 <- edges_value_vec20
interaction_data$edges_value_vec21 <- edges_value_vec21
interaction_data$edges_value_vec22 <- edges_value_vec22
interaction_data$edges_value_vec23<- edges_value_vec23
interaction_data$edges_value_vec24 <- edges_value_vec24
interaction_data$edges_value_vec25 <- edges_value_vec25
interaction_data$edges_value_vec26 <- edges_value_vec26
interaction_data$edges_value_vec27<- edges_value_vec27
interaction_data$edges_value_vec28 <- edges_value_vec28
interaction_data$edges_value_vec29<- edges_value_vec29
interaction_data$edges_value_vec30 <- edges_value_vec30
interaction_data$edges_value_vec31 <- edges_value_vec31
interaction_data$edges_value_vec32 <- edges_value_vec32
interaction_data$edges_value_vec33 <- edges_value_vec33
interaction_data$edges_value_vec34 <- edges_value_vec34
interaction_data$edges_value_vec35 <- edges_value_vec35
interaction_data$edges_value_vec36 <- edges_value_vec36
interaction_data$edges_value_vec37 <- edges_value_vec37
interaction_data$edges_value_vec38 <- edges_value_vec38
interaction_data$edges_value_vec39 <- edges_value_vec39
interaction_data$edges_value_vec40 <- edges_value_vec40
interaction_data$edges_value_vec41 <- edges_value_vec41
interaction_data$edges_value_vec42 <- edges_value_vec42
interaction_data$edges_value_vec43 <- edges_value_vec43
interaction_data$edges_value_vec44 <- edges_value_vec44

#calculate the mean of the values at each edge for each time point
interaction_data$avg_edge_values<-rowMeans(interaction_data[4:47])

#view interaction plot - by Group
interaction_data%>%
  group_by(Group,Timepoint) %>%
  summarise(s_mean=mean(avg_edge_values)) %>%
  ggplot(aes(y=s_mean,x=Group,colour=Timepoint,group=Timepoint))+
  geom_point()+geom_line()+
  scale_x_discrete(labels = c("1" = "C", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  ylab("Mean of Edges with a Significant Interaction") +
  theme_classic()
#view interaction plot - by Timepoint
interaction_data %>%
  group_by(Group,Timepoint) %>%
  summarise(s_mean=mean(avg_edge_values)) %>%
  ggplot(aes(y=s_mean,x=Timepoint,colour=Group,group=Group))+
  geom_point()+geom_line()+
  scale_color_discrete(labels = c("1" = "C", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  ylab("Mean of Edges with a Significant Interaction") +
  theme_classic()
#note that if these interaction plots are not working, just try restarting R (just exit and open again)






#------------------------------------------------------------------------------#


#w/ covariates (age and sex) for linear trend
#for the whole brain structural connectome - 38 significant edges for interaction
#load in data, copied from matlab
C_mean <- c(0.00265353780259572,
            -0.00656476339035109,
            -0.000277441054516318,
            0.00172244245682187,
            -0.000336195368531119,
            0.00597145915105930,
            0.000539535265853398,
            0.000924228410569799,
            -0.00622488441827231,
            0.0114569016700612,
            0.0103016945275209,
            0.0115475348715839,
            0.00326625486540041,
            0.00275689854136931,
            0.0185033684364213,
            0.000358652125748565,
            -0.000946199460819322,
            0.00629695641881915,
            0.00335461688884080,
            0.00139221412803477,
            0.0170022628572510) #n = 21

SCD_mean <- c(-0.00335529845135799,
              0.0103567018438050,
              0.00844810679819172,
              0.00270675667719352,
              -0.00213030156459130,
              0.000865236512665461,
              0.00237655910118666,
              0.00519368004282142,
              0.000770943079104066,
              -0.00446114538419999,
              -0.00281367004745854,
              -0.000961578285023861,
              0.00217633876786494,
              0.00465644895534084,
              0.00204944705494493,
              -0.000276102805186983,
              0.00558118814555344,
              0.00557986565843316,
              -0.000924319945345439,
              0.00525245683875717,
              -0.00132559147180554,
              0.00199210078283471,
              0.00702406686654072,
              0.00792112912739240,
              0.00275141581479378,
              0.00788394320569921,
              0.0114837851664699,
              0.000466125434692257,
              0.00178363237664615,
              -0.00586624556490831,
              0.00307130093856007,
              0.00319357580887845,
              0.0113380802690819,
              -0.00103747078232608,
              -0.00412034997118712,
              -0.000466269957287914,
              -3.14096544884594e-05,
              0.000471500976490366,
              0.00202115414307686,
              0.0116166336164860) #n = 40

aMCI_mean <- c(0.00875792070648110,
               0.000388216612687682,
               0.00309037355477867,
               -0.000300356316719899,
               -0.00252194257742920,
               0.00710558362454269,
               0.00880792565081406,
               -0.0108346955233797,
               0.00905194669440512,
               -0.00107854623319827,
               0.00166899115179908,
               0.00507057002143998,
               0.00357019030225512,
               0.00137646422872308,
               0.00683787467774972,
               0.00356485930791256,
               0.00613420553011260,
               -0.00514001587769819,
               -0.00107827256154722,
               0.00548268510667299,
               0.00199321517406824,
               0.0167605794644395,
               0.00219475404215274,
               0.00272958338381934,
               0.00308902281538994,
               -0.00195099328475901,
               0.0158502795221027,
               0.00592632415238352,
               -0.0251622965560533) #n = 29

mMCI_mean <- c(0.00667210744720650,
               0.0116972464842370,
               0.00156979229349245,
               0.00143608855342809,
               0.00235782182434532,
               0.000415648963853582,
               0.00138595237663486,
               0.00146466640388565,
               0.00484396787584111,
               0.00135938156493655,
               -0.000951562862579106,
               -0.00363873686961848,
               0.00986141940388822,
               -0.00285676738328928,
               -0.00325233243857697,
               -0.0164446551506377,
               -0.00128481720087648,
               -0.000116420115099004,
               0.00183353975973553) #n = 19

AD_mean <- c(-0.00856440766222751,
             -0.0110896128951514,
             -0.00182029501932189,
             -0.0128405545803739,
             -0.00716398232291320,
             -0.000929666029028751,
             -0.0213297742401921,
             -0.00808096235795423,
             -0.00316877204950230,
             -0.0198376437993570) #n = 10

#combine into one dataframe
all_groups_mean <-data.frame(matrix(ncol=2,nrow=119))
colnames(all_groups_mean) <- c('Group','change_in_FBC(post-pre)')
all_groups_mean$`change_in_FBC(post-pre)`<-(c(C_mean,SCD_mean,aMCI_mean,mMCI_mean,AD_mean))
C_number <- rep('1',21)
SCD_number <- rep('2',40)
aMCI_number <- rep('3',29)
mMCI_number <- rep('4',19)
AD_number <- rep('5',10)
all_groups_mean$Group<-(c(C_number,SCD_number,aMCI_number,mMCI_number,AD_number))
all_groups_mean$Group<-as.factor(all_groups_mean$Group)


#F-test  - is there a difference between the groups?
all_groups_mean_mod <- lm(`change_in_FBC(post-pre)` ~ Group, data = all_groups_mean)
anova(all_groups_mean_mod)
etaSquared(all_groups_mean_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_all_groups_mean_mod <- glht(all_groups_mean_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_all_groups_mean_mod)
confint(post_hoc_all_groups_mean_mod)
#effect size for sig. post hoc tests
  #for AD vs. SCD
  all_groups_mean_ADvSCD_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 5 | all_groups_mean$Group == 2)
  all_groups_mean_ADvSCD_whole_brain$Group <- droplevels(all_groups_mean_ADvSCD_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_ADvSCD_whole_brain) #this looks like Hedges' g? 
  #for AD vs. aMCI
  all_groups_mean_ADvaMCI_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 5 | all_groups_mean$Group == 3)
  all_groups_mean_ADvaMCI_whole_brain$Group <- droplevels(all_groups_mean_ADvaMCI_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_ADvaMCI_whole_brain) #this looks like Hedges' g? 
  #for AD vs. mMCI
  all_groups_mean_ADvmMCI_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 5 | all_groups_mean$Group == 4)
  all_groups_mean_ADvmMCI_whole_brain$Group <- droplevels(all_groups_mean_ADvmMCI_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_ADvmMCI_whole_brain) #this looks like Hedges' g? 
  #for C vs. AD
  all_groups_mean_CvAD_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 1 | all_groups_mean$Group == 5)
  all_groups_mean_CvAD_whole_brain$Group <- droplevels(all_groups_mean_CvAD_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_CvAD_whole_brain) #this looks like Hedges' g? 

#get confidence intervals for each group
group_95CI_model <- lm(`change_in_FBC(post-pre)` ~ Group-1, all_groups_mean)
CI_95_group <- confint(group_95CI_model,level=0.95)
CI_95_group<-as.data.frame(CI_95_group)
#https://stats.stackexchange.com/questions/210515/confidence-intervals-for-group-means-r
# #modify to get the correct the group CIs (intercept + beta estimates)
# CI_95_group$`2.5 %`[2]<-CI_95_group$`2.5 %`[2]+CI_95_group$`2.5 %`[1]
# CI_95_group$`2.5 %`[3]<-CI_95_group$`2.5 %`[3]+CI_95_group$`2.5 %`[1]
# CI_95_group$`2.5 %`[4]<-CI_95_group$`2.5 %`[4]+CI_95_group$`2.5 %`[1]
# CI_95_group$`2.5 %`[5]<-CI_95_group$`2.5 %`[5]+CI_95_group$`2.5 %`[1]
# CI_95_group$`97.5 %`[2]<-CI_95_group$`97.5 %`[2]+CI_95_group$`97.5 %`[1]
# CI_95_group$`97.5 %`[3]<-CI_95_group$`97.5 %`[3]+CI_95_group$`97.5 %`[1]
# CI_95_group$`97.5 %`[4]<-CI_95_group$`97.5 %`[4]+CI_95_group$`97.5 %`[1]
# CI_95_group$`97.5 %`[5]<-CI_95_group$`97.5 %`[5]+CI_95_group$`97.5 %`[1]

#plot data in R
mean_group_data <-data.frame(matrix(ncol=4,nrow=5))
colnames(mean_group_data) <- c('Group','average_change_in_FBC(post-pre)','ymin','ymax')
mean_group_data$`average_change_in_FBC(post-pre)`<-tapply(all_groups_mean$`change_in_FBC(post-pre)`,all_groups_mean$Group,mean)
mean_group_data$Group<-c('1','2','3','4','5')
mean_group_data$Group<-as.factor(mean_group_data$Group)
mean_group_data$ymin <-CI_95_group$`2.5 %`
mean_group_data$ymax <-CI_95_group$`97.5 %`

#bar graph
ggplot(mean_group_data, aes(x = Group, y = `average_change_in_FBC(post-pre)`, fill = Group)) + 
  geom_bar(position = "dodge",stat="identity")+
  geom_errorbar( aes(x=Group, ymin=ymin, ymax=ymax),width=0.4,alpha=0.9,size=1.3)+
  ylab("Difference scores in structural connectivity") +
  scale_x_discrete(labels = c("C","SCD","aMCI","mMCI","AD"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(size=18))+
  theme(axis.title=element_text(size=24))

#violin plot
ggplot(all_groups_mean, aes(x = Group, y = `change_in_FBC(post-pre)`, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = `change_in_FBC(post-pre)`, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Difference scores in structural connectivity") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=18))+
  theme(axis.title=element_text(size=24))#+
#coord_flip()


#add participant ID to dataframe - need to run some code from below for this
mod_dataframe <- DPRC_neuropsych_data[seq(1, nrow(DPRC_neuropsych_data), 2), ]
all_groups_mean$ParticipantID <- mod_dataframe$ParticipantID

#find max FBC edge value by group
result_max <- all_groups_mean %>% 
  group_by(Group) %>%
  filter(`change_in_FBC(post-pre)` == max(`change_in_FBC(post-pre)`))

# #find min FBC edge value by group
#   result_min <- dataframe_edge_compare %>% 
#   group_by(Group) %>%
#   filter(FBC_edge_value == min(FBC_edge_value))

#find min FBC edge value by group (non-zero)
result_min <- all_groups_mean %>% 
  group_by(Group) %>%
  filter(`change_in_FBC(post-pre)`== min(`change_in_FBC(post-pre)`[`change_in_FBC(post-pre)`>0]))




#show interaction plot (F0 mean and F2 mean plotted)

setwd('H:/ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/')
#for functional-diffusion connectome, longitudinal study (n = 119)
DPRC_neuropsych_data <- read.csv("fMRI_dMRI_longitudinal_DPRC_neuropsych_data_lined_up_valid_participants.csv")
#rename first column 
colnames(DPRC_neuropsych_data)[1] <-'ParticipantID'
#Add in longitudinal values for participants
Individual_number <- c(1:119, 1:119)
DPRC_neuropsych_data$Individual_number <- Individual_number
#convert variables
DPRC_neuropsych_data$ParticipantID <- as.factor(DPRC_neuropsych_data$ParticipantID)
DPRC_neuropsych_data$Group <- as.factor(DPRC_neuropsych_data$Group)
#add variables to a new dataframe
interaction_data <- subset(DPRC_neuropsych_data, select=c('ParticipantID','Group'))

#go into directory with F0 and F2 files
setwd('V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/dkFiles/weighted/fs_default_ordered/thresholded_connectomes/C&SCD')

#extract values from the sig. edges (e.g., 11) from each F0 and F2 file, and put into a dataframe
#read in all participant .csv files
all_csvs <- lapply(list.files(), read.csv, header=FALSE)

#extract edge from node 27 and 30 for each participant and put into a dataframe
edges_value_vec1 <- unlist(lapply(all_csvs, `[`, 2, 9))
edges_value_vec2 <- unlist(lapply(all_csvs, `[`, 2,	27))
edges_value_vec3 <- unlist(lapply(all_csvs, `[`, 3,	27))
edges_value_vec4 <- unlist(lapply(all_csvs, `[`, 9,	77))
edges_value_vec5 <- unlist(lapply(all_csvs, `[`, 10, 74))
edges_value_vec6 <- unlist(lapply(all_csvs, `[`, 10, 75))
edges_value_vec7 <- unlist(lapply(all_csvs, `[`, 11, 61))
edges_value_vec8 <- unlist(lapply(all_csvs, `[`, 13, 61))
edges_value_vec9 <- unlist(lapply(all_csvs, `[`, 14,23))
edges_value_vec10 <- unlist(lapply(all_csvs, `[`, 20,	40))
edges_value_vec11 <- unlist(lapply(all_csvs, `[`, 22,	50))
edges_value_vec12 <- unlist(lapply(all_csvs, `[`, 22,	55))
edges_value_vec13 <- unlist(lapply(all_csvs, `[`, 22,	57))
edges_value_vec14 <- unlist(lapply(all_csvs, `[`, 22,	70))
edges_value_vec15 <- unlist(lapply(all_csvs, `[`, 22,	77))
edges_value_vec16 <- unlist(lapply(all_csvs, `[`, 23,	71))
edges_value_vec17 <- unlist(lapply(all_csvs, `[`, 26,	27))
edges_value_vec18 <- unlist(lapply(all_csvs, `[`, 26,	79))
edges_value_vec19 <- unlist(lapply(all_csvs, `[`, 27,	63))
edges_value_vec20 <- unlist(lapply(all_csvs, `[`, 27,	79))
edges_value_vec21 <- unlist(lapply(all_csvs, `[`, 28,	38))
edges_value_vec22 <- unlist(lapply(all_csvs, `[`, 28,	55))
edges_value_vec23 <- unlist(lapply(all_csvs, `[`, 29,	73))
edges_value_vec24 <- unlist(lapply(all_csvs, `[`, 29,	74))
edges_value_vec25 <- unlist(lapply(all_csvs, `[`, 37,	61))
edges_value_vec26 <- unlist(lapply(all_csvs, `[`, 37,	64))
edges_value_vec27 <- unlist(lapply(all_csvs, `[`, 40,	45))
edges_value_vec28 <- unlist(lapply(all_csvs, `[`, 43,	59))
edges_value_vec29 <- unlist(lapply(all_csvs, `[`, 43,	61))
edges_value_vec30 <- unlist(lapply(all_csvs, `[`, 52,	79))
edges_value_vec31 <- unlist(lapply(all_csvs, `[`, 55,	56))
edges_value_vec32 <- unlist(lapply(all_csvs, `[`, 55,	70))
edges_value_vec33 <- unlist(lapply(all_csvs, `[`, 55,	78))
edges_value_vec34 <- unlist(lapply(all_csvs, `[`, 57,	58))
edges_value_vec35 <- unlist(lapply(all_csvs, `[`, 57,	69))
edges_value_vec36 <- unlist(lapply(all_csvs, `[`, 62,	64))
edges_value_vec37 <- unlist(lapply(all_csvs, `[`, 73,	82))
edges_value_vec38 <- unlist(lapply(all_csvs, `[`, 74,	82))


timepoint <- as.factor(rep(c('F0','F2'),119))

interaction_data$Timepoint <- timepoint 
interaction_data$edges_value_vec1 <- edges_value_vec1
interaction_data$edges_value_vec2 <- edges_value_vec2
interaction_data$edges_value_vec3 <- edges_value_vec3
interaction_data$edges_value_vec4 <- edges_value_vec4
interaction_data$edges_value_vec5 <- edges_value_vec5
interaction_data$edges_value_vec6 <- edges_value_vec6
interaction_data$edges_value_vec7 <- edges_value_vec7
interaction_data$edges_value_vec8 <- edges_value_vec8
interaction_data$edges_value_vec9 <- edges_value_vec9
interaction_data$edges_value_vec10 <- edges_value_vec10
interaction_data$edges_value_vec11 <- edges_value_vec11
interaction_data$edges_value_vec12 <- edges_value_vec12
interaction_data$edges_value_vec13 <- edges_value_vec13
interaction_data$edges_value_vec14<- edges_value_vec14
interaction_data$edges_value_vec15 <- edges_value_vec15
interaction_data$edges_value_vec16<- edges_value_vec16
interaction_data$edges_value_vec17<- edges_value_vec17
interaction_data$edges_value_vec18<- edges_value_vec18
interaction_data$edges_value_vec19 <- edges_value_vec19
interaction_data$edges_value_vec20 <- edges_value_vec20
interaction_data$edges_value_vec21 <- edges_value_vec21
interaction_data$edges_value_vec22 <- edges_value_vec22
interaction_data$edges_value_vec23<- edges_value_vec23
interaction_data$edges_value_vec24 <- edges_value_vec24
interaction_data$edges_value_vec25 <- edges_value_vec25
interaction_data$edges_value_vec26 <- edges_value_vec26
interaction_data$edges_value_vec27<- edges_value_vec27
interaction_data$edges_value_vec28 <- edges_value_vec28
interaction_data$edges_value_vec29<- edges_value_vec29
interaction_data$edges_value_vec30 <- edges_value_vec30
interaction_data$edges_value_vec31 <- edges_value_vec31
interaction_data$edges_value_vec32 <- edges_value_vec32
interaction_data$edges_value_vec33 <- edges_value_vec33
interaction_data$edges_value_vec34 <- edges_value_vec34
interaction_data$edges_value_vec35 <- edges_value_vec35
interaction_data$edges_value_vec36 <- edges_value_vec36
interaction_data$edges_value_vec37 <- edges_value_vec37
interaction_data$edges_value_vec38 <- edges_value_vec38

#calculate the mean of the values at each edge for each time point
interaction_data$avg_edge_values<-rowMeans(interaction_data[4:41])

#view interaction plot - by Group
interaction_data%>%
  group_by(Group,Timepoint) %>%
  summarise(s_mean=mean(avg_edge_values)) %>%
  ggplot(aes(y=s_mean,x=Group,colour=Timepoint,group=Timepoint))+
  geom_point()+geom_line()+
  scale_x_discrete(labels = c("1" = "C", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  ylab("Mean of Edges with a Significant Interaction") +
  theme_classic()
#view interaction plot - by Timepoint
interaction_data %>%
  group_by(Group,Timepoint) %>%
  summarise(s_mean=mean(avg_edge_values)) %>%
  ggplot(aes(y=s_mean,x=Timepoint,colour=Group,group=Group))+
  geom_point()+geom_line()+
  scale_color_discrete(labels = c("1" = "C", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  ylab("Mean of Edges with a Significant Interaction") +
  theme_classic()
#note that if these interaction plots are not working, just try restarting R (just exit and open again)





#------------------------------------------------------------------------------#
#for the FPN - 1 sig. edge from the interaction. (edge 413, nodes 277 & 264)
#load in data, copied from matlab
C_FPN_edge <- c(0.0548476831112877,
                -0.0737918599792347,
                0.241845791208232,
                0.281086493842325,
                0.0620371508855180,
                0.0700130048840030,
                0.541613186971485,
                0.439539806573493,
                0.000359019940259799,
                -0.0514109921284823,
                0.166512270631744,
                -0.00674132459989129,
                0.202746911222031,
                1.42702260702475,
                0.0180264992395908,
                -0.201224537144799,
                0.0292402573043198,
                0.171723482532596,
                0.0565895982448833,
                -0.0796482380871203,
                -0.251124921338215) #n = 21

SCD_FPN_edge <- c(-0.0738221863145670,
                  0.168964640186742,
                  0.0861096102674984,
                  0.0308044668937782,
                  -0.00654716909978298,
                  0.0665242683472012,
                  -0.560705182837781,
                  0.0181893129151229,
                  0.0292473752706022,
                  -0.0413093119712940,
                  -0.177254038607383,
                  0.294642274898383,
                  -0.0196277360561550,
                  0.294580870327078,
                  0.0639439458082949,
                  -0.0650793822311156,
                  0.0450062186011070,
                  0.431884998428502,
                  -0.153493380514991,
                  0.00576117989080333,
                  -0.0974117978524100,
                  0.0832971408040227,
                  -0.0894581703394819,
                  0.140681913875138,
                  0.0192827400890047,
                  0.0936979525477953,
                  0.00172479757470601,
                  -0.0113243618697259,
                  0.177031403226335,
                  -0.0502638361012447,
                  -0.0268671741014451,
                  0.0437572072188370,
                  0.0867023021115070,
                  -0.0901745778940147,
                  0.366840847195705,
                  -0.0241952056016852,
                  0.0196664302599878,
                  -0.341176043340943,
                  0.0583519207392218,
                  0.0744098213059021) #n = 40

aMCI_FPN_edge <- c(0.0354048157833929,
                   -0.208351870903656,
                   -0.0662569969052511,
                   0.0912499143498702,
                   0.0846444789711550,
                   0.0283313531938812,
                   -0.122528434677613,
                   0.210917619741797,
                   0.309935937779138,
                   -0.0855037019801297,
                   -0.691953628738151,
                   0.0228019138464064,
                   0.0900456411436156,
                   0.0484910633932820,
                   0.151845212229464,
                   0.120275170102829,
                   0.0971994218420950,
                   -0.120394927003093,
                   -0.141458750637434,
                   -0.0451906507246960,
                   0.0776714441086476,
                   0.147208013348318,
                   -0.0297398849545294,
                   0.168232961414800,
                   -0.00983132752029449,
                   0.223089661914498,
                   -0.105115827555295,
                   -0.199641543517829,
                   -0.160797811649492) #n = 29

mMCI_FPN_edge <- c(-0.0449648852147440,
                   -0.153939303635131,
                   0.336507506924342,
                   -0.0879352942574564,
                   -0.0995872454317870,
                   -0.0562027974684190,
                   0.309709148896969,
                   -0.0835073401773910,
                   -0.0850757984433304,
                   -0.000396642617176200,
                   -0.0616633414913532,
                   -0.108409083785915,
                   -0.114709040538779,
                   0.417350073327885,
                   0.00400814731226330,
                   -0.145097041355084,
                   -0.0686138456096144,
                   0.0771908916178990,
                   -0.174491037528450) #n = 19

AD_FPN_edge <- c(0.172506069247830,
                 -0.380480807866817,
                 -0.0178615354932872,
                 -0.253063619855387,
                 -0.347997611444715,
                 -0.0706081916159810,
                 -1.09850610482048,
                 -0.550188002474992,
                 -0.0155057492199916,
                 -0.0659540393607574) #n = 10

#combine into one dataframe
all_groups_FPN_edge <-data.frame(matrix(ncol=2,nrow=119))
colnames(all_groups_FPN_edge) <- c('Group','change_in_FBC(post-pre)')
all_groups_FPN_edge$`change_in_FBC(post-pre)`<-(c(C_FPN_edge,SCD_FPN_edge,aMCI_FPN_edge,mMCI_FPN_edge,AD_FPN_edge))
C_number <- rep('1',21)
SCD_number <- rep('2',40)
aMCI_number <- rep('3',29)
mMCI_number <- rep('4',19)
AD_number <- rep('5',10)
all_groups_FPN_edge$Group<-(c(C_number,SCD_number,aMCI_number,mMCI_number,AD_number))
all_groups_FPN_edge$Group<-as.factor(all_groups_FPN_edge$Group)


#F-test  - is there a difference between the groups?
all_groups_FPN_edge_mod <- lm(`change_in_FBC(post-pre)` ~ Group, data = all_groups_FPN_edge)
anova(all_groups_FPN_edge_mod)
etaSquared(all_groups_FPN_edge_mod)
#run pairwise comparisons (post-hoc Tukey), given that the F-test was significant. 
post_hoc_all_groups_FPN_edge_mod <- glht(all_groups_FPN_edge_mod, linfct = mcp(Group = "Tukey"))
summary(post_hoc_all_groups_FPN_edge_mod)
confint(post_hoc_all_groups_FPN_edge_mod)
#effect size for sig. post hoc tests
  #for AD vs. SCD
  all_groups_FPN_edge_ADvSCD_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 5 | all_groups_FPN_edge$Group == 2)
  all_groups_FPN_edge_ADvSCD_FPN$Group <- droplevels(all_groups_FPN_edge_ADvSCD_FPN$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_ADvSCD_FPN) #this looks like Hedges' g? 
  #for AD vs. aMCI
  all_groups_FPN_edge_ADvaMCI_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 5 | all_groups_FPN_edge$Group == 3)
  all_groups_FPN_edge_ADvaMCI_FPN$Group <- droplevels(all_groups_FPN_edge_ADvaMCI_FPN$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_ADvaMCI_FPN) #this looks like Hedges' g? 
  #for AD vs. mMCI
  all_groups_FPN_edge_ADvmMCI_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 5 | all_groups_FPN_edge$Group == 4)
  all_groups_FPN_edge_ADvmMCI_FPN$Group <- droplevels(all_groups_FPN_edge_ADvmMCI_FPN$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_ADvmMCI_FPN) #this looks like Hedges' g? 
  #for C vs. AD
  all_groups_FPN_edge_CvAD_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 1 | all_groups_FPN_edge$Group == 5)
  all_groups_FPN_edge_CvAD_FPN$Group <- droplevels(all_groups_FPN_edge_CvAD_FPN$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_CvAD_FPN) #this looks like Hedges' g? 
#get confidence intervals for each group
group_95CI_model <- lm(`change_in_FBC(post-pre)` ~ Group-1, all_groups_FPN_edge)
CI_95_group <- confint(group_95CI_model,level=0.95)
CI_95_group<-as.data.frame(CI_95_group)

#plot data in R
FPN_edge_group_data <-data.frame(matrix(ncol=4,nrow=5))
colnames(FPN_edge_group_data) <- c('Group','average_change_in_FBC(post-pre)','ymin','ymax')
FPN_edge_group_data$`average_change_in_FBC(post-pre)`<-tapply(all_groups_FPN_edge$`change_in_FBC(post-pre)`,all_groups_FPN_edge$Group,mean)
FPN_edge_group_data$Group<-c('1','2','3','4','5')
FPN_edge_group_data$Group<-as.factor(FPN_edge_group_data$Group)
FPN_edge_group_data$ymin <-CI_95_group$`2.5 %`
FPN_edge_group_data$ymax <-CI_95_group$`97.5 %`

#bar graph
ggplot(FPN_edge_group_data, aes(x = Group, y = `average_change_in_FBC(post-pre)`, fill = Group)) + 
  geom_bar(position = "dodge",stat="identity")+
  geom_errorbar( aes(x=Group, ymin=ymin, ymax=ymax),width=0.4,alpha=0.9,size=1.3)+
  ylab("Difference scores in structural connectivity in the FPN") +
  scale_x_discrete(labels = c("C","SCD","aMCI","mMCI","AD"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(size=18))+
  theme(axis.title=element_text(size=24))

#violin plot
ggplot(all_groups_FPN_edge, aes(x = Group, y = `change_in_FBC(post-pre)`, fill = Group)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = `change_in_FBC(post-pre)`, color = Group), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1, aes(colour = Group)) + 
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, aes(colour = Group)) + 
  xlab("Group") + 
  ylab("Difference scores in structural connectivity") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=18))+
  theme(axis.title=element_text(size=24))#+
#coord_flip()









































