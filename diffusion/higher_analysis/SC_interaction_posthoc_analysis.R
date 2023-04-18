#Post hoc tests for sig. interaction
#run test between groups with the avg. of all nodes

pacman::p_load(ggplot2, multcomp, lsr, dplyr)
#add any necessary sources: 
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R") #for raincloud graph


#for the whole brain structural connectome - 11 significant edges for interaction
#load in data, copied from matlab
C_mean <- c(-0.191808348594555,
            -0.121680380963005,
            -0.0493854491564611,
            -0.0136283508519648,
            -0.0869110590089623,
            -0.0885750102032992,
            -0.0597876379225890,
            -0.0442220392874910,
            -0.0536713651234083,
            -0.0835214370688043,
            -0.0700340948671511,
            -0.0241197085148515,
            -0.0407904991328877,
            -0.0668530627975919,
            -0.0333919233165108,
            -0.00744450221894429,
            0.00450404794591839,
            0.0369173295142060,
            -0.0523310404273705,
            -0.00110135502464539,
            -0.0986051314983299) #n = 21

SCD_mean <- c(0.0137444554800869,
              0.0310796359795206,
              -0.000102900599878620,
              0.0333080026057859,
              0.0176154497070329,
              0.0414593862557899,
              0.0220462504204187,
              0.0193898296466826,
              0.0302520482917501,
              -0.0114046520657980,
              0.00109443226522985,
              0.00742820071591585,
              0.0214960997540635,
              0.0469031892594489,
              -0.0112137177342014,
              -0.0180480188480967,
              0.0239300667985140,
              0.00949547125971047,
              0.0161391453011598,
              0.0131375370543122,
              -0.00420206460185554,
              0.0207919351961424,
              0.00176856926564427,
              0.0544700219292865,
              0.000835940231337294,
              0.0100605973983831,
              0.0263686617898364,
              0.0147669146449233,
              0.0501017880549760,
              -0.0163286826128571,
              0.0287480400836874,
              -0.00692953890353958,
              0.0255522017732768,
              -0.00241182736309177,
              0.0274069450645319,
              -0.00154116829793654,
              -0.0505597663643120,
              -0.0385128209039499,
              0.00263789746425298,
              0.0302058896535371) #n = 40

aMCI_mean <- c(0.0253484335605596,
               -0.000531417777342776,
               0.0228514717180981,
               0.0213552680117668,
               -0.0159751165189477,
               0.0300981207623250,
               0.0102808980833523,
               0.0129752951015591,
               0.0398733899658798,
               -0.00499176514443175,
               -0.0272377523140665,
               0.0254453581879410,
               0.0419129389132819,
               0.00846836859948537,
               -0.0355351865429595,
               -0.00326935827326877,
               0.0410094835966829,
               0.00102045476599876,
               0.0141684852380869,
               0.0220105608030669,
               -0.00137882752642046,
               -0.0689655952600806,
               -0.0243818560245826,
               -0.0714643070826755,
               0.0154235267623480,
               -0.00977623239159323,
               0.0145837908324535,
               0.00449996232519181,
               0.0662819775473663) #n = 29

mMCI_mean <- c(0.0420928947514569,
               0.0122445582540049,
               0.00874645315198615,
               0.0291409731419403,
               0.0118256600432048,
               -0.0177737640784815,
               0.0227339169894717,
               0.0230797184603011,
               -0.00324672404428247,
               0.0110674961727134,
               -0.0218720986860891,
               -0.0138473803990690,
               -0.0630310791399378,
               0.0123588868951274,
               -0.0109943576279138,
               0.194723641982899,
               -0.0280549172674878,
               0.00756000938486777,
               0.0217897052940245) #n = 19

AD_mean <- c(-0.0569352124962124,
             -0.0418788071681668,
             0.0147501515213921,
             -0.0608479956083731,
             -0.0387229085406119,
             0.00273860900285520,
             0.00884224823392592,
             -0.0180716835934308,
             0.00738629790407801,
             -0.0968871597964457) #n = 10

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
  #for C vs. SCD
  all_groups_mean_CvSCD_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 1 | all_groups_mean$Group == 2)
  all_groups_mean_CvSCD_whole_brain$Group <- droplevels(all_groups_mean_CvSCD_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_CvSCD_whole_brain) #this looks like Hedges' g? 
  #for C vs. aMCI
  all_groups_mean_CvaMCI_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 1 | all_groups_mean$Group == 3)
  all_groups_mean_CvaMCI_whole_brain$Group <- droplevels(all_groups_mean_CvaMCI_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_CvaMCI_whole_brain) #this looks like Hedges' g? 
  #for C vs. mMCI
  all_groups_mean_CvmMCI_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 1 | all_groups_mean$Group == 4)
  all_groups_mean_CvmMCI_whole_brain$Group <- droplevels(all_groups_mean_CvmMCI_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_CvmMCI_whole_brain) #this looks like Hedges' g? 
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

#for interaction w/o covariates - 11 edges
#(11,59)
#(17,34)
#(34,59)
#(34,73)
#(35,63)
#(35,64)
#(38,58)
#(46,54)
#(48,83)
#(63,68)
#(68,80)

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
edges_value_vec1 <- unlist(lapply(all_csvs, `[`, 11, 59))
edges_value_vec2 <- unlist(lapply(all_csvs, `[`, 17, 34))
edges_value_vec3 <- unlist(lapply(all_csvs, `[`, 34, 59))
edges_value_vec4 <- unlist(lapply(all_csvs, `[`, 34, 73))
edges_value_vec5 <- unlist(lapply(all_csvs, `[`, 35, 63))
edges_value_vec6 <- unlist(lapply(all_csvs, `[`, 35, 64))
edges_value_vec7 <- unlist(lapply(all_csvs, `[`, 38, 58))
edges_value_vec8 <- unlist(lapply(all_csvs, `[`, 46, 54))
edges_value_vec9 <- unlist(lapply(all_csvs, `[`, 48, 83))
edges_value_vec10 <- unlist(lapply(all_csvs, `[`, 63, 68))
edges_value_vec11 <- unlist(lapply(all_csvs, `[`, 68, 80))

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

#calculate the mean of the values at each edge for each time point
interaction_data$avg_edge_values<-rowMeans(interaction_data[4:14])

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
#for the FPN - 1 sig. edge from the interaction. (edge 2819)
#load in data, copied from matlab
C_FPN_edge <- c(-0.00838000237176770,
                -0.218321481662706,
                -0.287260029845800,
                0.316387643233174,
                -0.303845533214037,
                -0.325228543397035,
                -0.190618827319876,
                -0.165417539492150,
                -0.128140755375018,
                -1.40579355511463,
                -0.00786982357240690,
                -0.165060848542921,
                0.0576153677699300,
                -0.244415304231467,
                -0.548845696639556,
                -0.0620563642409932,
                -0.00189077459643640,
                -0.0606182427278316,
                -0.0846961023181695,
                0.00555055898396240,
                -1.10682364397049) #n = 21

SCD_FPN_edge <- c(0.0729712238058843,
                  0.235614983155697,
                  0.0149894347550110,
                  -0.197342205610026,
                  0.106943893924199,
                  0.00223911440504810,
                  0.0534770330005980,
                  -0.00270813174517826,
                  0.0340429299916920,
                  0.0324765943218540,
                  -0.185004235362829,
                  0.244846479087974,
                  0.0108736441256150,
                  0.120644070777947,
                  0.0442120737964812,
                  -0.000551322375187999,
                  0.120729265989684,
                  0.169397490091050,
                  -0.379862067146365,
                  0.00663546625502450,
                  -0.161660159723486,
                  -0.0722729484461830,
                  -0.00887858268265890,
                  0.0982552082163468,
                  0.200620286025617,
                  0.0693997330268200,
                  -0.0751970272608296,
                  -0.00106752878809170,
                  -0.0306869148652718,
                  0.105566137677252,
                  -0.0334451537481820,
                  -0.0344271851455580,
                  -0.00935635341490912,
                  -0.134236787062559,
                  0.00872578742560132,
                  -0.162273388651984,
                  -0.277471800534747,
                  0.0293307216673037,
                  0.0411618489522420,
                  -0.120389107675352) #n = 40

aMCI_FPN_edge <- c(-0.0657959405565528,
                   0.0306138455654569,
                   -0.0309377031564788,
                   0.0546477954314602,
                   -0.0649042584795920,
                   0.0104063691474202,
                   0.156174773167683,
                   0.118674337508808,
                   -0.198211896758892,
                   -0.0247236212985660,
                   0.509104460934156,
                   -0.00153287446141190,
                   0.0472050672148857,
                   -0.315645287051374,
                   -0.316091511662645,
                   -0.00571695834609900,
                   0.0653789204159391,
                   -0.0108017027583088,
                   -0.0785543116837240,
                   -0.117629173146636,
                   0.0395930898874580,
                   0.141185677979362,
                   -0.00731332086259614,
                   0.0892530632522082,
                   -0.144752699917162,
                   0.213711515467632,
                   -0.0285515246992857,
                   -0.0197031099478029,
                   -0.103269454925941) #n = 29

mMCI_FPN_edge <- c(0.0523355151299210,
                   -0.179153183987316,
                   -0.167569070600578,
                   -0.197554440484089,
                   -0.129159072142808,
                   0.0433401439905820,
                   0.0246849306414600,
                   -0.0303803403689459,
                   0.0466013850062540,
                   0.103373874194204,
                   -0.0513627365680190,
                   0.0882250868163440,
                   -0.00332064024189680,
                   0.238198693153140,
                   -0.0633886450917930,
                   -0.510367257903260,
                   0.0159748598690475,
                   -0.0195958419355160,
                   0.0675433653260371) #n = 19

AD_FPN_edge <- c(0.0886571988963100,
                 -0.0597062918800984,
                 -0.102281773626959,
                 -0.123979724298304,
                 -0.381072645816139,
                 0.131149151056772,
                 -0.136881095616093,
                 -0.214511247859225,
                 -0.0149765701153047,
                 -0.0364929456733247) #n = 10

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
  #for C vs. SCD
  all_groups_FPN_edge_CvSCD_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 1 | all_groups_FPN_edge$Group == 2)
  all_groups_FPN_edge_CvSCD_FPN$Group <- droplevels(all_groups_FPN_edge_CvSCD_FPN$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_CvSCD_FPN) #this looks like Hedges' g? 
  #for C vs. aMCI
  all_groups_FPN_edge_CvaMCI_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 1 | all_groups_FPN_edge$Group == 3)
  all_groups_FPN_edge_CvaMCI_FPN$Group <- droplevels(all_groups_FPN_edge_CvaMCI_FPN$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_CvaMCI_FPN) #this looks like Hedges' g? 
  #for C vs. mMCI
  all_groups_FPN_edge_CvmMCI_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 1 | all_groups_FPN_edge$Group == 4)
  all_groups_FPN_edge_CvmMCI_FPN$Group <- droplevels(all_groups_FPN_edge_CvmMCI_FPN$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_CvmMCI_FPN) #this looks like Hedges' g? 
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






#------------------------------------------------------------------------------#
#whole-brain w/ covariates (sex and age) - 8 significant edges for interaction
#load in data, copied from matlab
C_mean <- c(-0.248371695290698,
            -0.167391155612413,
            -0.0679826418984302,
            -0.0198056366853688,
            -0.115611352890815,
            -0.123950554310607,
            -0.0832782100104459,
            -0.0630505581887735,
            -0.0739250452023715,
            -0.0921418889927240,
            -0.0809108651967829,
            -0.0305563176606932,
            -0.0532870626130680,
            -0.0880887831678522,
            -0.0420683177374441,
            -0.0110820581351833,
            0.0103161855370633,
            0.0504317966682423,
            -0.0728311508852295,
            -0.00188480699029991,
            -0.127991830229605)

SCD_mean <- c(0.0204801856047645,
              0.0419642629105750,
              -0.00164266565446584,
              0.0418397866980180,
              0.0246669405471894,
              0.0584926457794382,
              0.0302905280862231,
              0.0254080559543197,
              0.0420779876277681,
              -0.0131696142286658,
              0.00198838891193341,
              0.0111353333709034,
              0.0264631558588438,
              0.0562730993286726,
              -0.0134110751095646,
              -0.0216885135720253,
              0.0293461205772282,
              0.0145230328819660,
              0.0228918063114099,
              0.0178510898366057,
              -0.00290598962741956,
              0.0280013216303794,
              0.00133638346559459,
              0.0711131925097906,
              -0.000781053119093845,
              0.0159771633202588,
              0.0252394332795549,
              0.0186738359362624,
              0.0649852403390787,
              -0.0206743526966003,
              0.0396010578755422,
              -0.00927226328808191,
              0.0343560959540022,
              -0.00711058824636712,
              0.0383699103706248,
              -0.000569872987429007,
              -0.0734587130558714,
              -0.0562358635638728,
              0.00311742815759446,
              0.0354471039340884)

aMCI_mean <- c(0.0339106965332727,
               -0.00133242390744334,
               0.0312892881640153,
               0.0274538780191741,
               -0.0223048401913647,
               0.0402560352338438,
               0.00856185602318961,
               0.0187154509571718,
               0.0508891778298432,
               -0.00311430676740224,
               -0.0353888896066620,
               0.0350905275312987,
               0.0504371196941543,
               0.0141480885001889,
               -0.0461387614127348,
               -0.00947840452191440,
               0.0536050344483595,
               0.000251764342740135,
               0.0190864945767873,
               0.0291213562295696,
               -0.00159071641172935,
               -0.0962423326978854,
               -0.0365401415721325,
               -0.102242141376960,
               0.0269204247524873,
               -0.0140417218840897,
               0.0174584153035776,
               0.00347770178062307,
               0.0936926707544605)

mMCI_mean <- c(0.0560281672570330,
               0.0118415959382753,
               0.0132446515144881,
               0.0418992845775574,
               0.0142596091397303,
               -0.0273895301844152,
               0.0292154917633872,
               0.0325178037805709,
               -0.00641064337835151,
               0.0146341591109027,
               -0.0305797143335357,
               -0.0191618182113161,
               -0.0864655679406938,
               0.0169731523751445,
               -0.0142690328634120,
               0.269839504235904,
               -0.0350994359896058,
               0.0108696820815041,
               0.0223547246248423)

AD_mean <- c(-0.0758914551943944,
             -0.0575569422794416,
             0.0193528366299086,
             -0.0734624373168590,
             -0.0548247930415497,
             0.00215846660989525,
             0.0151576431682536,
             -0.0234801405877591,
             0.0106969637677793,
             -0.125810318089730)

#whole-brain w/ covariates (sex and age) - 121 significant edges for interaction
#load in data, copied from matlab
# C_mean <- c(-0.613582645619834,
#             -0.697598376983471,
#             -0.588091644603306,
#             -0.294359283289256,
#             -0.677713938008265,
#             -1.74558689823967,
#             -0.587305395181818,
#             -0.795000292628099,
#             -0.385990802644628,
#             -1.12838248734711,
#             -1.01891058488430,
#             -0.992985053305785,
#             -0.705634105157025,
#             -1.22419598836364,
#             -0.178252014008265,
#             -0.543899800938843,
#             0.130737187008264,
#             0.374063609942149,
#             -0.405154994900827,
#             0.169640666604132,
#             -0.496411424190082) #n = 21
# 
# SCD_mean <- c(-0.0938946458595042,
#               0.250653020950413,
#               0.0273065820661157,
#               0.337227241768595,
#               0.458458778057851,
#               0.276189275251322,
#               0.323353418743802,
#               0.280988750016529,
#               0.0856556686611571,
#               0.219210181669422,
#               -0.336146665311653,
#               -0.0855922941074380,
#               -0.161452026289256,
#               0.155313884760331,
#               -0.0913223612727273,
#               -0.0920952220826446,
#               0.271231225242232,
#               0.426897542925620,
#               0.661135668446281,
#               0.297477431876033,
#               0.0572630246446281,
#               0.183927651140496,
#               0.129957493983471,
#               0.0964281041909091,
#               0.0315284849809918,
#               0.151766072132438,
#               0.100810708528926,
#               0.117487054643802,
#               0.364750631462810,
#               0.0414645132314049,
#               0.00899574061652892,
#               0.0967945870214876,
#               0.208897012255372,
#               0.161806984809917,
#               0.783382650479339,
#               0.0262649707438016,
#               0.0138874530578512,
#               0.0861041129752066,
#               -0.220197263644628,
#               -0.121422568132232) #n = 40
# 
# aMCI_mean <- c(0.210628506107438,
#                -0.373920363930579,
#                0.164651322413223,
#                0.0655767509421488,
#                -0.228044966082645,
#                0.163120419190083,
#                7.25656528925613e-05,
#                0.331683450330579,
#                0.167724036942149,
#                0.173174716579835,
#                -0.102794827256198,
#                0.218909892355372,
#                -0.101715906504132,
#                0.206652946776860,
#                0.548375626578512,
#                0.246558109363636,
#                0.165779569933884,
#                0.165851425636364,
#                -0.0198737240413223,
#                0.392255778436364,
#                0.138555096269421,
#                0.343144791842975,
#                -0.136496337016529,
#                0.499187919793388,
#                -0.344153816752066,
#                0.230567536917355,
#                -0.301081897793388,
#                -0.225200317677686,
#                0.0404093502727273) #n = 29
# 
# mMCI_mean <- c(0.383942841074380,
#                0.144847250099174,
#                -0.240137133793388,
#                0.419965928958678,
#                -0.138513236707438,
#                -0.195165875950413,
#                0.196993165776859,
#                0.0959134040247934,
#                0.0115247029421488,
#                0.0643792689719008,
#                -0.174659344760331,
#                -0.0469748975553719,
#                -0.0584582965950413,
#                -0.0209536217297521,
#                0.0225810810520661,
#                0.0487569821404959,
#                -0.0178478506289256,
#                0.0978584039504132,
#                -0.589114572966942) #n = 19
# 
# AD_mean <- c(-0.230920575685950,
#              -0.620028244338843,
#              -0.0999853773966943,
#              -0.760662974074380,
#              -0.898808092016529,
#              -0.307957704057851,
#              -0.285101124033058,
#              -0.0880323645206612,
#              0.154799139454545,
#              -0.877245958851240) #n = 10

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
  #for C vs. SCD
  all_groups_mean_CvSCD_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 1 | all_groups_mean$Group == 2)
  all_groups_mean_CvSCD_whole_brain$Group <- droplevels(all_groups_mean_CvSCD_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_CvSCD_whole_brain) #this looks like Hedges' g? 
  #for C vs. aMCI
  all_groups_mean_CvaMCI_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 1 | all_groups_mean$Group == 3)
  all_groups_mean_CvaMCI_whole_brain$Group <- droplevels(all_groups_mean_CvaMCI_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_CvaMCI_whole_brain) #this looks like Hedges' g? 
  #for C vs. mMCI
  all_groups_mean_CvmMCI_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 1 | all_groups_mean$Group == 4)
  all_groups_mean_CvmMCI_whole_brain$Group <- droplevels(all_groups_mean_CvmMCI_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_CvmMCI_whole_brain) #this looks like Hedges' g? 
  #for SCD vs.AD
  all_groups_mean_SCDvAD_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 2 | all_groups_mean$Group == 5)
  all_groups_mean_SCDvAD_whole_brain$Group <- droplevels(all_groups_mean_SCDvAD_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_SCDvAD_whole_brain) #this looks like Hedges' g? 
  #for aMCI vs.AD
  all_groups_mean_aMCIvAD_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 3 | all_groups_mean$Group == 5)
  all_groups_mean_aMCIvAD_whole_brain$Group <- droplevels(all_groups_mean_aMCIvAD_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_aMCIvAD_whole_brain) #this looks like Hedges' g? 
  #for mMCI vs.AD
  all_groups_mean_mMCIvAD_whole_brain <- subset(all_groups_mean, all_groups_mean$Group == 4 | all_groups_mean$Group == 5)
  all_groups_mean_mMCIvAD_whole_brain$Group <- droplevels(all_groups_mean_mMCIvAD_whole_brain$Group)
  cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_mean_mMCIvAD_whole_brain) #this looks like Hedges' g? 
  
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


#show interaction plot (F0 mean and F2 mean plotted)

#for interaction w/ age & sex covariates - 8 edges
#(17,34)
#(34,73)
#(35,63)
#(35,64)
#(38,58)
#(46,54)
#(48,83)
#(68,80)

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
edges_value_vec1 <- unlist(lapply(all_csvs, `[`, 17, 34))
edges_value_vec2 <- unlist(lapply(all_csvs, `[`, 34, 73))
edges_value_vec3 <- unlist(lapply(all_csvs, `[`, 35, 63))
edges_value_vec4 <- unlist(lapply(all_csvs, `[`, 35, 64))
edges_value_vec5 <- unlist(lapply(all_csvs, `[`, 38, 58))
edges_value_vec6 <- unlist(lapply(all_csvs, `[`, 46, 54))
edges_value_vec7 <- unlist(lapply(all_csvs, `[`, 48, 83))
edges_value_vec8 <- unlist(lapply(all_csvs, `[`, 68, 80))


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

#calculate the mean of the values at each edge for each time point
interaction_data$avg_edge_values<-rowMeans(interaction_data[4:11])

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
#for the FPN w/ covariates age and sex - 3 sig. edge from the interaction. (edge 1945 and 1184 (from contrast 1) and 1632 (from contrast 3)).
#Taking the avg. between these two sig. edges for each participant
#load in data, copied from matlab
C_FPN_edge <- c(-0.341098201175298,
                -0.741772286412839,
                -0.120802339305718,
                -0.268764015358204,
                -0.374934766959177,
                -1.00249796893690,
                -0.652053738026956,
                -0.724602380170321,
                -0.161931363339272,
                -0.568865180704971,
                -0.743761848717864,
                -0.468944417746632,
                -1.62851777740803,
                0.221253422940732,
                -0.574890528584275,
                -0.945255971647818,
                -0.206790703443417,
                0.746385191723509,
                0.808959748695291,
                0.00894736070643971,
                -1.20109611166221) #n = 21

SCD_FPN_edge <- c(-0.578501272576200,
                  0.247329495634578,
                  -0.0784714187258353,
                  0.526729629403462,
                  0.00636884765651218,
                  -0.157387938402222,
                  0.924054336676310,
                  0.555286437815952,
                  -0.0967748567037267,
                  0.0958926476893870,
                  -0.760679825507355,
                  0.498177829466847,
                  -0.444328945732542,
                  0.262958123826054,
                  0.190921635494896,
                  -0.483583613853200,
                  -0.224592513071998,
                  0.531477715119701,
                  0.659693418155873,
                  -0.268409726654045,
                  0.191991384678706,
                  0.0585849105659282,
                  0.522798622830550,
                  0.268714064890265,
                  -0.0189481813178984,
                  0.665571373011644,
                  0.0657102150185516,
                  0.00126099305600517,
                  0.274260770902060,
                  0.380274597783396,
                  0.137804275956817,
                  0.0615525761163274,
                  -0.223762155781438,
                  -0.208886662982315,
                  0.858415796664173,
                  0.638450526030940,
                  0.258344574319903,
                  0.0852551410763216,
                  0.429136072242131,
                  0.0409432285686512) #n = 40

aMCI_FPN_edge <- c(0.138538567465249,
                   -0.486402208729564,
                   0.0819381659764223,
                   -0.499149418092339,
                   -0.317745920040560,
                   0.444499905267739,
                   0.278352343709256,
                   -0.0343449450300685,
                   0.167441454207432,
                   -0.150293552230660,
                   -0.308120611038211,
                   0.0692854373948362,
                   0.158689451361446,
                   0.108547552180537,
                   0.423728500011088,
                   -0.175691368982451,
                   0.173567444481243,
                   0.752996058408102,
                   -0.238364103629937,
                   -0.0794133574009519,
                   0.0527583461949537,
                   0.207267301818602,
                   -0.292934073851607,
                   0.525478420290786,
                   -0.239673331073529,
                   0.285143249745028,
                   -0.143315940874621,
                   0.209705922492762,
                   -0.779373915488239) #n = 29

mMCI_FPN_edge <- c(0.518671543149862,
                   -0.255922913607426,
                   0.203260919753490,
                   0.162246334491434,
                   0.317325770195537,
                   -0.503925463154759,
                   -0.580363101887734,
                   0.121410415692810,
                   0.494266545154690,
                   -0.119261478400639,
                   -0.318164288400657,
                   -0.0427724597971900,
                   0.150627126175260,
                   -0.207995555203672,
                   0.168768648615951,
                   -0.215476302299387,
                   -0.151612752778414,
                   0.0392243242905308,
                   -0.568719238152663) #n = 19

AD_FPN_edge <- c(-0.216066351070047,
                 -0.613535669586493,
                 -0.234643160025299,
                 -0.894476515331026,
                 -0.786157760263108,
                 -0.120735222072723,
                 -0.0313952641893684,
                 0.175925865082613,
                 -0.518688376005600,
                 -0.634583992667660) #n = 10

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
#for C vs. SCD
all_groups_FPN_edge_CvSCD_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 1 | all_groups_FPN_edge$Group == 2)
all_groups_FPN_edge_CvSCD_FPN$Group <- droplevels(all_groups_FPN_edge_CvSCD_FPN$Group)
cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_CvSCD_FPN) #this looks like Hedges' g? 
#for C vs. aMCI
all_groups_FPN_edge_CvaMCI_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 1 | all_groups_FPN_edge$Group == 3)
all_groups_FPN_edge_CvaMCI_FPN$Group <- droplevels(all_groups_FPN_edge_CvaMCI_FPN$Group)
cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_CvaMCI_FPN) #this looks like Hedges' g? 
#for C vs. mMCI
all_groups_FPN_edge_CvmMCI_FPN <- subset(all_groups_FPN_edge, all_groups_FPN_edge$Group == 1 | all_groups_FPN_edge$Group == 4)
all_groups_FPN_edge_CvmMCI_FPN$Group <- droplevels(all_groups_FPN_edge_CvmMCI_FPN$Group)
cohensD(`change_in_FBC(post-pre)`~ Group, data = all_groups_FPN_edge_CvmMCI_FPN) #this looks like Hedges' g? 
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
  ylab("Difference scores in structural connectivity in the FPN") +
  scale_x_discrete(labels = c("1" = "Control", "2" = "SCD", "3" = "aMCI", "4" = "mMCI", "5" = "AD")) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=18))+
  theme(axis.title=element_text(size=24))#+
  #coord_flip()











