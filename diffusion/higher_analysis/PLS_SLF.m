%This script will run partial least squares (PLS) analysis on the data. The
%data is comprised of the superior longitudinal fasciculus (SLF) and the 
%neuropsychological assessment from the DPRC cohort. Imputation has also 
%been performed on the data through the Missing Data Imputation (MDI) 
%Toolbox for Matlab. PLS is run through a toolbox on Matlab, which has been 
%downloaded from Rotman Research Institue:
%(https://www.rotman-baycrest.on.ca/index.php?section=84)

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 20/05/22

clc;
clear all;
close all;

%1. Load data, choose directory, define variables, any add necessary pathways 

%cd C:\Users\ltah262\Desktop\PLS
%load SLF_data.mat; 

%excelFile = input('Please name the excel file you wish to use: ', 's');
%fileLocation = input(['Where is this file, ' excelFile  ', located? Name the directory: '], 's');

%shortcut: 
%excelFile = 'DPRC_neuropsych_data.xlsx';
%fileLocation = 'H:\ltah262\PhD\ExecutiveFunction\NeuroPsychAssessment\data\PLS\MDI';
%cd ([fileLocation]);

%read in excel file into matlab format (n = 229)
%full_data = readtable(excelFile);
%[num, txt, raw] = xlsread(excelFile, 'A:I');

%remove participants who are missing too much data (n = 211)
%edited_data = full_data([1:3,5:47,49:55,57:65,67:118,120:148,150:156,158:161,172:211,213:229],:);

%removed 3 participants: 
%sub-ADPRC0004F2
%sub-ADPRC0063F0
%sub-ADPRC0075F0 %ok now
%sub-ADPRC0088F0 %ok now
%sub-ADPRC0163F0
%sub-ADPRC0210F0 %ok now
%sub-ADPRC0222F0 %ok now
%sub-ADPRC0234F0 %ok now
%sub-ADPRC0235F0 %ok now
%sub-ADPRC0238F0 %ok now
%sub-ADPRC0240F0 %ok now
%sub-ADPRC0241F0 %ok now
%sub-ADPRC0243F0 %ok now
%sub-ADPRC0244F0 %ok now
%sub-ADPRC0248F0 %ok now
%sub-ADPRC0249F0 %ok now
%sub-ADPRCPL01F0  %ok now
%sub-DDPRC0013F0  %ok now

%for longitudinal dataset (n = 119), removed: 
%sub-ADPRC0091F2 (missing all exec. func data)
%sub-ADPRC0153F2 (missing all exec. func data)
%sub-CDPRC0019F2 (missing all exec. func data)
%sub-CDPRC0023F2 (missing all exec. func data)
%sub-DDPRC0012F2 (missing all exec. func data)
%note: keeping sub-ADPRC0163F0 (even though more than 3 neuropsych tests missing)








%names of behavioural data
%all 11 neuropsych scores
behav_11names = {'TMT-A','Colour Naming','Word Reading','Haytime1','TMT-B','Inhibition','Category Switching','HayTime2','HayTotError','Letter Fluency','Category Fluency'};
%all 12 neuropsych scores
%behav_names = {'TMT-A','TMT-B','Colour Naming','Word Reading','Inhibition','Letter Fluency','Category Fluency','Category Switching','Haytime1','HayTime2','HayCatAError','HayCatBError'};
% 4 processing speed tests
%ProcSpeed_names = {'TMT-A','Colour Naming','Word Reading','Haytime1'};
%3 processing speed names
ProcSpeed_3names = {'TMT-A','Colour Naming','Word Reading'};
ProcSpeed_name = {'Processing Speed'};
% 4 Inhibition test 
%Inhibition_names = {'TMT-B','Inhibition','Category Switching','Haytime2'};
% 2 generation tests
%Generation_names = {'Letter Fluency','Category Fluency'};
%3 Inhib Rate scores
%InhibRate_names = {'TMT-B/TMT-A','Inhibition/ColourNaming','Inhibition/WordReading'};
% 7 inhibition scores (+ Inhib Rate scores)
%Inhibition_7names = {'TMT-B','Inhibition','Category Switching','Haytime2','TMT-B/TMT-A','Inhibition/ColourNaming','Inhibition/WordReading'};
%5 inhib names
%inhib_5names = {'TMT-B/TMT-A','Inhibition/ColourNaming','Category Switching','HayTime2-HayTime1','HayTotalError'};
inhib_5names = {'TMT-B/TMT-A','Interference Effect','Category Switching','HayTime2-HayTime1','HayTotalError'};
inhib_name = {'Inhibition'};
%13 neuropsych scores (omit HayError scores, and include Inhib Rate scores)
%behav_13names = {'TMT-A','TMT-B','Colour Naming','Word Reading','Inhibition','Letter Fluency','Category Fluency','Category Switching','Haytime1','HayTime2','TMT-B/TMT-A','Inhibition/ColourNaming','Inhibition/WordReading'};
%13 neuropsych scores (omit HayError scores, and include Inhib Rate
%scores)- reordered Processing speed vs. Inhibition vs. inhibition
%behav_reordered_13names = {'TMT-A','Colour Naming','Word Reading','Haytime1','TMT-B','Inhibition','Category Switching','HayTime2','TMT-B/TMT-A','Inhibition/ColourNaming','Inhibition/WordReading','Letter Fluency','Category Fluency'};


%add in pathways for PLS
addpath 'C:/Program Files/MATLAB/PLS/plscmd'
addpath  'C:/Program Files/MATLAB/PLS/plsgui'
%add in pathways for MDI
addpath 'C:/Program Files/MATLAB/MDI'


%2. Perform imputation on data
fileLocation = 'H:\ltah262\PhD\ExecutiveFunction\NeuroPsychAssessment\data\PLS\MDI';
cd ([fileLocation]);
excelFile_all_behav = 'DPRC_neuropsych_data_only_switchedVars(n=226).xlsx';
%or, for connectome data (n = 224):
excelFile_all_behav = 'DPRC_neuropsych_data_only_switchedVars(n=223).xlsx';
%or, for data w/o AD group (n = 201)
excelFile_all_behav = 'DPRC_neuropsych_data_only_switchedVars(n=201).xlsx';


excelFile_all_behav = 'DPRC_neuropsych_data_F0_only_switched_Vars(n=119).xlsx';


%read in excel file into matlab format (n = 211)
all_behav_data = readtable(excelFile_all_behav);

raw_scores = table2array(all_behav_data(:,[1,3,5,7,9,11,13,15,17,19,21,23,25,26,27,31,33,39]));
zscores = table2array(all_behav_data(:,[2,4,6,8,10,12,14,16,18,20,22,24,28,29,30,32,34,40]));

%raw_scores_noHayerror_rateInhibScores = table2array(all_behav_data(:,[1,3,5,7,9,11,13,15,17,19,25,26,27]));


%mditoolbox(raw_scores, 'TSR', imputed_raw_scores2);
MDIgui

%imputation done for 16 participants:
%sub-ADPRC0001F0
%sub-ADPRC0002F0
%sub-ADPRC0003F0
%sub-ADPRC0005F0
%sub-ADPRC0007F0
%sub-ADPRC0040F0
%sub-ADPRC0047F0
%sub-ADPRC0075F0 
%sub-ADPRC0088F0 
%sub-ADPRC0101F0
%sub-ADPRCPL01F0
%sub-ADPRCPL02F0
%sub-CDPRC0013F0
%sub-CDPRC0041F0
%sub-DDPRC0001F0
%sub-DDPRC0013F0



%3. Run PLS on the imputted data

%combine the imputted neuropsych data with the SLF data
cd H:\ltah262\PhD\ExecutiveFunction\NeuroPsychAssessment\data\PLS\MDI;
all_data = readtable('SLF_data_tracts_of_interest(n=201).xlsx');

%extract SLF data from the imputted files:
FD_data = table2array(all_data(:,{'mn_FD_SLF1_L', 'mn_FD_SLF2_L', 'mn_FD_SLF3_L', 'mn_FD_SLF1_R', 'mn_FD_SLF2_R', 'mn_FD_SLF3_R'}));
FC_data = table2array(all_data(:,{'mn_FC_SLF1_L', 'mn_FC_SLF2_L', 'mn_FC_SLF3_L', 'mn_FC_SLF1_R', 'mn_FC_SLF2_R', 'mn_FC_SLF3_R'}));
FDC_data = table2array(all_data(:,{'mn_FDC_SLF1_L', 'mn_FDC_SLF2_L', 'mn_FDC_SLF3_L', 'mn_FDC_SLF1_R', 'mn_FDC_SLF2_R', 'mn_FDC_SLF3_R'}));

%FPN structural connectome data
SC_FPN_data = readtable('SC_small_FPN_data(n=224).xlsx'); %(for SC analysis, there are 227 available participants, so 227-3 = 224 total)
SC_FPN_data = table2array(SC_FPN_data(:,2:11));
%all SC FPN nodes (82 nodes; 3,321 edges)
SC_FPN_data = readtable('SC_big_FPN_data(n=224).xlsx'); %(for SC analysis, there are 227 available participants, so 227-3 = 224 total)
SC_FPN_data = table2array(SC_FPN_data(:,2:3322));

FC_FPN_data = readtable('FC_small_FPN_C2pre_data(n=119).xlsx'); %(for FC analysis, there are 119 available participants total)
FC_FPN_data = table2array(FC_FPN_data(:,2:11));
%all FC FPN nodes (82 nodes; 3,321 edges)
FC_FPN_data = readtable('FC_big_FPN_data(n=223).xlsx'); %(for FC analysis, there are 119 available participants total)
FC_FPN_data = table2array(FC_FPN_data(:,2:3322));


%choose data from specific cognitive domains
all_11_behav_raw_scores = MDItoolbox_results.X_imputed(:,[1,3,4,9,2,5,8,10,16,6,7]);
%%inhib_behav5 = MDItoolbox_results.X_imputed(:,[15,13,8,17,16]);
ProcSpeed_behav3 = MDItoolbox_results.X_imputed(:,[1,3,4]);
%take average of ProcSpeed z-scores for one value
ProcSpeed_means = mean(ProcSpeed_behav3, 2);
%all_12_behav_raw_scores = MDItoolbox_results.X_imputed(:,[1:12]);
%inhibition_raw_scores = MDItoolbox_results.X_imputed(:,[2,5,8,10]); 
%proc_speed_raw_scores = MDItoolbox_results.X_imputed(:,[1,3,4,9]);
%InhibRate_raw_scores = MDItoolbox_results.X_imputed(:,[13:15]); 
%inhibition_InhibRate_raw_scores = MDItoolbox_results.X_imputed(:,[2,5,8,10,13,14,15]); 
%behav_13_reordered_raw_scores = MDItoolbox_results.X_imputed(:,[1,3,4,9,2,5,8,10,13,14,15,6,7]);
inhib_behav5 = MDItoolbox_results.X_imputed(:,[15,18,8,17,16]);
%take average of inhib z-scores for one value
inhib_means = mean(inhib_behav5, 2);



%contrast definitions:
contrast_1 = [1];
contrast_3 = [1;1;1];
contrast_5 = [1;1;1;1;1];

%zscores
%all_12_behav_zscores = MDItoolbox_results.X_imputed(:,[1:12]);
%inhibition_zscores = MDItoolbox_results.X_imputed(:,[2,5,8,10]); 
%proc_speed_zscores = MDItoolbox_results.X_imputed(:,[1,3,4,9]);
%no Hayling Error scores and Inhib rate scores includes
%behav13_noHayError_InhibRate_raw_scores = MDItoolbox_results.X_imputed(:,[1:10,13:15]);
%inhibition_InhibRate_zscores = MDItoolbox_results.X_imputed(:,[2,5,8,10,13,14,15]); 
%InhibRate_zscores = MDItoolbox_results.X_imputed(:,[13:15]); 
%behav_13_reordered_zscores = MDItoolbox_results.X_imputed(:,[1,3,4,9,2,5,8,10,13,14,15,6,7]);
%all_11_behav_zscores = MDItoolbox_results.X_imputed(:,[1,3,4,9,2,5,8,10,16,6,7]);
%inhib_behav5_zscores = MDItoolbox_results.X_imputed(:,[15,13,8,17,16]);
%ProcSpeed_behav3_zscores = MDItoolbox_results.X_imputed(:,[1,3,4]);

%age residualize data (behavioural + SLF data):
fileLocation = 'H:\ltah262\PhD\ExecutiveFunction\NeuroPsychAssessment\data\PLS\MDI';
cd ([fileLocation]);

%function format: residualize(xvecs,datamat)
age_data = table2array(readtable('age_data_(n=226).xlsx'));
age_data = table2array(readtable('age_data_(n=201).xlsx'));
sex_data = table2array(readtable('sex_data_(n=226).xlsx'));

two_covars = table2array(readtable('two_covariates_(n=226).xlsx'));
two_covars = table2array(readtable('two_covariates_(n=201).xlsx'));
two_covars = table2array(readtable('two_covariates_F2(n=114).xlsx'));

age_data = table2array(readtable('age_F2_data(n=114).xlsx'));


resid_behav_data = residualize(age_data, inhib_behav5);
resid_FD_data = residualize(two_covars, FD_data);
resid_FC_data = residualize(two_covars, FC_data);
resid_FDC_data = residualize(two_covars, FDC_data);

resid_SC_FPN_data = residualize(two_covars, SC_FPN_data);
resid_FC_FPN_data = residualize(two_covars, FC_FPN_data);


%set up PLS options
datamat_list{1} = cat(1, resid_FD_data);
num_subj_lst = height(resid_FD_data);
num_cond = 1;
option.stacked_behavdata = inhib_behav5; 
option.method = 5; %non-rotated
option.num_perm = 5000;
option.num_boot = 5000;
option.clim = 95;
option.stacked_designdata = contrast_5; 
resultRotated = pls_analysis(datamat_list, num_subj_lst, num_cond, option); 

%save('FD_behav_PLS_raw_scores.mat', 'resultRotated'); 

%4. Make PLS figures
figure 
bar((resultRotated.boot_result.orig_corr)*-1); 
hold on 
for i = 1:length(resultRotated.boot_result.orig_corr) %what is this - bootstrap correlation value?
    plot([i i], [(resultRotated.boot_result.ulcorr(i)*-1) (resultRotated.boot_result.llcorr(i)*-1)],'-k', 'LineWidth',2); %95 percent CI - upper and lower bound
end
xticklabels(inhib_5names); 
xtickangle(45); 
xlabel('Neuropsychological Assessment');
ylabel('Bootstrap Ratio Correlation Value');
box off 
set(gcf,'color','w'); 
title(['LV p = ' num2str(resultRotated.perm_result.sprob)]);





%for Group contrast PLS: 
%define group names: 
Group_names = {'C','SCD','aMCI','mMCI','AD'};
Group_names = {'C','SCD','aMCI','mMCI'};


%For fibre metric data:
% cd 'H:\ltah262\PhD\ExecutiveFunction\NeuroPsychAssessment\data\PLS\MDI\group_contrast\FibreMetrics';
% %for FD:
% FD_C = table2array(readtable('C_FD.xlsx'));
% FD_SCD = table2array(readtable('SCD_FD.xlsx'));
% FD_aMCI = table2array(readtable('aMCI_FD.xlsx'));
% FD_mMCI = table2array(readtable('mMCI_FD.xlsx'));
% FD_AD = table2array(readtable('AD_FD.xlsx'));
% %for FC:
% FC_C = table2array(readtable('C_FC.xlsx'));
% FC_SCD = table2array(readtable('SCD_FC.xlsx'));
% FC_aMCI = table2array(readtable('aMCI_FC.xlsx'));
% FC_mMCI = table2array(readtable('mMCI_FC.xlsx'));
% FC_AD = table2array(readtable('AD_FC.xlsx'));
% %for FDC:
% FDC_C = table2array(readtable('C_FDC.xlsx'));
% FDC_SCD = table2array(readtable('SCD_FDC.xlsx'));
% FDC_aMCI = table2array(readtable('aMCI_FDC.xlsx'));
% FDC_mMCI = table2array(readtable('mMCI_FDC.xlsx'));
% FDC_AD = table2array(readtable('AD_FDC.xlsx'));



%For behavioural measures (executive function neuropsych data):
cd 'H:\ltah262\PhD\ExecutiveFunction\NeuroPsychAssessment\data\PLS\MDI\cross-sectional\dMRI\group_contrast\Behavioural-ExecFunc';
%for all raw variables (11 variables):
all_behav11_group_ordered = table2array(readtable('all_behav11(n=226).xlsx'));
inhib_behav5_group_ordered = table2array(readtable('inhib_behav5_zscores(n=201).xlsx'));
ProcSpeed_behav3_group_ordered = table2array(readtable('ProcSpeed_behav3_zscores(n=201).xlsx'));

%normalise the data
%age residualize data (behavioural + SLF data):
fileLocation = 'H:\ltah262\PhD\ExecutiveFunction\NeuroPsychAssessment\data\PLS\MDI';
cd ([fileLocation]);
%load in age + SLF and/or multiple covariates (ordered by group):
age_data = table2array(readtable('age_ordered_by_group(n=201).xlsx'));
SLF_data_ordered = table2array(readtable('SLF_ordered_by_group(n=201).xlsx'));
two_covars_ordered = table2array(readtable('two_covariates_ordered_by_group(n=201).xlsx'));


%non-residualised SLF data: 
FD_C = SLF_data_ordered(1:35,1:6);
FD_SCD = SLF_data_ordered(36:94,1:6);
FD_aMCI = SLF_data_ordered(95:149,1:6);
FD_mMCI = SLF_data_ordered(150:201,1:6);
FD_AD = SLF_data_ordered(202:226,1:6);

FC_C = SLF_data_ordered(1:35,7:12);
FC_SCD = SLF_data_ordered(36:94,7:12);
FC_aMCI = SLF_data_ordered(95:149,7:12);
FC_mMCI = SLF_data_ordered(150:201,7:12);
FC_AD = SLF_data_ordered(202:226,7:12);

FDC_C = SLF_data_ordered(1:35,13:18);
FDC_SCD = SLF_data_ordered(36:94,13:18);
FDC_aMCI = SLF_data_ordered(95:149,13:18);
FDC_mMCI = SLF_data_ordered(150:201,13:18);
FDC_AD = SLF_data_ordered(202:226,13:18);

%residualize by age - 
%function format: residualize(xvecs,datamat)
resid_all_behav11 = residualize(age_data, all_behav11_group_ordered);
resid_inhib_behav5_data = residualize(age_data, inhib_behav5_group_ordered);
resid_ProcSpeed_behav3_data = residualize(age_data, ProcSpeed_behav3_group_ordered);
resid_SLF_data = residualize(two_covars_ordered, SLF_data_ordered);
%resid_inhib_behav5_data = residualize(age_data, inhib_behav5_2);

resid_FD_C = resid_SLF_data(1:35,1:6);
resid_FD_SCD = resid_SLF_data(36:94,1:6);
resid_FD_aMCI = resid_SLF_data(95:149,1:6);
resid_FD_mMCI = resid_SLF_data(150:201,1:6);
resid_FD_AD = resid_SLF_data(202:226,1:6);

resid_FC_C = resid_SLF_data(1:35,7:12);
resid_FC_SCD = resid_SLF_data(36:94,7:12);
resid_FC_aMCI = resid_SLF_data(95:149,7:12);
resid_FC_mMCI = resid_SLF_data(150:201,7:12);
resid_FC_AD = resid_SLF_data(202:226,7:12);

resid_FDC_C = resid_SLF_data(1:35,13:18);
resid_FDC_SCD = resid_SLF_data(36:94,13:18);
resid_FDC_aMCI = resid_SLF_data(95:149,13:18);
resid_FDC_mMCI = resid_SLF_data(150:201,13:18);
resid_FDC_AD = resid_SLF_data(202:226,13:18);



% %for FPN SC:
% %For behavioural measures (executive function neuropsych data):
% cd 'H:\ltah262\PhD\ExecutiveFunction\NeuroPsychAssessment\data\PLS\MDI\group_contrast\Behavioural-ExecFunc';
% %for all raw variables (11 variables):
% all_behav11_group_ordered = table2array(readtable('all_behav11(n=224).xlsx'));
% inhib_behav5_group_ordered = table2array(readtable('inhib_behav5(n=224).xlsx'));
% ProcSpeed_behav3_group_ordered = table2array(readtable('ProcSpeed_behav3(n=224).xlsx'));
% 
% fileLocation = 'H:\ltah262\PhD\ExecutiveFunction\NeuroPsychAssessment\data\PLS\MDI\group_contrast\FBC';
% cd ([fileLocation]);
% %load in age + FBC (ordered by group):
% age_data = table2array(readtable('age_data_ordered_by_group(n=224).xlsx'));
% SC_FPN_data_ordered = table2array(readtable('SC_FPN_data_ordered_by_group(n=224).xlsx'));
% 
% %non-residualised SC_FPN data: 
% SC_FPN_C = SC_FPN_data_ordered(1:35,1:10);
% SC_FPN_SCD = SC_FPN_data_ordered(36:94,1:10);
% SC_FPN_aMCI = SC_FPN_data_ordered(95:148,1:10);
% SC_FPN_mMCI = SC_FPN_data_ordered(149:200,1:10);
% SC_FPN_AD = SC_FPN_data_ordered(201:224,1:10);
% 
% %residualize by age - 
% %function format: residualize(xvecs,datamat)
% resid_all_behav11 = residualize(age_data, all_behav11_group_ordered);
% resid_inhib_behav5_data = residualize(age_data, inhib_behav5_group_ordered);
% resid_ProcSpeed_behav3_data = residualize(age_data, ProcSpeed_behav3_group_ordered);
% resid_SC_FPN_data = residualize(age_data, SC_FPN_data_ordered);
% %resid_inhib_behav5_data = residualize(age_data, inhib_behav5_2);
% 
% resid_SC_FPN_C = resid_SC_FPN_data(1:35,1:10);
% resid_SC_FPN_SCD = resid_SC_FPN_data(36:94,1:10);
% resid_SC_FPN_aMCI = resid_SC_FPN_data(95:148,1:10);
% resid_SC_FPN_mMCI = resid_SC_FPN_data(149:200,1:10);
% resid_SC_FPN_AD = resid_SC_FPN_data(201:224,1:10);
% 
%for SC FPN
inhib_behav5_group_ordered = table2array(readtable('inhib_behav5_orderedF0(n=119).xlsx'));%F0
proc_speed3_group_ordered = table2array(readtable('proc_speed_behav3_orderedF0(n=119).xlsx'));
age_data = table2array(readtable('age_ordered_by_groupF0(n=119).xlsx'));
two_covars_ordered = table2array(readtable('two_covariates_ordered_by_groupF0(n=119).xlsx'));

inhib_behav5_group_ordered = table2array(readtable('inhib_behav5_orderedF2(n=114).xlsx'));%F2
proc_speed3_group_ordered = table2array(readtable('proc_speed_behav3_orderedF2(n=114).xlsx'));
age_data = table2array(readtable('age_ordered_by_groupF2(n=114).xlsx'));
two_covars_ordered = table2array(readtable('two_covariates_ordered_by_groupF2(n=114).xlsx'));

inhib_behav5_group_ordered = table2array(readtable('inhib_behav5_ordered_post-pre(diff)(n=114).xlsx'));%diff (post-pre)
proc_speed3_group_ordered = table2array(readtable('proc_speed_behav3_ordered_post-pre(diff)(n=114).xlsx'));
age_data = table2array(readtable('age_ordered_by_groupF2(n=114).xlsx'));
two_covars_ordered = table2array(readtable('two_covariates_ordered_by_groupF2(n=114).xlsx'));

%big FPN -
SC_big_FPN_dataF0_ordered = table2array(readtable('SC_big_FPN_data_ordered_by_groupF0(n=119).xlsx')); %F0
SC_big_FPN_dataF2_ordered = table2array(readtable('SC_big_FPN_data_ordered_by_groupF2(n=114).xlsx')); %F2
SC_big_FPN_data_post_prediff_ordered = table2array(readtable('SC_big_FPN_data_ordered_by_group_post-prediff(n=114).xlsx')); %diff (post-pre)

SC_big_FPN_data_ordered = table2array(readtable('SC_big_FPN_data_ordered_by_group(n=224).xlsx')); %F0, n = 224
inhib_behav5_group_ordered = table2array(readtable('inhib_behav5_ordered(n=224).xlsx'));%F0
proc_speed3_group_ordered = table2array(readtable('proc_speed_behav3_ordered(n=224).xlsx'));
age_data = table2array(readtable('age_data_ordered_by_group(n=224).xlsx'));
two_covars_ordered = table2array(readtable('two_covariates_ordered_by_group(n=224).xlsx'));

%non-residualised SC_big_FPN data: 
SC_FPN_C = SC_big_FPN_dataF0_ordered(1:21,1:3321); %F0
SC_FPN_SCD = SC_big_FPN_dataF0_ordered(22:61,1:3321);
SC_FPN_aMCI = SC_big_FPN_dataF0_ordered(62:90,1:3321);
SC_FPN_mMCI = SC_big_FPN_dataF0_ordered(91:109,1:3321);
SC_FPN_AD = SC_big_FPN_dataF0_ordered(110:119,1:3321);

SC_FPN_C = SC_big_FPN_dataF2_ordered(1:21,1:3321); %F2
SC_FPN_SCD = SC_big_FPN_dataF2_ordered(22:60,1:3321);
SC_FPN_aMCI = SC_big_FPN_dataF2_ordered(61:88,1:3321);
SC_FPN_mMCI = SC_big_FPN_dataF2_ordered(89:106,1:3321);
SC_FPN_AD = SC_big_FPN_dataF2_ordered(107:114,1:3321);

SC_FPN_C = SC_big_FPN_data_post_prediff_ordered(1:21,1:3321); %diff (post-pre)
SC_FPN_SCD = SC_big_FPN_data_post_prediff_ordered(22:60,1:3321);
SC_FPN_aMCI = SC_big_FPN_data_post_prediff_ordered(61:88,1:3321);
SC_FPN_mMCI = SC_big_FPN_data_post_prediff_ordered(89:106,1:3321);
SC_FPN_AD = SC_big_FPN_data_post_prediff_ordered(107:114,1:3321);


SC_FPN_C = SC_big_FPN_data_ordered(1:35,1:3321); %F0 (n = 224)
SC_FPN_SCD = SC_big_FPN_data_ordered(36:94,1:3321);
SC_FPN_aMCI = SC_big_FPN_data_ordered(95:148,1:3321);
SC_FPN_mMCI = SC_big_FPN_data_ordered(149:200,1:3321);
SC_FPN_AD = SC_big_FPN_data_ordered(201:224,1:3321);


%residualise behav data by age
resid_inhib_behav5_data = residualize(age_data, inhib_behav5_group_ordered);
resid_proc_speed_behav3_data = residualize(age_data, proc_speed3_group_ordered);

%residualise brain data by age and sex: 
resid_SC_FPN_data = residualize(two_covars_ordered, SC_big_FPN_dataF0_ordered); %F0
resid_SC_FPN_C = resid_SC_FPN_data(1:21,1:3321);
resid_SC_FPN_SCD = resid_SC_FPN_data(22:61,1:3321);
resid_SC_FPN_aMCI = resid_SC_FPN_data(62:90,1:3321);
resid_SC_FPN_mMCI = resid_SC_FPN_data(91:109,1:3321);
resid_SC_FPN_AD = resid_SC_FPN_data(110:119,1:3321);

resid_SC_FPN_data = residualize(two_covars_ordered, SC_big_FPN_dataF2_ordered); %F2
resid_SC_FPN_C = resid_SC_FPN_data(1:21,1:3321);
resid_SC_FPN_SCD = resid_SC_FPN_data(22:60,1:3321);
resid_SC_FPN_aMCI = resid_SC_FPN_data(61:88,1:3321);
resid_SC_FPN_mMCI = resid_SC_FPN_data(89:106,1:3321);
resid_SC_FPN_AD = resid_SC_FPN_data(107:114,1:3321);

resid_SC_FPN_data = residualize(two_covars_ordered, SC_big_FPN_data_post_prediff_ordered); %post-pre diff
resid_SC_FPN_C = resid_SC_FPN_data(1:21,1:3321);
resid_SC_FPN_SCD = resid_SC_FPN_data(22:60,1:3321);
resid_SC_FPN_aMCI = resid_SC_FPN_data(61:88,1:3321);
resid_SC_FPN_mMCI = resid_SC_FPN_data(89:106,1:3321);
resid_FC_FPN_AD = resid_SC_FPN_data(107:114,1:3321);

resid_SC_big_FPN_data_ordered = residualize(two_covars_ordered, SC_big_FPN_data_ordered); %F0 (n = 224)
resid_SC_FPN_C = resid_SC_big_FPN_data_ordered(1:35,1:3321); %F0 (n = 224)
resid_SC_FPN_SCD = resid_SC_big_FPN_data_ordered(36:94,1:3321);
resid_SC_FPN_aMCI = resid_SC_big_FPN_data_ordered(95:148,1:3321);
resid_SC_FPN_mMCI = resid_SC_big_FPN_data_ordered(149:200,1:3321);
resid_SC_FPN_AD = resid_SC_big_FPN_data_ordered(201:224,1:3321);

%small FPN
SC_small_FPN_dataF0_ordered = table2array(readtable('SC_small_FPN_data_ordered_by_groupF0(n=119).xlsx')); %F0 
SC_small_FPN_dataF2_ordered = table2array(readtable('SC_small_FPN_data_ordered_by_groupF2(n=114).xlsx')); %F2
SC_small_FPN_data_post_prediff_ordered = table2array(readtable('SC_small_FPN_data_ordered_by_group_post-prediff(n=114).xlsx')); %diff (post-pre)

%non-residualised FC_small_FPN data: 
SC_FPN_C = SC_small_FPN_dataF0_ordered(1:21,1:10); %F0
SC_FPN_SCD = SC_small_FPN_dataF0_ordered(22:61,1:10);
SC_FPN_aMCI = SC_small_FPN_dataF0_ordered(62:90,1:10);
SC_FPN_mMCI = SC_small_FPN_dataF0_ordered(91:109,1:10);
SC_FPN_AD = SC_small_FPN_dataF0_ordered(110:119,1:10);

SC_FPN_C = SC_small_FPN_dataF2_ordered(1:21,1:10); %F2
SC_FPN_SCD = SC_small_FPN_dataF2_ordered(22:60,1:10);
SC_FPN_aMCI = SC_small_FPN_dataF2_ordered(61:88,1:10);
SC_FPN_mMCI = SC_small_FPN_dataF2_ordered(89:106,1:10);
SC_FPN_AD = SC_small_FPN_dataF2_ordered(107:114,1:10);

SC_FPN_C = SC_small_FPN_data_post_prediff_ordered(1:21,1:10); %post-pre diff
SC_FPN_SCD = SC_small_FPN_data_post_prediff_ordered(22:60,1:10);
SC_FPN_aMCI = SC_small_FPN_data_post_prediff_ordered(61:88,1:10);
SC_FPN_mMCI = SC_small_FPN_data_post_prediff_ordered(89:106,1:10);
SC_FPN_AD = SC_small_FPN_data_post_prediff_ordered(107:114,1:10);

%residualise brain data by age and sex: 
resid_SC_FPN_data = residualize(two_covars_ordered, SC_small_FPN_dataF0_ordered); %F0
resid_SC_FPN_C = resid_SC_FPN_data(1:21,1:10);
resid_SC_FPN_SCD = resid_SC_FPN_data(22:61,1:10);
resid_SC_FPN_aMCI = resid_SC_FPN_data(62:90,1:10);
resid_SC_FPN_mMCI = resid_SC_FPN_data(91:109,1:10);
resid_SC_FPN_AD = resid_SC_FPN_data(110:119,1:10);

resid_SC_FPN_data = residualize(two_covars_ordered, SC_small_FPN_dataF2_ordered); %F2
resid_SC_FPN_C = resid_SC_FPN_data(1:21,1:10);
resid_SC_FPN_SCD = resid_SC_FPN_data(22:60,1:10);
resid_SC_FPN_aMCI = resid_SC_FPN_data(61:88,1:10);
resid_SC_FPN_mMCI = resid_SC_FPN_data(89:106,1:10);
resid_SC_FPN_AD = resid_SC_FPN_data(107:114,1:10);

resid_SC_FPN_data = residualize(two_covars_ordered, SC_small_FPN_data_post_prediff_ordered); %post-pre diff
resid_SC_FPN_C = resid_SC_FPN_data(1:21,1:10);
resid_SC_FPN_SCD = resid_SC_FPN_data(22:60,1:10);
resid_SC_FPN_aMCI = resid_SC_FPN_data(61:88,1:10);
resid_SC_FPN_mMCI = resid_SC_FPN_data(89:106,1:10);
resid_SC_FPN_AD = resid_SC_FPN_data(107:114,1:10);




%for fMRI FPN:
inhib_behav5_group_ordered = table2array(readtable('inhib_behav5_z-scores_ordered(n=223).xlsx'));%F0
proc_speed3_group_ordered = table2array(readtable('proc_speed_behav3_z-scores_ordered(n=223).xlsx'));
age_data = table2array(readtable('age_ordered_by_group(n=223).xlsx'));
two_covars_ordered = table2array(readtable('two_covariates_ordered_by_group(n=223).xlsx'));

inhib_behav5_group_ordered = table2array(readtable('inhib_behav5_orderedF2(n=114).xlsx'));%F2
proc_speed3_group_ordered = table2array(readtable('proc_speed_behav3_orderedF2(n=114).xlsx'));
age_data = table2array(readtable('age_ordered_by_groupF2(n=114).xlsx'));
two_covars_ordered = table2array(readtable('two_covariates_ordered_by_groupF2(n=114).xlsx'));

inhib_behav5_group_ordered = table2array(readtable('inhib_behav5_ordered_post-pre(diff)(n=114).xlsx'));%diff (post-pre)
proc_speed3_group_ordered = table2array(readtable('proc_speed_behav3_ordered_post-pre(diff)(n=114).xlsx'));
age_data = table2array(readtable('age_ordered_by_groupF2(n=114).xlsx'));
two_covars_ordered = table2array(readtable('two_covariates_ordered_by_groupF2(n=114).xlsx'));

%big FPN -
FC_big_FPN_dataF0_ordered = table2array(readtable('FC_big_FPN_data_ordered_by_groupF0(n=119).xlsx')); %F0
FC_big_FPN_dataF2_ordered = table2array(readtable('FC_big_FPN_data_ordered_by_groupF2(n=114).xlsx')); %F2
FC_big_FPN_data_post_prediff_ordered = table2array(readtable('FC_big_FPN_data_ordered_by_group_post-prediff(n=114).xlsx')); %diff (post-pre)

FC_big_FPN_data_ordered = table2array(readtable('FC_big_FPN_data_ordered_by_group(n=223).xlsx')); %F0

%non-residualised FC_big_FPN data: 
FC_FPN_C = FC_big_FPN_dataF0_ordered(1:21,1:3321); %F0
FC_FPN_SCD = FC_big_FPN_dataF0_ordered(22:61,1:3321);
FC_FPN_aMCI = FC_big_FPN_dataF0_ordered(62:90,1:3321);
FC_FPN_mMCI = FC_big_FPN_dataF0_ordered(91:109,1:3321);
FC_FPN_AD = FC_big_FPN_dataF0_ordered(110:119,1:3321);

FC_FPN_C = FC_big_FPN_dataF2_ordered(1:21,1:3321); %F2
FC_FPN_SCD = FC_big_FPN_dataF2_ordered(22:60,1:3321);
FC_FPN_aMCI = FC_big_FPN_dataF2_ordered(61:88,1:3321);
FC_FPN_mMCI = FC_big_FPN_dataF2_ordered(89:106,1:3321);
FC_FPN_AD = FC_big_FPN_dataF2_ordered(107:114,1:3321);

FC_FPN_C = FC_big_FPN_data_post_prediff_ordered(1:21,1:3321); %diff (post-pre)
FC_FPN_SCD = FC_big_FPN_data_post_prediff_ordered(22:60,1:3321);
FC_FPN_aMCI = FC_big_FPN_data_post_prediff_ordered(61:88,1:3321);
FC_FPN_mMCI = FC_big_FPN_data_post_prediff_ordered(89:106,1:3321);
FC_FPN_AD = FC_big_FPN_data_post_prediff_ordered(107:114,1:3321);

%non-residualised FC_FPN data (n = 223, F0): 
FC_FPN_C = FC_big_FPN_data_ordered(1:34,1:3321);
FC_FPN_SCD = FC_big_FPN_data_ordered(35:93,1:3321);
FC_FPN_aMCI = FC_big_FPN_data_ordered(94:147,1:3321);
FC_FPN_mMCI = FC_big_FPN_data_ordered(148:199,1:3321);
FC_FPN_AD = FC_big_FPN_data_ordered(200:223,1:3321);


%residualise behav data by age
resid_inhib_behav5_data = residualize(age_data, inhib_behav5_group_ordered);
resid_proc_speed_behav3_data = residualize(age_data, proc_speed3_group_ordered);

%residualise brain data by age and sex: 
resid_FC_FPN_data = residualize(two_covars_ordered, FC_big_FPN_dataF0_ordered); %F0
resid_FC_FPN_C = resid_FC_FPN_data(1:21,1:3321);
resid_FC_FPN_SCD = resid_FC_FPN_data(22:61,1:3321);
resid_FC_FPN_aMCI = resid_FC_FPN_data(62:90,1:3321);
resid_FC_FPN_mMCI = resid_FC_FPN_data(91:109,1:3321);
resid_FC_FPN_AD = resid_FC_FPN_data(110:119,1:3321);

resid_FC_FPN_data = residualize(two_covars_ordered, FC_big_FPN_dataF2_ordered); %F2
resid_FC_FPN_C = resid_FC_FPN_data(1:21,1:3321);
resid_FC_FPN_SCD = resid_FC_FPN_data(22:60,1:3321);
resid_FC_FPN_aMCI = resid_FC_FPN_data(61:88,1:3321);
resid_FC_FPN_mMCI = resid_FC_FPN_data(89:106,1:3321);
resid_FC_FPN_AD = resid_FC_FPN_data(107:114,1:3321);

resid_FC_FPN_data = residualize(two_covars_ordered, FC_big_FPN_data_post_prediff_ordered); %post-pre diff
resid_FC_FPN_C = resid_FC_FPN_data(1:21,1:3321);
resid_FC_FPN_SCD = resid_FC_FPN_data(22:60,1:3321);
resid_FC_FPN_aMCI = resid_FC_FPN_data(61:88,1:3321);
resid_FC_FPN_mMCI = resid_FC_FPN_data(89:106,1:3321);
resid_FC_FPN_AD = resid_FC_FPN_data(107:114,1:3321);

%residualised FC_FPN data (n = 223, F0):
resid_FC_big_FPN_data_ordered = residualize(two_covars_ordered, FC_big_FPN_data_ordered); %F0
resid_FC_FPN_C = resid_FC_big_FPN_data_ordered(1:34,1:3321);
resid_FC_FPN_SCD = resid_FC_big_FPN_data_ordered(35:93,1:3321);
resid_FC_FPN_aMCI = resid_FC_big_FPN_data_ordered(94:147,1:3321);
resid_FC_FPN_mMCI = resid_FC_big_FPN_data_ordered(148:199,1:3321);
resid_FC_FPN_AD = resid_FC_big_FPN_data_ordered(200:223,1:3321);

%small FPN
FC_small_FPN_dataF0_ordered = table2array(readtable('FC_small_FPN_data_ordered_by_groupF0(n=119).xlsx')); %F0 
FC_small_FPN_dataF2_ordered = table2array(readtable('FC_small_FPN_data_ordered_by_groupF2(n=114).xlsx')); %F2
FC_small_FPN_data_post_prediff_ordered = table2array(readtable('FC_small_FPN_data_ordered_by_group_post-prediff(n=114).xlsx')); %diff (post-pre)


%non-residualised FC_small_FPN data: 
FC_FPN_C = FC_small_FPN_dataF0_ordered(1:21,1:10); %F0
FC_FPN_SCD = FC_small_FPN_dataF0_ordered(22:61,1:10);
FC_FPN_aMCI = FC_small_FPN_dataF0_ordered(62:90,1:10);
FC_FPN_mMCI = FC_small_FPN_dataF0_ordered(91:109,1:10);
FC_FPN_AD = FC_small_FPN_dataF0_ordered(110:119,1:10);

FC_FPN_C = FC_small_FPN_dataF2_ordered(1:21,1:10); %F2
FC_FPN_SCD = FC_small_FPN_dataF2_ordered(22:60,1:10);
FC_FPN_aMCI = FC_small_FPN_dataF2_ordered(61:88,1:10);
FC_FPN_mMCI = FC_small_FPN_dataF2_ordered(89:106,1:10);
FC_FPN_AD = FC_small_FPN_dataF2_ordered(107:114,1:10);

FC_FPN_C = FC_small_FPN_data_post_prediff_ordered(1:21,1:10); %post-pre diff
FC_FPN_SCD = FC_small_FPN_data_post_prediff_ordered(22:60,1:10);
FC_FPN_aMCI = FC_small_FPN_data_post_prediff_ordered(61:88,1:10);
FC_FPN_mMCI = FC_small_FPN_data_post_prediff_ordered(89:106,1:10);
FC_FPN_AD = FC_small_FPN_data_post_prediff_ordered(107:114,1:10);


%residualise brain data by age and sex: 
resid_FC_FPN_data = residualize(two_covars_ordered, FC_small_FPN_dataF0_ordered); %F0
resid_FC_FPN_C = resid_FC_FPN_data(1:21,1:10);
resid_FC_FPN_SCD = resid_FC_FPN_data(22:61,1:10);
resid_FC_FPN_aMCI = resid_FC_FPN_data(62:90,1:10);
resid_FC_FPN_mMCI = resid_FC_FPN_data(91:109,1:10);
resid_FC_FPN_AD = resid_FC_FPN_data(110:119,1:10);

resid_FC_FPN_data = residualize(two_covars_ordered, FC_small_FPN_dataF2_ordered); %F2
resid_FC_FPN_C = resid_FC_FPN_data(1:21,1:10);
resid_FC_FPN_SCD = resid_FC_FPN_data(22:60,1:10);
resid_FC_FPN_aMCI = resid_FC_FPN_data(61:88,1:10);
resid_FC_FPN_mMCI = resid_FC_FPN_data(89:106,1:10);
resid_FC_FPN_AD = resid_FC_FPN_data(107:114,1:10);

resid_FC_FPN_data = residualize(two_covars_ordered, FC_small_FPN_data_post_prediff_ordered); %post-pre diff
resid_FC_FPN_C = resid_FC_FPN_data(1:21,1:10);
resid_FC_FPN_SCD = resid_FC_FPN_data(22:60,1:10);
resid_FC_FPN_aMCI = resid_FC_FPN_data(61:88,1:10);
resid_FC_FPN_mMCI = resid_FC_FPN_data(89:106,1:10);
resid_FC_FPN_AD = resid_FC_FPN_data(107:114,1:10);



% %Run behav PLS by separating groups and compare groups with correlation between neural x behav scores 
% datamat_list{1} = cat(1, resid_FD_C);
% datamat_list{2} = cat(1, resid_FD_SCD);
% datamat_list{3} = cat(1, resid_FD_aMCI);
% datamat_list{4} = cat(1, resid_FD_mMCI);
% datamat_list{5} = cat(1, resid_FD_AD);
% num_subj_lst = [length(resid_FD_C)  length(resid_FD_SCD) length(resid_FD_aMCI) length(resid_FD_mMCI) length(resid_FD_AD)];
% num_cond = 1;
% option.stacked_behavdata = resid_inhib_behav5_data; %data must be in group order (same order as the neural behaviour)
% option.method = 5; %non-rotated
% option.num_perm = 5000;
% option.num_boot = 5000;
% option.clim = 95;
% option.stacked_designdata = ones(width(resid_inhib_behav5_data)*length(Group_names),1)*-1; %should be # of tasks * # of groups
% resultRotated_group_nonrotated_behav_PLS = pls_analysis(datamat_list, num_subj_lst, num_cond,option);
% 
% 
% %make figure
% figure
% bar(resultRotated_group_nonrotated_behav_PLS.boot_result.orig_corr); 
% xticklabels(Group_names); 
% Xt=round(width(option.stacked_behavdata)/2):width(option.stacked_behavdata):length(option.stacked_designdata);
% set(gca,'XTick',Xt);
% hold on
% for i = 1:length(option.stacked_designdata)
%     plot([i i], [resultRotated_group_nonrotated_behav_PLS.boot_result.ulcorr(i,1) resultRotated_group_nonrotated_behav_PLS.boot_result.llcorr(i,1)],'-k', 'LineWidth', 2);
% end
% xlabel('Group');
% ylabel('Bootstrap Ratio Correlation Value');
% title(['LV p = ' num2str(resultRotated_group_nonrotated_behav_PLS.perm_result.sprob(1))]);
% box off 
% set(gcf,'color','w'); 


%alternative figure (by group)
% figure
% group_figure = transpose(reshape(resultRotated_group_nonrotated_behav_PLS.boot_result.orig_corr, 5, 5));
% bar(group_figure);
% xticklabels(Group_names); 


% %Run Group contrast PLS - mean-centred
% datamat_list{1} = cat(1, FD_C);
% datamat_list{2} = cat(1, FD_SCD);
% datamat_list{3} = cat(1, FD_aMCI);
% datamat_list{4} = cat(1, FD_mMCI);
% datamat_list{5} = cat(1, FD_AD);
% num_subj_lst = [length(FD_C)  length(FD_SCD) length(FD_aMCI) length(FD_mMCI) length(FD_AD)];
% num_cond = 1;
% %option.stacked_behavdata = all_12_behav_raw_scores; %is behaviour scores being used here??? no - this is b/c we selected option 1 for method(mean-centred)
% option.method = 1; %mean-centred
% option.num_perm = 5000;
% option.num_boot = 5000;
% option.clim = 95;
% option.stacked_designdata = [2;1;0;-1;-2];
% option.meancentering_type = 0;
% resultRotated_group_differences_meancentred = pls_analysis(datamat_list, num_subj_lst, num_cond,option);
% %save('FD_behav_PLS_raw_scores.mat', 'resultRotated');

%make figure for group contrasts
% figure
% bar(resultRotated.boot_result.orig_usc(:,1)); %the bootstrap brainscores - we are looking at the first column
% hold on
% for i = 1:5
% plot([i i], [resultRotated.boot_result.ulusc(i,1) resultRotated.boot_result.llusc(i,1)],'-k', 'LineWidth', 2);
% end
% xticklabels(Group_names); 
% ylabel('Bootstrap Ratio Correlation Value');
% title(['LV p = ' num2str(resultRotated.perm_result.sprob(1))]);




% %example from the paper - groups are not ordered correctly (e.g., AD, PD, NC)
% U = [0.41,-0.42,0.37,-0.02,-0.14,-0.71;-0.41,0.44,-0.40,0.00,-0.14,-0.68;-0.43,0.25,0.74,0.46,0.01,-0.02;-0.07,0.31,0.40,-0.86,0.08,-0.01;-0.44,-0.47,-0.06,-0.09,0.74,-0.15;0.53,0.51,0.00,0.20,0.64,-0.13];
% figure
% bar(U(:,6)); %the bootstrap salience score - we are looking at the first + second column column



%Run PLS with Splithalf Resampling:
%note that this will take much longer & more computational power, so it
%may be good to run on a vm instead


% %set up PLS options
% datamat_list{1} = cat(1, resid_FD_data);
% num_subj_lst = height(resid_FD_data);
% num_cond = 1;
% option.stacked_behavdata = resid_behav_data; 
% option.method = 3; %regular behaviour PLS
% option.num_perm = 1000;
% option.num_boot = 5000;
% option.clim = 95;
% option.stacked_designdata = contrast_5; 
% option.num_split = 1000;
% resultRotated3 = pls_analysis(datamat_list, num_subj_lst, num_cond,option); 
% 
% %save('FD_behav_PLS_raw_scores.mat', 'resultRotated'); 
% 
% %4. Make PLS figures
% figure 
% bar((resultRotated3.boot_result.orig_corr)*-1); 
% hold on 
% for i = 1:length(resultRotated3.boot_result.orig_corr) %what is this - bootstrap correlation value?
%     plot([i i], [(resultRotated3.boot_result.ulcorr(i)*-1) (resultRotated3.boot_result.llcorr(i)*-1)],'-k', 'LineWidth',2); %95 percent CI - upper and lower bound
% end
% xticklabels(inhib_5names2); 
% xtickangle(45); 
% xlabel('Neuropsychological Assessment');
% ylabel('Bootstrap Ratio Correlation Value');
% box off 
% set(gcf,'color','w'); 
% title(['LV p = ' num2str(resultRotated3.perm_result.sprob)]);






%Run behav PLS by separating groups and compare groups with correlation between neural x behav scores 
datamat_list{1} = cat(1, resid_FDC_C);
datamat_list{2} = cat(1, resid_FDC_SCD);
datamat_list{3} = cat(1, resid_FDC_aMCI);
datamat_list{4} = cat(1, resid_FDC_mMCI);
%datamat_list{5} = cat(1, resid_FD_AD);
%num_subj_lst = [height(resid_FD_C)  height(resid_FD_SCD) height(resid_FD_aMCI) height(resid_FD_mMCI) height(resid_FD_AD)];
num_subj_lst = [height(resid_FD_C)  height(resid_FD_SCD) height(resid_FD_aMCI) height(resid_FD_mMCI)];
num_cond = 1;
option.stacked_behavdata = inhib_behav5_group_ordered; %data must be in group order (same order as the neural behaviour)
option.method = 3; %regular behaviour PLS
option.num_perm = 5000;
option.num_boot = 5000;
option.clim = 95;
option.num_split = 100;
%option.stacked_designdata = ones(width(inhib_behav5_group_ordered)*length(Group_names),1)*-1; %should be # of tasks * # of groups
option.stacked_designdata = ones(width(inhib_behav5_group_ordered)*length(Group_names),1); %should be # of tasks * # of groups
%option.stacked_designdata = ones(width(ProcSpeed_behav3_group_ordered)*length(Group_names),1); %should be # of tasks * # of groups
resultRotated_group_rotated_behav_PLS = pls_analysis(datamat_list, num_subj_lst, num_cond,option);


%make figure
figure
bar(resultRotated_group_rotated_behav_PLS.boot_result.orig_corr); 
xticklabels(Group_names); 
Xt=round(width(option.stacked_behavdata)/2):width(option.stacked_behavdata):length(option.stacked_designdata);
set(gca,'XTick',Xt);
hold on
for i = 1:length(option.stacked_designdata)
    plot([i i], [resultRotated_group_rotated_behav_PLS.boot_result.ulcorr(i,1) resultRotated_group_rotated_behav_PLS.boot_result.llcorr(i,1)],'-k', 'LineWidth', 2);
end
xlabel('Group');
ylabel('Bootstrap Ratio Correlation Value');
title(['LV p = ' num2str(resultRotated_group_rotated_behav_PLS.perm_result.sprob(1))]);
box off 
set(gcf,'color','w'); 




