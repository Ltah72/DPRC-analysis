function CreateStudyMatrices(excelFile,fileLocation, startdir, groupname, period)
%Create your matrices (design_matrix.txt and contrast_matrix.txt) to run
%your statistics on them using the general linear model (GLM). These
%matrices have been customised to my own study design, in which I am
%comparing 5 groups from the DPRC cohort (HC, SCD, aMCI, mMCI, and AD),
%plus the specifcied covariates (age, sex, and ACE score). If you have your
%own groups for your own study design, then you can create a new function
%to create your matrices, or you can manually create your own matrices.

%This function will read off of the excel file which contains all of the
%DPRC participant information. Make sure that the participants that you
%will be analysing with the given pipeline matches up with the participants
%in the excel file.

%      Inputs (5): excelFile - Name of the excel file containing
%                               particpant info.
%                  fileLocation - directory of where the excelFile is
%                                 contained.
%                  startdir - directory where your source data is
%                              stored
%                  groupname - name of the group/study that you are
%                              analysing
%                  period - time period of the participant MRI scans

%       Outputs (none): You writing the matrices to a text file.


%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 05/09/20

%navigate to the most basic file directory and then to the file location
cd /;
cd ([fileLocation]);

%read in excel file into matlab format
[num, txt, raw] = xlsread(excelFile, 'A:I');

%define variables
Auckland = 'ADPRC';
Christchurch = 'CDPRC';
Dunedin = 'DDPRC';

%calculate average scores to input as covariates
sumAge = 0;
sumSex = 0;
sumClinicSite = 0;
sumACE = 0;
sumGroup = 0;
NumParticipants = 0;
for i = 2:228
    sumAge = raw{i,2} + sumAge;
    sumSex = raw{i,7} + sumSex;
    sumGroup = raw{i,5} + sumGroup;
    if contains(raw{i,9},Auckland)
        sumClinicSite = 0 + sumClinicSite;
    elseif contains(raw{i,9},Christchurch)
        sumClinicSite = 1 + sumClinicSite;
    elseif contains(raw{i,9},Dunedin)
        sumClinicSite = 2 + sumClinicSite;
    end
    %if isnan(raw{i,5}) ~= 1
     %   sumACE = raw{i,6} + sumACE;
    %end
    NumParticipants = NumParticipants + 1;
end
MeanAge = sumAge / NumParticipants;
MeanSex = sumSex / NumParticipants;
MeanClinicSite = sumClinicSite / NumParticipants;
%MeanACE = sumACE / (NumParticipants-7);
MeanGroup = sumGroup / NumParticipants;

%create your matrices for the study

%Create design matrix with the covariates. Added covariates currently are:
%age, sex, clinical site, and overall ACE-III score.
fid = fopen('design_matrix_group_diff_cov-age_sex.txt', 'w');
if (fid == -1)
    disp('Error in creating the text file')
else
    for i = 2:228
        if raw{i,5} ~= 0 && raw{i,5} ~= -1 && raw{i,5} ~= 6
            classification = raw{i,5};
            if (classification == 1)
                matrix_line = '1 0 0 0 0';
            elseif (classification == 2)
                matrix_line = '0 1 0 0 0';
            elseif (classification == 3)
                matrix_line = '0 0 1 0 0';
            elseif (classification == 4)
                matrix_line = '0 0 0 1 0';
            elseif (classification == 5)
                matrix_line = '0 0 0 0 1';
            end
            Norm_Age = raw{i,2} - MeanAge;
            Norm_Sex = raw{i,7} - MeanSex;
            Norm_Group = raw{i,5} - MeanGroup;
            if contains(raw{i,9},Auckland)
                Norm_ClinicSite = 0 - MeanClinicSite;
            elseif contains(raw{i,9},Christchurch)
                Norm_ClinicSite = 1 - MeanClinicSite;
            elseif contains(raw{i,9},Dunedin)
                Norm_ClinicSite = 2 - MeanClinicSite;
            end
            %if contains(raw{i,9},'F0')
             %   timepoint = '1';
            %elseif contains(raw{i,9},'F2')
             %   timepoint = '-1';
            %end 
            %Norm_ACE = raw{i,6} - MeanACE;
            %fprintf(fid, '%s', matrix_line);
            fprintf(fid, '%s %.2f %.2f', matrix_line, Norm_Age, Norm_Sex);
            %fprintf(fid, '%.2f %.2f %.2f', Norm_Group, Norm_Age, Norm_Sex);
            %fprintf(fid, '%s %.2f %.2f %.2f', matrix_line, Norm_Age, Norm_Sex, Norm_ClinicSite);
            %fprintf(fid, '%s %s', timepoint, matrix_line);
            %fprintf(fid, '%s %s %.2f %.2f', timepoint, matrix_line, Norm_Age, Norm_Sex);
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end

%create your associated contrast matrix file
fid2 = fopen('contrast_matrix.txt', 'w');
if (fid2 == -1)
    disp('Error in creating the text file')
else
    Comparison1 = '1 -1 0 0 0 0 0 0';
    Comparison2 = '1 0 -1 0 0 0 0 0';
    Comparison3 = '1 0 0 -1 0 0 0 0';
    Comparison4 = '1 0 0 0 -1 0 0 0';
    Comparison5 = '0 1 -1 0 0 0 0 0';
    Comparison6 = '0 1 0 -1 0 0 0 0';
    Comparison7 = '0 1 0 0 -1 0 0 0';
    Comparison8 = '0 0 1 -1 0 0 0 0';
    Comparison9 = '0 0 1 0 -1 0 0 0';
    Comparison10 = '0 0 0 1 -1 0 0 0';
    
    %Cov1Pos = '0 0 0 0 0 0 0 0 1';
    %Cov1Neg = '0 0 0 0 0 0 0 0 -1';
    
    fprintf(fid2, '%s', Comparison1);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison2);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison3);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison4);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison5);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison6);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison7);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison8);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison9);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison10);
    %fprintf(fid2, '\n');
    %fprintf(fid2, '%s', Cov1Pos);
    %fprintf(fid2, '\n');
    %fprintf(fid2, '%s', Cov1Neg);
    
    fclose(fid2);
end

%move both your matrices files into the stats_matrices folder.
movefile('design_matrix.txt', [startdir, '/derivatives/' period, '/diff_data/' groupname, '/stats_matrices']);
movefile('contrast_matrix.txt', [startdir, '/derivatives/' period, '/diff_data/' groupname, '/stats_matrices']);









%make study matrices for testing differences between the clinical sites

%navigate to the most basic file directory and then to the file location
excelFile = 'covariates-participants-lined-up_update.xlsx';

%read in excel file into matlab format
[num, txt, raw] = xlsread(excelFile, 'A:H');

%define variables
Auckland = 'ADPRC';
Christchurch = 'CDPRC';
Dunedin = 'DDPRC';

%calculate average scores to input as covariates
sumAge = 0;
sumSex = 0;
sumClinicSite = 0;
sumACE = 0;
sumGroup = 0;
NumParticipants = 0;
for i = 2:length(txt)
    sumAge = raw{i,7} + sumAge;
    sumSex = raw{i,5} + sumSex;
    sumGroup = raw{i,3} + sumGroup;
    if contains(raw{i,1},Auckland)
        sumClinicSite = 0 + sumClinicSite;
    elseif contains(raw{i,1},Christchurch)
        sumClinicSite = 1 + sumClinicSite;
    elseif contains(raw{i,1},Dunedin)
        sumClinicSite = 2 + sumClinicSite;
    end
    if isnan(raw{i,6}) ~= 1
        sumACE = raw{i,6} + sumACE;
    end
    NumParticipants = NumParticipants + 1;
end
MeanAge = sumAge / NumParticipants;
MeanSex = sumSex / NumParticipants;
MeanClinicSite = sumClinicSite / NumParticipants;
MeanACE = sumACE / (NumParticipants-7);
MeanGroup = sumGroup / NumParticipants;


%Create design matrix for testing the clinical site.
fid = fopen('design_matrix_clinsite_covar-groups.txt', 'w');
if (fid == -1)
    disp('Error in creating the text file')
else
    for i = 2:length(txt)
        if contains(raw{i,1},Auckland)
            matrix_line = '1 0 0';
        elseif contains(raw{i,1},Christchurch)
            matrix_line = '0 1 0';
        elseif contains(raw{i,1},Dunedin)
            matrix_line = '0 0 1';
        end
        Norm_Age = raw{i,7} - MeanAge;
        Norm_Sex = raw{i,5} - MeanSex;
        Norm_Group = raw{i,3} - MeanGroup;
        if contains(raw{i,1},Auckland)
            Norm_ClinicSite = 0 - MeanClinicSite;
        elseif contains(raw{i,1},Christchurch)
            Norm_ClinicSite = 1 - MeanClinicSite;
        elseif contains(raw{i,1},Dunedin)
            Norm_ClinicSite = 2 - MeanClinicSite;
        end
        fprintf(fid, '%s %.2f %.2f %.2f', matrix_line, Norm_Group);
        fprintf(fid, '\n');
    end
end
fclose(fid);












%Create design matrix with 3 groups (non-MCI, MCI, AD). Added covariates currently are:
%age, sex, clinical site, and overall ACE-III score.
fid = fopen('design_matrix_3groups_group_diff.txt', 'w');
if (fid == -1)
    disp('Error in creating the text file')
else
    for i = 2:length(txt)
        if raw{i,5} ~= 0 && raw{i,5} ~= -1 && raw{i,5} ~= 6
            classification = raw{i,5};
            if (classification == 1) || (classification == 2)
                matrix_line = '1 0 0';
            elseif (classification == 3) || (classification == 4)
                matrix_line = '0 1 0';
            elseif (classification == 5)
                matrix_line = '0 0 1';
            end
            %Norm_Age = raw{i,2} - MeanAge;
            %Norm_Sex = raw{i,6} - MeanSex;
            %Norm_Group = raw{i,5} - MeanGroup;
            %if contains(raw{i,9},Auckland)
             %   Norm_ClinicSite = 0 - MeanClinicSite;
            %elseif contains(raw{i,9},Christchurch)
            %    Norm_ClinicSite = 1 - MeanClinicSite;
            %elseif contains(raw{i,9},Dunedin)
             %   Norm_ClinicSite = 2 - MeanClinicSite;
            %end
            %Norm_ACE = raw{i,10} - MeanACE;
            fprintf(fid, '%s', matrix_line);
            %fprintf(fid, '%s %.2f', matrix_line, Norm_Age);
            %fprintf(fid, '%.2f %.2f %.2f', Norm_Group, Norm_Age, Norm_Sex);
            %fprintf(fid, '%s %.2f %.2f %.2f', matrix_line, Norm_Age, Norm_Sex, Norm_ClinicSite);
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end

%create your associated contrast matrix file (3 groups)
fid2 = fopen('contrast_matrix_3groups_group_diff.txt', 'w');
if (fid2 == -1)
    disp('Error in creating the text file')
else
    Comparison1 = '1 -1 0';
    Comparison2 = '1 0 -1';
    Comparison3 = '0 1 -1';
    
    fprintf(fid2, '%s', Comparison1);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison2);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Comparison3);
    
    
    fclose(fid2);
end






%matrices for correlation analysis


%navigate to the most basic file directory and then to the file location
fileLocation = 'H:ltah262/PhD/ExecutiveFunction/NeuroPsychAssessment/data/June/';

cd /;
cd ([fileLocation]);

%read in excel file into matlab format
excelNeuroPsychFile = 'DPRC_neuropsych_data_lined_up_valid_participants.xlsx';
[num, txt, raw] = xlsread(excelNeuroPsychFile, 'A:BG');

%define variables
Auckland = 'ADPRC';
Christchurch = 'CDPRC';
Dunedin = 'DDPRC';

%convert neuropsych variables to double type
for i = 2:length(txt)
    if isnan(raw{i,23}) ~= 1
        raw{i,23} = str2double(raw{i,23});
    end
    if isnan(raw{i,24}) ~= 1
        raw{i,24} = str2double(raw{i,24});
    end
    if isnan(raw{i,25}) ~= 1
        raw{i,25} = str2double(raw{i,25});
    end
    if isnan(raw{i,26}) ~= 1
        raw{i,26} = str2double(raw{i,26});
    end
    if isnan(raw{i,27}) ~= 1
        raw{i,27} = str2double(raw{i,27});
    end
    if isnan(raw{i,28}) ~= 1
        raw{i,28} = str2double(raw{i,28});
    end
    if isnan(raw{i,30}) ~= 1
        raw{i,30} = str2double(raw{i,30});
    end
    if isnan(raw{i,31}) ~= 1
        raw{i,31} = str2double(raw{i,31});
    end
    if isnan(raw{i,33}) ~= 1
        raw{i,33} = str2double(raw{i,33});
    end
    if isnan(raw{i,34}) ~= 1
        raw{i,34} = str2double(raw{i,34});
    end
    if isnan(raw{i,39}) ~= 1
        raw{i,39} = str2double(raw{i,39});
    end
    if isnan(raw{i,40}) ~= 1
        raw{i,40} = str2double(raw{i,40});
    end
    if isnan(raw{i,42}) ~= 1
        raw{i,42} = str2double(raw{i,42});
    end
    if isnan(raw{i,43}) ~= 1
        raw{i,43} = str2double(raw{i,43});
    end
    if isnan(raw{i,48}) ~= 1
        raw{i,48} = str2double(raw{i,48});
    end
    if isnan(raw{i,49}) ~= 1
        raw{i,49} = str2double(raw{i,49});
    end
    if isnan(raw{i,51}) ~= 1
        raw{i,51} = str2double(raw{i,51});
    end
    if isnan(raw{i,52}) ~= 1
        raw{i,52} = str2double(raw{i,52});
    end
    if isnan(raw{i,53}) ~= 1
        raw{i,53} = str2double(raw{i,53});
    end
    if isnan(raw{i,54}) ~= 1
        raw{i,54} = str2double(raw{i,54});
    end
    if isnan(raw{i,55}) ~= 1
        raw{i,55} = str2double(raw{i,55});
    end
    if isnan(raw{i,56}) ~= 1
        raw{i,56} = str2double(raw{i,56});
    end
    if isnan(raw{i,57}) ~= 1
        raw{i,57} = str2double(raw{i,57});
    end
    if isnan(raw{i,58}) ~= 1
        raw{i,58} = str2double(raw{i,58});
    end
end

%calculate average scores to input as variables of interest for correlation
%analysis
sumTrailsARaw = 0;
sumTrailsAZ = 0;
sumTrailsBRaw = 0;
sumTrailsBZ = 0;
sumColorNamingRaw = 0;
sumColorNamingZ = 0;
sumWordReadingRaw = 0;
sumWordReadingZ = 0;
sumInhibitionRaw = 0;
sumInhibitionZ = 0;
sumLetFluencyRaw = 0;
sumLetFluencyZ = 0;
sumCatFluencyRaw = 0;
sumCatFluencyZ = 0;
sumSwitchingRaw = 0;
sumSwitchingZ = 0;
sumHayBTime1Raw = 0;
sumHayBTime1Z = 0;
sumHayBTime2Raw = 0;
sumHayBTime2Z = 0;
sumHayBCatARaw = 0;
sumHayBCatAZ = 0;
sumHayBCatBRaw = 0;
sumHayBCatBZ = 0;

%demographic information (which may be used as covariates in the analysis).
sumAge = 0;
sumSex = 0;
sumClinicSite = 0;
sumACE = 0;
NumParticipants = 0;

NumParticipants_TrailsARaw = 0;
NumParticipants_TrailsAZ = 0;
NumParticipants_TrailsBRaw = 0;
NumParticipants_TrailsBZ = 0;
NumParticipants_ColorNamingRaw = 0;
NumParticipants_ColorNamingZ = 0;
NumParticipants_WordReadingRaw = 0;
NumParticipants_WordReadingZ = 0;
NumParticipants_InhibitionRaw = 0;
NumParticipants_InhibitionZ = 0;
NumParticipants_LetFluencyRaw = 0;
NumParticipants_LetFluencyZ = 0;
NumParticipants_CatFluencyRaw = 0;
NumParticipants_CatFluencyZ = 0;
NumParticipants_SwitchingRaw = 0;
NumParticipants_SwitchingZ = 0;
NumParticipants_HayBTime1Raw = 0;
NumParticipants_HayBTime1Z = 0;
NumParticipants_HayBTime2Raw = 0;
NumParticipants_HayBTime2Z = 0;
NumParticipants_HayBCatARaw = 0;
NumParticipants_HayBCatAZ = 0;
NumParticipants_HayBCatBRaw = 0;
NumParticipants_HayBCatBZ = 0;



for i = 2:length(txt)
    sumAge = raw{i,2} + sumAge;
    sumSex = raw{i,7} + sumSex;
    
    if isnan(raw{i,59}) ~= 1
        sumACE = raw{i,59} + sumACE;
    end

    if contains(raw{i,1},Auckland)
        sumClinicSite = 0 + sumClinicSite;
    elseif contains(raw{i,1},Christchurch)
        sumClinicSite = 1 + sumClinicSite;
    elseif contains(raw{i,1},Dunedin)
        sumClinicSite = 2 + sumClinicSite;
    end
    
    %to not include NaNs in the summation. 
    if isnan(raw{i,23}) ~= 1
        sumTrailsARaw = raw{i,23} + sumTrailsARaw;
    end
    if isnan(raw{i,24}) ~= 1
        sumTrailsAZ = raw{i,24} + sumTrailsAZ;
    end
    if isnan(raw{i,25}) ~= 1
        sumTrailsBRaw = raw{i,25} + sumTrailsBRaw;
    end
    if isnan(raw{i,26}) ~= 1
        sumTrailsBZ = raw{i,26} + sumTrailsBZ;
    end
    if isnan(raw{i,27}) ~= 1
        sumColorNamingRaw = raw{i,27} + sumColorNamingRaw;
    end
    if isnan(raw{i,28}) ~= 1
        sumColorNamingZ = raw{i,28} + sumColorNamingZ;
    end
    if isnan(raw{i,30}) ~= 1
        sumWordReadingRaw = raw{i,30} + sumWordReadingRaw;
    end
    if isnan(raw{i,31}) ~= 1
        sumWordReadingZ = raw{i,31} + sumWordReadingZ;
    end
    if isnan(raw{i,33}) ~= 1
        sumInhibitionRaw = raw{i,33} + sumInhibitionRaw;
    end
    if isnan(raw{i,34}) ~= 1
        sumInhibitionZ = raw{i,34} + sumInhibitionZ;
    end
    if isnan(raw{i,39}) ~= 1
        sumLetFluencyRaw = raw{i,39} + sumLetFluencyRaw;
    end
    if isnan(raw{i,40}) ~= 1
        sumLetFluencyZ = raw{i,40} + sumLetFluencyZ;
    end
    if isnan(raw{i,42}) ~= 1
        sumCatFluencyRaw = raw{i,42} + sumCatFluencyRaw;
    end
    if isnan(raw{i,43}) ~= 1
        sumCatFluencyZ = raw{i,43} + sumCatFluencyZ;
    end
    if isnan(raw{i,48}) ~= 1
        sumSwitchingRaw = raw{i,48} + sumSwitchingRaw;
    end
    if isnan(raw{i,49}) ~= 1
        sumSwitchingZ = raw{i,49} + sumSwitchingZ;
    end
    if isnan(raw{i,51}) ~= 1
        sumHayBTime1Raw = raw{i,51} + sumHayBTime1Raw;
    end
    if isnan(raw{i,52}) ~= 1
        sumHayBTime1Z = raw{i,52} + sumHayBTime1Z;
    end
    if isnan(raw{i,53}) ~= 1
        sumHayBTime2Raw = raw{i,53} + sumHayBTime2Raw;
    end
    if isnan(raw{i,54}) ~= 1
        sumHayBTime2Z = raw{i,54} + sumHayBTime2Z;
    end
    if isnan(raw{i,55}) ~= 1
        sumHayBCatARaw = raw{i,55} + sumHayBCatARaw;
    end
    if isnan(raw{i,56}) ~= 1
        sumHayBCatAZ = raw{i,56} + sumHayBCatAZ;
    end
    if isnan(raw{i,57}) ~= 1
        sumHayBCatBRaw = raw{i,57} + sumHayBCatBRaw;
    end
    if isnan(raw{i,58}) ~= 1
        sumHayBCatBZ = raw{i,58} + sumHayBCatBZ;
    end
    
    NumParticipants = NumParticipants + 1;
    
    %if any data is missing, then do not count and add this to the number
    %of current participants
    if isnan(raw{i,23}) ~= 1
        NumParticipants_TrailsARaw = NumParticipants_TrailsARaw + 1;
    end
    if isnan(raw{i,24}) ~= 1
        NumParticipants_TrailsAZ = NumParticipants_TrailsAZ + 1;
    end
    if isnan(raw{i,25}) ~= 1
        NumParticipants_TrailsBRaw = NumParticipants_TrailsBRaw + 1;
    end
    if isnan(raw{i,26}) ~= 1
        NumParticipants_TrailsBZ = NumParticipants_TrailsBZ + 1;
    end
    if isnan(raw{i,27}) ~= 1
        NumParticipants_ColorNamingRaw = NumParticipants_ColorNamingRaw + 1;
    end
    if isnan(raw{i,28}) ~= 1
        NumParticipants_ColorNamingZ = NumParticipants_ColorNamingZ + 1;
    end
    if isnan(raw{i,30}) ~= 1
        NumParticipants_WordReadingRaw = NumParticipants_WordReadingRaw + 1;
    end
    if isnan(raw{i,31}) ~= 1
        NumParticipants_WordReadingZ = NumParticipants_WordReadingZ + 1;
    end
    if isnan(raw{i,33}) ~= 1
        NumParticipants_InhibitionRaw = NumParticipants_InhibitionRaw + 1;
    end
    if isnan(raw{i,34}) ~= 1
        NumParticipants_InhibitionZ = NumParticipants_InhibitionZ + 1;
    end
    if isnan(raw{i,39}) ~= 1
        NumParticipants_LetFluencyRaw = NumParticipants_LetFluencyRaw + 1;
    end
    if isnan(raw{i,40}) ~= 1
        NumParticipants_LetFluencyZ = NumParticipants_LetFluencyZ + 1;
    end
    if isnan(raw{i,42}) ~= 1
        NumParticipants_CatFluencyRaw = NumParticipants_CatFluencyRaw + 1;
    end
    if isnan(raw{i,43}) ~= 1
        NumParticipants_CatFluencyZ = NumParticipants_CatFluencyZ + 1;
    end
    if isnan(raw{i,48}) ~= 1
        NumParticipants_SwitchingRaw = NumParticipants_SwitchingRaw + 1;
    end
    if isnan(raw{i,49}) ~= 1
        NumParticipants_SwitchingZ = NumParticipants_SwitchingZ + 1;
    end
    if isnan(raw{i,51}) ~= 1
        NumParticipants_HayBTime1Raw = NumParticipants_HayBTime1Raw + 1;
    end
    if isnan(raw{i,52}) ~= 1
        NumParticipants_HayBTime1Z = NumParticipants_HayBTime1Z + 1;
    end
    if isnan(raw{i,53}) ~= 1
        NumParticipants_HayBTime2Raw = NumParticipants_HayBTime2Raw + 1;
    end
    if isnan(raw{i,54}) ~= 1
        NumParticipants_HayBTime2Z = NumParticipants_HayBTime2Z + 1;
    end
    if isnan(raw{i,55}) ~= 1
        NumParticipants_HayBCatARaw = NumParticipants_HayBCatARaw + 1;
    end
    if isnan(raw{i,56}) ~= 1
        NumParticipants_HayBCatAZ = NumParticipants_HayBCatAZ + 1;
    end
    if isnan(raw{i,57}) ~= 1
        NumParticipants_HayBCatBRaw = NumParticipants_HayBCatBRaw + 1;
    end
    if isnan(raw{i,58}) ~= 1
        NumParticipants_HayBCatBZ = NumParticipants_HayBCatBZ + 1;
    end 
        
end

MeanAge = sumAge / NumParticipants;
MeanSex = sumSex / NumParticipants;
MeanClinicSite = sumClinicSite / NumParticipants;
MeanACE = (sumACE) / (NumParticipants - 3);
%MeanGroup = sumGroup / NumParticipants;
%MeanACE = 88.6039

MeanTrailsARaw = sumTrailsARaw / NumParticipants_TrailsARaw;
MeanTrailsAZ = sumTrailsAZ/ NumParticipants_TrailsAZ;
MeanTrailsBRaw = sumTrailsBRaw / NumParticipants_TrailsBRaw;
MeanTrailsBZ = sumTrailsBZ / NumParticipants_TrailsBZ;
MeanColorNamingRaw = sumColorNamingRaw / NumParticipants_ColorNamingRaw;
MeanColorNamingZ = sumColorNamingZ / NumParticipants_ColorNamingZ;
MeanWordReadingRaw = sumWordReadingRaw / NumParticipants_WordReadingRaw;
MeanWordReadingZ = sumWordReadingZ/ NumParticipants_WordReadingZ;
MeanInhibitionRaw = sumInhibitionRaw / NumParticipants_InhibitionRaw ;
MeanInhibitionZ = sumInhibitionZ / NumParticipants_InhibitionZ;
MeanLetFluencyRaw = sumLetFluencyRaw / NumParticipants_LetFluencyRaw;
MeanLetFluencyZ = sumLetFluencyZ / NumParticipants_LetFluencyZ;
MeanCatFluencyRaw = sumCatFluencyRaw / NumParticipants_CatFluencyRaw;
MeanCatFluencyZ = sumCatFluencyZ / NumParticipants_CatFluencyZ;
MeanSwitchingRaw = sumSwitchingRaw / NumParticipants_SwitchingRaw;
MeanSwitchingZ = sumSwitchingZ / NumParticipants_SwitchingZ;
MeanHayBTime1Raw = sumHayBTime1Raw / NumParticipants_HayBTime1Raw;
MeanHayBTime1Z = sumHayBTime1Z / NumParticipants_HayBTime1Z;
MeanHayBTime2Raw = sumHayBTime2Raw / NumParticipants_HayBTime2Raw;
MeanHayBTime2Z = sumHayBTime2Z / NumParticipants_HayBTime2Z;
MeanHayBCatARaw = sumHayBCatARaw / NumParticipants_HayBCatARaw;
MeanHayBCatAZ = sumHayBCatAZ / NumParticipants_HayBCatAZ;
MeanHayBCatBRaw = sumHayBCatBRaw / NumParticipants_HayBCatBRaw;
MeanHayBCatBZ = sumHayBCatBZ / NumParticipants_HayBCatBZ;


%create your matrices for the study

%Create design matrix with the variables of interest for correlation.
%Utilise the DPRC neuropscyhological assessment excel sheet. Also,
%added covariates currently are:
%age, sex, clinical site, and overall ACE-III score.

%-----------------for overall correlation matrices------------------------%

fid = fopen('design_matrix_overall_correlation_HayBCatBZ.txt', 'w');
if (fid == -1)
    disp('Error in creating the text file')
else
    for i = 2:length(txt)
        matrix_line = '1';
        Norm_Age = raw{i,5} - MeanAge;
        Norm_Sex = raw{i,7} - MeanSex;
        Norm_ACE = raw{i,59} - MeanACE;
        %Norm_Group = raw{i,9} - MeanGroup;
        if contains(raw{i,1},Auckland)
            Norm_ClinicSite = 0 - MeanClinicSite;
        elseif contains(raw{i,1},Christchurch)
            Norm_ClinicSite = 1 - MeanClinicSite;
        elseif contains(raw{i,1},Dunedin)
            Norm_ClinicSite = 2 - MeanClinicSite;
        end
        
        if isnan(raw{i,23}) ~= 1
            Norm_TrailsARaw = raw{i,23} - MeanTrailsARaw;
        elseif isnan(raw{i,23}) == 1
            Norm_TrailsARaw = NaN;
         end
        if isnan(raw{i,24}) ~= 1
            Norm_TrailsAZ = raw{i,24} - MeanTrailsAZ;
        elseif isnan(raw{i,24}) == 1
            Norm_TrailsAZ = NaN;
         end
        if isnan(raw{i,25}) ~= 1
            Norm_TrailsBRaw = raw{i,25} - MeanTrailsBRaw;
        elseif isnan(raw{i,25}) == 1
            Norm_TrailsBRaw = NaN;
         end
        if isnan(raw{i,26}) ~= 1
            Norm_TrailsBZ = raw{i,26} - MeanTrailsBZ;
        elseif isnan(raw{i,26}) == 1
            Norm_TrailsBZ = NaN;
         end
        if isnan(raw{i,27}) ~= 1
            Norm_ColorNamingRaw = raw{i,27} - MeanColorNamingRaw;
        elseif isnan(raw{i,27}) == 1
            Norm_ColorNamingRaw = NaN;
         end
        if isnan(raw{i,28}) ~= 1
            Norm_ColorNamingZ = raw{i,28} - MeanColorNamingZ;
        elseif isnan(raw{i,28}) == 1
            Norm_ColorNamingZ = NaN;
         end
        if isnan(raw{i,30}) ~= 1
            Norm_WordReadingRaw = raw{i,30} - MeanWordReadingRaw;
        elseif isnan(raw{i,30}) == 1
            Norm_WordReadingRaw = NaN;
         end
        if isnan(raw{i,31}) ~= 1
            Norm_WordReadingZ = raw{i,31} - MeanWordReadingZ;
        elseif isnan(raw{i,31}) == 1
            Norm_WordReadingZ = NaN;
         end
        if isnan(raw{i,33}) ~= 1
            Norm_InhibitionRaw = raw{i,33} - MeanInhibitionRaw;
        elseif isnan(raw{i,33}) == 1
            Norm_InhibitionRaw = NaN;
         end
        if isnan(raw{i,34}) ~= 1
            Norm_InhibitionZ = raw{i,34} - MeanInhibitionZ;
        elseif isnan(raw{i,34}) == 1
            Norm_InhibitionZ = NaN;
         end
        if isnan(raw{i,39}) ~= 1
            Norm_LetFluencyRaw = raw{i,39} - MeanLetFluencyRaw;
        elseif isnan(raw{i,39}) == 1
            Norm_LetFluencyRaw = NaN;
         end
        if isnan(raw{i,40}) ~= 1
            Norm_LetFluencyZ = raw{i,40} - MeanLetFluencyZ;
        elseif isnan(raw{i,40}) == 1
            Norm_LetFluencyZ = NaN;
         end
        if isnan(raw{i,42}) ~= 1
            Norm_CatFluencyRaw = raw{i,42} - MeanCatFluencyRaw;
        elseif isnan(raw{i,42}) == 1
            Norm_CatFluencyRaw = NaN;
         end
        if isnan(raw{i,43}) ~= 1
            Norm_CatFluencyZ = raw{i,43} - MeanCatFluencyZ;
        elseif isnan(raw{i,43}) == 1
            Norm_CatFluencyZ = NaN;
         end
        if isnan(raw{i,48}) ~= 1
            Norm_SwitchingRaw = raw{i,48} - MeanSwitchingRaw;
        elseif isnan(raw{i,48}) == 1
            Norm_SwitchingRaw = NaN;
         end
        if isnan(raw{i,49}) ~= 1
            Norm_SwitchingZ = raw{i,49} - MeanSwitchingZ;
        elseif isnan(raw{i,49}) == 1
            Norm_SwitchingZ = NaN;
         end
        if isnan(raw{i,51}) ~= 1
            Norm_HayBTime1Raw = raw{i,51} - MeanHayBTime1Raw;
        elseif isnan(raw{i,51}) == 1
            Norm_HayBTime1Raw = NaN;
         end
        if isnan(raw{i,52}) ~= 1
            Norm_HayBTime1Z = raw{i,52} - MeanHayBTime1Z;
        elseif isnan(raw{i,52}) == 1
            Norm_HayBTime1Z = NaN;
         end
        if isnan(raw{i,53}) ~= 1
            Norm_HayBTime2Raw = raw{i,53} - MeanHayBTime2Raw;
        elseif isnan(raw{i,53}) == 1
            Norm_HayBTime2Raw = NaN;
         end
        if isnan(raw{i,54}) ~= 1
            Norm_HayBTime2Z = raw{i,54} - MeanHayBTime2Z;
        elseif isnan(raw{i,54}) == 1
            Norm_HayBTime2Z = NaN;
         end
        if isnan(raw{i,55}) ~= 1
            Norm_HayBCatARaw = raw{i,55} - MeanHayBCatARaw;
        elseif isnan(raw{i,55}) == 1
            Norm_HayBCatARaw = NaN;
         end
        if isnan(raw{i,56}) ~= 1
            Norm_HayBCatAZ = raw{i,56} - MeanHayBCatAZ;
        elseif isnan(raw{i,56}) == 1
            Norm_HayBCatAZ = NaN;
         end
        if isnan(raw{i,57}) ~= 1
            Norm_HayBCatBRaw = raw{i,57} - MeanHayBCatBRaw;
        elseif isnan(raw{i,57}) == 1
            Norm_HayBCatBRaw = NaN;
         end
        if isnan(raw{i,58}) ~= 1
            Norm_HayBCatBZ = raw{i,58} - MeanHayBCatBZ;
        elseif isnan(raw{i,58}) == 1
            Norm_HayBCatBZ = NaN;
         end
        
        %Norm_ACE = raw{i,10} - MeanACE;
        %fprintf(fid, '%s %.2f', matrix_line, Norm_HayBCatBZ);
        fprintf(fid, '%s %.2f', matrix_line, Norm_HayBCatBZ);
        %fprintf(fid, '%.2f %.2f %.2f', Norm_Group, Norm_Age, Norm_Sex);
        %fprintf(fid, '%s %.2f %.2f %.2f', matrix_line, Norm_Age, Norm_Sex, Norm_ClinicSite);
        fprintf(fid, '\n');
    end
end
fclose(fid);


%create your associated contrast matrix file
fid2 = fopen('contrast_matrix_overall_correlation.txt', 'w');
if (fid2 == -1)
    disp('Error in creating the text file')
else
    
    Cov1Pos = '0 1';
    Cov1Neg = '0 -1';
    
    %with covariates
    %Cov1Pos_covar-age = '0 1 0';
    %Cov1Neg_covar-age = '0 -1 0';
    
    fprintf(fid2, '%s', Cov1Pos);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Cov1Neg);
    
    fclose(fid2);
end















%-------------for by-group interaction correlation matrices---------------%

fid = fopen('design_matrix_overall_correlation_test_ACE.txt', 'w');
if (fid == -1)
    disp('Error in creating the text file')
else
    for i = 2:length(txt)
        if (raw{i,5} == 1)
            matrix_line = '1 0 0 0 0';
        elseif (raw{i,5} == 2)
            matrix_line = '0 1 0 0 0';
        elseif (raw{i,5} == 3)
            matrix_line = '0 0 1 0 0';
        elseif (raw{i,5} == 4)
            matrix_line = '0 0 0 1 0';
        elseif (raw{i,5} == 5)
            matrix_line = '0 0 0 0 1';
        end
        Norm_Age = raw{i,5} - MeanAge;
        Norm_Sex = raw{i,7} - MeanSex;
        Norm_ACE = raw{i,59} - MeanACE;
        %Norm_Group = raw{i,9} - MeanGroup;
        if contains(raw{i,1},Auckland)
            Norm_ClinicSite = 0 - MeanClinicSite;
        elseif contains(raw{i,1},Christchurch)
            Norm_ClinicSite = 1 - MeanClinicSite;
        elseif contains(raw{i,1},Dunedin)
            Norm_ClinicSite = 2 - MeanClinicSite;
        end
        
        if isnan(raw{i,23}) ~= 1
            Norm_TrailsARaw = raw{i,23} - MeanTrailsARaw;
        elseif isnan(raw{i,23}) == 1
            Norm_TrailsARaw = [];
        end
        if isnan(raw{i,24}) ~= 1
            Norm_TrailsAZ = raw{i,24} - MeanTrailsAZ;
        elseif isnan(raw{i,24}) == 1
            Norm_TrailsAZ = [];
        end
        if isnan(raw{i,25}) ~= 1
            Norm_TrailsBRaw = raw{i,25} - MeanTrailsBRaw;
        elseif isnan(raw{i,25}) == 1
            Norm_TrailsBRaw = [];
        end
        if isnan(raw{i,26}) ~= 1
            Norm_TrailsBZ = raw{i,26} - MeanTrailsBZ;
        elseif isnan(raw{i,26}) == 1
            Norm_TrailsBZ = [];
        end
        if isnan(raw{i,27}) ~= 1
            Norm_ColorNamingRaw = raw{i,27} - MeanColorNamingRaw;
        elseif isnan(raw{i,27}) == 1
            Norm_ColorNamingRaw = [];
        end
        if isnan(raw{i,28}) ~= 1
            Norm_ColorNamingZ = raw{i,28} - MeanColorNamingZ;
        elseif isnan(raw{i,28}) == 1
            Norm_ColorNamingZ = [];
        end
        if isnan(raw{i,30}) ~= 1
            Norm_WordReadingRaw = raw{i,30} - MeanWordReadingRaw;
        elseif isnan(raw{i,30}) == 1
            Norm_WordReadingRaw = [];
        end
        if isnan(raw{i,31}) ~= 1
            Norm_WordReadingZ = raw{i,31} - MeanWordReadingZ;
        elseif isnan(raw{i,31}) == 1
            Norm_WordReadingZ = [];
        end
        if isnan(raw{i,33}) ~= 1
            Norm_InhibitionRaw = raw{i,33} - MeanInhibitionRaw;
        elseif isnan(raw{i,33}) == 1
            Norm_InhibitionRaw = [];
        end
        if isnan(raw{i,34}) ~= 1
            Norm_InhibitionZ = raw{i,34} - MeanInhibitionZ;
        elseif isnan(raw{i,34}) == 1
            Norm_InhibitionZ = [];
        end
        if isnan(raw{i,39}) ~= 1
            Norm_LetFluencyRaw = raw{i,39} - MeanLetFluencyRaw;
        elseif isnan(raw{i,39}) == 1
            Norm_LetFluencyRaw = [];
        end
        if isnan(raw{i,40}) ~= 1
            Norm_LetFluencyZ = raw{i,40} - MeanLetFluencyZ;
        elseif isnan(raw{i,40}) == 1
            Norm_LetFluencyZ = [];
        end
        if isnan(raw{i,42}) ~= 1
            Norm_CatFluencyRaw = raw{i,42} - MeanCatFluencyRaw;
        elseif isnan(raw{i,42}) == 1
            Norm_CatFluencyRaw = [];
        end
        if isnan(raw{i,43}) ~= 1
            Norm_CatFluencyZ = raw{i,43} - MeanCatFluencyZ;
        elseif isnan(raw{i,43}) == 1
            Norm_CatFluencyZ = [];
        end
        if isnan(raw{i,48}) ~= 1
            Norm_SwitchingRaw = raw{i,48} - MeanSwitchingRaw;
        elseif isnan(raw{i,48}) == 1
            Norm_SwitchingRaw = [];
        end
        if isnan(raw{i,49}) ~= 1
            Norm_SwitchingZ = raw{i,49} - MeanSwitchingZ;
        elseif isnan(raw{i,49}) == 1
            Norm_SwitchingZ = [];
        end
        if isnan(raw{i,51}) ~= 1
            Norm_HayBTime1Raw = raw{i,51} - MeanHayBTime1Raw;
        elseif isnan(raw{i,51}) == 1
            Norm_HayBTime1Raw = [];
        end
        if isnan(raw{i,52}) ~= 1
            Norm_HayBTime1Z = raw{i,52} - MeanHayBTime1Z;
        elseif isnan(raw{i,52}) == 1
            Norm_HayBTime1Z = [];
        end
        if isnan(raw{i,53}) ~= 1
            Norm_HayBTime2Raw = raw{i,53} - MeanHayBTime2Raw;
        elseif isnan(raw{i,53}) == 1
            Norm_HayBTime2Raw = [];
        end
        if isnan(raw{i,54}) ~= 1
            Norm_HayBTime2Z = raw{i,54} - MeanHayBTime2Z;
        elseif isnan(raw{i,54}) == 1
            Norm_HayBTime2Z = [];
        end
        if isnan(raw{i,55}) ~= 1
            Norm_HayBCatARaw = raw{i,55} - MeanHayBCatARaw;
        elseif isnan(raw{i,55}) == 1
            Norm_HayBCatARaw = [];
        end
        if isnan(raw{i,56}) ~= 1
            Norm_HayBCatAZ = raw{i,56} - MeanHayBCatAZ;
        elseif isnan(raw{i,56}) == 1
            Norm_HayBCatAZ = [];
        end
        if isnan(raw{i,57}) ~= 1
            Norm_HayBCatBRaw = raw{i,57} - MeanHayBCatBRaw;
        elseif isnan(raw{i,57}) == 1
            Norm_HayBCatBRaw = [];
        end
        if isnan(raw{i,58}) ~= 1
            Norm_HayBCatBZ = raw{i,58} - MeanHayBCatBZ;
        elseif isnan(raw{i,58}) == 1
            Norm_HayBCatBZ = [];
        end
        
        %Norm_ACE = raw{i,10} - MeanACE;
        %fprintf(fid, '%s %.2f', matrix_line, Norm_HayBCatBZ);
        fprintf(fid, '%s %.2f', matrix_line, Norm_ACE);
        %fprintf(fid, '%.2f %.2f %.2f', Norm_Group, Norm_Age, Norm_Sex);
        %fprintf(fid, '%s %.2f %.2f %.2f', matrix_line, Norm_Age, Norm_Sex, Norm_ClinicSite);
        fprintf(fid, '\n');
    end
end
fclose(fid);


%create your associated contrast matrix file
fid2 = fopen('contrast_matrix_overall_correlation.txt', 'w');
if (fid2 == -1)
    disp('Error in creating the text file')
else
    
    Cov1Pos = '0 0 0 0 0 1';
    Cov1Neg = '0 0 0 0 0 -1';
    
    %with covariates
    %Cov1Pos_covar-age = '0 0 0 0 0 1 0';
    %Cov1Neg_covar-age = '0 0 0 0 0 -1 0';
    
    fprintf(fid2, '%s', Cov1Pos);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Cov1Neg);
    
    fclose(fid2);
end









end

