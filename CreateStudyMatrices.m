function CreateStudyMatrices(excelFile,fileLocation, startdir, groupname)
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

%      Inputs (4): excelFile - Name of the excel file containing 
%                               particpant info. 
%                  fileLocation - directory of where the excelFile is
%                                 contained. 
%                  startdir - directory where your source data is
%                              stored
%                  groupname - name of the group/study that you are
%                              analysing
                           
%       Outputs (none): You writing the matrices to a text file. 


%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 05/09/20

%navigate to the most basic file directory and then to the file location
cd /;
cd ([fileLocation]);

%read in excel file into matlab format
[num, txt, raw] = xlsread(excelFile, 'A:K');
  
%calculate average scores to input as covariates
sumAge = 0;
sumSex = 0;
sumACE = 0;
NumParticipants = 0;
for i = 2:length(txt)
    if raw{i,9} ~= 0
        sumAge = raw{i,5} + sumAge;
        sumSex = raw{i,7} + sumSex;
        sumACE = raw{i,10} + sumACE; 
        NumParticipants = NumParticipants + 1;
    end 
end 
MeanAge = sumAge / NumParticipants;
MeanSex = sumSex / NumParticipants;
MeanACE = sumACE / NumParticipants;
 
%create your matrices for the study
 
%Create design matrix with the covariates. Added covariates currently are: 
%age, sex,and overall ACE-III score. 
fid = fopen('design_matrix.txt', 'w');
if (fid == -1)
    disp('Error in creating the text file')
else
    for i = 2:length(txt)
        if raw{i,9} ~= 0
            classification = raw{i,9};
            if (classification == 1)
                matrix_line = '1 0 0 0 0 ';
            elseif (classification == 2)
                matrix_line = '0 1 0 0 0 ';
            elseif (classification == 3)
                matrix_line = '0 0 1 0 0 ';
            elseif (classification == 4)
                matrix_line = '0 0 0 1 0 ';
            elseif (classification == 5)
                matrix_line = '0 0 0 0 1 ';
            end
            Norm_Age = raw{i,5} - MeanAge;
            Norm_Sex = raw{i,7} - MeanSex;
            Norm_ACE = raw{i,10} - MeanACE;
            fprintf(fid, '%s %.2f %.2f %.2f', matrix_line, Norm_Age, Norm_Sex, Norm_ACE);
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
    Group1 = '1 0 0 0 0 0 0 0';
    Group2 = '0 1 0 0 0 0 0 0';
    Group3 = '0 0 1 0 0 0 0 0';
    Group4 = '0 0 0 1 0 0 0 0';
    Group5 = '0 0 0 0 1 0 0 0';
    Cov1Pos = '0 0 0 0 0 1 0 0';
    Cov1Neg = '0 0 0 0 0 -1 0 0';
    Cov2Pos = '0 0 0 0 0 0 1 0';
    Cov2Neg = '0 0 0 0 0 0 -1 0';
    Cov3Pos = '0 0 0 0 0 0 0 1';
    Cov3Neg = '0 0 0 0 0 0 0 -1';
    
    fprintf(fid2, '%s', Group1);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Group2);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Group3);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Group4);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Group5);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Cov1Pos);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Cov1Neg);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Cov2Pos);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Cov2Neg);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Cov3Pos);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s', Cov3Neg);
    
    fclose(fid2);
end 
    
%move both your matrices files into the stats_matrices folder. 
movefile('design_matrix.txt', [startdir, '/derivatives/diff_data/' groupname, '/stats_matrices']);
movefile('contrast_matrix.txt', [startdir, '/derivatives/diff_data/' groupname, '/stats_matrices']);



end

