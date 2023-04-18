clc;
clear all;

%cd V:\Archive\NECTAR_data\LENORE\derivatives\groups\F2\diff_data\longitudinal\connectome\hcpmmpFiles\weighted\interaction_study(F2-F0)\groups\ANOVA_Test
cd V:\Archive\NECTAR_data\LENORE\derivatives\groups\F2\diff_data\longitudinal\connectome\dkFiles\weighted\fs_default_ordered\interaction_study(F2-F0)_C&SCD_thresholded\groups\ANOVA_test

%no covariates
load p_vals_post-pre.mat
load p_vals_post-pre_linear_trend.mat
%with covariates (age and sex)
load p_vals_post-pre_covariates.mat
load p_vals_post-pre_covariates_linear_trend.mat
%matrix = unnamed;

%convert table matrix
matrix=table2array(matrix);
%get lower triangle half of matrix
matrix = tril(matrix); 
[rows,cols,vals]=find(matrix>=.95);

cd ../

%get controls
cd C
matrices=dir('dk*.csv');
C_values = [];
for i = 1:length(cols) %length(cols)
    this_edges_C_values = [];
    for j = 1:length(matrices)
        subj_matrix=load(matrices(j).name);
        this_edges_C_values = [this_edges_C_values; subj_matrix(rows(i), cols(i))];
    end
    C_values = [C_values this_edges_C_values];
end

cd ../

%get SCD
cd SCD

matrices=dir('dk*.csv');
SCD_values = [];
for i = 1:length(cols) %length(cols)
    this_edges_SCD_values = [];
    for j = 1:length(matrices)
        subj_matrix=load(matrices(j).name);
        this_edges_SCD_values = [this_edges_SCD_values; subj_matrix(rows(i), cols(i))];
    end
    SCD_values = [SCD_values this_edges_SCD_values];
end

 
 cd ../
 %get aMCI
  cd aMCI
  
matrices=dir('dk*.csv');
aMCI_values = [];
for i = 1:length(cols) %length(cols)
    this_edges_aMCI_values = [];
    for j = 1:length(matrices)
        subj_matrix=load(matrices(j).name);
        this_edges_aMCI_values = [this_edges_aMCI_values; subj_matrix(rows(i), cols(i))];
    end
    aMCI_values = [aMCI_values this_edges_aMCI_values];
end
 cd ../
  
 %get mMCI
  cd mMCI
 
matrices=dir('dk*.csv');
mMCI_values = [];
for i = 1:length(cols) %length(cols)
    this_edges_mMCI_values = [];
    for j = 1:length(matrices)
        subj_matrix=load(matrices(j).name);
        this_edges_mMCI_values = [this_edges_mMCI_values; subj_matrix(rows(i), cols(i))];
    end
    mMCI_values = [mMCI_values this_edges_mMCI_values];
end

 
 cd ../
   %get AD
  cd AD

matrices=dir('dk*.csv');
AD_values = [];
for i = 1:length(cols) %length(cols)
    this_edges_AD_values = [];
    for j = 1:length(matrices)
        subj_matrix=load(matrices(j).name);
        this_edges_AD_values = [this_edges_AD_values; subj_matrix(rows(i), cols(i))];
    end
    AD_values = [AD_values this_edges_AD_values];
end
 

 
%display figure
 %figure
 %bar([mean(this_edges_C_values) mean(this_edges_SCD_values) mean(this_edges_aMCI_values) mean(this_edges_mMCI_values) mean(this_edges_AD_values)]) 
 
 
 %display all sig. edges with interaction(e.g., 11, 78 or 121 figures) (write loop)
 for i = 1:length(cols)
      %for i = 1:2

     figure
     bar([mean(C_values(:,i)) mean(SCD_values(:,i)) mean(aMCI_values(:,i)) mean(mMCI_values(:,i)) mean(AD_values(:,i))])
 end
 
 %display the mean of all edges per group
 figure
 bar([mean(C_values,'all') mean(SCD_values,'all') mean(aMCI_values,'all') mean(mMCI_values,'all') mean(AD_values,'all')])
 
 %get mean across columns 
 %e.g,. mean(C_values,2)
 
 x = mean(AD_values,2)
 
 
 
 
 
 
 
 
 %for FPN:
 %no covariates
%just need to extract one value - row 2819
 
 cd V:\Archive\NECTAR_data\LENORE\derivatives\groups\F2\diff_data\longitudinal\connectome\FPN_Files\weighted\interaction_study(F2-F0)\groups

%get controls
cd C
matrices=dir('FPN*.csv');
this_participants_C_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_C_values = [this_participants_C_values; subj_matrix(413,1)];
 end 
    
 cd ../
 %get SCD
cd SCD
matrices=dir('FPN*.csv');
this_participants_SCD_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_SCD_values = [this_participants_SCD_values; subj_matrix(413,1)];
 end
 
 cd ../
 %get aMCI
cd aMCI
matrices=dir('FPN*.csv');
this_participants_aMCI_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_aMCI_values = [this_participants_aMCI_values; subj_matrix(413,1)];
 end 
 
 cd ../
 %get mMCI
cd mMCI
matrices=dir('FPN*.csv');
this_participants_mMCI_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_mMCI_values = [this_participants_mMCI_values; subj_matrix(413,1)];
 end 
 
  cd ../
 %get AD
cd AD
matrices=dir('FPN*.csv');
this_participants_AD_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_AD_values = [this_participants_AD_values; subj_matrix(413,1)];
 end
 
 %plot bar graph for FPN
 figure
 bar([mean(this_participants_C_values,'all') mean(this_participants_SCD_values,'all') mean(this_participants_aMCI_values,'all') mean(this_participants_mMCI_values,'all') mean(this_participants_AD_values,'all')])
 

 
 
 %for FPN:
 %with covariates age and sex
%just need to extract two values - row 1632 and 1945 (take the avg. of
%these two values)
 
 cd V:\Archive\NECTAR_data\LENORE\derivatives\groups\F2\diff_data\longitudinal\connectome\FPN_Files\weighted\interaction_study(F2-F0)\groups

%get controls
cd C
matrices=dir('FPN*.csv');
this_participants_C_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_C_values = [this_participants_C_values; (subj_matrix(1945) + subj_matrix(1632) + subj_matrix(1184))/3];
 end 
    
 cd ../
 %get SCD
cd SCD
matrices=dir('FPN*.csv');
this_participants_SCD_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_SCD_values = [this_participants_SCD_values; (subj_matrix(1945) + subj_matrix(1632)+ subj_matrix(1184))/3];
 end
 
 cd ../
 %get aMCI
cd aMCI
matrices=dir('FPN*.csv');
this_participants_aMCI_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_aMCI_values = [this_participants_aMCI_values; (subj_matrix(1945) + subj_matrix(1632)+ subj_matrix(1184))/3];
 end 
 
 cd ../
 %get mMCI
cd mMCI
matrices=dir('FPN*.csv');
this_participants_mMCI_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_mMCI_values = [this_participants_mMCI_values; (subj_matrix(1945) + subj_matrix(1632)+ subj_matrix(1184))/3];
 end 
 
  cd ../
 %get AD
cd AD
matrices=dir('FPN*.csv');
this_participants_AD_values = [];
 for i = 1:length(matrices)
      subj_matrix=load(matrices(i).name);
      this_participants_AD_values = [this_participants_AD_values; (subj_matrix(1945) + subj_matrix(1632)+ subj_matrix(1184))/3];
 end
 
 %plot bar graph for FPN
 figure
 bar([mean(this_participants_C_values,'all') mean(this_participants_SCD_values,'all') mean(this_participants_aMCI_values,'all') mean(this_participants_mMCI_values,'all') mean(this_participants_AD_values,'all')])
 

 
 
 