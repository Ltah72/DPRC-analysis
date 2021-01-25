function CreateConnectomeTemplate(participants)
%This function will create a Connectome template file for use in order to 
%visualise the comparison across your group nodes. Specifically, the 
%template file will be FreeSurfer's parcellation file (aparc+aseg.mgz). 
%This file will come from the fsaverage directory from FreeSurfer 
%($SUBJECTS_DIR/fsaverage/mri). This function will also call on another 
%function, avg_transmat_T12dwi.m, which will generate an average 
%transformation matrix for the template file. 

%       Inputs(1): participants = list of participants being used for 
%                  analysis          
%       Output(none): You are creating the template file 
%                   (subject_average_hcpmmp1_parcels_coreg.mif).

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 25/01/21


%take and convert the fsaverage aparc+aseg.mgz file into the subject template .mif file. 
copyfile(['$SUBJECTS_DIR/fsaverage/mri/aparc+aseg.mgz'], ['fsaverage_aparc+aseg.mgz']);
unix(['mrconvert -datatype uint32 fsaverage_aparc+aseg.mgz subject_average_hcpmmp1.mif']);

%Replace the random integers of the hcpmmp1.mif file with integers that start at 1 and increase by 1.
unix(['labelconvert subject_average_hcpmmp1.mif /SOFTWARE/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_original.txt /SOFTWARE/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt subject_average_hcpmmp1_parcels_nocoreg.mif']);

%get average of transformation matrix (T1 to dwi matrix) in order to do registration in the next step
avg_transmat_T12dwi(participants);

%Register the ordered atlas-based volumetric parcellation to diffusion space.
unix(['mrtransform subject_average_hcpmmp1_parcels_nocoreg.mif -linear average_matrix_transformMrtrixT12Dwi.txt -inverse -datatype uint32 subject_average_hcpmmp1_parcels_coreg.mif']);


end

