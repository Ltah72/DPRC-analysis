%Use this script for participants who have failed 5ttgen (5 tissue-type
%parcellation). This will reorient the coregistered anatomical files (T1w +
%T2 FLAIR) to the diffusion file, using 'fslreorient2std'. From having run
%this script, you can then continue on after the 5ttgen function, in the
%main Connectome script to continue the analysis.  

%From our study sample cohort, I have found that this reorientation would 
%most often have to be done for people who have Alzheimer's disease (AD). 

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 08/01/20


%go into the participant analysed folder + choose participants
cd([startdir '/derivatives/' period, '/diff_data/', groupname, '/IN']);

%choose participants who's file you need to reorient
participants = uipickfiles;
datafile = '_acq_data_dwi';

%go into the connectome folder
cd([startdir '/derivatives/' period, '/diff_data/', groupname, '/connectome/']);

for i = 1:length(participants)
    
    [upper_path, PAR_NAME, ~] = fileparts(participants{1,i});
    
    %reorient both images to the correct, standard orientation
    unix(['fslreorient2std ' PAR_NAME, '_T2coreg.nii reoriented_' PAR_NAME, '_T2coreg.nii']);
    unix(['fslreorient2std ' PAR_NAME, '_T1coreg.nii reoriented_' PAR_NAME, '_T1coreg.nii']);
    
    %convert to .mif (mrtrix) format
    unix(['mrconvert reoriented_' PAR_NAME, '_T1coreg.nii.gz reoriented_' PAR_NAME '_T1coreg.mif']);
    unix(['flirt -in ' PAR_NAME '_task-rest_run-1_boldref.nii.gz -ref ref_b0_' PAR_NAME '.nii -out bold2b0_' PAR_NAME '.mat']);
    unix(['mrconvert reoriented_' PAR_NAME, '_T2coreg.nii.gz reoriented_' PAR_NAME '_T2coreg.mif']);
    
    %re-run 5ttgen with the reoriented anat images
    unix(['5ttgen fsl -mask brain_mask_', PAR_NAME, datafile, '.mif -t2 reoriented_' PAR_NAME, '_T2coreg.mif reoriented_' PAR_NAME '_T1coreg.mif 4ttimage_' PAR_NAME '.mif']);
    
end