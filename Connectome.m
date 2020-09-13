%This script will create a connectome to represent the number of
%streamlines connecting different parts of the brain. This script assumes 
%that the user has ran the previous pipelines (preprocessing and CSD), and 
%it will use that processed data for this pipeline. The tractography inputs 
%for this will be with using ACT and SIFT. The brain will be parcellated 
%into different regions or nodes using an atlas. 

%See MRtrix documentation: 
%https://mrtrix.readthedocs.io/en/latest/quantitative_structural_connectivity/structural_connectome.html

% Steps: 
% 1. Create a 5tt image for ACT
%   a. Perform bias field correction on t1w and t2 FLAIR images
%   b. Co-registration t1w and dwi image
%   c. Generate a '4tt image'
%   d. Edit in a pathological image to create a 5ttimage
%   e. Create a 'seed' boundary between the white and grey matter for
%   streamline generation
% 2. Create streamlines for tractography with ACT
% 3. Refine streamlines with SIFT
% 4. Parcellate brain regions/nodes with an atlas (e.g. FreeSurfer)
% 5. Create the connectome
%   a. Convert the labels to MRtrix format
%   b. Co-register the parcellation to the grey matter boundary
%   c. Create the whole-brain connectome


clc;
clear all;
close all;


%define/add pathways
%startdir = input('Please enter data directory:', 's');
startdir = '/data/USERS/LENORE';

%Script directory is defined, so that it can be added to path below:
%ScriptDirectory = input('Please enter script directory:', 's');
ScriptDirectory = '/data/USERS/LENORE/scripts/dprc/diffusion';

%should be the same groupname from what the user analysed in the previous scripts (e.g. CSD.m).  
groupname = input('Which pre-processed group / study do you want to continue to analyse?: ', 's');

%create a connectome folder to place all of the data and analysis in
mkdir([startdir,'/derivatives/diff_data/', groupname, '/connectome/']);

%go into the group/study analysis folder
cd([startdir '/derivatives/diff_data/', groupname]);

%copy necessary files into the connectome folder
%dwi preprocessed folder
preprocessed_dwi = fullfile('preprocessed_dwi');
connectome = fullfile('connectome');
copyfile(preprocessed_dwi, fullfile(connectome, 'preprocessed_dwi'));
%wmfod bias field corrected and intensity normalised data
copyfile ('wmfod_norm_*', [startdir,'/derivatives/diff_data/', groupname, '/connectome/']);


%Add your script and all necessary files (e.g. data, functions) to path
addpath(genpath(startdir));
addpath(genpath(ScriptDirectory));

%go into the connectome folder
cd([startdir '/derivatives/diff_data/', groupname, '/connectome']);

%define variables
participants = dir(fullfile('preprocessed_dwi', '*.mif'));
datafile = '_acq_data_dwi';

%create 5ttImageCheck text file with header line
cd([startdir '/derivatives/diff_data/', groupname, '/qc/']);
fid4 = fopen('5ttImageCheck.txt', 'w');
if (fid4 == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid4, '%s     %s', 'Participant', '5ttImage_status');
    fclose(fid4);
end




%-------------------------------------------------------------------------%
%Step 1: Create a 5tt image for ACT
for i = 1:length(participants)
    
    full_name = participants(i).name;
    PAR_NAME = full_name(1:15);
    
    %copy participant's raw t1w and t2 FLAIR images, and the preprocessed/segmented
    %t2 FLAIR images into main group folder
    copyfile ([startdir '/sourcedata/' PAR_NAME '/anat/', '*.nii'], [startdir,'/derivatives/diff_data/', groupname, '/connectome']);
    
    %a) perform bias field correct t1w and t2 FLAIR images
    %get brain mask from both images
    unix(['bet ' PAR_NAME, '_T1w.nii bet_' PAR_NAME, '_T1w.nii -m -f 0.2']);
    unix(['bet ' PAR_NAME, '_FLAIR.nii bet_' PAR_NAME, '_FLAIR.nii -m -f 0.2']);
    
    %bias field correct directly with ants programme
    unix(['N4BiasFieldCorrection -d 3 -x bet_' PAR_NAME '_T1w_mask.nii.gz -i ' PAR_NAME, '_T1w.nii -o ' PAR_NAME '_bfc_T1w.nii']);
    unix(['N4BiasFieldCorrection -d 3 -x bet_' PAR_NAME '_FLAIR_mask.nii.gz -i ' PAR_NAME, '_FLAIR.nii -o ' PAR_NAME '_bfc_FLAIR.nii']);
    
    %b) co-registration of t1w and t2 FLAIR to dwi image through fsl FLIRT
    %convert to nii format to run fsl's FLIRT - use the best B0 volume,
    %which will be the first vol in the dwi sequence (vol 0).
    unix(['mrconvert -coord 3 0 preprocessed_dwi/', PAR_NAME, datafile,'.mif ref_b0_', PAR_NAME, '.nii']);
    
    %linear registration with 6 dof and the tranformation matrix output
    unix(['flirt -in ', PAR_NAME, '_bfc_T1w.nii -ref ref_b0_', PAR_NAME, '.nii -dof 6 -out t1_flirt_' PAR_NAME, '.nii -omat transform_flirt_t12dwi.mat']);
    unix(['flirt -in ', PAR_NAME, '_bfc_FLAIR.nii -ref ref_b0_', PAR_NAME, '.nii -dof 6 -out t2_flirt_' PAR_NAME, '.nii -omat transform_flirt_t22dwi.mat']);
    
    %apply the linear transformation to the t1 and t2 image
    unix(['transformconvert transform_flirt_t12dwi.mat ' PAR_NAME, '_bfc_T1w.nii ref_b0_', PAR_NAME, '.nii flirt_import transform_mrtrix_t12dwi.txt']);
    unix(['transformconvert transform_flirt_t22dwi.mat ' PAR_NAME, '_bfc_FLAIR.nii ref_b0_', PAR_NAME, '.nii flirt_import transform_mrtrix_t22dwi.txt']);
    unix(['mrtransform ' PAR_NAME, '_bfc_T1w.nii -linear transform_mrtrix_t12dwi.txt ' PAR_NAME, '_T1coreg.nii']);
    unix(['mrtransform ' PAR_NAME, '_bfc_FLAIR.nii -linear transform_mrtrix_t12dwi.txt -template ' PAR_NAME, '_T1coreg.nii ' PAR_NAME, '_T2coreg.nii']);
    
    %convert to .mif (mrtrix) format
    unix(['mrconvert ' PAR_NAME, '_T1coreg.nii ' PAR_NAME '_T1coreg.mif']);
    unix(['mrconvert ' PAR_NAME, '_T2coreg.nii ' PAR_NAME '_T2coreg.mif']);
    
    %c) generate 5tt image with brain mask, t1, and t2 FLAIR
    unix(['5ttgen fsl -mask brain_mask/', PAR_NAME, datafile, '.mif -t2 ' PAR_NAME, '_T2coreg.mif ' PAR_NAME '_T1coreg.mif 4ttimage_' PAR_NAME '.mif']);
    
    %d) edit in the pathological tissue (WMH) to the 5tt image:
    %convert WMH lesion mask to .mif file
    unix(['mrconvert ples_lpa_mr', PAR_NAME, '_FLAIR.nii.gz ples_lpa_mr', PAR_NAME, '_FLAIR.mif']);
    
    %reslice the WMH lesion mask to match the 4tt image
    unix(['mrtransform ples_lpa_mr' PAR_NAME, '_FLAIR.mif -template 4ttimage_' PAR_NAME, '.mif ' PAR_NAME, '_WMH_mask_transformed.mif']);
    
    %add the image into the 5tt image
    unix(['5ttedit -path ', PAR_NAME, '_WMH_mask_transformed.mif 4ttimage_' PAR_NAME, '.mif 5ttimage_' PAR_NAME '.mif']);
    
    %check that your 5ttimage conforms to MRtrix's expected format
    FiveTTImageCheck(PAR_NAME);

    %e) Create a 'seed' boundary between the white and grey matter for
    %streamline generation
    unix(['5tt2gmwmi 5ttimage_' PAR_NAME '.mif gmwmSeed_mask_' PAR_NAME '.mif']);
   
   
    %---------------------------------------------------------------------%
    %Step 2: Generate streamline tracks with -act
    unix(['tckgen -act 5ttimage_' PAR_NAME '.mif -backtrack -seed_gmwmi gmwmSeed_mask_' PAR_NAME '.mif -maxlength 250 -cutoff 0.06 -select 10000000 wmfod_norm_' PAR_NAME '.mif tracks_10M_' PAR_NAME '.tck']);
    
    %if you want to view the tracks, make a smaller version of them:
    unix(['tckedit tracks_10M_' PAR_NAME '.tck -number 200k smallerTracks_200k_' PAR_NAME '.tck']);
    
   
    %---------------------------------------------------------------------%
    %Step 3: Refine streamlines with tcksift2
    unix(['tcksift2 -act 5ttimage_' PAR_NAME '.mif -out_mu sift_mu_' PAR_NAME '.txt -out_coeffs sift_coeffs_' PAR_NAME '.txt tracks_10M_' PAR_NAME '.tck wmfod_norm_' PAR_NAME '.mif sift_1M_' PAR_NAME '.txt']);

%     %if you want to see the .tck SIFT file, you can generate the image using tcksift 
%     unix(['tcksift -act 5ttimage_' PAR_NAME '.mif -term_number 1000000 tracks_10M_' PAR_NAME '.tck wmfod_norm_' PAR_NAME '.mif sift_1M_' PAR_NAME '.tck']);
    
%     %Compare the tracks with and without SIFT with 10k tracks (for
%     %visualation only)
%     unix(['tckedit tracks_10M_' PAR_NAME '.tck -number 10k SuperSmallTracks_10k_' PAR_NAME '.tck']);
%     unix(['tckedit sift_1M_' PAR_NAME '.tck -number 10k SuperSmallTracks_10k_' PAR_NAME '_sift.tck -force']);
% 
%     
    %---------------------------------------------------------------------%
    %Step 4: Parcellate brain regions/nodes with FreeSurfer
%     %FreeSurfer command below - run this command if you are not using the
%     %output from the fmriprep pipeline
%     unix(['recon-all -i ' PAR_NAME '_T1w.nii -s ' PAR_NAME ' -all']);
%     
      %get necessary freesurfer output from fmriprep


    %---------------------------------------------------------------------%
    %Step 5: Create the connectome
    
    %a) Convert labels to mrtrix format
    unix(['labelconvert ' PAR_NAME '/mri/aparc+aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt /data/SOFTWARE/mrtrix3/share/mrtrix3/labelconvert/fs_default.txt ' PAR_NAME '_parcels.mif']);
    
    %b) Co-register the parcellation to the grey matter boundary
    unix(['mrtransform ' PAR_NAME '_parcels.mif -interp nearest -linear transform_mrtrix_t12dwi.txt -inverse -datatype uint32 ' PAR_NAME '_parcels_coreg.mif']);
    
    %c) Create whole-brain connectome
    unix(['tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift_1M_' PAR_NAME '.txt tracks_10M_' PAR_NAME '.tck ' PAR_NAME '_parcels_coreg.mif ' PAR_NAME '_parcels_coreg.csv -out_assignment assignments_' PAR_NAME '_parcels_coreg.csv']);
    
%     %visualise the connectome matrix in matlab:
%     %cannot view on dementia vm, as their is no openGL for viewing gui. 
%     connectome = importdata([PAR_NAME '_parcels_coreg.csv']);
%     imagesc(connectome)




    
    
    
    
    
end