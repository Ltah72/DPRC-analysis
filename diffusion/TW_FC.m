%Conduct track-weighted functional connectivity (TW-FC) - combine a
%structural dwi data with a functional resting-state fMRI data, with a
%specific network (e.g. Frontoparietal Network (FPN).



clc;
clear all;
close all;

%choose time period
period = input('Which time period do you want to analyse, e.g. F0, F2, all, etc?: ', 's');

%ask user for what kind of group/study analysis they will conduct (e.g.
%cross-sectional, F0s, longitudinal, F0 vs. F2, etc.):
groupname = input('Please name the group/study analysis: ', 's');

%Choose directory where the diffusion data is (same level as sourcedata
%folder).
%startdir = input('Please enter data directory:', 's');
%shortcut for debugging purposes:
derivdir = '/data/USERS/LENORE/derivatives';

%Define script directory, so that it can be added to path below:
%ScriptDirectory = input('Please enter script directory:', 's');
ScriptDirectory = '/data/USERS/LENORE/scripts/dprc/diffusion';

%Define connectome directory, where most of the pre-required files lie.
%ConnectomeDir = input('Please enter connectome directory:', 's');
ConnectomeDir = ([derivdir, '/groups/' period '/diff_data/' groupname '/connectome/']);

%Define fmriprep derivatives directory.
%FmriprepDerDir = input('Please enter fmriprep derivatives directory:', 's');
FmriprepDerDir = ([derivdir, '/fmriprepped_data/derivatives/fmriprep/']);

%Define nuisance regressor (NR) removed derivatives directory.
%NRDir = input('Please enter func denoised derivatives directory:', 's');
NRDir = ([derivdir '/fMRI_denoised/']);

%Create TW-FC directory for holding output data
mkdir([derivdir,'/derivatives/' period '/diff_data/', groupname, '/TW-FC/']);

%Add your script and all necessary files (e.g. data, functions) to path
addpath(genpath(derivdir));
addpath(genpath(ScriptDirectory));

%go into your sourcedata folder, where all participant files will be
%there.
cd(['/data/sourcedata/' period]);

%choose participants to run TW-dFC
TWFC_participants = uipickfiles;

%go back into your derivatives folder.
cd([derivdir '/groups/' period '/diff_data/', groupname, '/TW-FC/']);

%copy any necessary files over to the TW-FC directory:
%copy the wmfod template file
copyfile ([derivdir '/groups/' period, '/diff_data/' groupname, '/template/wmfod_template.mif'], [derivdir,'/groups/' period, '/diff_data/', groupname, '/TW-FC/']);


for i = 1:length(TWFC_participants)
    
    [upper_path, PAR_NAME, ~] = fileparts(TWFC_participants{1,i});
    PAR_NAME = PAR_NAME(1:15);
    
    %1. Co-register t1w dwi image through fsl FLIRT
    %a) Linear registration with 6 dof and the transformation matrix output
    unix(['flirt -in ', FmriprepDerDir, PAR_NAME, '/anat/' PAR_NAME, '_space-MNI152NLin6Asym_desc-preproc_T1w.nii.gz -ref ' ConnectomeDir, 'ref_b0_', PAR_NAME, '.nii -dof 6 -out t1_flirt_' PAR_NAME, '.nii -omat transform_flirt_t12dwi_' PAR_NAME '.mat']);
    %b) Apply the linear transformation to the t1 and t2 image
    unix(['transformconvert transform_flirt_t12dwi_' PAR_NAME '.mat ' FmriprepDerDir, PAR_NAME, '/anat/' PAR_NAME, '_space-MNI152NLin6Asym_desc-preproc_T1w.nii.gz ' ConnectomeDir, 'ref_b0_', PAR_NAME, '.nii flirt_import transform_mrtrix_t12dwi_' PAR_NAME '.txt']);
    unix(['mrtransform ' FmriprepDerDir, PAR_NAME, '/anat/' PAR_NAME, '_space-MNI152NLin6Asym_desc-preproc_T1w.nii.gz -linear transform_mrtrix_t12dwi_' PAR_NAME '.txt ' PAR_NAME, '_T1coreg-MNI.nii']);
    %c) Convert the co-registered T1w-dwi nifti to .mif file
    unix(['mrconvert ' PAR_NAME, '_T1coreg-MNI.nii ' PAR_NAME '_T1coreg-MNI.mif -force']);
    
    
    %2. Co-register fMRI-dwi file
    %a) Create transformation matrix
    unix(['flirt -in ' FmriprepDerDir, PAR_NAME '/func/' PAR_NAME '_task-rest_run-1_space-MNI152NLin6Asym_boldref.nii.gz -ref ' ConnectomeDir 'ref_b0_' PAR_NAME '.nii -out flirt_bold_' PAR_NAME '.nii -omat bold2b0_' PAR_NAME '.mat']);
    %b) Generate coregistered image (fmri main bold sequence coregistered to B0):
    unix(['flirt -in ' NRDir, PAR_NAME '_task-rest_run-1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold_NR.nii.gz -applyxfm -init bold2b0_' PAR_NAME '.mat -paddingsize 0.0 -interp trilinear -ref ' ConnectomeDir 'ref_b0_' PAR_NAME '.nii -out bold_coreg-MNI_' PAR_NAME '.nii.gz']);
    %c) Convert the co-registered BOLD-dwi nifti to .mif file
    unix(['mrconvert bold_coreg-MNI_' PAR_NAME '.nii.gz bold_coreg-MNI_' PAR_NAME '.mif']);
    
    
    %3. Warp individual track files to group template
    %create template to subject warp (use non-linear transformation)
    %co-register subject and template
    unix(['mrregister ' ConnectomeDir 'wmfod_norm_' PAR_NAME '.mif wmfod_template.mif -nl_warp_full warp_full_' PAR_NAME '.mif']);
    %create template to subject warp
    unix(['warpconvert warp_full_' PAR_NAME '.mif warpfull2deformation w_t2s_' PAR_NAME '.mif -from 2 -template ' ConnectomeDir 'wmfod_norm_' PAR_NAME '.mif']);
    % warp subject tck to template space
    unix(['tcktransform ' ConnectomeDir 'sift_1M_' PAR_NAME '.tck w_t2s_' PAR_NAME '.mif sift_1M_subject_at_template_' PAR_NAME '.tck']);
    
    
    
    %4. %---------------------TW-FC with specific network---------------------%
    
    
    
    
    %5. %-------------------TW-FC without specific network--------------------%
    
    %a) Create the track-weighted dynamic functional connectivity (TW-dFC) per each individual
    unix(['tckdfc -dynamic rectangle 81 sift_1M_subject_at_template_' PAR_NAME '.tck -template ' PAR_NAME '_T1coreg-MNI.mif bold_coreg-MNI_' PAR_NAME '.mif tdfc_dynamic_' PAR_NAME '.mif']);
    %b) Create summary image per each individual - for 95% confidence interval
    %calculate average image
    unix(['mrmath tdfc_dynamic_' PAR_NAME '.mif mean tdfc_average_' PAR_NAME '.mif -axis 3']);
    %calculate STD (voxels across volumes)
    unix(['mrmath tdfc_dynamic_' PAR_NAME '.mif std tdfc_std_' PAR_NAME '.mif -axis 3']);
    %divide STD image by square root of n to get the SE
    unix(['mrcalc tdfc_std_' PAR_NAME '.mif 490 -sqrt -divide tdfc_SE_' PAR_NAME '.mif']);
    %multiply +1.96 (95% CI) to the SE image (upper margin of error)
    unix(['mrcalc tdfc_SE_' PAR_NAME '.mif 1.96 -mult tdfc_upper_ME_' PAR_NAME '.mif']);
    %add the mean to the upper margin of error to get the upper bound
    unix(['mrcalc tdfc_average_' PAR_NAME '.mif tdfc_upper_ME_' PAR_NAME '.mif -add tdfc_upper_bound_' PAR_NAME '.mif']);
    
    %do the same steps, but this time for the lower bound (-1.96)
    %multiply -1.96 (95% CI) to the SE image (lower margin of error)
    unix(['mrcalc tdfc_SE_' PAR_NAME '.mif 1.96 -neg -mult tdfc_lower_ME_' PAR_NAME '.mif']); %works
    %add the mean to the lower margin of error to get the lower bound
    unix(['mrcalc tdfc_average_' PAR_NAME '.mif tdfc_lower_ME_' PAR_NAME '.mif -add tdfc_lower_bound_' PAR_NAME '.mif']);
    %subtract upper bound - lower bound to get the CI 95%
    unix(['mrcalc tdfc_upper_bound_' PAR_NAME '.mif tdfc_lower_bound_' PAR_NAME '.mif -subtract tdfc_CI95_' PAR_NAME '.mif']);
    
    

    
    
end

%6. Calculate group mean TW-dFC template from the normalised individual TW-dFC maps
%change zeros into NaNs
unix(['mrgrid tdfc_dynamic_' PAR_NAME '.mif crop -mask upsampled_mask_mod_' PAR_NAME, '.mif -fill nan cropped2_tdfc_dynamic_low-res-2.1_win81_' PAR_NAME '.mif']); % ~ 1.02 GB -- half the size of its uncropped version! 

%4D image
unix(['mrmath tdfc_dynamic*.mif mean average_tdfc_dynamic.mif']);
%3D image
unix(['mrmath tdfc_CI95*.mif mean average_tdfc_CI95.mif']);

unix(['mrcalc tdfc_CI95*.mif -add 2 -div average_tdfc_CI95_mrcalc.mif']); %same as mrmath's mean

unix(['mrmath tdfc_CI95*.mif sum sum_tdfc_CI95.mif']); %add 









