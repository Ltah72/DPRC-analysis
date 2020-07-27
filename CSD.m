%Constrained spherical deconvolution (CSD) pipeline will conduct CSD to
%estimate response functions for fibre orientation distributions (FOD).
%Steps are followed through on the MRtrix manual here:
%https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html


% Steps:
% 1.Generate response function estimates using multi-tissue CSD
%   a.Co-registration t1w and dwi image
%   b.5tt image
%   c.Generate response function (odf.txt file)
% 2.Compute a group average response function of each tissue
% 3.Create FOD images (wm, gm, csf)
%   a.Upsample dwi image
%   b.Upsample brain mask images
%   c.FOD estimation (multi-tissue SD)
% 4.Joint bias field correction and intensity normalistion
% 5.Generate FOD population template
% 6.Register all subjects to FOD template
%   a.Register to FOD template
%   b.Warp mask to template
%   c.Compute intersection of masks in template space

%Author: Lenore Tahara-Eckl
%Date: 30/06/20


clc;
clear all;
close all;

%define/add pathways
startdir = '/data/USERS/LENORE';

%Script directory is defined, so that it can be added to path below:
ScriptDirectory = '/data/USERS/LENORE/scripts/dprc/diffusion';

%go to the directory where your scripts are.
cd(ScriptDirectory);

%should be a groupname from what the user analysed in the preproccessing script.
groupname = input('Which pre-processed group do you want to analyse?: ', 's');

%go to group folder
cd([startdir '/derivatives/diff_data/', groupname]);

%make a folder to hold the participants' data (response functions, fod
%images, brain masks)
mkdir([startdir,'/derivatives/diff_data/', groupname, '/responseFunctionWM/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/responseFunctionGM/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/responseFunctionCSF/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/wmFODimages/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/gmFODimages/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/csfFODimages/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/pop_temp_wmFOD']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/pop_temp_gmFOD/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/pop_temp_csfFOD/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/pop_temp_upsampled_masks/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/upsampled_masks/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/warpedMasks/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/fixel_directory/']);

%Add your script and all necessary files (e.g. data, functions) to path
addpath(genpath(startdir));
addpath(genpath(ScriptDirectory));

%define dwi datafile
datafile = '_acq_data_dwi';

%create 5ttImageCheck text file with header line
fid4 = fopen('5ttImageCheck.txt', 'w');
if (fid4 == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid4, '%s     %s', 'Participant', '5ttImage_status');
    fclose(fid4);
end

%choose which participants will comprise the population template (~40 total):
%message user
disp('Choose participants for population template (~40). You should have an even number of participants to represent each group.');
cd([startdir '/derivatives/diff_data/', groupname, '/preprocessed_dwi/']);
pop_rep_pts = uipickfiles;
cd ../;
fid4 = fopen('pop_temp_pts.txt', 'a+');
if (fid4 == -1)
    disp('Error in opening the text file.')
else 
    for i = 1:length(pop_rep_pts)
        [upper_path, POP_PAR, ~] = fileparts(pop_rep_pts{1,i});
        POP_PAR = POP_PAR(1:15);
        fprintf(fid4, '%s', POP_PAR);
        fprintf(fid4, '\n');
    end
    fclose(fid4);
end
 
%go back into group folder
cd([startdir '/derivatives/diff_data/', groupname]);

%record file numbers and names (participants) in the group folder name.
participants = dir(fullfile('preprocessed_dwi', '*.mif'));


%CSD steps begin below:

%-------------------------------------------------------------------------%
%Step 1: Generate response function estimates using multi-tissue CSD
for i = 1:length(participants)
    
    full_name = participants(i).name;
    PAR_NAME = full_name(1:15);
    
    %copy participant's raw t1w and t2 FLAIR images, and the preprocessed/segmented 
    %t2 FLAIR images into main group folder
    copyfile ([startdir '/derivatives/diff_data/' PAR_NAME '/anat/', '*.nii'], [startdir,'/derivatives/diff_data/', groupname]);
    
    %a) co-registration of t1w and t2 FLAIR to dwi image through fsl FLIRT
    %convert to nii format to run fsl's FLIRT - use the best B0 volume,
    %which will be the first vol in the dwi sequence (vol 0).
    unix(['mrconvert -coord 3 0 preprocessed_dwi/', PAR_NAME, datafile,'.mif ref_b0_', PAR_NAME, '.nii']);
    
    %linear registration with 6 dof and the tranformation matrix output
    unix(['flirt -in ', PAR_NAME, '_T1w.nii -ref ref_b0_', PAR_NAME, '.nii -dof 6 -out t1_flirt_' PAR_NAME, '.nii -omat transform_flirt_t12dwi.mat']);
    unix(['flirt -in ', PAR_NAME, '_FLAIR.nii -ref ref_b0_', PAR_NAME, '.nii -dof 6 -out t2_flirt_' PAR_NAME, '.nii -omat transform_flirt_t22dwi.mat']);

    %apply the linear transformation to the t1 and t2 image
    unix(['transformconvert transform_flirt_t12dwi.mat ' PAR_NAME, '_T1w.nii ref_b0_', PAR_NAME, '.nii flirt_import transform_mrtrix_t12dwi.mat']);
    unix(['transformconvert transform_flirt_t22dwi.mat ' PAR_NAME, '_FLAIR.nii ref_b0_', PAR_NAME, '.nii flirt_import transform_mrtrix_t22dwi.mat']);
    unix(['mrtransform ' PAR_NAME, '_T1w.nii -linear transform_mrtrix_t12dwi.mat ' PAR_NAME, '_T1coereg.nii']); 
    unix(['mrtransform ' PAR_NAME, '_FLAIR.nii -linear transform_mrtrix_t12dwi.mat -template ' PAR_NAME, '_T1coereg.nii ' PAR_NAME, '_T2coereg.nii']); 

    %convert to .mif (mrtrix) format
    unix(['mrconvert ' PAR_NAME, '_T1coereg.nii ' PAR_NAME '_T1coereg.mif']);
    unix(['mrconvert ' PAR_NAME, '_T2coereg.nii ' PAR_NAME '_T2coereg.mif']);
    
    %b) generate 5tt image with brain mask, t1, and t2 FLAIR
    unix(['5ttgen fsl -mask brain_mask/', PAR_NAME, datafile, '.mif -t2 ' PAR_NAME, '_T2coereg.mif ' PAR_NAME '_T1coereg.mif 4ttimage_' PAR_NAME '.mif']);
    
    %c) edit in the pathological tissue (WMH) to the 5tt image:
    %convert WMH lesion mask to .mif file
    unix(['mrconvert ples_lpa_mr', PAR_NAME, '_FLAIR.nii.gz ples_lpa_mr', PAR_NAME, '_FLAIR.mif -force']);
    
    %reslice the WMH lesion mask to match the 4tt image
    unix(['mrtransform ples_lpa_mr' PAR_NAME, '_FLAIR.mif -template 4ttimage_' PAR_NAME, '.mif ' PAR_NAME, '_WMH_mask_transformed.mif']); 

    %add the image into the 5tt image
    unix(['5ttedit -path ', PAR_NAME, '_WMH_mask_transformed.mif 4ttimage_' PAR_NAME, '.mif 5ttimage_' PAR_NAME '.mif']);
    
    %check that your 5ttimage conforms to MRtrix's expected format
    FiveTTImageCheck(PAR_NAME);
    
    %d) generate a response function per each subject - use multi-shell data option
    unix(['dwi2response msmt_5tt -mask brain_mask/', PAR_NAME, datafile, '.mif preprocessed_dwi/', PAR_NAME, datafile, '.mif 5ttimage_', PAR_NAME, '.mif ' PAR_NAME '_out_wm ' PAR_NAME '_out_gm ' PAR_NAME '_out_csf']);
    
    %place a copy of each response function (e.g. out_wm.txt) into a separate folder
    copyfile ([PAR_NAME '_out_wm'], [startdir,'/derivatives/diff_data/', groupname, '/responseFunctionWM/']);
    copyfile ([PAR_NAME '_out_gm'], [startdir,'/derivatives/diff_data/', groupname, '/responseFunctionGM/']);
    copyfile ([PAR_NAME '_out_csf'], [startdir,'/derivatives/diff_data/', groupname, '/responseFunctionCSF/']);
    
end

%-------------------------------------------------------------------------%
%Step 2: Compute a group average response function of each tissue
unix(['responsemean ' startdir, '/derivatives/diff_data/' groupname, '/responseFunctionWM/* group_average_responseWM.txt']);
unix(['responsemean ' startdir, '/derivatives/diff_data/' groupname, '/responseFunctionGM/* group_average_responseGM.txt']);
unix(['responsemean ' startdir, '/derivatives/diff_data/' groupname, '/responseFunctionCSF/* group_average_responseCSF.txt']);

%-------------------------------------------------------------------------%
%Step 3: Create FOD images (wm, gm, csf fod images)

for j = 1:length(participants)
    
    full_name = participants(j).name;
    PAR_NAME = full_name(1:15);
    
    %a) upsample DW images
    %regrid option (Gaussian smoothing)
    unix(['mrgrid -voxel 1.25 preprocessed_dwi/' PAR_NAME, datafile, '.mif regrid upsampled_dwi_' PAR_NAME '.mif']);

    %b) compute upsampled brain mask images (use bet)
    unix(['bet ref_b0_' PAR_NAME, '.nii upsampled_mask_' PAR_NAME, '.nii -m -f 0.2']);
    unix(['mrconvert upsampled_mask_' PAR_NAME, '_mask.nii.gz upsampled_mask_' PAR_NAME '.mif']);
    
    %copy and rename each of the upsampled brain mask into a folder
    copyfile (['upsampled_mask_' PAR_NAME '.mif'], [startdir,'/derivatives/diff_data/', groupname, '/upsampled_masks/']);
    movefile(['upsampled_masks/upsampled_mask_' PAR_NAME '.mif'], ['upsampled_masks/PRE_' PAR_NAME '.mif']);
   
    %c) fibre orientation distribution estimation (using group response
    %estimates for each tissue and the upsampled images)
    unix(['dwi2fod msmt_csd upsampled_dwi_' PAR_NAME '.mif group_average_responseWM.txt wmfod_' PAR_NAME '.mif group_average_responseGM.txt gmfod_' PAR_NAME '.mif group_average_responseCSF.txt csffod_' PAR_NAME '.mif -mask upsampled_mask_' PAR_NAME '.mif']);
    
   
    %---------------------------------------------------------------------%
    %Step 4: Joint bias field correction and intensity normalisation
    unix(['mtnormalise wmfod_' PAR_NAME '.mif wmfod_norm_' PAR_NAME '.mif gmfod_' PAR_NAME '.mif gmfod_norm_' PAR_NAME '.mif csffod_' PAR_NAME '.mif csffod_norm_' PAR_NAME '.mif -mask upsampled_mask_' PAR_NAME '.mif']);
    
    
    %copy and rename each of the FOD images (e.g. wmfod_PARNAME.mif) into a separate folder
    copyfile (['wmfod_norm_' PAR_NAME '.mif'], [startdir,'/derivatives/diff_data/', groupname, '/wmFODimages/']);
    movefile(['wmFODimages/wmfod_norm_' PAR_NAME '.mif'], ['wmFODimages/PRE_' PAR_NAME '.mif']);
    copyfile (['gmfod_norm_' PAR_NAME '.mif'], [startdir,'/derivatives/diff_data/', groupname, '/gmFODimages/']);
    movefile(['gmFODimages/gmfod_norm_' PAR_NAME '.mif'], ['gmFODimages/PRE_' PAR_NAME '.mif']);
    copyfile (['csffod_norm_' PAR_NAME '.mif'], [startdir,'/derivatives/diff_data/', groupname, '/csfFODimages/']);
    movefile(['csfFODimages/csffod_norm_' PAR_NAME '.mif'], ['csfFODimages/PRE_' PAR_NAME '.mif']);
    
end


%-------------------------------------------------------------------------%
%Step 5: Generate a study-specific unbiased FOD template - this step is 
%computationally expensive and takes a long time to process. Will pick out 
%the participants that you chose at the beginning. 
for i = 1:length(pop_rep_pts)
    [upper_path, POP_PAR, ~] = fileparts(pop_rep_pts{1,i});
    POP_PAR = POP_PAR(1:15);
    copyfile (['wmfod_norm_' POP_PAR '.mif'], [startdir,'/derivatives/diff_data/', groupname, '/pop_temp_wmFOD/']);
    copyfile (['gmfod_norm_' POP_PAR '.mif'], [startdir,'/derivatives/diff_data/', groupname, '/pop_temp_gmFOD/']);
    copyfile (['csffod_norm_' POP_PAR '.mif'], [startdir,'/derivatives/diff_data/', groupname, '/pop_temp_csfFOD/']);
    copyfile (['upsampled_mask_' POP_PAR, '.mif'], [startdir,'/derivatives/diff_data/', groupname, '/pop_temp_upsampled_masks/']);
end

unix(['population_template pop_temp_wmFOD wmfod_template.mif pop_temp_gmFOD gmfod_template.mif pop_temp_csfFOD csffod_template.mif -mask_dir pop_temp_upsampled_masks -force']);

%-------------------------------------------------------------------------%
%Step 6: Register all subjects to FOD template

for k = 1:length(participants)
    
    full_name = participants(k).name;
    PAR_NAME = full_name(1:15);
    
    %register all subject FOD images to the FOD template
    unix(['mrregister wmfod_' PAR_NAME '.mif -mask1 upsampled_mask_' PAR_NAME '.mif wmfod_template.mif -nl_warp ' PAR_NAME '_subject2template_warp.mif ' PAR_NAME '_template2subject_warp.mif']);
   
    %warp all subject masks into template space
    unix(['mrtransform upsampled_mask_' PAR_NAME '.mif -warp ' PAR_NAME '_subject2template_warp.mif -interp nearest ' PAR_NAME '_output_warped_mask.mif']);
    
    %place a copy of output warped mask (e.g. PAR_NAME_output_warped_mask.mif) into a folder
    copyfile ([PAR_NAME '_output_warped_mask.mif'], [startdir,'/derivatives/diff_data/', groupname, '/warpedMasks/']);
    
end

%Compute the intersection of all subject masks in template space
unix(['mrmath ' startdir, '/derivatives/diff_data/' groupname, '/warpedMasks/* min fixel_directory/template_mask_intersection.mif -datatype bit']);



