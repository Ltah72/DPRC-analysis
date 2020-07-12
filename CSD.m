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
ScriptDirectory = '/data/USERS/LENORE/scripts/dprc/DKI';

%should be a groupname from what the user analysed in the preproccessing script.  
groupname = input('Which pre-processed group do you want to analyse?: ', 's');

%make a folder to hold the participants' data (response functions, fod
%images, brain masks)
mkdir([startdir,'/test/derivatives/', groupname, '/responseFunctionWM/']);
mkdir([startdir,'/test/derivatives/', groupname, '/responseFunctionGM/']);
mkdir([startdir,'/test/derivatives/', groupname, '/responseFunctionCSF/']);
mkdir([startdir,'/test/derivatives/', groupname, '/wmFODimages/']);
mkdir([startdir,'/test/derivatives/', groupname, '/gmFODimages/']);
mkdir([startdir,'/test/derivatives/', groupname, '/csfFODimages/']);
mkdir([startdir,'/test/derivatives/', groupname, '/upsampled_masks/']);
mkdir([startdir,'/test/derivatives/', groupname, '/warpedMasks/']);
mkdir([startdir,'/test/derivatives/', groupname, '/fixel_directory/']);

%Add your script and all necessary files (e.g. data, functions) to path
addpath(genpath(startdir));
addpath(genpath(ScriptDirectory));

%define dwi datafile
datafile = '_acq_data_dwi';

%go into group folder
cd([startdir '/test/derivatives/', groupname]);

%record file numbers and names (participants) in the group folder name. 
participants = dir(fullfile('preprocessed_dwi', '*.mif'));


%CSD steps begin below: 

%-------------------------------------------------------------------------%
%Step 1: Generate response function estimates using multi-tissue CSD
for i = 1:length(participants)
    
    full_name = participants(i).name;
    PAR_NAME = full_name(1:15);
    
    %copy participant's t1w images into main group folder
    copyfile ([startdir '/test/derivatives/' PAR_NAME '/anat/', '*.nii'], [startdir,'/test/derivatives/', groupname]);
    
    t1w = ([PAR_NAME(1:13), '_T1w']);
    
    %copy participant's brain mask into main group folder (rename as w/
    %start file name brain_mask_). 
    %copyfile ([startdir '/test/derivatives/' groupname '/input_brain_mask/', PAR_NAME, datafile, '.mif'], [startdir,'/test/derivatives/', groupname '/', 'brain_mask_' PAR_NAME, datafile, '.mif'])
 
    %a) co-registration of t1w and dwi image through fsl FLIRT
    %convert to nii format to run fsl's FLIRT - use the best B0 volume, 
    %which will be the first vol in the dwi sequence (vol 0).
    unix(['mrconvert -coord 3 0 preprocessed_dwi/', PAR_NAME, datafile,'.mif AP_cleaned_', PAR_NAME, '.nii ']);
    
    unix(['flirt -in ', t1w, '.nii -ref ref_b0_', PAR_NAME, '.nii -dof 12 -out anat_flirt_' PAR_NAME, '.nii']);
    
    %convert to .mif (mrtrix) format
    unix(['mrconvert anat_flirt_' PAR_NAME, '.nii.gz anat_flirt_' PAR_NAME '.mif']);
    
    %b) generate 5tt image
    unix(['5ttgen fsl anat_flirt_' PAR_NAME '.mif 5ttimage_' PAR_NAME '.mif']);
    
    %c)generate a response function per each subject - use multi-shell data option
    unix(['dwi2response msmt_5tt preprocessed_dwi/', PAR_NAME, datafile, '.mif 5ttimage_', PAR_NAME, '.mif ' PAR_NAME '_out_wm ' PAR_NAME '_out_gm ' PAR_NAME '_out_csf']);
    
    %place a copy of each response function (e.g. out_wm.txt) into a separate folder
    copyfile ([PAR_NAME '_out_wm'], [startdir,'/test/derivatives/', groupname, '/responseFunctionWM/']);
    copyfile ([PAR_NAME '_out_gm'], [startdir,'/test/derivatives/', groupname, '/responseFunctionGM/']);
    copyfile ([PAR_NAME '_out_csf'], [startdir,'/test/derivatives/', groupname, '/responseFunctionCSF/']);
    
end

%-------------------------------------------------------------------------%
%Step 2: Compute a group average response function of each tissue
unix(['responsemean ' startdir, '/test/derivatives/' groupname, '/responseFunctionWM/* group_average_responseWM.txt']);  
unix(['responsemean ' startdir, '/test/derivatives/' groupname, '/responseFunctionGM/* group_average_responseGM.txt']);  
unix(['responsemean ' startdir, '/test/derivatives/' groupname, '/responseFunctionCSF/* group_average_responseCSF.txt']);  

%-------------------------------------------------------------------------%
%Step 3: Create FOD images (wm, gm, csf fod images)

for j = 1:length(participants)
 
    full_name = participants(j).name;
    PAR_NAME = full_name(1:15);
    
    %a) upsample DW images
    %regrid option (Gaussian smoothing)
    unix(['mrgrid -voxel 1.25 preprocessed_dwi/' PAR_NAME, datafile, '.mif regrid upsampled_dwi3_' PAR_NAME '.mif']);
   
    %b) compute upsampled brain mask images
    unix(['dwi2mask upsampled_dwi_' PAR_NAME '.mif upsampled_mask_' PAR_NAME, '.mif']);
    
    %copy and rename each of the upsampled brain mask into a folder
    copyfile (['upsampled_mask_' PAR_NAME '.mif'], [startdir,'/test/derivatives/', groupname, '/upsampled_masks/']);
    movefile(['upsampled_masks/upsampled_mask_' PAR_NAME '.mif'], ['upsampled_masks/PRE_' PAR_NAME '.mif']);
    
    %c) fibre orientation distribution estimation (using group response
    %estimates for each tissue and the upsampled images)
    unix(['dwi2fod msmt_csd upsampled_dwi_' PAR_NAME '.mif group_average_responseWM.txt wmfod_' PAR_NAME '.mif group_average_responseGM.txt gmfod_' PAR_NAME '.mif group_average_responseCSF.txt csffod_' PAR_NAME '.mif -mask upsampled_mask_' PAR_NAME '.mif']);
    
    
    %---------------------------------------------------------------------%
    %Step 4: Joint bias field correction and intensity normalisation
    unix(['mtnormalise wmfod_' PAR_NAME '.mif wmfod_norm_' PAR_NAME '.mif gmfod_' PAR_NAME '.mif gmfod_norm_' PAR_NAME '.mif csffod_' PAR_NAME '.mif csffod_norm_' PAR_NAME '.mif -mask upsampled_mask_' PAR_NAME '.mif']);
        
    
    %copy and rename each of the FOD images (e.g. wmfod_PARNAME.mif) into a separate folder
    copyfile (['wmfod_norm_' PAR_NAME '.mif'], [startdir,'/test/derivatives/', groupname, '/wmFODimages/']);
    movefile(['wmFODimages/wmfod_norm_' PAR_NAME '.mif'], ['wmFODimages/PRE_' PAR_NAME '.mif']);
    copyfile (['gmfod_norm_' PAR_NAME '.mif'], [startdir,'/test/derivatives/', groupname, '/gmFODimages/']);
    movefile(['gmFODimages/gmfod_norm_' PAR_NAME '.mif'], ['gmFODimages/PRE_' PAR_NAME '.mif']);
    copyfile (['csffod_norm_' PAR_NAME '.mif'], [startdir,'/test/derivatives/', groupname, '/csfFODimages/']);
    movefile(['csfFODimages/csffod_norm_' PAR_NAME '.mif'], ['csfFODimages/PRE_' PAR_NAME '.mif']);
    
end


%-------------------------------------------------------------------------%
%Step 5: Generate a study-specific unbiased FOD template - this step takes
%a long time to process
unix(['population_template wmFODimages wmfod_template.mif gmFODimages gmfod_template.mif csfFODimages csffod_template.mif -mask_dir upsampled_masks']);


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
    copyfile ([PAR_NAME '_output_warped_mask.mif'], [startdir,'/test/derivatives/', groupname, '/warpedMasks/']);
    
end 

%Compute the intersection of all subject masks in template space
unix(['mrmath ' startdir, '/test/derivatives/' groupname, '/warpedMasks/* min fixel_directory/template_mask_intersection.mif -datatype bit']);



