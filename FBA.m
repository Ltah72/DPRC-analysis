%Fixel-based analysis (FBA) pipeline, which is sequential to running the 
%CSD pipeline. The FBA pipeline will compute the diffusion fibre metrics,
%which are fibre density (FD), fibre cross-section (FC), and fibre density
%and cross-section (FDC). Statistical analysis will be performed on these
%metrics, and whole-brain fibre tractography will be conducted. 

%MRtrix manual: https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html
	
% Steps:
% 1. Compute a white matter template analysis fixel mask 
% 2. Estimate participants’ fixels and diffusion fibre metric 
%   a.Warp FOD images to template space
%   b.Segment FOD images to estimate fixels and their FD metric
%   c.Reorient fixels
%   d.Assign subject fixels to template fixels (compute FD)
%   e.Compute FC metric
%   f.Compute FDC metric
% 3. Perform whole-brain fibre tractography on the FOD template
% 4. Reduce biases in tractogram densities (using SIFT)
% 5. Generate fixel-fixel connectivity matrix
% 6. Smooth fixel data using fixel-fixel connectivity
% 7. Perform statistical analysis of FD, FC, and FDC
% 8. Visualise results (with mrview)
% 9. Display results with streamlines
%   a.Reduce number of streamlines to 200k
%   b.Create .tsf file (map fixels to streamlines)
%   c.Visualise .tsf in mrview
% 10. FBA post-statistical inference
%   a.Calculate whole-brain FBA metrics per each participant and put onto a
%   text file
%   b.Express the effect size relative to controls

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 30/06/20

clc;
clear all;
close all;

%define/add pathways
startdir = '/data/USERS/LENORE';

%Script directory is defined, so that it can be added to path below:
ScriptDirectory = '/data/USERS/LENORE/scripts/dprc/diffusion';

%should be the same groupname from what the user analysed in the CSD script.  
groupname = input('Which pre-processed group / study do you want to continue to analyse?: ', 's');

%make directories for AFD data
mkdir([startdir,'/derivatives/diff_data/', groupname, '/log_fc_data/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/fd_data/']);
mkdir([startdir,'/derivatives/diff_data/', groupname, '/fdc_data/']);

%Add your script and all necessary files (e.g. data, functions) to path
addpath(genpath(startdir));
addpath(genpath(ScriptDirectory));

%go into group folder
cd([startdir '/derivatives/diff_data/', groupname]);

%record file numbers and names (participants) in the group folder name. 
participants = dir(fullfile('preprocessed_dwi', '*.mif'));

%FBA steps begin below: 
%-------------------------------------------------------------------------%
%Step 1: Compute a white matter template analysis fixel mask
unix(['fod2fixel -mask fixel_directory/template_mask_intersection.mif -fmls_peak_value 0.06 wmfod_template.mif fixel_directory/fixel_mask']);


%-------------------------------------------------------------------------%
%Step 2: Estimate participants’ fixels and FBA metrics 
for i = 1:length(participants)
    
    full_name = participants(i).name;
    PAR_NAME = full_name(1:15);
    
    %a) Warp each participant's FOD images to template space.    
    unix(['mrtransform wmfod_norm_' PAR_NAME '.mif -warp ' PAR_NAME '_subject2template_warp.mif -reorient_fod no ' PAR_NAME '_fod_in_template_space_NOT_REORIENTED.mif']);
    
    %b) segment FOD images to estimate fixels and their fibre density (FD)    
    mkdir(['fixel_directory/' PAR_NAME]);
    %unix(['fod2fixel -mask analysis_voxel_mask.mif ' PAR_NAME '_fod_in_template_space_NOT_REORIENTED.mif fixel_directory/' PAR_NAME '/fod_in_template_space_NOT_REORIENTED -afd fd.mif']);
    unix(['fod2fixel -mask fixel_directory/template_mask_intersection.mif ' PAR_NAME '_fod_in_template_space_NOT_REORIENTED.mif fixel_directory/' PAR_NAME '/fixel_in_template_space_NOT_REORIENTED -afd fd.mif']);
    
    %c) reorient fixel orientations
    unix(['fixelreorient fixel_directory/' PAR_NAME '/fixel_in_template_space_NOT_REORIENTED ' PAR_NAME '_subject2template_warp.mif fixel_directory/' PAR_NAME '/fixel_in_template_space']);
    
    %d) assign subject fixels to template fixels (compute FD)
    unix(['fixelcorrespondence fixel_directory/' PAR_NAME '/fixel_in_template_space/fd.mif fixel_directory/fixel_mask fixel_directory/' PAR_NAME '/fd fd.mif']);
    
    %e) compute fibre cross-section (FC) metric
    unix(['warp2metric ' PAR_NAME '_subject2template_warp.mif -fc fixel_directory/fixel_mask fixel_directory/' PAR_NAME '/fc fc.mif']);
       
    %for group analysis FC - take the log of FC
    mkdir(['fixel_directory/' PAR_NAME '/log_fc']);
    copyfile (['fixel_directory/' PAR_NAME '/fc/index.mif'], ['fixel_directory/' PAR_NAME '/log_fc']);
    copyfile (['fixel_directory/' PAR_NAME '/fc/directions.mif'], ['fixel_directory/' PAR_NAME '/log_fc']);
    unix(['mrcalc fixel_directory/' PAR_NAME '/fc/fc.mif -log fixel_directory/' PAR_NAME '/log_fc/log_fc.mif']);
    
    %f) Compute a combined measure of fibre density and cross-section (FDC)
    mkdir(['fixel_directory/' PAR_NAME '/fdc']);
    copyfile (['fixel_directory/' PAR_NAME '/fc/index.mif'], ['fixel_directory/' PAR_NAME '/fdc']);
    copyfile (['fixel_directory/' PAR_NAME '/fc/directions.mif'], ['fixel_directory/' PAR_NAME '/fdc']);
    unix(['mrcalc fixel_directory/' PAR_NAME '/fd/fd.mif fixel_directory/' PAR_NAME '/fc/fc.mif -mult fixel_directory/' PAR_NAME '/fdc/fdc.mif']);
    
    %copy and rename each participant file into the 3 metric files (fd, fc, fdc) into separate folders
    copyfile (['fixel_directory/' PAR_NAME '/fd/fd.mif'], ['fd_data/' PAR_NAME '_fd.mif']);
    copyfile (['fixel_directory/' PAR_NAME '/log_fc/log_fc.mif'], ['log_fc_data/' PAR_NAME '_log_fc.mif']);
    copyfile (['fixel_directory/' PAR_NAME '/fdc/fdc.mif'], ['fdc_data/' PAR_NAME '_fdc.mif']);

end 

%copy 1 index and 1 directions file of each metric (fd, fc, and fdc) into 
%the 3 metric files (any 1 participant file is okay to use for this).
copyfile (['fixel_directory/' PAR_NAME '/fd/index.mif'], ['fd_data/']);
copyfile (['fixel_directory/' PAR_NAME '/fd/directions.mif'], ['fd_data/']);
copyfile (['fixel_directory/' PAR_NAME '/log_fc/index.mif'], ['log_fc_data/']);
copyfile (['fixel_directory/' PAR_NAME '/log_fc/directions.mif'], ['log_fc_data/']);
copyfile (['fixel_directory/' PAR_NAME '/fdc/index.mif'], ['fdc_data/']);
copyfile (['fixel_directory/' PAR_NAME '/fdc/directions.mif'], ['fdc_data/']);

%-------------------------------------------------------------------------%
%Step 3: Perform whole-brain fibre tractography on the FOD template -- this 
%step is computationally expensive and takes a long time to process
unix(['tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 wmfod_template.mif -seed_image fixel_directory/template_mask_intersection.mif -mask fixel_directory/template_mask_intersection.mif -select 20000000 -cutoff 0.06 tracks_20_million.tck']);

%-------------------------------------------------------------------------%
%Step 4:Reduce biases in tractogram densities (using SIFT)
unix(['tcksift tracks_20_million.tck wmfod_template.mif tracks_2_million_sift.tck -term_number 2000000']);
    
%-------------------------------------------------------------------------%
%Step 5: Generate fixel-fixel connectivity matrix
unix(['fixelconnectivity fixel_directory/fixel_mask/ tracks_2_million_sift.tck matrix/']);

%-------------------------------------------------------------------------%
%Step 6: Smooth fixel data using fixel-fixel connectivity 
unix(['fixelfilter fd_data smooth fd_smooth -matrix matrix/']);
unix(['fixelfilter log_fc_data smooth log_fc_smooth -matrix matrix/']);
unix(['fixelfilter fdc_data smooth fdc_smooth -matrix matrix/']);


%-------------------------------------------------------------------------%
%Step 7: Perform statistical analysis of FD, FC, and FDC -- this step is 
%computationally expensive and takes a long time to process

%create text files of a participant list to run statistics 
CreateParticipantFixelList;

%use the design matrix (e.g. design_matrix.txt) that corresponds to the participant list (participant
%groups - i.e. status)

%create contrast matrix to specify the tests and/or co-variates that you
%will include in your analysis (e.g. contrast_matrix.txt)

%Run statistical tests for each AFD metric
unix(['fixelcfestats fd_smooth/ files_fd.txt stats_matrices/design_matrix.txt stats_matrices/contrast_matrix.txt matrix/ stats_fd/']);
unix(['fixelcfestats log_fc_smooth/ files_log_fc.txt stats_matrices/design_matrix.txt stats_matrices/contrast_matrix.txt matrix/ stats_log_fc/']);
unix(['fixelcfestats fdc_smooth/ files_fdc.txt stats_matrices/design_matrix.txt stats_matrices/contrast_matrix.txt matrix/ stats_fdc/']);



%-------------------------------------------------------------------------%
%Step 8: Visualise results
%To view the results load the population FOD template image in mrview, and 
%overlay the fixel images using the vector plot tool. Note that p-value 
%images are saved as (1 - p-value). Therefore to visualise all results at a 
%threshold of p < 0.05, within the mrview fixel plot tool, apply a lower 
%threshold at a value of 0.95.



%-------------------------------------------------------------------------%
%Step 9: Display results with streamlines
%cropping streamlines from the template-derived whole-brain tractogram to include streamline points that correspond to significant fixels

%Reduce number of streamlines to 200,000
unix(['tckedit tracks_2_million_sift.tck -num 200000 tracks_200k_sift.tck']);

%map fixel values to streamline points and save them in a 'track scalar file' (.tsf)
unix(['fixel2tsf stats_fd/fwe_1mpvalue.mif tracks_200k_sift.tck fd_WholeBrainfwe_pvalue.tsf']);
unix(['fixel2tsf stats_log_fc/fwe_1mpvalue.mif tracks_200k_sift.tck log_fc_WholeBrainfwe_pvalue.tsf']);
unix(['fixel2tsf stats_fdc/fwe_1mpvalue.mif tracks_200k_sift.tck fdc_WholeBrainfwe_pvalue.tsf']);

%visualise track scalar files using the tractogram tool in MRview. First 
%load the streamlines (tracks_200k_sift.tck). Then to dynamically threshold 
%(remove) streamline point by p-value select the “Thresholds” dropdown and 
%select “Separate Scalar file” and set to 0.95.


%-------------------------------------------------------------------------%
%Step 10: FBA post-statistical inference

%a) Calculate whole-brain FBA metrics per each participant and put onto a text file. 
CreateWholeBrainFBAMetricFiles(participants);
%use this files as inputs into a statistical programme in order to 
%calculate the statistics of the group differences


%b) Expressing the effect size relative to controls
%for fd
unix(['mrcalc stats_fd/abs_effect.mif stats_fd/beta0.mif  -div 100 -mult stats_fd/percentage_effect.mif -force']);
%for fdc
unix(['mrcalc stats_fdc/abs_effect.mif stats_fdc/beta0.mif  -div 100 -mult stats_fdc/percentage_effect.mif']);
%for fc
unix(['mrcalc 1 1 stats_log_fc/abs_effect.mif -exp -div -sub stats_log_fc/percentage_effect.mif']);

    
    
    
    
