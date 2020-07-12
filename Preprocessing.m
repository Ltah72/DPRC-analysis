%Preprocessing script for DPRC diffusion data. Preprocessing is based upon
%MRtrix recommended pre-processing steps:
%https://mrtrix.readthedocs.io/en/0.3.16/workflows/DWI_preprocessing_for_quantitative_analysis.html
%And UKB data recommendation:
%Maximov & Westlye (2019).Towards an optimised processing pipeline for diffusion magnetic resonance imaging data: Effects of artefact corrections on diffusion metrics and their age associations in UK Biobank.

%All preprocessing steps are deployed through the MRtrix3 programme.
%Utilised neuroimaging programs are FSL, MRtrix3, and ANTs are for
%these steps. Assumes BIDS formatting and organisation. This script calls
%upon 6 preprocessing steps (as functions) which are, in order:

%Steps                                          Methodology - programme origin, algorithm/model, author
%1. Noise correction                            (denoising -- MP-PCA, Veraart et al., 2016)
%2. Gibbs ringing correction                    (local sub-voxel shift, Kellner et al., 2016)
%3. Field distortion                            (TOPUP -- FSL)
%4. Eddy current distortions                    (Eddy -- FSL)
%5  Estimate a brain mask
%6. Bias field correction                       (ANTs -- N4BiasFieldCorrection)


%Author: Lenore Tahara-Eckl
%Date: 30/06/20

clc;
clear all;
close all;


%Define script directory, so that it can be added to path below:
%ScriptDirectory = input('Please enter script directory:', 's');
ScriptDirectory = '/data/USERS/LENORE/scripts/dprc/DKI';

%go to the directory where your scripts are.
cd(ScriptDirectory);

%Choose directory where the diffusion data is (same level as sourcedata
%folder).
%startdir = input('Please enter data directory:', 's');

%shortcut for debugging purposes:
startdir = '/data/USERS/LENORE';


%ask user for what group they want to analyse:
groupname = input('Please name a group that you want to analyse: ', 's');

%make directories
mkdir([startdir,'/test/derivatives/', groupname, '/preprocessed_dwi/']);
mkdir([startdir,'/test/derivatives/', groupname, '/brain_mask/']);

%Add your script and all necessary files (e.g. data, functions) to path
addpath(genpath(startdir));
addpath(genpath(ScriptDirectory));

%define dwi datafile
datafile = '_acq_data_dwi';


%go into your source data folder, where all participant files will be
%there.
cd([startdir '/test/']);

%create BestB0 text file with header line
fid3 = fopen('BestB0.txt', 'w');
if (fid3 == -1)
    disp('Error in opening in one or both of the correlation score files.')
else
    fprintf(fid3, '%s     %s %s %s', 'Participant', 'B0_status', 'BU_used', 'BD_used');
    fclose(fid3);
end
%choose participants you want to include in the group.
subjects = uipickfiles;

%Preprocessing starts here:
for i = 1:length(subjects)
    
    [upper_path, PAR_NAME, ~] = fileparts(subjects{1,i});
    full_path = ([upper_path '/' PAR_NAME '/']);
    
    %PAR_NAME = deepest_dir;
    
    %source = ([startdir '/sourcedata/', PAR_NAME, '/dwi/']);
    source_dwi = ([startdir '/test/', PAR_NAME, '/dwi/']);
    source_anat = ([startdir '/test/', PAR_NAME, '/anat/']);
    
    %create derivatives folder per each participant
    %mkdir([startdir,'/derivatives/', PAR_NAME, '/dwi/']);
    mkdir([startdir,'/test/derivatives/', PAR_NAME, '/dwi/']);
    mkdir([startdir,'/test/derivatives/', PAR_NAME, '/anat/']);
    
    %copy data from source folder into participants derivatives folder
    copyfile ([source_dwi, '*'], [startdir,'/test/derivatives/', PAR_NAME, '/dwi/']);
    copyfile ([source_anat, '*'], [startdir,'/test/derivatives/', PAR_NAME, '/anat/']);
    
    %move into participants folder
    %cd([startdir '/derivatives/' PAR_NAME, '/dwi/']);
    cd([startdir '/test/derivatives/' PAR_NAME, '/dwi/']);
    
    %count number of BDs that the participant has.
    BDs = dir([PAR_NAME '*BD*.nii']);
    
    %combine all volumes together (main, BUs, BDs) for denoising and Gibbs
    %ringing correction.
    if length(BDs) == 1
        NumBDs = 1;
        unix(['mrcat ' PAR_NAME, datafile, '.nii ', PAR_NAME, '_acq_BU_dwi.nii ', PAR_NAME, '_acq_BD1_dwi.nii ', 'combined_', PAR_NAME, datafile, '.nii']); %106 volumes total
    elseif length(BDs) == 2
        NumBDs = 2;
        unix(['mrcat ' PAR_NAME, datafile, '.nii ', PAR_NAME, '_acq_BU_dwi.nii ', PAR_NAME, '_acq_BD1_dwi.nii ', PAR_NAME, '_acq_BD2_dwi.nii ' 'combined_', PAR_NAME, datafile, '.nii']); %107 volumes total
    elseif length(BDs) == 3
        NumBDs = 3;
        unix(['mrcat ' PAR_NAME, datafile, '.nii ', PAR_NAME, '_acq_BU_dwi.nii ', PAR_NAME, '_acq_BD1_dwi.nii ', PAR_NAME, '_acq_BD2_dwi.nii ', PAR_NAME, '_acq_BD3_dwi.nii ', 'combined_', PAR_NAME, datafile, '.nii']); %108 volumes total
    end
    
    %convert nifti to to .mif format for mrtrix3 processing
    unix(['mrconvert combined_', PAR_NAME, datafile,'.nii combined_', PAR_NAME, datafile,'.mif']);
    
    %---------------------------------------------------------------------%
    %Step 1: Denoising
    unix(['dwidenoise combined_', PAR_NAME, datafile,'.mif d', PAR_NAME, datafile,'.mif -noise noise_' PAR_NAME, datafile,'.mif']);
    
    %create a copy in NIFIT format, in case you want to examine in fsleyes
    unix(['mrconvert d', PAR_NAME, datafile,'.mif d', PAR_NAME, datafile, '.nii']);
    unix(['mrconvert noise_', PAR_NAME, datafile,'.mif noise_', PAR_NAME, datafile, '.nii']);
    
    %eyeball your residuals as part of quality control; compare raw data
    %with the denoised data. If denoising did a good job, there should be
    %little or no anatomy in the residual maps or the diffusion weighted
    %images (this does not apply to the intial B0).
    unix(['mrcalc combined_', PAR_NAME, datafile,'.mif d', PAR_NAME, datafile '.mif -subtract res_' PAR_NAME, datafile,'.mif']);
    
    %create a copy in NIFIT format
    unix(['mrconvert res_', PAR_NAME, datafile,'.mif res_', PAR_NAME, datafile, '.nii']);
    
    %---------------------------------------------------------------------%
    %Step 2: Gibbs ringing
    unix(['mrdegibbs d' PAR_NAME, datafile,'.mif gd' PAR_NAME, datafile,'.mif']);
    
    %create a copy in NIFIT format
    unix(['mrconvert gd', PAR_NAME, datafile,'.mif gd', PAR_NAME, datafile, '.nii']);
    
    %---------------------------------------------------------------------%
    %intermediate steps: edit gradient files & organising
    
    %edit gradient text files (bval & bvec files). We need to add on the
    %last B0 file (vol 106) to the dataset.
    LastB0AddOn(PAR_NAME, datafile);
    
    %separate data into main dwi data and p-a files.
    unix(['mrconvert -coord 3 0:105 gd' PAR_NAME, datafile,'.mif cgd' PAR_NAME, datafile,'.mif']);
    
    if length(BDs) == 1
        unix(['mrconvert -coord 3 106 gd' PAR_NAME, datafile,'.mif PA_' PAR_NAME, datafile,'.mif']);
    elseif length(BDs) == 2
        unix(['mrconvert -coord 3 106:107 gd' PAR_NAME, datafile,'.mif PA_' PAR_NAME, datafile,'.mif']);
    elseif length(BDs) == 3
        unix(['mrconvert -coord 3 106:108 gd' PAR_NAME, datafile,'.mif PA_' PAR_NAME, datafile,'.mif']);
    end
    
    %match up gradient files (bvals, bvecs) with main dwi dataset.
    unix(['mrconvert -fslgrad ', PAR_NAME, datafile,'.bvec ' PAR_NAME, datafile,'.bval cgd' PAR_NAME, datafile,'.mif bcgd', PAR_NAME, datafile,'.mif']);
    
    %create a copy in NIFTI format
    unix(['mrconvert bcgd', PAR_NAME, datafile,'.mif bcgd', PAR_NAME, datafile, '.nii']);
    
    %---------------------------------------------------------------------%
    %Step 3-4: Field distortion (topup, eddy)
    
    %extract all B0s from the dataset
    unix(['dwiextract -bzero bcgd' PAR_NAME, datafile,'.mif AP_' PAR_NAME, datafile,'.mif']);
    
    %place a-p and p-a images into one file
    unix(['mrcat AP_' PAR_NAME, datafile,'.mif ' 'PA_' PAR_NAME, datafile,'.mif ' 'allB0s_' PAR_NAME, datafile,'.mif']);
    
    %choose best pair using the 'BestB0' function
    BU_used = BestB0(PAR_NAME, datafile, NumBDs, startdir);
    
    %run topup, eddy (w/ -repol)
    if BU_used ~= 1
        %run eddy with gradient edit
        unix(['dwifslpreproc -fslgrad ', PAR_NAME, datafile,'.bvec ' PAR_NAME, datafile,'.bval bbcgd' PAR_NAME, datafile, '.mif ebbcgd' PAR_NAME, datafile, '.mif -rpe_pair -pe_dir AP -se_epi TUB0s_' PAR_NAME, datafile, '.mif -eddy_options " --repol --ol_nstd=4" -readout_time 0.07']);
    else
        %don't need to apply gradient edit with eddy
        unix(['dwifslpreproc bbcgd' PAR_NAME, datafile, '.mif ebbcgd' PAR_NAME, datafile, '.mif -rpe_pair -pe_dir AP -se_epi TUB0s_' PAR_NAME, datafile, '.mif -eddy_options " --repol --ol_nstd=4" -readout_time 0.07']);
    end
    %---------------------------------------------------------------------%
    %Step 5: Estimate a brain mask (initial)
    
    %estimate brain mask
    unix(['dwi2mask ebbcgd' PAR_NAME, datafile, '.mif initial_mask_' PAR_NAME, datafile, '.mif']);
    
    %create a copy in NIFIT format
    unix(['mrconvert ebbcgd', PAR_NAME, datafile,'.mif ebbcgd', PAR_NAME, datafile, '.nii']); %eddy corrected data
    unix(['mrconvert initial_mask_', PAR_NAME, datafile,'.mif initial_mask_', PAR_NAME, datafile, '.nii']); %brain mask
    
    
    %---------------------------------------------------------------------%
    %Step 6: Bias field correction (B0s)
    unix(['dwibiascorrect ants -mask initial_mask_' PAR_NAME, datafile, '.mif ebbcgd' PAR_NAME, datafile, '.mif febbcgd' PAR_NAME, datafile, '.mif']);
    
    %create a copy in NIFIT format
    unix(['mrconvert febbcgd', PAR_NAME, datafile,'.mif febbcgd', PAR_NAME, datafile, '.nii']);
    
    %generate a new brain mask after bias field correction has been
    %applied.
    unix(['dwi2mask -clean_scale 3 febbcgd' PAR_NAME, datafile, '.mif brain_mask_' PAR_NAME, datafile, '.mif']);
    
    %create a copy in NIFIT format
    unix(['mrconvert brain_mask_', PAR_NAME, datafile,'.mif brain_mask_', PAR_NAME, datafile, '.nii']);
    
    %copy dwi file over and rename
    copyfile(['febbcgd', PAR_NAME, datafile,'.mif'], [PAR_NAME, datafile,'.mif']);
    movefile([PAR_NAME, datafile, '.mif'], [startdir,'/test/derivatives/', groupname, '/preprocessed_dwi']);
    
    %copy brain mask file over
    copyfile(['brain_mask_', PAR_NAME, datafile,'.mif'], [PAR_NAME, datafile,'.mif']);
    movefile([PAR_NAME, datafile, '.mif'], [startdir,'/test/derivatives/', groupname, '/brain_mask']);
    
    
end


