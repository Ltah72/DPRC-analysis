%Run Fmriprep (version 20.2.1) via Docker on subjects. For the anatomical 
%files, this includes: bias field correction, skull-stripping, brain 
%surface reconstruction (recon-all,FreeSurfer), brain mask estimation (ANTs 
%+ FreeSurfer), spatial normalisation (antsRegistration, ANTs), brain 
%tissue segmentation (FAST, FSL). For the functional files, this includes:
%slice time correction (3dtshift, AFNI), motion correction (mcflirt, FSL),
%coregistration to T1w with 6 dof (bbregister, FreeSurfer), and motion
%correction transformations (antsApplyTransforms, ANTs). 

%If it hasn't been done already, you may need to edit the files for
%fmriprep (e.g. .json files), which I have created a script for this,
%called ModifyFilesForFmriprep.m.

%Esteban, O., Markiewicz, C. J., Blair, R. W., Moodie, C. A., Isik, A. I., Erramuzpe, A., … Gorgolewski, K. J. (2019). fMRIPrep: a robust preprocessing pipeline for functional MRI.
%https://fmriprep.org/en/stable/

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 10/01/20

%Define fmriprep directory, so that it used:
%FmriprepDirectory = input('Please enter fmriprep directory:', 's');
FmriprepDirectory = '/data/USERS/LENORE/fmriprep_test/';

%go into fmriprep sourcedata directory
cd([FmriprepDirectory 'sourcedata/']);

%specify participants
participants = uipickfiles;

for i = 1:length(participants)
    
    [upper_path, PAR_NAME, ~] = fileparts(participants{1,i});
    
    %---------------------------------------------------------------------%
    %Step 1: Preprocess anatomical images (T1w and T2 FLAIR images)
    
    %a) Gibbs ringing correction (mrdegibbs) on all anat files (T1w, T2w
    %blade, + T2w FLAIR). Use 'force' to overwrite the raw file. 
    cd ([PAR_NAME '/anat']);
    unix(['mrdegibbs ' PAR_NAME '_T1w.nii ' PAR_NAME '_T1w.nii -force']);
    unix(['mrdegibbs ' PAR_NAME '_T2w.nii ' PAR_NAME '_T2w.nii -force']);
    unix(['mrdegibbs ' PAR_NAME '_FLAIR.nii ' PAR_NAME '_FLAIR.nii -force']);
    
    cd([FmriprepDirectory 'sourcedata/');
    
    %b) Run fmriprep to preprocess anat + func images
    unix(['sudo /home/ubuntu/.local/bin/./fmriprep-docker ' FmriprepDirectory 'sourcedata ' FmriprepDirectory 'derivatives participant --participant-label ' PAR_NAME ' --output-spaces func --fs-license-file /SOFTWARE/freesurfer/license.txt']);

end