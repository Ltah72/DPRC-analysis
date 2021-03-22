%For running FreeSurfer's recon-all command. Inputs are the participant's
%anatomical files (T1w and T2w FLAIR), and the output folder will be in the
%FreeSurfer's subject directory ($SUBJECTS_DIR). 


%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 12/01/21


%go into participants' sourcedata directory
cd (['/data/sourcedata/' period]);

%choose participants to run recon-all on
participants = uipickfiles;

for i = 1:length(participants)
    
    [upper_path, PAR_NAME, ~] = fileparts(participants{1,i});
    
    %go into participant's anat folder in sourcedata directory 
    cd (['/data/sourcedata/' period '/' PAR_NAME '/anat/']);
    
    unix(['recon-all -i ' PAR_NAME '_T1w.nii -FLAIR ' PAR_NAME '_FLAIR.nii -s ' PAR_NAME ' -all -parallel -openmp 25']);

    %copy the recon-all output into the fmriprep derivatives directory
    unix(['sudo cp -r $SUBJECTS_DIR/' PAR_NAME ' /data/USERS/LENORE/fmriprepped_data/derivatives/freesurfer/']); 
end 