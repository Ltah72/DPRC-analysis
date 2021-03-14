function RunEddyQuad(startdir, ScriptDirectory, PAR_NAME, datafile, period)
%Run eddy_quad on all participants. This will generate a report of the
%quality control of eddy for all participants. It will create an 
%eddy_quad.qc folder, which will be used later to input into a group eddy 
%qc analysis (via eddy_squad). This will ultimately help us to determine 
%the outliers of the participants. 

%       
%Inputs(5): startdir = start directory that you defined in the script -
%           where the data will be stored.
%           ScriptDirectory = where you defined where all your diffusion
%           scripts are
%           PAR_NAME = participant ID - current participant for input
%           datafile = dwi definition file (i.e. '_acq_data_dwi')
%           period = time period of the participant MRI scans


%Output(none): eddy_quad will generate a folder per each participant, 
%              which will each hold the eddy qc files, needed for the group
%              qc analysis (eddy_squad) at the end of the script pipeline. 

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 21/08/20


%go into eddyqc folder and move and rename files (with same base name
%as eddy datafile to main participant folder
cd eddyqc;

%define variables
files = dir('eddy*');
%base name for eddy quad correction, as fsl requires this. 
newFilePart = 'eddy_quad.';

for i = 1:length(files)
    [~,filename,extension] = fileparts(files(i).name);
    %add in the base name at the beginning to all of the files
    movefile([filename extension], [newFilePart filename extension]);
end

%go back into main participant folder
cd([startdir '/derivatives/' period, '/diff_data/' PAR_NAME, '/dwi/']);

%make a copy and rename the eddy corrected datafile, and place into the
%eddyqc folder
copyfile(['ebbcgd', PAR_NAME, datafile,'.nii'], 'eddyqc/eddy_quad.nii');

%go into the eddyqc folder with all of the files
cd eddyqc;

%perform eddy quad on the participant
unix(['eddy_quad eddy_quad -idx ' ScriptDirectory '/files/index.txt -par ' ScriptDirectory '/files/acqparams.txt -m ' startdir '/derivatives/' period, '/diff_data/' PAR_NAME '/dwi/brain_mask_' PAR_NAME, datafile, '.nii -b ' startdir '/derivatives/' period, '/diff_data/' PAR_NAME '/dwi/' PAR_NAME, datafile, '.bval']);

%go back into main participant folder
cd([startdir '/derivatives/' period, '/diff_data/' PAR_NAME, '/dwi/']);

end

