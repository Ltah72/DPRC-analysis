%This is the tract of interest (TOI) pipeline (e.g. superior longitudinal
%fasciculus (SLF)) which will select certain TOIs using manually generated 
%ROI masks and data from the FBA pipeline. Unlike the other pipelines, this 
%pipeline will require more manual user input (e.g. going into the viewing 
%gui, mrview, with the tract(s) that you select. Also, unlike my other 
%scripted pipelines, there is no 'formal manual' for this from mrtrix, 
%and many of the steps are referenced by the mrtrix community forum and 
%logic of steps/analysis - see in my created manual for more details.

% Steps:
% 1. Compute track density image (TDI)
% 2. Create mask for the track of interest (TOI)
% 3. Edit your tracks with TOI mask as an input
% 4. Manually include/exclude fibres on your TOI if needed
% 5. Create a fixel mask with your TOI
% 6. Threshold the TOI mask
% 7. Compute fixel-based metrics against TOI mask per each participant 

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 04/08/20

clc;
clear all;
close all;


%define/add pathways
startdir = '/data/USERS/LENORE';

%Script directory is defined, so that it can be added to path below:
ScriptDirectory = '/data/USERS/LENORE/scripts/dprc/diffusion';

%should be the same groupname from what the user analysed in the FBA script.  
groupname = input('Which pre-processed group / study do you want to continue to analyse?: ', 's');

%Create TOI folder
TOI = input('Please name a tract of interest (TOI) that you want to analyse: ', 's');

mkdir([startdir,'/derivatives/diff_data/', groupname, '/' TOI]);

%go into group folder
cd([startdir '/derivatives/diff_data/', groupname]);

%define variables
participants = dir(fullfile('preprocessed_dwi', '*.mif'));

%Add your script and all necessary files (e.g. data, functions) to path
addpath(genpath(startdir));
addpath(genpath(ScriptDirectory));


%copy needed files into the TOI directory
copyfile ([startdir '/derivatives/diff_data/' groupname '/wmfod_template.mif'], [startdir,'/derivatives/diff_data/', groupname '/' TOI]);
copyfile ([startdir '/derivatives/diff_data/' groupname '/tracks_2_million_sift.tck'], [startdir,'/derivatives/diff_data/', groupname '/' TOI]);

%go into TOI directory
cd ([startdir,'/derivatives/diff_data/', groupname, '/' TOI]);


%-------------------------------------------------------------------------%
%Step 1: Create a track map (TDI) in order to see the tracts a bit better.

unix(['tckmap -template wmfod_template.mif -dec tracks_2_million_sift.tck tdi_image.mif']);


%-------------------------------------------------------------------------%
%Step 2: Create a mask of your TOI - manually select ROI(s) in mrview for 
%this. Use TDI image and other anatomical references as guides to do this.

%%%%% go into mrview and create a manual mask for your ROI/TOI %%%%%

% Once you created your ROI mask, move it into the current processing
%directory.

mask = input('Please provide the full name of the TOI mask (e.g. mask_corpus_callosum.mif): ', 's');


%-------------------------------------------------------------------------%
%Step 3: Extract TOI from the previously generated whole-brain SIFT 
%tractography file
unix(['tckedit tracks_2_million_sift.tck -include ' mask ' ' TOI '_sift.tck']);
%e.g. of full command:
%unix(['tckedit tracks_2_million_sift.tck -include mask_corpus_callosum.mif CorpusCallosum_sift.tck']);


%-------------------------------------------------------------------------%
%Step 4: Manual inclusion/exclusion of fibres (if needed)
%%%%% go into mrview for this if needed %%%%%, 
% Make your edits to the line below %

%unix(['tckedit ' TOI '_sift.tck -exclude -1.9,-21.5,-40.1,30 -exclude 10.3,-14.8,3.9,15 -exclude -12.8,-18.3,3.2,15 ' TOI '_sift_edit.tck -force']);

%e.g. of full command:
%unix(['tckedit CorpusCallosum_sift.tck -exclude -1.9,-21.5,-40.1,30 -exclude 10.3,-14.8,3.9,15 -exclude -12.8,-18.3,3.2,15 CorpusCallosum_sift_edit.tck']);


%-------------------------------------------------------------------------%
%Step 5: Create fixel mask with your TOI - make a fixel-wise TDI
ManualEdits = input('After viewing your TOI, did you make additional edits? (y or n): ', 's');

%with edits
if ManualEdits == 'y'
    unix(['tck2fixel ' TOI '_sift_edit.tck ' startdir,'/derivatives/diff_data/', groupname '/fixel_directory/fixel_mask output_' TOI '_TOI_fixel_directory ' TOI '_fixel_mask.mif']);
%no edits done    
elseif ManualEdits == 'n'
    unix(['tck2fixel ' TOI '_sift.tck ' startdir,'/derivatives/diff_data/', groupname '/fixel_directory/fixel_mask output_' TOI '_TOI_fixel_directory ' TOI '_fixel_mask.mif']);
end

%e.g. of full command:
%unix(['tck2fixel CorpusCallosum_sift_edit.tck ' startdir,'/derivatives/diff_data/', groupname '/fixel_directory/fixel_mask output_CorpusCallosum_TOI_fixel_directory CorpusCallosum_fixel_mask.mif']);


%-------------------------------------------------------------------------%
%Step 6: Threshold the TOI mask
unix(['mrthreshold -abs 0.95 output_' TOI '_TOI_fixel_directory/' TOI '_fixel_mask.mif output_' TOI '_TOI_fixel_directory/' TOI '_fixel_mask_thr.mif']);


%-------------------------------------------------------------------------%
%Step 7: Compute fixel-based metrics (FD, FC, FDC) with the TOI fixel mask per each participant.
CreateTOIFBAMetricFiles(participants, TOI, startdir, groupname);



