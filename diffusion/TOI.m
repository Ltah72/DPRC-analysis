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
% 2-3. Options: (1)Generate seedpoints (see below first script) (2)Create 
% mask for the track of interest (TOI)
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
derivdir = '/data/USERS/LENORE/derivatives';

%Script directory is defined, so that it can be added to path below:
ScriptDirectory = '/data/USERS/LENORE/scripts/dprc/diffusion';

%should be the same groupname from what the user analysed in the FBA script.  
groupname = input('Which pre-processed group / study do you want to continue to analyse?: ', 's');

%choose time period
period = input('Which time period do you want to analyse, e.g. F0, F2, all, etc?: ', 's');

%Create TOI folder
TOI = input('Please name a tract of interest (TOI) that you want to analyse: ', 's');

mkdir([derivdir,'/groups/' period, '/diff_data/', groupname, '/' TOI]);

%go into group folder
cd([derivdir '/groups/' period, '/diff_data/', groupname]);

%define variables
participants = dir(fullfile('preprocessed_dwi', '*.mif'));

%Add your script and all necessary files (e.g. data, functions) to path
addpath(genpath(derivdir));
addpath(genpath(ScriptDirectory));


%copy needed files into the TOI directory
copyfile ([derivdir '/groups/' period, '/diff_data/' groupname '/wmfod_template.mif'], [derivdir,'/groups/' period, '/diff_data/', groupname '/' TOI]);
copyfile ([derivdir '/groups/' period, '/diff_data/' groupname '/tracks_2_million_sift.tck'], [derivdir,'/groups/' period, '/diff_data/', groupname '/' TOI]);

%go into TOI directory
cd ([derivdir,'/groups/' period, '/diff_data/', groupname, '/' TOI]);


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
    unix(['tck2fixel ' TOI '_sift_edit.tck ' derivdir,'/groups/' period, '/diff_data/', groupname '/fixel_directory/fixel_mask output_' TOI '_TOI_fixel_directory ' TOI '_fixel_mask.mif']);
%no edits done    
elseif ManualEdits == 'n'
    unix(['tck2fixel ' TOI '_sift.tck ' derivdir,'/groups/' period, '/diff_data/', groupname '/fixel_directory/fixel_mask output_' TOI '_TOI_fixel_directory ' TOI '_fixel_mask.mif']);
end

%e.g. of full command:
%unix(['tck2fixel CorpusCallosum_sift_edit.tck ' startdir,'/derivatives/diff_data/', groupname '/fixel_directory/fixel_mask output_CorpusCallosum_TOI_fixel_directory CorpusCallosum_fixel_mask.mif']);


%-------------------------------------------------------------------------%
%Step 6: Threshold the TOI mask
unix(['mrthreshold -abs 0.95 output_' TOI '_TOI_fixel_directory/' TOI '_fixel_mask.mif output_' TOI '_TOI_fixel_directory/' TOI '_fixel_mask_thr.mif']);


%-------------------------------------------------------------------------%
%Step 7: Compute fixel-based metrics (FD, FC, FDC) with the TOI fixel mask per each participant.
CreateTOIFBAMetricFiles(participants, TOI, derivdir, groupname, period);







%%%% Test SLF script below: %%%%

%SLF tract parameters

%left SLF
unix(['tckgen wmfod_template.mif -select 10k -maxlen 250 -minlen 10 -angle 45 -power 1.0 -cutoff 0.15 -seed_sphere -34.19,-30.87,21.51,3 -exclude 50,92,102,2 -exclude -36.41,-71.2,16.21,2 -exclude -36.41,-45.49,37.64,2 -exclude -2.33,-46.72,5.07,3 -exclude 43,96,68,3 -exclude -16.59,-64.65,17.04,2 -exclude -16.54,-68.61,15.82,2 -exclude -29.79,-60.68,17.97,1 -exclude -27.44,-68.82,-2.9,1 SLF_track_L.tck -force']);

%right SLF
unix(['tckgen wmfod_template.mif -select 10k -maxlen 250 -minlen 10 -angle 45 -power 1.0 -cutoff 0.15 -seed_sphere 34.93,-29.5,23.28,3 -exclude 17.76,-13.72,31.91,2 -exclude 28.78,-22.39,12.63,2 SLF_track_R.tck -force']);


%merge the left and right tracks as one track file. 
unix(['tckedit SLF_track_L.tck SLF_track_R.tck SLF_track_whole.tck']);

%create fixel mask with your track file. 
unix(['tck2fixel SLF_track_whole.tck fixel_mask output_SLF_whole_fixel_mask SLF_whole_fixel_mask.mif'])

%threshold your mask (whole)
unix(['mrthreshold -abs 0.5 output_SLF_whole_fixel_mask/SLF_whole_fixel_mask.mif output_SLF_whole_fixel_mask/SLF_whole_fixel_mask_thr.mif']);



%%%Left SLF %%%
%create fixel mask with your track file (left SLF). 
unix(['tck2fixel SLF_track_L.tck fixel_mask output_SLF_track_L_fixel_mask SLF_track_L_fixel_mask.mif'])
%threshold your mask (left SLF)
unix(['mrthreshold -abs 0.5 output_SLF_track_L_fixel_mask/SLF_L_fixel_mask.mif output_SLF_track_L_fixel_mask/SLF_track_L_fixel_mask_thr.mif']);

%%%Right SLF %%%
%create fixel mask with your track file (left SLF). 
unix(['tck2fixel SLF_track_R.tck fixel_mask output_SLF_track_R_fixel_mask SLF_track_R_fixel_mask.mif'])
%threshold your mask (right SLF)
unix(['mrthreshold -abs 0.5 output_SLF_track_R_fixel_mask/SLF_R_fixel_mask.mif output_SLF_track_R_fixel_mask/SLF_track_R_fixel_mask_thr.mif']);


%Compute fixel-based metrics (FD, FC, FDC) with the TOI fixel mask per each participant.
CreateTOIFBAMetricFiles(participants, TOI, derivdir, groupname);




















