%This script will look at the fMRI timeseries between the superior 
%longitudinal fasciculus (SLF) track end points.


%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 18/01/22

clc;
clear all;
close all;

%set up & add pathways
%Define script directory, so that it can be added to path below:
%ScriptDirectory = input('Please enter script directory:', 's');
ScriptDirectory = '/data/USERS/LENORE/scripts/dprc/diffusion';
%Choose directory where the denoised fMRI data is.
%cd /data/USERS/LENORE
cd /home/ltah262/ressci202000017-dprc_diff_vis/NECTAR_data/LENORE/derivatives/fMRI_denoised/unzipped_files
addpath /data/SOFTWARE/Pls/Pls/plsgui
addpath /data/SOFTWARE/Pls/Pls/plscmd
addpath(genpath(ScriptDirectory));

%define constant variables
%fMRI image resolution
voxel_size = [2 2 2];
%origin point (centre) for the fMRI image. This can be found when loading the image in spm -
%this also remains the same for all of the images in MNI space
origin = [46 64 37];
%Frequency of the signal, defined in Hz (i.e. collect 1 fMRI image in 735 ms)
Fs = 1/.735;

%Create text files with column names:
%for timeseries file:
fid = fopen('Ptpt_time_series_values.txt', 'w');
if (fid == -1)
    disp('Error in opening the file.')
else
    %print the correlation and asssociated p-value to a textfile for each participant.
    fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s', 'ParticipantID', 'bp_ant_ts_SLF1_L', 'bp_post_ts_SLF1_L', 'bp_ant_ts_SLF2_L', 'bp_post_ts_SLF2_L', 'bp_ant_ts_SLF3_L', 'bp_post_ts_SLF3_L', 'bp_ant_ts_SLF1_R', 'bp_post_ts_SLF1_R', 'bp_ant_ts_SLF2_R', 'bp_post_ts_SLF2_R', 'bp_ant_ts_SLF3_R', 'bp_post_ts_SLF3_R');
    fprintf(fid, '\n');
end
fclose(fid);

%for correlation + p-value file:
fid = fopen('Ptpt_corr_pval.txt', 'w');
if (fid == -1)
    disp('Error in opening the file.')
else
    %print the correlation and asssociated p-value to a textfile for each participant.
    fprintf(fid, '%s   %s %s %s %s %s %s %s %s %s %s  %s %s', 'ParticipantID', 'r_SLF1_L', 'p_SLF1_L', 'r_SLF2_L', 'p_SLF2_L', 'r_SLF3_L', 'p_SLF3_L', 'r_SLF1_R', 'p_SLF1_R', 'r_SLF2_R', 'p_SLF2_R', 'r_SLF3_R', 'p_SLF3_R');
    fprintf(fid, '\n');
end
fclose(fid);



%choose participants for analysis 
participants = uipickfiles;

for i = 1:length(participants)
    
    [upper_path, PAR_NAME, ~] = fileparts(participants{1,i});
    PAR_NAME = PAR_NAME(1:15);
    
    %load in the clean + preprocessed fMRI image
    brain = load_nii([PAR_NAME '_task-rest_run-1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold_NR.nii']);
    
    %anterior white matter tract xyz coordinates (WorldSpace in mm) for each SLF tract
    %left SLF tract
    ant_xyzmm_SLF1_L = [-12.6 28.6 37.4];
    ant_xyz_SLF1_L = round(double(ant_xyzmm_SLF1_L./(ones(size(ant_xyzmm_SLF1_L,1),1)*voxel_size)+ones(size(ant_xyzmm_SLF1_L,1),1)*origin));
    ant_xyzmm_SLF2_L = [-28.8 21 33.4];
    ant_xyz_SLF2_L = round(double(ant_xyzmm_SLF2_L./(ones(size(ant_xyzmm_SLF2_L,1),1)*voxel_size)+ones(size(ant_xyzmm_SLF2_L,1),1)*origin));
    ant_xyzmm_SLF3_L = [-35.6 36.6 -0.6];
    ant_xyz_SLF3_L = round(double(ant_xyzmm_SLF3_L./(ones(size(ant_xyzmm_SLF3_L,1),1)*voxel_size)+ones(size(ant_xyzmm_SLF3_L,1),1)*origin));
    %right SLF tract
    ant_xyzmm_SLF1_R = [17.2 30.6 35];
    ant_xyz_SLF1_R = round(double(ant_xyzmm_SLF1_R./(ones(size(ant_xyzmm_SLF1_R,1),1)*voxel_size)+ones(size(ant_xyzmm_SLF1_R,1),1)*origin));
    ant_xyzmm_SLF2_R = [30, 25, 30]; %guessed coordinates for now for SLF2_R
    ant_xyz_SLF2_R = round(double(ant_xyzmm_SLF2_R./(ones(size(ant_xyzmm_SLF2_R,1),1)*voxel_size)+ones(size(ant_xyzmm_SLF2_R,1),1)*origin));
    ant_xyzmm_SLF3_R = [36 34.8 -6.2];
    ant_xyz_SLF3_R = round(double(ant_xyzmm_SLF3_R./(ones(size(ant_xyzmm_SLF3_R,1),1)*voxel_size)+ones(size(ant_xyzmm_SLF3_R,1),1)*origin));
    
    %posterior white matter tract xyz coordinates (WorldSpace in mm) for each SLF tract
    %left SLF tract
    post_xyzmm_SLF1_L = [-15.4 -72 38.6];
    post_xyz_SLF1_L = round(double(post_xyzmm_SLF1_L./(ones(size(post_xyzmm_SLF1_L,1),1)*voxel_size)+ones(size(post_xyzmm_SLF1_L,1),1)*origin));
    post_xyzmm_SLF2_L = [-35.6 -75.4 23.4];
    post_xyz_SLF2_L = round(double(post_xyzmm_SLF2_L./(ones(size(post_xyzmm_SLF2_L,1),1)*voxel_size)+ones(size(post_xyzmm_SLF2_L,1),1)*origin));
    post_xyzmm_SLF3_L = [-49.8 -43.8 38.4];
    post_xyz_SLF3_L = round(double(post_xyzmm_SLF3_L./(ones(size(post_xyzmm_SLF3_L,1),1)*voxel_size)+ones(size(post_xyzmm_SLF3_L,1),1)*origin));
    %right SLF tract
    post_xyzmm_SLF1_R = [14.6 -64.6 48.6];
    post_xyz_SLF1_R = round(double(post_xyzmm_SLF1_R./(ones(size(post_xyzmm_SLF1_R,1),1)*voxel_size)+ones(size(post_xyzmm_SLF1_R,1),1)*origin));
    post_xyzmm_SLF2_R = [30, -70, 28]; %guessed coordinates for now for SLF2_R
    post_xyz_SLF2_R = round(double(post_xyzmm_SLF2_R./(ones(size(post_xyzmm_SLF2_R,1),1)*voxel_size)+ones(size(post_xyzmm_SLF2_R,1),1)*origin));
    post_xyzmm_SLF3_R = [42.8 -46.9 30.1];
    post_xyz_SLF3_R = round(double(post_xyzmm_SLF3_R./(ones(size(post_xyzmm_SLF3_R,1),1)*voxel_size)+ones(size(post_xyzmm_SLF3_R,1),1)*origin));
    
    %fit data into matlab format
    %left tract
    ant_ts_SLF1_L = squeeze(brain.img(ant_xyz_SLF1_L(1), ant_xyz_SLF1_L(2), ant_xyz_SLF1_L(3), :)); %returns 0s...why?
    post_ts_SLF1_L = squeeze(brain.img(post_xyz_SLF1_L(1), post_xyz_SLF1_L(2), post_xyz_SLF1_L(3), :));
    ant_ts_SLF2_L = squeeze(brain.img(ant_xyz_SLF2_L(1), ant_xyz_SLF2_L(2), ant_xyz_SLF2_L(3), :));
    post_ts_SLF2_L = squeeze(brain.img(post_xyz_SLF2_L(1), post_xyz_SLF2_L(2), post_xyz_SLF2_L(3), :));
    ant_ts_SLF3_L = squeeze(brain.img(ant_xyz_SLF3_L(1), ant_xyz_SLF3_L(2), ant_xyz_SLF3_L(3), :));
    post_ts_SLF3_L = squeeze(brain.img(post_xyz_SLF3_L(1), post_xyz_SLF3_L(2), post_xyz_SLF3_L(3), :));
    %right tract
    ant_ts_SLF1_R = squeeze(brain.img(ant_xyz_SLF1_R(1), ant_xyz_SLF1_R(2), ant_xyz_SLF1_R(3), :));
    post_ts_SLF1_R = squeeze(brain.img(post_xyz_SLF1_R(1), post_xyz_SLF1_R(2), post_xyz_SLF1_R(3), :));
    ant_ts_SLF2_R = squeeze(brain.img(ant_xyz_SLF2_R(1), ant_xyz_SLF2_R(2), ant_xyz_SLF2_R(3), :));
    post_ts_SLF2_R = squeeze(brain.img(post_xyz_SLF2_R(1), post_xyz_SLF2_R(2), post_xyz_SLF2_R(3), :));
    ant_ts_SLF3_R = squeeze(brain.img(ant_xyz_SLF3_R(1), ant_xyz_SLF3_R(2), ant_xyz_SLF3_R(3), :));
    post_ts_SLF3_R = squeeze(brain.img(post_xyz_SLF3_R(1), post_xyz_SLF3_R(2), post_xyz_SLF3_R(3), :));
    
    %bandpass filter the data
    %left tract
    bp_ant_ts_SLF1_L = bandpass(ant_ts_SLF1_L,[.01, .1], Fs);
    bp_ant_ts_SLF2_L = bandpass(ant_ts_SLF2_L,[.01, .1], Fs);
    bp_ant_ts_SLF3_L = bandpass(ant_ts_SLF3_L,[.01, .1], Fs);
    bp_post_ts_SLF1_L = bandpass(post_ts_SLF1_L,[.01, .1], Fs);
    bp_post_ts_SLF2_L = bandpass(post_ts_SLF2_L,[.01, .1], Fs);
    bp_post_ts_SLF3_L = bandpass(post_ts_SLF3_L,[.01, .1], Fs);
    %right tract
    bp_ant_ts_SLF1_R = bandpass(ant_ts_SLF1_R,[.01, .1], Fs);
    bp_ant_ts_SLF2_R = bandpass(ant_ts_SLF2_R,[.01, .1], Fs);
    bp_ant_ts_SLF3_R = bandpass(ant_ts_SLF3_R,[.01, .1], Fs);
    bp_post_ts_SLF1_R = bandpass(post_ts_SLF1_R,[.01, .1], Fs);
    bp_post_ts_SLF2_R = bandpass(post_ts_SLF2_R,[.01, .1], Fs);
    bp_post_ts_SLF3_R = bandpass(post_ts_SLF3_R,[.01, .1], Fs);
    
    %put anterior + posterior values from timeseries per each participant (490
    %values per participant) on a file.
    %combine the bp timeseries into one matrix per participant
    pt_ts = [bp_ant_ts_SLF1_L, bp_post_ts_SLF1_L, bp_ant_ts_SLF2_L, bp_post_ts_SLF2_L, bp_ant_ts_SLF3_L, bp_post_ts_SLF3_L, bp_ant_ts_SLF1_R, bp_post_ts_SLF1_R, bp_ant_ts_SLF2_R, bp_post_ts_SLF2_R, bp_ant_ts_SLF3_R, bp_post_ts_SLF3_R];
    
    %write matrix to text file for timeseries
    writematrix(pt_ts, 'pt_ts.txt');
    
    %write matrix to text file for the current participant's ID
    current_pt = [repmat(PAR_NAME,490,1)];
    writematrix(current_pt, 'current_pt.txt');
    
    %combine the participant ID with their timeseries matrix file
    unix(['paste current_pt.txt pt_ts.txt > current_full_pt_ts.txt']);
    
    %append participant text files
    unix('cat current_full_pt_ts.txt >> Ptpt_time_series_values.txt');
    
    %get the correlation metric and associated p-value, per each tract
    %left tract
    [r_SLF1_L, p_SLF1_L] = corr(bp_ant_ts_SLF1_L, bp_post_ts_SLF1_L);
    [r_SLF2_L, p_SLF2_L] = corr(bp_ant_ts_SLF2_L, bp_post_ts_SLF2_L);
    [r_SLF3_L, p_SLF3_L] = corr(bp_ant_ts_SLF3_L, bp_post_ts_SLF3_L);
    %right tract
    [r_SLF1_R, p_SLF1_R] = corr(bp_ant_ts_SLF1_R, bp_post_ts_SLF1_R);
    [r_SLF2_R, p_SLF2_R] = corr(bp_ant_ts_SLF2_R, bp_post_ts_SLF2_R);
    [r_SLF3_R, p_SLF3_R] = corr(bp_ant_ts_SLF3_R, bp_post_ts_SLF3_R);
    
    %put correlation + pvalue on a file
    fid = fopen('Ptpt_corr_pval.txt', 'a+');
    if (fid == -1)
        disp('Error in opening the file.')
    else
        %print the correlation and asssociated p-value to a textfile for each participant.
        fprintf(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f', PAR_NAME, r_SLF1_L, p_SLF1_L, r_SLF2_L, p_SLF2_L, r_SLF3_L, p_SLF3_L, r_SLF1_R, p_SLF1_R, r_SLF2_R, p_SLF2_R, r_SLF3_R, p_SLF3_R);
        fprintf(fid, '\n');
    end
    fclose(fid);
end

%plot the correlation as a scatterplot
plot(bp_ant_ts_SLF1_L, bp_post_ts_SLF1_L, '.k', 'MarkerSize', 10)

%plot + compare time series (w/ and w/o bandpass filter)
plot(ant_ts_SLF1_L);
hold on
plot(bp_ant_ts_SLF1_L+mean(ant_ts_SLF1_L));

plot(bp_ant_ts_SLF1_L);
hold on
plot(bp_post_ts_SLF1_L);



