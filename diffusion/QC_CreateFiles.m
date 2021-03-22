function QC_CreateFiles(startdir, period)
%The purpose of this function is to create all of the collated quality 
%control (qc) files in the dwiqc directory. This directory will hold all 
%participants' data together on one file for qc checks/analysis.So far, 
%the data qc includes checking for the best B0 and checking for the
%total movement with the eddy correction. 

%Inputs (2): period = time period of the participant MRI scans
%            startdir = start directory that you defined in the script -


%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 18/08/20


cd ([startdir '/derivatives/' period '/diff_data/dwiqc/']);

%create BestB0 text file with header line
fid = fopen('BestB0.txt', 'w');
if (fid == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid, '%s     %s %s %s', 'Participant', 'B0_status', 'BU_used', 'BD_used');
    fclose(fid);
end


%create eddy_qc text file, for looking at movement, with header line
fid2 = fopen('eddyqc_movement_all_vols.txt', 'w');
if (fid2 == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid2, '%s     %s %s %s', 'Participant', 'mvmt_vol_1', 'mvmt_prev_vol Volume');
    fclose(fid2);
end


%create eddy_qc text file, for looking average movement per participant, with header line
fid3 = fopen('eddyqc_movement_average.txt', 'w');
if (fid3 == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid3, '%s     %s %s', 'Participant', 'mvmt_vol_1', 'mvmt_prev_vol');
    fclose(fid3);
end


%create eddy_qc text file for looking at outliers, with header line
fid4 = fopen('eddyqc_outliers.txt', 'w');
if (fid4 == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid4, '%s     %s %s', 'Participant', 'num_outliers', 'outlier_%');
    fclose(fid4);
end


%create eddy_qc text file for looking at SNR and CNR, with header line
fid5 = fopen('eddyqc_SNR&CNR.txt', 'w');
if (fid5 == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid5, '%s         %s  %s  %s', 'Participant', 'SNR_b0', 'CNR_b1000', 'CNR_b2000');
    fclose(fid5);
end


%create a list for running eddy_squad (group analysis)
fid6 = fopen('EddySQUADList.txt', 'w');
if (fid6 == -1)
    disp('Error in creating the text file.')
else
    fclose(fid6);
end 


end

