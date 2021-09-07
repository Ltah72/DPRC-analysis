function eddyqc_ToText(PAR_NAME, derivdir, period)
%This function will copy the information and values from the eddy qc files,
%such movement and outlier values from all participants, and will place
%these all into one text file for the user to collate together into a 
%visual graph.

%Inputs(3): PAR_NAME = participant ID - current participant for input
%           startdir = start directory that you defined in the script -
%           where the data will be stored.
%           period = time period of the participant MRI scans

%Output(none): You are writing values to the file.

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 18/08/20


%go into eddyqc folder of each participant
cd eddyqc;

%get data from the specific movement file
fid = fopen('eddy_quad.eddy_movement_rms', 'r');

if (fid == -1)
    disp('Error in opening the eddy_movement_rms file.')
else
    %read in the values from the text files
    eddy_movement = fscanf(fid, '%f');
    
    %put values into a separate matrix
    eddy_movement_c1 = eddy_movement(1:2:212);
    eddy_movement_c2 = eddy_movement(2:2:212);
    
    %take average movement from each column
    average_eddy_movement_c1 = mean(eddy_movement_c1);
    average_eddy_movement_c2 = mean(eddy_movement_c2);
    
    fclose(fid);
    
    %go into the main qc folder
    cd ([derivdir '/groups/' period '/diff_data/dwiqc/']);
    
    %edit / add the participant's values to the file
    fid2 = fopen('eddyqc_movement_all_vols.txt', 'a+');
    if (fid2 == -1)
        disp('Error in opening the eddyqc_movement_all_vols.txt file.')
    else
        for i = 1:106
            %print values on the eddyqc_movement_all_vols.txt file
            fprintf(fid2, '\n');
            fprintf(fid2, '%s     %f %f %.f', PAR_NAME, eddy_movement_c1(i,1), eddy_movement_c2(i,1), i);
        end
    end
    
    fclose(fid2);
    
    %record the average movement of all vol per participant on the file.
    fid3 = fopen('eddyqc_movement_average.txt', 'a+');
    if (fid3 == -1)
        disp('Error in opening the eddyqc_movement_average.txt file.')
    else
        %print values on the eddyqc_movement_average.txt file
        fprintf(fid3, '\n');
        fprintf(fid3, '%s     %f %f', PAR_NAME, average_eddy_movement_c1, average_eddy_movement_c2);
    end
    fclose(fid3);
end

%go back into the current participant's qc directory
cd([derivdir '/groups/' period '/diff_data/' PAR_NAME, '/dwi/eddyqc/']);

%get data from the specific outlier file. Skip the first header line and
%place all values into a matrix.
eddy_outlier_map_values = dlmread('eddy_quad.eddy_outlier_map', ' ', 1, 0);
%only include col 1-72 for the number of slices - there should not be a col
%73.
eddy_outlier_map_values = eddy_outlier_map_values(:,1:72);

%count the number of outliers which are present in each participant's dataset.
NumOutliers = sum(eddy_outlier_map_values(:)==1);

%discount the number of outliers in the B0 shell from the total number of
%outliers
outliers_B0 = 6*72;

%calculate the percentage of outliers in B1000 and B2000 shells (total
%number of outliers)
TotalNumberValues = numel(eddy_outlier_map_values)-outliers_B0;
PercentOutliers = (NumOutliers / TotalNumberValues) * 100;

%go into the main qc folder
cd ([derivdir '/groups/' period '/diff_data/dwiqc/']);

%print the result to the collated file (eddyqc_num_outliers.txt)
fid4 = fopen('eddyqc_outliers.txt', 'a+');
if (fid4 == -1)
    disp('Error in opening the eddyqc_outliers.txt file.')
else
    fprintf(fid4, '\n');
    fprintf(fid4, '%s     %.f %f', PAR_NAME, NumOutliers, PercentOutliers);
end
fclose(fid4);

%get the eddy SNR and CNR values and print the text to file
%go into participant's main dwi eddyqc file
cd([derivdir '/groups/' period '/diff_data/' PAR_NAME, '/dwi/eddyqc/eddy_quad.qc/']);

%read in .json file
fname = 'qc.json';
qcFile = jsondecode(fileread(fname));

%extract the SNR and CNR values from the qc file
SNRvalue = qcFile.qc_cnr_avg(1);
CNR_B1000value = qcFile.qc_cnr_avg(2,1);
CNR_B2000value = qcFile.qc_cnr_avg(3,1);

%go into the main qc folder
cd ([derivdir '/groups/' period '/diff_data/dwiqc/']);

%print values to the SNR/CNR text file
fid5 = fopen('eddyqc_SNR&CNR.txt', 'a+');
if (fid5 == -1)
    disp('Error in opening the eddyqc_SNR&CNR.txt file.')
else
    fprintf(fid5, '\n');
    fprintf(fid5, '%s     %f %f %f', PAR_NAME, SNRvalue, CNR_B1000value, CNR_B2000value);
end
fclose(fid5);


%go back into the current participant's directory
cd([derivdir '/groups/' period '/diff_data/' PAR_NAME, '/dwi/']);

end

