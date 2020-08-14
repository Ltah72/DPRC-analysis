function CreateTOIFBAMetricFiles(participants, TOI, startdir, groupname)
%This function will create text files which will have all participants
%tract of interest (TOI) FBA metrics (FD, log FC, and FDC). Specifically, 
%the mean, median, std, std_rv (standard error), min, max, and count will 
%be displayed per each metric, per each participant. This will be used for
%group comparisons.

%   Inputs(4) - participants: this will give the number of participants to 
%               calcaluate the metrics per each of them. 
%               TOI: your named tract of interest (TOI)
%               startdir: defined start directory by user
%               groupname: used as part of the directory for where the data
%               is stored. 


%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 03/08/2020


%create headers for each text file
%for FD
FDstr = (['FD_' TOI '_TOI.txt']);
fid9 = fopen(FDstr, 'w');
if (fid9 == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid9, '%s       %s          %s      %s  %s %s   %s', 'mn', 'md', 'std', 'std_rv', 'min', 'max', 'count');
    fprintf(fid9, '\n');
    fclose(fid9);
end

%for log FC
FCstr = (['FC_log_' TOI '_TOI.txt']);
fid10 = fopen(FCstr, 'w');
if (fid10 == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid10, '%s         %s           %s        %s   %s      %s     %s', 'mn', 'md', 'std', 'std_rv', 'min', 'max', 'count');
    fprintf(fid10, '\n');
    fclose(fid10);
end

%for FDC
FDCstr = (['FDC_' TOI '_TOI.txt']);
fid11 = fopen(FDCstr, 'w');
if (fid11 == -1)
    disp('Error in creating the text file.')
else
    fprintf(fid11, '%s       %s          %s      %s  %s %s   %s', 'mn', 'md', 'std', 'std_rv', 'min', 'max', 'count');
    fprintf(fid11, '\n');
    fclose(fid11);
end

%print results onto the text files
for i = 1:length(participants)
    
    full_name = participants(i).name;
    PAR_NAME = full_name(1:15);
    
    %for FD
    fid9 = fopen(FDstr, 'a+');
    %fid9 = fopen('FD_TOI.txt', 'a+');
    if (fid9 == -1)
        disp('Error in opening the text file.')
    else
        unix(['mrstats -mask output_' TOI '_TOI_fixel_directory/' TOI '_fixel_mask_thr.mif ' startdir '/derivatives/diff_data/' groupname '/fd_smooth/' PAR_NAME '_fd.mif -output mean -output median -output std -output std_rv -output min -output max -output count >> FD_' TOI '_TOI.txt']);
        fclose(fid9);
    end
    
    %for FC_log
    fid10 = fopen(FCstr, 'a+');
    if (fid10 == -1)
        disp('Error in opening the text file.')
    else
        unix(['mrstats -mask output_' TOI '_TOI_fixel_directory/' TOI '_fixel_mask_thr.mif ' startdir '/derivatives/diff_data/' groupname '/log_fc_smooth/' PAR_NAME '_log_fc.mif -output mean -output median -output std -output std_rv -output min -output max -output count >> FC_log_' TOI '_TOI.txt']);
        fclose(fid10);
    end
    
    %for FDC
    fid11 = fopen(FDCstr, 'a+');
    if (fid11 == -1)
        disp('Error in opening the text file.')
    else
        unix(['mrstats -mask output_' TOI '_TOI_fixel_directory/' TOI '_fixel_mask_thr.mif ' startdir '/derivatives/diff_data/' groupname '/fdc_smooth/' PAR_NAME '_fdc.mif -output mean -output median -output std -output std_rv -output min -output max -output count >> FDC_' TOI '_TOI.txt']);
        fclose(fid11);
    end
    
end
end


