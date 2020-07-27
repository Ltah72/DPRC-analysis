function CreateParticipantFixelList
%Create a list of all of the participants which will be used for the 
%fixel-based analysis. For the metrics of FD, FC, and FDC.   

%for fd list
fid1 = fopen('files_fd.txt', 'w');
files1 = dir('fd_smooth/*sub*');
fprintf(fid1, '%s\n', files1.name);
fclose(fid1);

%for log_fc list
fid2 = fopen('files_log_fc.txt', 'w');
files2 = dir('log_fc_smooth/*sub*');
fprintf(fid2, '%s\n', files2.name);
fclose(fid2);

%for fdc list
fid3 = fopen('files_fdc.txt', 'w');
files3 = dir('fdc_smooth/*sub*');
fprintf(fid3, '%s\n', files3.name);
fclose(fid3);

end

