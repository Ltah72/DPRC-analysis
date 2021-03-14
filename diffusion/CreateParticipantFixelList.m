function CreateParticipantFixelList
%Create a list of all of the participants which will be used for the 
%fixel-based analysis. For the metrics of FD, FC, and FDC.   


%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 15/07/2020

%for fd list
fid1 = fopen('files_fd.txt', 'w');
files1 = dir('template/fd_smooth/*sub*');
fprintf(fid1, '%s\n', files1.name);
fclose(fid1);

%for log_fc list
fid2 = fopen('files_log_fc.txt', 'w');
files2 = dir('template/log_fc_smooth/*sub*');
fprintf(fid2, '%s\n', files2.name);
fclose(fid2);

%for fdc list
fid3 = fopen('files_fdc.txt', 'w');
files3 = dir('template/fdc_smooth/*sub*');
fprintf(fid3, '%s\n', files3.name);
fclose(fid3);

end

