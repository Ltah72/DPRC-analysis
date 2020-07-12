function CreateParticipantFixelList()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%for fd list
fid1 = fopen('fd_smooth/files_fd.txt', 'w');
files1 = dir('fd_smooth/*sub*');
fprintf(fid1, '%s\n', files1.name);
fclose(fid1);

%for log_fc list
fid2 = fopen('log_fc_smooth/files_log_fc.txt', 'w');
files2 = dir('log_fc_smooth/*sub*');
fprintf(fid2, '%s\n', files2.name);
fclose(fid2);

%for fdc list
fid3 = fopen('fdc_smooth/files_fdc.txt', 'w');
files3 = dir('fdc_smooth/*sub*');
fprintf(fid3, '%s\n', files3.name);
fclose(fid3);

end

