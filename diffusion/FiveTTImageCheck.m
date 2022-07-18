function FiveTTImageCheck(PAR_NAME, startdir, period, groupname)
%This function will check each participant's 5tt image file in that it
%conforms to the expected five-tissue-type (5TT) format. It will print out
%all of the checking of the results onto a text file (5ttImage.txt). 

%   4 inputs PAR_NAME: current participant. 
%            startdir = start directory that you defined in the script -
%                       where the data will be stored.
%            period = time period of the participant MRI scans
%            groupname = name of the group you are analysing - used as 
%                        part of the directory for where the data is 
%                        stored.

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 25/07/2020


%check 5tt image
results = unix(['5ttcheck 5ttimage_' PAR_NAME, '.mif']);

%go back into the qc folder to write about the file
cd([startdir '/derivatives/' period, '/diff_data/', groupname, '/qc/']);

%open txt file
fid = fopen('5ttImageCheck.txt','a+');

%unix(['5ttcheck 5ttimage_' PAR_NAME, '.mif >> 5ttImageCheck.txt']);

if results == 0
    text2file = ('Input image checked OK');
elseif results == 1
    text2file = ('Error - check 5ttimage');
end

%print out results on the 5ttImageCheck.txt file.
fprintf(fid, '\n');
fprintf(fid, '%s %s', PAR_NAME, text2file);

fclose(fid);

%go back into connectome folder to continue processing
cd([startdir '/derivatives/' period, '/diff_data/', groupname, '/connectome/']);

end

