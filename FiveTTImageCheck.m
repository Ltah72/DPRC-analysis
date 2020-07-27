function FiveTTImageCheck(PAR_NAME)
%This function will check each participant's 5tt image file in that it
%conforms to the expected five-tissue-type (5TT) format. It will print out
%all of the checking of the results onto a text file (5ttImage.txt). 

%   1 input = PAR_NAME: current participant. 

%Author: Lenore Tahara-Eckl
%Date: 25/07/2020


fid5 = fopen('5ttImageCheck.txt','a+');

results = unix(['5ttcheck 5ttimage_' PAR_NAME, '.mif']);

if results == 0
    text2file = ('Input image checked OK');
elseif results == 1
    text2file = ('Error - check 5ttimage');
end

%print out results on the 5ttImageCheck.txt file.
fprintf(fid5, '\n');
fprintf(fid5, '%s %s', PAR_NAME, text2file);

fclose(fid5);

end

