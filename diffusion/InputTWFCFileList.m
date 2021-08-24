function InputTWFCFileList(TWFC_participants)
%Create a text file listing all of the participant track-weighted dynamic 
%functional connectivity file names. We will be looking the 3D files which
%display the confidence interval at 95%. 
%We will use this as an input for conducting the statistics on the
%whole-brain dynamic functional connnectivity (via MRtrix's mrclusterstats). 


%Inputs(1): participants = list of participants being used for analysis
           
%Output(none): You are writing values to the file.

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 17/08/21


fid = fopen('TWFCInput.txt', 'a+');
for i = 1:length(TWFC_participants)
    [~, PAR_NAME, ~] = fileparts(TWFC_participants{1,i});
    if (fid == -1)
        disp('Error in opening the file.')
    else
        %print the file name to a textfile for each participant. 
         fprintf(fid, '%s%s', 'tdfc_CI95_', PAR_NAME, '.mif');
         fprintf(fid, '\n');
    end
    
end
    
fclose(fid);



end
