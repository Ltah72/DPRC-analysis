function InputConnectomeFileList(participants)
%Create a text file listing all of the participant connectome file names.
%We will use this as an input for conducting the statistics on the
%whole-brain connectome (via MRtrix's connectomestats). 


%Inputs(1): participants = list of participants being used for analysis
           
%Output(none): You are writing values to the file.

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 27/12/20


fid = fopen('ConnectomeInput.txt', 'a+');
for i = 1:length(participants)
    [~, PAR_NAME, ~] = fileparts(participants{1,i});
    if (fid == -1)
        disp('Error in opening the file.')
    else
        %print the file name to a textfile for each participant. 
         fprintf(fid, '%s%s', 'hcpmmp1_', PAR_NAME, '.csv');
         fprintf(fid, '\n');
    end
    
end
    
fclose(fid);



end

