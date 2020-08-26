function RunEddySquad(subjects, startdir)
%This function will conduct 'eddy_squad', which will combine all
%participants' eddy qc data as a group study. Here, we will be able to view
%which participants seem like outliers. This will create a directory called
%'squad' located in your ([startdir /derivatives/diff_data/dwiqc/squad])
%folder. You can view the pdf file (group_qc.pdf) for the group summary 
%report, and the JASON file (group_db.json) for specific values. My script 
%already puts most of the values onto a text file (motion and outlier
%data), so if you want to extract more data from this file, go ahead. 
%Remember to reference the original authors (e.g. Bastiani et al., 2019) 
%who created this function!


%Inputs(2):     subjects = total number of subjects from your dataset
%               startdir = start directory that you defined in the script -
%               where the data will be stored.
%Output(none):  eddy_squad will generate a folder called 'squad' in the
%               dwiqc directory. This folder will contain several 
%               measurements from the group report of the given dataset.  


%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 21/08/20


%go into the main dwiqc folder
cd ([startdir '/derivatives/diff_data/dwiqc']);

fid = fopen('EddySQUADList.txt', 'a+');
for i = 1:length(subjects)
    [~, PAR_NAME, ~] = fileparts(subjects{1,i});
    if (fid == -1)
        disp('Error in opening the file.')
    else
        %print the file pathway to a textfile for each participant. 
         fprintf(fid, '%s%s%s%s', startdir, '/derivatives/diff_data/', PAR_NAME, '/dwi/eddyqc/eddy_quad.qc');
         fprintf(fid, '\n');
    end
    
end
    
fclose(fid);

%Run eddy_squad with the generated participant list. 
unix(['eddy_squad EddySQUADList.txt']);

end

