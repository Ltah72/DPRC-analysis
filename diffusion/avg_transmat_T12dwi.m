function avg_transmat_T12dwi(participants)
%The purpose of this function is to take the average transformation matrix
%(anatomical (T1w) to diffusion (dwi)), to then apply it to the subject
%average parcellation template (fsaverage). 

%       Inputs(1): participants = list of participants being used for 
%                  analysis          
%       Output(none): You are writing values to the file.

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 19/01/21

%preallocate matrix
SumMatrix = zeros(4:4);

%take data from each participants' matrix .txt file
for i = 1:length(participants)
    
    [~, PAR_NAME, ~] = fileparts(participants{1,i});
  
    %exclude the first row for reading in line
    %p_matrix  = dlmread('transform_mrtrix_t12dwi.txt', ' ',1,0);

    p_matrix  = dlmread(['transform_mrtrix_t12dwi_' PAR_NAME '.txt'], ' ',1,0);
    
    SumMatrix = SumMatrix + p_matrix;
end

%take average of all participant matrices
average_matrix_transformMrtrixT12Dwi = SumMatrix / length(participants);

%re-fill the last row as 0 0 0 1 (to not be affected by the previous calculations)
average_matrix_transformMrtrixT12Dwi(4,:) = [0,0,0,1];

%print average matrix to text file 
writematrix(average_matrix_transformMrtrixT12Dwi, 'average_matrix_transformMrtrixT12Dwi.txt', 'Delimiter', 'space')

end

