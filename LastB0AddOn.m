function LastB0AddOn(PAR_NAME, datafile)
%We are editing the gradient files in order to account for the last B0
%volume (vol 106) in the dwi series. bvecs and bvals files must have same
%number of diffusion directions as DW-image. Note that most of these B0
%values (both for the bval and bvec) should be around 0.

%           inputs (2) = participant ID, and datafile format
%           outputs (none) = we are editing the text files

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 15/04/2020

fid = fopen([PAR_NAME, datafile,'.bval'], 'r');
fid2 = fopen([PAR_NAME, datafile,'.bvec'], 'r');

if (fid == -1 || fid2 == -1)
    disp('Error in opening in one or both of the gradient files.')
else
    %read in the gradient text files
    data = fscanf(fid, '%d');
    fclose(fid);
    data2 = fscanf(fid2, '%f');
    fclose(fid2);
end

%edit the gradient text files-
%edit bval file
fid3 = fopen([PAR_NAME, datafile, '.bval'], 'w');
shell = data(1:105).';

%bval values should be the same as first B0
lastBval = data(1,1);
shell(1,106) = lastBval;

fprintf(fid3, '%d ', shell);

fclose(fid3);

%edit bvec file
fid4 = fopen([PAR_NAME, datafile,'.bvec'], 'w');
x = data2(1:105).';
y = data2(106:210).';
z = data2(211:315).';

%bvec values should be the same as first B0; note that most values should
%be around 0, anyways. 
lastBvecx = x(1,1); 
lastBvecy = y(1,1); 
lastBvecz = z(1,1); 
x(1,106) = lastBvecx; 
y(1,106) = lastBvecy; 
z(1,106) = lastBvecz; 


% lastBvec = data2(1,1) * -1; %make the default value positive
% x(1,106) = lastBvec * -1; %neg orientation value
% y(1,106) = lastBvec * -1; %neg orientation value
% z(1,106) = lastBvec; %pos orientation value

fprintf(fid4, '%f ', x);
fprintf(fid4, '\n');
fprintf(fid4, '%f ', y);
fprintf(fid4, '\n');
fprintf(fid4, '%f ', z);

fclose(fid4);

end

