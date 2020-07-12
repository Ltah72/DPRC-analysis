function GradientEdit_forBestB0(inputBU,PAR_NAME, datafile)
%Edit the gradient files (.bval and .bvec files) from the given BestB0 pair
%that has been chosen. This function is called by the BestB0.m function.
%Not the most neccessary thing, since we are assuming that all B0 (BUs)
%contain the same bval and bvec values. However, other f/mri scan machines
%may not necessarily follow this assumption, Jesper Andersson suggested to
%do this on the fsl forum here:
%https://www.jiscmail.ac.uk/cgi-bin/wa-jisc.exe?A2=ind2007&L=FSL&O=D&X=65845096BD798100FD&Y=ltah262%40aucklanduni.ac.nz&P=94911

% inputs (3) = BU number used for the best B0 pair for topup, participant
%              ID, and datafile format
% outputs (none) = we are editing the text files

%Author: Lenore Tahara-Eckl
%Date: 09/07/2020



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

%edit bval file
fid3 = fopen([PAR_NAME, datafile, '.bval'], 'w');
shell = data(1:106).';
Bval0 = data(1,1);

if inputBU == '1'
    %bval vol 21 should be switched with bval vol 0
    Bval21 = shell(1,22);
    shell(1,1) = Bval21;
    fprintf(fid3, '%d ', shell);
    
elseif inputBU == '2'
    %bval vol 42 should be switched with bval vol 0
    Bval42 = shell(1,43);
    shell(1,1) = Bval42;
    fprintf(fid3, '%d ', shell);
    
elseif inputBU == '3'
    %bval vol 63 should be switched with bval vol 0
    Bval63 = shell(1,64);
    shell(1,1) = Bval63;
    fprintf(fid3, '%d ', shell);
    
    
elseif inputBU == '4'
    %bval vol 84 should be switched with bval vol 0
    Bval84 = shell(1,85);
    shell(1,1) = Bval84;
    fprintf(fid3, '%d ', shell);
    
    
elseif inputBU == '5'
    %bval vol 105 should be switched with bval vol 0
    Bval105 = shell(1,106);
    shell(1,1) = Bval105;
    fprintf(fid3, '%d ', shell);
    
    
end

fclose(fid3);



%edit bvec file
fid4 = fopen([PAR_NAME, datafile,'.bvec'], 'w');
x = data2(1:106).';
y = data2(107:212).';
z = data2(213:318).';
Bvec0x = x(1,1);
Bvec0y = y(1,1);
Bvec0z = z(1,1);

if inputBU == '1'
    %bvec vol 22 should be switched with bvec vol 0
    Bvec21x = x(1,22);
    Bvec21y = y(1,22);
    Bvec21z = z(1,22);
    x(1,1) = Bvec21x;
    x(1,22) = Bvec0x;
    y(1,1) = Bvec21y;
    y(1,22) = Bvec0y;
    z(1,1) = Bvec21z;
    z(1,22) = Bvec0z;
    
elseif inputBU == '2'
    %bvec vol 42 should be switched with bvec vol 0
    Bvec42x = x(1,43);
    Bvec42y = y(1,43);
    Bvec42z = z(1,43);
    x(1,1) = Bvec42x;
    x(1,43) = Bvec0x;
    y(1,1) = Bvec42y;
    y(1,43) = Bvec0y;
    z(1,1) = Bvec42z;
    z(1,43) = Bvec0z;
    
    
elseif inputBU == '3'
    %bvec vol 63 should be switched with bvec vol 0
    Bvec63x = x(1,64);
    Bvec63y = y(1,64);
    Bvec63z = z(1,64);
    x(1,1) = Bvec63x;
    x(1,64) = Bvec0x;
    y(1,1) = Bvec63y;
    y(1,64) = Bvec0y;
    z(1,1) = Bvec63z;
    z(1,64) = Bvec0z;
    
    
elseif inputBU == '4'
    %bvec vol 84 should be switched with bvec vol 0
    Bvec84x = x(1,85);
    Bvec84y = y(1,85);
    Bvec84z = z(1,85);
    x(1,1) = Bvec84x;
    x(1,85) = Bvec0x;
    y(1,1) = Bvec84y;
    y(1,85) = Bvec0y;
    z(1,1) = Bvec84z;
    z(1,85) = Bvec0z;
    
elseif inputBU == '5'
    %bvec vol 105 should be switched with bval vol 0
    Bvec105x = x(1,106);
    Bvec105y = y(1,106);
    Bvec105z = z(1,106);
    x(1,1) = Bvec105x;
    x(1,106) = Bvec0x;
    y(1,1) = Bvec105y;
    y(1,106) = Bvec0y;
    z(1,1) = Bvec105z;
    z(1,106) = Bvec0z;
    
end

%print out the new gradient files the switched B0s. 
fprintf(fid4, '%f ', x);
fprintf(fid4, '\n');
fprintf(fid4, '%f ', y);
fprintf(fid4, '\n');
fprintf(fid4, '%f ', z);

fclose(fid4);



end

