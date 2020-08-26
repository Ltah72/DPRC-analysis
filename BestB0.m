function [BU_used] = BestB0(PAR_NAME, datafile, NumBDs, startdir)

%This function will calculate the 'BestB0' to find the most suitable pair
%of blip-up (BU or AP) and blip-down (BD or PA) volumes used for topup in
%the topup/eddy steps of the preprocessing pipeline. This is done by the
%Biobank. This function also calls upon another function, called
%GradientEdit_forBestB0.m, which will edit the gradient files (.bval and
%.bvec files).

%To do this, first all of the B0 images in both phase encoding directions
%(AP and PA) are aligned with one another with a rigid-body registration
%with 6 dof (using fsl's FLIRT tool). Next, the correlation is calculated
%between each of the b0 images to all of the others (fsl's FSLCC tool). The
%function checks that at least one of the B0 images has a cross-correlation
%greater than 0.95 - if not, then that participant is flagged, but will
%still use the highest (and earliest, if tied) cross-correlation B0 for it.
%If the first B0 has sufficient quality (a correlation of 0.98 or higher -
%Jesper's criterion), we would select this as the 'best B0 image' to use.
%If this is not the case, then the second B0 is checked, and so on, so
%forth. If none of the B0s have a higher correlation than 0.98, but there
%is at least one B0 greater than 0.95, then we would select for the highest
%correlation B0 (and earliest, if tied) for it. A text file is generated,
%which will show the participant ID, the B0 status, and the BU and BD
%number used in the sequence. In addition, the best BU volume, which is
%selected, will now be the first volume in the dwi sequence - the first
%volume and selected volume will switch places (if this is done).

%Inputs(4): PAR_NAME = participant ID
%           datafile = dwi datafile
%           NumBDs = number of BD files participant has (typically 3)
%           startdir = start directory that you defined in the script -
%           where the data will be stored.

%Outputs (1): BU number used in the main dwi dataset sequence

%Reference:
%Alfaro-Almagro, F., Jenkinson, M., Bangerter, N. K., Andersson, J. L. R.,
%Griffanti, L., Douaud, G., … Smith, S. M. (2018). Image processing and
%Quality Control for the first 10,000 brain imaging datasets from UK
%Biobank.


%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 01/07/2020

%convert AP and PA from .mif to .nii
unix(['mrconvert AP_' PAR_NAME, datafile,'.mif AP_' PAR_NAME, datafile,'.nii']);
unix(['mrconvert PA_' PAR_NAME, datafile,'.mif PA_' PAR_NAME, datafile,'.nii']);

%put into separate vols to prepare for registration

%for APs or BUs
unix(['mrconvert AP_' PAR_NAME, datafile '.nii BU1.nii -coord 3 0']);
unix(['mrconvert AP_' PAR_NAME, datafile '.nii BU2.nii -coord 3 1']);
unix(['mrconvert AP_' PAR_NAME, datafile '.nii BU3.nii -coord 3 2']);
unix(['mrconvert AP_' PAR_NAME, datafile '.nii BU4.nii -coord 3 3']);
unix(['mrconvert AP_' PAR_NAME, datafile '.nii BU5.nii -coord 3 4']);
unix(['mrconvert AP_' PAR_NAME, datafile '.nii BU6.nii -coord 3 5']);

%for PAs or BDs
if NumBDs == 1
    copyfile (['PA_' PAR_NAME, datafile '.nii'], ['BD1.nii']);
elseif NumBDs == 2
    unix(['mrconvert PA_' PAR_NAME, datafile '.nii BD1.nii -coord 3 0']);
    unix(['mrconvert PA_' PAR_NAME, datafile '.nii BD2.nii -coord 3 1']);
elseif NumBDs == 3
    unix(['mrconvert PA_' PAR_NAME, datafile '.nii BD1.nii -coord 3 0']);
    unix(['mrconvert PA_' PAR_NAME, datafile '.nii BD2.nii -coord 3 1']);
    unix(['mrconvert PA_' PAR_NAME, datafile '.nii BD3.nii -coord 3 2']);
end

%create the registered image for the BU
unix(['flirt -in BU2.nii -ref BU1.nii -dof 6 -out BU_Registered_' PAR_NAME, '.nii']);
unix(['flirt -in BU3.nii -ref BU_Registered_' PAR_NAME, '.nii -dof 6 -out BU_Registered_' PAR_NAME, '.nii']);
unix(['flirt -in BU4.nii -ref BU_Registered_' PAR_NAME, '.nii -dof 6 -out BU_Registered_' PAR_NAME, '.nii']);
unix(['flirt -in BU5.nii -ref BU_Registered_' PAR_NAME, '.nii -dof 6 -out BU_Registered_' PAR_NAME, '.nii']);
unix(['flirt -in BU6.nii -ref BU_Registered_' PAR_NAME, '.nii -dof 6 -out BU_Registered_' PAR_NAME, '.nii']);

%create the registered image for the BD
%note that if BD = 1, then just use the 1 BD image as the registered image.
if NumBDs == 1
    copyfile (['BD1.nii'], ['BD_Registered_' PAR_NAME '.nii']);
elseif NumBDs == 2
    unix(['flirt -in BD2.nii -ref BD1.nii -dof 6 -out BD_Registered_' PAR_NAME, '.nii']);
elseif NumBDs == 3
    unix(['flirt -in BD2.nii -ref BD1.nii -dof 6 -out BD_Registered_' PAR_NAME, '.nii']);
    unix(['flirt -in BD3.nii -ref BD_Registered_' PAR_NAME, '.nii -dof 6 -out BD_Registered_' PAR_NAME, '.nii']);
end

%-------------------------------------------------------------------------%

%run cross-correlations between every volume in one 4D data set with every
%volume in another (for investigating similarities in ICA outputs), and put
%into a text file.

%for BUs
unix(['fslcc BU_Registered_' PAR_NAME '.nii.gz AP_' PAR_NAME, datafile '.nii > BUScores.txt'])

%for BDs
if (NumBDs == 2 || NumBDs == 3)
    unix(['fslcc BD_Registered_' PAR_NAME '.nii.gz PA_' PAR_NAME, datafile '.nii > BDScores.txt']);
end

%read scores and put into a vector
fid = fopen('BUScores.txt', 'r');
if (fid == -1)
    disp('Error in opening the BU correlation score file (BUScores.txt).')
else
    %read in the scores from the text files
    BUCorrScore = fscanf(fid, '%f');
    %reorganise to show B0 number and corr score
    BUCorrScore = BUCorrScore([2,3; 5,6; 8,9; 11,12; 14,15; 17,18]);
    fclose(fid);
end

if (NumBDs == 2 || NumBDs == 3)
    fid2 = fopen('BDScores.txt', 'r');
    if (fid2 == -1)
        disp('Error in opening in the BD correlation score file (BDScores.txt).')
    else
        BDCorrScore = fscanf(fid2, '%f');
        
        if NumBDs == 2
        %reorganise to show B0 number and corr score
        BDCorrScore = BDCorrScore([2,3; 5,6]);
        fclose(fid2);
        elseif NumBDs == 3
        BDCorrScore = BDCorrScore([2,3; 5,6; 8,9]);
        fclose(fid2);
        end 
    end
end


%-------------------------------------------------------------------------%


%create a text file to note down the BU and BD used for the B0, and if any
%participants have been flagged. Also, choose the best B0 according to
%Jesper's criterion (see info above and/or in paper referenced)
cd([startdir '/derivatives/diff_data/dwiqc/']);

fid3 = fopen('BestB0.txt', 'a+');
if (fid3 == -1)
    disp('Error in opening in the BestB0.txt file.')
else
    cd([startdir '/derivatives/diff_data/' PAR_NAME, '/dwi/']);
end

%check that at least 1 B0 is > 0.95; if not, then flag the participant.
if NumBDs == 1
    
    if ((BUCorrScore(1,2) > 0.95) || (BUCorrScore(2,2) > 0.95) || (BUCorrScore(3,2) > 0.95) || (BUCorrScore(4,2) > 0.95) || (BUCorrScore(5,2) > 0.95) || (BUCorrScore(6,2) > 0.95)) || (BDCorrScore(1,2) > 0.95)
        B0_status = 'OK!'; 
    else
        B0_status = 'all < 0.95';
    end
    
    %use BU if greater than 0.98, in sequential order
    if (BUCorrScore(1,2) >= 0.98)
        BU_used = 1;
    elseif (BUCorrScore(2,2) >= 0.98)
        BU_used = 2;
    elseif (BUCorrScore(3,2) >= 0.98)
        BU_used = 3;
    elseif (BUCorrScore(4,2) >= 0.98)
        BU_used = 4;
    elseif (BUCorrScore(5,2) >= 0.98)
        BU_used = 5;
    elseif (BUCorrScore(6,2) >= 0.98)
        BU_used = 6;
    else
        %if no BU is greater than 0.98, then use the largest (and earliest,
        %tied) available BU.
        [BU_score, BU_used] = max(BUCorrScore(:,2));
    end
    
    BD_used = 1;
    
elseif NumBDs == 2
    
    if ((BUCorrScore(1,2) > 0.95) || (BUCorrScore(2,2) > 0.95) || (BUCorrScore(3,2) > 0.95) || (BUCorrScore(4,2) > 0.95) || (BUCorrScore(5,2) > 0.95) || (BUCorrScore(6,2) > 0.95)) || (BDCorrScore(1,2) > 0.95) || (BDCorrScore(2,2) > 0.95)
        B0_status = 'OK!';
    else
        B0_status = 'all < 0.95';
    end
    %use BU if greater than 0.98, in sequential order
    if (BUCorrScore(1,2) >= 0.98)
        BU_used = 1;
    elseif (BUCorrScore(2,2) >= 0.98)
        BU_used = 2;
    elseif (BUCorrScore(3,2) >= 0.98)
        BU_used = 3;
    elseif (BUCorrScore(4,2) >= 0.98)
        BU_used = 4;
    elseif (BUCorrScore(5,2) >= 0.98)
        BU_used = 5;
    elseif (BUCorrScore(6,2) >= 0.98)
        BU_used = 6;
    else
        %if no BU is greater than 0.98, then use the largest (and earliest,
        %tied) available BU.
        [BU_score, BU_used] = max(BUCorrScore(:,2)); 
    end
    
    %use BD if greater than 0.98, in sequential order
    if (BDCorrScore(1,2) >= 0.98)
        BD_used = 1;
    elseif (BDCorrScore(2,2) >= 0.98)
        BD_used = 2;
    else
        %if no BD is greater than 0.98, then use the largest (and earliest,
        %tied) available BD.
        [BD_score, BD_used] = max(BDCorrScore(:,2));
    end
    
elseif NumBDs == 3
    
    if ((BUCorrScore(1,2) > 0.95) || (BUCorrScore(2,2) > 0.95) || (BUCorrScore(3,2) > 0.95) || (BUCorrScore(4,2) > 0.95) || (BUCorrScore(5,2) > 0.95) || (BUCorrScore(6,2) > 0.95)) || (BDCorrScore(1,2) > 0.95) || (BDCorrScore(2,2) > 0.95) || (BDCorrScore(3,2) > 0.95)
        B0_status = 'OK!';
    else 
        B0_status = 'all < 0.95';
    end
    %use BU if greater than 0.98, in sequential order
    if (BUCorrScore(1,2) >= 0.98)
        BU_used = 1;
    elseif (BUCorrScore(2,2) >= 0.98)
        BU_used = 2;
    elseif (BUCorrScore(3,2) >= 0.98)
        BU_used = 3;
    elseif (BUCorrScore(4,2) >= 0.98)
        BU_used = 4;
    elseif (BUCorrScore(5,2) >= 0.98)
        BU_used = 5;
    elseif (BUCorrScore(6,2) >= 0.98)
        BU_used = 6;
    else
        %if no BU is greater than 0.98, then use the largest (and earliest,
        %tied) available BU.
        [BU_score, BU_used] = max(BUCorrScore(:,2));
    end
    %use BD if greater than 0.98, in sequential order
    if (BDCorrScore(1,2) >= 0.98)
        BD_used = 1;
    elseif (BDCorrScore(2,2) >= 0.98)
        BD_used = 2;
    elseif (BDCorrScore(3,2) >= 0.98)
        BD_used = 3;
    else
        %if no BD is greater than 0.98, then use the largest (and earliest,
        %tied) available BD.
        [BD_score, BD_used] = max(BDCorrScore(:,2));
    end
end

%print out results on the BestB0.txt file.
fprintf(fid3, '\n');
fprintf(fid3, '%s %s       %d       %d', PAR_NAME, B0_status, BU_used, BD_used);

fclose(fid3);

%create B0 pair to use for topup.
inputBU = num2str(BU_used-1);
inputBD = num2str(BD_used+5);

unix(['mrconvert -coord 3 ' inputBU, ',' inputBD, ' allB0s_' PAR_NAME, datafile,'.mif TUB0s_' PAR_NAME, datafile,'.mif']);
unix(['mrconvert TUB0s_' PAR_NAME, datafile '.mif TUB0s_' PAR_NAME, datafile '.nii']);

unix(['mrconvert allB0s_' PAR_NAME, datafile, '.mif allB0s_' PAR_NAME, datafile, '.nii']);
%-------------------------------------------------------------------------%

%Last thing you need to do, is change the main dwi sequence, so that the BU
%image used comes first in the sequence. And switch places with the other
%the other BU in the dwi sequence. You need to edit the gradient files as
%well for this (.bval and .bvec files).

%select the BU used
unix(['mrconvert TUB0s_' PAR_NAME, datafile '.mif newBU1.mif -coord 3 0']);
unix(['mrconvert newBU1.mif newBU1.nii']);

%no changes needed if the BestB0 for the BU was vol 0.
if inputBU == '0'
    copyfile (['bcgd' PAR_NAME, datafile '.mif'], ['bbcgd' PAR_NAME, datafile '.mif']);
elseif inputBU == '1'
    %edit the sequence
    unix(['mrconvert -coord 3 1:20 bcgd' PAR_NAME, datafile,'.mif beforeBestB0' PAR_NAME, datafile,'.mif']);
    unix(['mrcat newBU1.mif beforeBestB0' PAR_NAME, datafile,'.mif bbcgd-pre1' PAR_NAME, datafile,'.mif']);
    unix(['mrcat bbcgd-pre1' PAR_NAME, datafile,'.mif BU1.nii bbcgd-pre2' PAR_NAME, datafile,'.mif']);
    unix(['mrconvert -coord 3 22:105 bcgd' PAR_NAME, datafile,'.mif afterBestB0' PAR_NAME, datafile,'.mif']);
    unix(['mrcat bbcgd-pre2' PAR_NAME, datafile,'.mif afterBestB0' PAR_NAME, datafile,'.mif bbcgd' PAR_NAME, datafile,'.mif']);
    %make a nifti copy
    unix(['mrconvert bbcgd' PAR_NAME, datafile '.mif bbcgd' PAR_NAME, datafile '.nii']);
    
elseif inputBU == '2'
    unix(['mrconvert -coord 3 1:41 bcgd' PAR_NAME, datafile,'.mif beforeBestB0' PAR_NAME, datafile,'.mif']);
    unix(['mrcat newBU1.mif beforeBestB0' PAR_NAME, datafile,'.mif bbcgd-pre1' PAR_NAME, datafile,'.mif']);
    unix(['mrcat bbcgd-pre1' PAR_NAME, datafile,'.mif BU1.nii bbcgd-pre2' PAR_NAME, datafile,'.mif']);
    unix(['mrconvert -coord 3 43:105 bcgd' PAR_NAME, datafile,'.mif afterBestB0' PAR_NAME, datafile,'.mif']);
    unix(['mrcat bbcgd-pre2' PAR_NAME, datafile,'.mif afterBestB0' PAR_NAME, datafile,'.mif bbcgd' PAR_NAME, datafile,'.mif']);
    %make a nifti copy
    unix(['mrconvert bbcgd' PAR_NAME, datafile '.mif bbcgd' PAR_NAME, datafile '.nii']);
    
elseif inputBU == '3'
    unix(['mrconvert -coord 3 1:62 bcgd' PAR_NAME, datafile,'.mif beforeBestB0' PAR_NAME, datafile,'.mif']);
    unix(['mrcat newBU1.mif beforeBestB0' PAR_NAME, datafile,'.mif bbcgd-pre1' PAR_NAME, datafile,'.mif']);
    unix(['mrcat bbcgd-pre1' PAR_NAME, datafile,'.mif BU1.nii bbcgd-pre2' PAR_NAME, datafile,'.mif']);
    unix(['mrconvert -coord 3 64:105 bcgd' PAR_NAME, datafile,'.mif afterBestB0' PAR_NAME, datafile,'.mif']);
    unix(['mrcat bbcgd-pre2' PAR_NAME, datafile,'.mif afterBestB0' PAR_NAME, datafile,'.mif bbcgd' PAR_NAME, datafile,'.mif']);
    %make a nifti copy
    unix(['mrconvert bbcgd' PAR_NAME, datafile '.mif bbcgd' PAR_NAME, datafile '.nii']);
    
elseif inputBU == '4'
    unix(['mrconvert -coord 3 1:83 bcgd' PAR_NAME, datafile,'.mif beforeBestB0' PAR_NAME, datafile,'.mif']);
    unix(['mrcat newBU1.mif beforeBestB0' PAR_NAME, datafile,'.mif bbcgd-pre1' PAR_NAME, datafile,'.mif']);
    unix(['mrcat bbcgd-pre1' PAR_NAME, datafile,'.mif BU1.nii bbcgd-pre2' PAR_NAME, datafile,'.mif']);
    unix(['mrconvert -coord 3 85:105 bcgd' PAR_NAME, datafile,'.mif afterBestB0' PAR_NAME, datafile,'.mif']);
    unix(['mrcat bbcgd-pre2' PAR_NAME, datafile,'.mif afterBestB0' PAR_NAME, datafile,'.mif bbcgd' PAR_NAME, datafile,'.mif']);
    %make a nifti copy
    unix(['mrconvert bbcgd' PAR_NAME, datafile '.mif bbcgd' PAR_NAME, datafile '.nii']);
    
elseif inputBU == '5'
    unix(['mrconvert -coord 3 1:104 bcgd' PAR_NAME, datafile,'.mif beforeBestB0' PAR_NAME, datafile,'.mif']);
    unix(['mrcat newBU1.mif beforeBestB0' PAR_NAME, datafile,'.mif bbcgd-pre1' PAR_NAME, datafile,'.mif']);
    unix(['mrcat bbcgd-pre1' PAR_NAME, datafile,'.mif BU1.nii bbcgd' PAR_NAME, datafile,'.mif']);
    %make a nifti copy
    unix(['mrconvert bbcgd' PAR_NAME, datafile '.mif bbcgd' PAR_NAME, datafile '.nii']);
    
end

%edit the .bval and .bvec file - only do this is if the BU used for topup 
%is not the first one (vol 0) in the diffusion sequence. 
if BU_used ~= 1
    GradientEdit_forBestB0(inputBU, PAR_NAME, datafile);
end

end

