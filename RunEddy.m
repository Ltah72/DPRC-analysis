function RunEddy(BU_used, SlicetoVol_Corr, PAR_NAME, datafile, ScriptDirectory)
%Run topup/eddy correction. This function will account for any gradient
%changes which were applied in the previous steps (i.e. BestB0,
%GradientEdit_forBestB0). Also, Depending on the user's input, this
%function can utilise the GPU for eddy cuda to run the slice-to-volume
%motion correction.

%Inputs(4): 
%           BU_used = BU number used in the main dwi dataset sequence
%           SlicetoVol_Corr = user indicated whether they want to use the
%           slice-to-vol motion correction (mp_order)
%           PAR_NAME = participant ID
%           datafile = dwi datafile
%           ScriptDirectory = directory where the processing scripts are
%           contained. 

%Outputs (none): applying eddy and editing the diffusion file.

%Reference:
%Andersson, J. L. R., Graham, M. S., Drobnjak, I., Zhang, H., Filippini, 
%N., & Bastiani, M. (2017). Towards a comprehensive framework for movement 
%and distortion correction of diffusion MR images: Within volume movement. 

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 20/10/2020



if BU_used ~= 1 && SlicetoVol_Corr == 'y'
    %run eddy with gradient edit if you've switched the first BU, and the run mp_order correction
    unix(['dwifslpreproc -fslgrad ', PAR_NAME, datafile, '.bvec ' PAR_NAME, datafile, '.bval bbcgd' PAR_NAME, datafile, '.mif ebbcgd' PAR_NAME, datafile, '.mif -rpe_pair -pe_dir AP -se_epi TUB0s_' PAR_NAME, datafile, '.mif -eddy_mask brain_mask_' PAR_NAME, datafile, '.mif -eddy_options " --repol --ol_nstd=3 --ol_type=both --mporder=6 --s2v_niter=5 --cnr_maps --residuals" -eddy_slspec=' ScriptDirectory '/files/DPRC_slspec.txt -eddyqc_all eddyqc -readout_time 0.07']);
elseif BU_used ~= 1
    %run eddy with gradient edit if you've switched the first BU
    unix(['dwifslpreproc -fslgrad ', PAR_NAME, datafile, '.bvec ' PAR_NAME, datafile, '.bval bbcgd' PAR_NAME, datafile, '.mif ebbcgd' PAR_NAME, datafile, '.mif -rpe_pair -pe_dir AP -se_epi TUB0s_' PAR_NAME, datafile, '.mif -eddy_mask brain_mask_' PAR_NAME, datafile, '.mif -eddy_options " --repol --ol_nstd=3 --ol_type=both --mb=3 --cnr_maps --residuals" -eddyqc_all eddyqc -readout_time 0.07']);
elseif SlicetoVol_Corr == 'y'
    %run just the mp_order correction
    unix(['dwifslpreproc bbcgd' PAR_NAME, datafile, '.mif ebbcgd' PAR_NAME, datafile, '.mif -rpe_pair -pe_dir AP -se_epi TUB0s_' PAR_NAME, datafile, '.mif -eddy_mask brain_mask_' PAR_NAME, datafile, '.mif -eddy_options " --repol --ol_nstd=3 --ol_type=both --mporder=6 --s2v_niter=5 --cnr_maps --residuals" -eddy_slspec=' ScriptDirectory '/files/DPRC_slspec.txt -eddyqc_all eddyqc -readout_time 0.07']);
else
    %don't need to apply gradient edit with eddy, if you have not switched the first BU, and don't apply mp_order correction
    unix(['dwifslpreproc bbcgd' PAR_NAME, datafile, '.mif ebbcgd' PAR_NAME, datafile, '.mif -rpe_pair -pe_dir AP -se_epi TUB0s_' PAR_NAME, datafile, '.mif -eddy_mask brain_mask_' PAR_NAME, datafile, '.mif -eddy_options " --repol --ol_nstd=3 --ol_type=both --mb=3 --cnr_maps --residuals" -eddyqc_all eddyqc -readout_time 0.07']);
end




end

