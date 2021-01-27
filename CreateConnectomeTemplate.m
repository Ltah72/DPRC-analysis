function CreateConnectomeTemplate(participants)
%This function will create a Connectome template file for use in order to 
%visualise the comparison across your group nodes. Specifically, the 
%template file will be FreeSurfer's subject average file (fs_average folder). 
%This file will come from FreeSurfer's subject directory, 
%($SUBJECTS_DIR/fsaverage/mri). Essentially, we will be applying the same 
%modifications to the gorup file as the single participant, as done in the 
%main Connectome pipeline script. This function will also call on another 
%function, avg_transmat_T12dwi.m, which will generate an average 
%transformation matrix for the template file. 

%       Inputs(1): participants = list of participants being used for 
%                  analysis          
%       Output(none): You are creating the template file 
%                   (subject_average_hcpmmp1_parcels_coreg.mif).

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 25/01/21


%create a copy of the fsaverage file, called 'subject_average' within FreeSurfer's subject directory
copyfile(['$SUBJECTS_DIR/fsaverage/*'], ['$SUBJECTS_DIR/subject_average']);

%map the annotation file of the HCP MMP 1.0 atlas from fsaverage to your subject. Remember to do that for 
%both hemispheres.
unix(['mri_surf2surf --srcsubject fsaverage --trgsubject subject_average --hemi lh --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.hcpmmp1.annot --tval $SUBJECTS_DIR/subject_average/label/lh.hcpmmp1.annot']);
unix(['mri_surf2surf --srcsubject fsaverage --trgsubject subject_average --hemi rh --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.hcpmmp1.annot --tval $SUBJECTS_DIR/subject_average/label/rh.hcpmmp1.annot']);
 
%map the HCP MMP 1.0 annotations onto the volumetric image and add (FreeSurfer-specific) subcortical segmentation. 
%Convert the resulting file to .mif format (use datatype uint32, which is liked best by MRtrix).
unix(['mri_aparc2aseg --old-ribbon --s subject_average --annot hcpmmp1 --o subject_average_hcpmmp1.mgz']);
unix(['mrconvert -datatype uint32 subject_average_hcpmmp1.mgz subject_average_hcpmmp1.mif']);

%Replace the random integers of the hcpmmp1.mif file with integers that start at 1 and increase by 1.
unix(['labelconvert subject_average_hcpmmp1.mif /SOFTWARE/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_original.txt /SOFTWARE/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt subject_average_hcpmmp1_parcels_nocoreg.mif']);

%get average of transformation matrix (T1 to dwi matrix) in order to do registration in the next step
avg_transmat_T12dwi(participants);

%Register the ordered atlas-based volumetric parcellation to diffusion space.
unix(['mrtransform subject_average_hcpmmp1_parcels_nocoreg.mif -linear average_matrix_transformMrtrixT12Dwi.txt -inverse -datatype uint32 subject_average_hcpmmp1_parcels_coreg.mif']);


end

