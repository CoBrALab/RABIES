EPI=$1
anat_file=$2
mask=$3
filename_template=$4
reg_script=$5

ImageMath 3 null_mask.nii.gz ThresholdAtMean $EPI 0
ImageMath 3 thresh_mask.nii.gz ThresholdAtMean $EPI 2
N4BiasFieldCorrection -d 3 -i $EPI -b 20 -s 1 -c [100x100x100x100,1e-6] -w thresh_mask.nii.gz -x null_mask.nii.gz -o corrected.nii.gz

#iterative registration and bias correction
bash $reg_script corrected.nii.gz $anat_file $mask $filename_template
antsApplyTransforms -d 3 -i $mask -t [${filename_template}_output_0GenericAffine.mat,1] -r $EPI -o ${filename_template}_resampled_mask.nii.gz -n GenericLabel
N4BiasFieldCorrection -d 3 -i $EPI -b 20 -s 1 -c [100x100x100x100,1e-6] -w ${filename_template}_resampled_mask.nii.gz -x null_mask.nii.gz -o iter_corrected.nii.gz
