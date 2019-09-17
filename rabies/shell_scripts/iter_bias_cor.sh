EPI=$1
anat_file=$2
mask=$3
filename_template=$4
reg_script=$5

ResampleImage 3 $EPI resampled.nii.gz 0.4x0.4x0.4 0 4
ImageMath 3 null_mask.nii.gz ThresholdAtMean resampled.nii.gz 0
ImageMath 3 thresh_mask.nii.gz ThresholdAtMean resampled.nii.gz 2
N4BiasFieldCorrection -d 3 -i resampled.nii.gz -b 20 -s 1 -c [100x100x100x100,1e-6] -w thresh_mask.nii.gz -x null_mask.nii.gz -o corrected.nii.gz -v

#iterative registration and bias correction
ResampleImage 3 corrected.nii.gz resampled100.nii.gz 0.1x0.1x0.1 0 4

bash $reg_script resampled100.nii.gz $anat_file $mask $filename_template

antsApplyTransforms -d 3 -i $mask -t ${filename_template}_output_InverseComposite.h5 -r resampled.nii.gz -o ${filename_template}_resampled_mask.nii.gz --verbose -n GenericLabel

N4BiasFieldCorrection -d 3 -i resampled.nii.gz -b 20 -s 1 -c [100x100x100x100,1e-6] -w ${filename_template}_resampled_mask.nii.gz -x null_mask.nii.gz -o iter_corrected.nii.gz -v

ResampleImage 3 iter_corrected.nii.gz ${filename_template}_bias_cor.nii.gz 0.1x0.1x0.1 0 4
