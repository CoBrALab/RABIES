#!/bin/bash

#Robust registration method adapted by Dr. Joanes Grandjean for the mutlirat study

if [[ -n ${__mb_debug:-} ]]; then

  set -x

fi

set -euo pipefail


moving_pre=$1
movingmask=$2
fixed_pre=$3
fixedmask=$4
filename_template=$5
moving='scan_autobox.nii.gz'
fixed='template_autobox.nii.gz'

3dAutobox -input $moving_pre -prefix $moving
3dAutobox -input $fixed_pre -prefix $fixed


antsRegistration --dimensionality 3 \
  --output [${filename_template}_tmp_,${filename_template}_tmp_warped_image.nii.gz] \
  --initial-moving-transform [$fixed,$moving,1] --winsorize-image-intensities [0.1,0.995] \
  --transform Rigid[0.1] --metric MI[$fixed,$moving,1,32,Regular,0.25] --convergence [1000x1000x1000x500,1e-6,10] --shrink-factors 3x2x1x1 --smoothing-sigmas 3x2x1x0vox -v 0 -z 1


antsRegistration --dimensionality 3 \
  --output [${filename_template}_tmp2_,${filename_template}_tmp2_warped_image.nii.gz] \
  --transform Similarity[0.1] --metric MI[$fixed,${filename_template}_tmp_warped_image.nii.gz,1,32,Regular,0.25] --convergence [1000x1000x1000x500x100,1e-8,10] --shrink-factors 4x3x2x1x1 --smoothing-sigmas 4x2x2x1x0vox \
  --transform SyN[0.1,3,0] --metric CC[$fixed,${filename_template}_tmp_warped_image.nii.gz,1,4] --convergence [20x15x10,1e-6,10] \
  --shrink-factors 3x2x1 \
  --smoothing-sigmas 2x1x0vox -v 0 -z 1


ComposeMultiTransform 3 ${filename_template}_output_1Warp.nii.gz -R ${fixed_pre} ${filename_template}_tmp2_1Warp.nii.gz ${filename_template}_tmp2_0GenericAffine.mat
ComposeMultiTransform 3 ${filename_template}_output_1InverseWarp.nii.gz -R ${fixed_pre} -i ${filename_template}_tmp2_0GenericAffine.mat ${filename_template}_tmp2_1InverseWarp.nii.gz

cp ${filename_template}_tmp_0GenericAffine.mat ${filename_template}_output_0GenericAffine.mat

antsApplyTransforms -i ${moving_pre} -r ${fixed_pre} -t ${filename_template}_output_1Warp.nii.gz -t ${filename_template}_output_0GenericAffine.mat -o ${filename_template}_output_warped_image.nii.gz


rm $moving
rm $fixed
rm ${filename_template}_tmp_*
rm ${filename_template}_tmp2_*
