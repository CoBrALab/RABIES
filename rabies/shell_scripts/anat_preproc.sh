#!/bin/bash
#Script to do optimal preprocessing on in-vivo/ex-vivo structural scans
#usage:
#anat_preproc.sh input.nii output.nii

#Operations
# registers to template atlas
# gets approximate brain mask from atlas
# Bias field correction with N4

set -euo pipefail

tmpdir=$(mktemp -d)

input=$1
REGTARGET=$2
REGMASK=$3
output=$4
autoreg=$5

echo "Running anatomical preprocessing"
ImageMath 3 ${tmpdir}/fullmask.nii.gz ThresholdAtMean ${input} 0
ImageMath 3 ${tmpdir}/thresholdmask.nii.gz ThresholdAtMean ${input} 1.2
echo "First iteration of N4BiasFieldCorrection and denoising."
N4BiasFieldCorrection -d 3 -s 4 -i ${input} -b [20] -c [200x200x200,0.0] -w ${tmpdir}/thresholdmask.nii.gz -x ${tmpdir}/fullmask.nii.gz -o ${tmpdir}/N4.nii.gz
DenoiseImage -d 3 -i ${tmpdir}/N4.nii.gz -o ${tmpdir}/denoise.nii.gz

if [[ ${autoreg} == "True" ]]; then
  $RABIES/rabies/shell_scripts/antsRegistration_affine.sh ${tmpdir}/denoise.nii.gz ${REGTARGET} ${REGMASK} ${tmpdir}/trans
else
  echo "Rigid registration to template"
  antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 \
    --output ${tmpdir}/trans_output_ \
    --use-histogram-matching 0 \
    --initial-moving-transform [${REGTARGET},${tmpdir}/denoise.nii.gz,1] \
    --transform Rigid[0.1] \
      --metric Mattes[${REGTARGET},${tmpdir}/denoise.nii.gz,1,64,None] \
      --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 8x6x4 --smoothing-sigmas 4x3x2 --masks [NULL,NULL] \
    --transform Similarity[0.1] \
      --metric Mattes[${REGTARGET},${tmpdir}/denoise.nii.gz,1,64,None] \
      --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 6x4x2 --smoothing-sigmas 3x2x1 --masks [NULL,NULL] \
    --transform Affine[0.1] \
      --metric Mattes[${REGTARGET},${tmpdir}/denoise.nii.gz,1,64,None] \
      --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 6x4x2 --smoothing-sigmas 3x2x1 --masks [NULL,NULL] \
    --transform Affine[0.1] \
      --metric Mattes[${REGTARGET},${tmpdir}/denoise.nii.gz,1,64,None] \
      --convergence [2000x2000x50,1e-6,10,1] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0 --masks [${REGMASK},NULL]
fi


antsApplyTransforms -d 3 -i ${REGMASK} -o ${tmpdir}/mask.nii.gz -t [${tmpdir}/trans_output_0GenericAffine.mat,1] -r ${input} -n GenericLabel

echo "Second iteration of N4BiasFieldCorrection and denoising."
N4BiasFieldCorrection -d 3 -s 2 -i ${input} -b [20] -c [200x200x200x200,0.0] -w ${tmpdir}/mask.nii.gz -r 1 -x ${tmpdir}/fullmask.nii.gz -o ${tmpdir}/N4.nii.gz
DenoiseImage -d 3 -i ${tmpdir}/N4.nii.gz -o ${output}

cp ${tmpdir}/mask.nii.gz $(dirname ${output})/$(basename ${output} .nii.gz)_mask.nii.gz

rm -rf ${tmpdir}
