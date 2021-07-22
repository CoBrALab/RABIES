#!/bin/bash

# Script that doesn't conduct any registration and only provide null transform files
# Allows to keep consistent input/output nipype workflow, even if a step doesn't requiring registration in a particular case

if [[ -n ${__mb_debug:-} ]]; then

  set -x

fi

set -euo pipefail


moving=$1
movingmask=$2
fixed=$3
fixedmask=$4
filename_template=$5

antsRegistration --dimensionality 3 \
  --output [${filename_template}_output_,${filename_template}_output_warped_image.nii.gz] \
  --transform Rigid[0.1] --metric Mattes[$fixed,$moving,1,128,None] --convergence [0,1e-6,10] --shrink-factors 1 --smoothing-sigmas 1vox \
  --transform Affine[0.1] --metric Mattes[$fixed,$moving,1,128,None] --convergence [0,1e-6,10] --shrink-factors 1 --smoothing-sigmas 0vox --masks [$fixedmask,$movingmask] \
  --transform SyN[0.2,2,0] --metric CC[$fixed,$moving,1,4] --convergence [0,1e-6,10] \
    --shrink-factors 1 \
    --smoothing-sigmas 1 \
    -z 1
