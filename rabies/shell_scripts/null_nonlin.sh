#!/bin/bash

#Registration method based on the defaults of the antsRegistrationSyN.sh script from the main distro

if [[ -n ${__mb_debug:-} ]]; then

  set -x

fi

set -euo pipefail


moving=$1
fixed=$2
mask=$3
filename_template=$4

antsRegistration --dimensionality 3 \
  --output [${filename_template}_output_,${filename_template}_output_warped_image.nii.gz] \
  --initial-moving-transform [$fixed,$moving,2] \
  --transform Rigid[0.1] --metric Mattes[$fixed,$moving,1,128,None] --convergence [0,1e-6,10] --shrink-factors 1 --smoothing-sigmas 1vox \
  --transform Affine[0.1] --metric Mattes[$fixed,$moving,1,128,None] --convergence [0,1e-6,10] --shrink-factors 1 --smoothing-sigmas 0vox --masks [$mask] \
  --transform SyN[0.2,2,0] --metric CC[$fixed,$moving,1,4] --convergence [0,1e-6,10] \
    --shrink-factors 1 \
    --smoothing-sigmas 1 \
    -z 1
