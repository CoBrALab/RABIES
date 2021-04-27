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
  --initial-moving-transform [$fixed,$moving,1] \
  --transform Rigid[0.1] --metric Mattes[$fixed,$moving,1,128,None] --convergence [2025x2025x2025x2025x675,1e-6,10] --shrink-factors 8x7x6x5x4 --smoothing-sigmas 3.98448927075x3.4822628776x2.97928762436x2.47510701762x1.96879525311vox --masks [NULL,NULL] \
  --transform Affine[0.1] --metric Mattes[$fixed,$moving,1,128,None] --convergence [675x225x200x200x200,1e-6,10] --shrink-factors 4x3x2x1x1 --smoothing-sigmas 1.96879525311x1.45813399545x0.936031382318x0.355182697615x0vox --masks [$mask,NULL] \
  --transform SyN[0.2,2,0] --metric CC[$fixed,$moving,1,4] --convergence [2025x2025x200x50x25x25,1e-6,10] \
    --shrink-factors 6x4x3x2x1x1 \
    --smoothing-sigmas 5.48872979374x3.4822628776x2.47510701762x1.45813399545x0.936031382318x0 \
    --masks [NULL,NULL] -z 1
