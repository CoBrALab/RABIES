#!/bin/bash
#Script to do optimal preprocessing on in-vivo/ex-vivo structural scans
#Taken using the CIC Bruker 7T
#usage:
#mouse-preprocessing-v3.sh input.mnc output.mnc

#Operations
# registers to DSURQE atlas
# gets approximate brain mask from atlas
# Bias field correction with N4

set -euo pipefail

REGTARGET=${QUARANTINE_PATH}/resources/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/100um/DSURQE.mnc
REGMASK=${QUARANTINE_PATH}/resources/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/100um/DSURQE_mask.mnc

tmpdir=$(mktemp -d)

input=$1
output=$2

mincnorm -cutoff 0 -int -out_ceil 65535 -short ${input} ${tmpdir}/renorm.mnc

ThresholdImage 3 ${tmpdir}/renorm.mnc ${tmpdir}/thresholdmask.mnc Otsu 1
minccalc -unsigned -byte -expression '1' ${tmpdir}/renorm.mnc ${tmpdir}/fullmask.mnc

N4BiasFieldCorrection -d 3 -s 4 -i ${tmpdir}/renorm.mnc -b [20] -c [200x200x200,0.0] -w ${tmpdir}/thresholdmask.mnc -x ${tmpdir}/fullmask.mnc -o ${tmpdir}/N4.mnc --verbose
minc_anlm --mt $(nproc) ${tmpdir}/N4.mnc ${tmpdir}/denoise.mnc


antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --minc --verbose \
  --output ${tmpdir}/trans \
  --use-histogram-matching 0 \
  --initial-moving-transform [${REGTARGET},${tmpdir}/denoise.mnc,1] \
  --transform Rigid[0.1] \
    --metric Mattes[${REGTARGET},${tmpdir}/denoise.mnc,1,64,None] \
    --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 8x6x4 --smoothing-sigmas 4x3x2 --masks [NULL,NULL] \
  --transform Similarity[0.1] \
    --metric Mattes[${REGTARGET},${tmpdir}/denoise.mnc,1,64,None] \
    --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 6x4x2 --smoothing-sigmas 3x2x1 --masks [NULL,NULL] \
  --transform Affine[0.1] \
    --metric Mattes[${REGTARGET},${tmpdir}/denoise.mnc,1,64,None] \
    --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 6x4x2 --smoothing-sigmas 3x2x1 --masks [NULL,NULL] \
  --transform Affine[0.1] \
    --metric Mattes[${REGTARGET},${tmpdir}/denoise.mnc,1,64,None] \
    --convergence [2000x2000x50,1e-6,10,1] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0 --masks [${REGMASK},NULL]

antsApplyTransforms -d 3 -i ${REGMASK} -o ${tmpdir}/mask.mnc -t [${tmpdir}/trans0_GenericAffine.xfm,1] -r ${tmpdir}/renorm.mnc -n GenericLabel --verbose

N4BiasFieldCorrection -d 3 -s 2 -i ${tmpdir}/renorm.mnc -b [20] -c [200x200x200x200,0.0] -w ${tmpdir}/mask.mnc -r 1 -x ${tmpdir}/fullmask.mnc -o ${tmpdir}/N4.mnc --verbose
minc_anlm --mt $(nproc) ${tmpdir}/N4.mnc ${output}


cp ${tmpdir}/mask.mnc $(dirname ${output})/$(basename ${output} .mnc)_mask.mnc

rm -rf ${tmpdir}
