#!/bin/bash

# This is a generic script for the execution of a registration command based on calling ants_generate_iterations.py

set -euo pipefail
set -x

movingfile=$1
fixedfile=$2
movingmask=NULL
fixedmask=$3
_arg_outputbasename=$4
method=$5

fixed_minimum_resolution=$(python -c "print(min([abs(x) for x in [float(x) for x in \"$(PrintHeader ${fixedfile} 1)\".split(\"x\")]]))")
fixed_maximum_resolution=$(python -c "print(max([ a*b for a,b in zip([abs(x) for x in [float(x) for x in \"$(PrintHeader ${fixedfile} 1)\".split(\"x\")]],[abs(x) for x in [float(x) for x in \"$(PrintHeader ${fixedfile} 2)\".split(\"x\")]])]))")

#set a minimal number of slices for evaluating iteration parameters
ratio=$(python -c "print(int(${fixed_maximum_resolution} / ${fixed_minimum_resolution}))")
if (( 190 > $ratio )); then
  steps_rigid=$(ants_generate_iterations.py --min 0.1 --max 19.0 --output rigid)
  steps_affine=$(ants_generate_iterations.py --min 0.1 --max 19.0 --output multilevel-halving)
  steps_syn=$(ants_generate_iterations.py --min 0.1 --max 19.0)
else
  steps_rigid=$(ants_generate_iterations.py --min ${fixed_minimum_resolution} --max ${fixed_maximum_resolution} --output rigid)
  steps_affine=$(ants_generate_iterations.py --min ${fixed_minimum_resolution} --max ${fixed_maximum_resolution} --output multilevel-halving)
  steps_syn=$(ants_generate_iterations.py --min ${fixed_minimum_resolution} --max ${fixed_maximum_resolution})
fi

if [[ $method == "SyN" ]]; then
  antsRegistration --dimensionality 3 --verbose \
    --output [ ${_arg_outputbasename}_output_,${_arg_outputbasename}_output_warped_image.nii.gz ] \
    --use-histogram-matching 1 \
    --initial-moving-transform [ ${fixedfile},${movingfile},1 ] \
    $(eval echo ${steps_affine}) \
    --transform SyN[ 0.1,3,0 ] \
    --metric CC[ ${fixedfile},${movingfile},1,4,None ] \
    $(eval echo ${steps_syn}) \
    --masks [ ${fixedmask},${movingmask} ]
elif [[ $method == "Affine" ]]; then
  antsRegistration --dimensionality 3 --verbose \
    --output [ ${_arg_outputbasename}_output_,${_arg_outputbasename}_output_warped_image.nii.gz ] \
    --use-histogram-matching 1 \
    --initial-moving-transform [ ${fixedfile},${movingfile},1 ] \
    $(eval echo ${steps_affine})
elif [[ $method == "Rigid" ]]; then
  antsRegistration --dimensionality 3 --verbose \
    --output [ ${_arg_outputbasename}_output_,${_arg_outputbasename}_output_warped_image.nii.gz ] \
    --use-histogram-matching 1 \
    --initial-moving-transform [ ${fixedfile},${movingfile},1 ] \
    $(eval echo ${steps_rigid})
else
  echo "Registration method must be 'Rigid','Affine', or 'SyN'."
fi
