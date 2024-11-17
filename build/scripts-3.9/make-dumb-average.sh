#!/bin/bash

quick_com_align() {
  tmpdir2=$(mktemp -d)

  moving=$1
  fixed=$2
  output=$3

  antsRegistration --dimensionality 3 --verbose --minc --float 0 --output [ ${tmpdir2}/COM_ ] --use-histogram-matching 0 \
   --initial-moving-transform [${fixed},${moving},1] --transform Translation[0.5] \
   --metric Mattes[${fixed},${moving},1,32,None] -c 0 -s 0 -f 1

  antsApplyTransforms -d 3 -r ${fixed} -i ${moving} -t ${tmpdir2}/COM_0_GenericAffine.xfm -o ${output} --verbose

  rm -rf ${tmpdir2}

}


tmpdir=$(mktemp -d)

template=$1
average_type=$2
trim_percent=$3
shift
shift
shift
inputs="$@"

# get path to script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

AverageImages 3 ${template} 1 ${inputs}

for file in "$@"; do
  quick_com_align ${file} ${template} ${tmpdir}/$(basename ${file})
done

${SCRIPT_DIR}/modelbuild_averager.py -o ${template} --normalize --image_type image --method ${average_type} --trim_percent ${trim_percent} --file_list ${tmpdir}/*

for file in "$@"; do
  quick_com_align ${file} ${template} ${tmpdir}/$(basename ${file})
done

${SCRIPT_DIR}/modelbuild_averager.py -o ${template} --normalize --image_type image --method ${average_type} --file_list ${tmpdir}/*

rm -rf ${tmpdir}
