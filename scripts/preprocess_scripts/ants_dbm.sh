#!/bin/bash

FILE_PATH=$1
template_anat=$2
cluster_type=$3
walltime=$4
memory_request=$5
local_threads=$6

fixed_minimum_resolution=$(python -c "print(min([abs(x) for x in [float(x) for x in \"$(PrintHeader ${template_anat} 1)\".split(\"x\")]]))")
fixed_maximum_resolution=$(python -c "print(max([ a*b for a,b in zip([abs(x) for x in [float(x) for x in \"$(PrintHeader ${template_anat} 1)\".split(\"x\")]],[abs(x) for x in [float(x) for x in \"$(PrintHeader ${template_anat} 2)\".split(\"x\")]])]))")

#set a minimal number of slices for evaluating iteration parameters
ratio=$(python -c "print(int(${fixed_maximum_resolution} / ${fixed_minimum_resolution}))")
if (( 190 > $ratio )); then
  steps=$(ants_generate_iterations.py --min 0.1 --max 19.0 --output twolevel_dbm)
else
  steps=$(ants_generate_iterations.py --min ${fixed_minimum_resolution} --max ${fixed_maximum_resolution} --output twolevel_dbm)
fi


echo twolevel_dbm.py --rigid-model-target $template_anat \
  --no-N4 --skip-dbm --transform SyN --float --average-type normmean --gradient-step 0.25 --model-iterations 3 \
  --modelbuild-command rabies_antsMultivariateTemplateConstruction2.sh --cluster-type $cluster_type \
  --walltime $walltime --memory-request $memory_request --local-threads $local_threads \
  $(eval echo ${steps}) \
  1level $FILE_PATH > exec.sh

echo "Running the following commonspace registration:"
cat exec.sh
bash exec.sh
