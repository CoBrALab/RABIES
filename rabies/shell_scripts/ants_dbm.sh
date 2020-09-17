#!bin/bash

FILE_PATH=$1
template_anat=$2
cluster_type=$3
walltime=$4
memory_request=$5
local_threads=$6

echo twolevel_dbm.py --rigid-model-target $template_anat \
--no-N4 --transform SyN --float --average-type normmean --gradient-step 0.25 --model-iterations 3 \
--modelbuild-command $RABIES/rabies/shell_scripts/antsMultivariateTemplateConstruction2.sh --cluster-type $cluster_type \
--walltime $walltime --memory-request $memory_request --local-threads $local_threads \
1level $FILE_PATH > exec.sh
echo "Running the following commonspace registration:"
cat exec.sh
bash exec.sh
