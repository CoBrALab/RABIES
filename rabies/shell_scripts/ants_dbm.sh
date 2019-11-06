#!bin/bash

FILE_PATH=$1
info_csv=$2
template=$template_anat

mkdir -p ants_dbm
cd ants_dbm

echo twolevel_dbm.py --rigid-model-target $template_anat \
--no-N4 --transform SyN --float --average-type normmean --gradient-step 0.25 --model-iterations 3 \
--modelbuild-command antsMultivariateTemplateConstruction2.sh --cluster-type $ants_dbm_cluster_type \
--walltime $ants_dbm_walltime --memory-request $ants_dbm_memory_request --local-threads $ants_dbm_local_threads \
1level $FILE_PATH > exec.sh
echo "Running the following commonspace registration:"
cat exec.sh
bash exec.sh

cd ..
