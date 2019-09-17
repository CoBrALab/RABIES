#!bin/bash

FILE_PATH=$1

mkdir -p pydpiper
cd pydpiper

MBM.py --verbose --pipeline-name=mbm_atlasReg \
--num-executors 150 \
--lsq6-large-rotations \
--no-nuc \
--init-model /opt/quarantine/resources/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/100um/DSURQE.mnc \
--maget-atlas-library /opt/quarantine/resources/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/ex-vivo/ \
--maget-no-pairwise \
--maget-mask \
--common-space-model /opt/quarantine/resources/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/ex-vivo/DSURQE_40micron.mnc \
--run-maget \
--csv-file $FILE_PATH

cd ..

if test -f "complete.txt"
then
  echo "Transforms of brain mask and labels was previously completed. Using previous ouputs."
else
  echo "Running transform of subject anat mask and labels."

  mkdir -p minc_files

  #move brain mask and labels to each of subject space
  while read anat_file; do if test $anat_file == file
  then
    continue
  fi; echo $anat_file; file=$(basename $anat_file); IFS='_' read -r -a array <<< "$file"; sub=${array[0]}; ses=${array[1]:4:1}; \
  mkdir -p ../anat_datasink/anat_mask/_subject_id_${sub}/_session_${ses}/; mkdir -p ../anat_datasink/anat_labels/_subject_id_${sub}/_session_${ses}/; \
  voted_labels=pydpiper/mbm_atlasReg_processed/${sub}_ses-${ses}_anat_preproc/voted.mnc; lsq6_mask=pydpiper/mbm_atlasReg_atlases/DSURQE_40micron_mask/tmp/${sub}_ses-${ses}_anat_preproc_I_lsq6_max_mask.mnc; lsq6_transform=pydpiper/mbm_atlasReg_processed/${sub}_ses-${ses}_anat_preproc/transforms/${sub}_ses-${ses}_anat_preproc_lsq6.xfm; \
  echo "antsApplyTransforms -i $lsq6_mask -t [$lsq6_transform,1] -r $anat_file -o minc_files/${sub}_ses-${ses}_anat_mask.mnc -n GenericLabel; \
  antsApplyTransforms -i $voted_labels -t [$lsq6_transform,1] -r $anat_file -o minc_files/${sub}_ses-${ses}_anat_labels.mnc -n GenericLabel; \
  mnc2nii minc_files/${sub}_ses-${ses}_anat_mask.mnc ../anat_datasink/anat_mask/_subject_id_${sub}/_session_${ses}/${sub}_ses-${ses}_anat_mask.nii.gz; \
  mnc2nii minc_files/${sub}_ses-${ses}_anat_labels.mnc ../anat_datasink/anat_labels/_subject_id_${sub}/_session_${ses}/${sub}_ses-${ses}_anat_labels.nii.gz"; \
  done < $FILE_PATH > transform_jobs.sh

  qbatch -o "-sync y" transform_jobs.sh

  echo "Complete transforms." > complete.txt
fi
