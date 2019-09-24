#!bin/bash

FILE_PATH=$1
template=$template_anat
mask=$template_mask

mkdir -p ants_dbm
cd ants_dbm

echo twolevel_dbm.py --rigid-model-target $atlas_labels \
--no-N4 --transform SyN --float --average-type normmean --gradient-step 0.25 --model-iterations 3 \
--modelbuild-command antsMultivariateTemplateConstruction2.sh --cluster-type $ants_dbm_cluster_type \
--walltime $ants_dbm_walltime --memory-request $ants_dbm_memory_request --local-threads $ants_dbm_local_threads \
1level $FILE_PATH > exec.sh
echo "Running the following commonspace registration:"
cat exec.sh
bash exec.sh

cd ..


if test -f "template_reg/template_reg_InverseComposite.h5"
then
  echo "Registration previously run. Using previous transforms for further steps."
else
  echo "Running template registration."
  mkdir -p template_reg
  cd template_reg

  moving=../ants_dbm/output/secondlevel/secondlevel_template0.nii.gz

  antsRegistration --dimensionality 3 --verbose 1 \
    --output [template_reg_,template_reg_warped_image.nii.gz] \
    --initial-moving-transform [$template,$moving,1] \
    --transform Rigid[0.1] --metric Mattes[$template,$moving,1,32,None] --convergence [2025x2025x2025x2025x2025,1e-6,10] --shrink-factors 16x15x14x13x12 --smoothing-sigmas 7.99225592362x7.49173910041x6.99114831402x6.49046645078x5.98967067114vox --masks [NULL,NULL] \
    --transform Rigid[0.1] --metric Mattes[$template,$moving,1,64,None] --convergence [2025x2025x2025x2025x2025,1e-6,10] --shrink-factors 12x11x10x9x8 --smoothing-sigmas 5.98967067114x5.48872979374x4.98760009911x4.48621831264x3.98448927075vox --masks [NULL,NULL] \
    --transform Rigid[0.1] --metric Mattes[$template,$moving,1,128,None] --convergence [2025x2025x2025x2025x675,1e-6,10] --shrink-factors 8x7x6x5x4 --smoothing-sigmas 3.98448927075x3.4822628776x2.97928762436x2.47510701762x1.96879525311vox --masks [$mask,NULL] \
    --transform Rigid[0.1] --metric Mattes[$template,$moving,1,256,None] --convergence [675x225x200x200x200,1e-6,10] --shrink-factors 4x3x2x1x1 --smoothing-sigmas 1.96879525311x1.45813399545x0.936031382318x0.355182697615x0vox --masks [NULL,NULL] \
    --transform Similarity[0.1] --metric Mattes[$template,$moving,1,256,None] --convergence [675x225x200x200x200,1e-6,10] --shrink-factors 4x3x2x1x1 --smoothing-sigmas 1.96879525311x1.45813399545x0.936031382318x0.355182697615x0vox --masks [$mask,NULL] \
    --transform Affine[0.1] --metric Mattes[$template,$moving,1,256,None] --convergence [675x225x200x200x200,1e-6,10] --shrink-factors 4x3x2x1x1 --smoothing-sigmas 1.96879525311x1.45813399545x0.936031382318x0.355182697615x0vox --masks [$mask,NULL] \
    --transform SyN[0.2,2,0] --metric CC[$template,$moving,1,4] --convergence [2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x200x50x25x25,1e-6,10] \
      --shrink-factors 8x8x8x8x8x8x8x8x8x7x6x5x4x3x2x1x1 \
      --smoothing-sigmas 7.99225592362x7.49173910041x6.99114831402x6.49046645078x5.98967067114x5.48872979374x4.98760009911x4.48621831264x3.98448927075x3.4822628776x2.97928762436x2.47510701762x1.96879525311x1.45813399545x0.936031382318x0.355182697615x0 \
      --masks [NULL,NULL] -z 1 -a 1

  cd ..
fi

if test -f "complete.txt"
then
  echo "Transforms of brain mask and labels was previously completed. Using previous ouputs."
else
  echo "Running transform of subject anat mask and labels."

  #move brain mask and labels to each of subject space
  antsApplyTransforms -i $mask -t template_reg/template_reg_InverseComposite.h5 -r ants_dbm/output/secondlevel/secondlevel_template0.nii.gz -o template_reg/transformed_mask.nii.gz -n GenericLabel
  antsApplyTransforms -i $atlas_labels -t template_reg/template_reg_InverseComposite.h5 -r ants_dbm/output/secondlevel/secondlevel_template0.nii.gz -o template_reg/transformed_labels.nii.gz -n GenericLabel
  #also WM and CSF masks
  antsApplyTransforms -i $WM_mask -t template_reg/template_reg_InverseComposite.h5 -r ants_dbm/output/secondlevel/secondlevel_template0.nii.gz -o template_reg/transformed_WM_mask.nii.gz -n GenericLabel
  antsApplyTransforms -i $CSF_mask -t template_reg/template_reg_InverseComposite.h5 -r ants_dbm/output/secondlevel/secondlevel_template0.nii.gz -o template_reg/transformed_CSF_mask.nii.gz -n GenericLabel

  template_folder=../anat_datasink/commonspace_template/
  mkdir -p $template_folder
  cp template_reg/transformed_mask.nii.gz $template_folder
  cp template_reg/transformed_labels.nii.gz $template_folder
  cp ants_dbm/output/secondlevel/secondlevel_template0.nii.gz $template_folder/commonspace_template.nii.gz

  while read anat_file; do file=$(basename $anat_file); IFS='_' read -r -a array <<< "$file"; sub=${array[0]}; ses=${array[1]:4:1}; \
  sub_dir=../anat_datasink/commonspace_transforms/_subject_id_${sub}/_session_${ses}/ ; mkdir -p $sub_dir ; mkdir -p ../anat_datasink/anat_mask/_subject_id_${sub}/_session_${ses}/; mkdir -p ../anat_datasink/anat_labels/_subject_id_${sub}/_session_${ses}/; \
  mkdir -p ../anat_datasink/WM_mask/_subject_id_${sub}/_session_${ses}/; mkdir -p ../anat_datasink/CSF_mask/_subject_id_${sub}/_session_${ses}/; \
  inverse_warp=ants_dbm/output/secondlevel/secondlevel_${sub}_ses-${ses}_anat_preproc*InverseWarp.nii.gz; \
  str="$(basename $inverse_warp)"; delimiter=1InverseWarp; s=$str$delimiter; array=(); \
  while [[ $s ]]; do
      array+=( "${s%%"$delimiter"*}" );
      s=${s#*"$delimiter"};
  done; spec=${array[0]}; warp=ants_dbm/output/secondlevel/${spec}1Warp.nii.gz; affine=ants_dbm/output/secondlevel/${spec}0GenericAffine.mat; inverse_warp=ants_dbm/output/secondlevel/${spec}1InverseWarp.nii.gz; \
  echo "cp $affine $sub_dir/${sub}_ses-${ses}_Affine.mat; cp $warp $sub_dir/${sub}_ses-${ses}_Warp.nii.gz; cp $inverse_warp $sub_dir/${sub}_ses-${ses}_InverseWarp.nii.gz; \
  antsApplyTransforms -i template_reg/transformed_labels.nii.gz -t [$affine,1] -t $inverse_warp -r $anat_file -o ../anat_datasink/anat_labels/_subject_id_${sub}/_session_${ses}/${sub}_ses-${ses}_anat_labels.nii.gz -n GenericLabel; \
  antsApplyTransforms -i template_reg/transformed_mask.nii.gz -t [$affine,1] -t $inverse_warp -r $anat_file -o ../anat_datasink/anat_mask/_subject_id_${sub}/_session_${ses}/${sub}_ses-${ses}_anat_mask.nii.gz -n GenericLabel; \
  antsApplyTransforms -i template_reg/transformed_WM_mask.nii.gz -t [$affine,1] -t $inverse_warp -r $anat_file -o ../anat_datasink/WM_mask/_subject_id_${sub}/_session_${ses}/${sub}_ses-${ses}_WM_mask.nii.gz -n GenericLabel; \
  antsApplyTransforms -i template_reg/transformed_CSF_mask.nii.gz -t [$affine,1] -t $inverse_warp -r $anat_file -o ../anat_datasink/CSF_mask/_subject_id_${sub}/_session_${ses}/${sub}_ses-${ses}_CSF_mask.nii.gz -n GenericLabel"; \
  done < $FILE_PATH > transform_jobs.sh

  qbatch -b $ants_dbm_cluster_type -o "-sync y" transform_jobs.sh

  echo "Complete transforms." > complete.txt
fi
