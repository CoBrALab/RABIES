mkdir -p realigned_files

for file in anat_preprocessing/anat_preproc/*anat.mnc; do antsApplyTransforms -d 3 -i $file -t /data/chamal/projects/Gabriel_DG/utils_scripts/anat2EPI_realign/anat2EPI0GenericAffine.mat \
-r /data/chamal/projects/Gabriel_DG/utils_scripts/anat2EPI_realign/ref_space_img.mnc -o realigned_files/$(basename $file) --verbose -n BSpline; done

for file in anat_preprocessing/mask_labels_files/*.mnc; do antsApplyTransforms -d 3 -i $file -t /data/chamal/projects/Gabriel_DG/utils_scripts/anat2EPI_realign/anat2EPI0GenericAffine.mat \
-r /data/chamal/projects/Gabriel_DG/utils_scripts/anat2EPI_realign/ref_space_img.mnc -o realigned_files/$(basename $file) --verbose -n GenericLabel; done

for file in anat_preprocessing/WM_CSF_masks/*mask_eroded.mnc; do antsApplyTransforms -d 3 -i $file -t /data/chamal/projects/Gabriel_DG/utils_scripts/anat2EPI_realign/anat2EPI0GenericAffine.mat \
-r /data/chamal/projects/Gabriel_DG/utils_scripts/anat2EPI_realign/ref_space_img.mnc -o realigned_files/$(basename $file) --verbose -n GenericLabel; done
