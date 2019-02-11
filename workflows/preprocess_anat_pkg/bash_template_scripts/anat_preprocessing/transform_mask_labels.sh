mkdir -p mask_labels_files

for file in anat_preproc/*anat.mnc; do var=$(basename $file); \
antsApplyTransforms -d 3 -i pydpyper/mbm_atlasReg_atlases/DSURQE_40micron_mask/resampled/${var:0:8}preproc_anat_I_lsq6_lsq12_and_nlin-resampled_mask.mnc \
-t pydpyper/mbm_atlasReg_processed/${var:0:8}preproc_anat/transforms/${var:0:8}preproc_anat__concat_lsq6_I_lsq6_lsq12_and_nlin.xfm \
-r $file -o mask_labels_files/${var:0:8}mask_resample.mnc --verbose -n GenericLabel; done

for file in anat_preproc/*anat.mnc; do var=$(basename $file); \
antsApplyTransforms -d 3 -i pydpyper/mbm_atlasReg_processed/${var:0:8}preproc_anat/voted.mnc \
-t pydpyper/mbm_atlasReg_processed/${var:0:8}preproc_anat/transforms/${var:0:8}preproc_anat__concat_lsq6_I_lsq6_lsq12_and_nlin.xfm \
-r $file -o mask_labels_files/${var:0:8}labels_resample.mnc --verbose -n GenericLabel; done
