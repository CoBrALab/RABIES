source activate custom_nipype

mkdir -p WM_CSF_masks
for file in mask_labels_files/*labels_resample.mnc; do var=$(basename $file); python /data/chamal/projects/Gabriel_DG/utils_scripts/compute_masks.py \
$file WM_CSF_masks/${var:0:8} ; done

for file in WM_CSF_masks/*_mask.nii.gz; do /data/chamal/projects/Gabriel_DG/utils_scripts/erode_mask.sh $file WM_CSF_masks/$(basename $file .nii.gz)_eroded.mnc; done
