module load anaconda/5.0.1-python3 minc-toolkit/1.9.15 minc-stuffs/0.1.21^minc-toolkit-1.9.15 qbatch/git pydpiper/2.0.10 minc-toolkit-extras
module load qbatch

module load AFNI ANTs
source activate /home/cic/desgab/miniconda3/envs/custom_nipype
PYTHONPATH=$PYTHONPATH:/home/cic/desgab/Desktop/mouse_procfmri

bash nii2mnc.sh
cd anat_preprocessing
bash preproc_anat.sh
cd pydpyper
bash script_mbm.sh
cd ..
bash transform_mask_labels.sh
bash compute_WM_CSF_masks.sh
cd ..
bash realign_space.sh
bash mnc2nii.sh

for file in $PWD/nii_files/*main_EPI.nii.gz; \
do var=$(basename $file); echo python /home/cic/desgab/Desktop/mouse_procfmri/run_main.py $file ${var:0:7} \
$PWD/nii_files/${var:0:8}preproc_anat.nii.gz \
$PWD/nii_files/${var:0:8}mask_resample.nii.gz \
$PWD/nii_files/${var:0:8}labels_resample.nii.gz \
$PWD/nii_files/${var:0:8}WM_mask_eroded.nii.gz \
$PWD/nii_files/${var:0:8}CSF_mask_eroded.nii.gz \
/data/chamal/projects/Gabriel_DG/rodent_pipeline/outputs/iso_16_bits; done > run_fmri_pipeline.sh

module load qbatch
qbatch --ppj 4 run_fmri_pipeline.sh
