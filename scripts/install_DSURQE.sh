#!/bin/bash
out_dir=$1
tmpdir=$(mktemp -d)
mkdir -p $out_dir

# Download DSURQE template
curl -L --retry 5 -o ${tmpdir}/DSURQE_40micron_average.mnc http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_average.mnc
curl -L --retry 5 -o ${tmpdir}/DSURQE_40micron_labels.mnc http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_labels.mnc
curl -L --retry 5 -o ${tmpdir}/DSURQE_40micron_mask.mnc http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_mask.mnc
curl -L --retry 5 -o ${out_dir}/DSURQE_40micron_R_mapping.csv http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_R_mapping.csv

# convert to nifti
mnc2nii ${tmpdir}/DSURQE_40micron_average.mnc ${out_dir}/DSURQE_40micron_average.nii
mnc2nii ${tmpdir}/DSURQE_40micron_labels.mnc ${out_dir}/DSURQE_40micron_labels.nii
mnc2nii ${tmpdir}/DSURQE_40micron_mask.mnc ${out_dir}/DSURQE_40micron_mask.nii

gzip -f ${out_dir}/DSURQE_40micron_*.nii

# create regional masks
gen_DSURQE_masks.py ${out_dir}/DSURQE_40micron_labels.nii.gz ${out_dir}/DSURQE_40micron_R_mapping.csv ${out_dir}/DSURQE_40micron

# download additional supportive files
curl -L --retry 5 https://zenodo.org/record/4706153/files/melodic_IC.nii.gz -o ${out_dir}/melodic_IC.nii.gz
curl -L --retry 5 https://zenodo.org/record/4706153/files/vascular_mask.nii.gz -o ${out_dir}/vascular_mask.nii.gz

rm -rf ${tmpdir}
