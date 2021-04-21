#!/bin/bash
# Download DSURQE template
template_dir=$RABIES/template_files
curl -L --retry 5 -o $template_dir/DSURQE_40micron_average.nii http://repo.mouseimaging.ca/repo/DSURQE_40micron_nifti/DSURQE_40micron_average.nii
curl -L --retry 5 -o $template_dir/DSURQE_40micron_labels.nii http://repo.mouseimaging.ca/repo/DSURQE_40micron_nifti/DSURQE_40micron_labels.nii
#curl -L --retry 5 -o $template_dir/DSURQE_40micron_mask.nii http://repo.mouseimaging.ca/repo/DSURQE_40micron_nifti/DSURQE_40micron_mask.nii
curl -L --retry 5 -o $template_dir/DSURQE_40micron_R_mapping.csv http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_R_mapping.csv

# create a 100um version
ResampleImage 3 $template_dir/DSURQE_40micron_average.nii $template_dir/DSURQE_100micron_average.nii 0.1x0.1x0.1 0 4
#antsApplyTransforms -i $template_dir/DSURQE_40micron_mask.nii -o $template_dir/DSURQE_100micron_mask.nii -r $template_dir/DSURQE_100micron_average.nii -n GenericLabel
antsApplyTransforms -i $template_dir/DSURQE_40micron_labels.nii -o $template_dir/DSURQE_100micron_labels.nii -r $template_dir/DSURQE_100micron_average.nii -n GenericLabel

gzip -f $template_dir/*.nii

#convert templates to the RAS axis convention
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_average.nii.gz
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_mask.nii.gz
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_labels.nii.gz

# create regional masks
gen_DSURQE_masks.py $template_dir/DSURQE_100micron_labels.nii.gz $template_dir/DSURQE_40micron_R_mapping.csv $template_dir/DSURQE_100micron
