mkdir -p minc
mkdir -p nii_files

for file in ../raw_nii/*; do nii2mnc $file/*_7/*.nii.gz minc/$(basename $file)_anat.mnc; \
cp $file/*_5/*.nii.gz nii_files/$(basename $file)_main_EPI.nii.gz; \
cp $file/*_6/*.nii.gz nii_files/$(basename $file)_rev_EPI.nii.gz; done
