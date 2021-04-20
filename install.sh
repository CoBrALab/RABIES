#!/bin/bash
### setup RABIES execution and DSURQE atlas

#will install RABIES in HOME directory
export RABIES_VERSION=0.2.1-dev
export RABIES=$PWD

echo "# added by RABIES" >> $HOME/.bashrc
echo "export RABIES_VERSION=0.2.1-dev" >> $HOME/.bashrc
echo "export RABIES=$RABIES" >> $HOME/.bashrc
echo 'export PATH=$PATH:$RABIES/rabies/shell_scripts' >> $HOME/.bashrc

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
python $RABIES/gen_masks.py $template_dir/DSURQE_100micron_labels.nii.gz $template_dir/DSURQE_40micron_R_mapping.csv $template_dir/DSURQE_100micron

# install twolevel_ants_dbm
git clone https://github.com/CoBrALab/twolevel_ants_dbm.git $RABIES/twolevel_ants_dbm
cd $RABIES/twolevel_ants_dbm
git checkout 3623b5bf683473887b7d63c1a1b24c0362ed8396
cd $HOME
echo 'export PATH=$RABIES/twolevel_ants_dbm:$PATH' >> $HOME/.bashrc

# install twolevel_ants_dbm
git clone https://github.com/CoBrALab/minc-toolkit-extras.git $RABIES/minc-toolkit-extras
cd $RABIES/minc-toolkit-extras
git checkout 3dec9b81b0e59c7daa90ea1901469111b2374182
echo 'export PATH=$RABIES/minc-toolkit-extras:$PATH' >> $HOME/.bashrc
cd $HOME
