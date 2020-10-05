### setup RABIES execution and DSURQE atlas

#will install RABIES in HOME directory
export RABIES_VERSION=0.2.0-dev
export RABIES=$HOME/RABIES-${RABIES_VERSION}
export PYTHONPATH="${PYTHONPATH}:$RABIES"

cd $HOME
#curl -L --retry 5 -o $HOME/RABIES.tar.gz https://github.com/CoBrALab/RABIES/archive/${RABIES_VERSION}.tar.gz
#tar zxf $HOME/RABIES.tar.gz
#rm $HOME/RABIES.tar.gz
git clone https://github.com/CoBrALab/RABIES
mv RABIES $RABIES

#creates an executable script to execute rabies
mkdir -p $RABIES/bin
echo -e '#! /usr/bin/env python \nfrom rabies.run_main import execute_workflow \nexecute_workflow()' > $RABIES/bin/rabies
chmod +x $RABIES/bin/rabies
echo "# added by RABIES" >> $HOME/.bashrc
echo "export RABIES_VERSION=0.2.0-dev" >> $HOME/.bashrc
echo "export RABIES=$RABIES" >> $HOME/.bashrc
echo 'export PYTHONPATH="${PYTHONPATH}:$RABIES"' >> $HOME/.bashrc
echo 'export PATH=$PATH:$RABIES/bin' >> $HOME/.bashrc
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

# create WM and CSF masks
DSURQE_100micron_anat=$template_dir/DSURQE_100micron_average.nii.gz
DSURQE_100micron_mask=$template_dir/DSURQE_100micron_mask.nii.gz
DSURQE_100micron_labels=$template_dir/DSURQE_100micron_labels.nii.gz
csv_labels=$template_dir/DSURQE_40micron_R_mapping.csv
python $RABIES/gen_masks.py $DSURQE_100micron_labels $csv_labels $template_dir/DSURQE_100micron

#convert templates to the RAS axis convention
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_average.nii.gz
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_mask.nii.gz
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_labels.nii.gz
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_WM_mask.nii.gz
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_CSF_mask.nii.gz
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_eroded_WM_mask.nii.gz
python $RABIES/convert_to_RAS.py $template_dir/DSURQE_100micron_eroded_CSF_mask.nii.gz

# install twolevel_ants_dbm
git clone https://github.com/CoBrALab/twolevel_ants_dbm.git $RABIES/twolevel_ants_dbm
cd $RABIES/twolevel_ants_dbm
git checkout 0b40b07075be1d921a102b162e8713baeaf44484
cd $HOME
echo 'export PATH=$RABIES/twolevel_ants_dbm:$PATH' >> $HOME/.bashrc
