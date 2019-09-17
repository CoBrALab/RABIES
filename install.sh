### setup RABIES execution and DSURQE atlas

#will install RABIES in HOME directory
export RABIES=$HOME/RABIES
export PYTHONPATH="${PYTHONPATH}:$RABIES"

git clone https://github.com/Gab-D-G/RABIES $HOME/RABIES

mkdir -p $RABIES/bin
echo -e '#! /usr/bin/env python \nfrom rabies.run_main import execute_workflow \nexecute_workflow()' > $RABIES/bin/rabies
chmod +x $RABIES/bin/rabies
echo "# added by RABIES" >> $HOME/.bashrc
echo "export RABIES=$RABIES" >> $HOME/.bashrc
echo 'export PYTHONPATH="${PYTHONPATH}:$RABIES"' >> $HOME/.bashrc
echo 'export PATH=$PATH:$RABIES/bin' >> $HOME/.bashrc


# Download DSURQE template
mkdir -p $RABIES/DSURQE_atlas/minc
curl -L --retry 5 -o $RABIES/DSURQE_atlas/minc/DSURQE_40micron_average.mnc http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_average.mnc
curl -L --retry 5 -o $RABIES/DSURQE_atlas/minc/DSURQE_40micron_labels.mnc http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_labels.mnc
curl -L --retry 5 -o $RABIES/DSURQE_atlas/minc/DSURQE_40micron_mask.mnc http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_mask.mnc
curl -L --retry 5 -o $RABIES/DSURQE_atlas/DSURQE_40micron_R_mapping.csv http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_R_mapping.csv
# create a 100um version
ResampleImage 3 $RABIES/DSURQE_atlas/minc/DSURQE_40micron_average.mnc $RABIES/DSURQE_atlas/minc/DSURQE_100micron_average.mnc 0.1x0.1x0.1 0 4
antsApplyTransforms -i $RABIES/DSURQE_atlas/minc/DSURQE_40micron_mask.mnc -o $RABIES/DSURQE_atlas/minc/DSURQE_100micron_mask.mnc -r $RABIES/DSURQE_atlas/minc/DSURQE_100micron_average.mnc -n GenericLabel
antsApplyTransforms -i $RABIES/DSURQE_atlas/minc/DSURQE_40micron_labels.mnc -o $RABIES/DSURQE_atlas/minc/DSURQE_100micron_labels.mnc -r $RABIES/DSURQE_atlas/minc/DSURQE_100micron_average.mnc -n GenericLabel

# convert to nifti
mkdir -p $RABIES/DSURQE_atlas/nifti
for file in $RABIES/DSURQE_atlas/minc/*; do mnc2nii $file $RABIES/DSURQE_atlas/nifti/$(basename ${file::-4}).nii; done
gzip -f $RABIES/DSURQE_atlas/nifti/*.nii

# create WM and CSF masks
DSURQE_100micron_anat=$RABIES/DSURQE_atlas/nifti/DSURQE_100micron_average.nii.gz
DSURQE_100micron_mask=$RABIES/DSURQE_atlas/nifti/DSURQE_100micron_mask.nii.gz
DSURQE_100micron_labels=$RABIES/DSURQE_atlas/nifti/DSURQE_100micron_labels.nii.gz
csv_labels=$RABIES/DSURQE_atlas/DSURQE_40micron_R_mapping.csv
python $RABIES/gen_masks.py $DSURQE_100micron_labels $csv_labels $RABIES/DSURQE_atlas/nifti/DSURQE_100micron

# install twolevel_ants_dbm
git clone https://github.com/CobraLab/twolevel_ants_dbm $RABIES/twolevel_ants_dbm && \
echo 'export PATH=$RABIES/twolevel_ants_dbm/twolevel_dbm.py:$PATH' >> $HOME/.bashrc
