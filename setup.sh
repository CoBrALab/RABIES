### setup RABIES execution and DSURQE atlas

export RABIES="$( cd "$(dirname "$0")" ; pwd -P )"
export PYTHONPATH="${PYTHONPATH}:$RABIES"

mkdir -p $RABIES/bin
echo -e '#! /usr/bin/env python \nfrom rabies.run_main import execute_workflow \nexecute_workflow()' > $RABIES/bin/rabies
chmod +x $RABIES/bin/rabies
echo "# added by Mouse_fmriPype" >> $HOME/.bashrc
echo "export RABIES=$RABIES" >> $HOME/.bashrc
echo 'export PYTHONPATH="${PYTHONPATH}:$RABIES"' >> $HOME/.bashrc
echo 'export PATH=$PATH:$RABIES/bin' >> $HOME/.bashrc


# Download DSURQE template
mkdir -p $HOME/DSURQE_atlas/minc
curl -L --retry 5 -o $HOME/DSURQE_atlas/minc/DSURQE_40micron_average.mnc http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_average.mnc
curl -L --retry 5 -o $HOME/DSURQE_atlas/minc/DSURQE_40micron_labels.mnc http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_labels.mnc
curl -L --retry 5 -o $HOME/DSURQE_atlas/minc/DSURQE_40micron_mask.mnc http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_mask.mnc
curl -L --retry 5 -o $HOME/DSURQE_atlas/DSURQE_40micron_R_mapping.csv http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_R_mapping.csv
# create a 100um version
ResampleImage 3 $HOME/DSURQE_atlas/minc/DSURQE_40micron_average.mnc $HOME/DSURQE_atlas/minc/DSURQE_100micron_average.mnc 0.1x0.1x0.1 0 4
antsApplyTransforms -i $HOME/DSURQE_atlas/minc/DSURQE_40micron_mask.mnc -o $HOME/DSURQE_atlas/minc/DSURQE_100micron_mask.mnc -r $HOME/DSURQE_atlas/minc/DSURQE_100micron_average.mnc -n GenericLabel
antsApplyTransforms -i $HOME/DSURQE_atlas/minc/DSURQE_40micron_labels.mnc -o $HOME/DSURQE_atlas/minc/DSURQE_100micron_labels.mnc -r $HOME/DSURQE_atlas/minc/DSURQE_100micron_average.mnc -n GenericLabel

# convert to nifti
mkdir -p $HOME/DSURQE_atlas/nifti
for file in $HOME/DSURQE_atlas/minc/*; do mnc2nii $file $HOME/DSURQE_atlas/nifti/$(basename ${file::-4}).nii; done
gzip $HOME/DSURQE_atlas/nifti/*.nii

# create WM and CSF masks
DSURQE_100micron_anat=$HOME/DSURQE_atlas/nifti/DSURQE_100micron_average.nii.gz
DSURQE_100micron_mask=$HOME/DSURQE_atlas/nifti/DSURQE_100micron_mask.nii.gz
DSURQE_100micron_labels=$HOME/DSURQE_atlas/nifti/DSURQE_100micron_labels.nii.gz
csv_labels=$HOME/DSURQE_atlas/DSURQE_40micron_R_mapping.csv
python $RABIES/gen_masks.py $DSURQE_100micron_labels $csv_labels $HOME/DSURQE_atlas/nifti/DSURQE_100micron
