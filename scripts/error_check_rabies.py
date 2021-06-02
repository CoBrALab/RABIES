#! /usr/bin/env python

import SimpleITK as sitk
import os
import numpy as np
import tempfile
import shutil
import subprocess
from rabies.preprocess_pkg.utils import resample_image_spacing, copyInfo_4DImage

tmppath = tempfile.mkdtemp()
os.makedirs(tmppath+'/inputs')

if 'XDG_DATA_HOME' in os.environ.keys():
    rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
else:
    rabies_path = os.environ['HOME']+'/.local/share/rabies'

template="%s/DSURQE_40micron_average.nii.gz" % (rabies_path)
mask="%s/DSURQE_40micron_mask.nii.gz" % (rabies_path)

img = sitk.ReadImage(template)
spacing = (float(1), float(1), float(1)) # resample to 1mmx1mmx1mm
resampled = resample_image_spacing(sitk.ReadImage(template), spacing)
array = sitk.GetArrayFromImage(resampled)
array_4d = np.repeat(array[:,:,:,np.newaxis], 15, axis=3)
array_4d += np.random.normal(0,array_4d.mean()/100,array_4d.shape) # add gaussian noise
sitk.WriteImage(resampled, tmppath+'/inputs/sub-token_T1w.nii.gz')
sitk.WriteImage(sitk.GetImageFromArray(array_4d, isVector=False), tmppath+'/inputs/sub-token_bold.nii.gz')

resampled = resample_image_spacing(sitk.ReadImage(mask), spacing)
array = sitk.GetArrayFromImage(resampled).astype(bool).astype(int)
binarized = sitk.GetImageFromArray(array, isVector=False)
binarized.CopyInformation(resampled)
sitk.WriteImage(binarized, tmppath+'/inputs/token_mask.nii.gz')
array[:,:,:6] = 0
binarized = sitk.GetImageFromArray(array, isVector=False)
binarized.CopyInformation(resampled)
sitk.WriteImage(binarized, tmppath+'/inputs/token_mask_half.nii.gz')


sitk.WriteImage(copyInfo_4DImage(sitk.ReadImage(tmppath+'/inputs/sub-token_bold.nii.gz'), sitk.ReadImage(tmppath+'/inputs/sub-token_T1w.nii.gz'), sitk.ReadImage(tmppath+'/inputs/sub-token_bold.nii.gz')), tmppath+'/inputs/sub-token_bold.nii.gz')

command = "rabies preprocess %s/inputs %s/outputs --debug --anat_bias_cor_method disable --bold_bias_cor_method disable \
    --anat_template %s/inputs/sub-token_T1w.nii.gz --brain_mask %s/inputs/token_mask.nii.gz --WM_mask %s/inputs/token_mask.nii.gz --CSF_mask %s/inputs/token_mask.nii.gz --vascular_mask %s/inputs/token_mask.nii.gz --labels %s/inputs/token_mask.nii.gz \
    --coreg_script null_nonlin --template_reg_script null_nonlin --data_type int16 -e --detect_dummy \
    --tpattern seq" % (tmppath,tmppath,tmppath,tmppath,tmppath,tmppath,tmppath,tmppath)
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

command = "rabies preprocess %s/inputs %s/outputs --debug --anat_bias_cor_method disable --bold_bias_cor_method disable \
    --anat_template %s/inputs/sub-token_T1w.nii.gz --brain_mask %s/inputs/token_mask.nii.gz --WM_mask %s/inputs/token_mask_half.nii.gz --CSF_mask %s/inputs/token_mask_half.nii.gz --vascular_mask %s/inputs/token_mask_half.nii.gz --labels %s/inputs/token_mask.nii.gz \
    --coreg_script null_nonlin --template_reg_script null_nonlin --data_type int16 --HMC_option 0" % (tmppath,tmppath,tmppath,tmppath,tmppath,tmppath,tmppath,tmppath)
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

command = "rabies confound_regression %s/outputs %s/outputs --run_aroma --FD_censoring --DVARS_censoring" % (tmppath,tmppath)
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

command = "rabies confound_regression %s/outputs %s/outputs --conf_list mot_6 --smoothing_filter 0.3 --commonspace_bold" % (tmppath,tmppath)
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

command = "rabies data_diagnosis %s/outputs %s/outputs --dual_ICA 1" % (tmppath,tmppath)
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

command = "rabies analysis %s/outputs %s/outputs --dual_ICA 1" % (tmppath,tmppath)
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )


shutil.rmtree(tmppath)
