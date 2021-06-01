#! /usr/bin/env python

import SimpleITK as sitk
import os
import numpy as np
import tempfile
import shutil
import subprocess
from rabies.preprocess_pkg.utils import resample_image_spacing

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
sitk.WriteImage(resampled, tmppath+'/inputs/sub-token_T1w.nii.gz')
sitk.WriteImage(sitk.GetImageFromArray(array_4d, isVector=False), tmppath+'/inputs/sub-token_bold.nii.gz')
resampled = resample_image_spacing(sitk.ReadImage(mask), spacing)
sitk.WriteImage(resampled, tmppath+'/inputs/token_mask.nii.gz')

command = "rabies preprocess %s/inputs %s/outputs --debug --anat_bias_cor_method disable --bold_bias_cor_method disable \
    --coreg_script null_nonlin --template_reg_script null_nonlin --data_type int16 -e --detect_dummy \
    --HMC_transform Affine --tpattern seq" % (tmppath,tmppath)
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )


command = "rabies preprocess %s/inputs %s/outputs --debug --anat_bias_cor_method disable --bold_bias_cor_method disable \
    --coreg_script null_nonlin --template_reg_script null_nonlin --data_type int16" % (tmppath,tmppath)
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )



shutil.rmtree(tmppath)
