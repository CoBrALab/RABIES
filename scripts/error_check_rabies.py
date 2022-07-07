#! /usr/bin/env python

import SimpleITK as sitk
import os
import sys
import numpy as np
import tempfile
import shutil
import subprocess
from rabies.utils import resample_image_spacing, copyInfo_4DImage

if len(sys.argv) == 2:
    tmppath = sys.argv[1]
else:
    tmppath = tempfile.mkdtemp()

os.makedirs(tmppath+'/inputs', exist_ok=True)

if 'XDG_DATA_HOME' in os.environ.keys():
    rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
else:
    rabies_path = os.environ['HOME']+'/.local/share/rabies'

template = f"{rabies_path}/DSURQE_40micron_average.nii.gz"
mask = f"{rabies_path}/DSURQE_40micron_mask.nii.gz"

img = sitk.ReadImage(template)
spacing = (float(1), float(1), float(1))  # resample to 1mmx1mmx1mm
resampled = resample_image_spacing(sitk.ReadImage(template), spacing)
array = sitk.GetArrayFromImage(resampled)
array_4d = np.repeat(array[np.newaxis, :, :, :], 15, axis=0)
array_4d += np.random.normal(0, array_4d.mean()
                             / 100, array_4d.shape)  # add gaussian noise
sitk.WriteImage(resampled, tmppath+'/inputs/sub-token_T1w.nii.gz')
sitk.WriteImage(resampled, tmppath+'/inputs/sub-token2_T1w.nii.gz')
sitk.WriteImage(resampled, tmppath+'/inputs/sub-token3_T1w.nii.gz')
sitk.WriteImage(sitk.GetImageFromArray(array_4d, isVector=False),
                tmppath+'/inputs/sub-token_bold.nii.gz')

resampled = resample_image_spacing(sitk.ReadImage(mask), spacing)
array = sitk.GetArrayFromImage(resampled)
array[array < 1] = 0
array[array > 1] = 1
binarized = sitk.GetImageFromArray(array, isVector=False)
binarized.CopyInformation(resampled)
sitk.WriteImage(binarized, tmppath+'/inputs/token_mask.nii.gz')
array[:, :, :6] = 0
binarized = sitk.GetImageFromArray(array, isVector=False)
binarized.CopyInformation(resampled)
sitk.WriteImage(binarized, tmppath+'/inputs/token_mask_half.nii.gz')


sitk.WriteImage(copyInfo_4DImage(sitk.ReadImage(tmppath+'/inputs/sub-token_bold.nii.gz'), sitk.ReadImage(tmppath
                + '/inputs/sub-token_T1w.nii.gz'), sitk.ReadImage(tmppath+'/inputs/sub-token_bold.nii.gz')), tmppath+'/inputs/sub-token_bold.nii.gz')

command = f"rabies --verbose 1 preprocess {tmppath}/inputs {tmppath}/outputs --anat_inho_cor method=disable,otsu_thresh=2,multiotsu=false --bold_inho_cor method=disable,otsu_thresh=2,multiotsu=false \
    --anat_template {tmppath}/inputs/sub-token_T1w.nii.gz --brain_mask {tmppath}/inputs/token_mask.nii.gz --WM_mask {tmppath}/inputs/token_mask.nii.gz --CSF_mask {tmppath}/inputs/token_mask.nii.gz --vascular_mask {tmppath}/inputs/token_mask.nii.gz --labels {tmppath}/inputs/token_mask.nii.gz \
    --bold2anat_coreg registration=no_reg,masking=false,brain_extraction=false --commonspace_reg masking=false,brain_extraction=false,fast_commonspace=true,template_registration=no_reg --data_type int16 --bold_only --detect_dummy \
    --tpattern seq-z"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

os.remove(f'{tmppath}/outputs/rabies_preprocess.pkl')
command = f"rabies --verbose 1 preprocess {tmppath}/inputs {tmppath}/outputs --anat_inho_cor method=disable,otsu_thresh=2,multiotsu=false --bold_inho_cor method=disable,otsu_thresh=2,multiotsu=false \
    --anat_template {tmppath}/inputs/sub-token_T1w.nii.gz --brain_mask {tmppath}/inputs/token_mask.nii.gz --WM_mask {tmppath}/inputs/token_mask.nii.gz --CSF_mask {tmppath}/inputs/token_mask.nii.gz --vascular_mask {tmppath}/inputs/token_mask.nii.gz --labels {tmppath}/inputs/token_mask.nii.gz \
    --bold2anat_coreg registration=no_reg,masking=true,brain_extraction=true --commonspace_reg masking=true,brain_extraction=true,fast_commonspace=false,template_registration=no_reg --data_type int16  \
    --HMC_option 0"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

command = f"rabies --verbose 1 confound_correction {tmppath}/outputs {tmppath}/outputs --ica_aroma apply=true,dim=0,random_seed=1 --frame_censoring FD_censoring=true,FD_threshold=0.05,DVARS_censoring=true,minimum_timepoint=3 --nativespace_analysis"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

### Add subjects for the group analysis to run
sitk.WriteImage(sitk.GetImageFromArray(array_4d, isVector=False),
                tmppath+'/inputs/sub-token2_bold.nii.gz')
sitk.WriteImage(sitk.GetImageFromArray(array_4d, isVector=False),
                tmppath+'/inputs/sub-token3_bold.nii.gz')

os.remove(f'{tmppath}/outputs/rabies_preprocess.pkl')
command = f"rabies --verbose 1 preprocess {tmppath}/inputs {tmppath}/outputs --anat_inho_cor method=disable,otsu_thresh=2,multiotsu=false --bold_inho_cor method=disable,otsu_thresh=2,multiotsu=false \
    --anat_template {tmppath}/inputs/sub-token_T1w.nii.gz --brain_mask {tmppath}/inputs/token_mask.nii.gz --WM_mask {tmppath}/inputs/token_mask_half.nii.gz --CSF_mask {tmppath}/inputs/token_mask_half.nii.gz --vascular_mask {tmppath}/inputs/token_mask_half.nii.gz --labels {tmppath}/inputs/token_mask.nii.gz \
    --bold2anat_coreg registration=no_reg,masking=false,brain_extraction=false --commonspace_reg masking=false,brain_extraction=false,fast_commonspace=true,template_registration=no_reg --data_type int16  \
    --HMC_option 0"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

os.remove(f'{tmppath}/outputs/rabies_confound_correction.pkl')
command = f"rabies --verbose 1 confound_correction --read_datasink {tmppath}/outputs {tmppath}/outputs --conf_list mot_6 --smoothing_filter 0.3"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

command = f"rabies --verbose 1 analysis {tmppath}/outputs {tmppath}/outputs --DR_ICA --NPR_temporal_comp 1 --seed_list {tmppath}/inputs/token_mask_half.nii.gz"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

os.remove(f'{tmppath}/outputs/rabies_confound_correction.pkl')
command = f"rabies --verbose 1 confound_correction {tmppath}/outputs {tmppath}/outputs"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

os.remove(f'{tmppath}/outputs/rabies_analysis.pkl')
command = f"rabies --verbose 1 analysis {tmppath}/outputs {tmppath}/outputs --NPR_temporal_comp 1 --data_diagnosis --DR_ICA"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )


if not len(sys.argv) == 2:
    shutil.rmtree(tmppath)
