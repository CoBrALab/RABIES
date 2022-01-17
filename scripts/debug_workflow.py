#! /usr/bin/env python

import SimpleITK as sitk
import os
import sys
import numpy as np
import tempfile
import shutil
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


output_folder = f'{tmppath}/outputs'
args = [
        #'--debug',
        'preprocess', f'{tmppath}/inputs', output_folder,
        '--anat_inho_cor_method','disable', '--bold_inho_cor_method', 'disable', 
        '--coreg_script', 'NULL', '--atlas_reg_script', 'NULL', '--data_type', 'int16', '--bold_only', '--fast_commonspace',
        '--anat_template', f'{tmppath}/inputs/sub-token_T1w.nii.gz',
        '--brain_mask', f'{tmppath}/inputs/token_mask.nii.gz', 
        '--WM_mask', f'{tmppath}/inputs/token_mask.nii.gz',
        '--CSF_mask', f'{tmppath}/inputs/token_mask.nii.gz',
        '--vascular_mask', f'{tmppath}/inputs/token_mask.nii.gz', 
        '--labels', f'{tmppath}/inputs/token_mask.nii.gz',
        ]

from rabies.parser import get_parser

parser = get_parser()

opts = parser.parse_args(args)

if not os.path.isdir(output_folder):
    os.makedirs(output_folder)
        
from rabies.run_main import prep_logging, preprocess, confound_correction, analysis
log = prep_logging(opts, output_folder)
    
# print complete CLI command
args = 'CLI INPUTS: \n'
for arg in vars(opts):
    input = f'-> {arg} = {getattr(opts, arg)} \n'
    args += input
log.info(args)


if opts.rabies_step == 'preprocess':
    workflow = preprocess(opts, log)
elif opts.rabies_step == 'confound_correction':
    workflow = confound_correction(opts, log)
elif opts.rabies_step == 'analysis':
    workflow = analysis(opts, log)
else:
    parser.print_help()
workflow.base_dir = output_folder


log.info(f'Running workflow with {opts.plugin} plugin.')
# execute workflow, with plugin_args limiting the cluster load for parallel execution
graph_out = workflow.run(plugin=opts.plugin, plugin_args={'max_jobs': 50, 'dont_resubmit_completed_jobs': True,
                                              'n_procs': opts.local_threads, 'qsub_args': f'-pe smp {str(opts.min_proc)}'})


