#! /usr/bin/env python

from rabies.utils import generate_token_data

import tempfile
tmppath = tempfile.mkdtemp()

#### increase the number of scans generated to 3 if running group analysis
generate_token_data(tmppath, number_scans=1)

output_folder = f'{tmppath}/outputs'

#### HERE ARE SET THE DESIRED PARAMETERS FOR PREPROCESSING
args = [
        #'--debug',
        'preprocess', f'{tmppath}/inputs', output_folder,
        '--anat_inho_cor', 'method=disable,otsu_thresh=2,multiotsu=false', 
        '--bold_inho_cor', 'method=disable,otsu_thresh=2,multiotsu=false',
        '--bold2anat_coreg', 'registration=no_reg,masking=false,brain_extraction=false', 
        '--commonspace_reg', 'masking=false,brain_extraction=false,fast_commonspace=true,template_registration=no_reg', 
        '--data_type', 'int16', '--bold_only',
        '--anat_template', f'{tmppath}/inputs/sub-token1_T1w.nii.gz',
        '--brain_mask', f'{tmppath}/inputs/token_mask.nii.gz', 
        '--WM_mask', f'{tmppath}/inputs/token_mask.nii.gz',
        '--CSF_mask', f'{tmppath}/inputs/token_mask.nii.gz',
        '--vascular_mask', f'{tmppath}/inputs/token_mask.nii.gz', 
        '--labels', f'{tmppath}/inputs/token_mask.nii.gz',
        ]

from rabies.run_main import execute_workflow
execute_workflow(args=args)
