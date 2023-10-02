#! /usr/bin/env python

import os
import tempfile
import shutil
import subprocess
from rabies.utils import generate_token_data

import argparse
def get_parser():
    """Build parser object"""
    parser = argparse.ArgumentParser(
        description=
            "Parser to handle testing using token data.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--complete", dest='complete', action='store_true',
        help=
            "Run a complete testing of the pipeline."
        )
    parser.add_argument(
        '--output_dir', action='store', type=str,
        help=
            "Provide an output directory instead of using a temporary directory.\n"
            "This prevents the deletion of outputs.\n"
        )
    parser.add_argument(
        "--custom", dest='custom', type=str, default=None,
        help=
            "Provide a custom command to run, without specifying input/output directories.\n "
            "The path to token data and output folder is appended at the end of --custom, \n"
            "then the command is run, and the script is exited. The command must be provided \n"
            "as a string between ''.\n"
            "If the preprocessing stage is run, the following arguments are automatically \n"
            "provided to ensure compatibility with token data:\n"
            "   --anat_inho_cor method=disable,otsu_thresh=2,multiotsu=false --bold_inho_cor method=disable,otsu_thresh=2,multiotsu=false \ \n"
            "   --anat_template {tmppath}/inputs/sub-token1_T1w.nii.gz --brain_mask {tmppath}/inputs/token_mask.nii.gz  \ \n"
            "   --WM_mask {tmppath}/inputs/token_mask.nii.gz --CSF_mask {tmppath}/inputs/token_mask.nii.gz  \ \n"
            "   --vascular_mask {tmppath}/inputs/token_mask.nii.gz --labels {tmppath}/inputs/token_mask.nii.gz \ \n"
            "   --bold2anat_coreg registration=no_reg,masking=false,brain_extraction=false \ \n"
            "   --commonspace_reg masking=false,brain_extraction=false,fast_commonspace=true,template_registration=no_reg --data_type int16 \n"
        )
    
    return parser

parser = get_parser()
opts = parser.parse_args()

if opts.output_dir is None:
    tmppath = tempfile.mkdtemp()
else:
    tmppath = opts.output_dir

generate_token_data(tmppath, number_scans=3)

if not opts.custom is None:
    minimal_preproc = f"rabies --inclusion_ids {tmppath}/inputs/sub-token1_bold.nii.gz --verbose 1 preprocess {tmppath}/inputs {tmppath}/outputs --anat_inho_cor method=disable,otsu_thresh=2,multiotsu=false --bold_inho_cor method=disable,otsu_thresh=2,multiotsu=false \
        --anat_template {tmppath}/inputs/sub-token1_T1w.nii.gz --brain_mask {tmppath}/inputs/token_mask.nii.gz --WM_mask {tmppath}/inputs/token_mask.nii.gz --CSF_mask {tmppath}/inputs/token_mask.nii.gz --vascular_mask {tmppath}/inputs/token_mask.nii.gz --labels {tmppath}/inputs/token_mask.nii.gz \
        --bold2anat_coreg registration=no_reg,masking=false,brain_extraction=false --commonspace_reg masking=false,brain_extraction=false,fast_commonspace=true,template_registration=no_reg --data_type int16"
    minimal_cc = f"rabies --verbose 1 confound_correction {tmppath}/outputs {tmppath}/outputs"

    import sys
    command = opts.custom
    if 'preprocess' in command:
        command += f" --anat_inho_cor method=disable,otsu_thresh=2,multiotsu=false --bold_inho_cor method=disable,otsu_thresh=2,multiotsu=false \
    --anat_template {tmppath}/inputs/sub-token1_T1w.nii.gz --brain_mask {tmppath}/inputs/token_mask.nii.gz --WM_mask {tmppath}/inputs/token_mask.nii.gz --CSF_mask {tmppath}/inputs/token_mask.nii.gz --vascular_mask {tmppath}/inputs/token_mask.nii.gz --labels {tmppath}/inputs/token_mask.nii.gz \
    --bold2anat_coreg registration=no_reg,masking=false,brain_extraction=false --commonspace_reg masking=false,brain_extraction=false,fast_commonspace=true,template_registration=no_reg --data_type int16"
    if 'confound_correction' in command or 'analysis' in command:
        if not os.path.isfile(f'{tmppath}/outputs/rabies_preprocess_workflow.pkl'):
            # provide preprocess outputs to run cc stage
            process = subprocess.run(
                minimal_preproc,
                check=True,
                shell=True,
                )
        if 'analysis' in command:
            if not os.path.isfile(f'{tmppath}/outputs/rabies_confound_correction_workflow.pkl'):
                # provide cc outputs to run analysis stage
                process = subprocess.run(
                    minimal_cc,
                    check=True,
                    shell=True,
                    )

    command += f" {tmppath}/outputs {tmppath}/outputs"

    process = subprocess.run(
        command,
        check=True,
        shell=True,
        )
    sys.exit()

command = f"rabies --exclusion_ids {tmppath}/inputs/sub-token2_bold.nii.gz {tmppath}/inputs/sub-token3_bold.nii.gz --verbose 1 preprocess {tmppath}/inputs {tmppath}/outputs --anat_inho_cor method=disable,otsu_thresh=2,multiotsu=false --bold_inho_cor method=disable,otsu_thresh=2,multiotsu=false \
    --anat_template {tmppath}/inputs/sub-token1_T1w.nii.gz --brain_mask {tmppath}/inputs/token_mask.nii.gz --WM_mask {tmppath}/inputs/token_mask.nii.gz --CSF_mask {tmppath}/inputs/token_mask.nii.gz --vascular_mask {tmppath}/inputs/token_mask.nii.gz --labels {tmppath}/inputs/token_mask.nii.gz \
    --bold2anat_coreg registration=no_reg,masking=false,brain_extraction=false --commonspace_reg masking=false,brain_extraction=false,fast_commonspace=true,template_registration=no_reg --data_type int16 --bold_only --detect_dummy \
    --tpattern seq-z --apply_STC --voxelwise_motion --isotropic_HMC"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

command = f"rabies --inclusion_ids {tmppath}/inputs/sub-token1_bold.nii.gz --verbose 1 --force preprocess {tmppath}/inputs {tmppath}/outputs --anat_inho_cor method=disable,otsu_thresh=2,multiotsu=false --bold_inho_cor method=disable,otsu_thresh=2,multiotsu=false \
    --anat_template {tmppath}/inputs/sub-token1_T1w.nii.gz --brain_mask {tmppath}/inputs/token_mask.nii.gz --WM_mask {tmppath}/inputs/token_mask.nii.gz --CSF_mask {tmppath}/inputs/token_mask.nii.gz --vascular_mask {tmppath}/inputs/token_mask.nii.gz --labels {tmppath}/inputs/token_mask.nii.gz \
    --bold2anat_coreg registration=no_reg,masking=true,brain_extraction=true --commonspace_reg masking=true,brain_extraction=true,fast_commonspace=true,template_registration=no_reg --data_type int16  \
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

# rerunning confound correction without censoring, which removes all scans
os.remove(f'{tmppath}/outputs/rabies_confound_correction.pkl')
command = f"rabies --verbose 1 confound_correction {tmppath}/outputs {tmppath}/outputs --nativespace_analysis"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

# testing --data_diagnosis in native space
command = f"rabies --verbose 1 analysis {tmppath}/outputs {tmppath}/outputs --data_diagnosis"
process = subprocess.run(
    command,
    check=True,
    shell=True,
    )

if opts.complete:
    ###GROUP LEVEL, RUNNING ALL 3 SCANS
    shutil.rmtree(f'{tmppath}/outputs/')
    command = f"rabies --verbose 1 preprocess {tmppath}/inputs {tmppath}/outputs --anat_inho_cor method=disable,otsu_thresh=2,multiotsu=false --bold_inho_cor method=disable,otsu_thresh=2,multiotsu=false \
        --anat_template {tmppath}/inputs/sub-token1_T1w.nii.gz --brain_mask {tmppath}/inputs/token_mask.nii.gz --WM_mask {tmppath}/inputs/token_mask_half.nii.gz --CSF_mask {tmppath}/inputs/token_mask_half.nii.gz --vascular_mask {tmppath}/inputs/token_mask_half.nii.gz --labels {tmppath}/inputs/token_mask.nii.gz \
        --bold2anat_coreg registration=no_reg,masking=false,brain_extraction=false --commonspace_reg masking=false,brain_extraction=false,fast_commonspace=true,template_registration=no_reg --data_type int16  \
        --HMC_option 0"
    process = subprocess.run(
        command,
        check=True,
        shell=True,
        )

    command = f"rabies --verbose 1 confound_correction --read_datasink {tmppath}/outputs {tmppath}/outputs --conf_list mot_6 --smoothing_filter 0.3"
    process = subprocess.run(
        command,
        check=True,
        shell=True,
        )

    command = f"rabies --force --verbose 1 analysis {tmppath}/outputs {tmppath}/outputs --NPR_temporal_comp 1 --data_diagnosis --DR_ICA  --seed_list {tmppath}/inputs/token_mask_half.nii.gz"
    process = subprocess.run(
        command,
        check=True,
        shell=True,
        )

    # test group ICA
    command = f"rabies --force --verbose 1 analysis {tmppath}/outputs {tmppath}/outputs --group_ica apply=true,dim=0,random_seed=1"
    process = subprocess.run(
        command,
        check=True,
        shell=True,
        )

if opts.output_dir is None:
    shutil.rmtree(tmppath)
