# RABIES: Rodent Automated Bold Improvement of EPI Sequences.

![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/processing_schema.jpg)

## Command Line Interface
```
usage: rabies [-h] [--bids_input BIDS_INPUT] [-e BOLD_ONLY]
              [-c COMMONSPACE_METHOD] [-b BIAS_REG_SCRIPT] [-r COREG_SCRIPT]
              [-p PLUGIN] [--min_proc MIN_PROC] [-d DEBUG] [-v VERBOSE]
              [--isotropic_resampling ISOTROPIC_RESAMPLING]
              [--upsampling UPSAMPLING] [--data_type DATA_TYPE]
              [--cluster_type {local,sge,pbs,slurm}] [--walltime WALLTIME]
              [--memory_request MEMORY_REQUEST]
              [--local_threads LOCAL_THREADS]
              [--template_reg_script TEMPLATE_REG_SCRIPT] [--STC STC]
              [--TR TR] [--tpattern TPATTERN] [--anat_template ANAT_TEMPLATE]
              [--brain_mask BRAIN_MASK] [--WM_mask WM_MASK]
              [--CSF_mask CSF_MASK] [--labels LABELS]
              [--csv_labels CSV_LABELS]
              input_dir output_dir

RABIES performs preprocessing of rodent fMRI images. Can either run on
datasets that only contain EPI images, or both structural and EPI images.
Refer to the README documentation for the input folder structure.

positional arguments:
  input_dir             the root folder of the input data directory.
  output_dir            the output path to drop outputs from major
                        preprocessing steps.

optional arguments:
  -h, --help            show this help message and exit
  --bids_input BIDS_INPUT
                        If the provided input data folder is in the BIDS
                        format to use the BIDS reader.Note that all .nii
                        inputs will be converted to compressed .gz format.
                        (default: False)
  -e BOLD_ONLY, --bold_only BOLD_ONLY
                        preprocessing with only EPI scans. commonspace
                        registration and distortion correction is executed
                        through registration of the EPIs to a common template
                        atlas. (default: False)
  -c COMMONSPACE_METHOD, --commonspace_method COMMONSPACE_METHOD
                        specify either 'pydpiper' or 'ants_dbm' as common
                        space registration method. Pydpiper can only be
                        executed in parallel with SGE or PBS. ***pydpiper
                        option in development (default: ants_dbm)
  -b BIAS_REG_SCRIPT, --bias_reg_script BIAS_REG_SCRIPT
                        specify a registration script for iterative bias field
                        correction. 'default' is a rigid registration.
                        (default: Rigid)
  -r COREG_SCRIPT, --coreg_script COREG_SCRIPT
                        Specify EPI to anat coregistration script. Built-in
                        options include 'Rigid', 'Affine', 'SyN' (non-linear)
                        and 'light_SyN', but can specify a custom registration
                        script following the template script structure (see
                        RABIES/rabies/shell_scripts/ for template). (default:
                        SyN)
  -p PLUGIN, --plugin PLUGIN
                        Specify the nipype plugin for workflow execution.
                        Consult nipype plugin documentation for detailed
                        options. Linear, MultiProc, SGE and SGEGraph have been
                        tested. (default: Linear)
  --min_proc MIN_PROC   For parallel processing, specify the minimal number of
                        nodes to be assigned. (default: 1)
  -d DEBUG, --debug DEBUG
                        Run in debug mode. Default=False (default: False)
  -v VERBOSE, --verbose VERBOSE
                        Increase output verbosity. **doesn't do anything for
                        now. (default: False)

Options for the resampling of the EPI for motion realignment, susceptibility distortion correction and common space resampling::
  --isotropic_resampling ISOTROPIC_RESAMPLING
                        Whether to resample the EPI to an isotropic resolution
                        based on the lowest dimension. (default: False)
  --upsampling UPSAMPLING
                        Can specify a upsampling parameter to increase the EPI
                        resolution upon resampling and minimize information
                        lost from the transform interpolation. (default: 1.0)
  --data_type DATA_TYPE
                        Specify resampling data format to control for file
                        size. Can specify a numpy data type from https://docs.
                        scipy.org/doc/numpy/user/basics.types.html. (default:
                        float64)

cluster options if commonspace method is ants_dbm (taken from twolevel_dbm.py)::
  --cluster_type {local,sge,pbs,slurm}
                        Choose the type of cluster system to submit jobs to
                        (default: local)
  --walltime WALLTIME   Option for job submission specifying requested time
                        per pairwise registration. (default: 20:00:00)
  --memory_request MEMORY_REQUEST
                        Option for job submission specifying requested memory
                        per pairwise registration. (default: 8gb)
  --local_threads LOCAL_THREADS
                        For local execution, how many subject-wise modelbuilds
                        to run in parallel, defaults to number of CPUs
                        (default: 2)
  --template_reg_script TEMPLATE_REG_SCRIPT
                        Registration script that will be used for registration
                        of the generated template to the provided atlas for
                        masking and labeling. Can choose a predefined
                        registration script among Rigid,Affine,SyN or
                        light_SyN, or provide a custom script. (default: SyN)

Specify Slice Timing Correction info that is fed to AFNI 3dTshift.:
  --STC STC             Whether to run STC or not. (default: True)
  --TR TR               Specify repetition time (TR). (default: 1.0s)
  --tpattern TPATTERN   Specify if interleaved or sequential acquisition.
                        'alt' for interleaved, 'seq' for sequential. (default:
                        alt)

Template files.:
  --anat_template ANAT_TEMPLATE
                        Anatomical file for the commonspace template.
                        (default: /home/gabriel/RABIES/DSURQE_atlas/nifti/DSUR
                        QE_100micron_average.nii.gz)
  --brain_mask BRAIN_MASK
                        Brain mask for the template. (default: /home/gabriel/R
                        ABIES/DSURQE_atlas/nifti/DSURQE_100micron_mask.nii.gz)
  --WM_mask WM_MASK     White matter mask for the template. (default: /home/ga
                        briel/RABIES/DSURQE_atlas/nifti/DSURQE_100micron_erode
                        d_WM_mask.nii.gz)
  --CSF_mask CSF_MASK   CSF mask for the template. (default: /home/gabriel/RAB
                        IES/DSURQE_atlas/nifti/DSURQE_100micron_eroded_CSF_mas
                        k.nii.gz)
  --labels LABELS       Atlas file with anatomical labels. (default: /home/gab
                        riel/RABIES/DSURQE_atlas/nifti/DSURQE_100micron_labels
                        .nii.gz)
  --csv_labels CSV_LABELS
                        csv file with info on the labels. (default: /home/gabr
                        iel/RABIES/DSURQE_atlas/DSURQE_40micron_R_mapping.csv)


```

## Example execution command
```
  rabies nii_inputs/ rabies_outputs/ -e 0 -c ants_dbm -b default -r SyN -p SGEGraph -d 0 -v 0
```

# Input data folder structure
Input folder can follow either the BIDS structure (https://bids.neuroimaging.io/) or the following:
* input_folder/sub-subject_id/ses-#/func/sub-subject_id_ses-#_run-#_bold.nii.gz
* input_folder/sub-subject_id/ses-#/anat/sub-subject_id_ses-#_anat.nii.gz
* input_folder/data_info.csv (with columns header: group,subject_id,num_session,num_run)


<!DOCTYPE html>
<html>
<head>
 <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
 <meta name="Author" content="Made by 'tree'">
 <meta name="GENERATOR" content="$Version: $ tree v1.7.0 (c) 1996 - 2014 by Steve Baker, Thomas Moore, Francesc Rocher, Florian Sesser, Kyosuke Tokoro $">
 <title>Directory Tree</title>
 <style type="text/css">
  <!--
  BODY { font-family : ariel, monospace, sans-serif; }
  P { font-weight: normal; font-family : ariel, monospace, sans-serif; color: black; background-color: transparent;}
  B { font-weight: normal; color: black; background-color: transparent;}
  A:visited { font-weight : normal; text-decoration : none; background-color : transparent; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }
  A:link    { font-weight : normal; text-decoration : none; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }
  A:hover   { color : #000000; font-weight : normal; text-decoration : underline; background-color : yellow; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }
  A:active  { color : #000000; font-weight: normal; background-color : transparent; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }
  .VERSION { font-size: small; font-family : arial, sans-serif; }
  .NORM  { color: black;  background-color: transparent;}
  .FIFO  { color: purple; background-color: transparent;}
  .CHAR  { color: yellow; background-color: transparent;}
  .DIR   { color: blue;   background-color: transparent;}
  .BLOCK { color: yellow; background-color: transparent;}
  .LINK  { color: aqua;   background-color: transparent;}
  .SOCK  { color: fuchsia;background-color: transparent;}
  .EXEC  { color: green;  background-color: transparent;}
  -->
 </style>
</head>
<body>
  <h1>Example Directory Tree</h1><p>
	<a href="nii_input">nii_input</a><br>
	├── <a href="nii_input/data_info.csv">data_info.csv</a><br>
	├── <a href="nii_input/sub-jgrAesAWc11L/">sub-jgrAesAWc11L</a><br>
	│   ├── <a href="nii_input/sub-jgrAesAWc11L/ses-1/">ses-1</a><br>
	│   │   ├── <a href="nii_input/sub-jgrAesAWc11L/ses-1/anat/">anat</a><br>
	│   │   │   └── <a href="nii_input/sub-jgrAesAWc11L/ses-1/anat/sub-jgrAesAWc11L_ses-1_anat.nii.gz">sub-jgrAesAWc11L_ses-1_anat.nii.gz</a><br>
	│   │   └── <a href="nii_input/sub-jgrAesAWc11L/ses-1/func/">func</a><br>
	│   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11L/ses-1/func/sub-jgrAesAWc11L_ses-1_run-1_bold.nii.gz">sub-jgrAesAWc11L_ses-1_run-1_bold.nii.gz</a><br>
	│   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11L/ses-1/func/sub-jgrAesAWc11L_ses-1_run-2_bold.nii.gz">sub-jgrAesAWc11L_ses-1_run-2_bold.nii.gz</a><br>
	│   │   &nbsp;&nbsp;&nbsp; └── <a href="nii_input/sub-jgrAesAWc11L/ses-1/func/sub-jgrAesAWc11L_ses-1_run-3_bold.nii.gz">sub-jgrAesAWc11L_ses-1_run-3_bold.nii.gz</a><br>
	│   └── <a href="nii_input/sub-jgrAesAWc11L/ses-2/">ses-2</a><br>
	│   &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11L/ses-2/anat/">anat</a><br>
	│   &nbsp;&nbsp;&nbsp; │   └── <a href="nii_input/sub-jgrAesAWc11L/ses-2/anat/sub-jgrAesAWc11L_ses-2_anat.nii.gz">sub-jgrAesAWc11L_ses-2_anat.nii.gz</a><br>
	│   &nbsp;&nbsp;&nbsp; └── <a href="nii_input/sub-jgrAesAWc11L/ses-2/func/">func</a><br>
	│   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11L/ses-2/func/sub-jgrAesAWc11L_ses-2_run-1_bold.nii.gz">sub-jgrAesAWc11L_ses-2_run-1_bold.nii.gz</a><br>
	│   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11L/ses-2/func/sub-jgrAesAWc11L_ses-2_run-2_bold.nii.gz">sub-jgrAesAWc11L_ses-2_run-2_bold.nii.gz</a><br>
	│   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_input/sub-jgrAesAWc11L/ses-2/func/sub-jgrAesAWc11L_ses-2_run-3_bold.nii.gz">sub-jgrAesAWc11L_ses-2_run-3_bold.nii.gz</a><br>
	└── <a href="nii_input/sub-jgrAesAWc11R/">sub-jgrAesAWc11R</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11R/ses-1/">ses-1</a><br>
	&nbsp;&nbsp;&nbsp; │   ├── <a href="nii_input/sub-jgrAesAWc11R/ses-1/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; │   │   └── <a href="nii_input/sub-jgrAesAWc11R/ses-1/anat/sub-jgrAesAWc11R_ses-1_anat.nii.gz">sub-jgrAesAWc11R_ses-1_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   └── <a href="nii_input/sub-jgrAesAWc11R/ses-1/func/">func</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11R/ses-1/func/sub-jgrAesAWc11R_ses-1_run-1_bold.nii.gz">sub-jgrAesAWc11R_ses-1_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11R/ses-1/func/sub-jgrAesAWc11R_ses-1_run-2_bold.nii.gz">sub-jgrAesAWc11R_ses-1_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; └── <a href="nii_input/sub-jgrAesAWc11R/ses-1/func/sub-jgrAesAWc11R_ses-1_run-3_bold.nii.gz">sub-jgrAesAWc11R_ses-1_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; └── <a href="nii_input/sub-jgrAesAWc11R/ses-2/">ses-2</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11R/ses-2/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   └── <a href="nii_input/sub-jgrAesAWc11R/ses-2/anat/sub-jgrAesAWc11R_ses-2_anat.nii.gz">sub-jgrAesAWc11R_ses-2_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_input/sub-jgrAesAWc11R/ses-2/func/">func</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11R/ses-2/func/sub-jgrAesAWc11R_ses-2_run-1_bold.nii.gz">sub-jgrAesAWc11R_ses-2_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_input/sub-jgrAesAWc11R/ses-2/func/sub-jgrAesAWc11R_ses-2_run-2_bold.nii.gz">sub-jgrAesAWc11R_ses-2_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_input/sub-jgrAesAWc11R/ses-2/func/sub-jgrAesAWc11R_ses-2_run-3_bold.nii.gz">sub-jgrAesAWc11R_ses-2_run-3_bold.nii.gz</a><br>
	<br><br>
	</p>
	<p>

14 directories, 17 files
	<br><br>
	</p>
	<hr>
	<p class="VERSION">
		 tree v1.7.0 © 1996 - 2014 by Steve Baker and Thomas Moore <br>
		 HTML output hacked and copyleft © 1998 by Francesc Rocher <br>
		 JSON output hacked and copyleft © 2014 by Florian Sesser <br>
		 Charsets / OS/2 support © 2001 by Kyosuke Tokoro
	</p>
</body>
</html>

<br/>

## Docker Execution
After installing the container from the Dockerfile, can run RABIES interactively within a docker container
***note that currently the DSURQE atlas download fails within the container due to a required minc to nifti conversion, and thus, custom atlas files must be specified
```sh
    docker run -it --rm \
		-v /nii_inputs_path:/nii_inputs:ro \
		-v /outputs_path:/outputs \
		rabies /nii_inputs /outputs --further_execution_specifications
```
<br/>

# Managing outputs
Important outputs will be found in datasink folders:
* anat_datasink: Includes outputs specific to the anatomical workflow
* bold_datasink: Includes corrected EPI timeseries, EPI masks and key EPI outputs from the preprocessing workflow
* commonspace_datasink: Outputs from the common space registration
* transforms_datasink: Contains all transforms
* confounds_datasink: contains outputs to use for confound regression steps

<br/>

## Publications
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. Recurrent functional connectivity gradients identified along specific frequency bands of oscillatory coherence and across anesthesia protocols for mouse fMRI. Presented at Society for Neuroscience 2019, Chicago, IL
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. (2019) Dynamic functional connectivity properties are differentially affected by anesthesia protocols and compare across species. Presented at Organization for Human Brain Mapping 2019, Rome, Italy
* Gabriel Desrosiers-Gregoire, Daniel Gallino, Gabriel A. Devenyi, M. Mallar Chakravarty. (2019) Comparison of the BOLD-evoked response to hypercapnic challenge in mice anesthetized under isoflurane and dexmedetomidine. Presented at International Society for Magnetic Resonance in Medicine 2019, Montreal, QC

## References

* fmriprep - https://github.com/poldracklab/fmriprep
