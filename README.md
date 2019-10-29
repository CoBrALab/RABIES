# RABIES: Rodent Automated Bold Improvement of EPI Sequences.

![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/processing_schema.jpg)

## Command Line Interface
```
usage: rabies [-h] [-e BOLD_ONLY] [-c COMMONSPACE_METHOD] [-b BIAS_REG_SCRIPT]
              [-r COREG_SCRIPT] [-p PLUGIN] [-d DEBUG] [-v VERBOSE]
              [--cluster_type {local,sge,pbs,slurm}] [--walltime WALLTIME]
              [--memory_request MEMORY_REQUEST]
              [--local_threads LOCAL_THREADS] [--STC STC] [--TR TR]
              [--tpattern TPATTERN] [--anat_template ANAT_TEMPLATE]
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
                        options include 'Rigid', 'Affine' and 'SyN' (non-
                        linear), but can specify a custom registration script
                        following the template script structure (see
                        RABIES/rabies/shell_scripts/ for template). (default:
                        SyN)
  -p PLUGIN, --plugin PLUGIN
                        Specify the nipype plugin for workflow execution.
                        Consult nipype plugin documentation for detailed
                        options. Linear, MultiProc, SGE and SGEGraph have been
                        tested. (default: Linear)
  -d DEBUG, --debug DEBUG
                        Run in debug mode. Default=False (default: False)
  -v VERBOSE, --verbose VERBOSE
                        Increase output verbosity. **doesn't do anything for
                        now. (default: False)

cluster options if commonspace method is ants_dbm (taken from twolevel_dbm.py)::
  --cluster_type {local,sge,pbs,slurm}
                        Choose the type of cluster system to submit jobs to
                        (default: local)
  --walltime WALLTIME   Option for job submission specifying requested time
                        per pairwise registration. (default: 20:00:00)
  --memory_request MEMORY_REQUEST
                        Option for job submission specifying requested memory
                        per pairwise registration. (default: 8gb)
  --local_threads LOCAL_THREADS, -j LOCAL_THREADS
                        For local execution, how many subject-wise modelbuilds
                        to run in parallel, defaults to number of CPUs
                        (default: 8)

Specify Slice Timing Correction info that is fed to AFNI 3dTshift.:
  --STC STC             Whether to run STC or not. (default: True)
  --TR TR               Specify repetition time (TR). (default: 1.0s)
  --tpattern TPATTERN   Specify if interleaved or sequential acquisition.
                        'alt' for interleaved, 'seq' for sequential. (default:
                        alt)

Template files.:
  --anat_template ANAT_TEMPLATE
                        Anatomical file for the commonspace template.
                        (default: /home/cic/desgab/RABIES/DSURQE_atlas/nifti/D
                        SURQE_100micron_average.nii.gz)
  --brain_mask BRAIN_MASK
                        Brain mask for the template. (default: /home/cic/desga
                        b/RABIES/DSURQE_atlas/nifti/DSURQE_100micron_mask.nii.
                        gz)
  --WM_mask WM_MASK     White matter mask for the template. (default: /home/ci
                        c/desgab/RABIES/DSURQE_atlas/nifti/DSURQE_100micron_er
                        oded_WM_mask.nii.gz)
  --CSF_mask CSF_MASK   CSF mask for the template. (default: /home/cic/desgab/
                        RABIES/DSURQE_atlas/nifti/DSURQE_100micron_eroded_CSF_
                        mask.nii.gz)
  --labels LABELS       Atlas file with anatomical labels. (default: /home/cic
                        /desgab/RABIES/DSURQE_atlas/nifti/DSURQE_100micron_lab
                        els.nii.gz)
  --csv_labels CSV_LABELS
                        csv file with info on the labels. (default: /home/cic/
                        desgab/RABIES/DSURQE_atlas/DSURQE_40micron_R_mapping.c
                        sv)


```

## Example execution command
```
  rabies -e 0 -c ants_dbm -b default -r SyN -p SGEGraph -d 0 -v 0 nii_inputs/ rabies_outputs/
```

## Input data folder structure
* input_folder/subject_id/ses-#/bold/subject_id_ses-#_run-#_bold.nii.gz
* input_folder/subject_id/ses-#/anat/subject_id_ses-#_anat.nii.gz
* input_folder/data_info.csv (with columns header: group,subject_id,num_session,num_run)


<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
 <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
 <meta name="Author" content="Made by 'tree'">
 <meta name="GENERATOR" content="$Version: $ tree v1.6.0 (c) 1996 - 2011 by Steve Baker, Thomas Moore, Francesc Rocher, Kyosuke Tokoro $">
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
</head>
<body>
	<h1>Example Directory Tree</h1><p>
	<a href="nii_inputs/">nii_inputs/</a><br>
	└── <a href="nii_inputs//nii_inputs/">nii_inputs</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/data_info.csv">data_info.csv</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/">jgrAesAWc11L</a><br>
	&nbsp;&nbsp;&nbsp; │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/">ses-1</a><br>
	&nbsp;&nbsp;&nbsp; │   │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; │   │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/anat/jgrAesAWc11L_ses-1_anat.nii.gz">jgrAesAWc11L_ses-1_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/bold/jgrAesAWc11L_ses-1_run-1_bold.nii.gz">jgrAesAWc11L_ses-1_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/bold/jgrAesAWc11L_ses-1_run-2_bold.nii.gz">jgrAesAWc11L_ses-1_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/bold/jgrAesAWc11L_ses-1_run-3_bold.nii.gz">jgrAesAWc11L_ses-1_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/">ses-2</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/anat/jgrAesAWc11L_ses-2_anat.nii.gz">jgrAesAWc11L_ses-2_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/bold/jgrAesAWc11L_ses-2_run-1_bold.nii.gz">jgrAesAWc11L_ses-2_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/bold/jgrAesAWc11L_ses-2_run-2_bold.nii.gz">jgrAesAWc11L_ses-2_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/bold/jgrAesAWc11L_ses-2_run-3_bold.nii.gz">jgrAesAWc11L_ses-2_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/">jgrAesAWc11R</a><br>
	&nbsp;&nbsp;&nbsp; │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/">ses-1</a><br>
	&nbsp;&nbsp;&nbsp; │   │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; │   │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/anat/jgrAesAWc11R_ses-1_anat.nii.gz">jgrAesAWc11R_ses-1_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/bold/jgrAesAWc11R_ses-1_run-1_bold.nii.gz">jgrAesAWc11R_ses-1_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/bold/jgrAesAWc11R_ses-1_run-2_bold.nii.gz">jgrAesAWc11R_ses-1_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/bold/jgrAesAWc11R_ses-1_run-3_bold.nii.gz">jgrAesAWc11R_ses-1_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/">ses-2</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/anat/jgrAesAWc11R_ses-2_anat.nii.gz">jgrAesAWc11R_ses-2_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/bold/jgrAesAWc11R_ses-2_run-1_bold.nii.gz">jgrAesAWc11R_ses-2_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/bold/jgrAesAWc11R_ses-2_run-2_bold.nii.gz">jgrAesAWc11R_ses-2_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/bold/jgrAesAWc11R_ses-2_run-3_bold.nii.gz">jgrAesAWc11R_ses-2_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/">jgrAesAWc11R1L</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/">ses-1</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/anat/jgrAesAWc11R1L_ses-1_anat.nii.gz">jgrAesAWc11R1L_ses-1_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/bold/jgrAesAWc11R1L_ses-1_run-1_bold.nii.gz">jgrAesAWc11R1L_ses-1_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/bold/jgrAesAWc11R1L_ses-1_run-2_bold.nii.gz">jgrAesAWc11R1L_ses-1_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/bold/jgrAesAWc11R1L_ses-1_run-3_bold.nii.gz">jgrAesAWc11R1L_ses-1_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/">ses-2</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/anat/jgrAesAWc11R1L_ses-2_anat.nii.gz">jgrAesAWc11R1L_ses-2_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/bold/jgrAesAWc11R1L_ses-2_run-1_bold.nii.gz">jgrAesAWc11R1L_ses-2_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/bold/jgrAesAWc11R1L_ses-2_run-2_bold.nii.gz">jgrAesAWc11R1L_ses-2_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/bold/jgrAesAWc11R1L_ses-2_run-3_bold.nii.gz">jgrAesAWc11R1L_ses-2_run-3_bold.nii.gz</a><br>
	<br><br>
	</p>
	<p>

22 directories, 25 files
	<br><br>
	</p>
	<hr>
	<p class="VERSION">
		 tree v1.6.0 © 1996 - 2011 by Steve Baker and Thomas Moore <br>
		 HTML output hacked and copyleft © 1998 by Francesc Rocher <br>
		 Charsets / OS/2 support © 2001 by Kyosuke Tokoro
	</p>
</body>
</html>

<br/>

## Docker container
After installing the container from the Dockerfile, can run RABIES interactively within a docker container
```sh
    docker run -it --rm \
		-v /nii_inputs_path:/nii_inputs:ro \
		-v /outputs_path:/outputs \
		rabies
```
<br/>

## Publications
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. Recurrent functional connectivity gradients identified along specific frequency bands of oscillatory coherence and across anesthesia protocols for mouse fMRI. Presented at Society for Neuroscience 2019, Chicago, IL
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. (2019) Dynamic functional connectivity properties are differentially affected by anesthesia protocols and compare across species. Presented at Organization for Human Brain Mapping 2019, Rome, Italy
* Gabriel Desrosiers-Gregoire, Daniel Gallino, Gabriel A. Devenyi, M. Mallar Chakravarty. (2019) Comparison of the BOLD-evoked response to hypercapnic challenge in mice anesthetized under isoflurane and dexmedetomidine. Presented at International Society for Magnetic Resonance in Medicine 2019, Montreal, QC

## References

* fmriprep - https://github.com/poldracklab/fmriprep
* minc-toolkit v2 - https://github.com/BIC-MNI/minc-toolkit-v2
* minc-stuffs - https://github.com/Mouse-Imaging-Centre/minc-stuffs
* minc2-simple - https://github.com/vfonov/minc2-simple
* pydpiper - https://github.com/Mouse-Imaging-Centre/pydpiper
