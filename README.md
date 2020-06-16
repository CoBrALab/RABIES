# RABIES: Rodent Automated Bold Improvement of EPI Sequences.

![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/processing_schema.jpg)

## Command Line Interface
```
usage: rabies [-h] [-e] [--apply_despiking] [--apply_slice_mc]
              [-b BIAS_REG_SCRIPT] [-r COREG_SCRIPT] [-p PLUGIN]
              [--local_threads LOCAL_THREADS] [--min_proc MIN_PROC]
              [--data_type DATA_TYPE] [--debug]
              [--nativespace_resampling NATIVESPACE_RESAMPLING]
              [--commonspace_resampling COMMONSPACE_RESAMPLING]
              [--anatomical_resampling ANATOMICAL_RESAMPLING]
              [--cluster_type {local,sge,pbs,slurm}] [--walltime WALLTIME]
              [--memory_request MEMORY_REQUEST]
              [--template_reg_script TEMPLATE_REG_SCRIPT] [--detect_dummy]
              [--no_STC] [--TR TR] [--tpattern TPATTERN]
              [--anat_template ANAT_TEMPLATE] [--brain_mask BRAIN_MASK]
              [--WM_mask WM_MASK] [--CSF_mask CSF_MASK]
              [--vascular_mask VASCULAR_MASK] [--labels LABELS]
              bids_dir output_dir

RABIES performs preprocessing of rodent fMRI images. Can either run on
datasets that only contain EPI images, or both structural and EPI images.
Refer to the README documentation for the input folder structure.

positional arguments:
  bids_dir              the root folder of the BIDS-formated input data
                        directory.
  output_dir            the output path to drop outputs from major
                        preprocessing steps.

optional arguments:
  -h, --help            show this help message and exit
  -e, --bold_only       Apply preprocessing with only EPI scans. commonspace
                        registration is executed through registration of the
                        EPI-generated template from ants_dbm to the anatomical
                        template. (default: False)
  --apply_despiking     Whether to apply despiking of the EPI timeseries based
                        on AFNI's 3dDespike https://afni.nimh.nih.gov/pub/dist
                        /doc/program_help/3dDespike.html. (default: False)
  --apply_slice_mc      Whether to apply a slice-specific motion correction
                        after initial volumetric rigid correction. This second
                        motion correction can correct for interslice
                        misalignment resulting from within-TR motion.With this
                        option, motion corrections and the subsequent
                        resampling from registration are applied
                        sequentially,since the 2D slice registrations cannot
                        be concatenate with 3D transforms. (default: False)
  -b BIAS_REG_SCRIPT, --bias_reg_script BIAS_REG_SCRIPT
                        specify a registration script for iterative bias field
                        correction. This registration step consists of
                        aligning the volume with the commonspace template to
                        provide a brain mask and optimize the bias field
                        correction. (default: Rigid)
  -r COREG_SCRIPT, --coreg_script COREG_SCRIPT
                        Specify EPI to anat coregistration script. Built-in
                        options include 'Rigid', 'Affine', 'SyN' (non-linear)
                        and 'light_SyN', but can specify a custom registration
                        script following the template script structure (see
                        RABIES/rabies/shell_scripts/ for template). (default:
                        light_SyN)
  -p PLUGIN, --plugin PLUGIN
                        Specify the nipype plugin for workflow execution.
                        Consult nipype plugin documentation for detailed
                        options. Linear, MultiProc, SGE and SGEGraph have been
                        tested. (default: Linear)
  --local_threads LOCAL_THREADS
                        For local MultiProc execution, set the maximum number
                        of processors run in parallel, defaults to number of
                        CPUs. This option only applies to the MultiProc
                        execution plugin, otherwise it is set to 1. (default:
                        12)
  --min_proc MIN_PROC   For SGE parallel processing, specify the minimal
                        number of nodes to be assigned. (default: 1)
  --data_type DATA_TYPE
                        Specify data format outputs to control for file size
                        among 'int16','int32','float32' and 'float64'.
                        (default: float32)
  --debug               Run in debug mode. (default: False)

Options for the resampling of the EPI. Axis resampling specifications must follow the format 'dim1xdim2xdim3' (in mm) with the RAS axis convention (dim1=Right-Left, dim2=Anterior-Posterior, dim3=Superior-Inferior).:
  --nativespace_resampling NATIVESPACE_RESAMPLING
                        Can specify a resampling dimension for the nativespace
                        outputs. Must be of the form dim1xdim2xdim3 (in mm).
                        The original dimensions are conserved if 'origin' is
                        specified. (default: origin)
  --commonspace_resampling COMMONSPACE_RESAMPLING
                        Can specify a resampling dimension for the commonspace
                        outputs. Must be of the form dim1xdim2xdim3 (in mm).
                        The original dimensions are conserved if 'origin' is
                        specified.***this option specifies the resampling for
                        the --bold_only workflow (default: origin)
  --anatomical_resampling ANATOMICAL_RESAMPLING
                        To optimize the efficiency of registration, the
                        provided anatomical template is resampled based on the
                        provided input images. The dimension with the lowest
                        resolution among the provided anatomical images (EPI
                        images instead if --bold_only is True) is selected as
                        a basis for resampling the template to isotropic
                        resolution, if the provided resolution is lower than
                        the original resolution of the template.
                        Alternatively, the user can provide a custom
                        resampling dimension. This allows to accelerate
                        registration steps with minimal sampling dimensions.
                        (default: inputs_defined)

cluster options for running ants_dbm (options copied from twolevel_dbm.py)::
  --cluster_type {local,sge,pbs,slurm}
                        Choose the type of cluster system to submit jobs to
                        (default: local)
  --walltime WALLTIME   Option for job submission specifying requested time
                        per pairwise registration. (default: 20:00:00)
  --memory_request MEMORY_REQUEST
                        Option for job submission specifying requested memory
                        per pairwise registration. (default: 8gb)
  --template_reg_script TEMPLATE_REG_SCRIPT
                        Registration script that will be used for registration
                        of the generated template to the provided atlas for
                        masking and labeling. Can choose a predefined
                        registration script among Rigid,Affine,SyN or
                        light_SyN, or provide a custom script. (default:
                        light_SyN)

Options for the generation of EPI reference volume.:
  --detect_dummy        Detect and remove dummy volumes, and generate a
                        reference EPI based on these volumes if detected.
                        (default: False)

Specify Slice Timing Correction info that is fed to AFNI 3dTshift
    (https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html). The STC is applied in the
    anterior-posterior orientation, assuming slices were acquired in this direction.:
  --no_STC              Select this option to ignore the STC step. (default:
                        False)
  --TR TR               Specify repetition time (TR). (default: 1.0s)
  --tpattern TPATTERN   Specify if interleaved or sequential acquisition.
                        'alt' for interleaved, 'seq' for sequential. (default:
                        alt)

Template files.:
  --anat_template ANAT_TEMPLATE
                        Anatomical file for the commonspace template.
                        (default: /home/gabriel/RABIES-0.1.2-dev/template_file
                        s/DSURQE_40micron_average.nii.gz)
  --brain_mask BRAIN_MASK
                        Brain mask for the template. (default: /home/gabriel/R
                        ABIES-0.1.2-dev/template_files/DSURQE_100micron_mask.n
                        ii.gz)
  --WM_mask WM_MASK     White matter mask for the template. (default: /home/ga
                        briel/RABIES-0.1.2-dev/template_files/DSURQE_100micron
                        _eroded_WM_mask.nii.gz)
  --CSF_mask CSF_MASK   CSF mask for the template. (default: /home/gabriel/RAB
                        IES-0.1.2-dev/template_files/DSURQE_100micron_eroded_C
                        SF_mask.nii.gz)
  --vascular_mask VASCULAR_MASK
                        Can provide a mask of major blood vessels for
                        computing confound timeseries. The default mask was
                        generated by applying MELODIC ICA and selecting the
                        resulting component mapping onto major veins.
                        (Grandjean et al. 2020, NeuroImage; Beckmann et al.
                        2005) (default: /home/gabriel/RABIES-0.1.2-dev/templat
                        e_files/vascular_mask.nii.gz)
  --labels LABELS       Atlas file with anatomical labels. (default: /home/gab
                        riel/RABIES-0.1.2-dev/template_files/DSURQE_40micron_l
                        abels.nii.gz)
```

## Execution syntax
```
  rabies bids_inputs/ rabies_outputs/ -e -r light_SyN --template_reg_script light_SyN --TR 1.0s --cluster_type sge -p SGEGraph
```
### Running RABIES interactively within a container
Singularity execution
```sh
    singularity run -B /local_input_folder_path:/nii_inputs:ro \
    -B /local_output_folder_path:/rabies_out \
    /path_to_singularity_image/rabies.sif /nii_inputs /rabies_out \
    --rabies_execution_specifications
```
Docker execution
```sh
    docker run -it --rm \
		-v /local_input_folder_path:/nii_inputs:ro \
		-v /local_output_folder_path:/outputs \
		rabies /nii_inputs /outputs --further_execution_specifications
```

<br/>


# Input data folder structure
Input folder must follow the BIDS structure (https://bids.neuroimaging.io/), and must include the tags 'sub', 'ses' and 'run'. The 'ses' and 'run' tags must
have a numerical specification (i.e., ses-1, ses-2, ...)
* Example BOLD scan format: input_folder/sub-{subject_id}/ses-{session_number}/func/sub-{subject_id}_ses-{session_number}_task-rest_acq-EPI_run-{run_number}_bold.nii.gz
* Example Anatomical scan format: input_folder/sub-{subject_id}/ses-{session_number}/anat/sub-{subject_id}_ses-{session_number}_acq-FLASH_T1w.nii.gz

## Directory Tree of an example input folder

<br/>

<!DOCTYPE html>
<html>
<head>
 <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
 <meta name="Author" content="Made by 'tree'">
 <meta name="GENERATOR" content="$Version: $ tree v1.7.0 (c) 1996 - 2014 by Steve Baker, Thomas Moore, Francesc Rocher, Florian Sesser, Kyosuke Tokoro $">
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
	<h1>Directory Tree</h1><p>
	<a href="bids_input">bids_input</a><br>
	├── <a href="bids_input/sub-MFC067/">sub-MFC067</a><br>
	│   └── <a href="bids_input/sub-MFC067/ses-1/">ses-1</a><br>
	│   &nbsp;&nbsp;&nbsp; ├── <a href="bids_input/sub-MFC067/ses-1/anat/">anat</a><br>
	│   &nbsp;&nbsp;&nbsp; │   └── <a href="bids_input/sub-MFC067/ses-1/anat/sub-MFC067_ses-1_acq-FLASH_T1w.nii.gz">sub-MFC067_ses-1_acq-FLASH_T1w.nii.gz</a><br>
	│   &nbsp;&nbsp;&nbsp; └── <a href="bids_input/sub-MFC067/ses-1/func/">func</a><br>
	│   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="bids_input/sub-MFC067/ses-1/func/sub-MFC067_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz">sub-MFC067_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz</a><br>
	└── <a href="bids_input/sub-MFC068/">sub-MFC068</a><br>
	&nbsp;&nbsp;&nbsp; └── <a href="bids_input/sub-MFC068/ses-1/">ses-1</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="bids_input/sub-MFC068/ses-1/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   └── <a href="bids_input/sub-MFC068/ses-1/anat/sub-MFC068_ses-1_acq-FLASH_T1w.nii.gz">sub-MFC068_ses-1_acq-FLASH_T1w.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="bids_input/sub-MFC068/ses-1/func/">func</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="bids_input/sub-MFC068/ses-1/func/sub-MFC068_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz">sub-MFC068_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz</a><br>
	<br><br>
	</p>
	<p>

8 directories, 4 files
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

# Managing outputs
Important outputs will be found in datasink folders:
* anat_datasink: Includes outputs specific to the anatomical workflow
* bold_datasink: Includes corrected EPI timeseries (corrected_bold/ for native space and commonspace_bold/ for registered to commonspace), EPI masks and key EPI outputs from the preprocessing workflow
* commonspace_datasink: Outputs from the common space registration
* transforms_datasink: Contains all transforms
* confounds_datasink: contains outputs relevant for confound regression

## Recommendations for Quality Control (QC)
Registration overlaps and motion timecourses are presented in .png format in the rabies_out/QC_report directory:
* motion_trace: timecourses of the 6 rigid body parameters
* EPI2Anat: registration of the EPI to the anatomical image within subject
* Anat2Template: registration of the anatomical image to the dataset-generated template
* Template2Commonspace: registration of the dataset template to the provided commonspace template
The following image presents an example of the overlap for the EPI2Anat registration:
![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/sub-jgrAesMEDISOc11L_ses-1_run-1_EPI2Anat.png)
<br/>
For direct investigation of the output .nii files relevant for QC:
<br/>
* bias correction: can visualize if bias correction was correctly applied to correct intensity inhomogeneities for the anatomical scan (anat_datasink/anat_preproc/) and EPI reference image (bold_datasink/bias_cor_bold/)
* commonspace registration: verify that each anatomical image (commonspace_datasink/ants_dbm_outputs/ants_dbm/output/secondlevel/secondlevel_template0sub-*_ses-*_preproc0WarpedToTemplate.nii.gz) was properly realigned to the dataset-generated template (commonspace_datasink/ants_dbm_template/secondlevel_template0.nii.gz)
* template registration: verify that the dataset-generated template (commonspace_datasink/warped_template/secondlevel_template0_output_warped_image.nii.gz) was realigned properly to the provided commonspace template (--anat_template input)
* EPI_Coregistration: verify for each session that the bias field-corrected reference EPI (bold_datasink/bias_cor_bold_warped2anat/) was appropriately registered to the anatomical scan of that session (anat_datasink/anat_preproc/)

<br/>

## Detailed graph for the nipype workflow structure:
### main_wf:
![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/graph_orig.png)
### bold_main_wf:
![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/graph_bold_workflow.png)

## Publications
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. Recurrent functional connectivity gradients identified along specific frequency bands of oscillatory coherence and across anesthesia protocols for mouse fMRI. Presented at Society for Neuroscience 2019, Chicago, IL
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. (2019) Dynamic functional connectivity properties are differentially affected by anesthesia protocols and compare across species. Presented at Organization for Human Brain Mapping 2019, Rome, Italy
* Gabriel Desrosiers-Gregoire, Daniel Gallino, Gabriel A. Devenyi, M. Mallar Chakravarty. (2019) Comparison of the BOLD-evoked response to hypercapnic challenge in mice anesthetized under isoflurane and dexmedetomidine. Presented at International Society for Magnetic Resonance in Medicine 2019, Montreal, QC

## References

* fmriprep - https://github.com/poldracklab/fmriprep
