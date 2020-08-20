# RABIES: Rodent Automated Bold Improvement of EPI Sequences.

![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/processing_schema.jpg)

## Command Line Interface
```
usage: rabies [-h] [-p PLUGIN] [--local_threads LOCAL_THREADS]
              [--scale_min_memory SCALE_MIN_MEMORY] [--min_proc MIN_PROC]
              [--data_type {int16,int32,float32,float64}] [--debug]
              Processing step ...

RABIES performs processing of rodent fMRI images. Can either run on datasets
that only contain EPI images, or both structural and EPI images.

optional arguments:
  -h, --help            show this help message and exit

Commands:
  The RABIES workflow is seperated into three different processing steps:
  preprocessing, confound regression and analysis. Outputs from the
  preprocessing provides the inputs for the subsequent confound regression,
  and finally analysis.

  Processing step       Description
    preprocess          Conducts preprocessing on an input dataset in BIDS
                        format. Preprocessing includes realignment for motion,
                        correction for susceptibility distortions through non-
                        linear registration, registration to a commonspace
                        atlas and associated masks, as well as further options
                        (see --help).
    confound_regression
                        Flexible options for confound regression applied on
                        preprocessing outputs from RABIES. Smoothing is
                        applied first, followed by ICA-AROMA, detrending, then
                        regression of confound timeseries orthogonal to the
                        application of temporal filters (nilearn.clean_img,
                        Lindquist 2018), standardization of timeseries and
                        finally scrubbing. The corrections follow user
                        specifications.
    analysis            Optional analysis to conduct on cleaned timeseries.

Options for managing the execution of the workflow.:
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
  --scale_min_memory SCALE_MIN_MEMORY
                        For a parallel execution with MultiProc, the minimal
                        memory attributed to nodes can be scaled with this
                        multiplier to avoid memory crashes. (default: 1.0)
  --min_proc MIN_PROC   For SGE parallel processing, specify the minimal
                        number of nodes to be assigned to avoid memory
                        crashes. (default: 1)
  --data_type {int16,int32,float32,float64}
                        Specify data format outputs to control for file size.
                        (default: float32)
  --debug               Run in debug mode. (default: False)

```

# Preprocessing
```
usage: rabies preprocess [-h] [-e] [--disable_anat_preproc]
                         [--apply_despiking] [--apply_slice_mc]
                         [--detect_dummy] [--autoreg] [-r COREG_SCRIPT]
                         [--bias_reg_script BIAS_REG_SCRIPT]
                         [--template_reg_script TEMPLATE_REG_SCRIPT]
                         [--nativespace_resampling NATIVESPACE_RESAMPLING]
                         [--commonspace_resampling COMMONSPACE_RESAMPLING]
                         [--anatomical_resampling ANATOMICAL_RESAMPLING]
                         [--cluster_type {local,sge,pbs,slurm}]
                         [--walltime WALLTIME] [--TR TR] [--no_STC]
                         [--tpattern {alt,seq}]
                         [--anat_template ANAT_TEMPLATE]
                         [--brain_mask BRAIN_MASK] [--WM_mask WM_MASK]
                         [--CSF_mask CSF_MASK] [--vascular_mask VASCULAR_MASK]
                         [--labels LABELS]
                         bids_dir output_dir

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
                        template.
  --disable_anat_preproc
                        This option disables the preprocessing of anatomical
                        images before commonspace template generation.
  --apply_despiking     Whether to apply despiking of the EPI timeseries based
                        on AFNI's 3dDespike https://afni.nimh.nih.gov/pub/dist
                        /doc/program_help/3dDespike.html.
  --apply_slice_mc      Whether to apply a slice-specific motion correction
                        after initial volumetric rigid correction. This second
                        motion correction can correct for interslice
                        misalignment resulting from within-TR motion.With this
                        option, motion corrections and the subsequent
                        resampling from registration are applied
                        sequentially,since the 2D slice registrations cannot
                        be concatenate with 3D transforms.
  --detect_dummy        Detect and remove initial dummy volumes from the EPI,
                        and generate a reference EPI based on these volumes if
                        detected.Dummy volumes will be removed from the output
                        preprocessed EPI.

Options for the registration steps.:
  --autoreg             Choosing this option will conduct an adaptive
                        registration framework which will adjust parameters
                        according to the input images.This option overrides
                        other registration specifications.
  -r COREG_SCRIPT, --coreg_script COREG_SCRIPT
                        Specify EPI to anat coregistration script. Built-in
                        options include 'Rigid', 'Affine', 'autoreg_affine',
                        'autoreg_SyN', 'SyN' (non-linear), 'light_SyN', but
                        can specify a custom registration script following the
                        template script structure (see
                        RABIES/rabies/shell_scripts/ for template).
  --bias_reg_script BIAS_REG_SCRIPT
                        specify a registration script for iterative bias field
                        correction. This registration step consists of
                        aligning the volume with the commonspace template to
                        provide a brain mask and optimize the bias field
                        correction. The registration script options are the
                        same as --coreg_script.
  --template_reg_script TEMPLATE_REG_SCRIPT
                        Registration script that will be used for registration
                        of the generated dataset template to the provided
                        commonspace atlas for masking and labeling. Can choose
                        a predefined registration script among
                        Rigid,Affine,SyN or light_SyN, or provide a custom
                        script.

Options for the resampling of the EPI. Axis resampling specifications must follow the format 'dim1xdim2xdim3' (in mm) with the RAS axis convention (dim1=Right-Left, dim2=Anterior-Posterior, dim3=Superior-Inferior).:
  --nativespace_resampling NATIVESPACE_RESAMPLING
                        Can specify a resampling dimension for the nativespace
                        outputs. Must be of the form dim1xdim2xdim3 (in mm).
                        The original dimensions are conserved if 'origin' is
                        specified.
  --commonspace_resampling COMMONSPACE_RESAMPLING
                        Can specify a resampling dimension for the commonspace
                        outputs. Must be of the form dim1xdim2xdim3 (in mm).
                        The original dimensions are conserved if 'origin' is
                        specified.***this option specifies the resampling for
                        the --bold_only workflow
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

cluster options for running ants_dbm (options copied from twolevel_dbm.py)::
  --cluster_type {local,sge,pbs,slurm}
                        Choose the type of cluster system to submit jobs to
  --walltime WALLTIME   Option for job submission specifying requested time
                        per pairwise registration.

Specify Slice Timing Correction info that is fed to AFNI 3dTshift
    (https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html). The STC is applied in the
    anterior-posterior orientation, assuming slices were acquired in this direction.:
  --TR TR               Specify repetition time (TR) in seconds.
  --no_STC              Select this option to ignore the STC step.
  --tpattern {alt,seq}  Specify if interleaved or sequential acquisition.
                        'alt' for interleaved, 'seq' for sequential.

Provided commonspace atlas files.:
  --anat_template ANAT_TEMPLATE
                        Anatomical file for the commonspace template.
  --brain_mask BRAIN_MASK
                        Brain mask for the template.
  --WM_mask WM_MASK     White matter mask for the template.
  --CSF_mask CSF_MASK   CSF mask for the template.
  --vascular_mask VASCULAR_MASK
                        Can provide a mask of major blood vessels for
                        computing confound timeseries. The default mask was
                        generated by applying MELODIC ICA and selecting the
                        resulting component mapping onto major veins.
                        (Grandjean et al. 2020, NeuroImage; Beckmann et al.
                        2005)
  --labels LABELS       Atlas file with anatomical labels.
```

## Execution syntax
```sh
rabies -p SGEGraph preprocess bids_inputs/ rabies_outputs/ --autoreg--TR 1.0s --cluster_type sge
```
### Running RABIES interactively within a container
Singularity execution
```sh
singularity run -B /local_input_folder_path:/nii_inputs:ro \
-B /local_output_folder_path:/rabies_out \
/path_to_singularity_image/rabies.sif preprocess /nii_inputs /rabies_out \
--rabies_execution_specifications
```
Docker execution
```sh
docker run -it --rm \
-v /local_input_folder_path:/nii_inputs:ro \
-v /local_output_folder_path:/outputs \
rabies preprocess /nii_inputs /outputs --further_execution_specifications
```

<br/>


# Input data format
Input folder must follow the BIDS structure (https://bids.neuroimaging.io/). RABIES will iterate through subjects and search for all available functional scans with suffix 'bold' or 'cbv'.
If anatomical scans are used for preprocessing (--bold_only False), each functional scan will be matched to one corresponding anatomical scan with suffix 'T1w' or 'T2w' of the same subject/session.

## Directory Tree of an example input folder
* An example dataset for testing RABIES is available http://doi.org/10.5281/zenodo.3937697 with the following structure:

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
	<a href="test_dataset">test_dataset</a><br>
	├── <a href="test_dataset/1-s2.0-S1053811912010890-main.pdf">1-s2.0-S1053811912010890-main.pdf</a><br>
	├── <a href="test_dataset/2020-06-14_Statement.pdf">2020-06-14_Statement.pdf</a><br>
	├── <a href="test_dataset/2020-07-12_Statement.pdf">2020-07-12_Statement.pdf</a><br>
	├── <a href="test_dataset/640602122485.pdf">640602122485.pdf</a><br>
	├── <a href="test_dataset/Cahier%20de%20deuil%20pour%20p%C3%A8re.doc">Cahier de deuil pour père.doc</a><br>
	├── <a href="test_dataset/EpicInstaller-10.17.0-8e5decce9d5243acb9dd2b6a10722800.msi">EpicInstaller-10.17.0-8e5decce9d5243acb9dd2b6a10722800.msi</a><br>
	├── <a href="test_dataset/factory_image.iso">factory_image.iso</a><br>
	├── <a href="test_dataset/fmriprep.pdf">fmriprep.pdf</a><br>
	├── <a href="test_dataset/Fusar-Poli%20et%20al.%20-%202013%20-%20The%20Psychosis%20High-Risk%20State%20A%20Comprehensive%20State-of-the-Art%20Review%20HHS%20Public%20Access.pdf">Fusar-Poli et al. - 2013 - The Psychosis High-Risk State A Comprehensive State-of-the-Art Review HHS Public Access.pdf</a><br>
	├── <a href="test_dataset/In_Pursuit_of_the_Unknown_17_Equations_That_Change..._----_(Intro).pdf">In_Pursuit_of_the_Unknown_17_Equations_That_Change..._----_(Intro).pdf</a><br>
	├── <a href="test_dataset/lex_ai_matt_botvinick.mp3">lex_ai_matt_botvinick.mp3</a><br>
	├── <a href="test_dataset/mont_st-hilaire.pdf">mont_st-hilaire.pdf</a><br>
	├── <a href="test_dataset/Paperpile-receipt-07-01-2020.pdf">Paperpile-receipt-07-01-2020.pdf</a><br>
	├── <a href="test_dataset/Papers/">Papers</a><br>
	│   ├── <a href="test_dataset/Papers/1-s2.0-S1053811909011525-main.pdf">1-s2.0-S1053811909011525-main.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/Brunner%20et%20al.%202020.pdf">Brunner et al. 2020.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/Burt%20et%20al.%202020.pdf">Burt et al. 2020.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/Cold%20Spring%20Harb%20Perspect%20Biol-2009-Sansom-a002519.pdf">Cold Spring Harb Perspect Biol-2009-Sansom-a002519.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/fnsys_2010_00008.pdf">fnsys_2010_00008.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/Karl%20Zilles%202018.pdf">Karl Zilles 2018.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/Liu%20et%20al.%202018.pdf">Liu et al. 2018.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/Mandino%20et%20al.%202020.pdf">Mandino et al. 2020.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/Muschelli%20et%20al.%202014.pdf">Muschelli et al. 2014.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/Shirer%20et%20al.%202015.pdf">Shirer et al. 2015.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/Shrout-Fleiss-ICC.pdf">Shrout-Fleiss-ICC.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/temp_Abstract_mchakravarty_20200130060324449.pdf">temp_Abstract_mchakravarty_20200130060324449.pdf</a><br>
	│   ├── <a href="test_dataset/Papers/WuNordahl-Psychopathology.pdf">WuNordahl-Psychopathology.pdf</a><br>
	│   └── <a href="test_dataset/Papers/Yousefi%20et%20al.%202018.pdf">Yousefi et al. 2018.pdf</a><br>
	├── <a href="test_dataset/Pathfinder%202e%20-%20Core%20Rulebook.pdf">Pathfinder 2e - Core Rulebook.pdf</a><br>
	├── <a href="test_dataset/PCA%20tutorial.pdf">PCA tutorial.pdf</a><br>
	├── <a href="test_dataset/PodcastAddict_export_20190426_141903.opml">PodcastAddict_export_20190426_141903.opml</a><br>
	├── <a href="test_dataset/pour%20le%2012%20septembre.pdf">pour le 12 septembre.pdf</a><br>
	├── <a href="test_dataset/rat.zip">rat.zip</a><br>
	├── <a href="test_dataset/showDocument.do.pdf">showDocument.do.pdf</a><br>
	├── <a href="test_dataset/spreadsheet-analysis.csv">spreadsheet-analysis.csv</a><br>
	├── <a href="test_dataset/test_dataset/">test_dataset</a><br>
	│   ├── <a href="test_dataset/test_dataset/sub-MFC067/">sub-MFC067</a><br>
	│   │   └── <a href="test_dataset/test_dataset/sub-MFC067/ses-1/">ses-1</a><br>
	│   │   &nbsp;&nbsp;&nbsp; ├── <a href="test_dataset/test_dataset/sub-MFC067/ses-1/anat/">anat</a><br>
	│   │   &nbsp;&nbsp;&nbsp; │   └── <a href="test_dataset/test_dataset/sub-MFC067/ses-1/anat/sub-MFC067_ses-1_acq-FLASH_T1w.nii.gz">sub-MFC067_ses-1_acq-FLASH_T1w.nii.gz</a><br>
	│   │   &nbsp;&nbsp;&nbsp; └── <a href="test_dataset/test_dataset/sub-MFC067/ses-1/func/">func</a><br>
	│   │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="test_dataset/test_dataset/sub-MFC067/ses-1/func/sub-MFC067_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz">sub-MFC067_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz</a><br>
	│   └── <a href="test_dataset/test_dataset/sub-MFC068/">sub-MFC068</a><br>
	│   &nbsp;&nbsp;&nbsp; └── <a href="test_dataset/test_dataset/sub-MFC068/ses-1/">ses-1</a><br>
	│   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="test_dataset/test_dataset/sub-MFC068/ses-1/anat/">anat</a><br>
	│   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   └── <a href="test_dataset/test_dataset/sub-MFC068/ses-1/anat/sub-MFC068_ses-1_acq-FLASH_T1w.nii.gz">sub-MFC068_ses-1_acq-FLASH_T1w.nii.gz</a><br>
	│   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="test_dataset/test_dataset/sub-MFC068/ses-1/func/">func</a><br>
	│   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="test_dataset/test_dataset/sub-MFC068/ses-1/func/sub-MFC068_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz">sub-MFC068_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz</a><br>
	├── <a href="test_dataset/test_dataset.zip">test_dataset.zip</a><br>
	├── <a href="test_dataset/Thom%20Yorke%20%20Jonny%20Greenwood%20-%20Motion%20Picture%20Soundtrack%20(acoustic,%20april%2021st%201995,%20Netherlands).mp3">Thom Yorke  Jonny Greenwood - Motion Picture Soundtrack (acoustic, april 21st 1995, Netherlands).mp3</a><br>
	└── <a href="test_dataset/TR202277%20Mallar%20Chakravarty%20(1).pdf">TR202277 Mallar Chakravarty (1).pdf</a><br>
	<br><br>
	</p>
	<p>

10 directories, 41 files
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

# Confound Regression
```
usage: rabies confound_regression [-h] [--wf_name WF_NAME]
                                  [--commonspace_bold] [--TR TR]
                                  [--highpass HIGHPASS] [--lowpass LOWPASS]
                                  [--smoothing_filter SMOOTHING_FILTER]
                                  [--run_aroma] [--aroma_dim AROMA_DIM]
                                  [--conf_list [{WM_signal,CSF_signal,vascular_signal,global_signal,aCompCor,mot_6,mot_24,mean_FD} [{WM_signal,CSF_signal,vascular_signal,global_signal,aCompCor,mot_6,mot_24,mean_FD} ...]]]
                                  [--apply_scrubbing]
                                  [--scrubbing_threshold SCRUBBING_THRESHOLD]
                                  [--timeseries_interval TIMESERIES_INTERVAL]
                                  [--diagnosis_output]
                                  [--seed_list [SEED_LIST [SEED_LIST ...]]]
                                  preprocess_out output_dir

positional arguments:
  preprocess_out        path to RABIES preprocessing output directory with the
                        datasinks.
  output_dir            path to drop confound regression output datasink.

optional arguments:
  -h, --help            show this help message and exit
  --wf_name WF_NAME     Can specify a name for the workflow of this confound
                        regression run, to avoid potential overlaps with
                        previous runs (can be useful if investigating multiple
                        strategies).
  --commonspace_bold    If should run confound regression on the commonspace
                        bold output.
  --TR TR               Specify repetition time (TR) in seconds.
  --highpass HIGHPASS   Specify highpass filter frequency.
  --lowpass LOWPASS     Specify lowpass filter frequency.
  --smoothing_filter SMOOTHING_FILTER
                        Specify smoothing filter size in mm.
  --run_aroma           Whether to run ICA AROMA or not.
  --aroma_dim AROMA_DIM
                        Can specify a number of dimension for MELODIC.
  --conf_list [{WM_signal,CSF_signal,vascular_signal,global_signal,aCompCor,mot_6,mot_24,mean_FD} [{WM_signal,CSF_signal,vascular_signal,global_signal,aCompCor,mot_6,mot_24,mean_FD} ...]]
                        list of regressors.
  --apply_scrubbing     Whether to apply scrubbing or not. A temporal mask
                        will be generated based on the FD threshold. The
                        frames that exceed the given threshold together with 1
                        back and 2 forward frames will be masked out from the
                        data after the application of all other confound
                        regression steps (as in Power et al. 2012).
  --scrubbing_threshold SCRUBBING_THRESHOLD
                        Scrubbing threshold for the mean framewise
                        displacement in mm? (averaged across the brain mask)
                        to select corrupted volumes.
  --timeseries_interval TIMESERIES_INTERVAL
                        Specify a time interval in the timeseries to keep.
                        e.g. "0,80". By default all timeseries are kept.
  --diagnosis_output    Run a diagnosis for each image by computing melodic-
                        ICA on the corrected timeseries,and compute a tSNR map
                        from the input uncorrected image.
  --seed_list [SEED_LIST [SEED_LIST ...]]
                        Can provide a list of seed .nii images that will be
                        used to evaluate seed-based correlation maps during
                        data diagnosis.
```

# Analysis
```
usage: rabies analysis [-h] [--FC_matrix] [--ROI_type {parcellated,voxelwise}]
                       [--group_ICA] [--TR TR] [--dim DIM] [--DR_ICA]
                       [--IC_file IC_FILE]
                       confound_regression_out output_dir

positional arguments:
  confound_regression_out
                        path to RABIES confound regression output directory
                        with the datasink.
  output_dir            the output path to drop analysis outputs.

optional arguments:
  -h, --help            show this help message and exit

Options for performing a whole-brain timeseries correlation matrix analysis.:
  --FC_matrix           Choose this option to derive a whole-brain functional
                        connectivity matrix, based on the correlation of
                        regional timeseries for each subject cleaned
                        timeseries.
  --ROI_type {parcellated,voxelwise}
                        Define the types of ROI to extract regional timeseries
                        for correlation matrix analysis. Options are
                        'parcellated', in which case the atlas labels provided
                        for preprocessing are used as ROIs, or 'voxelwise', in
                        which case all voxel timeseries are cross-correlated.

Options for performing group-ICA using FSL's MELODIC on the whole dataset cleaned timeseries.Note that confound regression must have been conducted on commonspace outputs.:
  --group_ICA           Choose this option to conduct group-ICA.
  --TR TR               Specify repetition time (TR) in seconds.
  --dim DIM             You can specify the number of ICA components to be
                        derived. The default uses an automatic estimation.

Options for performing a dual regression analysis based on a previous group-ICA run from FSL's MELODIC. Note that confound regression must have been conducted on commonspace outputs.:
  --DR_ICA              Choose this option to conduct dual regression on each
                        subject cleaned timeseries.
  --IC_file IC_FILE     Option to provide a melodic_IC.nii.gz file with the
                        ICA components from a previous group-ICA run. If none
                        is provided, a group-ICA will be run with the dataset
                        cleaned timeseries.
```

## Publications
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. Recurrent functional connectivity gradients identified along specific frequency bands of oscillatory coherence and across anesthesia protocols for mouse fMRI. Presented at Society for Neuroscience 2019, Chicago, IL
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. (2019) Dynamic functional connectivity properties are differentially affected by anesthesia protocols and compare across species. Presented at Organization for Human Brain Mapping 2019, Rome, Italy
* Gabriel Desrosiers-Gregoire, Daniel Gallino, Gabriel A. Devenyi, M. Mallar Chakravarty. (2019) Comparison of the BOLD-evoked response to hypercapnic challenge in mice anesthetized under isoflurane and dexmedetomidine. Presented at International Society for Magnetic Resonance in Medicine 2019, Montreal, QC

## References

* fmriprep - https://github.com/poldracklab/fmriprep
