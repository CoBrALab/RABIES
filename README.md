# RABIES: Rodent Automated Bold Improvement of EPI Sequences.

![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/processing_schema.png)

# Installation

The software is available on PyPi https://pypi.org/project/rabies/.
You can install it on your local python environment with ```pip install rabies```.

For general uses, we recommend instead using a containerized installation with Singularity (https://singularity.lbl.gov) or Docker (https://www.docker.com),
which will contain all additional dependencies from *dependencies.txt* and ensure reproducible behavior.
See https://hub.docker.com/repository/docker/gabdesgreg/rabies/general for available versions.

* Install Singularity .sif image ```singularity build rabies.sif docker://gabdesgreg/rabies:tagname```
* Install Docker image: ```docker pull gabdesgreg/rabies:tagname```

# Command Line Interface (CLI)

The following section describes the --help outputs from the command line interface

## General CLI
<details><summary><b>Click to expand</b></summary>
<p>

```
usage: rabies [-h]
              [-p {Linear,MultiProc,SGE,SGEGraph,PBS,LSF,SLURM,SLURMGraph}]
              [--local_threads LOCAL_THREADS]
              [--scale_min_memory SCALE_MIN_MEMORY] [--min_proc MIN_PROC]
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
                        atlas and associated masks, evaluation of confounding
                        timecourses, and includes various execution options
                        (see --help).
    confound_regression
                        Different options for confound regression are
                        available to apply directly on preprocessing outputs
                        from RABIES. Detrending of timeseries is always
                        applied. Otherwise only selected confound regression
                        and denoising strategies are applied. The denoising
                        steps are applied in the following order: detrending
                        first then ICA-AROMA, then highpass/lowpass filtering
                        followed by temporal censoring, then regression of
                        confound timeseries, standardization of timeseries,
                        and finally smoothing.
    analysis            A few built-in resting-state functional connectivity
                        (FC) analysis options are provided to conduct rapid
                        analysis on the cleaned timeseries. The options
                        include seed-based FC, voxelwise or parcellated whole-
                        brain FC, group-ICA and dual regression.
    data_diagnosis      This workflow regroups a set of tools which allow to
                        establish the influence of confounding sources on FC
                        measures. We recommend to conduct a data diagnosis
                        from this workflow to complement FC analysis by
                        evaluating the effectiveness of the confound
                        correction strategies and ensure reliable subsequent
                        interpretations.

Options for managing the execution of the workflow.:
  -p {Linear,MultiProc,SGE,SGEGraph,PBS,LSF,SLURM,SLURMGraph}, --plugin {Linear,MultiProc,SGE,SGEGraph,PBS,LSF,SLURM,SLURMGraph}
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
```
</p>
</details>

## Preprocessing CLI
<details><summary><b>Click to expand</b></summary>
<p>

```
usage: rabies preprocess [-h] [-e]
                         [--bold_bias_cor_method {otsu_reg,thresh_reg,mouse-preprocessing-v5.sh}]
                         [--anat_bias_cor_method {otsu_reg,thresh_reg,mouse-preprocessing-v5.sh}]
                         [--disable_anat_preproc] [--apply_despiking]
                         [--HMC_option {intraSubjectBOLD,0,1,2,3}]
                         [--HMC_transform {Rigid,Affine}] [--apply_slice_mc]
                         [--detect_dummy]
                         [--data_type {int16,int32,float32,float64}] [--debug]
                         [--coreg_script COREG_SCRIPT]
                         [--template_reg_script TEMPLATE_REG_SCRIPT]
                         [--fast_commonspace]
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
                        template. (default: False)
  --bold_bias_cor_method {otsu_reg,thresh_reg,mouse-preprocessing-v5.sh}
                        Choose the algorithm for bias field correction of the
                        EPI before registration.otsu_reg will conduct an
                        initial serie of N4BiasFieldCorrection oriented by
                        Otsu masking method, followed by a rigid registration
                        to provide a brain mask orienting the final
                        correction.thresh_reg will instead use an initial
                        voxel intensity thresholding method for masking, and
                        will conduct a subsequent rigid registration to
                        provide a brain mask orienting the final correction.
                        (default: mouse-preprocessing-v5.sh)
  --anat_bias_cor_method {otsu_reg,thresh_reg,mouse-preprocessing-v5.sh}
                        Same as --bold_bias_cor_method but for the anatomical
                        image. (default: mouse-preprocessing-v5.sh)
  --disable_anat_preproc
                        This option disables the preprocessing of anatomical
                        images before commonspace template generation.
                        (default: False)
  --apply_despiking     Whether to apply despiking of the EPI timeseries based
                        on AFNI's 3dDespike https://afni.nimh.nih.gov/pub/dist
                        /doc/program_help/3dDespike.html. (default: False)
  --HMC_option {intraSubjectBOLD,0,1,2,3}
                        Select an option for head motion realignment among the
                        pre-built options from https://github.com/ANTsX/ANTsR/
                        blob/60eefd96fedd16bceb4703ecd2cd5730e6843807/R/ants_m
                        otion_estimation.R. (default: intraSubjectBOLD)
  --HMC_transform {Rigid,Affine}
                        Select between Rigid and Affine registration for head
                        motion realignment. (default: Rigid)
  --apply_slice_mc      Whether to apply a slice-specific motion correction
                        after initial volumetric rigid correction. This second
                        motion correction can correct for interslice
                        misalignment resulting from within-TR motion.With this
                        option, motion corrections and the subsequent
                        resampling from registration are applied
                        sequentially,since the 2D slice registrations cannot
                        be concatenate with 3D transforms. (default: False)
  --detect_dummy        Detect and remove initial dummy volumes from the EPI,
                        and generate a reference EPI based on these volumes if
                        detected.Dummy volumes will be removed from the output
                        preprocessed EPI. (default: False)
  --data_type {int16,int32,float32,float64}
                        Specify data format outputs to control for file size.
                        (default: float32)
  --debug               Run in debug mode. (default: False)

Options for the registration steps. Built-in options for selecting registration scripts include 'Rigid', 'Affine', 'SyN' (non-linear), 'light_SyN', 'heavy_SyN', 'multiRAT', but can specify a custom registration script following the template script structure (see RABIES/rabies/shell_scripts/ for template).'Rigid', 'Affine' and 'SyN' options rely on an adaptive registration framework which adapts parameters to the images dimensions:
  --coreg_script COREG_SCRIPT
                        Specify EPI to anat coregistration script. (default:
                        SyN)
  --template_reg_script TEMPLATE_REG_SCRIPT
                        Registration script that will be used for registration
                        of the generated dataset template to the provided
                        commonspace atlas for masking and labeling. (default:
                        SyN)
  --fast_commonspace    Choosing this option will skip the generation of a
                        dataset template, and instead, each anatomical scan
                        will be individually registered to the commonspace
                        template using the --template_reg_script.Note that
                        this option, although faster, is expected to reduce
                        the quality of commonspace registration. (default:
                        False)

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

Specify Slice Timing Correction info that is fed to AFNI 3dTshift
    (https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html). The STC is applied in the
    anterior-posterior orientation, assuming slices were acquired in this direction.:
  --TR TR               Specify repetition time (TR) in seconds. (default:
                        1.0s)
  --no_STC              Select this option to ignore the STC step. (default:
                        False)
  --tpattern {alt,seq}  Specify if interleaved or sequential acquisition.
                        'alt' for interleaved, 'seq' for sequential. (default:
                        alt)

Provided commonspace atlas files.:
  --anat_template ANAT_TEMPLATE
                        Anatomical file for the commonspace template.
                        (default: /home/gabriel/.local/share/rabies/DSURQE_40m
                        icron_average.nii.gz)
  --brain_mask BRAIN_MASK
                        Brain mask for the template. (default: /home/gabriel/.
                        local/share/rabies/DSURQE_40micron_mask.nii.gz)
  --WM_mask WM_MASK     White matter mask for the template. (default: /home/ga
                        briel/.local/share/rabies/DSURQE_40micron_eroded_WM_ma
                        sk.nii.gz)
  --CSF_mask CSF_MASK   CSF mask for the template. (default: /home/gabriel/.lo
                        cal/share/rabies/DSURQE_40micron_eroded_CSF_mask.nii.g
                        z)
  --vascular_mask VASCULAR_MASK
                        Can provide a mask of major blood vessels for
                        computing confound timeseries. The default mask was
                        generated by applying MELODIC ICA and selecting the
                        resulting component mapping onto major veins.
                        (Grandjean et al. 2020, NeuroImage; Beckmann et al.
                        2005) (default: /home/gabriel/.local/share/rabies/vasc
                        ular_mask.nii.gz)
  --labels LABELS       Atlas file with anatomical labels. (default: /home/gab
                        riel/.local/share/rabies/DSURQE_40micron_labels.nii.gz
                        )
```

</p>
</details>

## Confound Regression CLI
<details><summary><b>Click to expand</b></summary>
<p>

```
usage: rabies confound_regression [-h] [--wf_name WF_NAME]
                                  [--commonspace_bold] [--TR TR]
                                  [--highpass HIGHPASS] [--lowpass LOWPASS]
                                  [--smoothing_filter SMOOTHING_FILTER]
                                  [--run_aroma] [--aroma_dim AROMA_DIM]
                                  [--conf_list [{WM_signal,CSF_signal,vascular_signal,global_signal,aCompCor,mot_6,mot_24,mean_FD} [{WM_signal,CSF_signal,vascular_signal,global_signal,aCompCor,mot_6,mot_24,mean_FD} ...]]]
                                  [--FD_censoring]
                                  [--FD_threshold FD_THRESHOLD]
                                  [--DVARS_censoring] [--standardize]
                                  [--timeseries_interval TIMESERIES_INTERVAL]
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
                        strategies). (default: confound_regression_wf)
  --commonspace_bold    If should run confound regression on the commonspace
                        bold output. (default: False)
  --TR TR               Specify repetition time (TR) in seconds. (default:
                        1.0s)
  --highpass HIGHPASS   Specify highpass filter frequency. (default: None)
  --lowpass LOWPASS     Specify lowpass filter frequency. (default: None)
  --smoothing_filter SMOOTHING_FILTER
                        Specify smoothing filter size in mm. (default: None)
  --run_aroma           Whether to run ICA-AROMA or not. The classifier
                        implemented within RABIES is a slightly modified
                        version from the original (Pruim et al. 2015), with
                        parameters and masks adapted for rodent images.
                        (default: False)
  --aroma_dim AROMA_DIM
                        Can specify a number of dimension for the MELODIC run
                        before ICA-AROMA. (default: 0)
  --conf_list [{WM_signal,CSF_signal,vascular_signal,global_signal,aCompCor,mot_6,mot_24,mean_FD} [{WM_signal,CSF_signal,vascular_signal,global_signal,aCompCor,mot_6,mot_24,mean_FD} ...]]
                        list of nuisance regressors that will be applied on
                        voxel timeseries. mot_6 corresponds to the 6 rigid
                        body parameters, and mot_24 corresponds to the 6 rigid
                        parameters, their temporal derivative, and all 12
                        parameters squared (Friston et al. 1996). aCompCor
                        corresponds the timeseries of components from a PCA
                        conducted on the combined WM and CSF masks voxel
                        timeseries, including all components that together
                        explain 50 percent. of the variance, as in Muschelli
                        et al. 2014. (default: [])
  --FD_censoring        Whether to remove timepoints that exceed a framewise
                        displacement threshold, before other confound
                        regression steps. The frames that exceed the given
                        threshold together with 1 back and 4 forward frames
                        will be masked out (based on Power et al. 2012).
                        (default: False)
  --FD_threshold FD_THRESHOLD
                        Scrubbing threshold for the mean framewise
                        displacement in mm (averaged across the brain mask) to
                        select corrupted volumes. (default: 0.05)
  --DVARS_censoring     Whether to remove timepoints that present outlier
                        values on the DVARS metric (temporal derivative of
                        global signal). This option will create a distribution
                        of DVARS values which has no outlier values above or
                        below 2.5 standard deviations. (default: False)
  --standardize         Whether to standardize timeseries to unit variance.
                        (default: False)
  --timeseries_interval TIMESERIES_INTERVAL
                        Specify a time interval in the timeseries to keep.
                        e.g. "0,80". By default all timeseries are kept.
                        (default: all)
```

</p>
</details>

## Analysis CLI
<details><summary><b>Click to expand</b></summary>
<p>

```
usage: rabies analysis [-h] [--scan_list [SCAN_LIST [SCAN_LIST ...]]]
                       [--seed_list [SEED_LIST [SEED_LIST ...]]] [--FC_matrix]
                       [--ROI_type {parcellated,voxelwise}] [--group_ICA]
                       [--TR TR] [--dim DIM] [--DR_ICA] [--IC_file IC_FILE]
                       confound_regression_out output_dir

positional arguments:
  confound_regression_out
                        path to RABIES confound regression output directory
                        with the datasink.
  output_dir            the output path to drop analysis outputs.

optional arguments:
  -h, --help            show this help message and exit
  --scan_list [SCAN_LIST [SCAN_LIST ...]]
                        This option offers to run the analysis on a subset of
                        the scans.The scans selected are specified by
                        providing the full path to each EPI file from the
                        input BIDS folder.The list of scan can be specified
                        manually as a list of file name '--scan_list
                        scan1.nii.gz scan2.nii.gz ...' or the files can be
                        imbedded into a .txt file with one filename per row.By
                        default, 'all' will use all the scans previously
                        processed. (default: ['all'])
  --seed_list [SEED_LIST [SEED_LIST ...]]
                        Can provide a list of seed .nii images that will be
                        used to evaluate seed-based correlation maps based on
                        Pearson's r.Each seed must consist of a binary mask
                        representing the ROI in commonspace. (default: [])

Options for performing a whole-brain timeseries correlation matrix analysis.:
  --FC_matrix           Choose this option to derive a whole-brain functional
                        connectivity matrix, based on the Pearson's r
                        correlation of regional timeseries for each subject
                        cleaned timeseries. (default: False)
  --ROI_type {parcellated,voxelwise}
                        Define the types of ROI to extract regional timeseries
                        for correlation matrix analysis. Options are
                        'parcellated', in which case the atlas labels provided
                        for preprocessing are used as ROIs, or 'voxelwise', in
                        which case all voxel timeseries are cross-correlated.
                        (default: parcellated)

Options for performing group-ICA using FSL's MELODIC on the whole dataset cleaned timeseries.Note that confound regression must have been conducted on commonspace outputs.:
  --group_ICA           Choose this option to conduct group-ICA. (default:
                        False)
  --TR TR               Specify repetition time (TR) in seconds. (default:
                        1.0s)
  --dim DIM             You can specify the number of ICA components to be
                        derived. The default uses an automatic estimation.
                        (default: 0)

Options for performing a dual regression analysis based on a previous group-ICA run from FSL's MELODIC. Note that confound regression must have been conducted on commonspace outputs.:
  --DR_ICA              Choose this option to conduct dual regression on each
                        subject timeseries. This analysis will output the
                        spatial maps corresponding to the linear coefficients
                        from the second linear regression. See
                        rabies.analysis_pkg.analysis_functions.dual_regression
                        for the specific code. (default: False)
  --IC_file IC_FILE     Option to provide a melodic_IC.nii.gz file with the
                        ICA components from a previous group-ICA run. If none
                        is provided, a group-ICA will be run with the dataset
                        cleaned timeseries. (default: None)
```

</p>
</details>

## Data Diagnosis CLI
<details><summary><b>Click to expand</b></summary>
<p>

```
usage: rabies data_diagnosis [-h] [--scan_list [SCAN_LIST [SCAN_LIST ...]]]
                             [--dual_regression]
                             [--dual_convergence DUAL_CONVERGENCE]
                             [--prior_maps PRIOR_MAPS]
                             [--prior_bold_idx [PRIOR_BOLD_IDX [PRIOR_BOLD_IDX ...]]]
                             [--prior_confound_idx [PRIOR_CONFOUND_IDX [PRIOR_CONFOUND_IDX ...]]]
                             [--DSURQE_regions]
                             confound_regression_out output_dir

positional arguments:
  confound_regression_out
                        path to RABIES confound regression output directory
                        with the datasink.
  output_dir            the output path to drop data_diagnosis outputs.

optional arguments:
  -h, --help            show this help message and exit
  --scan_list [SCAN_LIST [SCAN_LIST ...]]
                        This option offers to run the analysis on a subset of
                        the scans.The scans selected are specified by
                        providing the full path to each EPI file from the
                        input BIDS folder.The list of scan can be specified
                        manually as a list of file name '--scan_list
                        scan1.nii.gz scan2.nii.gz ...' or the files can be
                        imbedded into a .txt file with one filename per row.By
                        default, 'all' will use all the scans previously
                        processed. (default: ['all'])
  --dual_regression     Whether to evaluate dual regression outputs. (default:
                        False)
  --dual_convergence DUAL_CONVERGENCE
                        Can specify a number of components to compute using a
                        dual convergence framework. (default: 0)
  --prior_maps PRIOR_MAPS
                        Provide a 4D nifti image with a series of spatial
                        priors representing common sources of signal (e.g. ICA
                        components from a group-ICA run).The default file
                        corresponds to a MELODIC run on a combined group of
                        anesthetized-ventilated mice with MEDISO and awake
                        mice.Confound regression consisted of highpass at 0.01
                        Hz, FD censoring at 0.03mm, DVARS censoring, and
                        mot_6,WM_signal,CSF_signal as regressors. (default:
                        /home/gabriel/.local/share/rabies/melodic_IC.nii.gz)
  --prior_bold_idx [PRIOR_BOLD_IDX [PRIOR_BOLD_IDX ...]]
  --prior_confound_idx [PRIOR_CONFOUND_IDX [PRIOR_CONFOUND_IDX ...]]
  --DSURQE_regions      (default: False)
```

</p>
</details>

# Step 1: Organize Input Data

Input folder must follow the BIDS structure (https://bids.neuroimaging.io/). RABIES will iterate through subjects and search for all available functional scans with suffix 'bold' or 'cbv'.
If anatomical scans are used for preprocessing (--bold_only False), each functional scan will be matched to one corresponding anatomical scan with suffix 'T1w' or 'T2w' of the same subject/session.
<br/>
<br/>
Mandatory BIDS specifications are:
* 'sub-{subject ID}' and 'ses-{session ID}' for both functional and anatomical images
* 'bold' or 'cbv' suffix for functional images
* 'T1w' or 'T2w' for anatomical images
* 'run-{run #}' is necessary for functional images if there are multiple scans per session

### Directory tree of the example dataset
* http://doi.org/10.5281/zenodo.3937697 has the following BIDS structure:

<!DOCTYPE html>
<html>
<head>
 <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
 <meta name="Author" content="Made by 'tree'">
 <meta name="GENERATOR" content="$Version: $ tree v1.7.0 (c) 1996 - 2014 by Steve Baker, Thomas Moore, Francesc Rocher, Florian Sesser, Kyosuke Tokoro $">
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
	<p>
	<a href="test_dataset">test_dataset</a><br>
	├── <a href="test_dataset/sub-MFC067/">sub-MFC067</a><br>
	│   └── <a href="test_dataset/sub-MFC067/ses-1/">ses-1</a><br>
	│   &nbsp;&nbsp;&nbsp; ├── <a href="test_dataset/sub-MFC067/ses-1/anat/">anat</a><br>
	│   &nbsp;&nbsp;&nbsp; │   └── <a href="test_dataset/sub-MFC067/ses-1/anat/sub-MFC067_ses-1_acq-FLASH_T1w.nii.gz">sub-MFC067_ses-1_acq-FLASH_T1w.nii.gz</a><br>
	│   &nbsp;&nbsp;&nbsp; └── <a href="test_dataset/sub-MFC067/ses-1/func/">func</a><br>
	│   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="test_dataset/sub-MFC067/ses-1/func/sub-MFC067_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz">sub-MFC067_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz</a><br>
	└── <a href="test_dataset/sub-MFC068/">sub-MFC068</a><br>
	&nbsp;&nbsp;&nbsp; └── <a href="test_dataset/sub-MFC068/ses-1/">ses-1</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="test_dataset/sub-MFC068/ses-1/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   └── <a href="test_dataset/sub-MFC068/ses-1/anat/sub-MFC068_ses-1_acq-FLASH_T1w.nii.gz">sub-MFC068_ses-1_acq-FLASH_T1w.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="test_dataset/sub-MFC068/ses-1/func/">func</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="test_dataset/sub-MFC068/ses-1/func/sub-MFC068_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz">sub-MFC068_ses-1_task-rest_acq-EPI_run-1_bold.nii.gz</a><br>
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

# Step 2: Execute RABIES
The following section describes the basic syntax to run RABIES with an example dataset available here http://doi.org/10.5281/zenodo.3937697

## Option A: How to run with local installation

<details><summary><b>Click to expand</b></summary>
<p>

**preprocess**
```sh
rabies -p MultiProc preprocess test_dataset/ preprocess_outputs/ --TR 1.0s --no_STC
```
First, this will run the minimal preprocessing step on the test dataset and store outputs into preprocess_outputs/ folder. The option -p MultiProc specifies to run the pipeline in parallel according to available local threads.
<br/>

**confound_regression**
```sh
rabies -p MultiProc confound_regression preprocess_outputs/ confound_regression_outputs/ --TR 1.0s --commonspace_bold --smoothing_filter 0.3 --conf_list WM_signal CSF_signal vascular_signal mot_6
```
Next, to conduct the modeling and regression of confounding sources, the confound_regression step can be run with custom options for denoising. In this case, we apply a highpass filtering at 0.01Hz, together with the voxelwise regression of the 6 rigid realignment parameters and the mean WM,CSF and vascular signal which are derived from masks provided along with the anatomical template. Finally, a smoothing filter 0.3mm is applied. We are running this on the commonspace outputs from preprocess (--commonspace_bold), since we will run analysis in commonspace in the next step.
<br/>

**analysis**
```sh
rabies -p MultiProc analysis confound_regression_outputs analysis_outputs/ --TR 1.0s --group_ICA --DR_ICA
```
Finally, RABIES has a few standard analysis options provided, which are specified in the Analysis documentation. In this example, we are going to run group independent component analysis (--group_ICA), using FSL's MELODIC function, followed by a dual regression (--DR_ICA) to back propagate the group components onto individual subjects.

</p>
</details>

## Option B: How to run with a container (Singularity and Docker)

<details><summary><b>Click to expand</b></summary>
<p>

Containers are independent computing environments which have their own dependencies installed to ensure consistent and reliable
execution of the software regardless of the user.
Singularity containers can be exported to remote high-performance computing platforms (e.g. computecanada).
<br/>
The main difference for the execution of a container consists in relating the paths for all relevant directories from the local
computer to the container's internal folders. This is done using -B for Singularity and -v for Docker. See below for examples:
<br/>

### Singularity execution

**preprocess**
```sh
singularity run -B $PWD/test_dataset:/test_dataset:ro \
-B $PWD/preprocess_outputs:/preprocess_outputs/ \
/path_to_singularity_image/rabies.sif -p MultiProc preprocess /test_dataset/ /preprocess_outputs/ --TR 1.0s --no_STC
```
<br/>

**confound_regression**
```sh
singularity run -B $PWD/test_dataset:/test_dataset:ro \
-B $PWD/preprocess_outputs:/preprocess_outputs/ \
-B $PWD/confound_regression_outputs:/confound_regression_outputs/ \
/path_to_singularity_image/rabies.sif -p MultiProc confound_regression /preprocess_outputs/ /confound_regression_outputs/ --TR 1.0s --highpass 0.01 --commonspace_bold --smoothing_filter 0.3 --conf_list WM_signal CSF_signal vascular_signal mot_6
```
Note here that the path to the dataset is still linked to the container with -B, even though it is not explicitely part of the inputs in the confound regression call. This is necessary since the paths used in the preprocess steps are still accessed in the background, and there will be an error if the paths are not kept consistent across processing steps.
<br/>

**analysis**
```sh
singularity run -B $PWD/test_dataset:/test_dataset:ro \
-B $PWD/preprocess_outputs:/preprocess_outputs/ \
-B $PWD/confound_regression_outputs:/confound_regression_outputs/ \
-B $PWD/analysis_outputs:/analysis_outputs/ \
/path_to_singularity_image/rabies.sif -p MultiProc analysis /confound_regression_outputs /analysis_outputs/ --TR 1.0s --group_ICA --DR_ICA
```
<br/>

### Docker execution
```sh
docker run -it --rm \
-v /local_input_folder_path:/nii_inputs:ro \
-v /local_output_folder_path:/outputs \
rabies preprocess /nii_inputs /outputs --further_execution_specifications
```
The docker execution has very similar syntax to the singularity execution, except that -B is replaced by -v, and some further specifications may be needed (e.g. -it, --rm).

</p>
</details>

# Step 3: Understand your outputs

In this section, there is a description for all the output files provided at each processing steps, together with recommandations for quality control.

## Preprocessing Outputs
<details><summary><b>Click to expand</b></summary>
<p>

Important outputs will be found in the datasink folders. All the different preprocessing outputs are found below:
- **anat_datasink**: Includes outputs specific to the anatomical preprocessing workflow
    - anat_preproc: preprocessed anatomical scans that are used for further registrations
    - anat_mask: brain mask in the anatomical native space
    - WM_mask: WM mask in the anatomical native space
    - CSF_mask: CSF mask in the anatomical native space
    - anat_labels: atlas labels in the anatomical native space
- **bold_datasink**: Includes corrected EPI timeseries (corrected_bold/ for native space and commonspace_bold/ for registered to commonspace), EPI masks and other key EPI outputs from the preprocessing workflow
    - input_bold: original raw EPI images used as inputs into the pipeline
    - corrected_bold: EPI timeseries after preprocessing in native space
    - corrected_bold_ref: reference 3D EPI image (temporal median) after correction
    - bold_brain_mask: brain mask in the corrected_bold space
    - bold_WM_mask: WM mask in the corrected_bold space
    - bold_CSF_mask: CSF mask in the corrected_bold space
    - bold_labels: atlas labels in the corrected_bold space
    - commonspace_bold: EPI timeseries after preprocessing in common space
    - commonspace_bold_mask: brain mask in the commonspace_bold space
    - commonspace_bold_WM_mask: WM mask in the commonspace_bold space
    - commonspace_bold_CSF_mask: CSF mask in the commonspace_bold space
    - commonspace_vascular_mask: vascular mask in the commonspace_bold space
    - commonspace_bold_labels: atlas labels in the commonspace_bold space
    - initial_bold_ref: initial reference 3D EPI image that was subsequently used for bias-field correction
    - bias_cor_bold: reference 3D EPI after bias-field correction which is then used for co-registration
    - bias_cor_bold_warped2anat: bias_cor_bold warped to the co-registration target anatomical image
- **commonspace_datasink**: Outputs from the common space registration
    - ants_dbm_template: the dataset template generated from the registration of anatomical images, using two-level ants dbm (https://github.com/CoBrALab/twolevel_ants_dbm), can be found here
    - warped_template: ants_dbm_template warped to the provided common space template after registration
    - ants_dbm_outputs: a complete output from the two-level ants dbm run for the generation of a dataset anatomical template
- **transforms_datasink**: Contains all transforms
    - affine_bold2anat: affine transforms from the EPI co-registration to the anatomical image
    - warp_bold2anat: non-linear transforms from the EPI co-registration to the anatomical image
    - inverse_warp_bold2anat: inverse of the non-linear transforms from the EPI co-registration to the anatomical image
    - anat_to_template_affine: affine transforms from the registration of the anatomical image to ants_dbm_template registration
    - anat_to_template_warp: non-linear transforms from the registration of the anatomical image to ants_dbm_template registration
    - anat_to_template_inverse_warp: inverse of the non-linear transforms from the registration of the anatomical image to ants_dbm_template
    - template_to_common_affine: affine transforms from the registration of the ants_dbm_template to the commonspace template
    - template_to_common_warp: non-linear transforms from the registration of the ants_dbm_template to the commonspace template
    - template_to_common_inverse_warp: inverse of the non-linear transforms from the registration of the ants_dbm_template to the commonspace template

- **confounds_datasink**: contains confounding features from the EPI that are relevant for subsequent confound regression
    - confounds_csv: a .csv file with the diverse potential confound timecourses. Includes up to 24 motion parameters (6 rigid parameters, their temporal derivative, and all 12 parameters squared; Friston et al. 1996), the global signal, the WM mask signal, the CSF mask signal, the vascular mask signal and aCompCor timecourses (Muschelli et al. 2014).
    - FD_csv: a .csv file with the timecourse of the voxelwise mean and maximal framewise displacement (FD) estimations
    - FD_voxelwise: a .nii image which contains FD timecourses for all single voxel
    - pos_voxelwise: a .nii image which contains the relative positioning timecourses for all single voxel

### Recommendations for Quality Control (QC)
Visual QC outputs in .png format are generate for several processing milestones. These can be found in the rabies_out/QC_report folder.
The milestones include:
* **template_files**: displays the overlap of the provided anatomical template with it's associated masks and labels.
* **anat_denoising**: represents the processing steps for the denoising of the anatomical image before main registration steps. It includes 1-the raw image, 2-the initial correction, 3-the overlap of the registered mask used for final correction, 4-final corrected output.
* **bold_denoising**: same as anat_denoising, but for the functional image.
* **EPI2Anat**: registration of the EPI to the anatomical image within subject
* **Anat2Template**: registration of the anatomical image to the dataset-generated template
* **Template2Commonspace**: registration of the dataset template to the provided commonspace template
* **temporal_denoising**: includes the timecourse of the head motion realignment parameters and associated framewise displacement. Also include spatial maps for temporal variability and tSNR.
The following image presents an example of the overlap for the EPI2Anat registration:
![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/sub-jgrAesMEDISOc11L_ses-1_run-1_EPI2Anat.png)

</p>
</details>

## Confound Regression Outputs

<details><summary><b>Click to expand</b></summary>
<p>

Important outputs from confound regression will be found in the confound_regression_datasink present in the provided output folder:
- **confound_regression_datasink**: Includes outputs specific to the anatomical preprocessing workflow
    - cleaned_timeseries: Resulting timeseries after the application of confound regression
    - VE_file: .pkl file which contains a dictionary vectors, where each vector corresponds to the voxelwise the variance explained (VE) from each regressor in the regression model
    - aroma_out: if --run_aroma is selected, the outputs from running ICA-AROMA will be saved, which includes the MELODIC ICA outputs and the component classification results
    - subject_melodic_ICA: if --diagnosis_output is activated, will contain the outputs from MELODIC ICA run on each individual scan
    - tSNR_map: if --diagnosis_output is activated, this will contain the tSNR map for each scan before confound regression

</p>
</details>

## Analysis Outputs

<details><summary><b>Click to expand</b></summary>
<p>

Important outputs from analysis will be found in the analysis_datasink present in the provided output folder:
- **analysis_datasink**: Includes outputs specific to the anatomical preprocessing workflow
    - group_ICA_dir: complete output from MELODIC ICA, which includes a HTML report for visualization
    - group_IC_file: MELODIC ICA output file with the ICA components
    - DR_data_file: dual regression outputs in the form of a .pkl file which contains a 2D numpy array of component number by voxel number
    - DR_nii_file: dual regression outputs in the form of a .nii file which contains all component 3D maps concatenated into a single .nii file, where the component numbers correspond to the provided template ICA file
    - matrix_data_file: .pkl file which contains a 2D numpy array representing the whole-brain correlation matrix. If using parcellation, the row/column ROI indices are in increasing number of the atlas label number
    - matrix_fig: .png file offered for visualization which represent the correlation matrix
    - seed_correlation_maps: nifti files with voxelwise correlation maps for all provided seeds for seed-based FC

</p>
</details>

# Acknowledgments

**Acknowledging RABIES:** We currently ask users to acknowledge the usage of this software by citing the Github page.

## Presentations
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Francesca Mandino, Joanes Grandjean, M. Mallar Chakravarty. The Convergence of Different fMRI Analysis Streams.
Organization for Human Brain Mapping 2020, Virtual Conference (06/2020)
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. Recurrent functional connectivity gradients identified along specific frequency bands of oscillatory coherence and across anesthesia protocols for mouse fMRI. Presented at Society for Neuroscience 2019, Chicago, IL
* Gabriel Desrosiers-Gregoire, Gabriel A. Devenyi, Joanes Grandjean, M. Mallar Chakravarty. (2019) Dynamic functional connectivity properties are differentially affected by anesthesia protocols and compare across species. Presented at Organization for Human Brain Mapping 2019, Rome, Italy
* Gabriel Desrosiers-Gregoire, Daniel Gallino, Gabriel A. Devenyi, M. Mallar Chakravarty. (2019) Comparison of the BOLD-evoked response to hypercapnic challenge in mice anesthetized under isoflurane and dexmedetomidine. Presented at International Society for Magnetic Resonance in Medicine 2019, Montreal, QC

## References and Influential Material
* fMRIPrep - https://github.com/poldracklab/fmriprep - Esteban, O., Markiewicz, C. J., Blair, R. W., Moodie, C. A., Isik, A. I., Erramuzpe, A., ... & Oya, H. (2019). fMRIPrep: a robust preprocessing pipeline for functional MRI. Nature methods, 16(1), 111-116.
* Two-level DBM - https://github.com/CoBrALab/twolevel_ants_dbm
* ICA-AROMA - https://github.com/maartenmennes/ICA-AROMA - Pruim, R. H., Mennes, M., van Rooij, D., Llera, A., Buitelaar, J. K., & Beckmann, C. F. (2015). ICA-AROMA: A robust ICA-based strategy for removing motion artifacts from fMRI data. Neuroimage, 112, 267-277.
* Grandjean, J., Canella, C., Anckaerts, C., Ayrancı, G., Bougacha, S., Bienert, T., ... & Garin, C. M. (2020). Common functional networks in the mouse brain revealed by multi-centre resting-state fMRI analysis. Neuroimage, 205, 116278.
* FSL MELODIC - https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MELODIC - Beckmann, C. F., Mackay, C. E., Filippini, N., & Smith, S. M. (2009). Group comparison of resting-state FMRI data using multi-subject ICA and dual regression. Neuroimage, 47(Suppl 1), S148.
* ANTs - https://github.com/ANTsX/ANTs
* AFNI - https://afni.nimh.nih.gov/
* Friston, K. J., Williams, S., Howard, R., Frackowiak, R. S., & Turner, R. (1996). Movement‐related effects in fMRI time‐series. Magnetic resonance in medicine, 35(3), 346-355.
* Muschelli, J., Nebel, M. B., Caffo, B. S., Barber, A. D., Pekar, J. J., & Mostofsky, S. H. (2014). Reduction of motion-related artifacts in resting state fMRI using aCompCor. Neuroimage, 96, 22-35.
* Hong, X., To, X. V., Teh, I., Soh, J. R., & Chuang, K. H. (2015). Evaluation of EPI distortion correction methods for quantitative MRI of the brain at high magnetic field. Magnetic resonance imaging, 33(9), 1098-1105.
