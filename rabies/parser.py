import os
import argparse
from pathlib import Path
import pathos.multiprocessing as multiprocessing  # Better multiprocessing

if 'XDG_DATA_HOME' in os.environ.keys():
    rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
else:
    rabies_path = os.environ['HOME']+'/.local/share/rabies'


def get_parser():
    """Build parser object"""
    parser = argparse.ArgumentParser(
        description=
            "RABIES performs multiple stages of rodent fMRI image processing, including preprocessing, \n"
            "confound correction, simple analyses and data quality assessment.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparsers = parser.add_subparsers(
        title='Processing options',
        description=
            "The RABIES workflow is seperated into three main processing stages: preprocessing, \n"
            "confound correction and analysis. Outputs from the preprocessing provide the inputs for\n"
            "the subsequent confound correction, and finally analysis.",
        help='Description',
        dest='rabies_stage',
        metavar='Processing stage')

    preprocess = subparsers.add_parser("preprocess",
        help=
            "\n"
            "Conducts preprocessing on an input dataset in BIDS format. Preprocessing includes \n"
            "motion realignment, susceptibility distortions correction through non-linear \n"
            "registration, alignment to commonspace, anatomical parcellation and evaluation of \n"
            "nuisance timecourses.\n"
            "\n",
        formatter_class=argparse.RawTextHelpFormatter)
    confound_correction = subparsers.add_parser("confound_correction",
        help=
            "\n"
            "Flexible options for confound correction are applied directly on preprocessing outputs\n"
            "from RABIES to derive cleaned timeseries. Various correction strategies, if selected, are\n"
            "applied in the following order, following best practices from human litterature:\n"
            "   #1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)\n"
            "   #2 - Linear detrending of fMRI timeseries and nuisance regressors\n"
            "   #3 - Apply ICA-AROMA.\n"
            "   #4 - If frequency filtering and frame censoring are applied, simulate data in censored\n" 
            "       timepoints using the Lomb-Scargle periodogram, as suggested in Power et al. (2014, \n"
            "       Neuroimage), for both the fMRI timeseries and nuisance regressors prior to filtering.\n"
            "   #5 - As recommended in Lindquist et al. (2019, Human brain mapping), make the nuisance \n"
            "       regressors orthogonal to the temporal frequency filter.\n"
            "   #6 - Apply highpass and/or lowpass filtering on the fMRI timeseries (with simulated \n"
            "       timepoints).\n"
            "   #7 - Re-apply the frame censoring mask onto filtered fMRI timeseries and nuisance \n"
            "       regressors, taking out the simulated timepoints. Edge artefacts from frequency \n"
            "       filtering can also be removed as recommended in Power et al. (2014, Neuroimage).\n"
            "   #8 - Apply confound regression using the selected nuisance regressors (see --conf_list\n" 
            "       options).\n"
            "   #9 - Standardize timeseries\n"
            "   #10 - Apply Gaussian spatial smoothing.\n"
            "\n",
        formatter_class=argparse.RawTextHelpFormatter)
    analysis = subparsers.add_parser("analysis",
        help=
            "\n"
            "Conduct simple resting-state functional connectivity (FC) analysis, or data quality\n" 
            "diagnosis, on cleaned timeseries after confound correction. Analysis options include\n"
            "seed-based FC, whole-brain FC matrix, group-ICA and dual regression. --data_diagnosis\n" 
            "computes features of data quality at the individual scan and group levels, as in \n"
            "Desrosiers-Gregoire et al. (in prep)\n" 
            "\n",
        formatter_class=argparse.RawTextHelpFormatter)

    ####Execution
    g_execution = parser.add_argument_group(
        title='Execution Options', 
        description=
            "Options for parallel execution and memory management."
        )
    g_execution.add_argument(
        "-p", "--plugin", default='Linear',
        choices=['Linear', 'MultiProc', 'SGE', 'SGEGraph',
                'PBS', 'LSF', 'SLURM', 'SLURMGraph'],
        help=
            "Specify the nipype plugin for workflow execution.\n"
            "Consult https://nipype.readthedocs.io/en/0.11.0/users/plugins.html for details.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_execution.add_argument(
        '--local_threads', type=int, default=multiprocessing.cpu_count(),
        help=
            "For --plugin MultiProc, set the maximum number of processors run in parallel.\n"
            "Defaults to number of CPUs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_execution.add_argument(
        "--scale_min_memory", type=float, default=1.0,
        help=
            "For --plugin MultiProc, set the memory scaling factor attributed to nodes during\n"
            "execution. Increase the scaling if memory crashes are reported.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_execution.add_argument(
        "--min_proc", type=int, default=1,
        help=
            "For --plugin SGE/SGEGraph, scale the number of nodes attributed to jobs to\n"
            "avoid memory crashes.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_execution.add_argument(
        "--verbose", type=int, default=1,
        help=
            "Set the verbose level. 0=WARNING, 1=INFO, 2 or above=DEBUG.\n"
            "(default: %(default)s)\n"
            "\n"
        )


    ####Preprocessing
    preprocess.add_argument(
        'bids_dir', action='store', type=Path,
        help=
            "The root folder of the BIDS-formated input data directory.\n"
            "\n"
        )
    preprocess.add_argument(
        'output_dir', action='store', type=Path,
        help=
            "the output path to drop outputs from major preprocessing steps.\n"
            "\n"
        )
    preprocess.add_argument(
        "--bold_only", dest='bold_only', action='store_true',
        help=
            "Apply preprocessing with only EPI scans. Commonspace registration is executed directly using\n"
            "the corrected EPI 3D reference images. The commonspace registration simultaneously applies\n"
            "distortion correction, this option will produce only commonspace outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--anat_autobox', dest='anat_autobox', action='store_true',
        help=
            "Crops out extra space around the brain on the structural image using AFNI's 3dAutobox\n"
            "https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutobox.html.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--bold_autobox', dest='bold_autobox', action='store_true',
        help=
            "Crops out extra space around the brain on the EPI image using AFNI's 3dAutobox\n"
            "https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutobox.html.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--apply_despiking', dest='apply_despiking', action='store_true',
        help=
            "Applies AFNI's 3dDespike https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        "--HMC_option", type=str, default='intraSubjectBOLD',
        choices=['intraSubjectBOLD', '0', '1', '2', '3'],
        help=
            "Select an option for head motion realignment among the pre-built options from\n"
            "https://github.com/ANTsX/ANTsR/blob/master/R/ants_motion_estimation.R.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--apply_slice_mc', dest='apply_slice_mc', action='store_true',
        help=
            "Whether to apply a slice-specific motion correction after initial volumetric HMC. This can \n"
            "correct for interslice misalignment resulting from within-TR motion. With this option, \n"
            "motion corrections and the subsequent resampling from registration are applied sequentially\n"
            "since the 2D slice registrations cannot be concatenate with 3D transforms. \n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--detect_dummy', dest='detect_dummy', action='store_true',
        help=
            "Detect and remove initial dummy volumes from the EPI, and generate a reference EPI based on\n"
            "these volumes if detected. Dummy volumes will be removed from the output preprocessed EPI.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        "--data_type", type=str, default='float32',
        choices=['int16', 'int32', 'float32', 'float64'],
        help=
            "Specify data format outputs to control for file size.\n"
            "(default: %(default)s)\n"
            "\n"
        )

    g_registration = preprocess.add_argument_group(
        title='Registration Options', 
        description=
            "Customize registration operations and troubleshoot registration failures.\n"
            "*** Rigid: conducts only rigid registration.\n"
            "*** Affine: conducts Rigid then Affine registration.\n"
            "*** SyN: conducts Rigid, Affine then non-linear registration.\n"
            "*** no_reg: skip registration.\n"
        )
    g_registration.add_argument(
        "--anat_inho_cor_method", type=str, default='SyN',
        choices=['Rigid', 'Affine', 'SyN', 'no_reg', 'N4_reg', 'disable'],
        help=
            "Select a registration type for masking during inhomogeneity correction of the structural \n"
            "image. \n"
            "*** N4_reg: previous correction script prior to version 0.3.1.\n"
            "*** disable: disables the inhomogeneity correction.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument("--bold_inho_cor_method", type=str, default='Rigid',
        choices=['Rigid', 'Affine', 'SyN', 'no_reg', 'N4_reg', 'disable'],
        help=
            "Select a registration type for masking during inhomogeneity correction of the EPI.\n"
            "*** N4_reg: previous correction script prior to version 0.3.1.\n"
            "*** disable: disables the inhomogeneity correction.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument("--bold_inho_cor_otsu", type=int, default=2,
        help=
            "The inhomogeneity correction script necessitates an initial correction with a Otsu\n"
            "masking strategy (prior to registration of an anatomical mask). This option sets the \n"
            "Otsu threshold level to capture the right intensity distribution.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        '--atlas_reg_script', type=str, default='SyN',
        choices=['Rigid', 'Affine', 'SyN', 'no_reg'],
        help=
            "Specify a registration script for alignment of the dataset-generated unbiased template \n"
            "to the commonspace atlas.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        "--coreg_script", type=str, default='SyN',
        choices=['Rigid', 'Affine', 'SyN', 'no_reg'],
        help=
            "Specify the registration script for cross-modal alignment between the EPI and structural\n"
            "images. This operation is responsible for correcting EPI susceptibility distortions.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        "--commonspace_masking", dest='commonspace_masking', action='store_true',
        help=
            "Combine masks derived from the inhomogeneity correction step to support registration \n"
            "during the generation of the unbiased template, and then during atlas registration. \n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        "--coreg_masking", dest='coreg_masking', action='store_true',
        help=
            "Use the mask from the EPI inhomogeneity correction step to support registration to the\n"
            "structural image.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        "--brain_extraction", dest='brain_extraction', action='store_true',
        help=
            "If using --commonspace_masking and/or --coreg_masking, this option will conduct brain\n"
            "extractions prior to registration based on the initial mask during inhomogeneity\n"
            "correction. This will enhance brain edge-matching, but requires good quality masks.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        "--fast_commonspace", dest='fast_commonspace', action='store_true',
        help=
            "Skip the generation of a dataset-generated unbiased template, and instead, register each\n"
            "anatomical scan independently directly onto the commonspace atlas, using the\n"
            "--atlas_reg_script registration. This option can be faster, but may decrease the quality\n"
            "of alignment between subjects."
            "(default: %(default)s)\n"
            "\n"
        )

    g_resampling = preprocess.add_argument_group(
        title='Resampling Options', 
        description=
            "The following options allow to resample the voxel dimensions for the preprocessed EPIs\n"
            "or for the anatomical images during registration.\n"
            "The resampling syntax must be 'dim1xdim2xdim3' (in mm), follwing the RAS axis convention\n"
            "(dim1=Right-Left, dim2=Anterior-Posterior, dim3=Superior-Inferior). If 'inputs_defined'\n"
            "is provided instead of axis dimensions, the original dimensions are preserved.\n"
        )
    g_resampling.add_argument(
        '--nativespace_resampling', type=str, default='inputs_defined',
        help=
            "Can specify a resampling dimension for the nativespace fMRI outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_resampling.add_argument(
        '--commonspace_resampling', type=str, default='inputs_defined',
        help=
            "Can specify a resampling dimension for the commonspace fMRI outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_resampling.add_argument(
        '--anatomical_resampling', type=str, default='inputs_defined',
        help=
            "\n"
            "This specifies resampling dimensions for the anatomical registration targets. By \n"
            "default, images are resampled to isotropic resolution based on the smallest dimension\n"
            "among the provided anatomical images (EPI images instead if --bold_only is True). \n"
            "Increasing voxel resampling size will increase registration speed at the cost of \n"
            "accuracy.\n"
            "(default: %(default)s)\n"
            "\n"
        )

    g_stc = preprocess.add_argument_group(
        title='STC Options', 
        description=
            "Specify Slice Timing Correction (STC) info that is fed to AFNI's 3dTshift\n"
            "(https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html). The STC is applied\n"
            "in the anterior-posterior orientation, and thus RABIES assumes slices were acquired in\n"
            "this direction.\n"
        )
    g_stc.add_argument(
        '--apply_STC', dest='apply_STC', action='store_true',
        help=
            "Select this option to apply the STC step.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_stc.add_argument(
        '--TR', type=str, default='auto',
        help=
            "Specify repetition time (TR) in seconds. (e.g. --TR 1.2)\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_stc.add_argument(
        '--tpattern', type=str, default='alt',
        choices=['alt', 'seq'],
        help=
            "Specify if interleaved ('alt') or sequential ('seq') acquisition.\n"
            "(default: %(default)s)\n"
            "\n"
        )

    g_atlas = preprocess.add_argument_group(
        title='Template Files', 
        description=
            "Specify commonspace template and associated mask/label files. By default, RABIES\n"
            "provides the mouse DSURQE atlas\n"
            "https://wiki.mouseimaging.ca/display/MICePub/Mouse+Brain+Atlases.\n"
        )
    g_atlas.add_argument(
        '--anat_template', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_average.nii.gz",
        help=
            "Anatomical file for the commonspace atlas.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--brain_mask', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_mask.nii.gz",
        help=
            "Brain mask aligned with the template.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--WM_mask', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_eroded_WM_mask.nii.gz",
        help=
            "White matter mask aligned with the template.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--CSF_mask', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_eroded_CSF_mask.nii.gz",
        help=
            "CSF mask aligned with the template.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--vascular_mask', action='store', type=Path,
        default=f"{rabies_path}/vascular_mask.nii.gz",
        help=
            "Can provide a mask of major blood vessels to compute associated nuisance timeseries.\n"
            "The default mask was generated by applying MELODIC ICA and selecting the resulting \n"
            "component mapping onto major brain vessels.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--labels', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_labels.nii.gz",
        help=
            "Labels file providing the atlas anatomical annotations.\n"
            "(default: %(default)s)\n"
            "\n"
        )



    ####Confound correction
    confound_correction.add_argument(
        'preprocess_out', action='store', type=Path,
        help=
            "path to RABIES preprocessing output directory.\n"
            "\n"
        )
    confound_correction.add_argument(
        'output_dir', action='store', type=Path,
        help=
            "path for confound correction output directory.\n"
            "\n"
        )
    confound_correction.add_argument(
        '--nativespace_analysis', dest='nativespace_analysis', action='store_true',
        help=
            "Conduct confound correction and analysis in native space.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--conf_list', type=str,
        nargs="*",  # 0 or more values expected => creates a list
        default=[],
        choices=["WM_signal", "CSF_signal", "vascular_signal",
                "global_signal", "aCompCor", "mot_6", "mot_24", "mean_FD"],
        help=
            "Select list of nuisance regressors that will be applied on voxel timeseries, i.e., confound\n"
            "regression.\n"
            "*** WM/CSF/vascular/global_signal: correspond to mean signal from WM/CSF/vascular/brain \n"
            "   masks.\n"
            "*** mot_6: 6 rigid head motion correction parameters.\n"
            "*** mot_24: mot_6 + their temporal derivative, then all 12 parameters squared, as in \n"
            "   Friston et al. (1996, Magnetic Resonance in Medicine).\n"
            "*** aCompCor: method from Muschelli et al. (2014, Neuroimage), where component timeseries\n"
            "   are obtained using PCA, conducted on the combined WM and CSF masks voxel timeseries. \n"
            "   Components adding up to 50 percent of the variance are included.\n"
            "*** mean_FD: the mean framewise displacement timecourse.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--FD_censoring', dest='FD_censoring', action='store_true', default=False,
        help=
            "Apply frame censoring based on a framewise displacement threshold (i.e.scrubbing).\n"
            "The frames that exceed the given threshold, together with 1 back and 2 forward frames\n"
            "will be masked out, as in Power et al. (2012, Neuroimage).\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--TR', type=str, default='auto',
        help=
            "Specify repetition time (TR) in seconds. (e.g. --TR 1.2)\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--highpass', type=float, default=None,
        help=
            "Specify highpass filter frequency.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--lowpass', type=float, default=None,
        help=
            "Specify lowpass filter frequency.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--edge_cutoff', type=float, default=0,
        help=
            "Specify the number of seconds to cut at beginning and end of acquisition if applying a\n"
            "frequency filter. Frequency filters generate edge effects at begining and end of the\n" 
            "timeseries. We recommend to cut those timepoints (around 30sec at both end for 0.01Hz \n" 
            "highpass.).\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--smoothing_filter', type=float, default=None,
        help=
            "Specify filter size in mm for spatial smoothing. Will apply nilearn's function \n"
            "https://nilearn.github.io/modules/generated/nilearn.image.smooth_img.html\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--FD_threshold', type=float, default=0.05,
        help=
            "--FD_censoring threshold in mm.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--DVARS_censoring', dest='DVARS_censoring', action='store_true',default=False,
        help=
            "Whether to remove timepoints that present outlier values on the DVARS metric (temporal\n"
            "derivative of global signal). This method will censor timepoints until the distribution\n" 
            "of DVARS values across time does not contain outliers values above or below 2.5 standard\n" 
            "deviations.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--minimum_timepoint', type=int,default=3,
        help=
            "Can set a minimum number of timepoints remaining after frame censoring. If the threshold\n" 
            "is not met, an empty file is generated and the scan is not considered in further steps.\n" 
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--standardize', dest='standardize', action='store_true',default=False,
        help=
            "Whether to standardize timeseries (z-scoring).\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--run_aroma', dest='run_aroma', action='store_true', default=False,
        help=
            "Whether to run ICA-AROMA or not. The original classifier (Pruim et al. 2015) was modified\n" 
            "to incorporate rodent-adapted masks and classification hyperparameters.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--aroma_dim', type=int, default=0,
        help=
            "Specify a pre-determined number of MELODIC components to derive for ICA-AROMA.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--aroma_random_seed', type=int, default=1,
        help=
            "For reproducibility, this option sets a fixed random seed for MELODIC.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--read_datasink', dest='read_datasink', action='store_true', default=False,
        help=
            "\n"
            "Choose this option to read preprocessing outputs from datasinks instead of the saved \n"
            "preprocessing workflow graph. This allows to run confound correction without having \n"
            "available RABIES preprocessing folders, but the targetted datasink folders must follow the\n"
            "structure of RABIES preprocessing.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--timeseries_interval', type=str, default='all',
        help=
            "Before confound correction, can crop the timeseries within a specific interval.\n"
            "e.g. '0,80' for timepoint 0 to 80.\n"
            "(default: %(default)s)\n"
            "\n"
        )


    ####Analysis
    analysis.add_argument(
        'confound_correction_out', action='store', type=Path,
        help=
            "path to RABIES confound correction output directory.\n"
            "\n"
        )
    analysis.add_argument(
        'output_dir', action='store', type=Path,
        help=
            "path for analysis outputs.\n"
            "\n"
        )
    analysis.add_argument(
        '--scan_list', type=str,
        nargs="*",  # 0 or more values expected => creates a list
        default=['all'],
        help=
            "This option offers to run the analysis on a subset of the scans. The scans are selected by\n"
            "providing the full path to the corresponding EPI file in the input BIDS folder. The list \n"
            "of scan can be specified manually as a list of file name '--scan_list scan1.nii.gz \n"
            "scan2.nii.gz ...' or the files can be imbedded into a .txt file with one filename per row.\n"
            "By default, 'all' will use all the scans previously processed.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--prior_maps', action='store', type=Path,
        default=f"{rabies_path}/melodic_IC.nii.gz",
        help=
            "Provide a 4D nifti image with a series of spatial priors representing common sources of\n"
            "signal (e.g. ICA components from a group-ICA run). This 4D prior map file will be used for \n"
            "Dual regression, Dual ICA and --data_diagnosis. The RABIES default corresponds to a MELODIC \n"
            "run on a combined group of anesthetized-ventilated and awake mice. Confound correction \n"
            "consisted of highpass at 0.01 Hz, FD censoring at 0.03mm, DVARS censoring, and \n"
            "mot_6,WM_signal,CSF_signal as regressors.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--prior_bold_idx', type=int,
        nargs="*",  # 0 or more values expected => creates a list
        default=[5, 12, 19],
        help=
            "Specify the indices for the priors corresponding to BOLD sources from --prior_maps. These will\n"
            "be fitted during Dual ICA and provide the BOLD components during --data_diagnosis.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--prior_confound_idx', type=int,
        nargs="*",  # 0 or more values expected => creates a list
        default=[0, 1, 2, 6, 7, 8, 9, 10, 11,
                    13, 14, 21, 22, 24, 26, 28, 29],
        help=
            "Specify the indices for the confound components from --prior_maps. This is pertinent for the\n" 
            "--data_diagnosis outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        "--data_diagnosis", dest='data_diagnosis', action='store_true',
        help=
            "This option carries out the spatiotemporal diagnosis as described in Desrosiers-Gregoire et al. \n"
            "The diagnosis generates key temporal and spatial features both at the scan level and the group\n"
            "level, allowing the identification of sources of confounds and data quality issues. We recommend \n"
            "using this data diagnosis workflow, more detailed in the publication, to improve the control for \n"
            "data quality issues and prevent the corruptions of analysis outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )

    analysis.add_argument(
        '--seed_list', type=str,
        nargs="*",  # 0 or more values expected => creates a list
        default=[],
        help=
            "Can provide a list of Nifti files providing a mask for an anatomical seed, which will be used\n"
            "to evaluate seed-based connectivity maps using on Pearson's r. Each seed must consist of \n"
            "a binary mask representing the ROI in commonspace.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--seed_prior_list', type=str,
        nargs="*",  # 0 or more values expected => creates a list
        default=[],
        help=
            "For analysis QC of seed-based FC during --data_diagnosis, prior network maps are required for \n"
            "each seed provided in --seed_list. Provide the list of prior files in matching order of the \n"
            "--seed_list arguments to match corresponding seed maps.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument("--FC_matrix", dest='FC_matrix', action='store_true',
        help=
            "Compute whole-brain connectivity matrices using Pearson's r between ROI timeseries.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        "--ROI_type", type=str, default='parcellated',
        choices=['parcellated', 'voxelwise'],
        help=
            "Define ROIs for --FC_matrix between 'parcellated' from the provided atlas during preprocessing,\n"
            "or 'voxelwise' to derive the correlations between every voxel."
            "(default: %(default)s)\n"
            "\n"
        )
    g_group_ICA = analysis.add_argument_group(
        title='Group ICA', 
        description=
            "Options for performing group-ICA using FSL's MELODIC on the whole dataset cleaned timeseries.\n"
            "Note that confound correction must have been conducted on commonspace outputs.\n"
        )
    g_group_ICA.add_argument(
        "--group_ICA", dest='group_ICA', action='store_true',
        help=
            "Perform group-ICA.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_group_ICA.add_argument(
        '--dim', type=int, default=0,
        help=
            "Derive a fixed number of ICA components during group-ICA. The default uses an automatic \n"
            "estimation.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_group_ICA.add_argument(
        '--melodic_random_seed', type=int, default=1,
        help=
            "For reproducibility, can manually set a random seed for MELODIC. \n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        "--DR_ICA", dest='DR_ICA', action='store_true',
        help=
            "Conduct dual regression on each subject timeseries, using the priors from --prior_maps. The\n"
            "linear coefficients from both the first and second regressions will be provided as outputs.\n"
            "Requires that confound correction was conducted on commonspace outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--dual_ICA', type=int, default=0,
        help=
            "Option for performing a Dual ICA. Specify how many subject-specific sources to compute \n"
            "during dual ICA. Dual ICA will provide a fit for each --prior_bold_idx from --prior_maps.\n"
            "(default: %(default)s)\n"
            "\n"
        )

    return parser
