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
            "confound regression and analysis. Outputs from the preprocessing provide the inputs for\n"
            "the subsequent confound regression, and finally analysis.",
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
    g_execution = parser.add_argument_group(title='Execution Options', 
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
    preprocess.add_argument('bids_dir', action='store', type=Path,
                            help="""
                            the root folder of the BIDS-formated input data directory.
                            """)
    preprocess.add_argument('output_dir', action='store', type=Path,
                            help="""
                            the output path to drop outputs from major preprocessing steps.
                            """)
    preprocess.add_argument("--bold_only", dest='bold_only', action='store_true',
                            help="""
                            Apply preprocessing with only EPI scans. Commonspace registration
                            is executed directly using the corrected EPI 3D reference images.
                            The commonspace registration simultaneously applies distortion
                            correction, this option will produce only commonspace outputs.
                            """)
    preprocess.add_argument('--anat_autobox', dest='anat_autobox', action='store_true',
                            help="""
                            Crops out extra space around the brain on the structural image using AFNI's 3dAutobox https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutobox.html.
                            """)
    preprocess.add_argument('--bold_autobox', dest='bold_autobox', action='store_true',
                            help="""
                            Crops out extra space around the brain on the EPI image using AFNI's 3dAutobox https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutobox.html.
                            """)
    preprocess.add_argument('--apply_despiking', dest='apply_despiking', action='store_true',
                            help="""
                            Applies AFNI's 3dDespike https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html.
                            """)
    preprocess.add_argument("--HMC_option", type=str, default='intraSubjectBOLD',
                            choices=['intraSubjectBOLD', '0', '1', '2', '3'],
                            help="""
                            Select an option for head motion realignment among the pre-built options
                            from https://github.com/ANTsX/ANTsR/blob/master/R/ants_motion_estimation.R.
                            """)
    preprocess.add_argument('--apply_slice_mc', dest='apply_slice_mc', action='store_true',
                            help="""
                            Whether to apply a slice-specific motion correction after initial volumetric HMC.
                            This can correct for interslice misalignment resulting from within-TR motion.
                            With this option, motion corrections and the subsequent resampling from registration are applied
                            sequentially, since the 2D slice registrations cannot be concatenate with 3D transforms.
                            """)
    preprocess.add_argument('--detect_dummy', dest='detect_dummy', action='store_true',
                            help="""
                            Detect and remove initial dummy volumes from the EPI, and generate
                            a reference EPI based on these volumes if detected.
                            Dummy volumes will be removed from the output preprocessed EPI.
                            """)
    preprocess.add_argument("--data_type", type=str, default='float32',
                            choices=['int16', 'int32', 'float32', 'float64'],
                            help="Specify data format outputs to control for file size.")

    g_registration = preprocess.add_argument_group(title='Registration Options', description="""
        Customize options for various registration steps. The three in-built registration script options provided
        consist of a 'Rigid', 'Affine' and 'SyN' (non-linear) registration. Other options also allow alternative
        workflows to troubleshoot registration failures.
        """)
    g_registration.add_argument("--bold_inho_cor_method", type=str, default='Rigid',
                            choices=['Rigid', 'Affine',
                                     'SyN', 'no_reg', 'N4_reg', 'disable'],
                            help="""
                            Select a registration type for masking during inhomogeneity correction of the EPI.
                            """)
    g_registration.add_argument("--bold_inho_cor_otsu", type=int, default=2,
                            help="""
                            """)
    g_registration.add_argument("--anat_inho_cor_method", type=str, default='SyN',
                            choices=['Rigid', 'Affine',
                                     'SyN', 'no_reg', 'N4_reg', 'disable'],
                            help="""
                            Select a registration type for masking during inhomogeneity correction of the structural image.
                            """)
    g_registration.add_argument(
        '--atlas_reg_script',
        type=str,
        default='SyN',
        choices=['Rigid', 'Affine', 'SyN', 'NULL'],
        help="""
        Specify a registration script for alignment of the unbiased dataset template to the atlas.
        """)
    g_registration.add_argument("--fast_commonspace", dest='fast_commonspace', action='store_true',
                                help="""
                                This option will skip the generation of a dataset template from https://github.com/CoBrALab/twolevel_ants_dbm.
                                Instead each anatomical scan will be individually registered to the
                                commonspace template using the --atlas_reg_script.
                                Note that this option, although faster, is expected to reduce the
                                quality of commonspace registration.
                                """)
    g_registration.add_argument("--commonspace_masking", dest='commonspace_masking', action='store_true',
                                help="""
                                If true, will use masks originating from the inhomogeneity correction step
                                to orient commonspace alignment.
                                """)
    g_registration.add_argument("--coreg_script", type=str, default='SyN',
                                choices=['Rigid', 'Affine', 'SyN', 'NULL'],
                                help="""
                                Specify EPI to anat coregistration script.
                                """)
    g_registration.add_argument("--coreg_masking", dest='coreg_masking', action='store_true',
                                help="""
                                If true, will use masks originating from the EPI inhomogeneity correction step
                                to orient alignment to the target anatomical image.
                                """)
    g_registration.add_argument("--brain_extraction", dest='brain_extraction', action='store_true',
                                help="""
                                If using masking during registration from --commonspace_masking/--coreg_masking, can
                                specify to extract the brains using the masks prior to registration, which will 
                                enhance brain edge-matching, but requires good quality masks.
                                """)

    g_resampling = preprocess.add_argument_group(title='Resampling Options', description="""
        The following options allow to customize the voxel dimensions for the preprocessed EPIs or for
        the anatomical images during registration.
        Axis resampling specifications must follow the format 'dim1xdim2xdim3' (in mm) with the RAS axis
        convention (dim1=Right-Left, dim2=Anterior-Posterior, dim3=Superior-Inferior).
        The original dimensions are conserved if 'inputs_defined' is specified.
        """)
    g_resampling.add_argument('--nativespace_resampling', type=str, default='inputs_defined',
                              help="""
                              Can specify a resampling dimension for the nativespace outputs.
                              """)
    g_resampling.add_argument('--commonspace_resampling', type=str, default='inputs_defined',
                              help="""
                              Can specify a resampling dimension for the commonspace outputs.
                              """)
    g_resampling.add_argument('--anatomical_resampling', type=str, default='inputs_defined',
                              help="""
                              Can specify resampling dimensions for the template files to optimize
                              registration efficiency. By defaults ('inputs_defined'), the resampling
                              dimension is estimated from the input images. The smallest dimension among
                              the anatomical images (EPI images instead if --bold_only is True) defines
                              the isotropic resolution for resampling.
                              """)

    g_stc = preprocess.add_argument_group(title='STC Options', description="""
        Specify Slice Timing Correction (STC) info that is fed to AFNI 3dTshift
        (https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html). The STC is applied in the
        anterior-posterior orientation, assuming slices were acquired in this direction.
        """)
    g_stc.add_argument('--TR', type=str, default='auto',
                       help="""
                       Specify repetition time (TR) in seconds. (e.g. --TR 1.2)
                       """)
    g_stc.add_argument('--apply_STC', dest='apply_STC', action='store_true',
                       help="""
                       Select this option to apply the STC step.
                       """)
    g_stc.add_argument('--tpattern', type=str, default='alt',
                       choices=['alt', 'seq'],
                       help="""
                       Specify if interleaved or sequential acquisition. 'alt' for interleaved, 'seq' for sequential.
                       """)

    g_atlas = preprocess.add_argument_group(title='Template Files', description="""
        Specify commonspace template and associated mask/label files.
        A mouse atlas is provided as default https://wiki.mouseimaging.ca/display/MICePub/Mouse+Brain+Atlases.
        """)
    g_atlas.add_argument('--anat_template', action='store', type=Path,
                         default=f"{rabies_path}/DSURQE_40micron_average.nii.gz",
                         help="""
                         Anatomical file for the commonspace template.
                         """)
    g_atlas.add_argument('--brain_mask', action='store', type=Path,
                         default=f"{rabies_path}/DSURQE_40micron_mask.nii.gz",
                         help="""
                         Brain mask for the template.
                         """)
    g_atlas.add_argument('--WM_mask', action='store', type=Path,
                         default=f"{rabies_path}/DSURQE_40micron_eroded_WM_mask.nii.gz",
                         help="""
                         White matter mask for the template.
                         """)
    g_atlas.add_argument('--CSF_mask', action='store', type=Path,
                         default=f"{rabies_path}/DSURQE_40micron_eroded_CSF_mask.nii.gz",
                         help="""
                         CSF mask for the template.
                         """)
    g_atlas.add_argument('--vascular_mask', action='store', type=Path,
                         default=f"{rabies_path}/vascular_mask.nii.gz",
                         help="""
                         Can provide a mask of major blood vessels for computing confound timeseries.
                         The default mask was generated by applying MELODIC ICA and selecting the resulting component
                         mapping onto major veins.
                         """)
    g_atlas.add_argument('--labels', action='store', type=Path,
                         default=f"{rabies_path}/DSURQE_40micron_labels.nii.gz",
                         help="""
                         Atlas file with anatomical labels.
                         """)



    ####Confound correction
    confound_correction.add_argument('preprocess_out', action='store', type=Path,
                                     help="""
                                     path to RABIES preprocessing output directory with the datasinks.
                                     """)
    confound_correction.add_argument('output_dir', action='store', type=Path,
                                     help="""
                                     path to drop confound regression output datasink.
                                     """)
    confound_correction.add_argument('--read_datasink', dest='read_datasink', action='store_true', default=False,
                                     help="""
                                     Choose this option to read preprocessing outputs from datasinks instead
                                     of the saved preprocessing workflow graph. This allows to run confound correction
                                     without having available RABIES preprocessing folders.
                                     Using this option, it is assumed that outputs in the datasink folders are in line
                                     with the expected outputs from preprocessing with RABIES.
                                     """)
    confound_correction.add_argument('--output_name', type=str, default='confound_correction_wf',
                                     help="""
                                     Creates a new output folder to store the workflow of this CR run, to avoid potential
                                     overlaps with previous runs (can be useful if investigating multiple strategies).
                                     """)
    confound_correction.add_argument('--nativespace_analysis', dest='nativespace_analysis', action='store_true',
                                     help="""
                                     Use to specify confound correction and analysis on native space outputs.
                                     """)
    confound_correction.add_argument('--TR', type=str, default='auto',
                                     help="""
                                     Specify repetition time (TR) in seconds. (e.g. --TR 1.2)
                                     """)
    confound_correction.add_argument('--highpass', type=float, default=None,
                                     help="""
                                     Specify highpass filter frequency.
                                     """)
    confound_correction.add_argument('--lowpass', type=float, default=None,
                                     help="""
                                     Specify lowpass filter frequency.
                                     """)
    confound_correction.add_argument('--edge_cutoff', type=float, default=0,
                                     help="""
                                     Specify number of seconds to cut at beginning and end of acquisition if applying 
                                     a frequency filter.
                                     Applying frequency filters generate edge effects at begining and end of the sequence. 
                                     We recommend to cut those timepoints (around 30sec at both end for 0.01Hz highpass.).
                                     """)
    confound_correction.add_argument('--smoothing_filter', type=float, default=None,
                                     help="""
                                     Specify spatial smoothing filter size in mm.
                                     Uses nilearn's function https://nilearn.github.io/modules/generated/nilearn.image.smooth_img.html
                                     """)
    confound_correction.add_argument('--run_aroma', dest='run_aroma', action='store_true', default=False,
                                     help="""
                                     Whether to run ICA-AROMA or not. The classifier implemented within RABIES
                                     is a slightly modified version from the original (Pruim et al. 2015),
                                     with parameters and masks adapted for rodent images.
                                     """)
    confound_correction.add_argument('--aroma_dim', type=int, default=0,
                                     help="""
                                     Can specify a fixed number of dimension for the MELODIC run before ICA-AROMA.
                                     """)
    confound_correction.add_argument('--aroma_random_seed', type=int, default=1,
                             help="""
                             For reproducibility, can manually set a random seed for MELODIC.
                             """)
    confound_correction.add_argument('--conf_list', type=str,
                                     nargs="*",  # 0 or more values expected => creates a list
                                     default=[],
                                     choices=["WM_signal", "CSF_signal", "vascular_signal",
                                              "global_signal", "aCompCor", "mot_6", "mot_24", "mean_FD"],
                                     help="""
                                     List of nuisance regressors that will be applied on voxel timeseries, i.e., confound regression.
                                     WM/CSF/vascular/global_signal: correspond to mean signal from WM/CSF/vascular/brain masks.
                                     mot_6: 6 rigid HMC parameters.
                                     mot_24: mot_6 + their temporal derivative, and all 12 parameters squared (Friston et al. 1996).
                                     aCompCor: corresponds to the timeseries of components from a PCA conducted on the combined
                                     WM and CSF masks voxel timeseries, including all components that together explain 50/100 of
                                     the variance, as in Muschelli et al. 2014.
                                     mean_FD: the mean framewise displacement timecourse
                                     """)
    confound_correction.add_argument('--FD_censoring', dest='FD_censoring', action='store_true', default=False,
                                     help="""
                                     Whether to remove timepoints that exceed a framewise displacement threshold.
                                     The frames that exceed the given threshold together with 1 back
                                     and 4 forward frames will be masked out (based on Power et al. 2012).
                                     """)
    confound_correction.add_argument('--FD_threshold', type=float, default=0.05,
                                     help="""
                                     Scrubbing threshold for the mean framewise displacement in mm (averaged across the brain
                                     mask) to select corrupted volumes.
                                     """)
    confound_correction.add_argument('--DVARS_censoring', dest='DVARS_censoring', action='store_true',default=False,
                                     help="""
                                     Whether to remove timepoints that present outlier values on the DVARS metric (temporal derivative
                                     of global signal).
                                     This will censor out timepoints until a distribution of DVARS values is obtained without outliers
                                     values above or below 2.5 standard deviations.
                                     """)
    confound_correction.add_argument('--minimum_timepoint', type=int,default=3,
                                     help="""
                                     Can select a threshold number of timepoint to remain after censoring, and return empty files for
                                     scans that don't pass threshold.
                                     """)
    confound_correction.add_argument('--standardize', dest='standardize', action='store_true',default=False,
                                     help="""
                                     Whether to standardize timeseries to unit variance.
                                     """)
    confound_correction.add_argument('--timeseries_interval', type=str, default='all',
                                     help="""
                                     Specify a time interval in the timeseries to keep. e.g. "0,80". By default all timeseries are kept.
                                     """)


    ####Analysis
    analysis.add_argument('confound_correction_out', action='store', type=Path,
                          help="""
                          path to RABIES confound regression output directory with the datasink.
                          """)
    analysis.add_argument('output_dir', action='store', type=Path,
                          help='path to drop analysis outputs.')
    analysis.add_argument('--output_name', type=str, default='analysis_wf',
                          help="""
                          Creates a new output folder to store the workflow of this analysis run, to avoid potential
                          overlaps with previous runs.
                          """)
    analysis.add_argument('--scan_list', type=str,
                          nargs="*",  # 0 or more values expected => creates a list
                          default=['all'],
                          help="""
                          This option offers to run the analysis on a subset of the scans.
                          The scans selected are specified by providing the full path to each EPI file from the input BIDS folder.
                          The list of scan can be specified manually as a list of file name '--scan_list scan1.nii.gz scan2.nii.gz ...'
                          or the files can be imbedded into a .txt file with one filename per row.
                          By default, 'all' will use all the scans previously processed.
                          """)
    analysis.add_argument('--seed_list', type=str,
                          nargs="*",  # 0 or more values expected => creates a list
                          default=[],
                          help="""
                          Can provide a list of seed .nii images that will be used to evaluate seed-based correlation maps
                          based on Pearson's r. Each seed must consist of a binary mask representing the ROI in commonspace.
                          """)
    analysis.add_argument('--seed_prior_list', type=str,
                          nargs="*",  # 0 or more values expected => creates a list
                          default=[],
                          help="""
                          For analysis QC of seed-based FC, prior network maps are required for each seed-FC provided in --seed_list.
                          Provide the list of prior files in matching order of the --seed_list inputs to match corresponding seed maps.
                          """)
    analysis.add_argument('--prior_maps', action='store', type=Path,
                            default=f"{rabies_path}/melodic_IC.nii.gz",
                            help="""
                            Provide a 4D nifti image with a series of spatial priors representing common sources of
                            signal (e.g. ICA components from a group-ICA run). This 4D prior map file will be used for 
                            Dual regression, Dual ICA and --data_diagnosis.
                            Default: Corresponds to a MELODIC run on a combined group of anesthetized-
                            ventilated and awake mice. Confound regression consisted of highpass at 0.01 Hz,
                            FD censoring at 0.03mm, DVARS censoring, and mot_6,WM_signal,CSF_signal as regressors.
                            """)
    analysis.add_argument('--prior_bold_idx', type=int,
                            nargs="*",  # 0 or more values expected => creates a list
                            default=[5, 12, 19],
                            help="""
                            Specify the indices for the priors to fit from --prior_maps. Only selected priors will be fitted 
                            for Dual ICA, and these priors will correspond to the BOLD components during --data_diagnosis.
                            """)
    analysis.add_argument('--prior_confound_idx', type=int,
                                nargs="*",  # 0 or more values expected => creates a list
                                default=[0, 1, 2, 6, 7, 8, 9, 10, 11,
                                         13, 14, 21, 22, 24, 26, 28, 29],
                                help="""
                                Specify the indices for the confound components from --prior_maps.
                                This is pertinent for the --data_diagnosis outputs.
                                """)
    analysis.add_argument("--data_diagnosis", dest='data_diagnosis', action='store_true',
                             help="""
                             This option carries out the spatiotemporal diagnosis as described in Desrosiers-Gregoire et al.
                             The diagnosis outputs key temporal and spatial features at the scan level allowing the
                             identification of sources of confounds in individual scans. A follow-up group-level correlation
                             between spatial features is also conducted to evaluate corruption of group-level outputs.
                             We recommend to conduct a data diagnosis from this workflow to complement FC analysis by
                             evaluating the intrinsic corruption of the dataset as well as the effectiveness of the confound
                             correction strategies.""")

    g_fc_matrix = analysis.add_argument_group(title='FC Matrix', description="""
        Options for performing a whole-brain timeseries correlation matrix analysis.
        """)
    g_fc_matrix.add_argument("--FC_matrix", dest='FC_matrix', action='store_true',
                             help="""
                             Choose this option to derive a whole-brain functional connectivity matrix, based on the
                             Pearson's r correlation of regional timeseries for each subject cleaned timeseries.
                             """)
    g_fc_matrix.add_argument("--ROI_type", type=str, default='parcellated',
                             choices=['parcellated', 'voxelwise'],
                             help="""
                             Define the types of ROI to extract regional timeseries for correlation matrix analysis.
                             Options are 'parcellated', in which case the atlas labels provided for preprocessing are used
                             as ROIs, or 'voxelwise', in which case all voxel timeseries are cross-correlated.
                             """)
    g_group_ICA = analysis.add_argument_group(title='Group ICA', description="""
        Options for performing group-ICA using FSL's MELODIC on the whole dataset cleaned timeseries.
        Note that confound regression must have been conducted on commonspace outputs.
        """)
    g_group_ICA.add_argument("--group_ICA", dest='group_ICA', action='store_true',
                             help="""
                             Choose this option to conduct group-ICA.
                             """)
    g_group_ICA.add_argument('--dim', type=int, default=0,
                             help="""
                             You can specify the number of ICA components to be derived. The default uses an automatic estimation.
                             """)
    g_group_ICA.add_argument('--melodic_random_seed', type=int, default=1,
                             help="""
                             For reproducibility, can manually set a random seed for MELODIC.
                             """)
    g_DR_ICA = analysis.add_argument_group(title='DR ICA', description="""
        Options for performing a dual regression analysis based on a previous group-ICA run from FSL's MELODIC.
        Note that confound regression must have been conducted on commonspace outputs.
        """)
    g_DR_ICA.add_argument("--DR_ICA", dest='DR_ICA', action='store_true',
                          help="""
                          Choose this option to conduct dual regression on each subject timeseries. This analysis will
                          output the spatial maps corresponding to the linear coefficients from the second linear
                          regression. See rabies.analysis_pkg.analysis_functions.dual_regression for the specific code.
                          """)
    g_dual_ICA = analysis.add_argument_group(title='Dual ICA', description="""
        Options for performing a Dual ICA.
        Need to provide the prior maps to fit --prior_maps, and the associated indices for the target components
        in --prior_bold_idx
        """)
    g_dual_ICA.add_argument('--dual_ICA', type=int, default=0,
                            help="""
                            Specify how many subject-specific sources to compute using dual ICA.
                            """)

    return parser
