import os
import sys
import pickle
import logging
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
        description="""RABIES performs processing of rodent fMRI images. Can either run
        on datasets that only contain EPI images, or both structural and EPI images.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(title='Commands',
                                       description='The RABIES workflow is seperated into three different processing steps: preprocessing, '
                                                   'confound regression and analysis. Outputs from the preprocessing provides the inputs for '
                                                   'the subsequent confound regression, and finally analysis.',
                                       help='Description',
                                       dest='rabies_step',
                                       metavar='Processing step')

    preprocess = subparsers.add_parser("preprocess",
                                       help="""Conducts preprocessing on an input dataset in BIDS format.
        Preprocessing includes realignment for motion, correction for susceptibility distortions through non-linear registration,
        registration to a commonspace atlas and associated masks, evaluation of confounding timecourses, and includes various
        execution options (see --help).
        """, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    confound_regression = subparsers.add_parser("confound_regression",
                                                help="""Different options for confound regression are available to apply directly
        on preprocessing outputs from RABIES. Detrending of timeseries is always applied. Otherwise only
        selected confound regression and denoising strategies are applied.
        The denoising steps are applied in the following order: detrending first then ICA-AROMA, then highpass/lowpass filtering followed by temporal censoring,
        then regression of confound timeseries, standardization of timeseries, and finally smoothing.
        """, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    analysis = subparsers.add_parser("analysis",
                                     help="""
        A few built-in resting-state functional connectivity (FC) analysis options are provided to conduct rapid analysis on the cleaned timeseries.
        The options include seed-based FC, voxelwise or parcellated whole-brain FC, group-ICA and dual regression.
        """, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    data_diagnosis = subparsers.add_parser("data_diagnosis",
                                     help="""
        This workflow regroups a set of tools which allow to establish the influence of confounding sources on FC measures. We recommend to
        conduct a data diagnosis from this workflow to complement FC analysis by evaluating the effectiveness of the confound correction strategies
        and ensure reliable subsequent interpretations.
        """, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    g_execution = parser.add_argument_group(
        "Options for managing the execution of the workflow.")
    g_execution.add_argument("-p", "--plugin", default='Linear',
                             choices=['Linear', 'MultiProc', 'SGE', 'SGEGraph', 'PBS', 'LSF', 'SLURM', 'SLURMGraph'],
                             help="Specify the nipype plugin for workflow execution. Consult nipype plugin documentation for detailed options."
                             " Linear, MultiProc, SGE and SGEGraph have been tested.")
    g_execution.add_argument('--local_threads', type=int, default=multiprocessing.cpu_count(),
                             help="""For local MultiProc execution, set the maximum number of processors run in parallel,
        defaults to number of CPUs. This option only applies to the MultiProc execution plugin, otherwise
        it is set to 1.""")
    g_execution.add_argument("--scale_min_memory", type=float, default=1.0,
                             help="For a parallel execution with MultiProc, the minimal memory attributed to nodes can be scaled with this multiplier to avoid memory crashes.")
    g_execution.add_argument("--min_proc", type=int, default=1,
                             help="For SGE parallel processing, specify the minimal number of nodes to be assigned to avoid memory crashes.")

    preprocess.add_argument('bids_dir', action='store', type=Path,
                            help='the root folder of the BIDS-formated input data directory.')
    preprocess.add_argument('output_dir', action='store', type=Path,
                            help='the output path to drop outputs from major preprocessing steps.')
    preprocess.add_argument("-e", "--bold_only", dest='bold_only', action='store_true',
                            help="Apply preprocessing with only EPI scans. commonspace registration"
                            " is executed through registration of the EPI-generated template from ants_dbm"
                            " to the anatomical template.")
    preprocess.add_argument("--bold_bias_cor_method", type=str, default='mouse-preprocessing-v5.sh',
                            choices=['otsu_reg', 'thresh_reg', 'mouse-preprocessing-v5.sh'],
                            help="Choose the algorithm for bias field correction of the EPI before registration."
                            "otsu_reg will conduct an initial serie of N4BiasFieldCorrection oriented by Otsu masking method, "
                            "followed by a rigid registration to provide a brain mask orienting the final correction."
                            "thresh_reg will instead use an initial voxel intensity thresholding method for masking, and will "
                            "conduct a subsequent rigid registration to provide a brain mask orienting the final correction.")
    preprocess.add_argument("--anat_bias_cor_method", type=str, default='mouse-preprocessing-v5.sh',
                            choices=['otsu_reg', 'thresh_reg', 'mouse-preprocessing-v5.sh'],
                            help="Same as --bold_bias_cor_method but for the anatomical image.")
    preprocess.add_argument("--disable_anat_preproc", dest='disable_anat_preproc', action='store_true',
                            help="This option disables the preprocessing of anatomical images before commonspace template generation.")
    preprocess.add_argument('--apply_despiking', dest='apply_despiking', action='store_true',
                            help="Whether to apply despiking of the EPI timeseries based on AFNI's "
                            "3dDespike https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html.")
    preprocess.add_argument("--HMC_option", type=str, default='intraSubjectBOLD',
                            choices=['intraSubjectBOLD', '0', '1', '2', '3'],
                            help="Select an option for head motion realignment among the pre-built options "
                            "from https://github.com/ANTsX/ANTsR/blob/60eefd96fedd16bceb4703ecd2cd5730e6843807/R/ants_motion_estimation.R.")
    preprocess.add_argument("--HMC_transform", type=str, default='Rigid',
                            choices=['Rigid', 'Affine'],
                            help="Select between Rigid and Affine registration for head motion realignment.")
    preprocess.add_argument('--apply_slice_mc', dest='apply_slice_mc', action='store_true',
                            help="Whether to apply a slice-specific motion correction after initial volumetric rigid correction. "
                            "This second motion correction can correct for interslice misalignment resulting from within-TR motion."
                            "With this option, motion corrections and the subsequent resampling from registration are applied sequentially,"
                            "since the 2D slice registrations cannot be concatenate with 3D transforms.")
    preprocess.add_argument('--detect_dummy', dest='detect_dummy', action='store_true',
                            help="Detect and remove initial dummy volumes from the EPI, and generate "
                            "a reference EPI based on these volumes if detected."
                            "Dummy volumes will be removed from the output preprocessed EPI.")
    preprocess.add_argument("--data_type", type=str, default='float32',
                            choices=['int16', 'int32', 'float32', 'float64'],
                            help="Specify data format outputs to control for file size.")
    preprocess.add_argument("--debug", dest='debug', action='store_true',
                            help="Run in debug mode.")

    g_registration = preprocess.add_argument_group(
        "Options for the registration steps. Built-in options for selecting registration scripts include 'Rigid', 'Affine', 'SyN' (non-linear), 'light_SyN', 'heavy_SyN', 'multiRAT', but"
        " can specify a custom registration script following the template script structure (see RABIES/rabies/shell_scripts/ for template)."
        "'Rigid', 'Affine' and 'SyN' options rely on an adaptive registration framework which adapts parameters to the images dimensions")
    g_registration.add_argument("--coreg_script", type=str, default='SyN',
                                help="Specify EPI to anat coregistration script.")
    g_registration.add_argument(
        '--template_reg_script',
        type=str,
        default='SyN',
        help="""Registration script that will be used for registration of the generated dataset
        template to the provided commonspace atlas for masking and labeling.""")
    g_registration.add_argument("--fast_commonspace", dest='fast_commonspace', action='store_true',
                                help="Choosing this option will skip the generation of a dataset template, and instead, "
                                "each anatomical scan will be individually registered to the commonspace template using the --template_reg_script."
                                "Note that this option, although faster, is expected to reduce the quality of commonspace registration.")

    g_resampling = preprocess.add_argument_group("Options for the resampling of the EPI. "
                                                 "Axis resampling specifications must follow the format 'dim1xdim2xdim3' (in mm) with the RAS axis convention (dim1=Right-Left, dim2=Anterior-Posterior, dim3=Superior-Inferior).")
    g_resampling.add_argument('--nativespace_resampling', type=str, default='origin',
                              help="Can specify a resampling dimension for the nativespace outputs. Must be of the form dim1xdim2xdim3 (in mm). The original dimensions are conserved "
                              "if 'origin' is specified.")
    g_resampling.add_argument('--commonspace_resampling', type=str, default='origin',
                              help="Can specify a resampling dimension for the commonspace outputs. Must be of the form dim1xdim2xdim3 (in mm). The original dimensions are conserved "
                              "if 'origin' is specified."
                              "***this option specifies the resampling for the --bold_only workflow")
    g_resampling.add_argument(
        '--anatomical_resampling', type=str, default='inputs_defined',
        help="""To optimize the efficiency of registration, the provided anatomical template is resampled based on the provided
        input images. The dimension with the lowest resolution among the provided anatomical images (EPI images instead if --bold_only is True)
        is selected as a basis for resampling the template to isotropic resolution, if the provided resolution is lower than the original
        resolution of the template. Alternatively, the user can provide a custom resampling dimension. This allows to accelerate
        registration steps with minimal sampling dimensions.""")

    g_ants_dbm = preprocess.add_argument_group(
        'cluster options for running ants_dbm (options copied from twolevel_dbm.py):')
    g_ants_dbm.add_argument(
        '--cluster_type',
        default="local",
        choices=["local", "sge", "pbs", "slurm"],
        help="Choose the type of cluster system to submit jobs to")
    g_ants_dbm.add_argument(
        '--walltime',
        default="20:00:00",
        help="""Option for job submission
        specifying requested time per pairwise registration.""")

    g_stc = preprocess.add_argument_group("""Specify Slice Timing Correction info that is fed to AFNI 3dTshift
    (https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html). The STC is applied in the
    anterior-posterior orientation, assuming slices were acquired in this direction.""")
    g_stc.add_argument('--TR', type=str, default='1.0s',
                       help="Specify repetition time (TR) in seconds.")
    g_stc.add_argument('--no_STC', dest='no_STC', action='store_true',
                       help="Select this option to ignore the STC step.")
    g_stc.add_argument('--tpattern', type=str, default='alt',
                       choices=['alt', 'seq'],
                       help="Specify if interleaved or sequential acquisition. 'alt' for interleaved, 'seq' for sequential.")

    g_atlas = preprocess.add_argument_group(
        'Provided commonspace atlas files.')
    g_atlas.add_argument('--anat_template', action='store', type=Path,
                         default="%s/DSURQE_40micron_average.nii.gz" % (
                             rabies_path),
                         help='Anatomical file for the commonspace template.')
    g_atlas.add_argument('--brain_mask', action='store', type=Path,
                         default="%s/DSURQE_40micron_mask.nii.gz" % (
                             rabies_path),
                         help='Brain mask for the template.')
    g_atlas.add_argument('--WM_mask', action='store', type=Path,
                         default="%s/DSURQE_40micron_eroded_WM_mask.nii.gz" % (
                             rabies_path),
                         help='White matter mask for the template.')
    g_atlas.add_argument('--CSF_mask', action='store', type=Path,
                         default="%s/DSURQE_40micron_eroded_CSF_mask.nii.gz" % (
                             rabies_path),
                         help='CSF mask for the template.')
    g_atlas.add_argument('--vascular_mask', action='store', type=Path,
                         default="%s/vascular_mask.nii.gz" % (
                             rabies_path),
                         help='Can provide a mask of major blood vessels for computing confound timeseries. The default mask was generated by applying MELODIC ICA and selecting the resulting component mapping onto major veins. (Grandjean et al. 2020, NeuroImage; Beckmann et al. 2005)')
    g_atlas.add_argument('--labels', action='store', type=Path,
                         default="%s/DSURQE_40micron_labels.nii.gz" % (
                             rabies_path),
                         help='Atlas file with anatomical labels.')

    confound_regression.add_argument('preprocess_out', action='store', type=Path,
                                     help='path to RABIES preprocessing output directory with the datasinks.')
    confound_regression.add_argument('output_dir', action='store', type=Path,
                                     help='path to drop confound regression output datasink.')
    confound_regression.add_argument('--wf_name', type=str, default='confound_regression_wf',
                                     help='Can specify a name for the workflow of this confound regression run, to avoid potential '
                                     'overlaps with previous runs (can be useful if investigating multiple strategies).')
    confound_regression.add_argument('--commonspace_bold', dest='commonspace_bold', action='store_true',
                                     help='If should run confound regression on the commonspace bold output.')
    confound_regression.add_argument('--TR', type=str, default='1.0s',
                                     help="Specify repetition time (TR) in seconds.")
    confound_regression.add_argument('--highpass', type=float, default=None,
                                     help='Specify highpass filter frequency.')
    confound_regression.add_argument('--lowpass', type=float, default=None,
                                     help='Specify lowpass filter frequency.')
    confound_regression.add_argument('--smoothing_filter', type=float, default=None,
                                     help='Specify smoothing filter size in mm.')
    confound_regression.add_argument('--run_aroma', dest='run_aroma', action='store_true',
                                     default=False,
                                     help="Whether to run ICA-AROMA or not. The classifier implemented within RABIES "
                                     "is a slightly modified version from the original (Pruim et al. 2015), with parameters and masks adapted for rodent images.")
    confound_regression.add_argument('--aroma_dim', type=int,
                                     default=0,
                                     help='Can specify a number of dimension for the MELODIC run before ICA-AROMA.')
    confound_regression.add_argument('--conf_list', type=str,
                                     nargs="*",  # 0 or more values expected => creates a list
                                     default=[],
                                     choices=["WM_signal", "CSF_signal", "vascular_signal",
                                              "global_signal", "aCompCor", "mot_6", "mot_24", "mean_FD"],
                                     help="list of nuisance regressors that will be applied on voxel timeseries. mot_6 corresponds to the 6 rigid body parameters, "
                                     "and mot_24 corresponds to the 6 rigid parameters, their temporal derivative, and all 12 parameters squared (Friston et al. 1996). "
                                     "aCompCor corresponds the timeseries of components from a PCA conducted on the combined WM and CSF masks voxel timeseries, including "
                                     "all components that together explain 50 percent. of the variance, as in Muschelli et al. 2014.")
    confound_regression.add_argument('--FD_censoring', dest='FD_censoring', action='store_true',
                                     default=False,
                                     help="""Whether to remove timepoints that exceed a framewise displacement threshold, before other confound regression steps.
                                     The frames that exceed the given threshold together with 1 back and 4 forward frames will be masked out (based on Power et al. 2012).""")
    confound_regression.add_argument('--FD_threshold', type=float,
                                     default=0.05,
                                     help='Scrubbing threshold for the mean framewise displacement in mm (averaged across the brain mask) to select corrupted volumes.')
    confound_regression.add_argument('--DVARS_censoring', dest='DVARS_censoring', action='store_true',
                                     default=False,
                                     help="""Whether to remove timepoints that present outlier values on the DVARS metric (temporal derivative of global signal).
                                     This option will create a distribution of DVARS values which has no outlier values above or below 2.5 standard deviations.""")
    confound_regression.add_argument('--minimum_timepoint', type=int,
                                     default=3,
                                     help='Can specify a minimal number of timepoint to remain in the timeseries after censoring, otherwise it returns empty files.')
    confound_regression.add_argument('--standardize', dest='standardize', action='store_true',
                                     default=False,
                                     help="""Whether to standardize timeseries to unit variance.""")
    confound_regression.add_argument('--timeseries_interval', type=str, default='all',
                                     help='Specify a time interval in the timeseries to keep. e.g. "0,80". By default all timeseries are kept.')

    analysis.add_argument('confound_regression_out', action='store', type=Path,
                          help='path to RABIES confound regression output directory with the datasink.')
    analysis.add_argument('output_dir', action='store', type=Path,
                          help='the output path to drop analysis outputs.')
    analysis.add_argument('--wf_name', type=str, default='analysis_wf',
                                     help='Can specify a name for the workflow of this analysis run, to avoid potential '
                                     'overlaps with previous runs.')
    analysis.add_argument('--scan_list', type=str,
                                     nargs="*",  # 0 or more values expected => creates a list
                                     default=['all'],
                                     help="This option offers to run the analysis on a subset of the scans."
                                     "The scans selected are specified by providing the full path to each EPI file from the input BIDS folder."
                                     "The list of scan can be specified manually as a list of file name '--scan_list scan1.nii.gz scan2.nii.gz ...' "
                                     "or the files can be imbedded into a .txt file with one filename per row."
                                     "By default, 'all' will use all the scans previously processed."
                                     )
    analysis.add_argument('--seed_list', type=str,
                                     nargs="*",  # 0 or more values expected => creates a list
                                     default=[],
                                     help="Can provide a list of seed .nii images that will be used to evaluate seed-based correlation maps based on Pearson's r."
                                     "Each seed must consist of a binary mask representing the ROI in commonspace.")
    g_fc_matrix = analysis.add_argument_group(
        'Options for performing a whole-brain timeseries correlation matrix analysis.')
    g_fc_matrix.add_argument("--FC_matrix", dest='FC_matrix', action='store_true',
                             help="Choose this option to derive a whole-brain functional connectivity matrix, based on the Pearson's r correlation of regional timeseries "
                             "for each subject cleaned timeseries.")
    g_fc_matrix.add_argument("--ROI_type", type=str, default='parcellated',
                             choices=['parcellated', 'voxelwise'],
                             help="Define the types of ROI to extract regional timeseries for correlation matrix analysis. "
                             "Options are 'parcellated', in which case the atlas labels provided for preprocessing are used as ROIs, or "
                             "'voxelwise', in which case all voxel timeseries are cross-correlated.")
    g_group_ICA = analysis.add_argument_group("Options for performing group-ICA using FSL's MELODIC on the whole dataset cleaned timeseries."
                                              "Note that confound regression must have been conducted on commonspace outputs.")
    g_group_ICA.add_argument("--group_ICA", dest='group_ICA', action='store_true',
                             help="Choose this option to conduct group-ICA.")
    g_group_ICA.add_argument('--TR', type=str, default='1.0s',
                             help="Specify repetition time (TR) in seconds.")
    g_group_ICA.add_argument('--dim', type=int, default=0,
                             help="You can specify the number of ICA components to be derived. The default uses an automatic estimation.")
    g_DR_ICA = analysis.add_argument_group("Options for performing a dual regression analysis based on a previous group-ICA run from FSL's MELODIC. "
                                           "Note that confound regression must have been conducted on commonspace outputs.")
    g_DR_ICA.add_argument("--DR_ICA", dest='DR_ICA', action='store_true',
                          help="Choose this option to conduct dual regression on each subject timeseries. This analysis will output the spatial "
                          "maps corresponding to the linear coefficients from the second linear regression. "
                          "See rabies.analysis_pkg.analysis_functions.dual_regression for the specific code.")
    g_DR_ICA.add_argument('--IC_file', action='store', type=Path,
                          default=None,
                          help="Option to provide a melodic_IC.nii.gz file with the ICA components from a previous group-ICA run. "
                          "If none is provided, a group-ICA will be run with the dataset cleaned timeseries.")

    data_diagnosis.add_argument('confound_regression_out', action='store', type=Path,
                          help='path to RABIES confound regression output directory with the datasink.')
    data_diagnosis.add_argument('output_dir', action='store', type=Path,
                          help='the output path to drop data_diagnosis outputs.')
    data_diagnosis.add_argument('--wf_name', type=str, default='data_diagnosis_wf',
                                     help='Can specify a name for the workflow of this data diagnosis run, to avoid potential '
                                     'overlaps with previous runs.')
    data_diagnosis.add_argument('--scan_list', type=str,
                                     nargs="*",  # 0 or more values expected => creates a list
                                     default=['all'],
                                     help="This option offers to run the analysis on a subset of the scans."
                                     "The scans selected are specified by providing the full path to each EPI file from the input BIDS folder."
                                     "The list of scan can be specified manually as a list of file name '--scan_list scan1.nii.gz scan2.nii.gz ...' "
                                     "or the files can be imbedded into a .txt file with one filename per row."
                                     "By default, 'all' will use all the scans previously processed."
                                     )
    data_diagnosis.add_argument('--dual_convergence', type=int, default=0,
                             help="Can specify a number of components to compute using a dual convergence framework.")
    data_diagnosis.add_argument('--prior_maps', action='store', type=Path,
                          default="%s/melodic_IC.nii.gz" % (
                              rabies_path),
                          help="Provide a 4D nifti image with a series of spatial priors representing common sources of signal (e.g. ICA components from a group-ICA run)."
                          "The default file corresponds to a MELODIC run on a combined group of anesthetized-ventilated mice with MEDISO and awake mice."
                          "Confound regression consisted of highpass at 0.01 Hz, FD censoring at 0.03mm, DVARS censoring, and mot_6,WM_signal,CSF_signal as regressors."
                          )
    data_diagnosis.add_argument('--prior_bold_idx', type=str,
                                     nargs="*",  # 0 or more values expected => creates a list
                                     default=[5,12,19],
                                     help="")
    data_diagnosis.add_argument('--prior_confound_idx', type=str,
                                     nargs="*",  # 0 or more values expected => creates a list
                                     default=[0,1,2,6,7,8,9,10,11,13,14,21,22,24,26,28,29],
                                     help="")
    data_diagnosis.add_argument('--DSURQE_regions', dest='DSURQE_regions', action='store_true',
                                     default=False,
                                     help=""" """)

    return parser


def execute_workflow():
    # generates the parser CLI and execute the workflow based on specified parameters.
    parser = get_parser()
    opts = parser.parse_args()

    output_folder = os.path.abspath(str(opts.output_dir))
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    # Shouldn't be using more than one thread if not MultiProc
    if not opts.plugin == 'MultiProc':
        opts.local_threads = 1

    # managing log info
    cli_file = '%s/rabies_%s.pkl' % (output_folder, opts.rabies_step, )
    with open(cli_file, 'wb') as handle:
        pickle.dump(opts, handle, protocol=pickle.HIGHEST_PROTOCOL)

    logging.basicConfig(filename='%s/rabies_%s.log' % (output_folder, opts.rabies_step, ), filemode='w',
                        format='%(asctime)s - %(levelname)s - %(message)s', level=os.environ.get("LOGLEVEL", "INFO"))
    log = logging.getLogger('root')

    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)

    # verify default template installation
    install_DSURQE(log)

    from .__version__ import __version__
    log.info('Running RABIES - version: '+__version__)

    # print complete CLI command
    args = 'CLI INPUTS: \n'
    for arg in vars(opts):
        input = '-> {arg} = {value} \n'.format(
            arg=arg, value=getattr(opts, arg))
        args += input
    log.info(args)

    if opts.rabies_step == 'preprocess':
        workflow = preprocess(opts, None, None, None, log)
    elif opts.rabies_step == 'confound_regression':
        workflow = confound_regression(opts, None, None, log)
    elif opts.rabies_step == 'analysis':
        workflow = analysis(opts, log)
    elif opts.rabies_step == 'data_diagnosis':
        workflow = data_diagnosis(opts, log)
    else:
        parser.print_help()

    try:
        log.info('Running workflow with %s plugin.' % opts.plugin)
        # execute workflow, with plugin_args limiting the cluster load for parallel execution
        workflow.run(plugin=opts.plugin, plugin_args={'max_jobs': 50, 'dont_resubmit_completed_jobs': True,
                                                      'n_procs': opts.local_threads, 'qsub_args': '-pe smp %s' % (str(opts.min_proc))})

    except Exception as e:
        log.critical('RABIES failed: %s', e)
        raise


def preprocess(opts, cr_opts, analysis_opts, data_diagnosis_opts, log):
    # Verify input and output directories
    data_dir_path = os.path.abspath(str(opts.bids_dir))
    output_folder = os.path.abspath(str(opts.output_dir))
    if not os.path.isdir(data_dir_path):
        raise ValueError("The provided BIDS data path doesn't exists.")
    else:
        # print the input data directory tree
        log.info("INPUT BIDS DATASET:  \n" + list_files(data_dir_path))

    import SimpleITK as sitk
    if str(opts.data_type) == 'int16':
        opts.data_type = sitk.sitkInt16
    elif str(opts.data_type) == 'int32':
        opts.data_type = sitk.sitkInt32
    elif str(opts.data_type) == 'float32':
        opts.data_type = sitk.sitkFloat32
    elif str(opts.data_type) == 'float64':
        opts.data_type = sitk.sitkFloat64
    else:
        raise ValueError('Invalid --data_type provided.')

    # template options
    # set OS paths to template and atlas files, and convert files to RAS convention if they aren't already
    from rabies.preprocess_pkg.utils import convert_to_RAS
    if not os.path.isfile(opts.anat_template):
        raise ValueError("--anat_template file %s doesn't exists." % (opts.anat_template))
    opts.template_anat = convert_to_RAS(
        str(opts.anat_template), output_folder+'/template_files')

    if not os.path.isfile(opts.brain_mask):
        raise ValueError("--brain_mask file %s doesn't exists." % (opts.brain_mask))
    opts.template_mask = convert_to_RAS(
        str(opts.brain_mask), output_folder+'/template_files')

    if not os.path.isfile(opts.WM_mask):
        raise ValueError("--WM_mask file %s doesn't exists." % (opts.WM_mask))
    opts.WM_mask = convert_to_RAS(
        str(opts.WM_mask), output_folder+'/template_files')

    if not os.path.isfile(opts.CSF_mask):
        raise ValueError("--CSF_mask file %s doesn't exists." % (opts.CSF_mask))
    opts.CSF_mask = convert_to_RAS(
        str(opts.CSF_mask), output_folder+'/template_files')

    if not os.path.isfile(opts.vascular_mask):
        raise ValueError("--vascular_mask file %s doesn't exists." % (opts.vascular_mask))
    opts.vascular_mask = convert_to_RAS(
        str(opts.vascular_mask), output_folder+'/template_files')

    if not os.path.isfile(opts.labels):
        raise ValueError("--labels file %s doesn't exists." % (opts.labels))
    opts.atlas_labels = convert_to_RAS(
        str(opts.labels), output_folder+'/template_files')

    from rabies.main_wf import init_main_wf
    workflow = init_main_wf(data_dir_path, output_folder,
                            opts, cr_opts=cr_opts, analysis_opts=analysis_opts, data_diagnosis_opts=data_diagnosis_opts)

    workflow.base_dir = output_folder

    # setting workflow options for debug mode
    if opts.debug:
        log.setLevel(os.environ.get("LOGLEVEL", "DEBUG"))

        # Change execution parameters
        workflow.config['execution'] = {'stop_on_first_crash': 'true',
                                        'remove_unnecessary_outputs': 'false',
                                        'keep_inputs': 'true'}

        # Change logging parameters
        workflow.config['logging'] = {'workflow_level': 'DEBUG',
                                      'filemanip_level': 'DEBUG',
                                      'interface_level': 'DEBUG',
                                      'utils_level': 'DEBUG'}
        log.debug('Debug ON')

    return workflow


def confound_regression(opts, analysis_opts, data_diagnosis_opts, log):

    cli_file = '%s/rabies_preprocess.pkl' % (opts.preprocess_out, )
    with open(cli_file, 'rb') as handle:
        preprocess_opts = pickle.load(handle)

    workflow = preprocess(preprocess_opts, opts, analysis_opts, data_diagnosis_opts, log)

    return workflow


def analysis(opts, log):

    cli_file = '%s/rabies_confound_regression.pkl' % (
        opts.confound_regression_out, )
    with open(cli_file, 'rb') as handle:
        confound_regression_opts = pickle.load(handle)

    workflow = confound_regression(confound_regression_opts, opts, None, log)

    return workflow


def data_diagnosis(opts, log):

    cli_file = '%s/rabies_confound_regression.pkl' % (
        opts.confound_regression_out, )
    with open(cli_file, 'rb') as handle:
        confound_regression_opts = pickle.load(handle)

    workflow = confound_regression(confound_regression_opts, None, opts, log)

    return workflow


def install_DSURQE(log):

    install=False
    # verifies whether default template files are installed and installs them otherwise
    if not os.path.isfile("%s/DSURQE_40micron_average.nii.gz" % (rabies_path)):
        install=True
    elif not os.path.isfile("%s/DSURQE_40micron_mask.nii.gz" % (rabies_path)):
        install=True
    elif not os.path.isfile("%s/DSURQE_40micron_eroded_WM_mask.nii.gz" % (rabies_path)):
        install=True
    elif not os.path.isfile("%s/DSURQE_40micron_eroded_CSF_mask.nii.gz" % (rabies_path)):
        install=True
    elif not os.path.isfile("%s/DSURQE_40micron_labels.nii.gz" % (rabies_path)):
        install=True
    elif not os.path.isfile("%s/vascular_mask.nii.gz" % (rabies_path)):
        install=True
    elif not os.path.isfile("%s/melodic_IC.nii.gz" % (rabies_path)):
        install=True
    if install:
        from rabies.preprocess_pkg.utils import run_command
        log.info("SOME FILES FROM THE DEFAULT TEMPLATE ARE MISSING. THEY WILL BE INSTALLED BEFORE FURTHER PROCESSING.")
        rc = run_command('install_DSURQE.sh %s' % (rabies_path), verbose = True)


def list_files(startpath):
    string = ''
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        string+='{}{}/ \n'.format(indent, os.path.basename(root))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            string+='{}{} \n'.format(subindent, f)
    return string
