import os
import pickle
import logging
import argparse
from pathlib import Path
import pathos.multiprocessing as multiprocessing  # Better multiprocessing


def get_parser():
    import rabies
    dir_path = os.path.dirname(os.path.realpath(rabies.__file__))

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
        on preprocessing outputs from RABIES. Only selected confound regression and denoising strategies are applied.
        The denoising steps are applied in the following order: ICA-AROMA first, followed by detrending, then
        regression of confound timeseries orthogonal to the application of temporal filters
        (nilearn.clean_img, Lindquist 2018), standardization of timeseries, scrubbing, and finally smoothing.
        """, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    analysis = subparsers.add_parser("analysis",
                                     help="""
        A few built-in resting-state functional connectivity (FC) analysis options are provided to conduct rapid analysis on the cleaned timeseries.
        The options include seed-based FC, voxelwise or parcellated whole-brain FC, group-ICA and dual regression.
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
    preprocess.add_argument("--bias_cor_method", type=str, default='otsu_reg',
                            choices=['otsu_reg', 'thresh_reg'],
                            help="Choose the algorithm for bias field correction of the EPI before registration."
                            "otsu_reg will conduct an initial serie of N4BiasFieldCorrection oriented by Otsu masking method, "
                            "followed by a rigid registration to provide a brain mask orienting the final correction."
                            "thresh_reg will instead use an initial voxel intensity thresholding method for masking, and will "
                            "conduct a subsequent rigid registration to provide a brain mask orienting the final correction.")
    preprocess.add_argument("--disable_anat_preproc", dest='disable_anat_preproc', action='store_true',
                            help="This option disables the preprocessing of anatomical images before commonspace template generation.")
    preprocess.add_argument('--apply_despiking', dest='apply_despiking', action='store_true',
                            help="Whether to apply despiking of the EPI timeseries based on AFNI's "
                            "3dDespike https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html.")
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
        "Options for the registration steps. Built-in options for selecting registration scripts include 'Rigid', 'Affine', 'autoreg_SyN', 'SyN' (non-linear), 'light_SyN', 'multiRAT', but"
        " can specify a custom registration script following the template script structure (see RABIES/rabies/shell_scripts/ for template).")
    g_registration.add_argument("--autoreg", dest='autoreg', action='store_true',
                                help="Choosing this option will conduct an adaptive registration framework which will adjust parameters according to the input images."
                                "This option overrides other registration specifications.")
    g_registration.add_argument("--coreg_script", type=str, default='autoreg_SyN',
                                help="Specify EPI to anat coregistration script.")
    g_registration.add_argument("--anat_reg_script", type=str, default='Rigid',
                                help="specify a registration script for the preprocessing of the anatomical images.")
    g_registration.add_argument("--bias_reg_script", type=str, default='Rigid',
                                help="specify a registration script for iterative bias field correction. This registration step"
                                " consists of aligning the volume with the commonspace template to provide"
                                " a brain mask and optimize the bias field correction.")
    g_registration.add_argument(
        '--template_reg_script',
        type=str,
        default='autoreg_SyN',
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
                         default="%s/../template_files/DSURQE_40micron_average.nii.gz" % (
                             dir_path),
                         help='Anatomical file for the commonspace template.')
    g_atlas.add_argument('--brain_mask', action='store', type=Path,
                         default="%s/../template_files/DSURQE_100micron_mask.nii.gz" % (
                             dir_path),
                         help='Brain mask for the template.')
    g_atlas.add_argument('--WM_mask', action='store', type=Path,
                         default="%s/../template_files/DSURQE_100micron_eroded_WM_mask.nii.gz" % (
                             dir_path),
                         help='White matter mask for the template.')
    g_atlas.add_argument('--CSF_mask', action='store', type=Path,
                         default="%s/../template_files/DSURQE_100micron_eroded_CSF_mask.nii.gz" % (
                             dir_path),
                         help='CSF mask for the template.')
    g_atlas.add_argument('--vascular_mask', action='store', type=Path,
                         default="%s/../template_files/vascular_mask.nii.gz" % (
                             dir_path),
                         help='Can provide a mask of major blood vessels for computing confound timeseries. The default mask was generated by applying MELODIC ICA and selecting the resulting component mapping onto major veins. (Grandjean et al. 2020, NeuroImage; Beckmann et al. 2005)')
    g_atlas.add_argument('--labels', action='store', type=Path,
                         default="%s/../template_files/DSURQE_40micron_labels.nii.gz" % (
                             dir_path),
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
    confound_regression.add_argument('--apply_scrubbing', dest='apply_scrubbing', action='store_true',
                                     default=False,
                                     help="""Whether to apply scrubbing or not. A temporal mask will be generated based on the FD threshold.
                        The frames that exceed the given threshold together with 1 back and 2 forward frames will be masked out
                        from the data after the application of all other confound regression steps (as in Power et al. 2012).""")
    confound_regression.add_argument('--scrubbing_threshold', type=float,
                                     default=0.05,
                                     help='Scrubbing threshold for the mean framewise displacement in mm (averaged across the brain mask) to select corrupted volumes.')
    confound_regression.add_argument('--timeseries_interval', type=str, default='all',
                                     help='Specify a time interval in the timeseries to keep. e.g. "0,80". By default all timeseries are kept.')
    confound_regression.add_argument('--diagnosis_output', dest='diagnosis_output', action='store_true',
                                     default=False,
                                     help="Run a diagnosis for each individual image by computing melodic-ICA on the corrected timeseries,"
                                     "and compute a tSNR map from the input uncorrected image.")

    analysis.add_argument('confound_regression_out', action='store', type=Path,
                          help='path to RABIES confound regression output directory with the datasink.')
    analysis.add_argument('output_dir', action='store', type=Path,
                          help='the output path to drop analysis outputs.')
    analysis.add_argument('--seed_list', type=str,
                                     nargs="*",  # 0 or more values expected => creates a list
                                     default=[],
                                     help="Can provide a list of seed .nii images that will be used to evaluate seed-based correlation maps."
                                     "Each seed must consist of a binary mask representing the ROI in commonspace.")
    g_fc_matrix = analysis.add_argument_group(
        'Options for performing a whole-brain timeseries correlation matrix analysis.')
    g_fc_matrix.add_argument("--FC_matrix", dest='FC_matrix', action='store_true',
                             help="Choose this option to derive a whole-brain functional connectivity matrix, based on the correlation of regional timeseries "
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
                          help="Choose this option to conduct dual regression on each subject cleaned timeseries.")
    g_DR_ICA.add_argument('--IC_file', action='store', type=Path,
                          default=None,
                          help="Option to provide a melodic_IC.nii.gz file with the ICA components from a previous group-ICA run. "
                          "If none is provided, a group-ICA will be run with the dataset cleaned timeseries.")

    return parser


def execute_workflow():
    # generates the parser CLI and execute the workflow based on specified parameters.
    parser = get_parser()
    opts = parser.parse_args()

    output_folder = os.path.abspath(str(opts.output_dir))
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    # options for workflow execution
    if not opts.plugin == 'MultiProc':
        opts.local_threads = 1

    # managing log info
    cli_file = '%s/rabies_%s.pkl' % (output_folder, opts.rabies_step, )
    with open(cli_file, 'wb') as handle:
        pickle.dump(opts, handle, protocol=pickle.HIGHEST_PROTOCOL)

    logging.basicConfig(filename='%s/rabies_%s.log' % (output_folder, opts.rabies_step, ), filemode='w',
                        format='%(asctime)s - %(levelname)s - %(message)s', level=os.environ.get("LOGLEVEL", "INFO"))
    log = logging.getLogger('root')

    from ._info import __version__
    log.info('Running RABIES - version: '+__version__)

    # print complete CLI command
    args = 'CLI INPUTS: \n'
    for arg in vars(opts):
        input = '-> {arg} = {value} \n'.format(
            arg=arg, value=getattr(opts, arg))
        args += input
    log.info(args)

    if opts.rabies_step == 'preprocess':
        workflow = preprocess(opts, None, None, log)
    elif opts.rabies_step == 'confound_regression':
        workflow = confound_regression(opts, None, log)
    elif opts.rabies_step == 'analysis':
        workflow = analysis(opts, log)
    else:
        parser.print_help()

    try:
        print('Running workflow with %s plugin.' % opts.plugin)
        # execute workflow, with plugin_args limiting the cluster load for parallel execution
        workflow.run(plugin=opts.plugin, plugin_args={'max_jobs': 50, 'dont_resubmit_completed_jobs': True,
                                                      'n_procs': opts.local_threads, 'qsub_args': '-pe smp %s' % (str(opts.min_proc))})

    except Exception as e:
        log.critical('RABIES failed: %s', e)
        raise


def preprocess(opts, cr_opts, analysis_opts, log):
    # obtain parser parameters
    data_dir_path = os.path.abspath(str(opts.bids_dir))
    output_folder = os.path.abspath(str(opts.output_dir))

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

    if opts.autoreg:
        opts.bias_reg_script = define_reg_script('Affine')
        opts.anat_reg_script = define_reg_script('Affine')
        opts.coreg_script = define_reg_script('autoreg_SyN')
        opts.template_reg_script = define_reg_script('autoreg_SyN')
    else:
        opts.bias_reg_script = define_reg_script(opts.bias_reg_script)
        opts.anat_reg_script = define_reg_script(opts.anat_reg_script)
        opts.coreg_script = define_reg_script(opts.coreg_script)
        opts.template_reg_script = define_reg_script(opts.template_reg_script)

    # template options
    # set OS paths to template and atlas files, and convert files to RAS convention if they aren't already
    from rabies.preprocess_pkg.utils import convert_to_RAS
    if not os.path.isfile(opts.anat_template):
        raise ValueError("--anat_template file doesn't exists.")
    opts.template_anat = convert_to_RAS(
        str(opts.anat_template), output_folder+'/template_files')

    if not os.path.isfile(opts.brain_mask):
        raise ValueError("--brain_mask file doesn't exists.")
    opts.template_mask = convert_to_RAS(
        str(opts.brain_mask), output_folder+'/template_files')

    if not os.path.isfile(opts.WM_mask):
        raise ValueError("--WM_mask file doesn't exists.")
    opts.WM_mask = convert_to_RAS(
        str(opts.WM_mask), output_folder+'/template_files')

    if not os.path.isfile(opts.CSF_mask):
        raise ValueError("--CSF_mask file doesn't exists.")
    opts.CSF_mask = convert_to_RAS(
        str(opts.CSF_mask), output_folder+'/template_files')

    if not os.path.isfile(opts.vascular_mask):
        raise ValueError("--vascular_mask file doesn't exists.")
    opts.vascular_mask = convert_to_RAS(
        str(opts.vascular_mask), output_folder+'/template_files')

    if not os.path.isfile(opts.labels):
        raise ValueError("--labels file doesn't exists.")
    opts.atlas_labels = convert_to_RAS(
        str(opts.labels), output_folder+'/template_files')

    from rabies.main_wf import init_main_wf
    workflow = init_main_wf(data_dir_path, output_folder,
                            opts, cr_opts=cr_opts, analysis_opts=analysis_opts)

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
        print('Debug ON')

    return workflow


def confound_regression(opts, analysis_opts, log):

    cli_file = '%s/rabies_preprocess.pkl' % (opts.preprocess_out, )
    with open(cli_file, 'rb') as handle:
        preprocess_opts = pickle.load(handle)

    workflow = preprocess(preprocess_opts, opts, analysis_opts, log)

    return workflow


def analysis(opts, log):

    cli_file = '%s/rabies_confound_regression.pkl' % (
        opts.confound_regression_out, )
    with open(cli_file, 'rb') as handle:
        confound_regression_opts = pickle.load(handle)

    workflow = confound_regression(confound_regression_opts, opts, log)

    return workflow


def define_reg_script(reg_option):
    import rabies
    dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
    if reg_option == 'SyN':
        reg_script = dir_path+'/shell_scripts/SyN_registration.sh'
    elif reg_option == 'Affine':
        reg_script = dir_path+'/shell_scripts/antsRegistration_affine.sh'
    elif reg_option == 'autoreg_SyN':
        reg_script = dir_path+'/shell_scripts/antsRegistration_affine_SyN.sh'
    elif reg_option == 'light_SyN':
        reg_script = dir_path+'/shell_scripts/light_SyN_registration.sh'
    elif reg_option == 'Rigid':
        reg_script = dir_path+'/shell_scripts/antsRegistration_rigid.sh'
    elif reg_option == 'multiRAT':
        reg_script = dir_path+'/shell_scripts/multiRAT_registration.sh'
    else:
        '''
        For user-provided antsRegistration command.
        '''
        if os.path.isfile(reg_option):
            reg_script = reg_option
        else:
            raise ValueError(
                'REGISTRATION ERROR: THE REG SCRIPT FILE DOES NOT EXISTS')
    return reg_script
