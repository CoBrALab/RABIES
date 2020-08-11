import os
import sys

import pickle
import logging
import argparse
from pathlib import Path
import pathos.multiprocessing as multiprocessing  # Better multiprocessing

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
        registration to a commonspace atlas and associated masks, as well as further options (see --help).
        """)
    confound_regression = subparsers.add_parser("confound_regression",
        help="""Flexible options for confound regression applied
        on preprocessing outputs from RABIES. Smoothing is applied first, followed by ICA-AROMA, detrending, then
        regression of confound timeseries orthogonal to the application of temporal filters
        (nilearn.clean_img, Lindquist 2018), standardization of timeseries and finally scrubbing. The corrections
        follow user specifications.
        """)
    analysis = subparsers.add_parser("analysis",
        help="""
        Optional analysis to conduct on cleaned timeseries.
        """)

    g_execution = parser.add_argument_group("Options for managing the execution of the workflow.")
    g_execution.add_argument("-p", "--plugin", type=str, default='Linear',
                        help="Specify the nipype plugin for workflow execution. Consult nipype plugin documentation for detailed options."
                             " Linear, MultiProc, SGE and SGEGraph have been tested.")
    g_execution.add_argument('--local_threads',type=int,default=multiprocessing.cpu_count(),
        help="""For local MultiProc execution, set the maximum number of processors run in parallel,
        defaults to number of CPUs. This option only applies to the MultiProc execution plugin, otherwise
        it is set to 1.""")
    g_execution.add_argument("--scale_min_memory", type=float, default=1.0,
                        help="For a parallel execution with MultiProc, the minimal memory attributed to nodes can be scaled with this multiplier to avoid memory crashes.")
    g_execution.add_argument("--min_proc", type=int, default=1,
                        help="For SGE parallel processing, specify the minimal number of nodes to be assigned to avoid memory crashes.")
    g_execution.add_argument("--data_type", type=str, default='float32',
                        help="Specify data format outputs to control for file size among 'int16','int32','float32' and 'float64'.")
    g_execution.add_argument("--debug", dest='debug', action='store_true',
                        help="Run in debug mode.")

    preprocess.add_argument('bids_dir', action='store', type=Path,
                        help='the root folder of the BIDS-formated input data directory.')
    preprocess.add_argument('output_dir', action='store', type=Path,
                        help='the output path to drop outputs from major preprocessing steps.')
    preprocess.add_argument("-e", "--bold_only", dest='bold_only', action='store_true',
                        help="Apply preprocessing with only EPI scans. commonspace registration"
                              " is executed through registration of the EPI-generated template from ants_dbm"
                              " to the anatomical template.")
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

    g_registration = preprocess.add_argument_group("Options for the registration steps.")
    g_registration.add_argument("--autoreg", dest='autoreg', action='store_true',
                        help="Choosing this option will conduct an adaptive registration framework which will adjust parameters according to the input images."
                        "This option overrides other registration specifications.")
    g_registration.add_argument("-r", "--coreg_script", type=str, default='light_SyN',
                        help="Specify EPI to anat coregistration script. Built-in options include 'Rigid', 'Affine', 'autoreg_affine', 'autoreg_SyN', 'SyN' (non-linear), 'light_SyN', but"
                        " can specify a custom registration script following the template script structure (see RABIES/rabies/shell_scripts/ for template).")
    g_registration.add_argument("--bias_reg_script", type=str, default='Rigid',
                        help="specify a registration script for iterative bias field correction. This registration step"
                        " consists of aligning the volume with the commonspace template to provide"
                        " a brain mask and optimize the bias field correction. The registration script options are the same as --coreg_script.")
    g_registration.add_argument(
        '--template_reg_script',
        type=str,
        default='light_SyN',
        help="""Registration script that will be used for registration of the generated dataset
        template to the provided commonspace atlas for masking and labeling. Can choose a predefined
        registration script among Rigid,Affine,SyN or light_SyN, or provide a custom script.""")

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
        '--anatomical_resampling',type=str,default='inputs_defined',
        help="""To optimize the efficiency of registration, the provided anatomical template is resampled based on the provided
        input images. The dimension with the lowest resolution among the provided anatomical images (EPI images instead if --bold_only is True)
        is selected as a basis for resampling the template to isotropic resolution, if the provided resolution is lower than the original
        resolution of the template. Alternatively, the user can provide a custom resampling dimension. This allows to accelerate
        registration steps with minimal sampling dimensions.""")

    g_ants_dbm = preprocess.add_argument_group('cluster options for running ants_dbm (options copied from twolevel_dbm.py):')
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
    g_ants_dbm.add_argument(
        '--memory_request',
        default="8gb",
        help="""Option for job submission
        specifying requested memory per pairwise registration.""")

    g_stc = preprocess.add_argument_group("""Specify Slice Timing Correction info that is fed to AFNI 3dTshift
    (https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html). The STC is applied in the
    anterior-posterior orientation, assuming slices were acquired in this direction.""")
    g_stc.add_argument('--TR', type=str, default='1.0s',
                        help="Specify repetition time (TR) in seconds.")
    g_stc.add_argument('--no_STC', dest='no_STC', action='store_true',
                        help="Select this option to ignore the STC step.")
    g_stc.add_argument('--tpattern', type=str, default='alt',
                        help="Specify if interleaved or sequential acquisition. 'alt' for interleaved, 'seq' for sequential.")

    g_atlas = preprocess.add_argument_group('Provided commonspace atlas files.')
    g_atlas.add_argument('--anat_template', action='store', type=Path,
                        default="%s/template_files/DSURQE_40micron_average.nii.gz" % (os.environ["RABIES"]),
                        help='Anatomical file for the commonspace template.')
    g_atlas.add_argument('--brain_mask', action='store', type=Path,
                        default="%s/template_files/DSURQE_100micron_mask.nii.gz" % (os.environ["RABIES"]),
                        help='Brain mask for the template.')
    g_atlas.add_argument('--WM_mask', action='store', type=Path,
                        default="%s/template_files/DSURQE_100micron_eroded_WM_mask.nii.gz" % (os.environ["RABIES"]),
                        help='White matter mask for the template.')
    g_atlas.add_argument('--CSF_mask', action='store', type=Path,
                        default="%s/template_files/DSURQE_100micron_eroded_CSF_mask.nii.gz" % (os.environ["RABIES"]),
                        help='CSF mask for the template.')
    g_atlas.add_argument('--vascular_mask', action='store', type=Path,
                        default="%s/template_files/vascular_mask.nii.gz" % (os.environ["RABIES"]),
                        help='Can provide a mask of major blood vessels for computing confound timeseries. The default mask was generated by applying MELODIC ICA and selecting the resulting component mapping onto major veins. (Grandjean et al. 2020, NeuroImage; Beckmann et al. 2005)')
    g_atlas.add_argument('--labels', action='store', type=Path,
                        default="%s/template_files/DSURQE_40micron_labels.nii.gz" % (os.environ["RABIES"]),
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
    confound_regression.add_argument('--smoothing_filter', type=float, default=0.3,
                        help='Specify smoothing filter size in mm.')
    confound_regression.add_argument('--run_aroma', dest='run_aroma', action='store_true',
                        default=False,
                        help='Whether to run ICA AROMA or not.')
    confound_regression.add_argument('--aroma_dim', type=int,
                        default=0,
                        help='Can specify a number of dimension for MELODIC.')
    confound_regression.add_argument('--conf_list', type=str,
                        nargs="*",  # 0 or more values expected => creates a list
                        default=[],
                        help='list of regressors. Possible options: WM_signal,CSF_signal,vascular_signal,aCompCor,global_signal,mot_6,mot_24, mean_FD')
    confound_regression.add_argument('--apply_scrubbing', dest='apply_scrubbing', action='store_true',
                        default=False,
                        help="""Whether to apply scrubbing or not. A temporal mask will be generated based on the FD threshold.
                        The frames that exceed the given threshold together with 1 back and 2 forward frames will be masked out
                        from the data after the application of all other confound regression steps (as in Power et al. 2012).""")
    confound_regression.add_argument('--scrubbing_threshold', type=float,
                        default=0.1,
                        help='Scrubbing threshold for the mean framewise displacement in mm? (averaged across the brain mask) to select corrupted volumes.')
    confound_regression.add_argument('--timeseries_interval', type=str, default='all',
                        help='Specify a time interval in the timeseries to keep. e.g. "0,80". By default all timeseries are kept.')
    confound_regression.add_argument('--diagnosis_output', dest='diagnosis_output', action='store_true',
                        default=False,
                        help="Run a diagnosis for each image by computing melodic-ICA on the corrected timeseries,"
                             "and compute a tSNR map from the input uncorrected image.")
    confound_regression.add_argument('--seed_list', type=str,
                        nargs="*",  # 0 or more values expected => creates a list
                        default=[],
                        help='Can provide a list of seed .nii images that will be used to evaluate seed-based correlation maps during data diagnosis.')

    analysis.add_argument('confound_regression_out', action='store', type=Path,
                        help='path to RABIES confound regression output directory with the datasink.')
    analysis.add_argument('output_dir', action='store', type=Path,
                        help='the output path to drop analysis outputs.')

    return parser


def execute_workflow():
    #generates the parser CLI and execute the workflow based on specified parameters.
    parser = get_parser()
    opts = parser.parse_args()

    output_folder=os.path.abspath(str(opts.output_dir))
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    #options for workflow execution
    plugin=opts.plugin
    os.environ["min_proc"]=str(opts.min_proc)
    if plugin=='MultiProc':
        os.environ["local_threads"]=str(opts.local_threads)
    else:
        os.environ["local_threads"]='1'
    os.environ["rabies_mem_scale"]=str(opts.scale_min_memory)
    import SimpleITK as sitk
    if str(opts.data_type)=='int16':
        os.environ["rabies_data_type"]=str(sitk.sitkInt16)
    elif str(opts.data_type)=='int32':
        os.environ["rabies_data_type"]=str(sitk.sitkInt32)
    elif str(opts.data_type)=='float32':
        os.environ["rabies_data_type"]=str(sitk.sitkFloat32)
    elif str(opts.data_type)=='float64':
        os.environ["rabies_data_type"]=str(sitk.sitkFloat64)
    else:
        raise ValueError('Invalid --data_type provided.')

    ###managing log info
    cli_file='%s/rabies_%s.pkl' % (output_folder, opts.rabies_step, )
    with open(cli_file, 'wb') as handle:
        pickle.dump(opts, handle, protocol=pickle.HIGHEST_PROTOCOL)

    logging.basicConfig(filename='%s/rabies_%s.log' % (output_folder, opts.rabies_step, ), filemode='w', format='%(asctime)s - %(levelname)s - %(message)s', level=os.environ.get("LOGLEVEL", "INFO"))
    log = logging.getLogger(__name__)

    from ._info import __version__
    log.info('Running RABIES - version: '+__version__)

    #print complete CLI command
    args='CLI INPUTS: \n'
    for arg in vars(opts):
        input='-> {arg} = {value} \n'.format(
            arg=arg, value=getattr(opts, arg))
        args+=input
    log.info(args)

    if opts.rabies_step == 'preprocess':
        workflow = preprocess(opts,None,None)
    elif opts.rabies_step == 'confound_regression':
        workflow = confound_regression(opts,None)
    elif opts.rabies_step == 'analysis':
        workflow = analysis(opts)
    else:
        parser.print_help()

    #setting workflow options for debug mode
    if opts.debug:
        # Change execution parameters
        workflow.config['execution'] = {'stop_on_first_crash' : 'true',
                                'remove_unnecessary_outputs': 'false',
                                'keep_inputs': 'true',
                                'log_directory' : os.getcwd()}

        # Change logging parameters
        workflow.config['logging'] = {'workflow_level' : 'DEBUG',
                                'filemanip_level' : 'DEBUG',
                                'interface_level' : 'DEBUG',
                                'utils_level' : 'DEBUG',
                                'log_to_file' : 'True',
                                'log_directory' : os.getcwd()}
        print('Debug ON')

    try:
        print('Running workflow with %s plugin.' % plugin)
        #execute workflow, with plugin_args limiting the cluster load for parallel execution
        workflow.run(plugin=plugin, plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'n_procs' : int(os.environ["local_threads"]), 'qsub_args': '-pe smp %s' % (os.environ["min_proc"])})

    except Exception as e:
        log.critical('RABIES failed: %s', e)
        raise


def preprocess(opts, cr_opts, analysis_opts):
    #obtain parser parameters
    output_folder=os.path.abspath(str(opts.output_dir))
    bold_preproc_only=opts.bold_only
    disable_anat_preproc=opts.disable_anat_preproc
    data_dir_path=os.path.abspath(str(opts.bids_dir))
    detect_dummy=opts.detect_dummy
    apply_despiking=opts.apply_despiking
    apply_slice_mc=opts.apply_slice_mc

    if opts.autoreg:
        bias_reg_script=define_reg_script('autoreg_affine')
    else:
        bias_reg_script=define_reg_script(opts.bias_reg_script)
    if opts.autoreg:
        coreg_script=define_reg_script('autoreg_SyN')
    else:
        coreg_script=define_reg_script(opts.coreg_script)
    if opts.autoreg:
        template_reg_script=define_reg_script('autoreg_SyN')
    else:
        template_reg_script=define_reg_script(opts.template_reg_script)

    #STC options
    stc_bool=not opts.no_STC
    stc_TR=opts.TR
    if opts.tpattern=="alt":
        stc_tpattern='alt-z'
    elif opts.tpattern=="seq":
        stc_tpattern='seq-z'
    else:
        raise ValueError('Invalid --tpattern provided.')

    #resampling options
    nativespace_resampling=opts.nativespace_resampling
    commonspace_resampling=opts.commonspace_resampling
    os.environ["anatomical_resampling"]=opts.anatomical_resampling

    #setting absolute paths for ants_dbm options options
    os.environ["ants_dbm_cluster_type"]=opts.cluster_type
    os.environ["ants_dbm_walltime"]=opts.walltime
    os.environ["ants_dbm_memory_request"]=opts.memory_request

    #template options
    # set OS paths to template and atlas files, and convert files to RAS convention if they aren't already
    from rabies.preprocess_pkg.utils import convert_to_RAS
    if not os.path.isfile(str(opts.anat_template)):
        raise ValueError("--anat_template file doesn't exists.")
    os.environ["template_anat"] = convert_to_RAS(str(opts.anat_template), os.environ["RABIES"]+'/template_files')

    if not os.path.isfile(str(opts.brain_mask)):
        raise ValueError("--brain_mask file doesn't exists.")
    os.environ["template_mask"] = convert_to_RAS(str(opts.brain_mask), os.environ["RABIES"]+'/template_files')

    if not os.path.isfile(str(opts.WM_mask)):
        raise ValueError("--WM_mask file doesn't exists.")
    os.environ["WM_mask"] = convert_to_RAS(str(opts.WM_mask), os.environ["RABIES"]+'/template_files')

    if not os.path.isfile(str(opts.CSF_mask)):
        raise ValueError("--CSF_mask file doesn't exists.")
    os.environ["CSF_mask"] = convert_to_RAS(str(opts.CSF_mask), os.environ["RABIES"]+'/template_files')

    if not os.path.isfile(str(opts.vascular_mask)):
        raise ValueError("--vascular_mask file doesn't exists.")
    os.environ["vascular_mask"] = convert_to_RAS(str(opts.vascular_mask), os.environ["RABIES"]+'/template_files')

    if not os.path.isfile(str(opts.labels)):
        raise ValueError("--labels file doesn't exists.")
    os.environ["atlas_labels"] = convert_to_RAS(str(opts.labels), os.environ["RABIES"]+'/template_files')

    if bold_preproc_only:
        from rabies.preprocess_pkg.bold_main_wf import init_EPIonly_bold_main_wf
        workflow = init_EPIonly_bold_main_wf(data_dir_path, output_folder, apply_despiking=apply_despiking, tr=stc_TR,
            tpattern=stc_tpattern, apply_STC=stc_bool, detect_dummy=detect_dummy, slice_mc=apply_slice_mc,
            bias_reg_script=bias_reg_script, coreg_script=coreg_script, template_reg_script=template_reg_script, commonspace_resampling=commonspace_resampling)
    elif not bold_preproc_only:
        from rabies.main_wf import init_unified_main_wf
        workflow = init_unified_main_wf(data_dir_path, output_folder, disable_anat_preproc=disable_anat_preproc, autoreg=opts.autoreg, apply_despiking=apply_despiking, tr=stc_TR,
            tpattern=stc_tpattern, detect_dummy=detect_dummy, slice_mc=apply_slice_mc, template_reg_script=template_reg_script, apply_STC=stc_bool,
            bias_reg_script=bias_reg_script, coreg_script=coreg_script, nativespace_resampling=nativespace_resampling,
            commonspace_resampling=commonspace_resampling, cr_opts=cr_opts, analysis_opts=analysis_opts)
    else:
        raise ValueError('bold_preproc_only must be true or false.')

    workflow.base_dir = output_folder

    return workflow

def confound_regression(opts, analysis_opts):

    cli_file='%s/rabies_preprocess.pkl' % (opts.preprocess_out, )
    with open(cli_file, 'rb') as handle:
        preprocess_opts=pickle.load(handle)

    workflow = preprocess(preprocess_opts, opts, analysis_opts)

    return workflow

def analysis(opts):

    cli_file='%s/rabies_confound_regression.pkl' % (opts.confound_regression_out, )
    with open(cli_file, 'rb') as handle:
        confound_regression_opts=pickle.load(handle)

    workflow = confound_regression(confound_regression_opts, opts)

    return workflow

def define_reg_script(reg_option):
    import rabies
    dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
    if reg_option=='SyN':
        reg_script=dir_path+'/shell_scripts/SyN_registration.sh'
    elif reg_option=='autoreg_affine':
        reg_script=dir_path+'/shell_scripts/antsRegistration_affine.sh'
    elif reg_option=='autoreg_SyN':
        reg_script=dir_path+'/shell_scripts/antsRegistration_affine_SyN.sh'
    elif reg_option=='light_SyN':
        reg_script=dir_path+'/shell_scripts/light_SyN_registration.sh'
    elif reg_option=='Affine':
        reg_script=dir_path+'/shell_scripts/Affine_registration.sh'
    elif reg_option=='Rigid':
        reg_script=dir_path+'/shell_scripts/Rigid_registration.sh'
    else:
        '''
        For user-provided antsRegistration command.
        '''
        if os.path.isfile(reg_option):
            reg_script=reg_option
        else:
            raise ValueError('REGISTRATION ERROR: THE REG SCRIPT FILE DOES NOT EXISTS')
    return reg_script

'''
def eval_memory(data_dir_path,native_spacing,commonspace_spacing,anat_resampling_spacing,scale_memory):
    import os
    import numpy as np
    import SimpleITK as sitk
    from bids.layout import BIDSLayout
    layout = BIDSLayout(data_dir_path, validate=False)
    #subject_list, session_iter, run_iter=prep_bids_iter(layout)
    func_list=layout.get(datatype='func', extension=['nii', 'nii.gz'], return_type='filename')
    anat_list=layout.get(datatype='anat', extension=['nii', 'nii.gz'], return_type='filename')

    #define the anatomical resampling
    if anat_resampling_spacing=='inputs_defined':
        img = sitk.ReadImage(anat_list[0], int(os.environ["rabies_data_type"]))
        low_dim=np.asarray(img.GetSpacing()[:3]).min()
        for file in anat_list[1:]:
            img = sitk.ReadImage(file, int(os.environ["rabies_data_type"]))
            new_low_dim=np.asarray(img.GetSpacing()[:3]).min()
            if new_low_dim<low_dim:
                low_dim=new_low_dim
        anat_resampling_spacing=(low_dim,low_dim,low_dim)
    else:
        shape=anat_resampling_spacing.split('x')
        anat_resampling_spacing=(float(shape[0]),float(shape[1]),float(shape[2]))

    def size_ratio(input_size,spacing,output_spacing):
        #get the resampled size
        sampling_ratio=np.asarray(spacing)/np.asarray(output_spacing)
        output_size = [int(input_size[0]*sampling_ratio[0]), int(input_size[1]*sampling_ratio[1]), int(input_size[2]*sampling_ratio[2])]
        size_ratio=np.array(output_size).prod()/np.array(input_size).prod()
        return size_ratio

    #determine the file sizes
    EPI_mb=0
    for file in func_list:
        image = sitk.ReadImage(file)
        bitpix=int(image.GetMetaData('bitpix'))
        file_size64 = os.path.getsize(file)*1e-6*64/bitpix #the size is scaled to a 64bitpix
        if file_size64>EPI_mb:
            EPI_mb=file_size64
            EPI_size=image.GetSize()
            EPI_spacing=image.GetSpacing()

    #if EPI_mb<100:
    #    EPI_mb=100

    anat_mb=0
    for file in anat_list:
        image = sitk.ReadImage(file)
        bitpix=int(image.GetMetaData('bitpix'))
        file_size64 = os.path.getsize(file)*1e-6*64/bitpix #the size is scaled to a 64bitpix
        if file_size64>anat_mb:
            anat_mb=file_size64
            anat_size=image.GetSize()
            anat_spacing=image.GetSpacing()

    #if anat_mb<5:
    #    anat_mb=5

    anat_resampled_mb=anat_mb*size_ratio(anat_size,anat_spacing,anat_resampling_spacing)
    if native_spacing=='origin':
        EPI_native_mb=EPI_mb
    else:
        shape=native_spacing.split('x')
        native_spacing=(float(shape[0]),float(shape[1]),float(shape[2]))
        EPI_native_mb=EPI_mb*size_ratio(EPI_size[:3],EPI_spacing[:3],native_spacing)
    if commonspace_spacing=='origin':
        EPI_commonspace_mb=EPI_mb
    else:
        shape=commonspace_spacing.split('x')
        commonspace_spacing=(float(shape[0]),float(shape[1]),float(shape[2]))
        EPI_commonspace_mb=EPI_mb*size_ratio(EPI_size[:3],EPI_spacing[:3],commonspace_spacing)
    low_dim=np.array(EPI_spacing[:3]).min()
    output_spacing=(low_dim,low_dim,low_dim)
    input_size=EPI_size[:3]
    spacing=EPI_spacing[:3]
    sampling_ratio=np.asarray(spacing)/np.asarray(output_spacing)
    output_size = [int(input_size[0]*sampling_ratio[0]), int(input_size[1]*sampling_ratio[1]), int(input_size[2]*sampling_ratio[2])]
    ratio=np.array(output_size).prod()/np.array(anat_size).prod()
    EPI_biascor_mb=anat_mb*ratio

    os.environ["anat_resampled_gb"]=str(scale_memory*round(anat_resampled_mb,2)/1000)
    os.environ["EPI_gb"]=str(scale_memory*round(EPI_mb,2)/1000)
    os.environ["EPI_native_gb"]=str(scale_memory*round(EPI_native_mb,2)/1000)
    os.environ["EPI_commonspace_gb"]=str(scale_memory*round(EPI_commonspace_mb,2)/1000)
    os.environ["EPI_biascor_gb"]=str(scale_memory*round(EPI_biascor_mb,2)/1000)

    print('STC mem '+os.environ["EPI_gb"])
    print('HMC mem '+os.environ["EPI_gb"])
    print('gen_ref mem '+os.environ["EPI_gb"])
    print('anat_preproc mem '+os.environ["anat_resampled_gb"])

    return {'anat_mb':round(anat_mb,2),'anat_resampled_mb':round(anat_resampled_mb,2),'EPI_mb':round(EPI_mb,2),'EPI_native_mb':round(EPI_native_mb,2),
            'EPI_commonspace_mb':round(EPI_commonspace_mb,2),'EPI_biascor_mb':round(EPI_biascor_mb,2)}
'''
