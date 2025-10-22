import os
import pickle
import SimpleITK as sitk
from nipype import logging, config
from .boilerplate import *
from .parser import get_parser,read_parser

# setting all default template files
if 'XDG_DATA_HOME' in os.environ.keys():
    rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
else:
    rabies_path = os.environ['HOME']+'/.local/share/rabies'

DSURQE_ANAT=f"{rabies_path}/DSURQE_40micron_average.nii.gz"
DSURQE_MASK=f"{rabies_path}/DSURQE_40micron_mask.nii.gz"
DSURQE_WM=f"{rabies_path}/DSURQE_40micron_eroded_WM_mask.nii.gz"
DSURQE_CSF=f"{rabies_path}/DSURQE_40micron_eroded_CSF_mask.nii.gz"
DSURQE_VASC=f"{rabies_path}/vascular_mask.nii.gz"
DSURQE_LABELS=f"{rabies_path}/DSURQE_40micron_labels.nii.gz"
DSURQE_MAPPINGS=f"{rabies_path}/DSURQE_40micron_R_mapping.csv"
DSURQE_ICA=f"{rabies_path}/melodic_IC.nii.gz"

EPICOMMON_ANAT=f"{rabies_path}/EPI_template.nii.gz"
EPICOMMON_MASK=f"{rabies_path}/EPI_brain_mask.nii.gz"
EPICOMMON_WM=f"{rabies_path}/EPI_WM_mask.nii.gz"
EPICOMMON_CSF=f"{rabies_path}/EPI_CSF_mask.nii.gz"
EPICOMMON_VASC=f"{rabies_path}/EPI_vascular_mask.nii.gz"
EPICOMMON_LABELS=f"{rabies_path}/EPI_labels.nii.gz"
EPICOMMON_ICA=f"{rabies_path}/melodic_IC_resampled.nii.gz"


def execute_workflow(args=None, return_workflow=False):
    # generates the parser CLI and execute the workflow based on specified parameters.
    parser = get_parser()
    opts = read_parser(parser, args)

    try: # convert the output path to absolute if not already the case
        opts.output_dir = os.path.abspath(str(opts.output_dir))
    except:
        parser.print_help()
        return

    if not os.path.isdir(opts.output_dir):
        os.makedirs(opts.output_dir)

    log = prep_logging(opts, opts.output_dir)

    # verify default template installation
    install_DSURQE(log)

    from .__version__ import __version__
    log.info('Running RABIES - version: '+__version__)

    # print complete CLI command
    args = 'CLI INPUTS: \n'
    for arg in vars(opts):
        input = f'-> {arg} = {getattr(opts, arg)} \n'
        args += input
    log.info(args)

    # inclusion/exclusion list are incompatible parameters
    if (not opts.inclusion_ids[0]=='all') and (not opts.exclusion_ids[0]=='none'):
        raise ValueError(f"""
           Either an inclusion list (--inclusion_ids) or exclusion list (--exclusion_ids)
           can be provided, not both.
           """)

    if opts.rabies_stage == 'preprocess':
        workflow = preprocess(opts, log)
    elif opts.rabies_stage == 'confound_correction':
        workflow = confound_correction(opts, log)
    elif opts.rabies_stage == 'analysis':
        workflow = analysis(opts, log)
    else:
        parser.print_help()
    workflow.base_dir = opts.output_dir

    # the cli parameters are saved after workflow has been prepared, since they have to be modified during workflow preparation
    cli_file = f'{opts.output_dir}/rabies_{opts.rabies_stage}.pkl'
    with open(cli_file, 'wb') as handle:
        pickle.dump(opts, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if return_workflow:
        return workflow
    try:
        log.info(f'Running workflow with {opts.plugin} plugin.')
        # execute workflow, with plugin_args limiting the cluster load for parallel execution
        graph_out = workflow.run(plugin=opts.plugin, plugin_args={'max_jobs': 50, 'dont_resubmit_completed_jobs': True,
                                                      'n_procs': opts.local_threads, 'qsub_args': f'-pe smp {str(opts.min_proc)}'})
        # save the workflow execution
        workflow_file = f'{opts.output_dir}/rabies_{opts.rabies_stage}_workflow.pkl'
        with open(workflow_file, 'wb') as handle:
            pickle.dump(graph_out, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    except Exception as e:
        log.critical(f'RABIES failed: {e}')
        raise


def prep_logging(opts, output_folder):
    cli_file = f'{output_folder}/rabies_{opts.rabies_stage}.pkl'
    if os.path.isfile(cli_file) and not opts.force:
        raise ValueError(f"""
            A previous run was indicated by the presence of {cli_file}.
            This can lead to inconsistencies between previous outputs and the log files.
            To prevent this, we recommend removing previous datasinks from the {opts.rabies_stage} 
            RABIES stage. To continue with your execution, the {cli_file} file must be  
            removed (use --force to automatically do so).
            """)

    # remove old versions of the log if already existing
    log_path = f'{output_folder}/rabies_{opts.rabies_stage}.log'
    if os.path.isfile(log_path):
        os.remove(log_path)

    config.update_config({'logging': {'log_directory': output_folder,
                                    'log_to_file': True}})

    # setting workflow logging level
    if opts.verbose==0:
        level="WARNING"
    elif opts.verbose==1:
        level="INFO"
    elif opts.verbose>=2:
        level="DEBUG"
        config.enable_debug_mode()
    else:
        raise ValueError(f"--verbose must be provided an integer of 0 or above. {opts.verbose} was provided instead.")

    # nipype has hard-coded 'nipype.log' filename; we rename it after it is created, and change the handlers
    logging.update_logging(config)
    os.rename(f'{output_folder}/pypeline.log', log_path)
    # change the handlers path to the desired file
    for logger in logging.loggers.keys():
        log = logging.getLogger(logger)
        handler = log.handlers[0]
        handler.baseFilename = log_path

    # set the defined level of verbose
    log = logging.getLogger('nipype.workflow')
    log.setLevel(level)
    log.debug('Debug ON')
    return log


def preprocess(opts, log):
    # convert the input path to absolute if not already the case
    opts.bids_dir = os.path.abspath(str(opts.bids_dir))

    if not os.path.isdir(opts.bids_dir):
        raise ValueError("The provided BIDS data path doesn't exists.")
    else:
        # print the input data directory tree
        log.info("INPUT BIDS DATASET:  \n" + list_files(opts.bids_dir))

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
    
    # if the default template is not used, then brain mask input is required, 
    # and other optional files are set to None if no input was provided 
    # to block downstream operations dependent on those inputs, but allow preprocessing nevertheless
    if not str(opts.anat_template)==DSURQE_ANAT:
        if str(opts.brain_mask)==DSURQE_MASK:
            raise ValueError("The default anatomical template was changed, but not the brain mask "
                             "- it is necessary to provide a new brain mask matching the template.")
        # make sure we have absolute paths
        opts.anat_template = os.path.abspath(opts.anat_template)
        opts.brain_mask = os.path.abspath(opts.brain_mask)
        
        for opt_key,default_file in zip(['WM_mask','CSF_mask','vascular_mask','labels'],
                                         [DSURQE_WM,DSURQE_CSF,DSURQE_VASC,DSURQE_LABELS]):
            
            opt_file = getattr(opts, opt_key)
            if str(opt_file)==default_file:
                opt_file=None
            else:
                opt_file = os.path.abspath(opt_file) # make sure we have absolute paths
            setattr(opts, opt_key, opt_file)

    # if --bold_only, the default atlas files change to EPI versions
    if opts.bold_only:
        for opt_key,default_file,EPI_file in zip(['anat_template','brain_mask','WM_mask','CSF_mask','vascular_mask','labels'],
                                         [DSURQE_ANAT, DSURQE_MASK, DSURQE_WM,DSURQE_CSF,DSURQE_VASC,DSURQE_LABELS],
                                         [EPICOMMON_ANAT, EPICOMMON_MASK, EPICOMMON_WM,EPICOMMON_CSF,EPICOMMON_VASC,EPICOMMON_LABELS]):
            opt_file = getattr(opts, opt_key)
            if str(opt_file)==default_file:
                setattr(opts, opt_key, EPI_file)
                log.info(f'With --bold_only, default --{opt_key} changed to {EPI_file}')

    # final check of template file formats
    from rabies.preprocess_pkg.utils import convert_to_RAS
    for opt_key,check_binary in zip(['anat_template', 'brain_mask', 'WM_mask','CSF_mask','vascular_mask','labels'],
                                    [False,True,True,True,True,False]):
        opt_file = getattr(opts, opt_key)
        if opt_file is not None: # some masks might be set to None
            if not os.path.isfile(opt_file):
                raise ValueError(f"--{opt_key} file {opt_file} doesn't exists.")
            opt_file = convert_to_RAS(
                str(opt_file), opts.output_dir+'/template_files')

            if check_binary:
                check_binary_masks(opt_file)
            check_template_overlap(opts.anat_template, opt_file)
            setattr(opts, opt_key, opt_file)

    if not opts.inherit_unbiased_template=='none':
        opts.inherit_unbiased_template = os.path.abspath(opts.inherit_unbiased_template)
        if not os.path.isdir(opts.inherit_unbiased_template):
            raise ValueError(f"--inherit_unbiased_template path {opts.inherit_unbiased_template} doesn't exists.")

    check_resampling_syntax(opts.nativespace_resampling)
    check_resampling_syntax(opts.commonspace_resampling)
    check_resampling_syntax(opts.anatomical_resampling)

    # write boilerplate
    boilerplate_file = f'{opts.output_dir}/boilerplate.txt'

    methods,ref_string = preprocess_boilerplate(opts)
    txt_boilerplate="#######PREPROCESSING\n\n"+methods+ref_string+'\n\n'
    with open(boilerplate_file, "w") as text_file:
        text_file.write(txt_boilerplate)

    from rabies.preprocess_pkg.main_wf import init_main_wf
    workflow = init_main_wf(opts.bids_dir, opts.output_dir, opts)

    return workflow


def confound_correction(opts, log):

    if opts.edge_cutoff == 0 and (opts.highpass is not None):
        log.warning(
            "\n############################################# WARNING\n"
            "Highpass filtering will be applied without removing timepoints at each edge "
            "of acquisition. This may introduce edge artefacts. We recommend removing "
            "~30sec at both end of the acquisition for a filter of 0.01Hz."
            "\n############################################# WARNING\n")

    cli_file = f'{opts.preprocess_out}/rabies_preprocess.pkl'
    with open(cli_file, 'rb') as handle:
        preprocess_opts = pickle.load(handle)

    boilerplate_file = f'{opts.output_dir}/boilerplate_confound_correction.txt'
    methods,ref_string = confound_correction_boilerplate(opts)
    txt_boilerplate="#######CONFOUND CORRECTION\n\n"+methods+ref_string+'\n\n'
    with open(boilerplate_file, "w") as text_file:
        text_file.write(txt_boilerplate)

    from rabies.confound_correction_pkg.main_wf import init_main_confound_correction_wf
    workflow = init_main_confound_correction_wf(preprocess_opts, opts)

    return workflow


def analysis(opts, log):

    cli_file = f'{opts.confound_correction_out}/rabies_confound_correction.pkl'
    with open(cli_file, 'rb') as handle:
        confound_correction_opts = pickle.load(handle)

    cli_file = f'{confound_correction_opts.preprocess_out}/rabies_preprocess.pkl'
    with open(cli_file, 'rb') as handle:
        preprocess_opts = pickle.load(handle)

    if preprocess_opts.bold_only:
        if str(opts.prior_maps)==DSURQE_ICA:
            file=EPICOMMON_ICA
            opts.prior_maps=file
            log.info('With --bold_only, default --prior_maps changed to '+file)


    from rabies.analysis_pkg.main_wf import init_main_analysis_wf
    workflow = init_main_analysis_wf(preprocess_opts, confound_correction_opts, opts)

    return workflow

def install_DSURQE(log):

    install = False
    # verifies whether default template files are installed and installs them otherwise
    for f in [DSURQE_ANAT,DSURQE_MASK,DSURQE_WM,DSURQE_CSF,DSURQE_VASC,DSURQE_LABELS,DSURQE_MAPPINGS,DSURQE_ICA,
              EPICOMMON_ANAT,EPICOMMON_MASK,EPICOMMON_WM,EPICOMMON_CSF,EPICOMMON_VASC,EPICOMMON_LABELS,EPICOMMON_ICA]:
        if not os.path.isfile(f):
            install = True

    if install:
        from rabies.preprocess_pkg.utils import run_command
        log.info(
            "SOME FILES FROM THE DEFAULT TEMPLATE ARE MISSING. THEY WILL BE INSTALLED BEFORE FURTHER PROCESSING.")
        rc,c_out = run_command(f'install_DSURQE.sh {rabies_path}', verbose=True)


def check_binary_masks(mask):
    img = sitk.ReadImage(mask)
    array = sitk.GetArrayFromImage(img)
    if ((array != 1)*(array != 0)).sum() > 0:
        raise ValueError(
            f"The file {mask} is not a binary mask. Non-binary masks cannot be processed.")


def check_template_overlap(template, mask):
    template_img = sitk.ReadImage(template)
    mask_img = sitk.ReadImage(mask)
    if not template_img.GetOrigin() == mask_img.GetOrigin() and template_img.GetDirection() == mask_img.GetDirection():
        raise ValueError(
            f"The file {mask} does not appear to overlap with provided template {template}.")


def check_resampling_syntax(resampling):
    if resampling == 'inputs_defined':
        return
    try:
        if not 'x' in resampling:
            raise
        shape = resampling.split('x')
        if not len(shape) == 3:
            raise
        spacing = (float(shape[0]), float(shape[1]), float(shape[2]))
    except:
        raise ValueError(
            f"Resampling {resampling} must follow the format 'dim1xdim2xdim3', e.g. '0.1x0.1x0.1', (in mm) following the RAS axis convention.")


def list_files(startpath):
    string = ''
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        string += f'{indent}{os.path.basename(root)}/ \n'
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            string += f'{subindent}{f} \n'
    return string
