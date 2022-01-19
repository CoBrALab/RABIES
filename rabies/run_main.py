import os
import pickle
import SimpleITK as sitk
from nipype import logging, config
from .boilerplate import *
from .parser import get_parser

if 'XDG_DATA_HOME' in os.environ.keys():
    rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
else:
    rabies_path = os.environ['HOME']+'/.local/share/rabies'


def execute_workflow():
    # generates the parser CLI and execute the workflow based on specified parameters.
    parser = get_parser()
    opts = parser.parse_args()

    try:
        output_folder = os.path.abspath(str(opts.output_dir))
    except:
        parser.print_help()
        return

    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    log = prep_logging(opts, output_folder)

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

    if opts.rabies_stage == 'preprocess':
        workflow = preprocess(opts, log)
    elif opts.rabies_stage == 'confound_correction':
        workflow = confound_correction(opts, log)
    elif opts.rabies_stage == 'analysis':
        workflow = analysis(opts, log)
    else:
        parser.print_help()
    workflow.base_dir = output_folder


    try:
        log.info(f'Running workflow with {opts.plugin} plugin.')
        # execute workflow, with plugin_args limiting the cluster load for parallel execution
        graph_out = workflow.run(plugin=opts.plugin, plugin_args={'max_jobs': 50, 'dont_resubmit_completed_jobs': True,
                                                      'n_procs': opts.local_threads, 'qsub_args': f'-pe smp {str(opts.min_proc)}'})
        # save the workflow execution
        workflow_file = f'{output_folder}/rabies_{opts.rabies_stage}_workflow.pkl'
        with open(workflow_file, 'wb') as handle:
            pickle.dump(graph_out, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    except Exception as e:
        log.critical(f'RABIES failed: {e}')
        raise


def prep_logging(opts, output_folder):
    cli_file = f'{output_folder}/rabies_{opts.rabies_stage}.pkl'
    if os.path.isfile(cli_file):
        raise ValueError(f"""
            A previous run was indicated by the presence of {cli_file}.
            This can lead to inconsistencies between previous outputs and the log files.
            To prevent this, you are required to manually remove {cli_file}, and we 
            recommend also removing previous datasinks from the {opts.rabies_stage} RABIES step.
            """)

    with open(cli_file, 'wb') as handle:
        pickle.dump(opts, handle, protocol=pickle.HIGHEST_PROTOCOL)

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
    elif opts.verbose<=2:
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
    # Verify input and output directories
    data_dir_path = os.path.abspath(str(opts.bids_dir))
    output_folder = os.path.abspath(str(opts.output_dir))
    if not os.path.isdir(data_dir_path):
        raise ValueError("The provided BIDS data path doesn't exists.")
    else:
        # print the input data directory tree
        log.info("INPUT BIDS DATASET:  \n" + list_files(data_dir_path))

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
    # if --bold_only, the default atlas files change to EPI versions
    if opts.bold_only:
        if str(opts.anat_template)==f"{rabies_path}/DSURQE_40micron_average.nii.gz":
            file=f"{rabies_path}/EPI_template.nii.gz"
            opts.anat_template=file
            log.info('With --bold_only, default --anat_template changed to '+file)
        if str(opts.brain_mask)==f"{rabies_path}/DSURQE_40micron_mask.nii.gz":
            file=f"{rabies_path}/EPI_brain_mask.nii.gz"
            opts.brain_mask=file
            log.info('With --bold_only, default --brain_mask changed to '+file)
        if str(opts.WM_mask)==f"{rabies_path}/DSURQE_40micron_eroded_WM_mask.nii.gz":
            file=f"{rabies_path}/EPI_WM_mask.nii.gz"
            opts.WM_mask=file
            log.info('With --bold_only, default --WM_mask changed to '+file)
        if str(opts.CSF_mask)==f"{rabies_path}/DSURQE_40micron_eroded_CSF_mask.nii.gz":
            file=f"{rabies_path}/EPI_CSF_mask.nii.gz"
            opts.CSF_mask=file
            log.info('With --bold_only, default --CSF_mask changed to '+file)
        if str(opts.vascular_mask)==f"{rabies_path}/vascular_mask.nii.gz":
            file=f"{rabies_path}/EPI_vascular_mask.nii.gz"
            opts.vascular_mask=file
            log.info('With --bold_only, default --vascular_mask changed to '+file)
        if str(opts.labels)==f"{rabies_path}/DSURQE_40micron_labels.nii.gz":
            file=f"{rabies_path}/EPI_labels.nii.gz"
            opts.labels=file
            log.info('With --bold_only, default --labels changed to '+file)

    # make sure we have absolute paths
    opts.anat_template = os.path.abspath(opts.anat_template)
    opts.brain_mask = os.path.abspath(opts.brain_mask)
    opts.WM_mask = os.path.abspath(opts.WM_mask)
    opts.CSF_mask = os.path.abspath(opts.CSF_mask)
    opts.vascular_mask = os.path.abspath(opts.vascular_mask)
    opts.labels = os.path.abspath(opts.labels)

    # To use brain extraction, make sure either --commonspace_masking or --coreg_masking is used
    if opts.brain_extraction:
        if not (opts.commonspace_masking or opts.coreg_masking):
            raise ValueError(f"To use --brain_extraction, you must select either --commonspace_masking or --coreg_masking.")

    # convert template files to RAS convention if they aren't already
    from rabies.preprocess_pkg.utils import convert_to_RAS
    if not os.path.isfile(opts.anat_template):
        raise ValueError(f"--anat_template file {opts.anat_template} doesn't exists.")
    opts.anat_template = convert_to_RAS(
        str(opts.anat_template), output_folder+'/template_files')

    if not os.path.isfile(opts.brain_mask):
        raise ValueError(f"--brain_mask file {opts.brain_mask} doesn't exists.")
    check_binary_masks(opts.brain_mask)
    opts.brain_mask = convert_to_RAS(
        str(opts.brain_mask), output_folder+'/template_files')
    check_template_overlap(opts.anat_template, opts.brain_mask)

    if not os.path.isfile(opts.WM_mask):
        raise ValueError(f"--WM_mask file {opts.WM_mask} doesn't exists.")
    check_binary_masks(opts.WM_mask)
    opts.WM_mask = convert_to_RAS(
        str(opts.WM_mask), output_folder+'/template_files')
    check_template_overlap(opts.anat_template, opts.WM_mask)

    if not os.path.isfile(opts.CSF_mask):
        raise ValueError(f"--CSF_mask file {opts.CSF_mask} doesn't exists.")
    check_binary_masks(opts.CSF_mask)
    opts.CSF_mask = convert_to_RAS(
        str(opts.CSF_mask), output_folder+'/template_files')
    check_template_overlap(opts.anat_template, opts.CSF_mask)

    if not os.path.isfile(opts.vascular_mask):
        raise ValueError(f"--vascular_mask file {opts.vascular_mask} doesn't exists.")
    check_binary_masks(opts.vascular_mask)
    opts.vascular_mask = convert_to_RAS(
        str(opts.vascular_mask), output_folder+'/template_files')
    check_template_overlap(opts.anat_template, opts.vascular_mask)

    if not os.path.isfile(opts.labels):
        raise ValueError(f"--labels file {opts.labels} doesn't exists.")
    opts.labels = convert_to_RAS(
        str(opts.labels), output_folder+'/template_files')
    check_template_overlap(opts.anat_template, opts.labels)

    check_resampling_syntax(opts.nativespace_resampling)
    check_resampling_syntax(opts.commonspace_resampling)
    check_resampling_syntax(opts.anatomical_resampling)

    # write boilerplate
    boilerplate_file = f'{output_folder}/boilerplate.txt'

    methods,ref_string = preprocess_boilerplate(opts)
    txt_boilerplate="#######PREPROCESSING\n\n"+methods+ref_string+'\n\n'
    with open(boilerplate_file, "w") as text_file:
        text_file.write(txt_boilerplate)

    from rabies.preprocess_pkg.main_wf import init_main_wf
    workflow = init_main_wf(data_dir_path, output_folder, opts)

    return workflow


def confound_correction(opts, log):

    if opts.edge_cutoff == 0 and ((opts.highpass is not None) or (opts.lowpass is not None)):
        log.warning(
            "\n############################################# WARNING\n"
            "Highpass and/or lowpass filtering will be applied without removing timepoints "
            "at the edge of acquisition. This may introduce edge artefacts not accounted for, and "
            "not recommended. We recommend removing ~30sec at both end of the acquisition."
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
        if str(opts.prior_maps)==f"{rabies_path}/melodic_IC.nii.gz":
            file=f"{rabies_path}/melodic_IC_resampled.nii.gz"
            opts.prior_maps=file
            log.info('With --bold_only, default --prior_maps changed to '+file)


    from rabies.analysis_pkg.main_wf import init_main_analysis_wf
    workflow = init_main_analysis_wf(preprocess_opts, confound_correction_opts, opts)

    return workflow

def install_DSURQE(log):

    install = False
    # verifies whether default template files are installed and installs them otherwise
    if not os.path.isfile(f"{rabies_path}/DSURQE_40micron_average.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/DSURQE_40micron_mask.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/DSURQE_40micron_eroded_WM_mask.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/DSURQE_40micron_eroded_CSF_mask.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/DSURQE_40micron_labels.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/vascular_mask.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/melodic_IC.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/EPI_template.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/EPI_brain_mask.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/EPI_WM_mask.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/EPI_CSF_mask.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/EPI_vascular_mask.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/EPI_labels.nii.gz"):
        install = True
    elif not os.path.isfile(f"{rabies_path}/melodic_IC_resampled.nii.gz"):
        install = True
    if install:
        from rabies.preprocess_pkg.utils import run_command
        log.info(
            "SOME FILES FROM THE DEFAULT TEMPLATE ARE MISSING. THEY WILL BE INSTALLED BEFORE FURTHER PROCESSING.")
        rc = run_command(f'install_DSURQE.sh {rabies_path}', verbose=True)


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
