import os
import pickle
import SimpleITK as sitk
from nipype import logging, config
import pathlib
from .boilerplate import *
from .parser import get_parser,read_parser
from .preprocess_pkg.utils import convert_to_RAS
from . import templates

def execute_workflow(args=None, return_workflow=False):
    # generates the parser CLI and execute the workflow based on specified parameters.
    parser = get_parser()
    opts = read_parser(parser, args)

    if opts.rabies_stage is None: # no processing stage was selected, print the help message instead
        parser.print_help()
        return

    if opts.rabies_stage == 'install': # installs template files, no workflow is executed
        install_template_sets(opts.template_set)
        return

    # convert all input paths to absolute paths
    for arg in vars(opts):
        attr = getattr(opts, arg)
        if isinstance(attr,pathlib.Path):
            # convert to string since some downstream code needs cannot process pathlib.Path as input
            setattr(opts, arg, str(attr.resolve()))

    if not os.path.isdir(opts.output_dir):
        os.makedirs(opts.output_dir)

    log = prep_logging(opts, opts.output_dir)


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
    
    ###ITK THREAD MANAGEMENT###
    if opts.num_ITK_threads=='optimal':
        num_ITK_threads=1 # until the number of scan is determined, this is set to 1 by default
        os.environ['RABIES_ITK_THREADS_STATE']='ON'
    elif opts.num_ITK_threads=='off':
        num_ITK_threads=1 # set to 1 to maximize nipype parallelization
        os.environ['RABIES_ITK_THREADS_STATE']='OFF' # this means run_command() won't set constraints on the number of threads
    else:
        num_ITK_threads=int(opts.num_ITK_threads)
        os.environ['RABIES_ITK_THREADS_STATE']='ON'
    # this is set as an os.environ variable instead of a nipype node input, to avoid re-running nodes when changing this parameter
    os.environ['RABIES_ITK_NUM_THREADS']=str(num_ITK_threads)

    if not opts.num_ITK_threads=='off':
        # python parallelization is entirely handled outside of SITK, so individual commands should inherit 1 CPU 
        sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)

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

    if str(opts.interpolation) == 'Linear':
        opts.interpolation_sitk = sitk.sitkLinear
    elif str(opts.interpolation) == 'BSpline3':
        opts.interpolation_sitk = sitk.sitkBSpline3
    elif str(opts.interpolation) == 'LanczosSinc':
        opts.interpolation_sitk = sitk.sitkLanczosWindowedSinc
    elif str(opts.interpolation) == 'CosineSinc':
        opts.interpolation_sitk = sitk.sitkCosineWindowedSinc
    else:
        raise ValueError('Invalid input for --interpolation.')

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


def image_extent(file):
    # largest physical dimension of an image, in mm
    img = sitk.ReadImage(file)
    return max([spacing*size for spacing,size
                in zip(img.GetSpacing()[:3], img.GetSize()[:3])])


def find_reference_scan(bids_dir, bold_only):
    # first input image found, used to compare the size of the data against the template
    modality = 'func' if bold_only else 'anat'
    scans = sorted(pathlib.Path(bids_dir).glob(f'sub-*/**/{modality}/*.nii*'))
    if len(scans)==0: # datasets without sessions, or without the expected modality folder
        scans = sorted(pathlib.Path(bids_dir).glob('sub-*/**/*.nii*'))
    return str(scans[0]) if len(scans)>0 else None


def check_template_scale(opts, log):
    # registering data to the template of another species does not fail, it produces
    # meaningless outputs, so the mismatch is caught before the workflow is built
    if opts.skip_scale_check:
        return

    scan = find_reference_scan(opts.bids_dir, opts.bold_only)
    if scan is None:
        log.warning("No input image was found to compare against the template size; "
                    "skipping the size check.")
        return

    extent = image_extent(scan)
    template_extents = {}
    for name in templates.TEMPLATE_SET_NAMES:
        template = templates.resolve(name, 'anat_template', bold_only=opts.bold_only)
        if os.path.isfile(template): # a set that is not installed cannot be compared
            template_extents[name] = image_extent(template)

    better = templates.scale_verdict(extent, template_extents, opts.template_set)
    if better is None:
        return

    raise ValueError(
        f"The input image {scan} is "
        f"{extent/template_extents[opts.template_set]:.1f} times the size of the "
        f"commonspace template of --template_set {opts.template_set}, but "
        f"{extent/template_extents[better]:.1f} times the size of the {better} one. "
        "This usually means the data comes from another species; consider "
        f"--template_set {better}. If the data and the selected template do match, "
        "re-run with --skip_scale_check.")


def check_inherited_template_set(opts):
    # --inherit_unbiased_template overrides the template files with those of a previous
    # run, so an explicitly selected set that disagrees with that run would be silently
    # discarded and the data registered to the wrong space
    cli_file = f'{opts.inherit_unbiased_template}/rabies_preprocess.pkl'
    if not os.path.isfile(cli_file):
        raise ValueError(f"--inherit_unbiased_template path {opts.inherit_unbiased_template} "
                         "does not contain a rabies_preprocess.pkl file.")
    with open(cli_file, 'rb') as handle:
        inherited_opts = pickle.load(handle)
    inherited_set = get_template_set(inherited_opts)

    if opts.explicit_template_set and not opts.template_set==inherited_set:
        raise ValueError(
            f"--template_set {opts.template_set} was selected, but "
            f"--inherit_unbiased_template inherits the template files of a run that used "
            f"the {inherited_set} set, and those files take precedence. Run with "
            f"--template_set {inherited_set}, or without --inherit_unbiased_template.")
    opts.template_set = inherited_set


def preprocess(opts, log):

    if not os.path.isdir(opts.bids_dir):
        raise ValueError("The provided BIDS data path doesn't exist.")
    else:
        # print the input data directory tree
        log.info("INPUT BIDS DATASET:  \n" + list_files(str(opts.bids_dir)))
    
    if not opts.inherit_unbiased_template=='none':
        opts.inherit_unbiased_template = os.path.abspath(opts.inherit_unbiased_template)
        if not os.path.isdir(opts.inherit_unbiased_template):
            raise ValueError(f"--inherit_unbiased_template path {opts.inherit_unbiased_template} doesn't exist.")
        # resolved before the template files, since the inherited run fixes which set is used
        check_inherited_template_set(opts)

    require_template_set(opts.template_set, log)
    templates.resolve_options(opts, log)
    check_template_scale(opts, log)

    # final check of template file formats
    for opt_key,check_binary in zip(['anat_template', 'brain_mask', 'WM_mask','CSF_mask','vascular_mask'],
                                    [False,True,True,True,True]):
        opt_file = getattr(opts, opt_key)
        if opt_file is not None: # some masks might be set to None
            if not os.path.isfile(opt_file):
                raise ValueError(f"--{opt_key} file {opt_file} doesn't exist.")
            opt_file = convert_to_RAS(
                str(opt_file), opts.output_dir+'/template_files')

            if check_binary:
                check_binary_masks(opt_file)
            check_template_overlap(opts.anat_template, opt_file)
            setattr(opts, opt_key, opt_file)

    check_resampling_syntax(opts.nativespace_resampling)
    check_resampling_syntax(opts.commonspace_resampling)
    check_resampling_syntax(opts.anatomical_resampling)

    # prepare nativespace/commonspace resampling arguments
    if opts.resampling_space=='both':
        opts.generate_commonspace = True
        opts.generate_nativespace = True
    elif opts.resampling_space=='common_only':
        opts.generate_commonspace = True
        opts.generate_nativespace = False
    elif opts.resampling_space=='native_only':
        opts.generate_commonspace = False
        opts.generate_nativespace = True

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

    if opts.conf_list:
        raise ValueError("--conf_list parameter is deprecated, and was replaced by --nuisance_regressors.")

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


def get_template_set(preprocess_opts):
    # runs preprocessed before --template_set existed carry no such attribute, and
    # were necessarily run with the mouse files that were the only ones available
    return getattr(preprocess_opts, 'template_set', 'mouse')


def analysis(opts, log):

    cli_file = f'{opts.confound_correction_out}/rabies_confound_correction.pkl'
    with open(cli_file, 'rb') as handle:
        confound_correction_opts = pickle.load(handle)

    cli_file = f'{confound_correction_opts.preprocess_out}/rabies_preprocess.pkl'
    with open(cli_file, 'rb') as handle:
        preprocess_opts = pickle.load(handle)

    labels_file = opts.ROI_labels_file
    # the template set is fixed at preprocessing; analysis files default to that same set
    # so that they are guaranteed to live in the commonspace the data was registered to
    template_set = get_template_set(preprocess_opts)
    bold_only = preprocess_opts.bold_only
    require_template_set(template_set, log)

    if labels_file is None:
        if getattr(preprocess_opts, 'custom_anat_template', False):
            # the set's atlas is not aligned with a template provided by the user
            opts.ROI_labels_file = None
        else:
            # files distributed with a template set are already aligned with its template,
            # so they are used as-is rather than converted and checked against it
            opts.ROI_labels_file = templates.resolve(template_set, 'labels', bold_only=bold_only)
        if opts.ROI_labels_file is None:
            # left as None so that no computation is attempted using the labels
            log.info("No labels file is available for the template used during preprocessing; "
                     "operations depending on --ROI_labels_file are disabled.")
        elif bold_only:
            log.info('With --bold_only, default --ROI_labels_file changed to '+opts.ROI_labels_file)
    else:
        if not os.path.isfile(labels_file):
            raise ValueError(f"--ROI_labels_file file {labels_file} doesn't exist.")
        # need to convert to RAS first, since the template file was also converted
        labels_file = convert_to_RAS(
            str(labels_file), preprocess_opts.output_dir+'/template_files')
        check_template_overlap(preprocess_opts.anat_template, labels_file)
        opts.ROI_labels_file = labels_file

    if opts.prior_maps is None:
        # the default prior maps are those of the template set used during preprocessing;
        # a set that provides none leaves them unset, and the analysis workflow raises if
        # an analysis requiring them was selected
        opts.prior_maps = templates.resolve(template_set, 'prior_maps', bold_only=bold_only)
        if opts.prior_maps is None:
            log.info(f"The {template_set} template set provides no ICA prior maps; "
                     "operations depending on --prior_maps are disabled.")
        elif bold_only:
            log.info('With --bold_only, default --prior_maps changed to '+opts.prior_maps)
    else:
        if not os.path.isfile(str(opts.prior_maps)):
            raise ValueError(f"--prior_maps file {opts.prior_maps} doesn't exist.")
        opts.prior_maps = os.path.abspath(str(opts.prior_maps))

    seed_file_list = []
    seed_name_list = []
    prebuilt_seeds = templates.seed_names(template_set)
    for seed in opts.seed_list:
        if seed in prebuilt_seeds:
            seed_file = templates.seed_file(template_set, seed, bold_only=bold_only)
            seed_name = seed
        else:
            seed_file = pathlib.Path(seed).resolve() # convert to absolute path
            seed_name = seed_file.name.rsplit(".nii")[0]
        if not os.path.isfile(seed_file):
            raise ValueError(
                f"Provide seed file path {seed_file} doesn't exist.")

        check_template_overlap(preprocess_opts.anat_template, seed_file)

        seed_file_list.append(seed_file)
        seed_name_list.append(seed_name)
    opts.seed_list = seed_file_list # override old input so that pre-built seeds are now files
    opts.seed_name_list = seed_name_list # store for developing iterables later

    from rabies.analysis_pkg.main_wf import init_main_analysis_wf
    workflow = init_main_analysis_wf(confound_correction_opts, opts)

    return workflow

def install_template_set(template_set):
    # downloads the files of a template set, unless they are already installed
    if len(templates.missing_files(template_set))==0:
        return False

    from rabies.utils import run_command
    script = templates.TEMPLATE_SETS[template_set]['install_script']
    rc,c_out = run_command(f'{script} {templates.rabies_path}', verbose=True)
    return True


def install_template_sets(template_set):
    # handles the `rabies install` stage, which runs without an output folder or workflow
    set_list = templates.TEMPLATE_SET_NAMES if template_set=='all' else [template_set]
    for name in set_list:
        print(f"Checking the {name} template set.")
        if install_template_set(name):
            print(f"Installed the {name} template set under {templates.rabies_path}.")
        else:
            print(f"The {name} template set is already installed.")
        missing = templates.missing_files(name)
        if len(missing)>0:
            raise ValueError(
                f"The {name} template set is still incomplete after installation. "
                f"Missing files: {missing}")


def require_template_set(template_set, log):
    # verifies that a template set is available, installing it if it is meant to be
    # installed on demand, and failing before the workflow is built otherwise
    missing = templates.missing_files(template_set)
    if len(missing)==0:
        return

    if templates.TEMPLATE_SETS[template_set]['auto_install']:
        log.info(
            f"SOME FILES FROM THE {template_set} TEMPLATE SET ARE MISSING. "
            "THEY WILL BE INSTALLED BEFORE FURTHER PROCESSING.")
        install_template_set(template_set)
        missing = templates.missing_files(template_set)
        if len(missing)==0:
            return

    raise ValueError(
        f"The {template_set} template set is not installed. Run `rabies install {template_set}` "
        "to download it, which must be done from a machine with network access before "
        f"running the pipeline. Missing files: {missing}")


def check_binary_masks(mask):
    img = sitk.ReadImage(mask)
    array = sitk.GetArrayFromImage(img)
    if ((array != 1)*(array != 0)).sum() > 0:
        raise ValueError(
            f"The file {mask} is not a binary mask. Non-binary masks cannot be processed.")


def check_template_overlap(template, mask):
    template_img = sitk.ReadImage(template)
    mask_img = sitk.ReadImage(mask)
    if not (template_img.GetOrigin() == mask_img.GetOrigin() and template_img.GetDirection() == mask_img.GetDirection()):
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
