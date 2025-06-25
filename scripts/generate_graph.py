import os
import pickle
import SimpleITK as sitk
from nipype import logging, config
from rabies.boilerplate import *
from rabies.parser import get_parser,read_parser
from rabies.run_main import *

if 'XDG_DATA_HOME' in os.environ.keys():
    rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
else:
    rabies_path = os.environ['HOME']+'/.local/share/rabies'


def execute_workflow(args=None):
    # generates the parser CLI and execute the workflow based on specified parameters.
    """
    Prepare and configure a RABIES workflow based on command-line arguments.
    
    Parses CLI arguments, validates output directory paths, initializes logging, verifies template installation, checks for incompatible inclusion/exclusion parameters, and sets up the specified workflow stage (preprocess, confound correction, or analysis). Saves the modified CLI options as a pickle file and returns the prepared workflow object.
    
    Parameters:
        args (list, optional): List of command-line arguments to override sys.argv. If None, uses sys.argv.
    
    Returns:
        workflow: The configured workflow object ready for execution.
    """
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

    return workflow


from rabies.utils import generate_token_data

import tempfile
tmppath = tempfile.mkdtemp()

#### increase the number of scans generated to 3 if running group analysis
generate_token_data(tmppath, number_scans=3)

output_folder = f'{tmppath}/outputs'

#### HERE ARE SET THE DESIRED PARAMETERS FOR PREPROCESSING
args = [
        f'--exclusion_ids',f'{tmppath}/inputs/sub-token1_bold.nii.gz',
        '-f',
        #'--debug',
        'preprocess', f'{tmppath}/inputs', output_folder,
        '--anat_inho_cor', 'method=disable,otsu_thresh=2,multiotsu=false', 
        '--bold_inho_cor', 'method=disable,otsu_thresh=2,multiotsu=false',
        '--bold2anat_coreg', 'registration=no_reg,masking=false,brain_extraction=false', 
        '--commonspace_reg', 'masking=false,brain_extraction=false,fast_commonspace=false,template_registration=SyN', 
        '--data_type', 'int16', 
        '--anat_template', f'{tmppath}/inputs/sub-token1_T1w.nii.gz',
        '--brain_mask', f'{tmppath}/inputs/token_mask.nii.gz', 
        '--WM_mask', f'{tmppath}/inputs/token_mask.nii.gz',
        '--CSF_mask', f'{tmppath}/inputs/token_mask.nii.gz',
        '--vascular_mask', f'{tmppath}/inputs/token_mask.nii.gz', 
        '--labels', f'{tmppath}/inputs/token_mask.nii.gz',
        ]

main_workflow = execute_workflow(args=args)

sub_wf = main_workflow.get_node('commonspace_reg_wf')

sub_wf.write_graph(graph2use='orig', dotfilename='./graph_orig.dot')


