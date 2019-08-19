import os
import sys
from mfp.main_wf import init_anat_init_wf, init_main_postPydpiper_wf

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from pathlib import Path

def get_parser():
    """Build parser object"""
    parser = ArgumentParser()

    parser.add_argument('input_dir', action='store', type=Path,
                        help='the root folder of the input data directory.')
    parser.add_argument('output_dir', action='store', type=Path,
                        help='the output path for the outcomes of preprocessing')
    parser.add_argument("-e", "--bold_only", type=bool, default=0,
                        help="preprocessing with only EPI scans. commonspace registration and distortion correction"
                              " is executed through registration of the EPIs to a common template atlas. Default=False")
    parser.add_argument("-c", "--commonspace_method", type=str, default='pydpiper',
                        help="specify either 'pydpiper' or 'ants_dbm' as common space registration method. Default=pydpiper")
    parser.add_argument("-b", "--bias_reg_script", type=str, default='default',
                        help="specify a registration script for iterative bias field correction. 'default' is a rigid registration. Default=default")
    parser.add_argument("-r", "--coreg_script", type=str, default='SyN',
                        help="Specify EPI to anat coregistration script. Built-in options include 'Rigid', 'Affine' and 'SyN' (non-linear), but"
                        " can specify a custom registration script following the template script structure (see Mouse_fmriPype/mfp/shell_scripts/ for template). Default=SyN")
    parser.add_argument("-p", "--plugin", type=str, default='Linear',
                        help="Specify the nipype plugin for workflow execution. Consult nipype plugin documentation for detailed options."
                             " Linear, MultiProc, SGE and SGEGraph have been tested. Default=Linear")
    parser.add_argument("-d", "--debug", type=bool, default=0,
                        help="Run in debug mode. Default=False")
    parser.add_argument("-v", "--verbose", type=bool, default=0,
                        help="Increase output verbosity. **doesn't do anything for now. Default=False")
    return parser


def execute_workflow():
    opts = get_parser().parse_args()

    bold_preproc_only=opts.bold_only
    bias_reg_script=opts.bias_reg_script
    coreg_script=opts.coreg_script
    commonspace_method=opts.commonspace_method
    data_dir_path=os.path.abspath(str(opts.input_dir))
    output_folder=os.path.abspath(str(opts.output_dir))
    plugin=opts.plugin

    data_csv=data_dir_path+'/data_info.csv'
    csv_labels='/data/chamal/projects/Gabriel_DG/atlases/DSURQE_atlas/labels/DSURQE_40micron_R_mapping.csv' #file with the id# of each label in the atlas to compute WM and CSF masks


    if bold_preproc_only:
        commonspace_transform=False
        compute_WM_CSF_masks=False
    elif not bold_preproc_only:
        anat_init_wf = init_anat_init_wf(data_csv, data_dir_path, output_folder, commonspace_method=commonspace_method)
        anat_init_wf.base_dir = output_folder

        dir_path = os.path.dirname(os.path.realpath(__file__))
        if commonspace_method=='pydpiper':
            model_script_path=dir_path+'/mfp/shell_scripts/pydpiper.sh'
            commonspace_transform=False
        elif commonspace_method=='ants_dbm':
            model_script_path=dir_path+'/mfp/shell_scripts/ants_dbm.sh'
            commonspace_transform=True
        else:
            raise ValueError('Invalid commonspace method.')

        commonspace_csv_file=output_folder+'/anat_init_wf/commonspace_prep/commonspace_input_files.csv'
    else:
        raise ValueError('bold_preproc_only must be true or false.')

    main_postPydpiper_wf = init_main_postPydpiper_wf(data_csv, data_dir_path, output_folder, bold_preproc_only=bold_preproc_only, csv_labels=csv_labels, bias_reg_script=bias_reg_script, coreg_script=coreg_script, commonspace_transform=commonspace_transform)
    main_postPydpiper_wf.base_dir = output_folder

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


    print('Running main workflow with %s plugin.' % plugin)
    if bold_preproc_only:
        main_postPydpiper_wf.run(plugin=plugin, plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})
    else:
        print('Running anat init.')
        anat_init_wf.run(plugin=plugin, plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})

        while(not os.path.isfile(commonspace_csv_file)):
            import time
            print('Anat init not finished, waiting 5min.')
            time.sleep(300)

        print('Running commonspace registration.')
        #run pydpiper
        out_dir=output_folder+'/commonspace/'
        os.system('mkdir -p %s' % (out_dir))
        cwd=os.getcwd()
        os.chdir(out_dir)
        os.system('bash %s %s' % (model_script_path,commonspace_csv_file))
        os.chdir(cwd)

        print('Running main workflow.')
        main_postPydpiper_wf.run(plugin=plugin, plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})
