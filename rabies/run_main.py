import os
import sys
from rabies.main_wf import init_anat_init_wf, init_main_postPydpiper_wf

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
    parser.add_argument("-b", "--bias_reg_script", type=str, default='Rigid',
                        help="specify a registration script for iterative bias field correction. 'default' is a rigid registration. Default=default")
    parser.add_argument("-r", "--coreg_script", type=str, default='SyN',
                        help="Specify EPI to anat coregistration script. Built-in options include 'Rigid', 'Affine' and 'SyN' (non-linear), but"
                        " can specify a custom registration script following the template script structure (see RABIES/rabies/shell_scripts/ for template). Default=SyN")
    parser.add_argument("-p", "--plugin", type=str, default='Linear',
                        help="Specify the nipype plugin for workflow execution. Consult nipype plugin documentation for detailed options."
                             " Linear, MultiProc, SGE and SGEGraph have been tested. Default=Linear")
    parser.add_argument("-d", "--debug", type=bool, default=0,
                        help="Run in debug mode. Default=False")
    parser.add_argument("-v", "--verbose", type=bool, default=0,
                        help="Increase output verbosity. **doesn't do anything for now. Default=False")

    g_stc = parser.add_argument_group('Specify Slice Timing Correction info that is fed to AFNI 3dTshift.')
    g_stc.add_argument('--STC', action='store_true', type=bool, default=True,
                        help="Whether to run STC or not. Default=True")
    g_stc.add_argument('--TR', action='store_true',
                        type=str, default='1.0s',
                        help="Anatomical file for the commonspace template. Default=1.0s")
    g_stc.add_argument('--tpattern', action='store_true',
                        type=str, default='alt',
                        help="Specify if interleaved or sequential acquisition. 'alt' for interleaved, 'seq' for sequential. Default=alt")

    g_template = parser.add_argument_group('Template files. ***under development, can only use the DSURQE atlas as default for now.')
    g_template.add_argument('--anat_template', action='store_true',
                        default=False,
                        help='Anatomical file for the commonspace template.')
    g_template.add_argument('--brain_mask', action='store_true',
                        default=False,
                        help='Brain mask for the template.')
    g_template.add_argument('--WM_mask', action='store_true',
                        default=False,
                        help='White matter mask for the template.')
    g_template.add_argument('--CSF_mask', action='store_true',
                        default=False,
                        help='CSF mask for the template.')
    g_template.add_argument('--labels', action='store_true',
                        default=False,
                        help='Atlas file with anatomical labels.')
    g_template.add_argument('--csv_labels', action='store_true',
                        default=False,
                        help='csv file with info on the labels.')
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

    #STC options
    stc_bool=opts.g_stc.STC
    stc_TR=opts.g_stc.TR
    if opts.g_stc.tpattern=="alt":
        stc_tpattern='altplus'
    elif opts.g_stc.tpattern=="seq":
        stc_tpattern='seqplus'
    else:
        raise ValueError('Invalid --tpattern provided.')

    #template options
    # set OS paths to template and atlas files
    os.environ["template_anat"] = "%s/DSURQE_atlas/nifti/DSURQE_100micron_average.nii.gz" % (os.environ["RABIES"])
    os.environ["template_mask"] = "%s/DSURQE_atlas/nifti/DSURQE_100micron_mask.nii.gz" % (os.environ["RABIES"])
    os.environ["template_anat_mnc"] = "%s/DSURQE_atlas/minc/DSURQE_100micron_average.mnc" % (os.environ["RABIES"])
    os.environ["template_mask_mnc"] = "%s/DSURQE_atlas/minc/DSURQE_100micron_mask.mnc" % (os.environ["RABIES"])
    os.environ["atlas_labels"] = "%s/DSURQE_atlas/nifti/DSURQE_100micron_labels.nii.gz" % (os.environ["RABIES"])
    os.environ["csv_labels"] = "%s/DSURQE_atlas/DSURQE_40micron_R_mapping.csv" % (os.environ["RABIES"])
    os.environ["WM_mask"] = "%s/DSURQE_atlas/nifti/DSURQE_100micron_eroded_WM_mask.nii.gz" % (os.environ["RABIES"])
    os.environ["CSF_mask"] = "%s/DSURQE_atlas/nifti/DSURQE_100micron_eroded_CSF_mask.nii.gz" % (os.environ["RABIES"])

    data_csv=data_dir_path+'/data_info.csv'
    csv_labels=os.environ["csv_labels"] #file with the id# of each label in the atlas to compute WM and CSF masks

    if bold_preproc_only:
        commonspace_transform=False
        compute_WM_CSF_masks=False
    elif not bold_preproc_only:
        anat_init_wf = init_anat_init_wf(data_csv, data_dir_path, output_folder, commonspace_method=commonspace_method)
        anat_init_wf.base_dir = output_folder

        dir_path = os.path.dirname(os.path.realpath(__file__))
        if commonspace_method=='pydpiper':
            model_script_path=dir_path+'/rabies/shell_scripts/pydpiper.sh'
            commonspace_transform=False
        elif commonspace_method=='ants_dbm':
            model_script_path=dir_path+'/rabies/shell_scripts/ants_dbm.sh'
            commonspace_transform=True
        else:
            raise ValueError('Invalid commonspace method.')

        commonspace_csv_file=output_folder+'/anat_init_wf/commonspace_prep/commonspace_input_files.csv'
    else:
        raise ValueError('bold_preproc_only must be true or false.')

    main_postPydpiper_wf = init_main_postPydpiper_wf(data_csv, data_dir_path, output_folder, apply_STC=stc_bool tr=tr, tpattern=tpattern, bold_preproc_only=bold_preproc_only, csv_labels=csv_labels, bias_reg_script=bias_reg_script, coreg_script=coreg_script, commonspace_transform=commonspace_transform)
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
