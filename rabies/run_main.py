import os
import sys
from rabies.main_wf import init_anat_init_wf, init_main_postcommonspace_wf

import argparse
from pathlib import Path
import pathos.multiprocessing as multiprocessing  # Better multiprocessing

def get_parser():
    """Build parser object"""
    parser = argparse.ArgumentParser(
        description="""RABIES performs preprocessing of rodent fMRI images. Can either run
        on datasets that only contain EPI images, or both structural and EPI images. Refer
        to the README documentation for the input folder structure.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_dir', action='store', type=Path,
                        help='the root folder of the input data directory.')
    parser.add_argument('output_dir', action='store', type=Path,
                        help='the output path to drop outputs from major preprocessing steps.')
    parser.add_argument("-e", "--bold_only", type=bool, default=False,
                        help="preprocessing with only EPI scans. commonspace registration and distortion correction"
                              " is executed through registration of the EPIs to a common template atlas.")
    parser.add_argument("-c", "--commonspace_method", type=str, default='ants_dbm',
                        help="specify either 'pydpiper' or 'ants_dbm' as common space registration method. Pydpiper can only be "
                        "executed in parallel with SGE or PBS. ***pydpiper option in development")
    parser.add_argument("-b", "--bias_reg_script", type=str, default='Rigid',
                        help="specify a registration script for iterative bias field correction. 'default' is a rigid registration.")
    parser.add_argument("-r", "--coreg_script", type=str, default='SyN',
                        help="Specify EPI to anat coregistration script. Built-in options include 'Rigid', 'Affine' and 'SyN' (non-linear), but"
                        " can specify a custom registration script following the template script structure (see RABIES/rabies/shell_scripts/ for template).")
    parser.add_argument("-p", "--plugin", type=str, default='Linear',
                        help="Specify the nipype plugin for workflow execution. Consult nipype plugin documentation for detailed options."
                             " Linear, MultiProc, SGE and SGEGraph have been tested.")
    parser.add_argument("-d", "--debug", type=bool, default=False,
                        help="Run in debug mode. Default=False")
    parser.add_argument("-v", "--verbose", type=bool, default=False,
                        help="Increase output verbosity. **doesn't do anything for now.")

    g_ants_dbm = parser.add_argument_group('cluster options if commonspace method is ants_dbm (taken from twolevel_dbm.py):')
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
    g_ants_dbm.add_argument(
        '--local_threads',
        '-j',
        type=int,
        default=multiprocessing.cpu_count(),
        help="""For local execution, how many subject-wise modelbuilds to run in parallel,
        defaults to number of CPUs""")


    g_stc = parser.add_argument_group('Specify Slice Timing Correction info that is fed to AFNI 3dTshift.')
    g_stc.add_argument('--STC', type=bool, default=True,
                        help="Whether to run STC or not.")
    g_stc.add_argument('--TR', type=str, default='1.0s',
                        help="Specify repetition time (TR).")
    g_stc.add_argument('--tpattern', type=str, default='alt',
                        help="Specify if interleaved or sequential acquisition. 'alt' for interleaved, 'seq' for sequential.")

    g_template = parser.add_argument_group('Template files.')
    g_template.add_argument('--anat_template', action='store', type=Path,
                        default="%s/DSURQE_atlas/nifti/DSURQE_100micron_average.nii.gz" % (os.environ["RABIES"]),
                        help='Anatomical file for the commonspace template.')
    g_template.add_argument('--brain_mask', action='store', type=Path,
                        default="%s/DSURQE_atlas/nifti/DSURQE_100micron_mask.nii.gz" % (os.environ["RABIES"]),
                        help='Brain mask for the template.')
    g_template.add_argument('--WM_mask', action='store', type=Path,
                        default="%s/DSURQE_atlas/nifti/DSURQE_100micron_eroded_WM_mask.nii.gz" % (os.environ["RABIES"]),
                        help='White matter mask for the template.')
    g_template.add_argument('--CSF_mask', action='store', type=Path,
                        default="%s/DSURQE_atlas/nifti/DSURQE_100micron_eroded_CSF_mask.nii.gz" % (os.environ["RABIES"]),
                        help='CSF mask for the template.')
    g_template.add_argument('--labels', action='store', type=Path,
                        default="%s/DSURQE_atlas/nifti/DSURQE_100micron_labels.nii.gz" % (os.environ["RABIES"]),
                        help='Atlas file with anatomical labels.')
    g_template.add_argument('--csv_labels', action='store', type=Path,
                        default="%s/DSURQE_atlas/DSURQE_40micron_R_mapping.csv" % (os.environ["RABIES"]),
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
    stc_bool=opts.STC
    stc_TR=opts.TR
    if opts.tpattern=="alt":
        stc_tpattern='altplus'
    elif opts.tpattern=="seq":
        stc_tpattern='seqplus'
    else:
        raise ValueError('Invalid --tpattern provided.')

    #setting absolute paths for ants_dbm options options
    os.environ["ants_dbm_cluster_type"]=opts.cluster_type
    os.environ["ants_dbm_walltime"]=opts.walltime
    os.environ["ants_dbm_memory_request"]=opts.memory_request
    os.environ["ants_dbm_local_threads"]=str(opts.local_threads)

    #template options
    # set OS paths to template and atlas files
    os.environ["template_anat"] = str(opts.anat_template)
    if not os.path.isfile(os.environ["template_anat"]):
        raise ValueError("--anat_template file doesn't exists.")

    os.environ["template_mask"] = str(opts.brain_mask)
    if not os.path.isfile(os.environ["template_mask"]):
        raise ValueError("--brain_mask file doesn't exists.")

    os.environ["WM_mask"] = str(opts.WM_mask)
    if not os.path.isfile(os.environ["WM_mask"]):
        raise ValueError("--WM_mask file doesn't exists.")

    os.environ["CSF_mask"] = str(opts.CSF_mask)
    if not os.path.isfile(os.environ["CSF_mask"]):
        raise ValueError("--CSF_mask file doesn't exists.")

    os.environ["atlas_labels"] = str(opts.labels)
    if not os.path.isfile(os.environ["atlas_labels"]):
        raise ValueError("--labels file doesn't exists.")

    os.environ["csv_labels"] = str(opts.csv_labels)
    if not os.path.isfile(os.environ["csv_labels"]):
        raise ValueError("--csv_labels file doesn't exists.")

    os.environ["template_anat_mnc"] = "%s/DSURQE_atlas/minc/DSURQE_100micron_average.mnc" % (os.environ["RABIES"])
    os.environ["template_mask_mnc"] = "%s/DSURQE_atlas/minc/DSURQE_100micron_mask.mnc" % (os.environ["RABIES"])

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
            model_script_path=dir_path+'/shell_scripts/pydpiper.sh'
            commonspace_transform=False
        elif commonspace_method=='ants_dbm':
            model_script_path=dir_path+'/shell_scripts/ants_dbm.sh'
            commonspace_transform=True
        else:
            raise ValueError('Invalid commonspace method.')

        commonspace_csv_file=output_folder+'/anat_init_wf/commonspace_prep/commonspace_input_files.csv'
        commonspace_info_csv=output_folder+'/anat_init_wf/commonspace_prep/commonspace_info.csv'
    else:
        raise ValueError('bold_preproc_only must be true or false.')

    main_postcommonspace_wf = init_main_postcommonspace_wf(data_csv, data_dir_path, output_folder, apply_STC=stc_bool, tr=stc_TR, tpattern=stc_tpattern, bold_preproc_only=bold_preproc_only, csv_labels=csv_labels, bias_reg_script=bias_reg_script, coreg_script=coreg_script, commonspace_transform=commonspace_transform)
    main_postcommonspace_wf.base_dir = output_folder

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
        main_postcommonspace_wf.run(plugin=plugin, plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})
    else:
        print('Running anat init.')
        anat_init_wf.run(plugin=plugin, plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})

        while(not os.path.isfile(commonspace_csv_file)):
            import time
            print('Anat init not finished, waiting 5min.')
            time.sleep(300)

        print('Running commonspace registration.')
        #run commonspace
        out_dir=output_folder+'/commonspace/'
        os.system('mkdir -p %s' % (out_dir))
        cwd=os.getcwd()
        os.chdir(out_dir)
        os.system('bash %s %s %s' % (model_script_path,commonspace_csv_file,commonspace_info_csv))
        os.chdir(cwd)

        print('Running main workflow.')
        main_postcommonspace_wf.run(plugin=plugin, plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})
