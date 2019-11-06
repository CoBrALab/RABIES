import os
import sys

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
                        help="Specify EPI to anat coregistration script. Built-in options include 'Rigid', 'Affine', 'SyN' (non-linear) and 'light_SyN', but"
                        " can specify a custom registration script following the template script structure (see RABIES/rabies/shell_scripts/ for template).")
    parser.add_argument("-p", "--plugin", type=str, default='Linear',
                        help="Specify the nipype plugin for workflow execution. Consult nipype plugin documentation for detailed options."
                             " Linear, MultiProc, SGE and SGEGraph have been tested.")
    parser.add_argument("--min_proc", type=int, default=1,
                        help="For parallel processing, specify the minimal number of nodes to be assigned.")
    parser.add_argument("-d", "--debug", type=bool, default=False,
                        help="Run in debug mode. Default=False")
    parser.add_argument("-v", "--verbose", type=bool, default=False,
                        help="Increase output verbosity. **doesn't do anything for now.")

    g_resampling = parser.add_argument_group('Options for the resampling of the EPI for motion realignment, susceptibility distortion correction and common space resampling:')
    g_resampling.add_argument('--isotropic_resampling', type=bool, default=False,
                        help="Whether to resample the EPI to an isotropic resolution based on the lowest dimension.")
    g_resampling.add_argument('--upsampling', type=float, default=1.0,
                        help="Can specify a upsampling parameter to increase the EPI resolution upon resampling and minimize information lost from the transform interpolation.")
    g_resampling.add_argument("--data_type", type=str, default='float64',
                        help="Specify resampling data format to control for file size. Can specify a numpy data type from https://docs.scipy.org/doc/numpy/user/basics.types.html.")

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
        type=int,
        default=multiprocessing.cpu_count(),
        help="""For local execution, how many subject-wise modelbuilds to run in parallel,
        defaults to number of CPUs""")
    g_ants_dbm.add_argument(
        '--template_reg_script',
        type=str,
        default='SyN',
        help="""Registration script that will be used for registration of the generated
        template to the provided atlas for masking and labeling. Can choose a predefined
        registration script among Rigid,Affine,SyN or light_SyN, or provide a custom script.""")


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
    #generates the parser CLI and execute the workflow based on specified parameters.
    opts = get_parser().parse_args()

    #obtain parser parameters
    bold_preproc_only=opts.bold_only
    bias_reg_script=opts.bias_reg_script
    coreg_script=define_reg_script(opts.coreg_script)
    commonspace_method=opts.commonspace_method
    data_dir_path=os.path.abspath(str(opts.input_dir))
    output_folder=os.path.abspath(str(opts.output_dir))
    plugin=opts.plugin
    os.environ["min_proc"]=str(opts.min_proc)

    #STC options
    stc_bool=opts.STC
    stc_TR=opts.TR
    if opts.tpattern=="alt":
        stc_tpattern='altplus'
    elif opts.tpattern=="seq":
        stc_tpattern='seqplus'
    else:
        raise ValueError('Invalid --tpattern provided.')

    #resampling options
    isotropic_resampling=opts.isotropic_resampling
    upsampling=opts.upsampling
    resampling_data_type=opts.data_type

    #setting absolute paths for ants_dbm options options
    os.environ["ants_dbm_cluster_type"]=opts.cluster_type
    os.environ["ants_dbm_walltime"]=opts.walltime
    os.environ["ants_dbm_memory_request"]=opts.memory_request
    os.environ["ants_dbm_local_threads"]=str(opts.local_threads)
    template_reg_option=opts.template_reg_script
    template_reg_script=define_reg_script(template_reg_option)

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

    data_csv=data_dir_path+'/data_info.csv' #this will be eventually replaced

    if bold_preproc_only:
        from rabies.preprocess_bold_pkg.bold_main_wf import init_EPIonly_bold_main_wf
        workflow = init_EPIonly_bold_main_wf(data_dir_path, data_csv, output_folder, tr=stc_TR, tpattern=stc_tpattern, apply_STC=stc_bool, bias_reg_script=bias_reg_script, coreg_script=coreg_script, template_reg_script=template_reg_script, isotropic_resampling=isotropic_resampling, upsampling=upsampling, resampling_data_type=resampling_data_type)
    elif not bold_preproc_only:
        from rabies.main_wf import init_unified_main_wf
        workflow = init_unified_main_wf(data_dir_path, data_csv, output_folder, tr=stc_TR, tpattern=stc_tpattern, commonspace_method=commonspace_method, template_reg_script=template_reg_script, apply_STC=stc_bool, bias_reg_script=bias_reg_script, coreg_script=coreg_script, isotropic_resampling=isotropic_resampling, upsampling=upsampling, resampling_data_type=resampling_data_type)
    else:
        raise ValueError('bold_preproc_only must be true or false.')

    workflow.base_dir = output_folder

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

    print('Running main workflow with %s plugin.' % plugin)
    #execute workflow, with plugin_args limiting the cluster load for parallel execution
    workflow.run(plugin=plugin, plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp %s' % (os.environ["min_proc"])})


def define_reg_script(reg_option):
    import rabies
    dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
    if reg_option=='SyN':
        reg_script=dir_path+'/shell_scripts/SyN_registration.sh'
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
