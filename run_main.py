'''
Pipeline Execution.
'''

import os
import sys
from mfp.main_wf import init_anat_init_wf, init_main_postPydpiper_wf

data_dir_path=os.path.abspath(sys.argv[1])
output_folder=os.path.abspath(sys.argv[2])
bold_preproc_only=str(sys.argv[3])
commonspace_method=sys.argv[4]
bias_reg_script=sys.argv[5]
coreg_script=sys.argv[6]
plugin=str(sys.argv[7])
debug=str(sys.argv[8])
atlas_file=sys.argv[9]

data_csv=data_dir_path+'/data_info.csv'
csv_labels='/data/chamal/projects/Gabriel_DG/atlases/DSURQE_atlas/labels/DSURQE_40micron_R_mapping.csv' #file with the id# of each label in the atlas to compute WM and CSF masks

if bold_preproc_only=='true':
    bold_preproc_only=True
    commonspace_transform=False
    compute_WM_CSF_masks=False
elif bold_preproc_only=='false':
    bold_preproc_only=False

    anat_init_wf = init_anat_init_wf(data_csv, data_dir_path, output_folder, commonspace_method=commonspace_method)
    anat_init_wf.base_dir = output_folder

    dir_path = os.path.dirname(os.path.realpath(__file__))
    if commonspace_method=='pydpiper':
        model_script_path=dir_path+'/mfp/shell_scripts/pydpiper.sh'
        commonspace_transform=False
        compute_WM_CSF_masks=True
    elif commonspace_method=='ants_dbm':
        model_script_path=dir_path+'/mfp/shell_scripts/ants_dbm.sh'
        commonspace_transform=True
        compute_WM_CSF_masks=False
    else:
        raise ValueError('Invalid commonspace method.')

    commonspace_csv_file=output_folder+'/anat_init_wf/commonspace_prep/commonspace_input_files.csv'

else:
    raise ValueError('bold_preproc_only must be true or false.')

main_postPydpiper_wf = init_main_postPydpiper_wf(data_csv, data_dir_path, output_folder, bold_preproc_only=bold_preproc_only, compute_WM_CSF_masks=compute_WM_CSF_masks, csv_labels=csv_labels, bias_reg_script=bias_reg_script, coreg_script=coreg_script, commonspace_transform=commonspace_transform)
main_postPydpiper_wf.base_dir = output_folder

if debug=='true':
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


if bold_preproc_only:
    if plugin=='local':
        print('Running main workflow locally.')
        main_postPydpiper_wf.run()
    if plugin=='parallel':
        print('Running main workflow in parallel.')
        main_postPydpiper_wf.run(plugin='SGEGraph', plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})

elif plugin=='local':
    print('Running locally.')
    print('Running anat init.')
    anat_init_wf.run()

    print('Running commonspace registration.')
    #run pydpiper
    out_dir=output_folder+'/commonspace/'
    os.system('mkdir -p %s' % (out_dir))
    cwd=os.getcwd()
    os.chdir(out_dir)
    os.system('bash %s %s %s' % (model_script_path,commonspace_csv_file, atlas_file))
    os.chdir(cwd)

    print('Running main workflow locally.')
    main_postPydpiper_wf.run()

elif plugin=='parallel':
    print('Running in parallel.')
    anat_init_wf.run(plugin='SGEGraph', plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})


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
    os.system('bash %s %s %s' % (model_script_path,commonspace_csv_file, atlas_file))
    os.chdir(cwd)

    print('Running main workflow in parallel.')
    main_postPydpiper_wf.run(plugin='SGEGraph', plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})
