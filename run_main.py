import os
import sys

data_dir_path=os.path.abspath(sys.argv[1])
output_folder=os.path.abspath(sys.argv[2])
bold_preproc_only=str(sys.argv[3])
csv_labels=os.path.abspath(sys.argv[4]) #file with the id# of each label in the atlas to compute WM and CSF masks
mbm_script=sys.argv[5]
bias_reg_script=sys.argv[6]
coreg_script=sys.argv[7]
plugin=str(sys.argv[8])
debug=str(sys.argv[9])

if bold_preproc_only=='true':
    bold_preproc_only=True
elif bold_preproc_only=='false':
    bold_preproc_only=False
else:
    raise ValueError('bold_preproc_only must be true or false.')


data_csv=data_dir_path+'/data_info.csv'


from mfp.main_wf import init_main_wf
workflow = init_main_wf(data_csv, data_dir_path, output_folder, bold_preproc_only=bold_preproc_only, keep_pydpiper_transforms=True, csv_labels=csv_labels, bias_reg_script=bias_reg_script, coreg_script=coreg_script, mbm_script=mbm_script)


#Specify the base directory for the working directory
workflow.base_dir = output_folder

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

if plugin=='local':
    print('Running locally.')
    workflow.run()
if plugin=='multiproc':
    print('Running locally with multiproc.')
    workflow.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
if plugin=='parallel':
    print('Running in parallel.')
    workflow.run(plugin='SGEGraph', plugin_args = {'max_jobs':50,'dont_resubmit_completed_jobs': True, 'qsub_args': '-pe smp 1'})
