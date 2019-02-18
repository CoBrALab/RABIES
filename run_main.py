import os
import sys

data_dir_path=os.path.abspath(sys.argv[1])
output_folder=os.path.abspath(sys.argv[2])
bias_reg_script=os.path.abspath(sys.argv[3])
coreg_script=os.path.abspath(sys.argv[4])
plugin=str(sys.argv[5])
debug=bool(sys.argv[6])

data_csv=data_dir_path+'/data_info.csv'
csv_labels='/opt/quarantine/resources/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/mappings/DSURQE_40micron_R_mapping.csv'

from mfp.main_wf import init_main_wf
workflow = init_main_wf(data_csv, data_dir_path, output_folder, TR=1.2, keep_pydpiper_transforms=True, anat_only=False, csv_labels=csv_labels, bias_reg_script=bias_reg_script, coreg_script=coreg_script, mbm_script='default')


#Specify the base directory for the working directory
workflow.base_dir = output_folder

print(plugin)
print(debug)
#workflow.run()
#workflow.run(plugin='SGEGraph', plugin_args = {'dont_resubmit_completed_jobs': True, qsub_args: '-l h_vmem=5G'})

'''
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
'''
