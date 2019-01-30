import sys
import os
sys.path.append('/home/cic/desgab/Desktop/mouse_procfmri')

bold_file=os.path.abspath(sys.argv[1])
name=sys.argv[2]
anat=sys.argv[3]
brain_mask=sys.argv[4]
labels=sys.argv[5]
WM_mask=sys.argv[6]
CSF_mask=sys.argv[7]
output_folder=sys.argv[8]

from preprocess_bold_pkg.main_wf import init_func_preproc_wf


workflow = init_func_preproc_wf(bold_file=bold_file, use_syn=True, TR=1.2, iterative_N4=True,
				apply_GSR=True, apply_STC=True, motioncorr_24params=False, name=name)


#set inputs

workflow.inputs.inputnode.bold_file = bold_file
workflow.inputs.inputnode.anat_preproc = anat
workflow.inputs.inputnode.anat_mask = brain_mask
workflow.inputs.inputnode.anat_labels = labels
workflow.inputs.inputnode.CSF_mask = CSF_mask
workflow.inputs.inputnode.WM_mask = WM_mask

#Specify the base directory for the working directory
workflow.base_dir = output_folder

#from nipype import config, logging
#config.enable_debug_mode()
#workflow.run(plugin='SGE')
workflow.run()
#workflow.run(plugin='MultiProc')
#config the execution
'''
from nipype import config, logging

config_dict={'execution': {'stop_on_first_crash': 'true',
						   'plugin': 'SGEGraph',
						   'remove_unnecessary_outputs': 'false',
                           'keep_inputs': 'false',
                           'crashdump_dir': '/data/chamal/projects/Gabriel_DG/debug/crash_folder',
                           'stop_on_first_rerun': 'false',
                           'create_report': 'true'},
             'logging': {'workflow_level': 'INFO',
                         'log_directory': '/data/chamal/projects/Gabriel_DG/debug/log_folder',
                         'interface_level': 'INFO',
                         'log_to_file': 'true'}}

config.update_config(config_dict)
logging.update_logging(config)

config.enable_debug_mode()
workflow.run(plugin='SGEGraph')
#workflow.run(plugin='SGEGraph', plugin_args = {'dont_resubmit_completed_jobs': True})plugin_args=dict(qsub_args='--ppj 4 --mem 10')
'''
