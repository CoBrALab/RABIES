import sys
import os

bold_file=os.path.abspath(sys.argv[1])
name=sys.argv[2]
anat=sys.argv[3]
brain_mask=sys.argv[4]
labels=sys.argv[5]
WM_mask=sys.argv[6]
CSF_mask=sys.argv[7]
output_folder=sys.argv[8]

from preprocess_bold_pkg.main_wf import init_func_preproc_wf


workflow = init_func_preproc_wf(bold_file=bold_file, use_syn=True, TR='1.2s', iterative_N4=True,
				apply_GSR=True, name=name)


#set inputs

workflow.inputs.inputnode.bold_file = bold_file
workflow.inputs.inputnode.anat_preproc = anat
workflow.inputs.inputnode.anat_mask = brain_mask
workflow.inputs.inputnode.anat_labels = labels
workflow.inputs.inputnode.CSF_mask = CSF_mask
workflow.inputs.inputnode.WM_mask = WM_mask

#Specify the base directory for the working directory
workflow.base_dir = output_folder

workflow.run(plugin='SGEGraph', plugin_args = {'dont_resubmit_completed_jobs': True})
