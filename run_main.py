import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path)

data_dir_path=os.path.abspath(sys.argv[1])
output_folder=os.path.abspath(sys.argv[2])
data_csv=data_dir_path+'/data_info.csv'
csv_labels='/opt/quarantine/resources/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/mappings/DSURQE_40micron_R_mapping.csv'

from workflows.main_wf import init_main_wf
workflow = init_main_wf(data_csv, data_dir_path, output_folder, TR=1.2, keep_pydpiper_transforms=True, anat_only=False, csv_labels=csv_labels, mbm_script='default')


#Specify the base directory for the working directory
workflow.base_dir = output_folder

#workflow.run(plugin='SGE')
#workflow.run()
workflow.run(plugin='MultiProc')
#workflow.run(plugin='SGEGraph')
