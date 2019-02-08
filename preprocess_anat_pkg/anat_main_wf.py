import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from .anat_preproc import init_anat_preproc_wf
from .pydpyper import init_pydpyper_wf
from .anat_mask_prep import init_anat_mask_prep_wf
from os.path import join as opj
from nipype.interfaces.io import SelectFiles, DataSink

from nipype.interfaces.utility import Function, IdentityInterface
def anat_main_wf(subject_csv, data_dir_path, output_folder, model_script_path='default', name='anat_main_wf'):

    workflow = pe.Workflow(name=name)
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_preproc', 'input_anat']),
        name='outputnode')


    import pandas as pd
    subject_list=pd.read_csv(subject_csv)['subject_id'].values.tolist()

    infosource = pe.Node(IdentityInterface(fields=['subject_id']),
                      name="infosource")
    infosource.iterables = [('subject_id', subject_list)]

    anat_file = opj('{subject_id}', 'anat', '{subject_id}_anat.nii.gz')
    bold_file = opj('{subject_id}', 'bold', '{subject_id}_bold.nii.gz')

    templates = {'anat': anat_file, 'bold': bold_file}

    selectfiles = pe.Node(SelectFiles(templates,
                                   base_directory=data_dir_path),
                       name="selectfiles")

    # Datasink - creates output folder for important outputs
    input_datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="datasink"),
                    name="input_datasink")

    anat_preproc_wf = init_anat_preproc_wf()

    # Datasink - creates output folder for important outputs
    pydpyper_datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="datasink"),
                    name="pydpyper_datasink")

    pydpyper_prep = pe.JoinNode(Function(input_names=['file_list'],
                              output_names=['csv_file'],
                              function=pydpyper_prep_func),
                     name='pydpyper_prep',
                     joinsource=infosource,
                     joinfield=['file_list'])


    if model_script_path=='default':
        dir_path = os.path.dirname(os.path.realpath(__file__))
        pydpyper_wf=init_pydpyper_wf(model_script_path=dir_path+'/utils_scripts/mbm_template.sh')
    else:
        pydpyper_wf=init_pydpyper_wf(model_script_path=model_script_path)


    labels = opj('mbm_atlasReg_processed','{subject_id}_preproc_anat','voted.mnc')
    nlin_mask = opj('mbm_atlasReg_processed','DSURQE_40micron_mask','resampled','{subject_id}_preproc_anat_I_lsq6_lsq12_and_nlin-resampled_mask.mnc')

    pydpyper_templates = {'labels': labels, 'nlin_mask': nlin_mask}

    pydpyper_selectfiles = pe.Node(SelectFiles(pydpyper_templates),
                       name="pydpyper_selectfiles")


    anat_mask_prep_wf=init_anat_mask_prep_wf(subject_csv=subject_csv)

    workflow.connect([
        (infosource, selectfiles, [("subject_id", "subject_id")]),
        (selectfiles, outputnode, [("anat", "input_anat"),
                                 ("bold", "input_bold")]),
        (selectfiles, input_datasink, [("anat", "input_datasink")]),
        (selectfiles, anat_preproc_wf, [("anat", "inputnode.anat_file")]),
        (anat_preproc_wf, pydpyper_prep, [("outputnode.preproc_anat", "file_list")]),
        (anat_preproc_wf, pydpyper_datasink, [("outputnode.preproc_anat", "pydpyper_inputs")]),
        (pydpyper_prep, pydpyper_wf, [("csv_file", "inputnode.csv_file")]),
        (infosource, pydpyper_selectfiles, [("subject_id", "subject_id")]),
        (pydpyper_wf, pydpyper_selectfiles, [("outputnode.pydpyper_directory", "base_directory")]),
    ])


    return workflow


def pydpyper_prep_func(file_list):
    import os
    import pandas as pd
    cwd = os.getcwd()
    csv_path=cwd+'/pydpyper_input_files.csv'
    anat_id=[]
    for file in file_list:
        anat_id.append(file)
    df = pd.DataFrame(data={"file": anat_id})
    df.to_csv(csv_path, sep=',',index=False)
    return csv_path
