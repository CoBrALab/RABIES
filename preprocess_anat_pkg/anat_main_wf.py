import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from .anat_preproc import init_anat_preproc_wf
from .pydpyper import init_pydpyper_wf
from .anat_mask_prep import init_anat_mask_prep_wf
from os.path import join as opj
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.interfaces.utility import Function


from .test import init_test_wf
from .bold_iter_test import init_bold_iter_wf

def anat_main_wf(data_csv, data_dir_path, output_folder, model_script_path='default', name='anat_main_wf'):

    workflow = pe.Workflow(name=name)
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_preproc', 'input_anat']),
        name='outputnode')


    #Load all information about the dataset form the csv file
    import pandas as pd
    data_df=pd.read_csv(data_csv, sep=',')
    subject_list=data_df['subject_id'].values.tolist()
    session_list=data_df['session'].values.tolist()
    run_list=data_df['run'].values.tolist()

    #create a infosource node for each information about the files
    infosource = pe.Node(niu.IdentityInterface(fields=['subject_id', 'session', 'run']),
                      name="infosource")
    infosource.iterables = [('subject_id', subject_list), ('session', session_list), ('run', run_list)]


    anat_file = opj('{subject_id}', 'ses-{session_num}', 'anat', '{subject_id}_ses-{session_num}_anat.nii.gz')

    templates = {'anat': anat_file}

    anat_selectfiles = pe.Node(SelectFiles(templates,
                                   base_directory=data_dir_path),
                       name="anat_selectfiles")

    # Datasink - creates output folder for important outputs
    datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="datasink"),
                    name="input_datasink")

    anat_preproc_wf = init_anat_preproc_wf()

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
    nlin_transform = opj('mbm_atlasReg_processed','{subject_id}_preproc_anat','transforms','{subject_id}_preproc_anat__concat_lsq6_I_lsq6_lsq12_and_nlin.xfm')

    pydpyper_templates = {'labels': labels, 'nlin_mask': nlin_mask, 'nlin_transform': nlin_transform}

    pydpyper_selectfiles = pe.Node(SelectFiles(pydpyper_templates),
                       name="pydpyper_selectfiles")

    anat_mask_prep_wf=init_anat_mask_prep_wf()
    bold_iter_wf=init_bold_iter_wf(data_dir_path=data_dir_path)

#    test_wf=init_test_wf()

    workflow.connect([
        (infosource, anat_selectfiles, [
            ("subject_id", "subject_id"),
            ("session", "session_num"),
            ]),
        (anat_selectfiles, datasink, [
            ("anat", "input_anat"),
            ]),
        (infosource, bold_iter_wf, [
            ("subject_id", "inputnode.subject_id"),
            ("session", "inputnode.session"),
            ("run", "inputnode.run"),
            ]),
        (bold_iter_wf, datasink, [
            ("outputnode.bold_file", "input_bold"),
            ]),
        (anat_selectfiles, anat_preproc_wf, [("anat", "inputnode.anat_file")]),
        (anat_preproc_wf, pydpyper_prep, [("outputnode.preproc_anat", "file_list")]),
        (anat_preproc_wf, datasink, [("outputnode.preproc_anat", "pydpyper_inputs")]),
        (pydpyper_prep, pydpyper_wf, [("csv_file", "inputnode.csv_file")]),
        (infosource, pydpyper_selectfiles, [("subject_id", "subject_id")]),
        (pydpyper_wf, pydpyper_selectfiles, [("outputnode.pydpyper_directory", "base_directory")]),
        (infosource, anat_mask_prep_wf, [("subject_id", "inputnode.subject_id")]),
        (anat_preproc_wf, anat_mask_prep_wf, [("outputnode.preproc_anat", "inputnode.anat_preproc")]),
        (pydpyper_selectfiles, anat_mask_prep_wf, [
            ("labels", "inputnode.labels"),
            ("nlin_mask", "inputnode.nlin_mask"),
            ("nlin_transform", "inputnode.nlin_transform"),
            ]),
    ])


    return workflow

'''
        (selectfiles, test_wf, [("anat", "inputnode.anat")]),
        (infosource, test_wf, [("subject_id", "inputnode.subject_id")]),

        (pydpyper_selectfiles, anat_mask_prep_wf, [
            ("labels", "inputnode.labels"),
            ("nlin_mask", "inputnode.nlin_mask"),
            ("nlin_transform", "inputnode.nlin_transform"),
            ]),
'''

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
