import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from .preprocess_anat_pkg.anat_main_wf import init_anat_main_wf
from .preprocess_bold_pkg.bold_main_wf import init_bold_main_wf
from os.path import join as opj
from nipype.interfaces.io import SelectFiles, DataSink


def init_main_wf(data_csv, data_dir_path, output_folder, TR, csv_labels, mbm_script='default', name='main_wf'):

    workflow = pe.Workflow(name=name)
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_preproc', 'input_anat']),
        name='outputnode')

    import pandas as pd
    data_df=pd.read_csv(data_csv, sep=',')
    subject_list=data_df['subject_id'].values.tolist()
    session_list=data_df['session'].values.tolist()
    run_list=data_df['run'].values.tolist()

    #create a list of list of session numbers for each subject
    session_iter={}
    for i in range(len(subject_list)):
        session_iter[subject_list[i]] = list(range(1,int(session_list[i])+1))

    #create a dictionary with list of bold run numbers for each subject
    run_iter={}
    for i in range(len(subject_list)):
        run_iter[subject_list[i]] = list(range(1,int(run_list[i])+1))

    #create a infosource node for each information about the files
    infosource = pe.Node(niu.IdentityInterface(fields=['subject_id']),
                      name="infosource")
    infosource.iterables = [('subject_id', subject_list)]


    anat_file = opj('{subject_id}', 'ses-{session}', 'anat', '{subject_id}_ses-{session}_anat.nii.gz')
    anat_selectfiles = pe.Node(SelectFiles({'anat': anat_file},
                                   base_directory=data_dir_path),
                       name="anat_selectfiles")
    anat_selectfiles.itersource = ('infosource', 'subject_id')
    anat_selectfiles.iterables = [('session', session_iter)]

    bold_file = opj('{subject_id}', 'ses-{session}', 'bold', '{subject_id}_ses-{session}_run-{run}_bold.nii.gz')
    bold_selectfiles = pe.Node(SelectFiles({'bold': bold_file},
                                   base_directory=data_dir_path),
                       name="bold_selectfiles")
    bold_selectfiles.itersource = ('infosource', 'subject_id')
    bold_selectfiles.iterables = [('session', session_iter), ('run', run_iter)]


    # Datasink - creates output folder for important outputs
    datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="datasink"),
                    name="input_datasink")


    #connect the select files and iterables
    workflow.connect([
        (infosource, anat_selectfiles, [
            ("subject_id", "subject_id"),
            ]),
        (anat_selectfiles, datasink, [
            ("anat", "input_anat"),
            ]),
        (infosource, bold_selectfiles, [
            ("subject_id", "subject_id"),
            ]),
        (bold_selectfiles, datasink, [
            ("bold", "input_bold"),
            ]),
    ])

    #####LINKING THE MAIN WORKFLOWS
    anat_main_wf=init_anat_main_wf(csv_labels=csv_labels, mbm_script=mbm_script)

    workflow.connect([
        (anat_selectfiles, anat_main_wf, [
            ("anat", "inputnode.anat_file"),
            ]),
        (anat_main_wf, datasink, [
            ("outputnode.pydpyper_inputs", "pydpyper_inputs"),
            ("outputnode.anat_preproc", "anat_preproc"),
            ("outputnode.anat_mask", "anat_mask"),
            ("outputnode.WM_mask", "WM_mask"),
            ("outputnode.CSF_mask", "CSF_mask"),
            ("outputnode.anat_labels", "anat_labels"),
            ]),
    ])

    bold_main_wf=init_bold_main_wf(TR=TR, use_syn=True, apply_STC=True, iterative_N4=True, motioncorr_24params=True, apply_GSR=False)

    workflow.connect([
        (bold_selectfiles, bold_main_wf, [
            ("bold", "inputnode.bold_file"),
            ]),
        (anat_main_wf, bold_main_wf, [
            ("outputnode.anat_preproc", "inputnode.anat_preproc"),
            ("outputnode.anat_mask", "inputnode.anat_mask"),
            ("outputnode.WM_mask", "inputnode.WM_mask"),
            ("outputnode.CSF_mask", "inputnode.CSF_mask"),
            ("outputnode.anat_labels", "inputnode.anat_labels"),
            ]),
    ])

    return workflow
