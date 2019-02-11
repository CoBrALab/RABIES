import os
from os.path import join as opj
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from .preprocess_anat_pkg.anat_preproc import init_anat_preproc_wf
from .preprocess_anat_pkg.pydpiper import init_pydpiper_wf
from .preprocess_anat_pkg.anat_mask_prep import init_anat_mask_prep_wf
from .preprocess_bold_pkg.bold_main_wf import init_bold_main_wf
from nipype.interfaces.io import SelectFiles, DataSink

from nipype.interfaces.utility import Function


def init_main_wf(data_csv, data_dir_path, output_folder, TR, csv_labels, anat_only=True, mbm_script='default', name='main_wf'):

    workflow = pe.Workflow(name=name)


    import pandas as pd
    data_df=pd.read_csv(data_csv, sep=',')
    subject_list=data_df['subject_id'].values.tolist()
    session_list=data_df['num_session'].values.tolist()
    run_list=data_df['num_run'].values.tolist()

    #create a dictionary with list of bold session numbers for each subject
    session_iter={}
    for i in range(len(subject_list)):
        session_iter[subject_list[i]] = list(range(1,int(session_list[i])+1))


    #create a dictionary with list of bold run numbers for each subject
    run_iter={}
    for i in range(len(subject_list)):
        run_iter[subject_list[i]] = list(range(1,int(run_list[i])+1))

    infosource = pe.Node(niu.IdentityInterface(fields=['subject_id']),
                      name="infosource")
    infosource.iterables = [('subject_id', subject_list)]


    anat_file = opj('{subject_id}', 'ses-{session}', 'anat', '{subject_id}_ses-{session}_anat.nii.gz')
    anat_selectfiles = pe.Node(SelectFiles({'anat': anat_file},
                                   base_directory=data_dir_path),
                       name="anat_selectfiles")
    anat_selectfiles.itersource = ('infosource', 'subject_id')
    anat_selectfiles.iterables = [('session', session_iter)]


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
    ])

    ###############ANAT PREPROCESSING WORKFLOW
    file_info = pe.Node(Function(input_names=['file_path'],
                              output_names=['subject_id', 'session'],
                              function=file_reader),
                     name='file_info')

    anat_preproc_wf = init_anat_preproc_wf()

    anat2nii = pe.Node(Function(input_names=['mnc_file'],
                              output_names=['nii_file'],
                              function=mnc2nii),
                     name='anat2nii')

    pydpiper_prep = pe.JoinNode(Function(input_names=['file_list'],
                              output_names=['csv_file'],
                              function=pydpiper_prep_func),
                     name='pydpiper_prep',
                     joinsource='anat_preproc_wf',
                     joinfield=['file_list'])

    if mbm_script=='default':
        dir_path = os.path.dirname(os.path.realpath(__file__))
        pydpiper_wf=init_pydpiper_wf(model_script_path=dir_path+'/preprocess_anat_pkg/shell_scripts/mbm_template.sh')
    else:
        pydpiper_wf=init_pydpiper_wf(model_script_path=model_script_path)


    labels = opj('mbm_atlasReg_processed','{subject_id}_ses-{session}_anat_preproc','voted.mnc')
    nlin_mask = opj('mbm_atlasReg_processed','DSURQE_40micron_mask','resampled','{subject_id}_ses-{session}_anat_preproc_I_lsq6_lsq12_and_nlin-resampled_mask.mnc')
    nlin_transform = opj('mbm_atlasReg_processed','{subject_id}_ses-{session}_anat_preproc','transforms','{subject_id}_ses-{session}_anat_preproc__concat_lsq6_I_lsq6_lsq12_and_nlin.xfm')

    pydpiper_templates = {'labels': labels, 'nlin_mask': nlin_mask, 'nlin_transform': nlin_transform}

    pydpiper_selectfiles = pe.Node(SelectFiles(pydpiper_templates),
                       name="pydpiper_selectfiles")

    anat_mask_prep_wf=init_anat_mask_prep_wf(csv_labels=csv_labels)


    workflow.connect([
        (anat_selectfiles, file_info, [("anat", "file_path")]),
        (anat_selectfiles, anat_preproc_wf, [("anat", "inputnode.anat_file")]),
        (anat_preproc_wf, datasink, [("outputnode.preproc_anat", "pydpiper_inputs")]),
        (anat_preproc_wf, anat2nii, [("outputnode.preproc_anat", "mnc_file")]),
        (anat2nii, datasink, [("nii_file", 'anat_preproc')]),
        (anat_preproc_wf, pydpiper_prep, [("outputnode.preproc_anat", "file_list")]),
        (pydpiper_prep, pydpiper_wf, [("csv_file", "inputnode.csv_file")]),
        (file_info, pydpiper_selectfiles, [
            ("subject_id", "subject_id"),
            ("session", "session"),
            ]),
        (pydpiper_wf, pydpiper_selectfiles, [("outputnode.pydpiper_directory", "base_directory")]),
        (file_info, anat_mask_prep_wf, [
            ("subject_id", "inputnode.subject_id"),
            ("session", "inputnode.session")]),
        (anat_preproc_wf, anat_mask_prep_wf, [("outputnode.preproc_anat", "inputnode.anat_preproc")]),
        (pydpiper_selectfiles, anat_mask_prep_wf, [
            ("labels", "inputnode.labels"),
            ("nlin_mask", "inputnode.nlin_mask"),
            ("nlin_transform", "inputnode.nlin_transform"),
            ]),
        (anat_mask_prep_wf, datasink, [
            ("outputnode.resampled_mask", "anat_mask"),
            ("outputnode.resampled_labels", "anat_labels"),
            ("outputnode.eroded_WM_mask", "WM_mask"),
            ("outputnode.eroded_CSF_mask", "CSF_mask"),
            ]),
    ])


    ########BOLD PREPROCESSING WORKFLOW
    if not anat_only:
        bold_main_wf=init_bold_main_wf(data_dir_path=data_dir_path, run_iter=run_iter, TR=TR, use_syn=True, apply_STC=True, iterative_N4=True, motioncorr_24params=True, apply_GSR=False)

        bold_file = opj('{subject_id}', 'ses-{session}', 'bold', '{subject_id}_ses-{session}_run-{run}_bold.nii.gz')
        bold_selectfiles = pe.Node(SelectFiles({'bold': bold_file},
                                       base_directory=data_dir_path),
                           name="bold_selectfiles")
        bold_selectfiles.itersource = ('infosource', 'subject_id')
        bold_selectfiles.iterables = [('run', run_iter)]


        workflow.connect([
            (file_info, bold_selectfiles, [
                ("subject_id", "subject_id"),
                ("session", "session"),
                ]),
            (bold_selectfiles, datasink, [
                ("bold", "input_bold"),
                ]),
            (bold_selectfiles, bold_main_wf, [
                ("bold", "inputnode.bold"),
                ]),
            (anat2nii, bold_main_wf, [("nii_file", 'inputnode.anat_preproc')]),
            (anat_mask_prep_wf, bold_main_wf, [
                ("outputnode.resampled_mask", "inputnode.anat_mask"),
                ("outputnode.resampled_labels", "inputnode.labels"),
                ("outputnode.eroded_WM_mask", "inputnode.WM_mask"),
                ("outputnode.eroded_CSF_mask", "inputnode.CSF_mask"),
                ]),
        ])

    return workflow


def file_reader(file_path):
    import os
    subject_id=os.path.basename(file_path).split('_ses-')[0]
    session=os.path.basename(file_path).split('_ses-')[1][0]
    return [subject_id, session]


def pydpiper_prep_func(file_list):
    import os
    import pandas as pd
    cwd = os.getcwd()
    csv_path=cwd+'/pydpiper_input_files.csv'
    anat_id=[]
    for file in file_list:
        anat_id.append(file)
    df = pd.DataFrame(data={"file": anat_id})
    df.to_csv(csv_path, sep=',',index=False)
    return csv_path

def mnc2nii(mnc_file):
    import os
    cwd = os.getcwd()
    basename=os.path.basename(mnc_file).split('.')[0]
    os.system('mnc2nii %s %s/%s.nii' % (mnc_file,cwd,basename))
    os.system('gzip *.nii')
    return '%s/%s.nii.gz' % (cwd,basename)
