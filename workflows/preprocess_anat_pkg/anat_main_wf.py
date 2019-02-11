import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from .anat_preproc import init_anat_preproc_wf
from .pydpyper import init_pydpyper_wf
from .anat_mask_prep import init_anat_mask_prep_wf
from os.path import join as opj
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.interfaces.utility import Function


def init_anat_main_wf(csv_labels, mbm_script='default', name='anat_main_wf'):

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_file']),
        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_preproc', 'pydpyper_inputs', 'anat_mask',
                'anat_labels', 'WM_mask', 'CSF_mask']),
        name='outputnode')

    file_info = pe.Node(Function(input_names=['file_path'],
                              output_names=['subject_id', 'session'],
                              function=file_reader),
                     name='file_info')

    anat_preproc_wf = init_anat_preproc_wf()

    anat2nii = pe.Node(Function(input_names=['mnc_file'],
                              output_names=['nii_file'],
                              function=mnc2nii),
                     name='anat2nii')

    pydpyper_prep = pe.JoinNode(Function(input_names=['file_list'],
                              output_names=['csv_file'],
                              function=pydpyper_prep_func),
                     name='pydpyper_prep',
                     joinsource=anat_preproc_wf,
                     joinfield=['file_list'])

    if mbm_script=='default':
        dir_path = os.path.dirname(os.path.realpath(__file__))
        pydpyper_wf=init_pydpyper_wf(model_script_path=dir_path+'/utils_scripts/mbm_template.sh')
    else:
        pydpyper_wf=init_pydpyper_wf(model_script_path=model_script_path)


    labels = opj('mbm_atlasReg_processed','{subject_id}_ses-{session}_anat_preproc','voted.mnc')
    nlin_mask = opj('mbm_atlasReg_processed','DSURQE_40micron_mask','resampled','{subject_id}_ses-{session}_anat_preproc_I_lsq6_lsq12_and_nlin-resampled_mask.mnc')
    nlin_transform = opj('mbm_atlasReg_processed','{subject_id}_ses-{session}_anat_preproc','transforms','{subject_id}_ses-{session}_anat_preproc__concat_lsq6_I_lsq6_lsq12_and_nlin.xfm')

    pydpyper_templates = {'labels': labels, 'nlin_mask': nlin_mask, 'nlin_transform': nlin_transform}

    pydpyper_selectfiles = pe.Node(SelectFiles(pydpyper_templates),
                       name="pydpyper_selectfiles")

    anat_mask_prep_wf=init_anat_mask_prep_wf(csv_labels=csv_labels)

    #ANAT PREPROCESSING WORKFLOW
    workflow.connect([
        (inputnode, file_info, [("anat_file", "file_path")]),
        (inputnode, anat_preproc_wf, [("anat_file", "inputnode.anat_file")]),
        (anat_preproc_wf, pydpyper_prep, [("outputnode.preproc_anat", "file_list")]),
        (anat_preproc_wf, outputnode, [("outputnode.preproc_anat", "pydpyper_inputs")]),
        (anat_preproc_wf, anat2nii, [("outputnode.preproc_anat", "mnc_file")]),
        (anat2nii, outputnode, [("nii_file", 'anat_preproc')]),
        (pydpyper_prep, pydpyper_wf, [("csv_file", "inputnode.csv_file")]),
        (file_info, pydpyper_selectfiles, [
            ("subject_id", "subject_id"),
            ("session", "session"),
            ]),
        (pydpyper_wf, pydpyper_selectfiles, [("outputnode.pydpyper_directory", "base_directory")]),
        (file_info, anat_mask_prep_wf, [
            ("subject_id", "inputnode.subject_id"),
            ("session", "inputnode.session")]),
        (anat_preproc_wf, anat_mask_prep_wf, [("outputnode.preproc_anat", "inputnode.anat_preproc")]),
        (pydpyper_selectfiles, anat_mask_prep_wf, [
            ("labels", "inputnode.labels"),
            ("nlin_mask", "inputnode.nlin_mask"),
            ("nlin_transform", "inputnode.nlin_transform"),
            ]),
        (anat_mask_prep_wf, outputnode, [
            ("outputnode.resampled_mask", "anat_mask"),
            ("outputnode.resampled_labels", "anat_labels"),
            ("outputnode.eroded_WM_mask", "WM_mask"),
            ("outputnode.eroded_CSF_mask", "CSF_mask"),
            ]),
    ])

    return workflow

def file_reader(file_path):
    import os
    subject_id=os.path.basename(file_path).split('_ses-')[0]
    session=os.path.basename(file_path).split('_ses-')[1][0]
    return [subject_id, session]


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

def mnc2nii(mnc_file):
    import os
    cwd = os.getcwd()
    basename=os.path.basename(mnc_file).split('.')[0]
    os.system('mnc2nii %s %s/%s.nii' % (mnc_file,cwd,basename))
    os.system('gzip *.nii')
    return '%s/%s.nii.gz' % (cwd,basename)
