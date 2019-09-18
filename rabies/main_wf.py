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

def file_reader(file_path):
    import os
    subject_id=os.path.basename(file_path).split('_ses-')[0]
    session=os.path.basename(file_path).split('_ses-')[1][0]
    return [subject_id, session]


def commonspace_prep_func(file_list, commonspace_method):
    print(commonspace_method)
    import os
    import pandas as pd
    cwd = os.getcwd()
    csv_path=cwd+'/commonspace_input_files.csv'
    anat_id=[]
    for sub_file_list in file_list:
        for file in sub_file_list:
            anat_id.append(file)
    if commonspace_method=='pydpiper':
        df = pd.DataFrame(data={"file": anat_id})
        df.to_csv(csv_path, sep=',',index=False)
    elif commonspace_method=='ants_dbm':
        df = pd.DataFrame(data=anat_id)
        df.to_csv(csv_path, header=False, sep=',',index=False)
    return csv_path

def mnc2nii(mnc_file):
    import os
    cwd = os.getcwd()
    basename=os.path.basename(mnc_file).split('.')[0]
    os.system('mnc2nii %s %s/%s.nii' % (mnc_file,cwd,basename))
    os.system('gzip *.nii')
    return '%s/%s.nii.gz' % (cwd,basename)

def move_pydpiper_transforms(file_list_buffer, pydpiper_directory, transform_csv, output_folder):
    import os
    os.system('mkdir -p %s/anat_datasink/pydpiper_transforms/' % (output_folder))
    os.system('mv %s %s/anat_datasink/pydpiper_transforms/' % (transform_csv, output_folder))
    os.system('mv %s/mbm_atlasReg_processed %s/anat_datasink/pydpiper_transforms/' % (pydpiper_directory, output_folder))
    return None


def init_anat_init_wf(data_csv, data_dir_path, output_folder, commonspace_method='pydpiper', name='anat_init_wf'):

    workflow = pe.Workflow(name=name)
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_preproc','csv_file']),
        name='outputnode')

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
                             container="anat_datasink"),
                    name="anat_datasink")

    #connect the select files and iterables
    workflow.connect([
        (infosource, anat_selectfiles, [
            ("subject_id", "subject_id"),
            ]),
        (anat_selectfiles, datasink, [
            ("anat", "input_anat"),
            ]),
    ])

    anat_preproc_wf = init_anat_preproc_wf()

    anat2nii = pe.Node(Function(input_names=['mnc_file'],
                              output_names=['nii_file'],
                              function=mnc2nii),
                     name='anat2nii')

    joinnode_buffer = pe.JoinNode(niu.IdentityInterface(fields=['file_list']),
                     name='joinnode_buffer',
                     joinsource='anat_selectfiles',
                     joinfield=['file_list'])

    commonspace_prep = pe.JoinNode(Function(input_names=['file_list', 'commonspace_method'],
                              output_names=['csv_file'],
                              function=commonspace_prep_func),
                     name='commonspace_prep',
                     joinsource='infosource',
                     joinfield=['file_list'])
    commonspace_prep.inputs.commonspace_method=commonspace_method

    #connect the anat workflow
    workflow.connect([
        (anat_selectfiles, anat_preproc_wf, [("anat", "inputnode.anat_file")]),
        (anat_preproc_wf, datasink, [("outputnode.preproc_anat", "pydpiper_inputs")]),
        (anat_preproc_wf, anat2nii, [("outputnode.preproc_anat", "mnc_file")]),
        (anat2nii, datasink, [("nii_file", 'anat_preproc')]),
        (anat_preproc_wf, outputnode, [("outputnode.preproc_anat", "anat_preproc")]),
        (joinnode_buffer, commonspace_prep, [("file_list", "file_list")]),
        (commonspace_prep, outputnode, [("csv_file", "csv_file")]),
    ])

    if commonspace_method=='pydpiper':
        workflow.connect([
            (anat_preproc_wf, joinnode_buffer, [("outputnode.preproc_anat", "file_list")]),
        ])
    elif commonspace_method=='ants_dbm':
        workflow.connect([
            (anat2nii, joinnode_buffer, [("nii_file", "file_list")]),
        ])

    return workflow


def init_main_postPydpiper_wf(data_csv, data_dir_path, output_folder, csv_labels, tr, tpattern, apply_STC=True, bold_preproc_only=False, commonspace_transform=False,compute_WM_CSF_masks=False,
                bias_cor_script='Default', bias_reg_script='Rigid', coreg_script='SyN', aCompCor_method='50%', name='main_wf'):

    workflow = pe.Workflow(name=name)
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_mask', 'anat_labels', 'WM_mask', 'CSF_mask','initial_bold_ref', 'bias_cor_bold', 'confounds_csv', 'itk_bold_to_anat', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv',
                'itk_anat_to_bold','boldref_warped2anat', 'native_corrected_bold', 'corrected_ref_bold', 'bold_brain_mask', 'bold_WM_mask', 'bold_CSF_mask', 'bold_labels', 'commonspace_bold']),
        name='outputnode')

    if not bold_preproc_only:
        ###############ANAT PREPROCESSING WORKFLOW
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

        # Datasink - creates output folder for important outputs
        datasink = pe.Node(DataSink(base_directory=output_folder,
                                 container="anat_datasink"),
                        name="anat_datasink")

        file_info = pe.Node(Function(input_names=['file_path'],
                                  output_names=['subject_id', 'session'],
                                  function=file_reader),
                         name='file_info')

        anat2nii = pe.Node(Function(input_names=['mnc_file'],
                                  output_names=['nii_file'],
                                  function=mnc2nii),
                         name='anat2nii')

        anat_nii = output_folder+'/'+opj('anat_datasink','anat_preproc','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_anat_preproc.nii.gz')
        labels = output_folder+'/'+opj('anat_datasink','anat_labels','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_anat_labels.nii.gz')
        mask = output_folder+'/'+opj('anat_datasink','anat_mask','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_anat_mask.nii.gz')

        if commonspace_transform:
            commonspace_affine = output_folder+'/'+opj('anat_datasink','commonspace_transforms','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_Affine.mat')
            commonspace_warp = output_folder+'/'+opj('anat_datasink','commonspace_transforms','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_Warp.nii.gz')
            commonspace_anat_template = output_folder+'/'+opj('anat_datasink','commonspace_template','commonspace_template.nii.gz')
            if not compute_WM_CSF_masks:
                WM_mask = output_folder+'/'+opj('anat_datasink','WM_mask','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_WM_mask.nii.gz')
                CSF_mask = output_folder+'/'+opj('anat_datasink','CSF_mask','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_CSF_mask.nii.gz')
                commonspace_templates = {'anat_nii':anat_nii,'labels': labels, 'mask': mask, 'WM_mask': WM_mask, 'CSF_mask': CSF_mask, 'commonspace_affine': commonspace_affine, 'commonspace_warp': commonspace_warp, 'commonspace_anat_template':commonspace_anat_template}
            else:
                commonspace_templates = {'anat_nii':anat_nii,'labels': labels, 'mask': mask, 'commonspace_affine': commonspace_affine, 'commonspace_warp': commonspace_warp, 'commonspace_anat_template':commonspace_anat_template}

        else:
            commonspace_templates = {'anat_nii':anat_nii,'labels': labels, 'mask': mask}

        commonspace_selectfiles = pe.Node(SelectFiles(commonspace_templates),
                           name="commonspace_selectfiles")
        commonspace_selectfiles.itersource = ('infosource', 'subject_id')
        commonspace_selectfiles.iterables = [('session', session_iter)]
        commonspace_selectfiles.base_directory=output_folder


        if compute_WM_CSF_masks:
            anat_mask_prep_wf=init_anat_mask_prep_wf(csv_labels=csv_labels)

            #connect the anat workflow
            workflow.connect([
                (infosource, commonspace_selectfiles, [
                    ("subject_id", "subject_id"),
                    ]),
                (commonspace_selectfiles, file_info, [("anat_nii", "file_path")]),
                (file_info, anat_mask_prep_wf, [
                    ("subject_id", "inputnode.subject_id"),
                    ("session", "inputnode.session")]),
                (commonspace_selectfiles, anat_mask_prep_wf, [
                    ("labels", "inputnode.labels"),
                    ]),
                (anat_mask_prep_wf, datasink, [
                    ("outputnode.eroded_WM_mask", "WM_mask"),
                    ("outputnode.eroded_CSF_mask", "CSF_mask"),
                    ]),
                (anat_mask_prep_wf, outputnode, [
                    ("outputnode.eroded_WM_mask", "WM_mask"),
                    ("outputnode.eroded_CSF_mask", "CSF_mask"),
                    ]),
            ])
        else:
            #connect the anat workflow
            workflow.connect([
                (infosource, commonspace_selectfiles, [
                    ("subject_id", "subject_id"),
                    ]),
                (commonspace_selectfiles, file_info, [("anat_nii", "file_path")]),
                (commonspace_selectfiles, outputnode, [
                    ("WM_mask", "WM_mask"),
                    ("CSF_mask", "CSF_mask"),
                    ]),
            ])


        ########BOLD PREPROCESSING WORKFLOW
        bold_main_wf=init_bold_main_wf(tr=tr, tpattern=tpattern, apply_STC=apply_STC, data_dir_path=data_dir_path, bias_cor_script=bias_cor_script, bias_reg_script=bias_reg_script, coreg_script=coreg_script, commonspace_transform=commonspace_transform, SyN_SDC=True, apply_STC=True, aCompCor_method=aCompCor_method)

        bold_file = opj('{subject_id}', 'ses-{session}', 'bold', '{subject_id}_ses-{session}_run-{run}_bold.nii.gz')
        bold_selectfiles = pe.Node(SelectFiles({'bold': bold_file},
                                       base_directory=data_dir_path),
                           name="bold_selectfiles")
        bold_selectfiles.itersource = ('infosource', 'subject_id')
        bold_selectfiles.iterables = [('run', run_iter)]

        # Datasink - creates output folder for important outputs
        bold_datasink = pe.Node(DataSink(base_directory=output_folder,
                                 container="bold_datasink"),
                        name="bold_datasink")


        if compute_WM_CSF_masks:
            workflow.connect([
                (anat_mask_prep_wf, bold_main_wf, [
                    ("outputnode.eroded_WM_mask", "inputnode.WM_mask"),
                    ("outputnode.eroded_CSF_mask", "inputnode.CSF_mask"),
                    ]),
            ])
        else:
            workflow.connect([
                (commonspace_selectfiles, bold_main_wf, [
                    ("WM_mask", "inputnode.WM_mask"),
                    ("CSF_mask", "inputnode.CSF_mask"),
                    ]),
            ])


        workflow.connect([
            (file_info, bold_selectfiles, [
                ("subject_id", "subject_id"),
                ("session", "session"),
                ]),
            (bold_selectfiles, bold_datasink, [
                ("bold", "input_bold"),
                ]),
            (bold_selectfiles, bold_main_wf, [
                ("bold", "inputnode.bold"),
                ]),
            (commonspace_selectfiles, bold_main_wf, [
                ("anat_nii", 'inputnode.anat_preproc'),
                ("labels", 'inputnode.labels'),
                ("mask", 'inputnode.anat_mask'),
                ]),
            (bold_main_wf, outputnode, [
                ("outputnode.bold_ref", "initial_bold_ref"),
                ("outputnode.corrected_EPI", "bias_cor_bold"),
                ("outputnode.EPI_brain_mask", "bold_brain_mask"),
                ("outputnode.EPI_WM_mask", "bold_WM_mask"),
                ("outputnode.EPI_CSF_mask", "bold_CSF_mask"),
                ("outputnode.EPI_labels", "bold_labels"),
                ("outputnode.confounds_csv", "confounds_csv"),
                ("outputnode.FD_voxelwise", "FD_voxelwise"),
                ("outputnode.pos_voxelwise", "pos_voxelwise"),
                ("outputnode.FD_csv", "FD_csv"),
                ("outputnode.itk_bold_to_anat", "itk_bold_to_anat"),
                ("outputnode.itk_anat_to_bold", "itk_anat_to_bold"),
                ("outputnode.output_warped_bold", "boldref_warped2anat"),
                ("outputnode.resampled_bold", "native_corrected_bold"),
                ("outputnode.resampled_ref_bold", "corrected_ref_bold"),
                ]),
            (outputnode, bold_datasink, [
                ("initial_bold_ref","initial_bold_ref"), #inspect initial bold ref
                ("bias_cor_bold","bias_cor_bold"), #inspect bias correction
                ("bold_brain_mask","bold_brain_mask"), #get the EPI labels
                ("bold_WM_mask","bold_WM_mask"), #get the EPI labels
                ("bold_CSF_mask","bold_CSF_mask"), #get the EPI labels
                ("bold_labels","bold_labels"), #get the EPI labels
                ("confounds_csv", "confounds_csv"), #confounds file
                ("FD_voxelwise", "FD_voxelwise"),
                ("pos_voxelwise", "pos_voxelwise"),
                ("FD_csv", "FD_csv"),
                ("itk_bold_to_anat", "itk_bold_to_anat"),
                ("itk_anat_to_bold", "itk_anat_to_bold"),
                ("boldref_warped2anat","boldref_warped2anat"), #warped EPI to anat
                ("native_corrected_bold", "native_corrected_bold"), #resampled EPI after motion realignment and SDC
                ("corrected_ref_bold", "corrected_ref_bold"), #resampled EPI after motion realignment and SDC
                ]),
        ])
        if commonspace_transform:
            workflow.connect([
                (commonspace_selectfiles, bold_main_wf, [
                    ('commonspace_affine', 'inputnode.commonspace_affine'),
                    ('commonspace_warp', 'inputnode.commonspace_warp'),
                    ('commonspace_anat_template', 'inputnode.commonspace_template'),
                    ]),
                (bold_main_wf, outputnode, [
                    ("outputnode.commonspace_bold", "commonspace_bold"),
                    ]),
                (outputnode, bold_datasink, [
                    ("commonspace_bold", "commonspace_corrected_bold"), #resampled EPI after motion realignment and SDC
                    ]),
            ])


    elif bold_preproc_only:
        print("BOLD preproc only!")
        import pandas as pd
        data_df=pd.read_csv(data_csv, sep=',')
        subject_list=data_df['subject_id'].values.tolist()

        infosource = pe.Node(niu.IdentityInterface(fields=['subject_id']),
                          name="infosource")
        infosource.iterables = [('subject_id', subject_list)]

        # Datasink - creates output folder for important outputs
        bold_datasink = pe.Node(DataSink(base_directory=output_folder,
                                 container="bold_datasink"),
                        name="bold_datasink")

        bold_main_wf=init_bold_main_wf(bold_preproc_only=True, tr=tr, tpattern=tpattern, apply_STC=apply_STC, data_csv=data_csv, data_dir_path=data_dir_path, bias_reg_script=bias_reg_script, coreg_script=coreg_script, SyN_SDC=True, apply_STC=True, aCompCor_method=aCompCor_method)

        workflow.connect([
            (infosource, bold_main_wf, [
                ("subject_id", "inputnode.subject_id"),
                ]),
            (bold_main_wf, outputnode, [
                ("outputnode.bold_ref", "initial_bold_ref"),
                ("outputnode.corrected_EPI", "bias_cor_bold"),
                ("outputnode.EPI_brain_mask", "bold_brain_mask"),
                ("outputnode.EPI_WM_mask", "bold_WM_mask"),
                ("outputnode.EPI_CSF_mask", "bold_CSF_mask"),
                ("outputnode.EPI_labels", "bold_labels"),
                ("outputnode.confounds_csv", "confounds_csv"),
                ("outputnode.FD_voxelwise", "FD_voxelwise"),
                ("outputnode.pos_voxelwise", "pos_voxelwise"),
                ("outputnode.FD_csv", "FD_csv"),
                ("outputnode.itk_bold_to_anat", "itk_bold_to_anat"),
                ("outputnode.itk_anat_to_bold", "itk_anat_to_bold"),
                ("outputnode.output_warped_bold", "boldref_warped2anat"),
                ("outputnode.resampled_bold", "native_corrected_bold"),
                ("outputnode.resampled_ref_bold", "corrected_ref_bold"),
                ]),
            (outputnode, bold_datasink, [
                ("initial_bold_ref","initial_bold_ref"), #inspect initial bold ref
                ("bias_cor_bold","bias_cor_bold"), #inspect bias correction
                ("bold_brain_mask","bold_brain_mask"), #get the EPI labels
                ("bold_WM_mask","bold_WM_mask"), #get the EPI labels
                ("bold_CSF_mask","bold_CSF_mask"), #get the EPI labels
                ("bold_labels","bold_labels"), #get the EPI labels
                ("confounds_csv", "confounds_csv"), #confounds file
                ("FD_voxelwise", "FD_voxelwise"),
                ("pos_voxelwise", "pos_voxelwise"),
                ("FD_csv", "FD_csv"),
                ("itk_bold_to_anat", "itk_bold_to_anat"),
                ("itk_anat_to_bold", "itk_anat_to_bold"),
                ("boldref_warped2anat","boldref_warped2anat"), #warped EPI to anat
                ("native_corrected_bold", "native_corrected_bold"), #resampled EPI after motion realignment and SDC
                ("corrected_ref_bold", "corrected_ref_bold"), #resampled EPI after motion realignment and SDC
                ]),
        ])


    return workflow
