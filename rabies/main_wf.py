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

def init_unified_main_wf(data_dir_path, data_csv, output_folder, tr, tpattern, commonspace_method='pydpiper', apply_STC=True,
                bias_cor_script='Default', bias_reg_script='Rigid', coreg_script='SyN', aCompCor_method='50%', name='main_wf'):

    workflow = pe.Workflow(name=name)

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_preproc','anat_mask', 'anat_labels', 'WM_mask', 'CSF_mask','initial_bold_ref', 'bias_cor_bold', 'confounds_csv', 'itk_bold_to_anat', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv',
                'itk_anat_to_bold','boldref_warped2anat', 'native_corrected_bold', 'corrected_ref_bold', 'bold_brain_mask', 'bold_WM_mask', 'bold_CSF_mask', 'bold_labels', 'commonspace_bold', 'commonspace_mask', 'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_labels']),
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

    ####setting up all iterables
    infosub_id = pe.Node(niu.IdentityInterface(fields=['subject_id']),
                      name="infosub_id")
    infosub_id.iterables = [('subject_id', subject_list)]

    infosession = pe.Node(niu.IdentityInterface(fields=['session', 'subject_id']),
                      name="infosession")
    infosession.itersource = ('infosub_id', 'subject_id')
    infosession.iterables = [('session', session_iter)]

    inforun = pe.Node(niu.IdentityInterface(fields=['run', 'subject_id']),
                      name="inforun")
    inforun.itersource = ('infosub_id', 'subject_id')
    inforun.iterables = [('run', run_iter)]


    anat_file = opj('{subject_id}', 'ses-{session}', 'anat', '{subject_id}_ses-{session}_anat.nii.gz')
    anat_selectfiles = pe.Node(SelectFiles({'anat': anat_file},
                                   base_directory=data_dir_path),
                       name="anat_selectfiles")

    anat_preproc_wf = init_anat_preproc_wf()

    anat2nii = pe.Node(Function(input_names=['mnc_file'],
                              output_names=['nii_file'],
                              function=mnc2nii),
                     name='anat2nii')

    bold_file = opj('{subject_id}', 'ses-{session}', 'bold', '{subject_id}_ses-{session}_run-{run}_bold.nii.gz')
    bold_selectfiles = pe.Node(SelectFiles({'bold': bold_file},
                                   base_directory=data_dir_path),
                       name="bold_selectfiles")

    # Datasink - creates output folder for important outputs
    datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="datasink"),
                    name="datasink")

    #####setting up commonspace registration within the workflow
    joinnode_session = pe.JoinNode(niu.IdentityInterface(fields=['file_list']),
                     name='joinnode_session',
                     joinsource='infosession',
                     joinfield=['file_list'])

    joinnode_sub_id = pe.JoinNode(niu.IdentityInterface(fields=['file_list']),
                     name='joinnode_sub_id',
                     joinsource='infosub_id',
                     joinfield=['file_list'])


    #connect iterables, joinnodes and SelectFiles
    workflow.connect([
        (infosub_id, infosession, [
            ("subject_id", "subject_id"),
            ]),
        (infosub_id, inforun, [
            ("subject_id", "subject_id"),
            ]),
        (infosub_id, anat_selectfiles, [
            ("subject_id", "subject_id"),
            ]),
        (infosession, anat_selectfiles, [
            ("session", "session")
            ]),
        (joinnode_session, joinnode_sub_id, [
            ("file_list", "file_list"),
            ]),
        (infosub_id, bold_selectfiles, [
            ("subject_id", "subject_id"),
            ]),
        (infosession, bold_selectfiles, [
            ("session", "session")
            ]),
        (inforun, bold_selectfiles, [
            ("run", "run")
            ]),
        ])


    if commonspace_method=='pydpiper':
        workflow.connect([
            (anat_preproc_wf, joinnode_session, [("outputnode.preproc_anat", "file_list")]),
        ])
    elif commonspace_method=='ants_dbm':

        commonspace_reg = pe.Node(Function(input_names=['file_list', 'output_folder'],
                                  output_names=['ants_dbm_common', 'common_to_template_transform', 'template_to_common_transform'],
                                  function=commonspace_reg_function),
                         name='commonspace_reg')
        commonspace_reg.inputs.output_folder = output_folder+'/datasink/'

        #setting SelectFiles for the commonspace registration
        ants_dbm_inverse_warp = output_folder+'/'+opj('datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_{subject_id}_ses-{session}_anat_preproc*1InverseWarp.nii.gz')
        ants_dbm_warp = output_folder+'/'+opj('datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_{subject_id}_ses-{session}_anat_preproc*1Warp.nii.gz')
        ants_dbm_affine = output_folder+'/'+opj('datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_{subject_id}_ses-{session}_anat_preproc*0GenericAffine.mat')
        ants_dbm_common_anat = output_folder+'/'+opj('datasink','ants_dbm_outputs','{ants_dbm_common}')
        common_to_template_transform = output_folder+'/'+opj('datasink','ants_dbm_outputs','{common_to_template_transform}')
        template_to_common_transform = output_folder+'/'+opj('datasink','ants_dbm_outputs','{template_to_common_transform}')

        commonspace_templates = {'ants_dbm_inverse_warp':ants_dbm_inverse_warp,'ants_dbm_warp': ants_dbm_warp, 'ants_dbm_affine': ants_dbm_affine, 'ants_dbm_common_anat': ants_dbm_common_anat, 'common_to_template_transform': common_to_template_transform, 'template_to_common_transform':template_to_common_transform}

        commonspace_selectfiles = pe.Node(SelectFiles(commonspace_templates),
                       name="commonspace_selectfiles")

        def transform_masks(reference_image,ants_dbm_inverse_warp, ants_dbm_affine,template_to_common_transform):
            import os
            cwd = os.getcwd()
            subject_id=os.path.basename(reference_image).split('_ses-')[0]
            session=os.path.basename(reference_image).split('_ses-')[1][0]
            filename_template = '%s_ses-%s' % (subject_id, session)
            input_image=os.environ["template_mask"]
            brain_mask='%s/%s_ses-%s_%s' % (cwd, filename_template, 'anat_mask.nii.gz')
            os.system('antsApplyTransforms -d 3 -i %s -t %s -t [%s,1] -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,ants_dbm_inverse_warp, ants_dbm_affine,template_to_common_transform,reference_image,brain_mask,))
            input_image=os.environ["WM_mask"]
            WM_mask='%s/%s_ses-%s_%s' % (cwd, filename_template, 'WM_mask.nii.gz')
            os.system('antsApplyTransforms -d 3 -i %s -t %s -t [%s,1] -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,ants_dbm_inverse_warp, ants_dbm_affine,template_to_common_transform,reference_image,WM_mask,))
            input_image=os.environ["CSF_mask"]
            CSF_mask='%s/%s_ses-%s_%s' % (cwd, filename_template, 'CSF_mask.nii.gz')
            os.system('antsApplyTransforms -d 3 -i %s -t %s -t [%s,1] -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,ants_dbm_inverse_warp, ants_dbm_affine,template_to_common_transform,reference_image,CSF_mask,))
            input_image=os.environ["atlas_labels"]
            anat_labels='%s/%s_ses-%s_%s' % (cwd, filename_template, 'atlas_labels.nii.gz')
            os.system('antsApplyTransforms -d 3 -i %s -t %s -t [%s,1] -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,ants_dbm_inverse_warp, ants_dbm_affine,template_to_common_transform,reference_image,anat_labels,))

            return brain_mask, WM_mask, CSF_mask, anat_labels

        transform_masks = pe.Node(Function(input_names=['reference_image','ants_dbm_inverse_warp', 'ants_dbm_affine','template_to_common_transform'],
                                  output_names=['brain_mask', 'WM_mask', 'CSF_mask', 'anat_labels'],
                                  function=transform_masks),
                         name='transform_masks')

        def commonspace_transforms(ants_dbm_warp, ants_dbm_affine,common_to_template_transform):
            return [ants_dbm_warp, ants_dbm_affine,common_to_template_transform],[0,0,0] #transforms_list,inverses
        commonspace_transforms_prep = pe.Node(Function(input_names=['ants_dbm_warp','ants_dbm_affine','common_to_template_transform'],
                                  output_names=['transforms_list','inverses'],
                                  function=commonspace_transforms),
                         name='commonspace_transforms_prep')

        bold_main_wf=init_bold_main_wf(tr=tr, tpattern=tpattern, apply_STC=apply_STC, data_dir_path=data_dir_path, bias_cor_script=bias_cor_script, bias_reg_script=bias_reg_script, coreg_script=coreg_script, commonspace_transform=True, SyN_SDC=True, aCompCor_method=aCompCor_method)

        workflow.connect([
            (anat2nii, joinnode_session, [("nii_file", "file_list")]),
            (joinnode_sub_id, commonspace_reg, [
                ("file_list", "file_list"),
                ]),
            (infosub_id, commonspace_selectfiles, [
                ("subject_id", "sub"),
                ]),
            (infosession, commonspace_selectfiles, [
                ("session", "ses")
                ]),
            (commonspace_reg, commonspace_selectfiles, [
                ("ants_dbm_common", "ants_dbm_common"),
                ("common_to_template_transform", "common_to_template_transform"),
                ("template_to_common_transform", "template_to_common_transform"),
                ]),
            (commonspace_reg, datasink, [
                ("ants_dbm_common", "ants_dbm_common"),
                ("common_to_template_transform", "common_to_template_transform"),
                ("template_to_common_transform", "template_to_common_transform"),
                ]),
            (commonspace_selectfiles, transform_masks, [
                ("template_to_common_transform", "template_to_common_transform"),
                ("ants_dbm_affine", "ants_dbm_affine"),
                ("ants_dbm_inverse_warp", "ants_dbm_inverse_warp"),
                ]),
            (anat2nii, transform_masks, [
                ("nii_file", "reference_image"),
                ]),
            (anat2nii, bold_main_wf, [
                ("nii_file", "inputnode.anat_preproc"),
                ]),
            (transform_masks, bold_main_wf, [
                ("anat_labels", 'inputnode.labels'),
                ("brain_mask", 'inputnode.anat_mask'),
                ("WM_mask", "inputnode.WM_mask"),
                ("CSF_mask", "inputnode.CSF_mask"),
                ]),
            (commonspace_selectfiles, commonspace_transforms_prep, [
                ("common_to_template_transform", "common_to_template_transform"),
                ("ants_dbm_affine", "ants_dbm_affine"),
                ("ants_dbm_warp", "ants_dbm_warp"),
                ]),
            (commonspace_selectfiles, datasink, [
                ("ants_dbm_affine", "ants_dbm_affine"),
                ("ants_dbm_warp", "ants_dbm_warp"),
                ]),
            (commonspace_transforms_prep, bold_main_wf, [
                ("transforms_list", "inputnode.commonspace_transforms_list"),
                ("inverses", "inputnode.commonspace_inverses"),
                ]),
            (bold_main_wf, outputnode, [
                ("outputnode.commonspace_bold", "commonspace_bold"),
                ]),
            (outputnode, datasink, [
                ("commonspace_bold", "commonspace_bold"), #resampled EPI after motion realignment and SDC
                ("commonspace_mask", "commonspace_bold_mask"),
                ("commonspace_WM_mask", "commonspace_bold_WM_mask"),
                ("commonspace_CSF_mask", "commonspace_bold_CSF_mask"),
                ("commonspace_labels", "commonspace_bold_labels"),
                ]),
            ])

    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (anat_selectfiles, anat_preproc_wf, [("anat", "inputnode.anat_file")]),
        (anat_preproc_wf, anat2nii, [("outputnode.preproc_anat", "mnc_file")]),
        (anat2nii, datasink, [("nii_file", 'anat_preproc')]),
        (anat_preproc_wf, outputnode, [("outputnode.preproc_anat", "anat_preproc")]),
        (bold_selectfiles, datasink, [
            ("bold", "input_bold"),
            ]),
        (bold_selectfiles, bold_main_wf, [
            ("bold", "inputnode.bold"),
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
        (outputnode, datasink, [
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

def mnc2nii(mnc_file):
    import os
    cwd = os.getcwd()
    basename=os.path.basename(mnc_file).split('.')[0]
    os.system('mnc2nii %s %s/%s.nii' % (mnc_file,cwd,basename))
    os.system('gzip *.nii')
    return '%s/%s.nii.gz' % (cwd,basename)


def commonspace_reg_function(file_list, output_folder):
    import os
    import pandas as pd
    cwd = os.getcwd()
    csv_path=cwd+'/commonspace_input_files.csv'
    files=[]
    for ses_file in file_list:
        for file in ses_file:
            files.append(file)
    df = pd.DataFrame(data=files)
    df.to_csv(csv_path, header=False, sep=',',index=False)

    model_script_path = os.environ["RABIES"]+ '/rabies/shell_scripts/ants_dbm.sh'
    print('Running commonspace registration.')
    os.system('bash %s %s' % (model_script_path,csv_path))


    template_folder=output_folder+'/ants_dbm_outputs/'
    os.system('mkdir -p %s' % (template_folder,))
    os.system('mv * %s' % (template_folder,))

    #ants dbm outputs
    ants_dbm_common = '/ants_dbm/output/secondlevel/secondlevel_template0.nii.gz'
    common_to_template_transform = '/template_reg/template_reg_Composite.h5'
    template_to_common_transform = '/template_reg/template_reg_InverseComposite.h5'

    return ants_dbm_common, common_to_template_transform, template_to_common_transform


'''
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
    info_csv_path=cwd+'/commonspace_info.csv'
    anat_id=[]
    sub_id=[]
    session=[]
    for sub_file_list in file_list:
        for file in sub_file_list:
            anat_id.append(file)
            sub_id.append(os.path.basename(file).split('_ses-')[0])
            session.append(os.path.basename(file).split('_ses-')[1][0])
    if commonspace_method=='pydpiper':
        df = pd.DataFrame(data={"file": anat_id})
        df.to_csv(csv_path, sep=',',index=False)
    elif commonspace_method=='ants_dbm':
        df = pd.DataFrame(data=anat_id)
        df.to_csv(csv_path, header=False, sep=',',index=False)
        df = pd.DataFrame(data=anat_id)
        df['sub_id']=sub_id
        df['session']=session
        df.to_csv(info_csv_path, header=False, sep=',',index=False)
    return csv_path, info_csv_path

def move_commonspace_transforms(file_list_buffer, commonspace_directory, transform_csv, output_folder):
    import os
    os.system('mkdir -p %s/anat_datasink/commonspace_transforms/' % (output_folder))
    os.system('mv %s %s/anat_datasink/commonspace_transforms/' % (transform_csv, output_folder))
    os.system('mv %s/mbm_atlasReg_processed %s/anat_datasink/commonspace_transforms/' % (commonspace_directory, output_folder))
    return None


def init_anat_init_wf(data_csv, data_dir_path, output_folder, commonspace_method='pydpiper', name='anat_init_wf'):

    workflow = pe.Workflow(name=name)
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_preproc','csv_file', 'info_csv']),
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
                              output_names=['csv_file', 'info_csv'],
                              function=commonspace_prep_func),
                     name='commonspace_prep',
                     joinsource='infosource',
                     joinfield=['file_list'])
    commonspace_prep.inputs.commonspace_method=commonspace_method

    #connect the anat workflow
    workflow.connect([
        (anat_selectfiles, anat_preproc_wf, [("anat", "inputnode.anat_file")]),
        (anat_preproc_wf, datasink, [("outputnode.preproc_anat", "commonspace_inputs")]),
        (anat_preproc_wf, anat2nii, [("outputnode.preproc_anat", "mnc_file")]),
        (anat2nii, datasink, [("nii_file", 'anat_preproc')]),
        (anat_preproc_wf, outputnode, [("outputnode.preproc_anat", "anat_preproc")]),
        (joinnode_buffer, commonspace_prep, [("file_list", "file_list")]),
        (commonspace_prep, outputnode, [("csv_file", "csv_file"),("info_csv", "info_csv")]),
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


def init_main_postcommonspace_wf(data_csv, data_dir_path, output_folder, csv_labels, tr, tpattern, apply_STC=True, bold_preproc_only=False, commonspace_transform=False,compute_WM_CSF_masks=False,
                bias_cor_script='Default', bias_reg_script='Rigid', coreg_script='SyN', aCompCor_method='50%', name='main_wf'):

    workflow = pe.Workflow(name=name)
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_mask', 'anat_labels', 'WM_mask', 'CSF_mask','initial_bold_ref', 'bias_cor_bold', 'confounds_csv', 'itk_bold_to_anat', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv',
                'itk_anat_to_bold','boldref_warped2anat', 'native_corrected_bold', 'corrected_ref_bold', 'bold_brain_mask', 'bold_WM_mask', 'bold_CSF_mask', 'bold_labels', 'commonspace_bold']),
        name='outputnode')

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

    file_info = pe.Node(Function(input_names=['file_path'],
                              output_names=['subject_id', 'session'],
                              function=file_reader),
                     name='file_info')


    # Datasink - creates output folder for important outputs
    datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="anat_datasink"),
                    name="anat_datasink")

    anat2nii = pe.Node(Function(input_names=['mnc_file'],
                              output_names=['nii_file'],
                              function=mnc2nii),
                     name='anat2nii')

    anat_nii = output_folder+'/'+opj('anat_datasink','anat_preproc','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_anat_preproc.nii.gz')
    labels = output_folder+'/'+opj('anat_datasink','anat_labels','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_anat_labels.nii.gz')
    mask = output_folder+'/'+opj('anat_datasink','anat_mask','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_anat_mask.nii.gz')

    commonspace_affine = output_folder+'/'+opj('anat_datasink','commonspace_transforms','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_Affine.mat')
    commonspace_warp = output_folder+'/'+opj('anat_datasink','commonspace_transforms','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_Warp.nii.gz')
    commonspace_anat_template = output_folder+'/'+opj('anat_datasink','commonspace_template','commonspace_template.nii.gz')
    if not compute_WM_CSF_masks:
        WM_mask = output_folder+'/'+opj('anat_datasink','WM_mask','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_WM_mask.nii.gz')
        CSF_mask = output_folder+'/'+opj('anat_datasink','CSF_mask','_subject_id_{subject_id}','_session_{session}','{subject_id}_ses-{session}_CSF_mask.nii.gz')
        commonspace_templates = {'anat_nii':anat_nii,'labels': labels, 'mask': mask, 'WM_mask': WM_mask, 'CSF_mask': CSF_mask, 'commonspace_affine': commonspace_affine, 'commonspace_warp': commonspace_warp, 'commonspace_anat_template':commonspace_anat_template}
    else:
        commonspace_templates = {'anat_nii':anat_nii,'labels': labels, 'mask': mask, 'commonspace_affine': commonspace_affine, 'commonspace_warp': commonspace_warp, 'commonspace_anat_template':commonspace_anat_template}

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
    bold_main_wf=init_bold_main_wf(tr=tr, tpattern=tpattern, apply_STC=apply_STC, data_dir_path=data_dir_path, bias_cor_script=bias_cor_script, bias_reg_script=bias_reg_script, coreg_script=coreg_script, commonspace_transform=commonspace_transform, SyN_SDC=True, aCompCor_method=aCompCor_method)

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


    return workflow
'''
