import os
from os.path import join as opj

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import SelectFiles

from .hmc import init_bold_hmc_wf
from .utils import init_bold_reference_wf
from .resampling import init_bold_preproc_trans_wf, init_bold_commonspace_trans_wf
from .stc import init_bold_stc_wf
from .sdc import init_sdc_wf
from .bias_correction import bias_correction_wf
from .registration import init_bold_reg_wf
from .confounds import init_bold_confs_wf

from nipype.interfaces.utility import Function

def init_bold_main_wf(data_dir_path, tr='1.0s', tpattern='altplus', apply_STC=True, data_csv=None, bias_cor_script='Default', bias_reg_script='Rigid', coreg_script='SyN', SyN_SDC=True, commonspace_transform=False,
                        aCompCor_method='50%', name='bold_main_wf'):

    """
    This workflow controls the functional preprocessing stages of the pipeline.


    **Parameters**

        use_syn : bool
            Use ANTs SyN-based susceptibility distortion correction (SDC) during
            EPI to anat coregistration.

    **Inputs**

        bold_file
            BOLD series NIfTI file
        reversed_bold_file
            EPI acquired with the reversed phase encoding direction to apply topup distortion correction
        anat_preproc
            Bias-corrected structural template image
        anat_mask
            Mask of the preprocessed anat
        anat_labels
            Labels derived from atlas registration.


    **Outputs**

        native_bold
            Preprocessed BOLD series, resampled to BOLD native space
        bold_anat
            BOLD series, resampled to anatw space
        bold_mask_anat
            BOLD series mask in anatw space
        bold_template
            BOLD series, resampled to template space
        bold_mask_mni
            BOLD series mask in template space
        confounds
            TSV of confounds

    """

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['subject_id', 'bold', 'anat_preproc', 'anat_mask', 'WM_mask', 'CSF_mask', 'labels', 'commonspace_transforms_list', 'commonspace_inverses']),
                      name="inputnode")

    outputnode = pe.Node(niu.IdentityInterface(
                fields=['input_bold', 'bold_ref', 'skip_vols','motcorr_params', 'corrected_EPI', 'output_warped_bold', 'itk_bold_to_anat', 'itk_anat_to_bold', 'commonspace_bold', 'commonspace_mask', 'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_labels',
                        'resampled_bold', 'resampled_ref_bold', 'EPI_brain_mask', 'EPI_WM_mask', 'EPI_CSF_mask', 'EPI_labels', 'confounds_csv', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv']),
                name='outputnode')

    bold_reference_wf = init_bold_reference_wf()
    bias_cor_wf = bias_correction_wf(bias_cor_script=bias_cor_script, bias_reg_script=bias_reg_script)

    if apply_STC:
        bold_stc_wf = init_bold_stc_wf(tr=tr, tpattern=tpattern)

    # BOLD buffer: an identity used as a pointer to the STC data for further use.
    boldbuffer = pe.Node(niu.IdentityInterface(fields=['bold_file']), name='boldbuffer')

    # HMC on the BOLD
    bold_hmc_wf = init_bold_hmc_wf(name='bold_hmc_wf')

    bold_reg_wf = init_bold_reg_wf(coreg_script=coreg_script)

    def SyN_coreg_transforms_prep(itk_bold_to_anat):
        return [itk_bold_to_anat],[0] #transforms_list,inverses
    transforms_prep = pe.Node(Function(input_names=['itk_bold_to_anat'],
                              output_names=['transforms_list','inverses'],
                              function=SyN_coreg_transforms_prep),
                     name='transforms_prep')

    # Apply transforms in 1 shot
    bold_bold_trans_wf = init_bold_preproc_trans_wf(
        name='bold_bold_trans_wf'
    )

    bold_confs_wf = init_bold_confs_wf(SyN_SDC=SyN_SDC, aCompCor_method=aCompCor_method, name="bold_confs_wf")


    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (inputnode, bold_reference_wf, [('bold', 'inputnode.bold_file')]),
        (inputnode, bias_cor_wf, [
            ('anat_preproc', 'inputnode.anat'),
            ('anat_mask', 'inputnode.anat_mask'),
            ]),
        (inputnode, bold_reg_wf, [
            ('anat_preproc', 'inputnode.anat_preproc'),
            ('anat_mask', 'inputnode.anat_mask')]),
        (inputnode, bold_bold_trans_wf, [('bold', 'inputnode.name_source')]),
        (inputnode, bold_confs_wf, [('anat_mask', 'inputnode.t1_mask'),
            ('WM_mask', 'inputnode.WM_mask'),
            ('CSF_mask', 'inputnode.CSF_mask'),
            ('labels', 'inputnode.t1_labels'),
            ]),
        (bold_reference_wf, bias_cor_wf, [
            ('outputnode.ref_image', 'inputnode.ref_EPI')]),
        (bold_reference_wf, bold_hmc_wf, [
            ('outputnode.ref_image', 'inputnode.ref_image'),
            ('outputnode.bold_file', 'inputnode.bold_file')]),
        (bold_hmc_wf, outputnode, [
            ('outputnode.motcorr_params', 'motcorr_params')]),
        (bold_reference_wf, outputnode, [
            ('outputnode.ref_image', 'bold_ref')]),
        (bias_cor_wf, bold_reg_wf, [
              ('outputnode.corrected_EPI', 'inputnode.ref_bold_brain')]),
        (bias_cor_wf, outputnode, [
              ('outputnode.corrected_EPI', 'corrected_EPI')]),
        (bold_reg_wf, outputnode, [
            ('outputnode.itk_bold_to_anat', 'itk_bold_to_anat'),
            ('outputnode.itk_anat_to_bold', 'itk_anat_to_bold'),
            ('outputnode.output_warped_bold', 'output_warped_bold'),
            ]),
        (bold_reg_wf, transforms_prep, [
            ('outputnode.itk_bold_to_anat', 'itk_bold_to_anat'),
            ]),
        (transforms_prep, bold_bold_trans_wf, [
            ('transforms_list', 'inputnode.transforms_list'),
            ('inverses', 'inputnode.inverses'),
            ]),
        (bold_reg_wf, bold_bold_trans_wf, [
            ('outputnode.output_warped_bold', 'inputnode.ref_file')]),
        (boldbuffer, bold_bold_trans_wf, [('bold_file', 'inputnode.bold_file')]),
        (bold_hmc_wf, bold_bold_trans_wf, [('outputnode.motcorr_params', 'inputnode.motcorr_params')]),
        (bold_bold_trans_wf, outputnode, [
            ('outputnode.bold_ref', 'resampled_ref_bold'),
            ('outputnode.bold', 'resampled_bold'),
            ]),
        (bold_bold_trans_wf, bold_confs_wf, [('outputnode.bold', 'inputnode.bold'),
            ('outputnode.bold_ref', 'inputnode.ref_bold'),
            ]),
        (bold_hmc_wf, bold_confs_wf, [('outputnode.motcorr_params', 'inputnode.movpar_file'),
            ]),
        (bold_confs_wf, outputnode, [
            ('outputnode.brain_mask', 'EPI_brain_mask'),
            ('outputnode.WM_mask', 'EPI_WM_mask'),
            ('outputnode.CSF_mask', 'EPI_CSF_mask'),
            ('outputnode.EPI_labels', 'EPI_labels'),
            ('outputnode.confounds_csv', 'confounds_csv'),
            ('outputnode.FD_csv', 'FD_csv'),
            ('outputnode.FD_voxelwise', 'FD_voxelwise'),
            ('outputnode.pos_voxelwise', 'pos_voxelwise'),
            ]),
        ])

    if apply_STC:
        workflow.connect([
            (bold_reference_wf, bold_stc_wf, [
                ('outputnode.skip_vols', 'inputnode.skip_vols'),
                ('outputnode.bold_file', 'inputnode.bold_file')]),
            (bold_stc_wf, boldbuffer, [('outputnode.stc_file', 'bold_file')]),
            ])
    else:
        workflow.connect([
            (bold_reference_wf, boldbuffer, [
                ('outputnode.bold_file', 'bold_file')]),
            ])

    if commonspace_transform:
        bold_commonspace_trans_wf = init_bold_commonspace_trans_wf(
            name='bold_commonspace_trans_wf'
        )

        workflow.connect([
            (inputnode, bold_commonspace_trans_wf, [
                ('bold', 'inputnode.name_source'),
                ('commonspace_transforms_list', 'inputnode.transforms_list'),
                ('commonspace_inverses', 'inputnode.inverses'),
                ]),
            (bold_commonspace_trans_wf, outputnode, [
                ('outputnode.bold', 'commonspace_bold'),
                ('outputnode.brain_mask', 'commonspace_mask'),
                ('outputnode.WM_mask', 'commonspace_WM_mask'),
                ('outputnode.CSF_mask', 'commonspace_CSF_mask'),
                ('outputnode.labels', 'commonspace_labels'),
                ]),
        ])

    return workflow


#alternative workflow for EPI-only commonspace
def init_EPIonly_bold_main_wf(data_dir_path, data_csv, output_folder, tr='1.0s', tpattern='altplus', apply_STC=True, bias_cor_script='Default', bias_reg_script='Rigid', coreg_script='SyN', SyN_SDC=True,
                        aCompCor_method='50%', name='bold_main_wf'):
    from nipype.interfaces.io import SelectFiles, DataSink

    print("BOLD preproc only!")

    workflow = pe.Workflow(name=name)

    outputnode = pe.Node(niu.IdentityInterface(
                fields=['input_bold', 'bold_ref', 'skip_vols','motcorr_params', 'corrected_EPI', 'ants_dbm_inverse_warp', 'ants_dbm_warp', 'ants_dbm_affine', 'ants_dbm_common_anat', 'common_to_template_transform', 'template_to_common_transform',
                        'resampled_bold', 'resampled_ref_bold', 'EPI_brain_mask', 'EPI_WM_mask', 'EPI_CSF_mask', 'EPI_labels', 'confounds_csv', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv']),
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

    bold_file = opj('{subject_id}', 'ses-{session}', 'bold', '{subject_id}_ses-{session}_run-{run}_bold.nii.gz')
    bold_selectfiles = pe.Node(SelectFiles({'bold': bold_file},
                                   base_directory=data_dir_path),
                       name="bold_selectfiles")


    # Datasink - creates output folder for important outputs
    datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="datasink"),
                    name="datasink")

    bold_reference_wf = init_bold_reference_wf()
    bias_cor_wf = bias_correction_wf(bias_cor_script=bias_cor_script, bias_reg_script=bias_reg_script)
    bias_cor_wf.inputs.inputnode.anat=os.environ["template_anat"]
    bias_cor_wf.inputs.inputnode.anat_mask=os.environ["template_mask"]

    if apply_STC:
        bold_stc_wf = init_bold_stc_wf(tr=tr, tpattern=tpattern)

    # BOLD buffer: an identity used as a pointer to the STC data for further use.
    boldbuffer = pe.Node(niu.IdentityInterface(fields=['bold_file']), name='boldbuffer')

    # HMC on the BOLD
    bold_hmc_wf = init_bold_hmc_wf(name='bold_hmc_wf')

    # Apply transforms in 1 shot
    bold_bold_trans_wf = init_bold_preproc_trans_wf(
        name='bold_bold_trans_wf'
    )
    #the EPIs are resampled to the template common space to correct for susceptibility distortion based on non-linear registration
    bold_bold_trans_wf.inputs.inputnode.ref_file = os.environ["template_anat"]


    def to_commonspace_transforms_prep(common_to_template_transform, ants_dbm_warp, ants_dbm_affine):
        #simply list transforms in the proper order
        return [common_to_template_transform, ants_dbm_warp, ants_dbm_affine],[0,0,0] #transforms_list,inverses

    transforms_prep = pe.Node(Function(input_names=['common_to_template_transform', 'ants_dbm_warp', 'ants_dbm_affine'],
                              output_names=['transforms_list','inverses'],
                              function=to_commonspace_transforms_prep),
                     name='transforms_prep')


    bold_confs_wf = init_bold_confs_wf(SyN_SDC=SyN_SDC, aCompCor_method=aCompCor_method, name="bold_confs_wf")
    #give the template masks and labels, which will be assign to every subject scan after corrections, since
    #all scans will be in common space after SDC
    bold_confs_wf.inputs.inputnode.t1_mask = os.environ["template_mask"]
    bold_confs_wf.inputs.inputnode.WM_mask = os.environ["WM_mask"]
    bold_confs_wf.inputs.inputnode.CSF_mask = os.environ["CSF_mask"]
    bold_confs_wf.inputs.inputnode.t1_labels = os.environ["atlas_labels"]

    #####setting up commonspace registration within the workflow
    joinnode_run = pe.JoinNode(niu.IdentityInterface(fields=['file_list']),
                     name='joinnode_run',
                     joinsource='inforun',
                     joinfield=['file_list'])

    joinnode_session = pe.JoinNode(niu.IdentityInterface(fields=['file_list']),
                     name='joinnode_session',
                     joinsource='infosession',
                     joinfield=['file_list'])

    joinnode_sub_id = pe.JoinNode(niu.IdentityInterface(fields=['file_list']),
                     name='joinnode_sub_id',
                     joinsource='infosub_id',
                     joinfield=['file_list'])

    commonspace_reg = pe.Node(Function(input_names=['file_list', 'output_folder'],
                              output_names=['ants_dbm_common', 'common_to_template_transform', 'template_to_common_transform'],
                              function=commonspace_reg_function),
                     name='commonspace_reg')
    commonspace_reg.inputs.output_folder = output_folder+'/datasink/'

    #setting SelectFiles for the commonspace registration
    ants_dbm_inverse_warp = output_folder+'/'+opj('datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_{sub}_ses-{ses}_run-{run}_bias_cor*1InverseWarp.nii.gz')
    ants_dbm_warp = output_folder+'/'+opj('datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_{sub}_ses-{ses}_run-{run}_bias_cor*1Warp.nii.gz')
    ants_dbm_affine = output_folder+'/'+opj('datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_{sub}_ses-{ses}_run-{run}_bias_cor*0GenericAffine.mat')
    ants_dbm_common_anat = output_folder+'/'+opj('datasink','ants_dbm_outputs','{ants_dbm_common}')
    common_to_template_transform = output_folder+'/'+opj('datasink','ants_dbm_outputs','{common_to_template_transform}')
    template_to_common_transform = output_folder+'/'+opj('datasink','ants_dbm_outputs','{template_to_common_transform}')

    commonspace_templates = {'ants_dbm_inverse_warp':ants_dbm_inverse_warp,'ants_dbm_warp': ants_dbm_warp, 'ants_dbm_affine': ants_dbm_affine, 'ants_dbm_common_anat': ants_dbm_common_anat, 'common_to_template_transform': common_to_template_transform, 'template_to_common_transform':template_to_common_transform}

    commonspace_selectfiles = pe.Node(SelectFiles(commonspace_templates),
                   name="commonspace_selectfiles")


    #connect iterables, joinnodes and SelectFiles
    workflow.connect([
        (infosub_id, infosession, [
            ("subject_id", "subject_id"),
            ]),
        (infosub_id, inforun, [
            ("subject_id", "subject_id"),
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
        (bias_cor_wf, joinnode_run, [
            ("outputnode.corrected_EPI", "file_list"),
            ]),
        (joinnode_run, joinnode_session, [
            ("file_list", "file_list"),
            ]),
        (joinnode_session, joinnode_sub_id, [
            ("file_list", "file_list"),
            ]),
        (joinnode_sub_id, commonspace_reg, [
            ("file_list", "file_list"),
            ]),
        (infosub_id, commonspace_selectfiles, [
            ("subject_id", "sub"),
            ]),
        (infosession, commonspace_selectfiles, [
            ("session", "ses")
            ]),
        (inforun, commonspace_selectfiles, [
            ("run", "run")
            ]),
        (commonspace_reg, commonspace_selectfiles, [
            ("ants_dbm_common", "ants_dbm_common"),
            ("common_to_template_transform", "common_to_template_transform"),
            ("template_to_common_transform", "template_to_common_transform"),
            ]),
        ])

    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (bold_selectfiles, bold_reference_wf, [('bold', 'inputnode.bold_file')]),
        (bold_selectfiles, bold_bold_trans_wf, [('bold', 'inputnode.name_source')]),
        (bold_reference_wf, bias_cor_wf, [
            ('outputnode.ref_image', 'inputnode.ref_EPI')
            ]),
        (bold_reference_wf, bold_hmc_wf, [
            ('outputnode.ref_image', 'inputnode.ref_image'),
            ('outputnode.bold_file', 'inputnode.bold_file')]),
        (boldbuffer, bold_bold_trans_wf, [('bold_file', 'inputnode.bold_file')]),
        (bold_hmc_wf, bold_bold_trans_wf, [('outputnode.motcorr_params', 'inputnode.motcorr_params')]),
        (commonspace_selectfiles, transforms_prep, [
            ('common_to_template_transform', 'common_to_template_transform'),
            ('ants_dbm_warp', 'ants_dbm_warp'),
            ('ants_dbm_affine', 'ants_dbm_affine'),
            ]),
        (transforms_prep, bold_bold_trans_wf, [
            ('transforms_list', 'inputnode.transforms_list'),
            ('inverses', 'inputnode.inverses'),
            ]),
        (bold_bold_trans_wf, bold_confs_wf, [('outputnode.bold', 'inputnode.bold'),
            ('outputnode.bold_ref', 'inputnode.ref_bold'),
            ]),
        (bold_hmc_wf, bold_confs_wf, [('outputnode.motcorr_params', 'inputnode.movpar_file'),
            ]),
        ])

    if apply_STC:
        workflow.connect([
            (bold_reference_wf, bold_stc_wf, [
                ('outputnode.skip_vols', 'inputnode.skip_vols'),
                ('outputnode.bold_file', 'inputnode.bold_file')]),
            (bold_stc_wf, boldbuffer, [('outputnode.stc_file', 'bold_file')]),
            ])
    else:
        workflow.connect([
            (bold_reference_wf, boldbuffer, [
                ('outputnode.bold_file', 'bold_file')]),
            ])

    # set outputs
    workflow.connect([
        (bold_hmc_wf, outputnode, [
            ('outputnode.motcorr_params', 'motcorr_params')]),
        (bold_reference_wf, outputnode, [
            ('outputnode.ref_image', 'bold_ref')]),
        (bias_cor_wf, outputnode, [
              ('outputnode.corrected_EPI', 'corrected_EPI')]),
        (bold_bold_trans_wf, outputnode, [
            ('outputnode.bold_ref', 'resampled_ref_bold'),
            ('outputnode.bold', 'resampled_bold'),
            ]),
        (bold_confs_wf, outputnode, [
            ('outputnode.brain_mask', 'EPI_brain_mask'),
            ('outputnode.WM_mask', 'EPI_WM_mask'),
            ('outputnode.CSF_mask', 'EPI_CSF_mask'),
            ('outputnode.EPI_labels', 'EPI_labels'),
            ('outputnode.confounds_csv', 'confounds_csv'),
            ('outputnode.FD_csv', 'FD_csv'),
            ('outputnode.FD_voxelwise', 'FD_voxelwise'),
            ('outputnode.pos_voxelwise', 'pos_voxelwise'),
            ]),
        (commonspace_selectfiles, outputnode, [
            ('common_to_template_transform', 'common_to_template_transform'),
            ('template_to_common_transform', 'template_to_common_transform'),
            ('ants_dbm_warp', 'ants_dbm_warp'),
            ('ants_dbm_inverse_warp', 'ants_dbm_inverse_warp'),
            ('ants_dbm_affine', 'ants_dbm_affine'),
            ('ants_dbm_common_anat', 'ants_dbm_common_anat'),
            ]),
        ])

    return workflow

def commonspace_reg_function(file_list, output_folder):
    import os
    import pandas as pd
    cwd = os.getcwd()
    csv_path=cwd+'/commonspace_input_files.csv'
    files=[]
    for sub_file_list in file_list:
        for ses_file in sub_file_list:
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
