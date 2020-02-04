import os
from os.path import join as opj

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, afni
from nipype.interfaces.io import SelectFiles

from .hmc import init_bold_hmc_wf
from .utils import init_bold_reference_wf, BIDSDataGraber, prep_bids_iter, convert_to_RAS
from .resampling import init_bold_preproc_trans_wf, init_bold_commonspace_trans_wf
from .stc import init_bold_stc_wf
from .bias_correction import bias_correction_wf
from .registration import init_bold_reg_wf, run_antsRegistration
from .confounds import init_bold_confs_wf
from nipype.interfaces.utility import Function

def init_bold_main_wf(data_dir_path, apply_despiking=False, tr='1.0s', tpattern='altplus', apply_STC=True, detect_dummy=False, slice_mc=False, bias_reg_script='Rigid', coreg_script='SyN',
                        nativespace_resampling='origin', commonspace_resampling='origin', aCompCor_method='50%', name='bold_main_wf'):
    """
    This workflow controls the functional preprocessing stages of the pipeline when both
    functional and anatomical images are provided.

    **Parameters**

        data_dir_path
            Path to the input data directory with proper input folder structure.
        tr
            repetition time for the EPI
        tpattern
            specification for the within TR slice acquisition method. The input is fed to AFNI's 3dTshift
        apply_STC
            whether to apply slice timing correction (STC) or not
        bias_reg_script
            path to registration script that will be applied for bias field correction. The script must
            follow the template structure of registration scripts in shell_scripts/.
            Default is set to 'Rigid' registration.
        coreg_script
            path to registration script for EPI to anat coregistraion. The script must
            follow the template structure of registration scripts in shell_scripts/.
            Default is set to 'SyN' registration.
        nativespace_resampling
            Specified dimensions for the resampling of the corrected EPI in native space.
        commonspace_resampling
            Specified dimensions for the resampling of the corrected EPI in common space.

    **Inputs**

        bold
            Input BOLD series NIfTI file
        anat_preproc
            Preprocessed anatomical image after bias field correction and denoising
        anat_mask
            Brain mask inherited from the common space registration
        WM_mask
            Eroded WM mask inherited from the common space registration
        CSF_mask
            Eroded CSF mask inherited from the common space registration
        labels
            Anatomical labels inherited from the common space registration
        commonspace_transforms_list
            list of transforms to be applied to resample to commonspace
        commonspace_inverses
            Specification for the application of inverse affine transform for
            the provided commonspace transforms

    **Outputs**

        input_bold
            The provided input BOLD file
        bold_ref
            Initial EPI median volume subsequently used as 3D reference EPI volume
        skip_vols
            Initial saturated volumes detected through the generation of the bold_ref
        motcorr_params
            motion parameters file provided from antsMotionCorr
        corrected_EPI
            3D reference EPI volume after bias field correction
        itk_bold_to_anat
            Composite transforms from the EPI space to the anatomical space
        itk_anat_to_bold
            Composite transforms from the anatomical space to the EPI space
        output_warped_bold
            Bias field corrected 3D EPI volume warped to the anatomical space
        resampled_bold
            Original BOLD timeseries resampled through motion realignment and
            susceptibility distortion correction based on registration to the
            anatomical image
        resampled_ref_bold
            3D median EPI volume from the resampled native BOLD timeseries
        confounds_csv
            .csv file with measured confound timecourses, including global signal,
            WM signal, CSF signal, 6 rigid body motion parameters + their first
            temporal derivate + the 12 parameters squared (24 motion parameters),
            and aCompCorr timecourses
        FD_voxelwise
            Voxelwise framewise displacement (FD) measures that can be integrated
            to future confound regression.
            These measures are computed from antsMotionCorrStats.
        pos_voxelwise
            Voxel distancing across time based on rigid body movement parameters,
            which can be integrated for a voxelwise motion regression
            These measures are computed from antsMotionCorrStats.
        FD_csv
            .csv file with global framewise displacement (FD) measures
        EPI_brain_mask
            EPI brain mask for resampled bold
        EPI_WM_mask
            EPI WM mask for resampled bold
        EPI_CSF_mask
            EPI CSF mask for resampled bold
        EPI_labels
            EPI anatomical labels for resampled bold
        commonspace_bold
            Motion and SDC-corrected EPI timeseries resampled into common space
            by applying transforms from the anatomical common space registration
        commonspace_mask
            EPI brain mask for commonspace bold
        commonspace_WM_mask
            EPI WM mask for commonspace bold
        commonspace_CSF_mask
            EPI CSF mask for commonspace bold
        commonspace_labels
            EPI anatomical labels for commonspace bold
    """

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['subject_id', 'bold', 'anat_preproc', 'anat_mask', 'WM_mask', 'CSF_mask', 'vascular_mask', 'labels', 'template_to_common_affine', 'template_to_common_warp','anat_to_template_affine','anat_to_template_warp']),
                      name="inputnode")

    outputnode = pe.Node(niu.IdentityInterface(
                fields=['input_bold', 'bold_ref', 'skip_vols','motcorr_params', 'corrected_EPI', 'output_warped_bold', 'affine_bold2anat', 'warp_bold2anat', 'inverse_warp_bold2anat','resampled_bold', 'resampled_ref_bold', 'EPI_brain_mask', 'EPI_WM_mask', 'EPI_CSF_mask', 'EPI_labels',
                        'confounds_csv', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv', 'commonspace_bold', 'commonspace_mask', 'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_labels']),
                name='outputnode')

    bold_reference_wf = init_bold_reference_wf(detect_dummy=detect_dummy)
    bias_cor_wf = bias_correction_wf(bias_reg_script=bias_reg_script)

    if apply_STC:
        bold_stc_wf = init_bold_stc_wf(tr=tr, tpattern=tpattern)

    # BOLD buffer: an identity used as a pointer to the STC data for further use.
    boldbuffer = pe.Node(niu.IdentityInterface(fields=['bold_file']), name='boldbuffer')

    # HMC on the BOLD
    bold_hmc_wf = init_bold_hmc_wf(slice_mc=slice_mc, name='bold_hmc_wf')

    bold_reg_wf = init_bold_reg_wf(coreg_script=coreg_script)

    def SyN_coreg_transforms_prep(warp_bold2anat,affine_bold2anat):
        return [warp_bold2anat,affine_bold2anat],[0,0] #transforms_list,inverses
    transforms_prep = pe.Node(Function(input_names=['warp_bold2anat','affine_bold2anat'],
                              output_names=['transforms_list','inverses'],
                              function=SyN_coreg_transforms_prep),
                     name='transforms_prep')

    # Apply transforms in 1 shot
    bold_bold_trans_wf = init_bold_preproc_trans_wf(resampling_dim=nativespace_resampling, name='bold_bold_trans_wf')

    bold_confs_wf = init_bold_confs_wf(aCompCor_method=aCompCor_method, name="bold_confs_wf")

    if apply_despiking:
        despike = pe.Node(
            afni.Despike(outputtype='NIFTI_GZ'),
            name='despike')
        workflow.connect([
            (inputnode, despike, [('bold', 'in_file')]),
            (despike, bold_reference_wf, [('out_file', 'inputnode.bold_file')]),
            ])
    else:
        workflow.connect([
            (inputnode, bold_reference_wf, [('bold', 'inputnode.bold_file')]),
            ])


    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
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
            ('vascular_mask', 'inputnode.vascular_mask'),
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
            ('outputnode.affine_bold2anat', 'affine_bold2anat'),
            ('outputnode.warp_bold2anat', 'warp_bold2anat'),
            ('outputnode.inverse_warp_bold2anat', 'inverse_warp_bold2anat'),
            ('outputnode.output_warped_bold', 'output_warped_bold'),
            ]),
        (bold_reg_wf, transforms_prep, [
            ('outputnode.affine_bold2anat', 'affine_bold2anat'),
            ('outputnode.warp_bold2anat', 'warp_bold2anat'),
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

    def commonspace_transforms( template_to_common_warp, template_to_common_affine,anat_to_template_warp, anat_to_template_affine, warp_bold2anat, affine_bold2anat):
        return [template_to_common_warp, template_to_common_affine,anat_to_template_warp, anat_to_template_affine, warp_bold2anat, affine_bold2anat],[0,0,0,0,0,0] #transforms_list,inverses
    commonspace_transforms_prep = pe.Node(Function(input_names=['template_to_common_warp', 'template_to_common_affine','anat_to_template_warp','anat_to_template_affine', 'warp_bold2anat', 'affine_bold2anat'],
                              output_names=['transforms_list','inverses'],
                              function=commonspace_transforms),
                     name='commonspace_transforms_prep')

    bold_commonspace_trans_wf = init_bold_commonspace_trans_wf(resampling_dim=commonspace_resampling, name='bold_commonspace_trans_wf')

    workflow.connect([
        (inputnode, commonspace_transforms_prep, [
            ("template_to_common_affine", "template_to_common_affine"),
            ("template_to_common_warp", "template_to_common_warp"),
            ("anat_to_template_affine", "anat_to_template_affine"),
            ("anat_to_template_warp", "anat_to_template_warp"),
            ]),
        (bold_reg_wf, commonspace_transforms_prep, [
            ('outputnode.affine_bold2anat', 'affine_bold2anat'),
            ('outputnode.warp_bold2anat', 'warp_bold2anat'),
            ]),
        (commonspace_transforms_prep, bold_commonspace_trans_wf, [
            ('transforms_list', 'inputnode.transforms_list'),
            ('inverses', 'inputnode.inverses'),
            ]),
        (boldbuffer, bold_commonspace_trans_wf, [('bold_file', 'inputnode.bold_file')]),
        (bold_hmc_wf, bold_commonspace_trans_wf, [('outputnode.motcorr_params', 'inputnode.motcorr_params')]),
        (inputnode, bold_commonspace_trans_wf, [
            ('bold', 'inputnode.name_source'),
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


def init_EPIonly_bold_main_wf(data_dir_path, data_csv, output_folder, bids_input=False, apply_despiking=False, tr='1.0s', tpattern='altplus', apply_STC=True, detect_dummy=False, slice_mc=False, bias_reg_script='Rigid', coreg_script='SyN', template_reg_script=None,
                        commonspace_resampling='origin', aCompCor_method='50%', name='bold_main_wf'):
    """
    This is an alternative workflow for EPI-only preprocessing, inluding commonspace
    registration based on the generation of a EPI template from the provided sample
    and registration of that template to the provided external anatomical template for
    masking and labeling.

    **Parameters**

        data_dir_path
            Path to the input data directory with proper input folder structure.
        data_csv
            csv file specifying subject id and number of sessions and runs
        output_folder
            path to output folder for the workflow and datasink
        bids_input
            specify if the provided input folder is in a BIDS format to use BIDS reader
        tr
            repetition time for the EPI
        tpattern
            specification for the within TR slice acquisition method. The input is fed to AFNI's 3dTshift
        apply_STC
            whether to apply slice timing correction (STC) or not
        bias_reg_script
            path to registration script that will be applied for bias field correction. The script must
            follow the template structure of registration scripts in shell_scripts/.
            Default is set to 'Rigid' registration.
        coreg_script
            path to registration script for EPI to anat coregistraion. The script must
            follow the template structure of registration scripts in shell_scripts/.
            Default is set to 'SyN' registration.
        commonspace_resampling
            Specified dimensions for the resampling of the corrected EPI in common space.

    **Outputs**

        input_bold
            The provided input BOLD file
        bold_ref
            Initial EPI median volume subsequently used as 3D reference EPI volume
        skip_vols
            Initial saturated volumes detected through the generation of the bold_ref
        motcorr_params
            motion parameters file provided from antsMotionCorr
        corrected_EPI
            3D reference EPI volume after bias field correction
        ants_dbm_affine
            Affine transforms from the EPI subject space to the EPI template
            space
        ants_dbm_warp
            Non-linear transforms from the EPI subject space
            to the EPI template space
        ants_dbm_inverse_warp
            Inverse for the non-linear transforms from the EPI subject space
            to the EPI template space
        ants_dbm_template
            EPI template generated from ants_dbm
        common_to_template_transform
            Inverse composite transforms from the registration of the EPI
            template to the anatomical template.
        template_to_common_transform
            Composite transforms from the registration of the EPI template to
            the anatomical template.
        warped_template
            ants_dbm template registered to the atlas template
        resampled_bold
            Original BOLD timeseries resampled through motion realignment and
            resampling to the anatomical common space, simultaneously correcting
            for susceptibility distortion correction
        resampled_ref_bold
            3D median EPI volume from the resampled BOLD timeseries
        confounds_csv
            .csv file with measured confound timecourses, including global signal,
            WM signal, CSF signal, 6 rigid body motion parameters + their first
            temporal derivate + the 12 parameters squared (24 motion parameters),
            and aCompCorr timecourses
        FD_voxelwise
            Voxelwise framewise displacement (FD) measures that can be integrated
            to future confound regression.
            These measures are computed from antsMotionCorrStats.
        pos_voxelwise
            Voxel distancing across time based on rigid body movement parameters,
            which can be integrated for a voxelwise motion regression
            These measures are computed from antsMotionCorrStats.
        FD_csv
            .csv file with global framewise displacement (FD) measures
        EPI_brain_mask
            EPI brain mask for resampled bold
        EPI_WM_mask
            EPI WM mask for resampled bold
        EPI_CSF_mask
            EPI CSF mask for resampled bold
        EPI_labels
            EPI anatomical labels for resampled bold
    """

    from nipype.interfaces.io import SelectFiles, DataSink

    print("BOLD preproc only!")

    workflow = pe.Workflow(name=name)

    outputnode = pe.Node(niu.IdentityInterface(
                fields=['input_bold', 'bold_ref', 'skip_vols','motcorr_params', 'corrected_EPI', 'ants_dbm_inverse_warp', 'ants_dbm_warp', 'ants_dbm_affine', 'ants_dbm_template', 'template_to_common_affine','template_to_common_warp','template_to_common_inverse_warp','warped_template',
                        'resampled_bold', 'resampled_ref_bold', 'EPI_brain_mask', 'EPI_WM_mask', 'EPI_CSF_mask', 'EPI_labels', 'confounds_csv', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv']),
                name='outputnode')

    if bids_input:
        #with BIDS input data
        from bids.layout import BIDSLayout
        layout = BIDSLayout(data_dir_path)
        subject_list, session_iter, run_iter=prep_bids_iter(layout)
        #set SelectFiles nodes
        bold_selectfiles = pe.Node(BIDSDataGraber(bids_dir=data_dir_path, datatype='func'), name='bold_selectfiles')
    else:
        import pandas as pd
        data_df=pd.read_csv(data_csv, sep=',', dtype=str)
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

        bold_file = opj('sub-{subject_id}', 'ses-{session}', 'func', 'sub-{subject_id}_ses-{session}_run-{run}_bold.nii.gz')
        bold_selectfiles = pe.Node(SelectFiles({'out_file': bold_file},
                                       base_directory=data_dir_path),
                           name="bold_selectfiles")


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


    # Datasink - creates output folder for important outputs
    bold_datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="bold_datasink"),
                    name="bold_datasink")

    commonspace_datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="commonspace_datasink"),
                    name="commonspace_datasink")

    transforms_datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="transforms_datasink"),
                    name="transforms_datasink")

    confounds_datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="confounds_datasink"),
                    name="confounds_datasink")

    #node to conver input image to consistent RAS orientation
    convert_to_RAS_node = pe.Node(Function(input_names=['img_file'],
                              output_names=['RAS_file'],
                              function=convert_to_RAS),
                     name='convert_to_RAS')

    bold_reference_wf = init_bold_reference_wf(detect_dummy=detect_dummy)
    bias_cor_wf = bias_correction_wf(bias_reg_script=bias_reg_script)
    bias_cor_wf.inputs.inputnode.anat=os.environ["template_anat"]
    bias_cor_wf.inputs.inputnode.anat_mask=os.environ["template_mask"]

    if apply_STC:
        bold_stc_wf = init_bold_stc_wf(tr=tr, tpattern=tpattern)

    # BOLD buffer: an identity used as a pointer to the STC data for further use.
    boldbuffer = pe.Node(niu.IdentityInterface(fields=['bold_file']), name='boldbuffer')

    # HMC on the BOLD
    bold_hmc_wf = init_bold_hmc_wf(slice_mc=slice_mc, name='bold_hmc_wf')

    # Apply transforms in 1 shot
    bold_bold_trans_wf = init_bold_preproc_trans_wf(resampling_dim=commonspace_resampling, name='bold_bold_trans_wf')
    #the EPIs are resampled to the template common space to correct for susceptibility distortion based on non-linear registration
    bold_bold_trans_wf.inputs.inputnode.ref_file = os.environ["template_anat"]


    def to_commonspace_transforms_prep(template_to_common_warp, template_to_common_affine, ants_dbm_warp, ants_dbm_affine):
        #simply list transforms in the proper order
        return [template_to_common_warp, template_to_common_affine, ants_dbm_warp, ants_dbm_affine],[0,0,0,0] #transforms_list,inverses

    transforms_prep = pe.Node(Function(input_names=['template_to_common_warp', 'template_to_common_affine', 'ants_dbm_warp', 'ants_dbm_affine'],
                              output_names=['transforms_list','inverses'],
                              function=to_commonspace_transforms_prep),
                     name='transforms_prep')


    bold_confs_wf = init_bold_confs_wf(aCompCor_method=aCompCor_method, name="bold_confs_wf")
    #give the template masks and labels, which will be assign to every subject scan after corrections, since
    #all scans will be in common space after SDC
    bold_confs_wf.inputs.inputnode.t1_mask = os.environ["template_mask"]
    bold_confs_wf.inputs.inputnode.WM_mask = os.environ["WM_mask"]
    bold_confs_wf.inputs.inputnode.CSF_mask = os.environ["CSF_mask"]
    bold_confs_wf.inputs.inputnode.vascular_mask = os.environ["vascular_mask"]
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

    #calculate the number of bold scans that will be registered
    num_bold=0
    for sub in subject_list:
        num_bold+=(len(session_iter[sub])*len(run_iter[sub]))
    if int(os.environ["local_threads"])<num_bold:
        num_bold=int(os.environ["local_threads"])

    commonspace_reg = pe.Node(Function(input_names=['file_list', 'output_folder'],
                              output_names=['ants_dbm_template'],
                              function=commonspace_reg_function),
                     name='commonspace_reg', n_procs=num_bold, mem_gb=1*num_bold)
    commonspace_reg.inputs.output_folder = output_folder+'/commonspace_datasink/'

    #execute the registration of the generate anatomical template with the provided atlas for labeling and masking
    template_reg = pe.Node(Function(input_names=['reg_script', 'moving_image', 'fixed_image', 'anat_mask'],
                              output_names=['affine', 'warp', 'inverse_warp', 'warped_image'],
                              function=run_antsRegistration),
                     name='template_reg', mem_gb=3)
    template_reg.plugin_args = {'qsub_args': '-pe smp %s' % (str(3*int(os.environ["min_proc"]))), 'overwrite': True}
    template_reg.inputs.fixed_image = os.environ["template_anat"]
    template_reg.inputs.anat_mask = os.environ["template_mask"]
    template_reg.inputs.reg_script = template_reg_script

    #setting SelectFiles for the commonspace registration
    ants_dbm_inverse_warp = output_folder+'/'+opj('commonspace_datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_sub-{sub}_ses-{ses}_run-{run}_bias_cor*1InverseWarp.nii.gz')
    ants_dbm_warp = output_folder+'/'+opj('commonspace_datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_sub-{sub}_ses-{ses}_run-{run}_bias_cor*1Warp.nii.gz')
    ants_dbm_affine = output_folder+'/'+opj('commonspace_datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_sub-{sub}_ses-{ses}_run-{run}_bias_cor*0GenericAffine.mat')
    template_to_common_affine = '/'+opj('{template_to_common_affine}')
    template_to_common_warp = '/'+opj('{template_to_common_warp}')
    template_to_common_inverse_warp = '/'+opj('{template_to_common_inverse_warp}')

    commonspace_templates = {'ants_dbm_inverse_warp':ants_dbm_inverse_warp,'ants_dbm_warp': ants_dbm_warp, 'ants_dbm_affine': ants_dbm_affine, 'template_to_common_affine':template_to_common_affine,'template_to_common_warp': template_to_common_warp, 'template_to_common_inverse_warp':template_to_common_inverse_warp}

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
        (commonspace_reg, template_reg, [
            ("ants_dbm_template", "moving_image"),
            ]),
        (template_reg, commonspace_selectfiles, [
            ("affine", "template_to_common_affine"),
            ("warp", "template_to_common_warp"),
            ("inverse_warp", "template_to_common_inverse_warp"),
            ]),
        ])

    #convert file to RAS orientation before any operation in the workflow
    workflow.connect([
        (bold_selectfiles, convert_to_RAS_node, [('out_file', 'img_file')]),
        ])
    if apply_despiking:
        despike = pe.Node(
            afni.Despike(outputtype='NIFTI_GZ'),
            name='despike')
        workflow.connect([
            (convert_to_RAS_node, despike, [('RAS_file', 'in_file')]),
            (despike, bold_reference_wf, [('out_file', 'inputnode.bold_file')]),
            ])
    else:
        workflow.connect([
            (convert_to_RAS_node, bold_reference_wf, [('RAS_file', 'inputnode.bold_file')]),
            ])


    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (convert_to_RAS_node, bold_bold_trans_wf, [('RAS_file', 'inputnode.name_source')]),
        (bold_reference_wf, bias_cor_wf, [
            ('outputnode.ref_image', 'inputnode.ref_EPI')
            ]),
        (bold_reference_wf, bold_hmc_wf, [
            ('outputnode.ref_image', 'inputnode.ref_image'),
            ('outputnode.bold_file', 'inputnode.bold_file')]),
        (boldbuffer, bold_bold_trans_wf, [('bold_file', 'inputnode.bold_file')]),
        (bold_hmc_wf, bold_bold_trans_wf, [('outputnode.motcorr_params', 'inputnode.motcorr_params')]),
        (commonspace_selectfiles, transforms_prep, [
            ("template_to_common_warp", "template_to_common_warp"),
            ("template_to_common_affine", "template_to_common_affine"),
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
        (template_reg, outputnode, [
            ("affine", "template_to_common_affine"),
            ("warp", "template_to_common_warp"),
            ("inverse_warp", "template_to_common_inverse_warp"),
            ("warped_image", "warped_template"),
            ]),
        (commonspace_selectfiles, outputnode, [
            ('ants_dbm_warp', 'ants_dbm_warp'),
            ('ants_dbm_inverse_warp', 'ants_dbm_inverse_warp'),
            ('ants_dbm_affine', 'ants_dbm_affine'),
            ]),
        (commonspace_reg, outputnode, [
            ("ants_dbm_template", "ants_dbm_template"),
            ]),
        (commonspace_reg, commonspace_datasink, [
            ("ants_dbm_template", "ants_dbm_template"),
            ]),
        (outputnode, bold_datasink, [
            ("bold_ref","initial_bold_ref"), #inspect initial bold ref
            ("corrected_EPI","bias_cor_bold"), #inspect bias correction
            ("EPI_brain_mask","bold_brain_mask"), #get the EPI labels
            ("EPI_WM_mask","bold_WM_mask"), #get the EPI labels
            ("EPI_CSF_mask","bold_CSF_mask"), #get the EPI labels
            ("EPI_labels","bold_labels"), #get the EPI labels
            ("resampled_bold", "corrected_bold"), #resampled EPI after motion realignment and SDC
            ("resampled_ref_bold", "corrected_bold_ref"), #resampled EPI after motion realignment and SDC
            ]),
        (outputnode, confounds_datasink, [
            ("confounds_csv", "confounds_csv"), #confounds file
            ("FD_voxelwise", "FD_voxelwise"),
            ("pos_voxelwise", "pos_voxelwise"),
            ("FD_csv", "FD_csv"),
            ]),
        (outputnode, transforms_datasink, [
            ('ants_dbm_warp', 'ants_dbm_warp'),
            ('ants_dbm_inverse_warp', 'ants_dbm_inverse_warp'),
            ('ants_dbm_affine', 'ants_dbm_affine'),
            ]),
        (outputnode, transforms_datasink, [
            ("template_to_common_affine", "template_to_common_affine"),
            ("template_to_common_warp", "template_to_common_warp"),
            ("template_to_common_inverse_warp", "template_to_common_inverse_warp"),
            ]),
        (outputnode, commonspace_datasink, [
            ("warped_template", "warped_template"),
            ]),
        ])

    return workflow


def commonspace_reg_function(file_list, output_folder):
    import os
    import numpy as np
    import pandas as pd
    #create a csv file of the input image list
    cwd = os.getcwd()
    csv_path=cwd+'/commonspace_input_files.csv'

    import itertools
    merged = list(itertools.chain.from_iterable(file_list))
    df = pd.DataFrame(data=merged)
    df.to_csv(csv_path, header=False, sep=',',index=False)

    model_script_path = os.environ["RABIES"]+ '/rabies/shell_scripts/ants_dbm.sh'

    print('Running commonspace registration.')
    command='bash %s %s' % (model_script_path,csv_path)
    if os.system(command) != 0:
        raise ValueError('Error in running commonspace registration with: '+command)

    #copy all outputs to provided output folder to prevent deletion of the files after the node has run
    template_folder=output_folder+'/ants_dbm_outputs/'
    command='rm -r %s' % (template_folder,)
    if os.system(command) != 0:
        raise ValueError('Error in '+command)
    command='mkdir %s' % (template_folder,)
    if os.system(command) != 0:
        raise ValueError('Error in '+command)
    command='cp -r * %s' % (template_folder,)
    if os.system(command) != 0:
        raise ValueError('Error in '+command)

    ###verify that all outputs are present
    #ants dbm outputs
    ants_dbm_template = template_folder+'/ants_dbm/output/secondlevel/secondlevel_template0.nii.gz'
    if not os.path.isfile(ants_dbm_template):
        raise ValueError(ants_dbm_template+" doesn't exists.")

    i=0
    for file in merged:
        file=str(file)
        filename_template=os.path.basename(file).split('.')[0]
        anat_to_template_inverse_warp = '%s/ants_dbm/output/secondlevel/secondlevel_%s%s1InverseWarp.nii.gz' % (template_folder,filename_template,str(i),)
        if not os.path.isfile(anat_to_template_inverse_warp):
            raise ValueError(anat_to_template_inverse_warp+" file doesn't exists.")
        anat_to_template_warp = '%s/ants_dbm/output/secondlevel/secondlevel_%s%s1Warp.nii.gz' % (template_folder,filename_template,str(i),)
        if not os.path.isfile(anat_to_template_warp):
            raise ValueError(anat_to_template_warp+" file doesn't exists.")
        anat_to_template_affine = '%s/ants_dbm/output/secondlevel/secondlevel_%s%s0GenericAffine.mat' % (template_folder,filename_template,str(i),)
        if not os.path.isfile(anat_to_template_affine):
            raise ValueError(anat_to_template_affine+" file doesn't exists.")
        i+=1

    return ants_dbm_template
