import os
from os.path import join as opj
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from .preprocess_anat_pkg.anat_preproc import init_anat_preproc_wf
from .preprocess_anat_pkg.anat_mask_prep import init_anat_mask_prep_wf
from .preprocess_bold_pkg.bold_main_wf import init_bold_main_wf, commonspace_reg_function
from nipype.interfaces.io import SelectFiles, DataSink

from nipype.interfaces.utility import Function

def init_unified_main_wf(data_dir_path, data_csv, output_folder, tr, tpattern, apply_STC=True, commonspace_method='pydpiper',
                bias_reg_script='Rigid', coreg_script='SyN', isotropic_resampling=False, upsampling=1.0, name='main_wf'):
    '''
    This workflow includes complete anatomical and BOLD preprocessing within a single workflow.

    **Parameters**

        data_dir_path
            Path to the input data directory with proper input folder structure.
        data_csv
            csv file specifying subject id and number of sessions and runs
        output_folder
            path to output folder for the workflow and datasink
        tr
            repetition time for the EPI
        tpattern
            specification for the within TR slice acquisition method. The input is fed to AFNI's 3dTshift
        apply_STC
            whether to apply slice timing correction (STC) or not
        commonspace_method
            specified method for common space registration between 'ants_dbm' and 'pydpiper'
        bias_reg_script
            path to registration script that will be applied for bias field correction. The script must
            follow the template structure of registration scripts in shell_scripts/.
            Default is set to 'Rigid' registration.
        coreg_script
            path to registration script for EPI to anat coregistraion. The script must
            follow the template structure of registration scripts in shell_scripts/.
            Default is set to 'SyN' registration.

    **Outputs**
        fields=['anat_preproc','anat_mask', 'anat_labels', 'WM_mask', 'CSF_mask','initial_bold_ref', 'bias_cor_bold', 'confounds_csv', 'itk_bold_to_anat', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv',
                'itk_anat_to_bold','boldref_warped2anat', 'native_corrected_bold', 'corrected_ref_bold', 'bold_brain_mask', 'bold_WM_mask', 'bold_CSF_mask', 'bold_labels', 'commonspace_bold', 'commonspace_mask', 'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_labels']),

    '''

    workflow = pe.Workflow(name=name)

    #set output node
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_preproc','anat_mask', 'anat_labels', 'WM_mask', 'CSF_mask','initial_bold_ref', 'bias_cor_bold', 'confounds_csv', 'itk_bold_to_anat', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv',
                'itk_anat_to_bold','boldref_warped2anat', 'native_corrected_bold', 'corrected_ref_bold', 'bold_brain_mask', 'bold_WM_mask', 'bold_CSF_mask', 'bold_labels', 'commonspace_bold', 'commonspace_mask', 'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_labels']),
        name='outputnode')


    #read the data_info csv
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

    # Datasink - creates output folder for important outputs
    datasink = pe.Node(DataSink(base_directory=output_folder,
                             container="datasink"),
                    name="datasink")

    #set SelectFiles nodes
    anat_file = opj('{subject_id}', 'ses-{session}', 'anat', '{subject_id}_ses-{session}_anat.nii.gz')
    anat_selectfiles = pe.Node(SelectFiles({'anat': anat_file},
                                   base_directory=data_dir_path),
                       name="anat_selectfiles")

    bold_file = opj('{subject_id}', 'ses-{session}', 'bold', '{subject_id}_ses-{session}_run-{run}_bold.nii.gz')
    bold_selectfiles = pe.Node(SelectFiles({'bold': bold_file},
                                   base_directory=data_dir_path),
                       name="bold_selectfiles")

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

    #setting anat preprocessing nodes
    anat_preproc_wf = init_anat_preproc_wf()

    if commonspace_method=='pydpiper':
        raise ValueError("Depreciated. Use 'ants_dbm' instead as common space registration method.")
    elif commonspace_method=='ants_dbm':

        commonspace_reg = pe.Node(Function(input_names=['file_list', 'output_folder'],
                                  output_names=['ants_dbm_template', 'common_to_template_transform', 'template_to_common_transform'],
                                  function=commonspace_reg_function),
                         name='commonspace_reg')
        commonspace_reg.inputs.output_folder = output_folder+'/datasink/'

        #setting SelectFiles for the commonspace registration
        anat_to_template_inverse_warp = output_folder+'/'+opj('datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_{subject_id}_ses-{session}_anat_preproc*1InverseWarp.nii.gz')
        anat_to_template_warp = output_folder+'/'+opj('datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_{subject_id}_ses-{session}_anat_preproc*1Warp.nii.gz')
        anat_to_template_affine = output_folder+'/'+opj('datasink','ants_dbm_outputs','ants_dbm','output','secondlevel','secondlevel_{subject_id}_ses-{session}_anat_preproc*0GenericAffine.mat')
        ants_dbm_template_anat = output_folder+'/'+opj('datasink','ants_dbm_outputs','{ants_dbm_template}')
        common_to_template_transform = output_folder+'/'+opj('datasink','ants_dbm_outputs','{common_to_template_transform}')
        template_to_common_transform = output_folder+'/'+opj('datasink','ants_dbm_outputs','{template_to_common_transform}')

        commonspace_templates = {'anat_to_template_inverse_warp':anat_to_template_inverse_warp,'anat_to_template_warp': anat_to_template_warp, 'anat_to_template_affine': anat_to_template_affine, 'ants_dbm_template_anat': ants_dbm_template_anat, 'common_to_template_transform': common_to_template_transform, 'template_to_common_transform':template_to_common_transform}

        commonspace_selectfiles = pe.Node(SelectFiles(commonspace_templates),
                       name="commonspace_selectfiles")

        def transform_masks(reference_image,anat_to_template_inverse_warp, anat_to_template_affine,common_to_template_transform):
            import os
            cwd = os.getcwd()
            subject_id=os.path.basename(reference_image).split('_ses-')[0]
            session=os.path.basename(reference_image).split('_ses-')[1][0]
            filename_template = '%s_ses-%s' % (subject_id, session)
            input_image=os.environ["template_mask"]
            brain_mask='%s/%s_%s' % (cwd, filename_template, 'anat_mask.nii.gz')
            os.system('antsApplyTransforms -d 3 -i %s -t %s -t [%s,1] -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,anat_to_template_inverse_warp, anat_to_template_affine,common_to_template_transform,reference_image,brain_mask,))
            input_image=os.environ["WM_mask"]
            WM_mask='%s/%s_%s' % (cwd, filename_template, 'WM_mask.nii.gz')
            os.system('antsApplyTransforms -d 3 -i %s -t %s -t [%s,1] -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,anat_to_template_inverse_warp, anat_to_template_affine,common_to_template_transform,reference_image,WM_mask,))
            input_image=os.environ["CSF_mask"]
            CSF_mask='%s/%s_%s' % (cwd, filename_template, 'CSF_mask.nii.gz')
            os.system('antsApplyTransforms -d 3 -i %s -t %s -t [%s,1] -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,anat_to_template_inverse_warp, anat_to_template_affine,common_to_template_transform,reference_image,CSF_mask,))
            input_image=os.environ["atlas_labels"]
            anat_labels='%s/%s_%s' % (cwd, filename_template, 'atlas_labels.nii.gz')
            os.system('antsApplyTransforms -d 3 -i %s -t %s -t [%s,1] -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,anat_to_template_inverse_warp, anat_to_template_affine,common_to_template_transform,reference_image,anat_labels,))

            return brain_mask, WM_mask, CSF_mask, anat_labels

        transform_masks = pe.Node(Function(input_names=['reference_image','anat_to_template_inverse_warp', 'anat_to_template_affine','common_to_template_transform'],
                                  output_names=['brain_mask', 'WM_mask', 'CSF_mask', 'anat_labels'],
                                  function=transform_masks),
                         name='transform_masks')

        def commonspace_transforms(template_to_common_transform,anat_to_template_warp, anat_to_template_affine):
            return [template_to_common_transform,anat_to_template_warp, anat_to_template_affine],[0,0,0] #transforms_list,inverses
        commonspace_transforms_prep = pe.Node(Function(input_names=['template_to_common_transform','anat_to_template_warp','anat_to_template_affine'],
                                  output_names=['transforms_list','inverses'],
                                  function=commonspace_transforms),
                         name='commonspace_transforms_prep')

        bold_main_wf=init_bold_main_wf(tr=tr, tpattern=tpattern, apply_STC=apply_STC, data_dir_path=data_dir_path, bias_reg_script=bias_reg_script, coreg_script=coreg_script, commonspace_transform=True, SyN_SDC=True, isotropic_resampling=isotropic_resampling, upsampling=upsampling)

        workflow.connect([
            (anat_preproc_wf, joinnode_session, [("outputnode.preproc_anat", "file_list")]),
            (joinnode_sub_id, commonspace_reg, [
                ("file_list", "file_list"),
                ]),
            (infosub_id, commonspace_selectfiles, [
                ("subject_id", "subject_id"),
                ]),
            (infosession, commonspace_selectfiles, [
                ("session", "session")
                ]),
            (commonspace_reg, commonspace_selectfiles, [
                ("ants_dbm_template", "ants_dbm_template"),
                ("common_to_template_transform", "common_to_template_transform"),
                ("template_to_common_transform", "template_to_common_transform"),
                ]),
            (commonspace_reg, datasink, [
                ("ants_dbm_template", "ants_dbm_template"),
                ("common_to_template_transform", "common_to_template_transform"),
                ("template_to_common_transform", "template_to_common_transform"),
                ]),
            (commonspace_selectfiles, transform_masks, [
                ("common_to_template_transform", "common_to_template_transform"),
                ("anat_to_template_affine", "anat_to_template_affine"),
                ("anat_to_template_inverse_warp", "anat_to_template_inverse_warp"),
                ]),
            (anat_preproc_wf, transform_masks, [
                ("outputnode.preproc_anat", "reference_image"),
                ]),
            (anat_preproc_wf, bold_main_wf, [
                ("outputnode.preproc_anat", "inputnode.anat_preproc"),
                ]),
            (transform_masks, bold_main_wf, [
                ("anat_labels", 'inputnode.labels'),
                ("brain_mask", 'inputnode.anat_mask'),
                ("WM_mask", "inputnode.WM_mask"),
                ("CSF_mask", "inputnode.CSF_mask"),
                ]),
            (commonspace_selectfiles, commonspace_transforms_prep, [
                ("template_to_common_transform", "template_to_common_transform"),
                ("anat_to_template_affine", "anat_to_template_affine"),
                ("anat_to_template_warp", "anat_to_template_warp"),
                ]),
            (commonspace_selectfiles, datasink, [
                ("anat_to_template_affine", "anat_to_template_affine"),
                ("anat_to_template_warp", "anat_to_template_warp"),
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
        (anat_preproc_wf, datasink, [("outputnode.preproc_anat", "anat_preproc")]),
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
