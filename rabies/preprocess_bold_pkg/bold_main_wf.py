"""
Orchestrating the BOLD-preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""

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
                        aCompCor_method='50%', bold_preproc_only=False, name='bold_main_wf'):

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

    inputnode = pe.Node(niu.IdentityInterface(fields=['subject_id', 'bold', 'anat_preproc', 'anat_mask', 'WM_mask', 'CSF_mask', 'labels', 'commonspace_affine', 'commonspace_warp', 'commonspace_template']),
                      name="inputnode")

    outputnode = pe.Node(niu.IdentityInterface(
                fields=['input_bold', 'bold_ref', 'skip_vols','motcorr_params', 'corrected_EPI', 'output_warped_bold', 'itk_bold_to_anat', 'itk_anat_to_bold', 'commonspace_bold',
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

    # Apply transforms in 1 shot
    bold_bold_trans_wf = init_bold_preproc_trans_wf(
        use_fieldwarp=True,
        name='bold_bold_trans_wf'
    )

    bold_confs_wf = init_bold_confs_wf(SyN_SDC=SyN_SDC, aCompCor_method=aCompCor_method, name="bold_confs_wf")

    if bold_preproc_only:
        #takes as input a csv file with the path to the bold file, the paths to the anatomical template, the brain/CSF/WM masks, and the associated label file
        #the workflow will iterate simply through each subject

        sub_setup = pe.Node(Function(input_names=['subject_id', 'data_csv'],
                                  output_names=['bold', 'anat_preproc', 'anat_mask', 'WM_mask', 'CSF_mask', 'labels'],
                                  function=get_sub_files),
                         name='sub_setup')
        sub_setup.inputs.data_csv=data_csv

        workflow.connect([
            (inputnode, sub_setup, [('subject_id', 'subject_id')]),
            (sub_setup, bold_reference_wf, [('bold', 'inputnode.bold_file')]),
            (sub_setup, bias_cor_wf, [
                ('anat_preproc', 'inputnode.anat'),
                ('anat_mask', 'inputnode.anat_mask'),
                ]),
            (sub_setup, bold_reg_wf, [
                ('anat_preproc', 'inputnode.anat_preproc'),
                ('anat_mask', 'inputnode.anat_mask')]),
            (sub_setup, bold_bold_trans_wf, [('bold', 'inputnode.name_source')]),
            (sub_setup, bold_confs_wf, [('anat_mask', 'inputnode.t1_mask'),
                ('WM_mask', 'inputnode.WM_mask'),
                ('CSF_mask', 'inputnode.CSF_mask'),
                ('labels', 'inputnode.t1_labels'),
                ]),
            ])
    else:
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
            ])


    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
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
        (bold_reg_wf, bold_bold_trans_wf, [('outputnode.itk_bold_to_anat', 'inputnode.fieldwarp')]),
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

    workflow.connect([
        (bold_reg_wf, bold_bold_trans_wf, [
            ('outputnode.output_warped_bold', 'inputnode.ref_file')]),
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
        bold_commonspace_trans_wf=init_bold_commonspace_trans_wf()
        workflow.connect([
            (inputnode, bold_commonspace_trans_wf, [
                ('commonspace_affine', 'inputnode.commonspace_affine'),
                ('commonspace_warp', 'inputnode.commonspace_warp'),
                ('commonspace_template', 'inputnode.commonspace_template'),
                ]),
            (bold_bold_trans_wf, bold_commonspace_trans_wf, [
                ('outputnode.bold', 'inputnode.bold_file'),
                ]),
            (bold_commonspace_trans_wf, outputnode, [
                ('outputnode.commonspace_bold', 'commonspace_bold'),
                ]),
        ])

    return workflow



def get_sub_files(subject_id, data_csv):
    import os
    import pandas as pd
    import numpy as np
    data_df=pd.read_csv(data_csv, sep=',')
    subject_array=np.asarray(data_df['subject_id'].values)
    idx=int(np.where(subject_array==subject_id)[0])

    if os.path.isfile(data_df['bold_file'].values[idx]):
        return [data_df['bold_file'].values[idx],os.environ["template_anat"],os.environ["template_mask"],os.environ["WM_mask"],
                    os.environ["CSF_mask"],os.environ["atlas_labels"]]
    else:
        raise ValueError('Input bold file %s does not exist.' % (data_df['bold_file'].values[idx]))

'''
    if SyN_SDC:
        workflow.connect([
            (inputnode, bold_bold_trans_wf, [
                ('anat_preproc', 'inputnode.ref_file')]),
            ])
    else:
        workflow.connect([
            (bold_reference_wf, bold_bold_trans_wf, [
                ('outputnode.ref_image', 'inputnode.ref_file')]),
            ])
'''
