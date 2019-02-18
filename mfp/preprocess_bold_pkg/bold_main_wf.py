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
from .resampling import init_bold_preproc_trans_wf
from .stc import init_bold_stc_wf
from .sdc import init_sdc_wf
from .bias_correction import bias_correction_wf
from .registration import init_bold_reg_wf
from .confounds import init_bold_confs_wf

def init_bold_main_wf(data_dir_path, TR, run_iter=None, anat_files_csv=None, bias_reg_script='Rigid', coreg_script='SyN', SyN_SDC=True, apply_STC=False, iterative_N4=True,
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

    inputnode = pe.Node(niu.IdentityInterface(fields=['subject_id', 'session', 'run_iter', 'bold', 'anat_preproc', 'anat_mask', 'WM_mask', 'CSF_mask', 'labels']),
                      name="inputnode")

    outputnode = pe.Node(niu.IdentityInterface(
                fields=['input_bold', 'bold_ref', 'skip_vols','hmc_xforms', 'corrected_EPI', 'output_warped_bold', 'itk_bold_to_anat', 'itk_anat_to_bold',
                        'resampled_bold', 'resampled_ref_bold', 'hmc_movpar_file', 'EPI_labels', 'confounds_csv']),
                name='outputnode')

    '''


    if bold_preproc_only:
        import pandas as pd
        data_df=pd.read_csv(anat_files_csv, sep=',')
        subject_list=data_df['subject_id'].values.tolist()
        session=data_df['session'].values.tolist()
        run_list=data_df['run_num'].values.tolist()
        anat_preproc_list=data_df['anat_preproc'].values.tolist()
        anat_mask_list=data_df['anat_mask'].values.tolist()
        WM_mask_list=data_df['WM_mask'].values.tolist()
        CSF_mask_list=data_df['CSF_mask'].values.tolist()
        labels_list=data_df['labels'].values.tolist()

        #create a dictionary with list of bold run numbers for each subject
        run_iter={}
        for i in range(len(subject_list)):
            run_iter[subject_list[i]] = list(range(1,int(run_list[i])+1))

        inputnode.iterables = [('subject_id', subject_list), ('session', session), ('anat_preproc', anat_preproc_list),
                                ('anat_mask', anat_mask_list), ('WM_mask', WM_mask_list), ('CSF_mask', CSF_mask_list), ('labels', labels_list)]
        inputnode.synchronize = True
    '''


    bold_reference_wf = init_bold_reference_wf()
    bias_cor_wf = bias_correction_wf(iterative=iterative_N4, bias_reg_script=bias_reg_script)

    if apply_STC:
        bold_stc_wf = init_bold_stc_wf(TR=TR)

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
    bold_confs_wf = init_bold_confs_wf(SyN_SDC=SyN_SDC, aCompCor_method=aCompCor_method, TR=TR, name="bold_confs_wf")


    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (inputnode, bold_reference_wf, [('bold', 'inputnode.bold_file')]),
        (bold_reference_wf, bias_cor_wf, [
            ('outputnode.ref_image', 'inputnode.ref_EPI')]),
        (inputnode, bias_cor_wf, [
            ('anat_preproc', 'inputnode.anat'),
            ('anat_mask', 'inputnode.anat_mask'),
            ]),
        (bold_reference_wf, bold_hmc_wf, [
            ('outputnode.ref_image', 'inputnode.ref_image'),
            ('outputnode.bold_file', 'inputnode.bold_file')]),
        (bold_hmc_wf, outputnode, [
            ('outputnode.xforms', 'hmc_xforms'),
            ('outputnode.movpar_file', 'hmc_movpar_file')]),
        (bold_reference_wf, outputnode, [
            ('outputnode.ref_image', 'bold_ref')]),
        (inputnode, bold_reg_wf, [
            ('anat_preproc', 'inputnode.anat_preproc'),
            ('anat_mask', 'inputnode.anat_mask')]),
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
        (inputnode, bold_bold_trans_wf, [('bold', 'inputnode.name_source')]),
        (bold_hmc_wf, bold_bold_trans_wf, [('outputnode.xforms', 'inputnode.hmc_xforms')]),
        (bold_bold_trans_wf, outputnode, [
            ('outputnode.bold_ref', 'resampled_ref_bold'),
            ('outputnode.bold', 'resampled_bold'),
            ]),
        (inputnode, bold_confs_wf, [('anat_mask', 'inputnode.t1_mask'),
            ('WM_mask', 'inputnode.WM_mask'),
            ('CSF_mask', 'inputnode.CSF_mask'),
            ('labels', 'inputnode.t1_labels'),
            ]),
        (bold_bold_trans_wf, bold_confs_wf, [('outputnode.bold', 'inputnode.bold'),
            ('outputnode.bold_ref', 'inputnode.ref_bold'),
            ]),
        (bold_hmc_wf, bold_confs_wf, [('outputnode.movpar_file', 'inputnode.movpar_file'),
            ]),
        (bold_confs_wf, outputnode, [
            ('outputnode.EPI_labels', 'EPI_labels'),
            ('outputnode.confounds_csv', 'confounds_csv'),
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


    return workflow
