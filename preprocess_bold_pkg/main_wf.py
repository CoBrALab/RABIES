"""
Orchestrating the BOLD-preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""

import os

import nibabel as nb

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from .hmc import init_bold_hmc_wf
from .utils import init_bold_reference_wf
from .resampling import init_bold_preproc_trans_wf
from .stc import init_bold_stc_wf
from .sdc import init_sdc_wf
from .bias_correction import bias_correction_wf
from .registration import init_bold_reg_wf
from .confounds import init_bold_confs_wf

def init_func_preproc_wf(bold_file, use_syn, TR, apply_STC=False, iterative_N4=True, motioncorr_24params=False, apply_GSR=False, name='main_wf'):

    """
    This workflow controls the functional preprocessing stages of the pipeline.


    **Parameters**

        bold_file
            BOLD series NIfTI file
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

    '''setting the workflow'''
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file', 'reversed_bold_file', 'anat_preproc', 'anat_mask',
                'anat_labels', 'WM_mask', 'CSF_mask']),
        name='inputnode')


    outputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file', 'bold_ref', 'skip_vols', 'hmc_xforms', 'output_warped_bold', 'itk_bold_to_anat', 'itk_anat_to_bold',
                'resampled_bold', 'resampled_ref_bold', 'hmc_movpar_file', 'cleaned_bold', 'GSR_cleaned_bold', 'EPI_labels', 'confounds_csv']),
        name='outputnode')


    bold_reference_wf = init_bold_reference_wf()
    bias_cor_wf = bias_correction_wf(iterative=iterative_N4)

    if apply_STC:
        bold_stc_wf = init_bold_stc_wf(TR=(str(TR)+'s'))

    # BOLD buffer: an identity used as a pointer to the STC data for further use.
    boldbuffer = pe.Node(niu.IdentityInterface(fields=['bold_file']), name='boldbuffer')

    # HMC on the BOLD
    bold_hmc_wf = init_bold_hmc_wf(name='bold_hmc_wf')

    bold_reg_wf = init_bold_reg_wf(SyN_reg=use_syn)

    # Apply transforms in 1 shot
    bold_bold_trans_wf = init_bold_preproc_trans_wf(
        use_fieldwarp=True,
        name='bold_bold_trans_wf'
    )

    bold_confs_wf = init_bold_confs_wf(apply_GSR=apply_GSR, TR=TR, motioncorr_24params=motioncorr_24params, name="bold_confs_wf")


    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (inputnode, bold_reference_wf, [('bold_file', 'inputnode.bold_file')]),
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
        (inputnode, outputnode, [('bold_file', 'bold_file')]),
        (inputnode, bold_reg_wf, [
            ('anat_preproc', 'inputnode.anat_preproc'),
            ('anat_mask', 'inputnode.anat_mask')]),
        (bold_reg_wf, outputnode, [
            ('outputnode.itk_bold_to_anat', 'itk_bold_to_anat'),
            ('outputnode.itk_anat_to_bold', 'itk_anat_to_bold'),
            ('outputnode.output_warped_bold', 'output_warped_bold'),
            ]),
        (boldbuffer, bold_bold_trans_wf, [('bold_file', 'inputnode.bold_file')]),
        (inputnode, bold_bold_trans_wf, [('bold_file', 'inputnode.name_source')]),
        (bold_hmc_wf, bold_bold_trans_wf, [('outputnode.xforms', 'inputnode.hmc_xforms')]),
        (bold_bold_trans_wf, outputnode, [
            ('outputnode.bold_ref', 'resampled_ref_bold'),
            ('outputnode.bold', 'resampled_bold'),
            ]),
        (inputnode, bold_confs_wf, [('anat_mask', 'inputnode.t1_mask'),
            ('WM_mask', 'inputnode.WM_mask'),
            ('CSF_mask', 'inputnode.CSF_mask'),
            ('anat_labels', 'inputnode.t1_labels'),
            ]),
        (bold_bold_trans_wf, bold_confs_wf, [('outputnode.bold', 'inputnode.bold'),
            ('outputnode.bold_ref', 'inputnode.ref_bold'),
            ]),
        (bold_hmc_wf, bold_confs_wf, [('outputnode.movpar_file', 'inputnode.movpar_file'),
            ]),
        (bold_confs_wf, outputnode, [
            ('outputnode.cleaned_bold', 'cleaned_bold'),
            ('outputnode.GSR_cleaned_bold', 'GSR_cleaned_bold'),
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

    workflow.connect([
        (bias_cor_wf, bold_reg_wf, [
              ('outputnode.corrected_EPI', 'inputnode.ref_bold_brain')]),
        ])


    workflow.connect([
        (bold_reg_wf, bold_bold_trans_wf, [('outputnode.itk_bold_to_anat', 'inputnode.fieldwarp')]),
        ])

    return workflow


def _get_series_len(bold_fname):
    from .utils import _get_vols_to_discard
    img = nb.load(bold_fname)
    if len(img.shape) < 4:
        return 1

    skip_vols = _get_vols_to_discard(img)

    return img.shape[3] - skip_vols
