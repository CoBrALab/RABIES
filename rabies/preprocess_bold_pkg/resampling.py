from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from .utils import slice_applyTransforms, init_bold_reference_wf, Merge
from nipype.interfaces.utility import Function

def init_bold_preproc_trans_wf(name='bold_preproc_trans_wf',
                               use_fieldwarp=True):
    """
    This workflow resamples the input fMRI in its native (original)
    space in a "single shot" from the original BOLD series.


    **Parameters**


        name : str
            Name of workflow (default: ``bold_mni_trans_wf``)
        use_fieldwarp : bool
            Include SDC warp in single-shot transform from BOLD to MNI

    **Inputs**

        bold_file
            Individual 3D volumes, not motion corrected
        bold_mask
            Skull-stripping mask of reference image
        name_source
            BOLD series NIfTI file
            Used to recover original information lost during processing
        hmc_xforms
            List of affine transforms aligning each volume to ``ref_image`` in ITK format
        fieldwarp
            a :abbr:`DFM (displacements field map)` in ITK format

    **Outputs**

        bold
            BOLD series, resampled in native space, including all preprocessing
        bold_mask

        bold_ref
            BOLD reference image: an average-like 3D image of the time-series
        bold_ref_brain
            Same as ``bold_ref``, but once the brain mask has been applied

    """
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=[
        'name_source', 'bold_file', 'motcorr_params', 'fieldwarp','ref_file']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['bold', 'bold_ref', 'bold_ref_brain']),
        name='outputnode')


    bold_transform = pe.Node(
        slice_applyTransforms(use_fieldwarp=use_fieldwarp), name='bold_transform')

    merge = pe.Node(Merge(), name='merge')

    # Generate a new BOLD reference
    bold_reference_wf = init_bold_reference_wf()

    workflow.connect([
        (inputnode, merge, [('name_source', 'header_source')]),
        (inputnode, bold_transform, [
            ('bold_file', 'in_file'),
            ('motcorr_params', 'motcorr_params'),
            ]),
        (bold_transform, merge, [('out_files', 'in_files')]),
        (merge, bold_reference_wf, [('out_file', 'inputnode.bold_file')]),
        (merge, outputnode, [('out_file', 'bold')]),
        (bold_reference_wf, outputnode, [
            ('outputnode.ref_image', 'bold_ref')]),
    ])


    if use_fieldwarp:
        workflow.connect([
            (inputnode, bold_transform, [('fieldwarp', 'fieldwarp')]),
            (inputnode, bold_transform, [('ref_file', 'ref_file')]),
        ])

    return workflow


def init_bold_commonspace_trans_wf(name='bold_commonspace_trans_wf'):
    """
    Apply transforms of the EPI to commonspace, resampling the resolution to the target template.
    """
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=[
        'bold_file', 'commonspace_affine', 'commonspace_warp', 'commonspace_template']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['commonspace_bold']),
        name='outputnode')

    transform_commonspace = pe.Node(Function(input_names=["reference_image",'input_image','commonspace_affine', 'commonspace_warp'],
                              output_names=["output_image"],
                              function=apply_transform),
                     name='transform_commonspace')

    workflow.connect([
        (inputnode, transform_commonspace, [
            ('bold_file', 'input_image'),
            ('commonspace_template', 'reference_image'),
            ('commonspace_affine', 'commonspace_affine'),
            ('commonspace_warp', 'commonspace_warp'),
            ]),
        (transform_commonspace, outputnode, [
            ('output_image', 'commonspace_bold'),
            ]),
    ])

    return workflow


def apply_transform(reference_image,input_image,commonspace_affine,commonspace_warp):
    import os
    subject_id=os.path.basename(input_image).split('_ses-')[0]
    session=os.path.basename(input_image).split('_ses-')[1][0]
    run=os.path.basename(input_image).split('_run-')[1][0]
    cwd = os.getcwd()
    output_image='%s/%s_ses-%s_run-%s_commonspace.nii.gz' % (cwd, subject_id, session, run)
    os.system('antsApplyTransforms -d 3 -i %s -t %s -t %s -r %s -o %s -e 3 --verbose' % (input_image,commonspace_warp,commonspace_affine,reference_image,output_image,))

    return output_image
