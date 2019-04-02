from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from .utils import applyTransforms, init_bold_reference_wf, Merge


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
        'name_source', 'bold_file', 'hmc_xforms', 'fieldwarp','ref_file']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['bold', 'bold_ref', 'bold_ref_brain']),
        name='outputnode')


    bold_transform = pe.Node(
        applyTransforms(use_fieldwarp=use_fieldwarp), name='bold_transform')


    merge = pe.Node(Merge(), name='merge')

    # Generate a new BOLD reference
    bold_reference_wf = init_bold_reference_wf()

    workflow.connect([
        (inputnode, merge, [('name_source', 'header_source')]),
        (inputnode, bold_transform, [
            ('bold_file', 'in_file'),
            ('hmc_xforms', 'xforms'),
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
