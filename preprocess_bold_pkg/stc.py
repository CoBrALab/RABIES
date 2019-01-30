# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Slice-Timing Correction (STC) of BOLD images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_bold_stc_wf

"""
from nipype import logging
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, afni
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, SimpleInterface
)

# pylint: disable=R0914
def init_bold_stc_wf(TR, name='bold_stc_wf'):
    """
    This workflow performs :abbr:`STC (slice-timing correction)` over the input
    :abbr:`BOLD (blood-oxygen-level dependent)` image.

    .. workflow::
        :graph2use: orig
        :simple_form: yes

        from fmriprep.workflows.bold import init_bold_stc_wf
        wf = init_bold_stc_wf(
            metadata={"RepetitionTime": 2.0,
                      "SliceTiming": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]},
            )

    **Parameters**

        name : str
            Name of workflow (default: ``bold_stc_wf``)

    **Inputs**

        bold_file
            BOLD series NIfTI file
        skip_vols
            Number of non-steady-state volumes detected at beginning of ``bold_file``

    **Outputs**

        stc_file
            Slice-timing corrected BOLD series NIfTI file

    """
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['bold_file', 'skip_vols']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['stc_file']), name='outputnode')


    # It would be good to fingerprint memory use of afni.TShift
    slice_timing_correction = pe.Node(
        afni.TShift(outputtype='NIFTI_GZ', tr=TR),
        name='slice_timing_correction')

    copy_xform = pe.Node(CopyXForm(), name='copy_xform')

    workflow.connect([
        (inputnode, slice_timing_correction, [('bold_file', 'in_file'),
                                              ('skip_vols', 'ignore')]),
        (slice_timing_correction, copy_xform, [('out_file', 'in_file')]),
        (inputnode, copy_xform, [('bold_file', 'hdr_file')]),
        (copy_xform, outputnode, [('out_file', 'stc_file')]),
    ])

    return workflow




from nipype.utils.filemanip import fname_presuffix
import shutil


class CopyXFormInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='the file we get the data from')
    hdr_file = File(exists=True, mandatory=True, desc='the file we get the header from')


class CopyXFormOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='written file path')


class CopyXForm(SimpleInterface):
    """
    Copy the x-form matrices from `hdr_file` to `out_file`.
    """
    input_spec = CopyXFormInputSpec
    output_spec = CopyXFormOutputSpec

    def _run_interface(self, runtime):
        out_name = fname_presuffix(self.inputs.in_file,
                                   suffix='_xform',
                                   newpath=runtime.cwd)
        # Copy and replace header
        shutil.copy(self.inputs.in_file, out_name)
        _copyxform(self.inputs.hdr_file, out_name)
        self._results['out_file'] = out_name
        return runtime



def _copyxform(ref_image, out_image, message=None):
    # Read in reference and output
    # Use mmap=False because we will be overwriting the output image
    import nibabel as nb
    import numpy as np
    resampled = nb.load(out_image, mmap=False)
    orig = nb.load(ref_image)

    if not np.allclose(orig.affine, resampled.affine):
        LOG.debug(
            'Affines of input and reference images do not match, '
            'FMRIPREP will set the reference image headers. '
            'Please, check that the x-form matrices of the input dataset'
            'are correct and manually verify the alignment of results.')

    # Copy xform infos
    qform, qform_code = orig.header.get_qform(coded=True)
    sform, sform_code = orig.header.get_sform(coded=True)
    header = resampled.header.copy()
    header.set_qform(qform, int(qform_code))
    header.set_sform(sform, int(sform_code))
    header['descrip'] = 'xform matrices modified by %s.' % (message or '(unknown)')

    newimg = resampled.__class__(resampled.get_data(), orig.affine, header)
    newimg.to_filename(out_image)
