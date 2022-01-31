from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

from .bold_ref import init_bold_reference_wf
from rabies.utils import slice_applyTransforms, Merge

def init_bold_preproc_trans_wf(opts, resampling_dim, name='bold_native_trans_wf'):
    """
    This workflow carries out the resampling of the original EPI timeseries into preprocessed timeseries.
    This is accomplished by applying at each frame a combined transform which accounts for previously estimated 
    motion correction and susceptibility distortion correction, together with the alignment to common space if
    the outputs are desired in common space. All transforms are concatenated into a single resampling operation
    to mitigate interpolation effects from repeated resampling.
    This workflow also carries the resampling of brain masks and labels from the reference atlas onto the 
    preprocessed EPI timeseries.

    Command line interface parameters:
        Resampling Options:
            The following options allow to resample the voxel dimensions for the preprocessed EPIs
            or for the anatomical images during registration.
            The resampling syntax must be 'dim1xdim2xdim3' (in mm), follwing the RAS axis convention
            (dim1=Right-Left, dim2=Anterior-Posterior, dim3=Superior-Inferior). If 'inputs_defined'
            is provided instead of axis dimensions, the original dimensions are preserved.

        --nativespace_resampling NATIVESPACE_RESAMPLING
                            Can specify a resampling dimension for the nativespace fMRI outputs.
                            (default: inputs_defined)
                            
        --commonspace_resampling COMMONSPACE_RESAMPLING
                            Can specify a resampling dimension for the commonspace fMRI outputs.
                            (default: inputs_defined)

    Workflow:
        parameters
            opts: command line interface parameters
            resampling_dim: specify the desired output voxel dimensions after resampling

        inputs
            name_source: a reference file for naming the output
            bold_file: the EPI timeseries to resample
            motcorr_params: the motion correction parameters
            transforms_list: a list of transforms to apply onto EPI timeseries, including 
                susceptibility distortion correction and resampling to common space
            inverses: a list specifying whether the inverse affine transforms should be 
                applied in transforms_list
            ref_file: a reference image in the targetted space for resampling. Should be the structural 
                image from the same session if outputs are in native space, or the atlas template for
                outputs in common space
            mask_transforms_list: the list of transforms to apply onto the atlas parcellations
                to overlap with the EPI
            mask_inverses: a list specifying whether the inverse affine transforms should be 
                applied in mask_transforms_list

        outputs
            bold: the preprocessed EPI timeseries
            bold_ref: a volumetric 3D EPI generated from the preprocessed timeseries
            brain_mask: the brain mask resampled onto preprocessed EPI timeseries
            WM_mask: the WM mask resampled onto preprocessed EPI timeseries
            CSF_mask: the CSF mask resampled onto preprocessed EPI timeseries
            vascular_mask: the vascular mask resampled onto preprocessed EPI timeseries
            labels: the atlas labels resampled onto preprocessed EPI timeseries
    """

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=[
        'name_source', 'bold_file', 'motcorr_params', 'transforms_list', 'inverses', 'ref_file',
        'mask_transforms_list', 'mask_inverses']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['bold', 'bold_ref', 'brain_mask', 'WM_mask', 'CSF_mask', 'vascular_mask', 'labels']),
        name='outputnode')

    bold_transform = pe.Node(slice_applyTransforms(
        rabies_data_type=opts.data_type), name='bold_transform', mem_gb=1*opts.scale_min_memory)
    bold_transform.inputs.apply_motcorr = (not opts.apply_slice_mc)
    bold_transform.inputs.resampling_dim = resampling_dim

    merge = pe.Node(Merge(rabies_data_type=opts.data_type, clip_negative=True), name='merge', mem_gb=4*opts.scale_min_memory)
    merge.plugin_args = {
        'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}

    # Generate a new BOLD reference
    bold_reference_wf = init_bold_reference_wf(opts=opts)

    WM_mask_to_EPI = pe.Node(MaskEPI(), name='WM_mask_EPI')
    WM_mask_to_EPI.inputs.name_spec = 'EPI_WM_mask'
    WM_mask_to_EPI.inputs.mask = str(opts.WM_mask)

    CSF_mask_to_EPI = pe.Node(MaskEPI(), name='CSF_mask_EPI')
    CSF_mask_to_EPI.inputs.name_spec = 'EPI_CSF_mask'
    CSF_mask_to_EPI.inputs.mask = str(opts.CSF_mask)

    vascular_mask_to_EPI = pe.Node(MaskEPI(), name='vascular_mask_EPI')
    vascular_mask_to_EPI.inputs.name_spec = 'EPI_vascular_mask'
    vascular_mask_to_EPI.inputs.mask = str(opts.vascular_mask)

    brain_mask_to_EPI = pe.Node(MaskEPI(), name='Brain_mask_EPI')
    brain_mask_to_EPI.inputs.name_spec = 'EPI_brain_mask'
    brain_mask_to_EPI.inputs.mask = str(opts.brain_mask)

    propagate_labels = pe.Node(MaskEPI(), name='prop_labels_EPI')
    propagate_labels.inputs.name_spec = 'EPI_anat_labels'
    propagate_labels.inputs.mask = str(opts.labels)

    workflow.connect([
        (inputnode, merge, [('name_source', 'header_source')]),
        (inputnode, bold_transform, [
            ('bold_file', 'in_file'),
            ('motcorr_params', 'motcorr_params'),
            ('transforms_list', 'transforms'),
            ('inverses', 'inverses'),
            ('ref_file', 'ref_file')
            ]),
        (bold_transform, merge, [('out_files', 'in_files')]),
        (merge, bold_reference_wf, [('out_file', 'inputnode.bold_file')]),
        (merge, outputnode, [('out_file', 'bold')]),
        (inputnode, brain_mask_to_EPI, [
            ('name_source', 'name_source'),
            ('mask_transforms_list', 'transforms'),
            ('mask_inverses', 'inverses'),
            ]),
        (bold_reference_wf, brain_mask_to_EPI, [
            ('outputnode.ref_image', 'ref_EPI')]),
        (brain_mask_to_EPI, outputnode, [
            ('EPI_mask', 'brain_mask')]),
        (inputnode, WM_mask_to_EPI, [
            ('name_source', 'name_source'),
            ('mask_transforms_list', 'transforms'),
            ('mask_inverses', 'inverses'),
            ]),
        (bold_reference_wf, WM_mask_to_EPI, [
            ('outputnode.ref_image', 'ref_EPI')]),
        (WM_mask_to_EPI, outputnode, [
            ('EPI_mask', 'WM_mask')]),
        (inputnode, CSF_mask_to_EPI, [
            ('name_source', 'name_source'),
            ('mask_transforms_list', 'transforms'),
            ('mask_inverses', 'inverses'),
            ]),
        (bold_reference_wf, CSF_mask_to_EPI, [
            ('outputnode.ref_image', 'ref_EPI')]),
        (CSF_mask_to_EPI, outputnode, [
            ('EPI_mask', 'CSF_mask')]),
        (inputnode, vascular_mask_to_EPI, [
            ('name_source', 'name_source'),
            ('mask_transforms_list', 'transforms'),
            ('mask_inverses', 'inverses'),
            ]),
        (bold_reference_wf, vascular_mask_to_EPI, [
            ('outputnode.ref_image', 'ref_EPI')]),
        (vascular_mask_to_EPI, outputnode, [
            ('EPI_mask', 'vascular_mask')]),
        (inputnode, propagate_labels, [
            ('name_source', 'name_source'),
            ('mask_transforms_list', 'transforms'),
            ('mask_inverses', 'inverses'),
            ]),
        (bold_reference_wf, propagate_labels, [
            ('outputnode.ref_image', 'ref_EPI')]),
        (propagate_labels, outputnode, [
            ('EPI_mask', 'labels')]),
        (bold_reference_wf, outputnode, [
            ('outputnode.ref_image', 'bold_ref')]),
    ])

    return workflow

class MaskEPIInputSpec(BaseInterfaceInputSpec):
    mask = File(exists=True, mandatory=True,
                desc="Mask to transfer to EPI space.")
    ref_EPI = File(exists=True, mandatory=True,
                   desc="Motion-realigned and SDC-corrected reference 3D EPI.")
    transforms = traits.List(desc="List of transforms to apply to every volume.")
    inverses = traits.List(
        desc="Define whether some transforms must be inverse, with a boolean list where true defines inverse e.g.[0,1,0]")
    name_spec = traits.Str(desc="Specify the name of the mask.")
    name_source = File(exists=True, mandatory=True,
                       desc='Reference BOLD file for naming the output.')


class MaskEPIOutputSpec(TraitedSpec):
    EPI_mask = traits.File(desc="The generated EPI mask.")


class MaskEPI(BaseInterface):

    input_spec = MaskEPIInputSpec
    output_spec = MaskEPIOutputSpec

    def _run_interface(self, runtime):
        import os
        import SimpleITK as sitk
        from rabies.utils import exec_applyTransforms

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.name_source).name.rsplit(".nii")

        if self.inputs.name_spec is None:
            new_mask_path = os.path.abspath(
                f'{filename_split[0]}_EPI_mask.nii.gz')
        else:
            new_mask_path = os.path.abspath(f'{filename_split[0]}_{self.inputs.name_spec}.nii.gz')

        exec_applyTransforms(self.inputs.transforms, self.inputs.inverses, self.inputs.mask, self.inputs.ref_EPI, new_mask_path, mask=True)
        sitk.WriteImage(sitk.ReadImage(
            new_mask_path, sitk.sitkInt16), new_mask_path)

        setattr(self, 'EPI_mask', new_mask_path)
        return runtime

    def _list_outputs(self):
        return {'EPI_mask': getattr(self, 'EPI_mask')}
