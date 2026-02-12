from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
import os
from .bold_ref import init_bold_reference_wf
from ..utils import ResampleVolumes,ResampleMask

def init_bold_preproc_trans_wf(opts, resampling_dim, name='bold_native_trans_wf'):
    # bold_resampling_head_start
    """
    This workflow carries out the resampling of the original EPI timeseries into preprocessed timeseries.
    This is accomplished by applying at each frame a combined transform which accounts for previously estimated 
    motion correction and susceptibility distortion correction, together with the alignment to common space (the
    exact combination of transforms depends on which anatomical preprocessed timeseries are resampled into). 
    All transforms are concatenated into a single resampling operation to mitigate interpolation effects from 
    repeated resampling. 

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

        outputs
            bold: the preprocessed EPI timeseries
            bold_ref: a volumetric 3D EPI generated from the preprocessed timeseries
    """
    # bold_resampling_head_end

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=[
        'name_source', 'bold_file', 'motcorr_params', 'transforms_list', 'inverses', 'ref_file']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['bold', 'bold_ref']),
        name='outputnode')

    n_procs=int(os.environ['RABIES_ITK_NUM_THREADS'])
    bold_transform = pe.Node(ResampleVolumes(
        rabies_data_type=opts.data_type, clip_negative=True), 
        name='bold_transform', mem_gb=4*opts.scale_min_memory, n_procs=n_procs)
    bold_transform.inputs.resampling_dim = resampling_dim
    bold_transform.inputs.interpolation = opts.interpolation_sitk
    bold_transform.plugin_args = {
        'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}
    if opts.no_HMC:
        bold_transform.inputs.apply_motcorr = False
    else:
        bold_transform.inputs.apply_motcorr = True

    # Generate a new BOLD reference
    bold_reference_wf = init_bold_reference_wf(opts=opts)

    workflow.connect([
        (inputnode, bold_transform, [
            ('bold_file', 'in_file'),
            ('motcorr_params', 'motcorr_params'),
            ('transforms_list', 'transforms_3d_files'),
            ('inverses', 'inverses_3d'),
            ('ref_file', 'ref_file'),
            ('name_source', 'name_source'),
            ]),
        (bold_transform, bold_reference_wf, [('resampled_file', 'inputnode.bold_file')]),
        (bold_transform, outputnode, [('resampled_file', 'bold')]),
        (bold_reference_wf, outputnode, [
            ('outputnode.ref_image', 'bold_ref')]),
    ])

    return workflow


def init_mask_preproc_trans_wf(opts, name='mask_native_trans_wf'):
    # mask_resampling_head_start
    """
    This workflow carries the resampling of brain masks from the reference 
    commonspace onto the preprocessed EPI timeseries.

    Workflow:
        parameters
            opts: command line interface parameters

        inputs
            name_source: a reference file for naming the output
            ref_file: a reference image in the targetted space for resampling. Should be the structural 
                image from the same session if outputs are in native space, or the atlas template for
                outputs in common space
            mask_transforms_list: the list of transforms to apply onto the atlas parcellations
                to overlap with the EPI
            mask_inverses: a list specifying whether the inverse affine transforms should be 
                applied in mask_transforms_list

        outputs
            brain_mask: the brain mask resampled onto preprocessed EPI timeseries
            WM_mask: the WM mask resampled onto preprocessed EPI timeseries
            CSF_mask: the CSF mask resampled onto preprocessed EPI timeseries
            vascular_mask: the vascular mask resampled onto preprocessed EPI timeseries
    """
    # mask_resampling_head_end

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=[
        'name_source', 'ref_file',
        'mask_transforms_list', 'mask_inverses'
        ]),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['brain_mask', 'WM_mask', 'CSF_mask', 'vascular_mask']),
        name='outputnode')

    # integrate a node to resample each mask, only if the mask exists
    for opt_key in ['brain_mask', 'WM_mask','CSF_mask','vascular_mask']:
        opt_file = getattr(opts, opt_key)
        if opt_file is not None:
            mask_to_EPI = pe.Node(ResampleMask(), name=opt_key+'_resample')
            mask_to_EPI.inputs.name_suffix = opt_key+'_resampled'
            mask_to_EPI.inputs.mask_file = str(opt_file)

            workflow.connect([
                (inputnode, mask_to_EPI, [
                    ('name_source', 'name_source'),
                    ('mask_transforms_list', 'transforms'),
                    ('mask_inverses', 'inverses'),
                    ('ref_file', 'ref_file'),
                    ]),
                (mask_to_EPI, outputnode, [
                    ('resampled_file', opt_key)]),
            ])

    return workflow
