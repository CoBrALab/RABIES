from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from .utils import SliceMotionCorrection, antsMotionCorr

def init_bold_hmc_wf(opts, name='bold_hmc_wf'):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['bold_file', 'ref_image']),
                        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['motcorr_params', 'slice_corrected_bold']),
        name='outputnode')

    # Head motion correction (hmc)
    motion_estimation = pe.Node(antsMotionCorr(prebuilt_option=opts.HMC_option,transform_type='Rigid', second=False, rabies_data_type=opts.data_type),
                         name='ants_MC', mem_gb=1.1*opts.scale_min_memory)
    motion_estimation.plugin_args = {
        'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}

    workflow.connect([
        (inputnode, motion_estimation, [('ref_image', 'ref_file'),
                                        ('bold_file', 'in_file')]),
        (motion_estimation, outputnode, [
         ('csv_params', 'motcorr_params')]),
    ])

    if opts.apply_slice_mc:
        slice_mc_n_procs = int(opts.local_threads/4)+1
        slice_mc_node = pe.Node(SliceMotionCorrection(n_procs=slice_mc_n_procs),
                                name='slice_mc', mem_gb=1*slice_mc_n_procs, n_procs=slice_mc_n_procs)
        slice_mc_node.plugin_args = {
            'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}

        # conducting a volumetric realignment before slice-specific mc to correct for larger head translations and rotations
        workflow.connect([
            (inputnode, slice_mc_node, [('ref_image', 'ref_file'),
                                        ('bold_file', 'name_source')]),
            (motion_estimation, slice_mc_node, [
             ('mc_corrected_bold', 'in_file')]),
            (slice_mc_node, outputnode, [
             ('mc_corrected_bold', 'slice_corrected_bold')]),
        ])

    return workflow
