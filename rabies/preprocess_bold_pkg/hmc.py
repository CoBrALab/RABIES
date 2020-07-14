from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from .utils import SliceMotionCorrection


def init_bold_hmc_wf(slice_mc=False, name='bold_hmc_wf'):
    """
    This workflow estimates the motion parameters to perform HMC over the BOLD image.

    **Parameters**

        name : str
            Name of workflow (default: ``bold_hmc_wf``)

    **Inputs**

        bold_file
            BOLD series NIfTI file
        ref_image
            Reference image to which BOLD series is motion corrected

    **Outputs**

        movpar_file
            CSV file with antsMotionCorr motion parameters
    """
    import os

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['bold_file', 'ref_image']),
                        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['motcorr_params', 'slice_corrected_bold']),
        name='outputnode')

    # Head motion correction (hmc)
    motion_estimation = pe.Node(EstimateMotion(), name='ants_MC', mem_gb=1.1*float(os.environ["rabies_mem_scale"]))
    motion_estimation.plugin_args = {'qsub_args': '-pe smp %s' % (str(3*int(os.environ["min_proc"]))), 'overwrite': True}


    workflow.connect([
        (inputnode, motion_estimation, [('ref_image', 'ref_file'),
                              ('bold_file', 'in_file')]),
        (motion_estimation, outputnode, [('motcorr_params', 'motcorr_params')]),
    ])

    if slice_mc:
        slice_mc_n_procs=int(int(os.environ["local_threads"])/4)+1
        slice_mc_node = pe.Node(SliceMotionCorrection(n_procs=slice_mc_n_procs), name='slice_mc', mem_gb=1*slice_mc_n_procs, n_procs=slice_mc_n_procs)
        slice_mc_node.plugin_args = {'qsub_args': '-pe smp %s' % (str(3*int(os.environ["min_proc"]))), 'overwrite': True}

        #conducting a volumetric realignment before slice-specific mc to correct for larger head translations and rotations
        workflow.connect([
            (inputnode, slice_mc_node, [('ref_image', 'ref_file'),
                                        ('bold_file', 'name_source')]),
            (motion_estimation, slice_mc_node, [('mc_corrected_bold', 'in_file')]),
            (slice_mc_node, outputnode, [('mc_corrected_bold', 'slice_corrected_bold')]),
        ])

    return workflow



class EstimateMotionInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="4D EPI file")
    ref_file = File(exists=True, mandatory=True, desc="Reference image to which timeseries are realigned for motion estimation")

class EstimateMotionOutputSpec(TraitedSpec):
    motcorr_params = File(exists=True, desc="Motion estimation derived from antsMotionCorr")
    mc_corrected_bold = File(exists=True, desc="motion corrected time series")


class EstimateMotion(BaseInterface):
    """
    Runs ants motion correction interface and returns the motion estimation
    """

    input_spec = EstimateMotionInputSpec
    output_spec = EstimateMotionOutputSpec

    def _run_interface(self, runtime):
        import os
        from .utils import antsMotionCorr
        res = antsMotionCorr(in_file=self.inputs.in_file, ref_file=self.inputs.ref_file, second=False).run()
        csv_params = os.path.abspath(res.outputs.csv_params)

        setattr(self, 'csv_params', csv_params)
        setattr(self, 'mc_corrected_bold', os.path.abspath(res.outputs.mc_corrected_bold))

        return runtime

    def _list_outputs(self):
        return {'mc_corrected_bold': getattr(self, 'mc_corrected_bold'),
                'motcorr_params': getattr(self, 'csv_params')}
