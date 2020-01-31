from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)



def init_bold_hmc_wf(name='bold_hmc_wf'):
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
        niu.IdentityInterface(fields=['motcorr_params']),
        name='outputnode')

    # Head motion correction (hmc)
    motion_estimation = pe.Node(EstimateMotion(), name='ants_MC', mem_gb=3)
    motion_estimation.plugin_args = {'qsub_args': '-pe smp %s' % (str(3*int(os.environ["min_proc"]))), 'overwrite': True}


    workflow.connect([
        (inputnode, motion_estimation, [('ref_image', 'ref_file'),
                              ('bold_file', 'in_file')]),
        (motion_estimation, outputnode, [('motcorr_params', 'motcorr_params')]),
    ])

    return workflow



class EstimateMotionInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="4D EPI file")
    ref_file = File(exists=True, mandatory=True, desc="Reference image to which timeseries are realigned for motion estimation")

class EstimateMotionOutputSpec(TraitedSpec):
    motcorr_params = File(exists=True, desc="Motion estimation derived from antsMotionCorr")
    csv_params = File(exists=True, desc="CSV file with motion parameters")


class EstimateMotion(BaseInterface):
    """
    Runs ants motion correction interface and returns the motion estimation
    """

    input_spec = EstimateMotionInputSpec
    output_spec = EstimateMotionOutputSpec

    def _run_interface(self, runtime):
        import os
        import nibabel as nb
        from .utils import antsMotionCorr
        res = antsMotionCorr(in_file=self.inputs.in_file, ref_file=self.inputs.ref_file, second=False).run()
        csv_params = os.path.abspath(res.outputs.csv_params)

        setattr(self, 'csv_params', csv_params)

        return runtime

    def _list_outputs(self):
        return {'motcorr_params': getattr(self, 'csv_params')}
