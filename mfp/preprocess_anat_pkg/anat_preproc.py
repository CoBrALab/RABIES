from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def init_anat_preproc_wf(name='anat_preproc_wf'):
    '''
    This workflow execute required preprocessing steps for anatomical scans at
    the single subject levels.
    This includes:
    1-conversion of anat to minc format
    2-single scan anatomical preprocessing (N4 bias correction, brain mask, ...)
    with mouse-preprocessing-v3.sh script.
    The inputs for pydpiper are then ready.
    '''

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['anat_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['preproc_anat']), name='outputnode')

    anat_preproc = pe.Node(AnatPreproc(), name='Anat_Preproc')


    workflow.connect([
        (inputnode, anat_preproc, [("anat_file", "nii_anat")]),
        (anat_preproc, outputnode, [("preproc_anat", "preproc_anat")]),
    ])

    return workflow



class AnatPreprocInputSpec(BaseInterfaceInputSpec):
    nii_anat = File(exists=True, mandatory=True, desc="Anatomical image to preprocess")

class AnatPreprocOutputSpec(TraitedSpec):
    preproc_anat = File(exists=True, desc="Preprocessed anatomical image.")

class AnatPreproc(BaseInterface):

    input_spec = AnatPreprocInputSpec
    output_spec = AnatPreprocOutputSpec

    def _run_interface(self, runtime):
        import os
        cwd = os.getcwd()
        os.system('mkdir -p %s/mnc_anat/' % (cwd,))
        os.system('mkdir -p %s/anat_preproc/' % (cwd,))
        anat_file=os.path.basename(self.inputs.nii_anat).split('.')[0]
        os.system('nii2mnc %s %s/mnc_anat/%s.mnc' % (self.inputs.nii_anat,cwd,anat_file))
        dir_path = os.path.dirname(os.path.realpath(__file__))
        os.system('bash %s/../shell_scripts/mouse-preprocessing-v3.sh %s/mnc_anat/%s.mnc %s/anat_preproc/%s_preproc.mnc' % (dir_path,cwd,anat_file,cwd,anat_file))

        setattr(self, 'preproc_anat', '%s/anat_preproc/%s_preproc.mnc' % (cwd,anat_file))
        return runtime

    def _list_outputs(self):
        return {'preproc_anat': getattr(self, 'preproc_anat')}
