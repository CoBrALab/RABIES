from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def init_pydpiper_wf(model_script_path, name='pydpiper_wf'):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['csv_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['transform_csv', 'pydpiper_directory']), name='outputnode')

    pydpiper = pe.Node(Pydpiper(model_script_path=model_script_path),
                 name="pydpiper")


    workflow.connect([
        (inputnode, pydpiper, [("csv_file", "csv_file")]),
        (pydpiper, outputnode, [("transform_csv", "transform_csv"),
                                ("pydpiper_directory", "pydpiper_directory")]),
    ])

    return workflow


class PydpiperInputSpec(BaseInterfaceInputSpec):
    csv_file = File(exists=True, mandatory=True, desc="Anatomical image to preprocess")
    model_script_path = File(exists=True, mandatory=True, desc="Anatomical image to preprocess")

class PydpiperOutputSpec(TraitedSpec):
    transform_csv = File(exists=True, desc="CSV with path info.")
    pydpiper_directory = traits.Str(exists=True, desc="Path to pydpiper directory.")

class Pydpiper(BaseInterface):

    input_spec = PydpiperInputSpec
    output_spec = PydpiperOutputSpec

    def _run_interface(self, runtime):
        import os
        cwd = os.getcwd()
        os.system('bash %s %s' % (self.inputs.model_script_path,self.inputs.csv_file,))
        setattr(self, 'transform_csv', '%s/transforms.csv' % (cwd,))
        setattr(self, 'pydpiper_directory', cwd)
        return runtime

    def _list_outputs(self):
        return {'transform_csv': getattr(self, 'transform_csv'),
                'pydpiper_directory': getattr(self, 'pydpiper_directory')}
