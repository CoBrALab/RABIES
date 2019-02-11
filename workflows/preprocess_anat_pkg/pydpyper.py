from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def init_pydpyper_wf(model_script_path, name='pydpyper_wf'):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['csv_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['transform_csv', 'pydpyper_directory']), name='outputnode')

    pydpyper = pe.Node(Pydpyper(model_script_path=model_script_path),
                 name="pydpyper")


    workflow.connect([
        (inputnode, pydpyper, [("csv_file", "csv_file")]),
        (pydpyper, outputnode, [("transform_csv", "transform_csv"),
                                ("pydpyper_directory", "pydpyper_directory")]),
    ])

    return workflow


class PydpyperInputSpec(BaseInterfaceInputSpec):
    csv_file = File(exists=True, mandatory=True, desc="Anatomical image to preprocess")
    model_script_path = File(exists=True, mandatory=True, desc="Anatomical image to preprocess")

class PydpyperOutputSpec(TraitedSpec):
    transform_csv = File(exists=True, desc="CSV with path info.")
    pydpyper_directory = traits.Str(exists=True, desc="Path to pydpyper directory.")

class Pydpyper(BaseInterface):

    input_spec = PydpyperInputSpec
    output_spec = PydpyperOutputSpec

    def _run_interface(self, runtime):
        import os
        cwd = os.getcwd()
        os.system('mkdir -p %s/pydpyper_run/' % (cwd,))
        os.system('bash %s %s %s/pydpyper_run/' % (self.inputs.model_script_path,self.inputs.csv_file,cwd,))
        setattr(self, 'transform_csv', '%s/pydpyper_run/transforms.csv' % (cwd,))
        setattr(self, 'pydpyper_directory', '%s/pydpyper_run/' % (cwd,))
        return runtime

    def _list_outputs(self):
        return {'transform_csv': getattr(self, 'transform_csv'),
                'pydpyper_directory': getattr(self, 'pydpyper_directory')}
