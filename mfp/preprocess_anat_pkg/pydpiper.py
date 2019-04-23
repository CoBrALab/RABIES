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
        from os.path import join as opj
        import pandas as pd
        cwd = os.getcwd()

        #verify that all required outputs are present, otherwise try to run pydpiper again
        run=True
        while(run):
            #run pydpiper
            os.system('bash %s %s' % (self.inputs.model_script_path,self.inputs.csv_file,))
            run=False

            files=pd.read_csv(self.inputs.csv_file)['file']
            for file in files:
                subject_id=os.path.basename(file).split('_ses-')[0]
                session=os.path.basename(file).split('_ses-')[1][0]
                filename_template = '%s_ses-%s' % (subject_id, session)
                labels = opj(cwd,'mbm_atlasReg_processed',filename_template+'_anat_preproc','voted.mnc')
                lsq6_mask = opj(cwd,'mbm_atlasReg_atlases','DSURQE_40micron_mask','tmp',filename_template+'_anat_preproc_I_lsq6_max_mask.mnc')
                lsq6_transform = opj(cwd,'mbm_atlasReg_processed',filename_template+'_anat_preproc','transforms',filename_template+'_anat_preproc_lsq6.xfm')
                if not os.path.isfile(labels) or not os.path.isfile(lsq6_mask) or not os.path.isfile(lsq6_transform):
                    print('MISSING OUTPUT FILES. PYDPIPER IS EXECUTED AGAIN.')
                    run=True

        setattr(self, 'transform_csv', '%s/transforms.csv' % (cwd,))
        setattr(self, 'pydpiper_directory', cwd)
        return runtime

    def _list_outputs(self):
        return {'transform_csv': getattr(self, 'transform_csv'),
                'pydpiper_directory': getattr(self, 'pydpiper_directory')}
