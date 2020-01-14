from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def init_anat_preproc_wf(name='anat_preproc_wf'):
    '''
    This workflow executes anatomical preprocessing based on anat_preproc.sh,
    which includes initial N4 bias field correction and Adaptive
    Non-Local Means Denoising (Manjon et al. 2010), followed by rigid
    registration to a template atlas to obtain brain mask to then compute an
    optimized N4 correction and denoising.
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
        import numpy as np
        import nibabel as nb
        from nibabel import processing
        from rabies.preprocess_bold_pkg.utils import resample_image

        cwd = os.getcwd()
        out_dir='%s/anat_preproc/' % (cwd,)
        os.system('mkdir -p %s' % (out_dir,))
        anat_file=os.path.basename(self.inputs.nii_anat).split('.')[0]

        #resample the anatomical image to isotropic rez if it is not already the case to facilitate registration
        dim=nb.load(self.inputs.nii_anat).header.get_zooms()
        low_dim=np.asarray(dim).min()
        if not (dim==low_dim).sum()==3:
            print('ANAT IMAGE NOT ISOTROPIC. WILL RESAMPLE TO ISOTROPIC RESOLUTION BASED ON LOWEST AVAILABLE RESOLUTION.')
            input_anat=out_dir+anat_file+'_resampled.nii.gz'
            processing.resample_to_output(nb.load(self.inputs.nii_anat), voxel_sizes=(low_dim,low_dim,low_dim), order=4).to_filename(input_anat)
        else:
            input_anat=self.inputs.nii_anat

        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_anat='%s%s_preproc.nii.gz' % (out_dir,anat_file)
        os.system('bash %s/../shell_scripts/anat_preproc.sh %s %s' % (dir_path,input_anat,output_anat))

        resample_image(nb.load(output_anat), os.environ["rabies_data_type"]).to_filename(output_anat)

        setattr(self, 'preproc_anat', output_anat)
        return runtime

    def _list_outputs(self):
        return {'preproc_anat': getattr(self, 'preproc_anat')}
