from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def init_anat_preproc_wf(disable_anat_preproc=False, name='anat_preproc_wf'):
    '''
    This workflow executes anatomical preprocessing based on anat_preproc.sh,
    which includes initial N4 bias field correction and Adaptive
    Non-Local Means Denoising (Manjon et al. 2010), followed by rigid
    registration to a template atlas to obtain brain mask to then compute an
    optimized N4 correction and denoising.
    '''

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['anat_file', 'template_anat']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['preproc_anat']), name='outputnode')

    anat_preproc = pe.Node(AnatPreproc(disable_anat_preproc=disable_anat_preproc), name='Anat_Preproc')


    workflow.connect([
        (inputnode, anat_preproc, [("anat_file", "nii_anat"), ("template_anat", "template_anat")]),
        (anat_preproc, outputnode, [("preproc_anat", "preproc_anat")]),
    ])

    return workflow


class AnatPreprocInputSpec(BaseInterfaceInputSpec):
    nii_anat = File(exists=True, mandatory=True, desc="Anatomical image to preprocess")
    template_anat = File(exists=True, mandatory=True, desc="anatomical template for registration.")
    disable_anat_preproc = traits.Bool(desc="If anatomical preprocessing is disabled, then only copy the input to a new file named _preproc.nii.gz.")

class AnatPreprocOutputSpec(TraitedSpec):
    preproc_anat = File(exists=True, desc="Preprocessed anatomical image.")

class AnatPreproc(BaseInterface):

    input_spec = AnatPreprocInputSpec
    output_spec = AnatPreprocOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import SimpleITK as sitk
        from rabies.preprocess_bold_pkg.utils import resample_image_spacing

        cwd = os.getcwd()
        out_dir='%s/anat_preproc/' % (cwd,)
        command='mkdir -p %s' % (out_dir,)
        if os.system(command) != 0:
            raise ValueError('Error in '+command)
        anat_file=os.path.basename(self.inputs.nii_anat).split('.')[0]
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_anat='%s%s_preproc.nii.gz' % (out_dir,anat_file)

        #resample the anatomical image to the resolution of the provided template
        anat_image=sitk.ReadImage(self.inputs.nii_anat, int(os.environ["rabies_data_type"]))
        anat_dim=anat_image.GetSpacing()

        template_image=sitk.ReadImage(self.inputs.template_anat, int(os.environ["rabies_data_type"]))
        template_dim=template_image.GetSpacing()
        if not (np.array(anat_dim)==np.array(template_dim)).sum()==3:
            print('Anat image will be resampled to the template resolution.')
            resampled_anat=resample_image_spacing(anat_image,template_dim)
            input_anat=out_dir+anat_file+'_resampled.nii.gz'
            sitk.WriteImage(resampled_anat, input_anat)
        else:
            input_anat=self.inputs.nii_anat

        if self.inputs.disable_anat_preproc:
            #resample image to specified data format
            sitk.WriteImage(sitk.ReadImage(input_anat, int(os.environ["rabies_data_type"])), output_anat)
        else:
            command='bash %s/../shell_scripts/anat_preproc.sh %s %s %s' % (dir_path,input_anat,self.inputs.template_anat, output_anat)
            if os.system(command) != 0:
                raise ValueError('Error in '+command)

            #resample image to specified data format
            sitk.WriteImage(sitk.ReadImage(output_anat, int(os.environ["rabies_data_type"])), output_anat)

        setattr(self, 'preproc_anat', output_anat)
        return runtime

    def _list_outputs(self):
        return {'preproc_anat': getattr(self, 'preproc_anat')}
