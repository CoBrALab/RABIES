from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def init_anat_preproc_wf(reg_script, disable_anat_preproc=False, rabies_data_type=8, rabies_mem_scale=1.0, name='anat_preproc_wf'):
    '''
    This workflow executes anatomical preprocessing based on anat_preproc.sh,
    which includes initial N4 bias field correction and Adaptive
    Non-Local Means Denoising (Manjon et al. 2010), followed by rigid
    registration to a template atlas to obtain brain mask to then compute an
    optimized N4 correction and denoising.
    '''

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['anat_file', 'template_anat', 'template_mask']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['preproc_anat']), name='outputnode')

    anat_preproc = pe.Node(AnatPreproc(reg_script=reg_script, disable_anat_preproc=disable_anat_preproc, rabies_data_type=rabies_data_type),
                           name='Anat_Preproc', mem_gb=0.6*rabies_mem_scale)

    workflow.connect([
        (inputnode, anat_preproc, [
            ("anat_file", "nii_anat"),
            ("template_anat", "template_anat"),
            ("template_mask", "template_mask"),
            ]),
        (anat_preproc, outputnode, [("preproc_anat", "preproc_anat")]),
    ])

    return workflow


class AnatPreprocInputSpec(BaseInterfaceInputSpec):
    nii_anat = File(exists=True, mandatory=True,
                    desc="Anatomical image to preprocess")
    template_anat = File(exists=True, mandatory=True,
                         desc="anatomical template for registration.")
    template_mask = File(exists=True, mandatory=True,
                         desc="The brain mask of the anatomical template.")
    disable_anat_preproc = traits.Bool(
        desc="If anatomical preprocessing is disabled, then only copy the input to a new file named _preproc.nii.gz.")
    reg_script = traits.Str(exists=True, mandatory=True,
                            desc="Specifying the script to use for registration.")
    rabies_data_type = traits.Int(mandatory=True,
        desc="Integer specifying SimpleITK data type.")


class AnatPreprocOutputSpec(TraitedSpec):
    preproc_anat = File(exists=True, desc="Preprocessed anatomical image.")


class AnatPreproc(BaseInterface):

    input_spec = AnatPreprocInputSpec
    output_spec = AnatPreprocOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import SimpleITK as sitk
        from rabies.preprocess_pkg.utils import resample_image_spacing, run_command

        cwd = os.getcwd()
        out_dir = '%s/anat_preproc/' % (cwd,)
        command = 'mkdir -p %s' % (out_dir,)
        rc = run_command(command)

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(self.inputs.nii_anat).name.rsplit(".nii")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_anat = '%s%s_preproc.nii.gz' % (out_dir, filename_split[0],)

        # resample the anatomical image to the resolution of the provided template
        anat_image = sitk.ReadImage(
            self.inputs.nii_anat, self.inputs.rabies_data_type)
        anat_dim = anat_image.GetSpacing()

        template_image = sitk.ReadImage(
            self.inputs.template_anat, self.inputs.rabies_data_type)
        template_dim = template_image.GetSpacing()
        if not (np.array(anat_dim) == np.array(template_dim)).sum() == 3:
            print('Anat image will be resampled to the template resolution.')
            resampled_anat = resample_image_spacing(anat_image, template_dim)
            input_anat = out_dir+filename_split[0]+'_resampled.nii.gz'
            sitk.WriteImage(resampled_anat, input_anat)
        else:
            input_anat = self.inputs.nii_anat

        if self.inputs.disable_anat_preproc:
            # resample image to specified data format
            sitk.WriteImage(sitk.ReadImage(input_anat, self.inputs.rabies_data_type), output_anat)
        else:
            command = 'bash %s/../shell_scripts/anat_preproc.sh %s %s %s %s %s' % (
                dir_path, input_anat, self.inputs.template_anat, self.inputs.template_mask, output_anat, self.inputs.reg_script)
            rc = run_command(command)

            # resample image to specified data format
            sitk.WriteImage(sitk.ReadImage(output_anat, self.inputs.rabies_data_type), output_anat)

        setattr(self, 'preproc_anat', output_anat)
        return runtime

    def _list_outputs(self):
        return {'preproc_anat': getattr(self, 'preproc_anat')}
