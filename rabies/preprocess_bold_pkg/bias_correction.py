import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from nipype.interfaces.base import CommandLine
from nipype.interfaces.ants import N4BiasFieldCorrection
from nipype.interfaces.ants.resampling import ApplyTransforms


def bias_correction_wf(bias_reg_script='Rigid', name='bias_correction_wf'):

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['ref_EPI', 'anat', 'anat_mask', 'name_source']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['corrected_EPI', 'resampled_mask', 'warped_EPI']),
        name='outputnode')


    bias_correction = pe.Node(EPIBiasCorrection(reg_script=bias_reg_script), name='bias_correction', mem_gb=0.3*float(os.environ["rabies_mem_scale"]))

    workflow.connect([
        (inputnode, bias_correction, [('ref_EPI', 'input_ref_EPI'),
                                      ('anat', 'anat'),
                                      ('anat_mask', 'anat_mask'),
                                      ('name_source', 'name_source'),
                                      ]),
        (bias_correction, outputnode, [('corrected_EPI', 'corrected_EPI'),
                                      ('warped_EPI', 'warped_EPI'),
                                      ('resampled_mask', 'resampled_mask')
                                      ]),
    ])


    return workflow



class EPIBiasCorrectionInputSpec(BaseInterfaceInputSpec):
    input_ref_EPI = File(exists=True, mandatory=True, desc="The input 3D ref EPI to correct for bias fields")
    anat = File(exists=True, mandatory=True, desc="Anatomical reference image for registration")
    anat_mask = File(exists=True, mandatory=True, desc="Brain mask for the anatomical image")
    reg_script = traits.Str(exists=True, mandatory=True, desc="Specifying the script to use for registration.")
    name_source = File(exists=True, mandatory=True, desc='Reference BOLD file for naming the output.')

class EPIBiasCorrectionOutputSpec(TraitedSpec):
    corrected_EPI = File(exists=True, desc="input ref EPI corrected for bias fields")
    warped_EPI = File(desc="output warped image from antsRegistration")
    resampled_mask = File(exists=True, desc="resampled EPI mask after registration")


class EPIBiasCorrection(BaseInterface):
    '''
    This interfaces will first apply N4BiasFieldCorrection to a 3D EPI by using an anatomical mask derived from a corresponding anatomical
    image, and will then register the image to the anatomical image to improve the localization of the mask
    for the bias field correction, which will then be applied once more with the relocated mask.
    '''

    input_spec = EPIBiasCorrectionInputSpec
    output_spec = EPIBiasCorrectionOutputSpec

    def _run_interface(self, runtime):
        import os
        import subprocess
        import numpy as np
        import SimpleITK as sitk

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(self.inputs.name_source).name.rsplit(".nii")

        import rabies
        dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
        bias_cor_script_path=dir_path+'/shell_scripts/iter_bias_cor.sh'

        if os.path.isfile(self.inputs.reg_script):
            reg_script_path=self.inputs.reg_script
        else:
            msg='THE BIASCOR PATH %s DOES NOT EXISTS' % self.inputs.reg_script
            raise ValueError(msg)

        cwd=os.getcwd()
        warped_image='%s/%s_output_warped_image.nii.gz' % (cwd, filename_split[0])
        resampled_mask='%s/%s_resampled_mask.nii.gz' % (cwd, filename_split[0])
        biascor_EPI='%s/%s_bias_cor.nii.gz' % (cwd, filename_split[0],)

        #resample to isotropic resolution based on lowest dimension
        input_ref_EPI=sitk.ReadImage(self.inputs.input_ref_EPI, int(os.environ["rabies_data_type"]))
        dim=input_ref_EPI.GetSpacing()
        low_dim=np.asarray(dim).min()
        from rabies.preprocess_bold_pkg.utils import resample_image_spacing
        sitk.WriteImage(resample_image_spacing(input_ref_EPI, (low_dim,low_dim,low_dim)), cwd+'/resampled.nii.gz')

        command='bash %s %s %s %s %s %s' % (bias_cor_script_path,cwd+'/resampled.nii.gz', self.inputs.anat, self.inputs.anat_mask, filename_split[0], reg_script_path)
        from rabies.preprocess_bold_pkg.utils import run_command
        rc = run_command(command)

        #resample to anatomical image resolution
        dim=sitk.ReadImage(self.inputs.anat, int(os.environ["rabies_data_type"])).GetSpacing()
        low_dim=np.asarray(dim).min()
        sitk.WriteImage(resample_image_spacing(sitk.ReadImage(cwd+'/iter_corrected.nii.gz', int(os.environ["rabies_data_type"])), (low_dim,low_dim,low_dim)), biascor_EPI)

        sitk.WriteImage(sitk.ReadImage(biascor_EPI, int(os.environ["rabies_data_type"])), biascor_EPI)
        sitk.WriteImage(sitk.ReadImage(warped_image, int(os.environ["rabies_data_type"])), warped_image)
        sitk.WriteImage(sitk.ReadImage(resampled_mask, int(os.environ["rabies_data_type"])), resampled_mask)

        setattr(self, 'corrected_EPI', biascor_EPI)
        setattr(self, 'warped_EPI', warped_image)
        setattr(self, 'resampled_mask', resampled_mask)

        return runtime

    def _list_outputs(self):
        return {'corrected_EPI': getattr(self, 'corrected_EPI'),
                'warped_EPI': getattr(self, 'warped_EPI'),
                'resampled_mask': getattr(self, 'resampled_mask')}
