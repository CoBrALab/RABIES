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


def bias_correction_wf(bias_cor_script='Default', bias_reg_script='Rigid', name='bias_correction_wf'):

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['ref_EPI', 'anat', 'anat_mask']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['corrected_EPI', 'resampled_mask', 'warped_EPI']),
        name='outputnode')


    bias_correction = pe.Node(EPIBiasCorrection(bias_cor_script=bias_cor_script, reg_script=bias_reg_script), name='bias_correction')

    workflow.connect([
        (inputnode, bias_correction, [('ref_EPI', 'input_ref_EPI'),
                                      ('anat', 'anat'),
                                      ('anat_mask', 'anat_mask'),
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
    bias_cor_script = traits.Str(exists=True, mandatory=True, desc="Specifying the script to use for registration.")
    reg_script = traits.Str(exists=True, mandatory=True, desc="Specifying the script to use for registration.")

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
        import numpy as np
        import nibabel as nb
        from nibabel import processing

        subject_id=os.path.basename(self.inputs.input_ref_EPI).split('_ses-')[0]
        session=os.path.basename(self.inputs.input_ref_EPI).split('_ses-')[1][0]
        run=os.path.basename(self.inputs.input_ref_EPI).split('_run-')[1][0]
        filename_template = '%s_ses-%s_run-%s' % (subject_id, session, run)

        if self.inputs.bias_cor_script=='Default':
            import rabies
            dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
            bias_cor_script_path=dir_path+'/shell_scripts/iter_bias_cor.sh'
        else:
            '''
            For user-provided bias correction command.
            '''
            if os.path.isfile(self.inputs.bias_cor_script):
                bias_cor_script_path=self.inputs.bias_cor_script
            else:
                msg='THE BIASCOR PATH %s DOES NOT EXISTS' % self.inputs.bias_cor_script
                raise ValueError(msg)

        if self.inputs.reg_script=='Rigid':
            import rabies
            dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
            reg_script_path=dir_path+'/shell_scripts/Rigid_registration.sh'
        else:
            '''
            For user-provided reg script.
            '''
            if os.path.isfile(self.inputs.reg_script):
                reg_script_path=self.inputs.reg_script
            else:
                msg='THE BIASCOR PATH %s DOES NOT EXISTS' % self.inputs.reg_script
                raise ValueError(msg)

        cwd=os.getcwd()
        warped_image='%s/%s_output_warped_image.nii.gz' % (cwd, filename_template)
        resampled_mask='%s/%s_resampled_mask.nii.gz' % (cwd, filename_template)
        biascor_EPI='%s/%s_bias_cor.nii.gz' % (cwd, filename_template)

        #resample to isotropic resolution based on lowest dimension
        dim=nb.load(self.inputs.input_ref_EPI).header.get_zooms()
        low_dim=np.asarray(dim).min()
        processing.resample_to_output(nb.load(self.inputs.input_ref_EPI), voxel_sizes=(low_dim,low_dim,low_dim), order=4).to_filename(cwd+'/resampled.nii.gz')

        command='bash %s %s %s %s %s %s' % (bias_cor_script_path,cwd+'/resampled.nii.gz', self.inputs.anat, self.inputs.anat_mask, filename_template, reg_script_path)
        if os.system(command) != 0:
            raise ValueError('Error in '+command)

        #resample to anatomical image resolution
        dim=nb.load(self.inputs.anat).header.get_zooms()
        low_dim=np.asarray(dim).min()
        processing.resample_to_output(nb.load(cwd+'/iter_corrected.nii.gz'), voxel_sizes=(low_dim,low_dim,low_dim), order=4).to_filename(biascor_EPI)

        from .utils import resample_image
        resample_image(nb.load(biascor_EPI), os.environ["rabies_data_type"]).to_filename(biascor_EPI)
        resample_image(nb.load(warped_image), os.environ["rabies_data_type"]).to_filename(warped_image)
        resample_image(nb.load(resampled_mask), os.environ["rabies_data_type"]).to_filename(resampled_mask)

        setattr(self, 'corrected_EPI', biascor_EPI)
        setattr(self, 'warped_EPI', warped_image)
        setattr(self, 'resampled_mask', resampled_mask)

        return runtime

    def _list_outputs(self):
        return {'corrected_EPI': getattr(self, 'corrected_EPI'),
                'warped_EPI': getattr(self, 'warped_EPI'),
                'resampled_mask': getattr(self, 'resampled_mask')}
