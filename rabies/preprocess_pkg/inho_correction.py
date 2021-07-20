from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

def init_inho_correction_wf(opts, image_type, inho_cor_method='Affine', name='inho_correction_wf'):

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['target_img', 'anat_ref', 'anat_mask', 'name_source']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['corrected', 'denoise_mask', 'init_denoise']),
        name='outputnode')

    anat_preproc = pe.Node(InhoCorrection(image_type=image_type, inho_cor_method=inho_cor_method, rabies_data_type=opts.data_type),
                           name='InhoCorrection', mem_gb=0.6*opts.scale_min_memory)

    workflow.connect([
        (inputnode, anat_preproc, [
            ("target_img", "target_img"),
            ("anat_ref", "anat_ref"),
            ("anat_mask", "anat_mask"),
            ("name_source", "name_source"),
            ]),
        (anat_preproc, outputnode, [
            ("corrected", "corrected"),
            ("init_denoise", "init_denoise"),
            ("denoise_mask", "denoise_mask"),
            ]),
    ])

    return workflow

class InhoCorrectionInputSpec(BaseInterfaceInputSpec):
    target_img = File(exists=True, mandatory=True,
                    desc="Anatomical image to preprocess")
    anat_ref = File(exists=True, mandatory=True,
                         desc="anatomical template for registration.")
    anat_mask = File(exists=True, mandatory=True,
                         desc="The brain mask of the anatomical template.")
    name_source = File(exists=True, mandatory=True,
                       desc='Reference BOLD file for naming the output.')
    image_type = traits.Str(
        desc="Between 'EPI' or 'structural'.")
    inho_cor_method = traits.Str(
        desc="Option for inhomogeneity correction.")
    rabies_data_type = traits.Int(mandatory=True,
        desc="Integer specifying SimpleITK data type.")


class InhoCorrectionOutputSpec(TraitedSpec):
    corrected = File(exists=True, desc="Preprocessed anatomical image.")
    denoise_mask = File(
        exists=True, desc="resampled mask after registration")
    init_denoise = File(
        exists=True, desc="Initial correction before registration.")


class InhoCorrection(BaseInterface):

    input_spec = InhoCorrectionInputSpec
    output_spec = InhoCorrectionOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import SimpleITK as sitk
        from rabies.preprocess_pkg.utils import resample_image_spacing, run_command

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.name_source).name.rsplit(".nii")
        cwd = os.getcwd()

        # resample the anatomical image to the resolution of the provided template
        target_img = sitk.ReadImage(
            self.inputs.target_img, self.inputs.rabies_data_type)
        anat_dim = target_img.GetSpacing()

        template_image = sitk.ReadImage(
            self.inputs.anat_ref, self.inputs.rabies_data_type)
        template_dim = template_image.GetSpacing()
        if not (np.array(anat_dim) == np.array(template_dim)).sum() == 3:
            import logging
            log = logging.getLogger('root')
            log.debug('Anat image will be resampled to the template resolution.')
            resampled = resample_image_spacing(target_img, template_dim)
            target_img = f'{cwd}/{filename_split[0]}_resampled.nii.gz'
            sitk.WriteImage(resampled, target_img)
        else:
            target_img = self.inputs.target_img

        if self.inputs.inho_cor_method=='disable':
            # outputs correspond to the inputs
            corrected=target_img
            init_denoise=corrected
            resampled_mask=self.inputs.anat_mask
        elif self.inputs.inho_cor_method in ['Rigid','Affine','SyN']:
            corrected = f'{cwd}/{filename_split[0]}_inho_cor.nii.gz'
            if self.inputs.image_type=='EPI':
                processing_script='EPI-preprocessing.sh'
            elif self.inputs.image_type=='structural':
                processing_script='structural-preprocessing.sh'
            else:
                raise ValueError(f"Image type must be 'EPI' or 'structural', {self.inputs.image_type}")
            command = f'{processing_script} {target_img} {corrected} {self.inputs.anat_ref} {self.inputs.anat_mask} {self.inputs.inho_cor_method}'
            rc = run_command(command)

            resampled_mask = corrected.split('.nii.gz')[0]+'_mask.nii.gz'
            init_denoise = corrected.split('.nii.gz')[0]+'_init_denoise.nii.gz'
        else:
            raise ValueError("Wrong inho_cor_method.")

        # resample image to specified data format
        sitk.WriteImage(sitk.ReadImage(corrected, self.inputs.rabies_data_type), corrected)

        setattr(self, 'corrected', corrected)
        setattr(self, 'init_denoise', init_denoise)
        setattr(self, 'denoise_mask', resampled_mask)
        return runtime

    def _list_outputs(self):
        return {'corrected': getattr(self, 'corrected'),
                'init_denoise': getattr(self, 'init_denoise'),
                'denoise_mask': getattr(self, 'denoise_mask')}
