import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def bias_correction_wf(bias_cor_method='otsu_reg', rabies_data_type=8, rabies_mem_scale=1.0, name='bias_correction_wf'):

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['ref_EPI', 'anat', 'anat_mask', 'name_source']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['corrected_EPI', 'resampled_mask', 'warped_EPI']),
        name='outputnode')

    if bias_cor_method=='otsu_reg':
        bias_correction = pe.Node(OtsuEPIBiasCorrection(rabies_data_type=rabies_data_type),
                                  name='bias_correction', mem_gb=0.3*rabies_mem_scale)

    elif bias_cor_method=='thresh_reg':
        bias_correction = pe.Node(EPIBiasCorrection(rabies_data_type=rabies_data_type),
                                  name='bias_correction', mem_gb=0.3*rabies_mem_scale)
    else:
        raise ValueError("Wrong --bias_cor_method.")


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


class OtsuEPIBiasCorrectionInputSpec(BaseInterfaceInputSpec):
    input_ref_EPI = File(exists=True, mandatory=True,
                         desc="The input 3D ref EPI to correct for bias fields")
    anat = File(exists=True, mandatory=True,
                desc="Anatomical reference image for registration")
    anat_mask = File(exists=True, mandatory=True,
                     desc="Brain mask for the anatomical image")
    name_source = File(exists=True, mandatory=True,
                       desc='Reference BOLD file for naming the output.')
    rabies_data_type = traits.Int(mandatory=True,
        desc="Integer specifying SimpleITK data type.")


class OtsuEPIBiasCorrectionOutputSpec(TraitedSpec):
    corrected_EPI = File(
        exists=True, desc="input ref EPI corrected for bias fields")
    warped_EPI = File(desc="output warped image from antsRegistration")
    resampled_mask = File(
        exists=True, desc="resampled EPI mask after registration")


class OtsuEPIBiasCorrection(BaseInterface):
    '''
    This interfaces will use multiple iterations of otsu masking combined a rigid registration
    to evaluate a final bias field correction for the 3D reference EPI.
    '''

    input_spec = OtsuEPIBiasCorrectionInputSpec
    output_spec = OtsuEPIBiasCorrectionOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import SimpleITK as sitk

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.name_source).name.rsplit(".nii")

        import rabies
        from rabies.preprocess_pkg.utils import run_command, resample_image_spacing
        dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
        reg_script_path=dir_path+'/shell_scripts/antsRegistration_rigid.sh'

        cwd = os.getcwd()
        warped_image = '%s/%s_output_warped_image.nii.gz' % (
            cwd, filename_split[0])
        resampled = '%s/%s_resampled.nii.gz' % (
            cwd, filename_split[0])
        resampled_mask = '%s/%s_resampled_mask.nii.gz' % (
            cwd, filename_split[0])
        biascor_EPI = '%s/%s_bias_cor.nii.gz' % (cwd, filename_split[0],)

        # resample to isotropic resolution based on lowest dimension
        input_ref_EPI_img = sitk.ReadImage(
            self.inputs.input_ref_EPI, self.inputs.rabies_data_type)
        dim = input_ref_EPI_img.GetSpacing()
        low_dim = np.asarray(dim).min()
        sitk.WriteImage(resample_image_spacing(
            input_ref_EPI_img, (low_dim, low_dim, low_dim)), resampled)

        # the -b will be rounded up to the nearest multiple of 10 of the image largest dimension
        largest_dim = (np.array(input_ref_EPI_img.GetSize())*np.array(input_ref_EPI_img.GetSpacing())).max()
        b_value = int(np.ceil(largest_dim/10)*10)

        bias_cor_input = resampled
        otsu_bias_cor(target=bias_cor_input, otsu_ref=bias_cor_input, out_name='corrected_iter1.nii.gz', b_value=b_value, n_iter=100)
        otsu_bias_cor(target=bias_cor_input, otsu_ref='corrected_iter1.nii.gz', out_name='corrected_iter2.nii.gz', b_value=b_value, n_iter=100)

        command = 'bash %s %s %s %s %s' % (reg_script_path, 'corrected_iter2.nii.gz', self.inputs.anat, self.inputs.anat_mask, filename_split[0],)
        rc = run_command(command)

        command = 'antsApplyTransforms -d 3 -i %s -t [%s_output_0GenericAffine.mat,1] -r %s -o %s -n GenericLabel' % (self.inputs.anat_mask,filename_split[0], 'corrected_iter2.nii.gz',resampled_mask)
        rc = run_command(command)

        otsu_bias_cor(target=bias_cor_input, otsu_ref='corrected_iter2.nii.gz', out_name=cwd+'/final_otsu.nii.gz', b_value=b_value, mask=resampled_mask, n_iter=75)

        # resample to anatomical image resolution
        dim = sitk.ReadImage(self.inputs.anat, self.inputs.rabies_data_type).GetSpacing()
        low_dim = np.asarray(dim).min()
        sitk.WriteImage(resample_image_spacing(sitk.ReadImage(cwd+'/final_otsu.nii.gz',
                                                              self.inputs.rabies_data_type), (low_dim, low_dim, low_dim)), biascor_EPI)

        sitk.WriteImage(sitk.ReadImage(biascor_EPI, self.inputs.rabies_data_type), biascor_EPI)
        sitk.WriteImage(sitk.ReadImage(warped_image, self.inputs.rabies_data_type), warped_image)
        sitk.WriteImage(sitk.ReadImage(resampled_mask, self.inputs.rabies_data_type), resampled_mask)

        setattr(self, 'corrected_EPI', biascor_EPI)
        setattr(self, 'warped_EPI', warped_image)
        setattr(self, 'resampled_mask', resampled_mask)

        return runtime

    def _list_outputs(self):
        return {'corrected_EPI': getattr(self, 'corrected_EPI'),
                'warped_EPI': getattr(self, 'warped_EPI'),
                'resampled_mask': getattr(self, 'resampled_mask')}

def otsu_bias_cor(target, otsu_ref, out_name, b_value, mask=None, n_iter=100):
    import SimpleITK as sitk
    from rabies.preprocess_pkg.utils import run_command
    command = 'ImageMath 3 null_mask.nii.gz ThresholdAtMean %s 0' % (otsu_ref)
    rc = run_command(command)
    command = 'ThresholdImage 3 %s otsu_weight.nii.gz Otsu 4' % (otsu_ref)
    rc = run_command(command)

    otsu_img = sitk.ReadImage(
        'otsu_weight.nii.gz', sitk.sitkUInt8)
    otsu_array = sitk.GetArrayFromImage(otsu_img)

    if mask is not None:
        resampled_mask_img = sitk.ReadImage(
            mask, sitk.sitkUInt8)
        resampled_mask_array = sitk.GetArrayFromImage(resampled_mask_img)

        otsu_array = otsu_array*resampled_mask_array

    combined_mask=(otsu_array==1.0)+(otsu_array==2.0)+(otsu_array==3.0)
    mask_img=sitk.GetImageFromArray(combined_mask.astype('uint8'), isVector=False)
    mask_img.CopyInformation(otsu_img)
    sitk.WriteImage(mask_img, 'mask123.nii.gz')

    combined_mask=(otsu_array==2.0)+(otsu_array==3.0)+(otsu_array==4.0)
    mask_img=sitk.GetImageFromArray(combined_mask.astype('uint8'), isVector=False)
    mask_img.CopyInformation(otsu_img)
    sitk.WriteImage(mask_img, 'mask234.nii.gz')

    combined_mask=(otsu_array==1.0)+(otsu_array==2.0)+(otsu_array==3.0)+(otsu_array==4.0)
    mask_img=sitk.GetImageFromArray(combined_mask.astype('uint8'), isVector=False)
    mask_img.CopyInformation(otsu_img)
    sitk.WriteImage(mask_img, 'mask1234.nii.gz')

    command = 'N4BiasFieldCorrection -d 3 -i %s -b %s -s 1 -c [%sx%sx%s,0.0] -w mask123.nii.gz -x null_mask.nii.gz -o corrected2.nii.gz' % (target, str(b_value), str(n_iter),str(n_iter),str(n_iter),)
    rc = run_command(command)

    command = 'N4BiasFieldCorrection -d 3 -i corrected2.nii.gz -b %s -s 1 -c [%sx%sx%s,0.0] -w mask234.nii.gz -x null_mask.nii.gz -o corrected3.nii.gz' % (str(b_value), str(n_iter),str(n_iter),str(n_iter),)
    rc = run_command(command)

    command = 'N4BiasFieldCorrection -d 3 -i corrected3.nii.gz -b %s -s 1 -c [%sx%sx%s,0.0] -w mask1234.nii.gz -x null_mask.nii.gz -o %s' % (str(b_value), str(n_iter),str(n_iter),str(n_iter),out_name,)
    rc = run_command(command)

class EPIBiasCorrectionInputSpec(BaseInterfaceInputSpec):
    input_ref_EPI = File(exists=True, mandatory=True,
                         desc="The input 3D ref EPI to correct for bias fields")
    anat = File(exists=True, mandatory=True,
                desc="Anatomical reference image for registration")
    anat_mask = File(exists=True, mandatory=True,
                     desc="Brain mask for the anatomical image")
    name_source = File(exists=True, mandatory=True,
                       desc='Reference BOLD file for naming the output.')
    rabies_data_type = traits.Int(mandatory=True,
        desc="Integer specifying SimpleITK data type.")


class EPIBiasCorrectionOutputSpec(TraitedSpec):
    corrected_EPI = File(
        exists=True, desc="input ref EPI corrected for bias fields")
    warped_EPI = File(desc="output warped image from antsRegistration")
    resampled_mask = File(
        exists=True, desc="resampled EPI mask after registration")


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
        import SimpleITK as sitk

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.name_source).name.rsplit(".nii")

        import rabies
        from rabies.preprocess_pkg.utils import run_command, resample_image_spacing
        dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
        reg_script_path=dir_path+'/shell_scripts/antsRegistration_rigid.sh'
        bias_cor_script_path = dir_path+'/shell_scripts/iter_bias_cor.sh'

        cwd = os.getcwd()
        warped_image = '%s/%s_output_warped_image.nii.gz' % (
            cwd, filename_split[0])
        resampled_mask = '%s/%s_resampled_mask.nii.gz' % (
            cwd, filename_split[0])
        biascor_EPI = '%s/%s_bias_cor.nii.gz' % (cwd, filename_split[0],)

        # resample to isotropic resolution based on lowest dimension
        input_ref_EPI = sitk.ReadImage(
            self.inputs.input_ref_EPI, self.inputs.rabies_data_type)
        dim = input_ref_EPI.GetSpacing()
        low_dim = np.asarray(dim).min()
        from rabies.preprocess_pkg.utils import resample_image_spacing
        sitk.WriteImage(resample_image_spacing(
            input_ref_EPI, (low_dim, low_dim, low_dim)), cwd+'/resampled.nii.gz')

        command = 'bash %s %s %s %s %s %s' % (bias_cor_script_path, self.inputs.input_ref_EPI,
                                              self.inputs.anat, self.inputs.anat_mask, filename_split[0], reg_script_path)
        from rabies.preprocess_pkg.utils import run_command
        rc = run_command(command)

        # resample to anatomical image resolution
        dim = sitk.ReadImage(self.inputs.anat, self.inputs.rabies_data_type).GetSpacing()
        low_dim = np.asarray(dim).min()
        sitk.WriteImage(resample_image_spacing(sitk.ReadImage(cwd+'/iter_corrected.nii.gz',
                                                              self.inputs.rabies_data_type), (low_dim, low_dim, low_dim)), biascor_EPI)

        sitk.WriteImage(sitk.ReadImage(biascor_EPI, self.inputs.rabies_data_type), biascor_EPI)
        sitk.WriteImage(sitk.ReadImage(warped_image, self.inputs.rabies_data_type), warped_image)
        sitk.WriteImage(sitk.ReadImage(resampled_mask, self.inputs.rabies_data_type), resampled_mask)

        setattr(self, 'corrected_EPI', biascor_EPI)
        setattr(self, 'warped_EPI', warped_image)
        setattr(self, 'resampled_mask', resampled_mask)

        return runtime

    def _list_outputs(self):
        return {'corrected_EPI': getattr(self, 'corrected_EPI'),
                'warped_EPI': getattr(self, 'warped_EPI'),
                'resampled_mask': getattr(self, 'resampled_mask')}
