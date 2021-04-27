from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def bias_correction_wf(opts, name='bias_correction_wf'):

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['ref_EPI', 'anat', 'anat_mask', 'name_source']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['corrected_EPI', 'denoise_mask', 'init_denoise']),
        name='outputnode')

    if opts.bold_bias_cor_method=='otsu_reg':
        bias_correction = pe.Node(OtsuEPIBiasCorrection(rabies_data_type=opts.data_type),
                                  name='bias_correction', mem_gb=0.3*opts.scale_min_memory)

    elif opts.bold_bias_cor_method=='thresh_reg':
        bias_correction = pe.Node(ThreshBiasCorrection(rabies_data_type=opts.data_type),
                                  name='bias_correction', mem_gb=0.3*opts.scale_min_memory)
    elif opts.bold_bias_cor_method=='mouse-preprocessing-v5.sh':
        from nipype.interfaces.utility import Function

        def mouse_preprocess_func(input_ref_EPI,anat,anat_mask,name_source):
            import os

            import pathlib  # Better path manipulation
            filename_split = pathlib.Path(
                name_source).name.rsplit(".nii")
            cwd = os.getcwd()
            corrected_EPI = '%s/%s_bias_cor.nii.gz' % (cwd, filename_split[0],)

            import rabies
            dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
            from rabies.preprocess_pkg.utils import run_command
            command = 'rabies-mouse-preprocessing-v5.sh %s %s %s %s' % (input_ref_EPI, corrected_EPI,anat,anat_mask)
            rc = run_command(command)
            denoise_mask = corrected_EPI.split('.nii.gz')[0]+'_mask.nii.gz'
            init_denoise = corrected_EPI.split('.nii.gz')[0]+'_init_denoise.nii.gz'
            return corrected_EPI, denoise_mask, init_denoise

        bias_correction = pe.Node(Function(input_names=['input_ref_EPI','anat','anat_mask','name_source'],
                                                    output_names=['corrected_EPI', 'denoise_mask', 'init_denoise'],
                                                    function=mouse_preprocess_func),
                                           name='bold_bias_correction')
    else:
        raise ValueError("Wrong --bold_bias_cor_method.")

    workflow.connect([
        (inputnode, bias_correction, [('ref_EPI', 'input_ref_EPI'),
                                      ('anat', 'anat'),
                                      ('anat_mask', 'anat_mask'),
                                      ('name_source', 'name_source'),
                                      ]),
        (bias_correction, outputnode, [('corrected_EPI', 'corrected_EPI'),
                                       ('denoise_mask', 'denoise_mask'),
                                       ('init_denoise', 'init_denoise'),
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
    denoise_mask = File(
        exists=True, desc="resampled mask after registration")
    init_denoise = File(
        exists=True, desc="Initial correction before registration.")


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

        from rabies.preprocess_pkg.utils import run_command, resample_image_spacing
        from rabies.preprocess_pkg.registration import run_antsRegistration

        cwd = os.getcwd()
        resampled = '%s/%s_resampled.nii.gz' % (
            cwd, filename_split[0])
        resampled_mask = '%s/%s_resampled_mask.nii.gz' % (
            cwd, filename_split[0])
        biascor_EPI = '%s/%s_bias_cor.nii.gz' % (cwd, filename_split[0],)

        input_ref_EPI_img = sitk.ReadImage(
            self.inputs.input_ref_EPI, self.inputs.rabies_data_type)

        # the -b will be rounded up to the nearest multiple of 10 of the image largest dimension
        largest_dim = (np.array(input_ref_EPI_img.GetSize())*np.array(input_ref_EPI_img.GetSpacing())).max()
        b_value = int(np.ceil(largest_dim/10)*10)

        bias_cor_input = self.inputs.input_ref_EPI
        otsu_bias_cor(target=bias_cor_input, otsu_ref=bias_cor_input, out_name=cwd+'/corrected_iter1.nii.gz', b_value=b_value)
        otsu_bias_cor(target=bias_cor_input, otsu_ref=cwd+'/corrected_iter1.nii.gz', out_name=cwd+'/corrected_iter2.nii.gz', b_value=b_value)

        [affine, warp, inverse_warp, warped_image] = run_antsRegistration(reg_method='Rigid', moving_image=cwd+'/corrected_iter2.nii.gz', fixed_image=self.inputs.anat, anat_mask=self.inputs.anat_mask)

        command = 'antsApplyTransforms -d 3 -i %s -t [%s,1] -r %s -o %s -n GenericLabel' % (self.inputs.anat_mask, affine, cwd+'/corrected_iter2.nii.gz',resampled_mask)
        rc = run_command(command)

        otsu_bias_cor(target=bias_cor_input, otsu_ref=cwd+'/corrected_iter2.nii.gz', out_name=cwd+'/final_otsu.nii.gz', b_value=b_value, mask=resampled_mask)

        # resample to anatomical image resolution
        dim = sitk.ReadImage(self.inputs.anat, self.inputs.rabies_data_type).GetSpacing()
        low_dim = np.asarray(dim).min()
        sitk.WriteImage(resample_image_spacing(sitk.ReadImage(cwd+'/final_otsu.nii.gz',
                                                              self.inputs.rabies_data_type), (low_dim, low_dim, low_dim)), biascor_EPI)

        sitk.WriteImage(sitk.ReadImage(cwd+'/corrected_iter2.nii.gz', self.inputs.rabies_data_type), cwd+'/corrected_iter2.nii.gz')
        sitk.WriteImage(sitk.ReadImage(biascor_EPI, self.inputs.rabies_data_type), biascor_EPI)
        sitk.WriteImage(sitk.ReadImage(warped_image, self.inputs.rabies_data_type), warped_image)
        sitk.WriteImage(sitk.ReadImage(resampled_mask, self.inputs.rabies_data_type), resampled_mask)

        setattr(self, 'init_denoise', cwd+'/corrected_iter2.nii.gz')
        setattr(self, 'corrected_EPI', biascor_EPI)
        setattr(self, 'warped_EPI', warped_image)
        setattr(self, 'denoise_mask', resampled_mask)

        return runtime

    def _list_outputs(self):
        return {'corrected_EPI': getattr(self, 'corrected_EPI'),
                'warped_EPI': getattr(self, 'warped_EPI'),
                'init_denoise': getattr(self, 'init_denoise'),
                'denoise_mask': getattr(self, 'denoise_mask')}

def otsu_bias_cor(target, otsu_ref, out_name, b_value, mask=None, n_iter=200):
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

    combined_mask=(otsu_array==1.0)+(otsu_array==2.0)
    mask_img=sitk.GetImageFromArray(combined_mask.astype('uint8'), isVector=False)
    mask_img.CopyInformation(otsu_img)
    sitk.WriteImage(mask_img, 'mask12.nii.gz')

    combined_mask=(otsu_array==3.0)+(otsu_array==4.0)
    mask_img=sitk.GetImageFromArray(combined_mask.astype('uint8'), isVector=False)
    mask_img.CopyInformation(otsu_img)
    sitk.WriteImage(mask_img, 'mask34.nii.gz')

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

    command = 'N4BiasFieldCorrection -d 3 -i %s -b %s -s 1 -c [%sx%sx%s,1e-4] -w mask12.nii.gz -x null_mask.nii.gz -o corrected1.nii.gz' % (target, str(b_value), str(n_iter),str(n_iter),str(n_iter),)
    rc = run_command(command)

    command = 'N4BiasFieldCorrection -d 3 -i corrected1.nii.gz -b %s -s 1 -c [%sx%sx%s,1e-4] -w mask34.nii.gz -x null_mask.nii.gz -o corrected2.nii.gz' % (str(b_value), str(n_iter),str(n_iter),str(n_iter),)
    rc = run_command(command)

    command = 'N4BiasFieldCorrection -d 3 -i corrected2.nii.gz -b %s -s 1 -c [%sx%sx%s,1e-4] -w mask123.nii.gz -x null_mask.nii.gz -o corrected3.nii.gz' % (str(b_value), str(n_iter),str(n_iter),str(n_iter),)
    rc = run_command(command)

    command = 'N4BiasFieldCorrection -d 3 -i corrected3.nii.gz -b %s -s 1 -c [%sx%sx%s,1e-4] -w mask234.nii.gz -x null_mask.nii.gz -o corrected4.nii.gz' % (str(b_value), str(n_iter),str(n_iter),str(n_iter),)
    rc = run_command(command)

    command = 'N4BiasFieldCorrection -d 3 -i corrected4.nii.gz -b %s -s 1 -c [%sx%sx%s,1e-4] -w mask1234.nii.gz -x null_mask.nii.gz -o %s' % (str(b_value), str(n_iter),str(n_iter),str(n_iter),out_name,)
    rc = run_command(command)

class ThreshBiasCorrectionInputSpec(BaseInterfaceInputSpec):
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


class ThreshBiasCorrectionOutputSpec(TraitedSpec):
    corrected_EPI = File(
        exists=True, desc="input ref EPI corrected for bias fields")
    warped_EPI = File(desc="output warped image from antsRegistration")
    denoise_mask = File(
        exists=True, desc="resampled mask after registration")
    init_denoise = File(
        exists=True, desc="Initial correction before registration.")


class ThreshBiasCorrection(BaseInterface):
    '''
    This interfaces will first apply N4BiasFieldCorrection to a 3D EPI by using an anatomical mask derived from a corresponding anatomical
    image, and will then register the image to the anatomical image to improve the localization of the mask
    for the bias field correction, which will then be applied once more with the relocated mask.
    '''

    input_spec = ThreshBiasCorrectionInputSpec
    output_spec = ThreshBiasCorrectionOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import SimpleITK as sitk

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.name_source).name.rsplit(".nii")

        from rabies.preprocess_pkg.utils import run_command, resample_image_spacing
        from rabies.preprocess_pkg.registration import run_antsRegistration

        cwd = os.getcwd()
        resampled_mask = '%s/%s_resampled_mask.nii.gz' % (
            cwd, filename_split[0])
        biascor_EPI = '%s/%s_bias_cor.nii.gz' % (cwd, filename_split[0],)

        command = 'ImageMath 3 null_mask.nii.gz ThresholdAtMean %s 0' % (self.inputs.input_ref_EPI)
        rc = run_command(command)
        command = 'ImageMath 3 thresh_mask.nii.gz ThresholdAtMean %s 2' % (self.inputs.input_ref_EPI)
        rc = run_command(command)

        command = 'N4BiasFieldCorrection -d 3 -s 4 -i %s -b 20 -s 1 -c [100x100x100x100,1e-6] -w thresh_mask.nii.gz -x null_mask.nii.gz -o corrected.nii.gz' % (self.inputs.input_ref_EPI)
        rc = run_command(command)
        command = 'DenoiseImage -d 3 -i corrected.nii.gz -o denoise.nii.gz'
        rc = run_command(command)

        [affine, warp, inverse_warp, warped_image] = run_antsRegistration(reg_method='Rigid', moving_image=cwd+'/denoise.nii.gz', fixed_image=self.inputs.anat, anat_mask=self.inputs.anat_mask)

        command = 'antsApplyTransforms -d 3 -i %s -t [%s,1] -r %s -o %s -n GenericLabel' % (self.inputs.anat_mask, affine, self.inputs.input_ref_EPI,resampled_mask)
        rc = run_command(command)

        command = 'N4BiasFieldCorrection -d 3 -s 2 -i %s -b 20 -s 1 -c [100x100x100x100,1e-6] -w %s -x null_mask.nii.gz -o %s' % (self.inputs.input_ref_EPI, resampled_mask,cwd+'/iter_corrected.nii.gz')
        rc = run_command(command)
        command = 'DenoiseImage -d 3 -i iter_corrected.nii.gz -o iter_corrected_denoise.nii.gz'
        rc = run_command(command)


        # resample to anatomical image resolution
        dim = sitk.ReadImage(self.inputs.anat, self.inputs.rabies_data_type).GetSpacing()
        low_dim = np.asarray(dim).min()
        sitk.WriteImage(resample_image_spacing(sitk.ReadImage(cwd+'/iter_corrected_denoise.nii.gz',
                                                              self.inputs.rabies_data_type), (low_dim, low_dim, low_dim)), biascor_EPI)

        sitk.WriteImage(sitk.ReadImage(biascor_EPI, self.inputs.rabies_data_type), biascor_EPI)
        sitk.WriteImage(sitk.ReadImage(warped_image, self.inputs.rabies_data_type), warped_image)
        sitk.WriteImage(sitk.ReadImage(resampled_mask, self.inputs.rabies_data_type), resampled_mask)

        setattr(self, 'corrected_EPI', biascor_EPI)
        setattr(self, 'warped_EPI', warped_image)
        setattr(self, 'denoise_mask', resampled_mask)
        setattr(self, 'init_denoise', cwd+'/denoise.nii.gz')

        return runtime

    def _list_outputs(self):
        return {'corrected_EPI': getattr(self, 'corrected_EPI'),
                'warped_EPI': getattr(self, 'warped_EPI'),
                'init_denoise': getattr(self, 'init_denoise'),
                'denoise_mask': getattr(self, 'denoise_mask')}


def init_anat_preproc_wf(opts, name='anat_preproc_wf'):
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
        fields=['anat_preproc','init_denoise', 'denoise_mask']), name='outputnode')

    anat_preproc = pe.Node(AnatPreproc(bias_cor_method=opts.anat_bias_cor_method, disable_anat_preproc=opts.disable_anat_preproc, rabies_data_type=opts.data_type),
                           name='Anat_Preproc', mem_gb=0.6*opts.scale_min_memory)

    workflow.connect([
        (inputnode, anat_preproc, [
            ("anat_file", "nii_anat"),
            ("template_anat", "template_anat"),
            ("template_mask", "template_mask"),
            ]),
        (anat_preproc, outputnode, [
            ("anat_preproc", "anat_preproc"),
            ("init_denoise", "init_denoise"),
            ("denoise_mask", "denoise_mask"),
            ]),
    ])

    return workflow


class AnatPreprocInputSpec(BaseInterfaceInputSpec):
    nii_anat = File(exists=True, mandatory=True,
                    desc="Anatomical image to preprocess")
    template_anat = File(exists=True, mandatory=True,
                         desc="anatomical template for registration.")
    template_mask = File(exists=True, mandatory=True,
                         desc="The brain mask of the anatomical template.")
    bias_cor_method = traits.Str(
        desc="Option for bias correction.")
    disable_anat_preproc = traits.Bool(
        desc="If anatomical preprocessing is disabled, then only copy the input to a new file named _preproc.nii.gz.")
    rabies_data_type = traits.Int(mandatory=True,
        desc="Integer specifying SimpleITK data type.")


class AnatPreprocOutputSpec(TraitedSpec):
    anat_preproc = File(exists=True, desc="Preprocessed anatomical image.")
    denoise_mask = File(
        exists=True, desc="resampled mask after registration")
    init_denoise = File(
        exists=True, desc="Initial correction before registration.")


class AnatPreproc(BaseInterface):

    input_spec = AnatPreprocInputSpec
    output_spec = AnatPreprocOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import SimpleITK as sitk
        from rabies.preprocess_pkg.utils import resample_image_spacing, run_command

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(self.inputs.nii_anat).name.rsplit(".nii")
        cwd = os.getcwd()
        output_anat = '%s/%s_preproc.nii.gz' % (cwd, filename_split[0],)

        # resample the anatomical image to the resolution of the provided template
        anat_image = sitk.ReadImage(
            self.inputs.nii_anat, self.inputs.rabies_data_type)
        anat_dim = anat_image.GetSpacing()

        template_image = sitk.ReadImage(
            self.inputs.template_anat, self.inputs.rabies_data_type)
        template_dim = template_image.GetSpacing()
        if not (np.array(anat_dim) == np.array(template_dim)).sum() == 3:
            import logging
            log = logging.getLogger('root')
            log.debug('Anat image will be resampled to the template resolution.')
            resampled_anat = resample_image_spacing(anat_image, template_dim)
            input_anat = '%s/%s_resampled.nii.gz' % (cwd, filename_split[0],)
            sitk.WriteImage(resampled_anat, input_anat)
        else:
            input_anat = self.inputs.nii_anat

        if self.inputs.disable_anat_preproc:
            # resample image to specified data format
            sitk.WriteImage(sitk.ReadImage(input_anat, self.inputs.rabies_data_type), output_anat)
            init_denoise=output_anat
            resampled_mask=self.inputs.template_mask
            init_denoise=cwd+'/denoise.nii.gz'
            resampled_mask=cwd+'/resampled_mask.nii.gz'
        else:
            if self.inputs.bias_cor_method=='otsu_reg':
                bias_correction = OtsuEPIBiasCorrection(input_ref_EPI = input_anat,
                    anat = self.inputs.template_anat,
                    anat_mask = self.inputs.template_mask,
                    name_source = self.inputs.nii_anat,
                    rabies_data_type=self.inputs.rabies_data_type)
                out = bias_correction.run()
                output_anat = out.outputs.corrected_EPI
                resampled_mask = out.outputs.denoise_mask
                init_denoise = out.outputs.init_denoise
            elif self.inputs.bias_cor_method=='thresh_reg':
                bias_correction = ThreshBiasCorrection(input_ref_EPI = input_anat,
                    anat = self.inputs.template_anat,
                    anat_mask = self.inputs.template_mask,
                    name_source = self.inputs.nii_anat,
                    rabies_data_type=self.inputs.rabies_data_type)
                out = bias_correction.run()
                output_anat = out.outputs.corrected_EPI
                resampled_mask = out.outputs.denoise_mask
                init_denoise = out.outputs.init_denoise

            elif self.inputs.bias_cor_method=='mouse-preprocessing-v5.sh':
                command = 'rabies-mouse-preprocessing-v5.sh %s %s %s %s' % (input_anat, output_anat, self.inputs.template_anat,self.inputs.template_mask)
                rc = run_command(command)

                resampled_mask = output_anat.split('.nii.gz')[0]+'_mask.nii.gz'
                init_denoise = output_anat.split('.nii.gz')[0]+'_init_denoise.nii.gz'

            else:
                raise ValueError("Wrong --anat_bias_cor_method.")

            # resample image to specified data format
            sitk.WriteImage(sitk.ReadImage(output_anat, self.inputs.rabies_data_type), output_anat)

        setattr(self, 'anat_preproc', output_anat)
        setattr(self, 'init_denoise', init_denoise)
        setattr(self, 'denoise_mask', resampled_mask)
        return runtime

    def _list_outputs(self):
        return {'anat_preproc': getattr(self, 'anat_preproc'),
                'init_denoise': getattr(self, 'init_denoise'),
                'denoise_mask': getattr(self, 'denoise_mask')}
