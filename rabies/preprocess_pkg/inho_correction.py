from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

def init_inho_correction_wf(opts, image_type, name='inho_correction_wf'):
    """
    Corrects an input 3D image for intensity inhomogeneities. The image is denoised with non-local mean 
    denoising (Manjón et al., 2010) followed by iterative correction for intensity inhomogeneities (Sled 
    et al., 1998). Initial masking is achieved via intensity thresholding, giving an initial correction of 
    the image, and a registration is then conducted to register a brain mask for a final round of correction.

    References:
        Manjón, J. V., Coupé, P., Martí-Bonmatí, L., Collins, D. L., & Robles, M. (2010). Adaptive non-local means 
            denoising of MR images with spatially varying noise levels. Journal of Magnetic Resonance Imaging: 
            JMRI, 31(1), 192–203.
        Sled, J. G., Zijdenbos, A. P., & Evans, A. C. (1998). A nonparametric method for automatic correction of 
            intensity nonuniformity in MRI data. IEEE Transactions on Medical Imaging, 17(1), 87–97.            

    Command line interface parameters:
        --anat_inho_cor_method {Rigid,Affine,SyN,no_reg,N4_reg,disable}
                                Select a registration type for masking during inhomogeneity correction of the structural 
                                image. 
                                *** N4_reg: previous correction script prior to version 0.3.1.
                                *** disable: disables the inhomogeneity correction.
                                (default: SyN)
                                
        --bold_inho_cor_method {Rigid,Affine,SyN,no_reg,N4_reg,disable}
                                Select a registration type for masking during inhomogeneity correction of the EPI.
                                *** N4_reg: previous correction script prior to version 0.3.1.
                                *** disable: disables the inhomogeneity correction.
                                (default: Rigid)
                                
        --bold_inho_cor_otsu BOLD_INHO_COR_OTSU
                                The inhomogeneity correction script necessitates an initial correction with a Otsu
                                masking strategy (prior to registration of an anatomical mask). This option sets the 
                                Otsu threshold level to capture the right intensity distribution.
                                (default: 2)

    Workflow:
        parameters
            opts: command line interface parameters
            image_type: between 'EPI' and 'structural'. Defines which script to run depending on 
                image type

        inputs
            target_img: the image to correct
            anat_ref: the registration target with a brain mask
            anat_mask: the brain mask of the registration target
            name_source: reference file for naming purpose

        outputs
            corrected: the output image after the final correction
            denoise_mask: the brain mask resampled on the corrected image
            init_denoise: the image after a first round of correction
    """

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['target_img', 'anat_ref', 'anat_mask', 'name_source']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['corrected', 'denoise_mask', 'init_denoise']),
        name='outputnode')

    if image_type=='EPI':
        inho_cor_method=opts.bold_inho_cor_method
    elif image_type=='structural':
        inho_cor_method=opts.anat_inho_cor_method
    else:
        raise
    anat_preproc = pe.Node(InhoCorrection(image_type=image_type, inho_cor_method=inho_cor_method, otsu_threshold=opts.bold_inho_cor_otsu, rabies_data_type=opts.data_type),
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
    otsu_threshold = traits.Int(
        desc="")
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
        from rabies.utils import resample_image_spacing, run_command

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.name_source).name.rsplit(".nii")
        cwd = os.getcwd()

        # resample the anatomical image to the resolution of the provided template
        target_img = sitk.ReadImage(
            self.inputs.target_img, self.inputs.rabies_data_type)
        if not target_img.GetDimension()==3:
            raise ValueError(f"Input image {self.inputs.target_img} is not 3-dimensional.")
        anat_dim = target_img.GetSpacing()

        template_image = sitk.ReadImage(
            self.inputs.anat_ref, self.inputs.rabies_data_type)
        template_dim = template_image.GetSpacing()
        if not (np.array(anat_dim) == np.array(template_dim)).sum() == 3:
            from nipype import logging
            log = logging.getLogger('nipype.workflow')
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
        elif self.inputs.inho_cor_method in ['Rigid','Affine','SyN', 'no_reg']:
            corrected = f'{cwd}/{filename_split[0]}_inho_cor.nii.gz'
            if self.inputs.image_type=='EPI':
                processing_script='EPI-preprocessing.sh'
            elif self.inputs.image_type=='structural':
                processing_script='structural-preprocessing.sh'
            else:
                raise ValueError(f"Image type must be 'EPI' or 'structural', {self.inputs.image_type}")
            command = f'{processing_script} {target_img} {corrected} {self.inputs.anat_ref} {self.inputs.anat_mask} {self.inputs.inho_cor_method} {str(self.inputs.otsu_threshold)}'
            rc = run_command(command)

            resampled_mask = corrected.split('.nii.gz')[0]+'_mask.nii.gz'
            init_denoise = corrected.split('.nii.gz')[0]+'_init_denoise.nii.gz'
        elif self.inputs.inho_cor_method in ['N4_reg']:
            if self.inputs.image_type=='EPI':

                bias_correction = OtsuEPIBiasCorrection(
                    input_ref_EPI=target_img, anat=self.inputs.anat_ref, anat_mask=self.inputs.anat_mask, 
                    name_source=self.inputs.name_source, rabies_data_type=self.inputs.rabies_data_type)
                out = bias_correction.run()
                corrected = out.outputs.corrected_EPI
                resampled_mask = out.outputs.denoise_mask
                init_denoise = out.outputs.init_denoise

            elif self.inputs.image_type=='structural':
                corrected = f'{cwd}/{filename_split[0]}_inho_cor.nii.gz'
                input_anat = target_img

                command = 'ImageMath 3 null_mask.nii.gz ThresholdAtMean %s 0' % (input_anat)
                rc = run_command(command)
                command = 'ImageMath 3 thresh_mask.nii.gz ThresholdAtMean %s 1.2' % (input_anat)
                rc = run_command(command)

                command = 'N4BiasFieldCorrection -d 3 -s 4 -i %s -b [20] -c [200x200x200,0.0] -w thresh_mask.nii.gz -x null_mask.nii.gz -o N4.nii.gz' % (input_anat)
                rc = run_command(command)
                command = 'DenoiseImage -d 3 -i N4.nii.gz -o denoise.nii.gz'
                rc = run_command(command)

                from rabies.preprocess_pkg.registration import run_antsRegistration
                [affine, warp, inverse_warp, warped_image] = run_antsRegistration(reg_method='Affine', moving_image=input_anat, fixed_image=self.inputs.anat_ref, fixed_mask=self.inputs.anat_mask)

                command = 'antsApplyTransforms -d 3 -i %s -t [%s,1] -r %s -o resampled_mask.nii.gz -n GenericLabel' % (self.inputs.anat_mask, affine, input_anat)
                rc = run_command(command)

                command = 'N4BiasFieldCorrection -d 3 -s 2 -i %s -b [20] -c [200x200x200x200,0.0] -w resampled_mask.nii.gz -r 1 -x null_mask.nii.gz -o N4.nii.gz' % (input_anat)
                rc = run_command(command)
                command = 'DenoiseImage -d 3 -i N4.nii.gz -o %s' % (corrected)
                rc = run_command(command)

                # resample image to specified data format
                sitk.WriteImage(sitk.ReadImage(corrected, self.inputs.rabies_data_type), corrected)
                init_denoise=cwd+'/denoise.nii.gz'
                resampled_mask=cwd+'/resampled_mask.nii.gz'

            else:
                raise ValueError(f"Image type must be 'EPI' or 'structural', {self.inputs.image_type}")

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

        from rabies.utils import run_command, resample_image_spacing
        from rabies.preprocess_pkg.registration import run_antsRegistration

        cwd = os.getcwd()
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
        #sitk.WriteImage(resample_image_spacing(
        #    input_ref_EPI_img, (low_dim, low_dim, low_dim)), resampled)

        # the -b will be rounded up to the nearest multiple of 10 of the image largest dimension
        largest_dim = (np.array(input_ref_EPI_img.GetSize())*np.array(input_ref_EPI_img.GetSpacing())).max()
        b_value = int(np.ceil(largest_dim/10)*10)

        bias_cor_input = self.inputs.input_ref_EPI
        otsu_bias_cor(target=bias_cor_input, otsu_ref=bias_cor_input, out_name=cwd+'/corrected_iter1.nii.gz', b_value=b_value)
        otsu_bias_cor(target=bias_cor_input, otsu_ref=cwd+'/corrected_iter1.nii.gz', out_name=cwd+'/corrected_iter2.nii.gz', b_value=b_value)

        [affine, warp, inverse_warp, warped_image] = run_antsRegistration(reg_method='Rigid', moving_image=cwd+'/corrected_iter2.nii.gz', fixed_image=self.inputs.anat, fixed_mask=self.inputs.anat_mask)

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
    from rabies.utils import run_command
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
