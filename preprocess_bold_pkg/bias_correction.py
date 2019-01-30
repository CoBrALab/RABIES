import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from nipype.interfaces.ants import N4BiasFieldCorrection
from nipype.interfaces.ants.resampling import ApplyTransforms


def bias_correction_wf(iterative=False, name='bias_correction_wf'):

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['ref_EPI', 'anat', 'anat_mask']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['corrected_EPI', 'resampled_mask', 'warped_EPI']),
        name='outputnode')


    bias_correction = pe.Node(EPIBiasCorrection(use_thresh_mask=True, iterative_registration=iterative), name='bias_correction')

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
    iterative_registration = traits.Bool(mandatory=True, desc="Use the iterative algorithm with registration and brain masks")
    use_thresh_mask = traits.Bool(mandatory=True, desc="Use threshold masks to mask the brain for correction")

class EPIBiasCorrectionOutputSpec(TraitedSpec):
    corrected_EPI = File(exists=True, desc="input ref EPI corrected for bias fields")
    warped_EPI = File(exists=True, desc="output warped image from antsRegistration")
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
        from nipype.interfaces.base import CommandLine
        from nipype.interfaces.ants import N4BiasFieldCorrection
        from .registration import setup_antsCoRegistration
        from nipype.interfaces.ants.resampling import ApplyTransforms

        null_mask = os.path.abspath('tmp/null_mask.nii.gz')
        thresh_mask = os.path.abspath('tmp/thresh_mask.nii.gz')
        resample_EPI = os.path.abspath('tmp/resampled.nii.gz')
        resample_100iso_EPI = os.path.abspath('tmp/resampled_100iso.nii.gz')
        n4_corrected = os.path.abspath('tmp/n4_corrected.nii.gz')
        resampled_mask = os.path.abspath('iteration2/EPIMask_resample.nii.gz')

        if self.inputs.use_thresh_mask:
            os.makedirs('tmp', exist_ok=True)

            resample = CommandLine('ResampleImage', args='3 ' + self.inputs.input_ref_EPI + ' ' + resample_EPI + ' 0.4x0.4x0.4 [BSpline]')
            resample.run()

            gen_null_mask = CommandLine('ImageMath', args='3 ' + null_mask + ' ThresholdAtMean ' + resample_EPI + ' 0', terminal_output='stream')
            gen_null_mask.run()

            #thresholding at 2%
            gen_thresh_mask = CommandLine('ImageMath', args='3 ' + thresh_mask + ' ThresholdAtMean ' + resample_EPI + ' 2', terminal_output='stream')
            gen_thresh_mask.run()

            n4_correct = N4BiasFieldCorrection(dimension=3, copy_header=True, bspline_fitting_distance=20,
                                n_iterations=[200,200,200,200], convergence_threshold=1e-6 ,shrink_factor=1, weight_image=thresh_mask, mask_image=null_mask,
                                output_image=n4_corrected)
            n4_correct.inputs.input_image=resample_EPI
            n4_res = n4_correct.run()
            print('Executed bias field correction with threshold mask.')

        else:
            os.makedirs('tmp', exist_ok=True)

            resample = CommandLine('ResampleImage', args='3 ' + self.inputs.input_ref_EPI + ' tmp/resample.nii.gz 0.4x0.4x0.4 [BSpline]')
            resample.run()

            n4_correct = N4BiasFieldCorrection(dimension=3, copy_header=True, bspline_fitting_distance=20,
                                n_iterations=[200,200,200,200], convergence_threshold=1e-6 ,shrink_factor=1,
                                output_image=n4_corrected)
            n4_correct.inputs.input_image=resample_EPI
            n4_res = n4_correct.run()

        if self.inputs.iterative_registration:

            resample = CommandLine('ResampleImage', args='3 ' + n4_corrected + ' ' + resample_100iso_EPI + ' 0.1x0.1x0.1 [BSpline]')
            resample.run()

            reg = setup_antsCoRegistration(reg_type='Rigid', moving_image=resample_100iso_EPI, fixed_image=self.inputs.anat)
            reg.inputs.moving_image_masks = ['NULL']
            reg.inputs.fixed_image_masks = ['NULL']
            res=reg.run()

            warped_EPI=res.outputs.warped_image

            os.makedirs('iteration2', exist_ok=True)

            trans = ApplyTransforms(dimension=3, input_image=self.inputs.anat_mask, transforms=res.outputs.inverse_composite_transform, reference_image=resample_EPI, output_image=resampled_mask)
            trans.run()

            gen_mask = CommandLine('ImageMath', args='3 iteration2/null_mask.nii.gz ThresholdAtMean ' + resample_EPI + ' 0', terminal_output='stream')
            gen_mask.run()

            n4_correct = N4BiasFieldCorrection(dimension=3, copy_header=True, bspline_fitting_distance=20,
                                n_iterations=[200,200,200,200], convergence_threshold=1e-6 ,shrink_factor=1, weight_image=resampled_mask, mask_image='iteration2/null_mask.nii.gz',
                                output_image='iteration2/corrected.nii.gz')
            n4_correct.inputs.input_image=resample_EPI
            n4_res = n4_correct.run()
            print('Executed bias field correction through registration with the structural image.')

        resample = CommandLine('ResampleImage', args='3 ' + n4_res.outputs.output_image + ' ' + resample_100iso_EPI + ' 0.1x0.1x0.1 [BSpline]')
        resample.run()

        setattr(self, 'corrected_EPI', resample_100iso_EPI)
        setattr(self, 'warped_EPI', warped_EPI)
        setattr(self, 'resampled_mask', resampled_mask)

        return runtime

    def _list_outputs(self):
        return {'corrected_EPI': getattr(self, 'corrected_EPI'),
                'warped_EPI': getattr(self, 'warped_EPI'),
                'resampled_mask': getattr(self, 'resampled_mask')}
