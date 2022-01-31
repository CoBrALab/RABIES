import os
import pathlib
import SimpleITK as sitk
import numpy as np
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from rabies.utils import copyInfo_4DImage, copyInfo_3DImage, run_command
from .hmc import antsMotionCorr

def init_bold_reference_wf(opts, name='gen_bold_ref'):
    """
    The 4D raw EPI file is used to generate a representative volumetric 3D EPI. This volume later becomes the target for 
    motion realignment and the estimation of susceptibility distortions through registration to the structural image. 
    Two iterations of motion realignment to an initial median of the volumes are conducted, then a trimmed mean is 
    computed on the realignment volumes, ignoring 5% extreme, and this average becomes the reference image. The final
    image is then corrected using non-local means denoising (Manjón et al., 2010).

    References:
        Manjón, J. V., Coupé, P., Martí-Bonmatí, L., Collins, D. L., & Robles, M. (2010). Adaptive non-local means 
            denoising of MR images with spatially varying noise levels. Journal of Magnetic Resonance Imaging: 
            JMRI, 31(1), 192–203.

    Command line interface parameters:
        --detect_dummy        Detect and remove initial dummy volumes from the EPI, and generate a reference EPI based on
                            these volumes if detected. Dummy volumes will be removed from the output preprocessed EPI.
                            (default: False)

    Workflow:
        parameters
            opts: command line interface parameters

        inputs
            bold_file: Nifti file with EPI timeseries

        outputs
            ref_image: the reference EPI volume
            bold_file: the input EPI timeseries, but after removing dummy volumes if --detect_dummy is selected
    """

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['bold_file', 'ref_image']),
        name='outputnode')

    gen_ref = pe.Node(EstimateReferenceImage(HMC_option=opts.HMC_option, detect_dummy=opts.detect_dummy, rabies_data_type=opts.data_type),
                      name='gen_ref', mem_gb=2*opts.scale_min_memory)
    gen_ref.plugin_args = {
        'qsub_args': f'-pe smp {str(2*opts.min_proc)}', 'overwrite': True}

    workflow.connect([
        (inputnode, gen_ref, [('bold_file', 'in_file')]),
        (gen_ref, outputnode, [('ref_image', 'ref_image'),
                               ('bold_file', 'bold_file')]),
    ])

    return workflow


class EstimateReferenceImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="4D EPI file")
    HMC_option = traits.Str(desc="Select one of the pre-built options from https://github.com/ANTsX/ANTsR/blob/60eefd96fedd16bceb4703ecd2cd5730e6843807/R/ants_motion_estimation.R")
    detect_dummy = traits.Bool(
        desc="specify if should detect and remove dummy scans, and use these volumes as reference image.")
    rabies_data_type = traits.Int(mandatory=True,
                                  desc="Integer specifying SimpleITK data type.")


class EstimateReferenceImageOutputSpec(TraitedSpec):
    ref_image = File(exists=True, desc="3D reference image")
    bold_file = File(
        exists=True, desc="Input bold file without dummy volumes if detect_dummy is True.")


class EstimateReferenceImage(BaseInterface):
    """
    Given a 4D EPI file, estimates an optimal reference image that could be later
    used for motion estimation and coregistration purposes. If the detect_dummy
    option is selected, it will use detected anat saturated volumes (non-steady
    state). Otherwise, a median of a subset of motion corrected volumes is used.
    In the later case, a first median is extracted from the raw data and used as
    reference for motion correction, then a new median image is extracted from
    the corrected series, and the process is repeated one more time to generate
    a final image reference image.
    The final image is corrected using non-local means denoising. (Manjon et al. Journal of Magnetic Resonance Imaging, June 2010.)
    """

    input_spec = EstimateReferenceImageInputSpec
    output_spec = EstimateReferenceImageOutputSpec

    def _run_interface(self, runtime):

        from nipype import logging
        log = logging.getLogger('nipype.workflow')

        in_nii = sitk.ReadImage(self.inputs.in_file,
                                self.inputs.rabies_data_type)
        if not in_nii.GetDimension()==4:
            raise ValueError(f"Input image {self.inputs.in_file} is not 4-dimensional.")

        data_slice = sitk.GetArrayFromImage(in_nii)[:50, :, :, :]

        n_volumes_to_discard = _get_vols_to_discard(in_nii)

        filename_split = pathlib.Path(self.inputs.in_file).name.rsplit(".nii")
        out_ref_fname = os.path.abspath(
            f'{filename_split[0]}_bold_ref.nii.gz')

        if (not n_volumes_to_discard == 0) and self.inputs.detect_dummy:
            log.info("Detected "+str(n_volumes_to_discard)
                  + " dummy scans. Taking the median of these volumes as reference EPI.")
            median_image_data = np.median(
                data_slice[:n_volumes_to_discard, :, :, :], axis=0)

            out_bold_file = os.path.abspath(
                f'{filename_split[0]}_cropped_dummy.nii.gz')
            img_array = sitk.GetArrayFromImage(
                in_nii)[n_volumes_to_discard:, :, :, :]

            image_4d = copyInfo_4DImage(sitk.GetImageFromArray(
                img_array, isVector=False), in_nii, in_nii)
            sitk.WriteImage(image_4d, out_bold_file)

        else:
            out_bold_file = self.inputs.in_file

            n_volumes_to_discard = 0
            if self.inputs.detect_dummy:
                log.info(
                    "Detected no dummy scans. Generating the ref EPI based on multiple volumes.")
            # if no dummy scans, will generate a median from a subset of max 100
            # slices of the time series
            if in_nii.GetSize()[-1] > 100:
                slice_fname = os.path.abspath("slice.nii.gz")
                image_4d = copyInfo_4DImage(sitk.GetImageFromArray(
                    data_slice[20:100, :, :, :], isVector=False), in_nii, in_nii)
                sitk.WriteImage(image_4d, slice_fname)
                median_fname = os.path.abspath("median.nii.gz")
                image_3d = copyInfo_3DImage(sitk.GetImageFromArray(
                    np.median(data_slice[20:100, :, :, :], axis=0), isVector=False), in_nii)
                sitk.WriteImage(image_3d, median_fname)
            else:
                slice_fname = self.inputs.in_file
                median_fname = os.path.abspath("median.nii.gz")
                image_3d = copyInfo_3DImage(sitk.GetImageFromArray(
                    np.median(data_slice, axis=0), isVector=False), in_nii)
                sitk.WriteImage(image_3d, median_fname)

            # First iteration to generate reference image.
            res = antsMotionCorr(in_file=slice_fname,
                                 ref_file=median_fname, prebuilt_option=self.inputs.HMC_option, transform_type='Rigid', second=False, rabies_data_type=self.inputs.rabies_data_type).run()

            median = np.median(sitk.GetArrayFromImage(sitk.ReadImage(
                res.outputs.mc_corrected_bold, self.inputs.rabies_data_type)), axis=0)
            tmp_median_fname = os.path.abspath("tmp_median.nii.gz")
            image_3d = copyInfo_3DImage(
                sitk.GetImageFromArray(median, isVector=False), in_nii)
            sitk.WriteImage(image_3d, tmp_median_fname)

            # Second iteration to generate reference image.
            res = antsMotionCorr(in_file=slice_fname,
                                 ref_file=tmp_median_fname, prebuilt_option=self.inputs.HMC_option, transform_type='Rigid', second=True,  rabies_data_type=self.inputs.rabies_data_type).run()

            # evaluate a trimmed mean instead of a median, trimming the 5% extreme values
            from scipy import stats
            median_image_data = stats.trim_mean(sitk.GetArrayFromImage(sitk.ReadImage(
                res.outputs.mc_corrected_bold, self.inputs.rabies_data_type)), 0.05, axis=0)

        # median_image_data is a 3D array of the median image, so creates a new nii image
        # saves it
        image_3d = copyInfo_3DImage(sitk.GetImageFromArray(
            median_image_data, isVector=False), in_nii)
        sitk.WriteImage(image_3d, out_ref_fname)

        # denoise the resulting reference image through non-local mean denoising
        # Denoising reference image.
        command = f'DenoiseImage -d 3 -i {out_ref_fname} -o {out_ref_fname}'
        rc = run_command(command)

        setattr(self, 'ref_image', out_ref_fname)
        setattr(self, 'bold_file', out_bold_file)

        return runtime

    def _list_outputs(self):
        return {'ref_image': getattr(self, 'ref_image'),
                'bold_file': getattr(self, 'bold_file')}


def _get_vols_to_discard(img):
    '''
    Takes a nifti file, extracts the mean signal of the first 50 volumes and computes which are outliers.
    is_outlier function: computes Modified Z-Scores (https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm) to determine which volumes are outliers.
    '''
    from nipype.algorithms.confounds import is_outlier
    data_slice = sitk.GetArrayFromImage(img)[:50, :, :, :]
    global_signal = data_slice.mean(axis=-1).mean(axis=-1).mean(axis=-1)
    return is_outlier(global_signal)
