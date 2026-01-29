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
from rabies.utils import copyInfo_3DImage, run_command
from rabies.simpleitk_timeseries_motion_correction.motion import framewise_register_pair
from rabies.simpleitk_timeseries_motion_correction.apply_transforms import framewise_resample_volume

def init_bold_reference_wf(opts, name='gen_bold_ref'):
    # gen_bold_ref_head_start
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
    """
    # gen_bold_ref_head_end

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['ref_image']),
        name='outputnode')

    n_procs=int(os.environ['RABIES_ITK_NUM_THREADS'])
    gen_ref = pe.Node(EstimateReferenceImage(HMC_level=opts.HMC_level, detect_dummy=opts.detect_dummy, rabies_data_type=opts.data_type),
                      name='gen_ref', mem_gb=2*opts.scale_min_memory, n_procs=n_procs)
    gen_ref.plugin_args = {
        'qsub_args': f'-pe smp {str(2*opts.min_proc)}', 'overwrite': True}

    workflow.connect([
        (inputnode, gen_ref, [
            ('bold_file', 'in_file'),
            ]),
        (gen_ref, outputnode, [('ref_image', 'ref_image')]),
    ])

    return workflow


class EstimateReferenceImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="4D EPI file")
    HMC_level = traits.Int(desc="Level for motion correction.")
    detect_dummy = traits.Bool(
        desc="specify if should detect and remove dummy scans, and use these volumes as reference image.")
    rabies_data_type = traits.Int(mandatory=True,
                                  desc="Integer specifying SimpleITK data type.")

class EstimateReferenceImageOutputSpec(TraitedSpec):
    ref_image = File(exists=True, desc="3D reference image")

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

        n_procs = int(os.environ['RABIES_ITK_NUM_THREADS']) if "RABIES_ITK_NUM_THREADS" in os.environ else os.cpu_count() # default to number of CPUs

        filename_split = pathlib.Path(self.inputs.in_file).name.rsplit(".nii")
        out_ref_fname = os.path.abspath(
            f'{filename_split[0]}_bold_ref.nii.gz')

        in_nii = sitk.ReadImage(self.inputs.in_file,
                                self.inputs.rabies_data_type)
        if not in_nii.GetDimension()==4:
            raise ValueError(f"Input image {self.inputs.in_file} is not 4-dimensional.")

        dummy_ref = False
        if self.inputs.detect_dummy:
            n_volumes_to_discard = _get_vols_to_discard(in_nii)
            if (not n_volumes_to_discard == 0):
                log.info("Detected "+str(n_volumes_to_discard)
                    + " dummy scans. Taking the mean of these volumes as reference EPI.")
                trimean_array = np.quantile(sitk.GetArrayFromImage(in_nii[:,:,:,:n_volumes_to_discard]), (0.2, 0.5, 0.8), axis=0).mean(axis=0)
                ref_3d = copyInfo_3DImage(sitk.GetImageFromArray(
                    trimean_array, isVector=False), in_nii)
                dummy_ref = True
            else:
                log.info(
                    "Detected no dummy scans. Generating the ref EPI based on multiple volumes.")
                
        if not dummy_ref:
            num_timepoints = in_nii.GetSize()[-1]
            if num_timepoints > 50:
                # select a set of 50 frames spread uniformally across time to avoid temporal biases
                subset_idx = np.linspace(0,num_timepoints-1,50).astype(int)
                subset_img_4d = sitk.JoinSeries([in_nii[:,:,:,int(i)] for i in subset_idx])
            else:
                subset_img_4d = in_nii

            trimean_array = np.quantile(sitk.GetArrayFromImage(subset_img_4d), (0.2, 0.5, 0.8), axis=0).mean(axis=0)
            ref_3d = copyInfo_3DImage(
                sitk.GetImageFromArray(trimean_array, isVector=False), in_nii)

            for round in range(2): # conducted 2 rounds of motion correction and re-calculation of the 3D reference
                transforms = framewise_register_pair(subset_img_4d, ref_3d, level=self.inputs.HMC_level, interpolation=sitk.sitkBSpline5, max_workers=n_procs)
                resampled_img = framewise_resample_volume(subset_img_4d, ref_3d, transforms, interpolation=sitk.sitkBSpline5, max_workers=n_procs)

                trimean_array = np.quantile(sitk.GetArrayFromImage(resampled_img), (0.2, 0.5, 0.8), axis=0).mean(axis=0)
                ref_3d = copyInfo_3DImage(
                    sitk.GetImageFromArray(trimean_array, isVector=False), in_nii)

        # saves the finalized average
        sitk.WriteImage(ref_3d, out_ref_fname)

        # denoise the resulting reference image through non-local mean denoising
        # Denoising reference image.
        command = f'DenoiseImage -d 3 -i {out_ref_fname} -o {out_ref_fname}'
        rc,c_out = run_command(command)

        setattr(self, 'ref_image', out_ref_fname)

        return runtime

    def _list_outputs(self):
        return {'ref_image': getattr(self, 'ref_image')}


def _get_vols_to_discard(img):
    '''
    Takes a nifti file, extracts the mean signal of the first 50 volumes and computes which are outliers.
    is_outlier function: computes Modified Z-Scores (https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm) to determine which volumes are outliers.
    '''
    from nipype.algorithms.confounds import is_outlier
    data_slice = sitk.GetArrayFromImage(img)[:50, :, :, :]
    global_signal = data_slice.mean(axis=-1).mean(axis=-1).mean(axis=-1)
    return is_outlier(global_signal)
