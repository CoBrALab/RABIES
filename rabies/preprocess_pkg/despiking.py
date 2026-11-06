import os
import pathlib
import numpy as np
import SimpleITK as sitk
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


class DespikingInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="4D EPI timeseries")
    smoothing_sigma = traits.Any(None, usedefault=True,
        desc="Standard deviation (in mm) of the in-plane Gaussian smoothing applied before spike detection. "
             "If None, no smoothing is applied.")
    threshold = traits.Float(4.0, usedefault=True,
        desc="Voxel-frames with a spike score at or above this threshold (in units of each voxel's standard "
             "deviation, as estimated by AFNI's 3dDespike) are flagged as spikes.")
    video_dir = traits.Any(None, usedefault=True,
        desc="If provided, a video of the timeseries before and after despiking is written in this folder.")
    name_source = traits.Any(None, usedefault=True,
        desc="Reference file for naming the video. If None, in_file is used.")
    fps = traits.Float(10, usedefault=True, desc="Frames per second for the video.")


class DespikingOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="Despiked 4D EPI timeseries")
    spike_flags = File(exists=True, desc="4D mask of the voxel-frames flagged as spikes (1=spike)")


class Despiking(BaseInterface):
    """
    Detects temporal spikes voxelwise, and replaces them by linear interpolation in time.

    1. Optionally, each frame is smoothed in-plane with a Gaussian of standard deviation smoothing_sigma (in mm).
       There is no smoothing across slices (the axis with the largest voxel dimension, or the superior-inferior
       axis for isotropic voxels) nor across time.
    2. AFNI's 3dDespike (https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html) is run on the
       (smoothed) timeseries, only to obtain the spike score of each voxel-frame: its deviation from a smooth curve
       fitted to the voxel timecourse, in units of the voxel's standard deviation.
    3. Voxel-frames with a score at or above the threshold are flagged.
    4. In the original (unsmoothed) timeseries, the flagged voxel-frames are replaced by linear interpolation in
       time between the nearest unflagged frames of the same voxel (at the ends of the scan, the nearest
       unflagged value is copied). All other voxel-frames are unchanged.

    Default behavior (smoothing_sigma=None, threshold=4.0): a voxelwise spike correction close to AFNI's 3dDespike
    defaults settings using the -localedit option. AFNI flags voxel-frames beyond 4 standard deviations (its default
    cut c2, at which it caps values) and, with the -localedit option, replaces them with neighboring values, as
    done here.

    Alternative recipe (smoothing_sigma=1, threshold=2.0): smoothed spike detection, developed for functional
    ultrasound (fUSI) data, where flashes are spatially coherent signal increases across a slice or a blob of
    tissue. Voxelwise despiking faces a tradeoff with head motion correction: a brief motion displaces tissue
    edges, which produces voxelwise intensity changes that look exactly like spikes. A stringent voxelwise
    threshold thus removes flashes but also erases the image signature of brief real motion, which the motion
    correction can no longer detect. Detecting spikes on spatially smoothed data resolves this tradeoff:
    smoothing averages out the local intensity changes caused by small tissue displacements, while spatially
    coherent artefacts survive, so a lower threshold can be used to remove artefacts aggressively while
    preserving the signature of real motion. Since flagged voxel-frames are repaired from the original data,
    voxel by voxel, the output is not blurred.

    Optionally (video_dir), a video of the timeseries before and after despiking is generated, with a single
    intensity scale across frames and rows (1st to 99th percentiles) so that the removed spikes remain visible.
    """

    input_spec = DespikingInputSpec
    output_spec = DespikingOutputSpec

    def _run_interface(self, runtime):
        from nipype import logging
        log = logging.getLogger('nipype.workflow')
        from rabies.utils import run_command

        n_procs = int(os.environ['RABIES_ITK_NUM_THREADS']) if "RABIES_ITK_NUM_THREADS" in os.environ else os.cpu_count()
        sigma = self.inputs.smoothing_sigma
        if sigma is not None and not sigma > 0:
            raise ValueError(f"smoothing_sigma must be None or above 0, {sigma} was provided.")

        split = pathlib.Path(self.inputs.in_file).name.rsplit(".nii")[0]
        out_file = os.path.abspath(f"{split}_despike.nii.gz")
        flags_file = os.path.abspath(f"{split}_spike_flags.nii.gz")

        img = sitk.ReadImage(self.inputs.in_file, sitk.sitkFloat32)
        if not img.GetDimension() == 4:
            raise ValueError(f"Input image {self.inputs.in_file} is not 4-dimensional.")
        n_t = img.GetSize()[3]
        if n_t < 15:
            raise ValueError(f"Despiking requires at least 15 frames, but {self.inputs.in_file} has {n_t}.")
        raw = sitk.GetArrayFromImage(img)  # (t, z, y, x)

        # 1. optional in-plane smoothing, then 3dDespike on the (smoothed) timeseries to obtain the spike score
        if sigma is None:
            detect_file = self.inputs.in_file
        else:
            spacing_zyx = img.GetSpacing()[:3][::-1]
            smoothed = _smooth_inplane(raw, spacing_zyx, sigma, n_procs)
            smoothed_img = sitk.GetImageFromArray(smoothed, isVector=False)
            smoothed_img.CopyInformation(img)
            del smoothed
            detect_file = os.path.abspath(f"{split}_smoothed.nii.gz")
            sitk.WriteImage(smoothed_img, detect_file)

        score_file = os.path.abspath(f"{split}_spike_score.nii.gz")
        afni_out = os.path.abspath(f"{split}_afni_despike.nii.gz")
        # -nomask: score all voxels (by default 3dDespike zeroes voxels outside its automask)
        # -ssave: saves the spike score of each voxel-frame (unsigned, at 0.1 precision)
        command = f'3dDespike -nomask -ssave {score_file} -prefix {afni_out} {detect_file}'
        rc, c_out = run_command(command)
        score = sitk.GetArrayFromImage(sitk.ReadImage(score_file, sitk.sitkFloat32))
        # only the spike score is used; the files generated for detection are removed
        for f in [afni_out, score_file] + ([detect_file] if sigma is not None else []):
            os.remove(f)

        # 2. flag spikes, and replace them by linear interpolation in time from the original timeseries
        flags = score >= self.inputs.threshold - 1e-6  # the score is saved at 0.1 precision
        log.info(f"Despiking: {flags.mean()*100:.3f}% of voxel-frames flagged as spikes.")
        corrected = _interpolate_flagged(raw.reshape(n_t, -1), flags.reshape(n_t, -1)).reshape(raw.shape)

        out_img = sitk.GetImageFromArray(corrected, isVector=False)
        out_img.CopyInformation(img)
        sitk.WriteImage(out_img, out_file)
        flags_img = sitk.GetImageFromArray(flags.astype(np.uint8), isVector=False)
        flags_img.CopyInformation(img)
        sitk.WriteImage(flags_img, flags_file)

        # 3. optional video before/after despiking
        if self.inputs.video_dir is not None:
            from simpleitk_timeseries_motion_correction.create_animation import main as create_animation
            name_source = self.inputs.in_file if self.inputs.name_source is None else self.inputs.name_source
            video_name = pathlib.Path(name_source).name.rsplit(".nii")[0]
            os.makedirs(self.inputs.video_dir, exist_ok=True)
            # a single intensity range across frames and rows, so that spikes remain visible
            create_animation(input_img=img,
                             output_file=os.path.join(self.inputs.video_dir, f'{video_name}_despiking.webp'),
                             additional_input_imgs=[out_img],
                             labels=['Before despiking', 'After despiking'],
                             scale=2.0,
                             fps=self.inputs.fps,
                             percentile_range=(1, 99),
                             framewise_scaling=False)

        setattr(self, 'out_file', out_file)
        setattr(self, 'spike_flags', flags_file)
        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file'),
                'spike_flags': getattr(self, 'spike_flags')}


def _smooth_inplane(arr, spacing_zyx, sigma_mm, n_procs):
    '''
    Gaussian smoothing of each frame of a (t, z, y, x) array, with standard deviation sigma_mm, in-plane only:
    there is no smoothing along time, nor along the slice axis (the axis with the largest voxel dimension; for
    isotropic voxels, z, the superior-inferior axis).
    '''
    import concurrent.futures
    from scipy.ndimage import gaussian_filter
    slice_axis = int(np.argmax(spacing_zyx))  # argmax returns the first maximum, so z in case of ties
    sigma_vox = [sigma_mm / s for s in spacing_zyx]
    sigma_vox[slice_axis] = 0
    out = np.empty_like(arr)
    chunks = np.array_split(np.arange(arr.shape[0]), max(1, min(n_procs, arr.shape[0])))

    def _run(idx):
        if len(idx) > 0:
            out[idx] = gaussian_filter(arr[idx], sigma=(0,) + tuple(sigma_vox), mode='nearest')
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_procs) as executor:
        list(executor.map(_run, chunks))
    return out


def _interpolate_flagged(data, flags):
    '''
    Replace flagged entries of (t, n) data by linear interpolation in time between each column's nearest unflagged
    entries (np.interp copies the nearest one at the ends); columns without any unflagged entry are left unchanged.
    Returns a new array.
    '''
    out = data.copy()
    frames = np.arange(data.shape[0])
    for i in np.where(flags.any(axis=0) & ~flags.all(axis=0))[0]:
        bad = flags[:, i]
        out[bad, i] = np.interp(frames[bad], frames[~bad], data[~bad, i])
    return out
