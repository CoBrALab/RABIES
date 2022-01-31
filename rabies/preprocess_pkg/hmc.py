import os
import pathlib
import SimpleITK as sitk
import numpy as np
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from rabies.utils import run_command
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

def init_bold_hmc_wf(opts, name='bold_hmc_wf'):
    """
    This workflow estimates motion during fMRI acquisition. To do so, each EPI frame is registered to a volumetric 
    target reference image with a rigid registration using ANTs' antsMotionCorr algorithm (Avants et al., 2009). 
    The resulting 3 translation and 3 rotation realignment parameters are saved for all frames into an output CSV file.

    References:
        Avants, B. B., Tustison, N., & Song, G. (2009). Advanced normalization tools (ANTS). The Insight Journal, 2, 1–35.

    Command line interface parameters:
        --HMC_option {intraSubjectBOLD,0,1,2,3}
                                Select an option for head motion realignment among the pre-built options from
                                https://github.com/ANTsX/ANTsR/blob/master/R/ants_motion_estimation.R.
                                (default: intraSubjectBOLD)
                                
        --apply_slice_mc      Whether to apply a slice-specific motion correction after initial volumetric HMC. This can 
                                correct for interslice misalignment resulting from within-TR motion. With this option, 
                                motion corrections and the subsequent resampling from registration are applied sequentially
                                since the 2D slice registrations cannot be concatenate with 3D transforms. 
                                (default: False)

    Workflow:
        parameters
            opts: command line interface parameters

        inputs
            bold_file: Nifti file with EPI timeseries to realign
            ref_image: the 3D image target for realignment

        outputs
            motcorr_params: CSV file which contains all translation and rotation parameters
            slice_corrected_bold: if using the experimental method --apply_slice_mc, these are the EPI frames 
                after both rigid and then slice-specific realignment
    """

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['bold_file', 'ref_image']),
                        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['motcorr_params', 'slice_corrected_bold']),
        name='outputnode')

    # Head motion correction (hmc)
    motion_estimation = pe.Node(antsMotionCorr(prebuilt_option=opts.HMC_option,transform_type='Rigid', second=False, rabies_data_type=opts.data_type),
                         name='ants_MC', mem_gb=1.1*opts.scale_min_memory)
    motion_estimation.plugin_args = {
        'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}

    workflow.connect([
        (inputnode, motion_estimation, [('ref_image', 'ref_file'),
                                        ('bold_file', 'in_file')]),
        (motion_estimation, outputnode, [
         ('csv_params', 'motcorr_params')]),
    ])

    if opts.apply_slice_mc:
        slice_mc_n_procs = int(opts.local_threads/4)+1
        slice_mc_node = pe.Node(SliceMotionCorrection(n_procs=slice_mc_n_procs),
                                name='slice_mc', mem_gb=1*slice_mc_n_procs, n_procs=slice_mc_n_procs)
        slice_mc_node.plugin_args = {
            'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}

        # conducting a volumetric realignment before slice-specific mc to correct for larger head translations and rotations
        workflow.connect([
            (inputnode, slice_mc_node, [('ref_image', 'ref_file'),
                                        ('bold_file', 'name_source')]),
            (motion_estimation, slice_mc_node, [
             ('mc_corrected_bold', 'in_file')]),
            (slice_mc_node, outputnode, [
             ('mc_corrected_bold', 'slice_corrected_bold')]),
        ])

    return workflow


class antsMotionCorrInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='input BOLD time series')
    ref_file = File(exists=True, mandatory=True,
                    desc='ref file to realignment time series')
    prebuilt_option = traits.Str(desc="Select one of the pre-built options from https://github.com/ANTsX/ANTsR/blob/60eefd96fedd16bceb4703ecd2cd5730e6843807/R/ants_motion_estimation.R")
    transform_type = traits.Str(desc="Specify between Rigid and Affine transform.")
    second = traits.Bool(desc="specify if it is the second iteration")
    rabies_data_type = traits.Int(mandatory=True,
                                  desc="Integer specifying SimpleITK data type.")


class antsMotionCorrOutputSpec(TraitedSpec):
    mc_corrected_bold = File(exists=True, desc="motion corrected time series")
    avg_image = File(
        exists=True, desc="average image of the motion corrected time series")
    csv_params = File(
        exists=True, desc="csv files with the 6-parameters rigid body transformations")


class antsMotionCorr(BaseInterface):
    """
    This interface performs motion realignment using antsMotionCorr function. It takes a reference volume to which
    EPI volumes from the input 4D file are realigned based on a Rigid registration.
    """

    input_spec = antsMotionCorrInputSpec
    output_spec = antsMotionCorrOutputSpec

    def _run_interface(self, runtime):
 
        # change the name of the first iteration directory to prevent overlap of files with second iteration
        if self.inputs.second:
            command = 'mv ants_mc_tmp first_ants_mc_tmp'
            rc = run_command(command)

        # make a tmp directory to store the files
        os.makedirs('ants_mc_tmp', exist_ok=True)

        # adaptation from https://github.com/ANTsX/ANTsR/blob/60eefd96fedd16bceb4703ecd2cd5730e6843807/R/ants_motion_estimation.R
        moving = self.inputs.in_file
        fixed = self.inputs.ref_file
        txtype = self.inputs.transform_type
        moreaccurate = self.inputs.prebuilt_option
        verbose = 0

        if txtype not in ['Rigid', 'Affine']:
            raise ValueError("Wrong transform type provided.")
        if not moreaccurate=="intraSubjectBOLD":
            moreaccurate=int(moreaccurate)
            if moreaccurate not in [0, 1, 2, 3]:
                raise ValueError("Wrong pre-built option provided.")

        img = sitk.ReadImage(self.inputs.in_file, self.inputs.rabies_data_type)
        n = img.GetSize()[3]
        if (n > 10):
            n = 10
        mibins = 20

        if (moreaccurate == 3):
            # check the size of the lowest dimension, and make sure that the first shrinking factor allow for at least 4 slices
            shrinking_factor = 4
            low_dim = np.asarray(img.GetSize()[:3]).min()
            if shrinking_factor > int(low_dim/4):
                shrinking_factor = int(low_dim/4)
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} ] -t {txtype}[0.1,3,0] -i 100x50x30 -s 2x1x0 -f {str(shrinking_factor)}x2x1 -u 1 -e 1 -l 1 -n {str(n)} -v {str(verbose)}"

        elif (moreaccurate == 2):
            # check the size of the lowest dimension, and make sure that the first shrinking factor allow for at least 4 slices
            shrinking_factor = 4
            img = sitk.ReadImage(self.inputs.in_file, self.inputs.rabies_data_type)
            low_dim = np.asarray(img.GetSize()[:3]).min()
            if shrinking_factor > int(low_dim/4):
                shrinking_factor = int(low_dim/4)
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} , regular, 0.25 ] -t {txtype}[ 0.1 ] -i 100x50x30 -s 2x1x0 -f {str(shrinking_factor)}x2x1 -u 1 -e 1 -l 1 -n {str(n)} -v {str(verbose)}"

        elif (moreaccurate == "intraSubjectBOLD"):
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} , regular, 0.2 ] -t {txtype}[ 0.25 ] -i 50x20 -s 1x0 -f 2x1 -u 1 -e 1 -l 1 -n {str(n)} -v {str(verbose)}"

        elif (moreaccurate == 1):
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} , regular, 0.25 ] -t {txtype}[ 0.1 ] -i 100 -s 0 -f 1 -u 1 -e 1 -l 1 -n {str(n)} -v {str(verbose)}"

        elif (moreaccurate == 0):
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} , regular, 0.02 ] -t {txtype}[ 0.1 ] -i 3 -s 0 -f 1 -u 1 -e 1 -l 1 -n {str(n)} -v {str(verbose)}"
        else:
            raise ValueError("Wrong moreaccurate provided.")
        rc = run_command(command)

        setattr(self, 'csv_params', os.path.abspath('ants_mc_tmp/motcorrMOCOparams.csv'))
        setattr(self, 'mc_corrected_bold', os.path.abspath('ants_mc_tmp/motcorr.nii.gz'))
        setattr(self, 'avg_image', os.path.abspath('ants_mc_tmp/motcorr_avg.nii.gz'))

        return runtime

    def _list_outputs(self):
        return {'mc_corrected_bold': getattr(self, 'mc_corrected_bold'),
                'csv_params': getattr(self, 'csv_params'),
                'avg_image': getattr(self, 'avg_image')}


def register_slice(fixed_image, moving_image):
    # function for 2D registration
    dimension = 2
    initial_transform = sitk.Transform(dimension, sitk.sitkIdentity)

    registration_method = sitk.ImageRegistrationMethod()

    # Similarity metric settings.
    registration_method.SetMetricAsMattesMutualInformation(
        numberOfHistogramBins=20)
    registration_method.SetMetricSamplingStrategy(registration_method.NONE)
    # registration_method.SetMetricSamplingPercentage(0.01)

    registration_method.SetInterpolator(sitk.sitkBSpline)

    # Optimizer settings.
    registration_method.SetOptimizerAsGradientDescent(
        learningRate=0.05, numberOfIterations=100, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    # registration_method.SetOptimizerScalesFromPhysicalShift()

    # Setup for the multi-resolution framework.
    registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[4, 2, 1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2, 1, 0])
    # registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    # Don't optimize in-place, we would possibly like to run this cell multiple times.
    registration_method.SetInitialTransform(initial_transform, inPlace=False)

    final_transform = registration_method.Execute(sitk.Cast(fixed_image, sitk.sitkFloat32),
                                                  sitk.Cast(moving_image, sitk.sitkFloat32))
    return final_transform


def slice_specific_registration(i, ref_file, timeseries_file):
    from nipype import logging
    log = logging.getLogger('nipype.workflow')

    log.debug('Slice-specific correction on volume '+str(i+1))
    ref_image = sitk.ReadImage(ref_file, sitk.sitkFloat32)
    timeseries_image = sitk.ReadImage(timeseries_file, sitk.sitkFloat32)
    volume_array = sitk.GetArrayFromImage(timeseries_image)[i, :, :, :]

    for j in range(volume_array.shape[1]):
        slice_array = volume_array[:, j, :]
        if slice_array.sum()==0:
            continue
        moving_image = sitk.GetImageFromArray(slice_array)
        fixed_image = ref_image[:, j, :]
        moving_image.CopyInformation(fixed_image)

        final_transform = register_slice(fixed_image, moving_image)
        moving_resampled = sitk.Resample(moving_image, fixed_image, final_transform,
                                         sitk.sitkBSplineResamplerOrder4, 0.0, moving_image.GetPixelID())

        resampled_slice = sitk.GetArrayFromImage(moving_resampled)
        volume_array[:, j, :] = resampled_slice
    return [i, volume_array]


class SliceMotionCorrectionInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='input BOLD time series')
    ref_file = File(exists=True, mandatory=True,
                    desc='ref file to realignment time series')
    name_source = File(exists=True, mandatory=True,
                       desc='Reference BOLD file for naming the output.')
    n_procs = traits.Int(exists=True, mandatory=True,
                         desc="Number of processors available to run in parallel.")


class SliceMotionCorrectionOutputSpec(TraitedSpec):
    mc_corrected_bold = File(exists=True, desc="motion corrected time series")


class SliceMotionCorrection(BaseInterface):
    """
    This interface performs slice-specific motion realignment of coronal slices to correct for interslice
    misalignment issues that arise from within-TR motion. It relies on 2D Rigid registration to the
    reference 3D EPI volume provided.
    """

    input_spec = SliceMotionCorrectionInputSpec
    output_spec = SliceMotionCorrectionOutputSpec

    def _run_interface(self, runtime):

        import multiprocessing as mp

        timeseries_image = sitk.ReadImage(
            self.inputs.in_file, sitk.sitkFloat32)

        pool = mp.Pool(processes=self.inputs.n_procs)
        results = [pool.apply_async(slice_specific_registration, args=(
            i, self.inputs.ref_file, self.inputs.in_file)) for i in range(timeseries_image.GetSize()[3])]
        results = [p.get() for p in results]
        # enforce proper order of the slices
        results.sort()
        results = [r[1] for r in results]

        timeseries_array = sitk.GetArrayFromImage(timeseries_image)
        for i in range(timeseries_image.GetSize()[3]):
            timeseries_array[i, :, :, :] = results[i]

        # clip potential negative values
        timeseries_array[(timeseries_array < 0).astype(bool)] = 0
        resampled_timeseries = sitk.GetImageFromArray(
            timeseries_array, isVector=False)
        resampled_timeseries.CopyInformation(timeseries_image)

        split = pathlib.Path(self.inputs.name_source).name.rsplit(".nii")
        out_name = os.path.abspath(split[0]+'_slice_mc.nii.gz')
        sitk.WriteImage(resampled_timeseries, out_name)

        setattr(self, 'mc_corrected_bold', out_name)

        return runtime

    def _list_outputs(self):
        return {'mc_corrected_bold': getattr(self, 'mc_corrected_bold')}
