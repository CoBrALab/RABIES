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
    # hmc_wf_head_start
    """
    This workflow estimates motion during fMRI acquisition. To do so, each EPI frame is registered to a volumetric 
    target reference image with a rigid registration using ANTs' antsMotionCorr algorithm (Avants et al., 2009). 
    This results in the measurement of 3 Euler angles in radians and 3 translations in mm (from ITK's 
    Euler3DTransform https://itk.org/Doxygen/html/classitk_1_1Euler3DTransform.html) at each time frame, which are 
    then stored into an output CSV file.

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
    # hmc_wf_head_end

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
            rc,c_out = run_command(command)

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
        if not (moreaccurate=="intraSubjectBOLD" or moreaccurate=="optim"):
            moreaccurate=int(moreaccurate)
            if moreaccurate not in [0, 1, 2, 3]:
                raise ValueError("Wrong pre-built option provided.")

        img = sitk.ReadImage(moving, self.inputs.rabies_data_type)
        ref_img = sitk.ReadImage(fixed, self.inputs.rabies_data_type)

        n = img.GetSize()[3]
        if (n > 10):
            n = 10
        mibins = 20

        # check the size of the lowest dimension, and make sure that the first shrinking factor allow for at least 4 slices
        shrinking_factor = 4
        low_dim = np.asarray(ref_img.GetSize()[:3]).min()
        if shrinking_factor > int(low_dim/4):
            shrinking_factor = int(low_dim/4)

        if (moreaccurate == 3):
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} ] -t {txtype}[0.1,3,0] -i 100x50x30 -s 2x1x0mm -f {str(shrinking_factor)}x2x1 -u 1 -e 1 -n {str(n)} -v {str(verbose)}"

        elif (moreaccurate == 2):
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} , regular, 0.25 ] -t {txtype}[ 0.1 ] -i 100x50x30 -s 2x1x0mm -f {str(shrinking_factor)}x2x1 -u 1 -e 1 -n {str(n)} -v {str(verbose)}"

        elif (moreaccurate == "intraSubjectBOLD"):
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} , regular, 0.2 ] -t {txtype}[ 0.25 ] -i 50x20 -s 1x0mm -f 2x1 -u 1 -e 1 -n {str(n)} -v {str(verbose)}"

        elif (moreaccurate == 1):
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} , regular, 0.25 ] -t {txtype}[ 0.1 ] -i 100 -s 0 -f 1 -u 1 -e 1 -n {str(n)} -v {str(verbose)}"

        elif (moreaccurate == 0):
            command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                {moving} , 1 , {str(mibins)} , regular, 0.02 ] -t {txtype}[ 0.1 ] -i 3 -s 0 -f 1 -u 1 -e 1 -n {str(n)} -v {str(verbose)}"
        
        elif (moreaccurate == "optim"):

            # generate sensible smoothing coefficients based on image dimensions
            low_dim = np.asarray(ref_img.GetSpacing()[:3]).min()
            largest_dim = (np.array(ref_img.GetSize()[:3])*np.array(ref_img.GetSpacing()[:3])).max()

            command=f'ants_generate_iterations.py --min {low_dim} --max {largest_dim}'
            import subprocess
            from nipype import logging
            log = logging.getLogger('nipype.workflow')
            log.debug('Running: '+command)
            try:
                process = subprocess.run(
                    command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    check=True,
                    shell=True,
                    )
            except Exception as e:
                log.warning(e.output.decode("utf-8"))
                raise
            out = process.stdout.decode("utf-8")
            log.debug(out)

            s = out.split('--smoothing-sigmas ')[-1].split('mm')[0].split('x')[-2:] # taking the last 2 smoothing sigmas
            f = out.split('--shrink-factors ')[-1].split(' ')[0].split('x')[-2:] # taking shrink factors

            if len(s)==1:
                command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                    {moving} , 1 , {str(mibins)} , Regular, 0.25, 1 ] -t {txtype}[ 0.1 ] -i 20 -s 0.0mm -f 1 -u 1 -e 1 -n {str(n)} -v {str(verbose)}"
            elif len(s)>1:
                command = f"antsMotionCorr -d 3 -o [ants_mc_tmp/motcorr,ants_mc_tmp/motcorr.nii.gz,ants_mc_tmp/motcorr_avg.nii.gz] -m MI[ {fixed} , \
                    {moving} , 1 , {str(mibins)} , Regular, 0.25, 1 ] -t {txtype}[ 0.1 ] -i 50x20 -s {s[0]}x{s[1]}mm -f {f[0]}x{f[1]} -u 1 -e 1 -n {str(n)} -v {str(verbose)}"
            else:
                raise ValueError("No smoothing coefficient was found.")
        else:
            raise ValueError("Wrong moreaccurate provided.")
        rc,c_out = run_command(command)

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

    registration_method.SetInterpolator(sitk.sitkLinear)

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
                                         sitk.sitkLinear, 0.0, moving_image.GetPixelID())

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


class EstimateMotionParamsInputSpec(BaseInterfaceInputSpec):
    motcorr_params = File(exists=True, mandatory=True,
                       desc="CSV file with the 6 rigid body parameters")
    raw_bold = File(exists=True, mandatory=True,
                      desc="Raw EPI before resampling.")
    raw_brain_mask = File(exists=True, mandatory=True,
                      desc="Brain mask of the raw EPI.")


class EstimateMotionParamsOutputSpec(TraitedSpec):
    motion_params_csv = traits.File(desc="CSV file of motion parameters")
    FD_csv = traits.File(desc="CSV file with global framewise displacement.")
    FD_voxelwise = traits.File(
        desc=".nii file with voxelwise framewise displacement.")
    pos_voxelwise = traits.File(desc=".nii file with voxelwise Positioning.")


class EstimateMotionParams(BaseInterface):
    # motion_param_head_start
    """
    This interface generates estimations of absolute displacement and framewise displacement, together with
    the expansion of the 6 motion parameters to include derivatives and squared parameters (Friston 24).
    Absolute and framewise displacement are computed within antsMotionCorrStats as follows:
        1. For each timepoint, the 3 Euler rotations and translations are converted to an affine matrix
        2. For each voxel within a brain mask representing the referential space post-motion realignment,
           the inverse transform is applied to generate a point pre-motion realignment.
        3. Absolute displacement is computed as the distance between the referential point post-correction 
           and the point pre-correction generated from the affine. For framewise displacement, the 
           distance is measured between the pre-correction points generated from the current and the 
           next timeframes. Distance is measured in mm with the Euclidean distance.
        4. From the distance measurements, voxelwise 4D timeseries are generated, and for framewise
           displacement, the mean and max displacement at each timeframe is stored in a CSV file.
    """
    # motion_param_head_end

    input_spec = EstimateMotionParamsInputSpec
    output_spec = EstimateMotionParamsOutputSpec

    def _run_interface(self, runtime):
        import numpy as np
        import os
        from rabies.utils import run_command
        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(self.inputs.raw_bold).name.rsplit(".nii")

        # generate a .nii file representing the positioning or framewise displacement for each voxel within the brain_mask
        # first the voxelwise positioning map
        command = f'antsMotionCorrStats -m {self.inputs.motcorr_params} -o {filename_split[0]}_pos_file.csv -x {self.inputs.raw_brain_mask} \
                    -d {self.inputs.raw_bold}'
        rc,c_out = run_command(command)
        pos_voxelwise = os.path.abspath(
            f"{filename_split[0]}_pos_file.nii.gz")

        # then the voxelwise framewise displacement map
        command = f'antsMotionCorrStats -m {self.inputs.motcorr_params} -o {filename_split[0]}_FD_file.csv -x {self.inputs.raw_brain_mask} \
                    -d {self.inputs.raw_bold} -f 1'
        rc,c_out = run_command(command)

        FD_csv = os.path.abspath(f"{filename_split[0]}_FD_file.csv")
        FD_voxelwise = os.path.abspath(f"{filename_split[0]}_FD_file.nii.gz")

        motion_24,motion_24_header = motion_24_params(self.inputs.motcorr_params)
        # write into a .csv
        import pandas as pd
        df = pd.DataFrame(motion_24)
        df.columns = motion_24_header
        motion_params_csv = os.path.abspath(f"{filename_split[0]}_motion_params.csv")
        df.to_csv(motion_params_csv)

        setattr(self, 'FD_csv', FD_csv)
        setattr(self, 'motion_params_csv', motion_params_csv)
        setattr(self, 'FD_voxelwise', FD_voxelwise)
        setattr(self, 'pos_voxelwise', pos_voxelwise)

        return runtime

    def _list_outputs(self):
        return {'motion_params_csv': getattr(self, 'motion_params_csv'),
                'FD_csv': getattr(self, 'FD_csv'),
                'pos_voxelwise': getattr(self, 'pos_voxelwise'),
                'FD_voxelwise': getattr(self, 'FD_voxelwise')}


def motion_24_params(movpar_csv):
    '''
    motioncorr_24params: 6 head motion parameters, their temporal derivative, and the 12 corresponding squared items (Friston et al. 1996, Magn. Reson. Med.)
    '''
    import numpy as np
    rigid_params = extract_rigid_movpar(movpar_csv)
    rotations = rigid_params[:,:3] # rotations are listed first
    translations = rigid_params[:,3:] # translations are last

    movpar = np.zeros([np.size(rigid_params, 0), 24])
    movpar[:, :3] = translations # add rotations first
    movpar[:, 3:6] = rotations # rotations second
    for i in range(6):
        # Compute temporal derivative as difference between two neighboring points
        movpar[0, 6+i] = 0
        movpar[1:, 6+i] = movpar[1:, i]-movpar[:-1, i]
        # add the squared coefficients
        movpar[:, 12+i] = movpar[:, i]**2
        movpar[:, 18+i] = movpar[:, 6+i]**2

    motion_24_header = ['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3', 'mov1_der', 'mov2_der', 'mov3_der', 'rot1_der', 'rot2_der', 'rot3_der',
                    'mov1^2', 'mov2^2', 'mov3^2', 'rot1^2', 'rot2^2', 'rot3^2', 'mov1_der^2', 'mov2_der^2', 'mov3_der^2', 'rot1_der^2', 'rot2_der^2', 'rot3_der^2']

    return movpar,motion_24_header


def extract_rigid_movpar(movpar_csv):
    import numpy as np
    import csv
    temp = []
    with open(movpar_csv) as csvfile:
        motcorr = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in motcorr:
            temp.append(row)
    movpar = np.zeros([(len(temp)-1), 6])
    j = 0
    for row in temp[1:]:
        for i in range(2, len(row)):
            movpar[j, i-2] = float(row[i])
        j = j+1
    return movpar
