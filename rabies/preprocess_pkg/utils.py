import os
import SimpleITK as sitk
import numpy as np
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, InputMultiPath, BaseInterface
)


def prep_bids_iter(layout, bold_only=False):
    '''
    This function takes as input a BIDSLayout, and generates iteration lists
    for managing the workflow's iterables depending on whether --bold_only is
    selected or not.
    '''
    import pathlib

    scan_info = []
    split_name = []
    scan_list = []
    run_iter = {}
    bold_scan_list = []

    subject_list = layout.get_subject()
    if len(subject_list) == 0:
        raise ValueError(
            "No subject information could be retrieved from the BIDS directory. The 'sub-' specification is mandatory.")
    if not bold_only:
        anat_bids = layout.get(subject=subject_list, suffix=[
                               'T2w', 'T1w'], extension=['nii', 'nii.gz'])
        if len(anat_bids) == 0:
            raise ValueError(
                "No anatomical file with the suffix 'T2w' or 'T1w' were found among the BIDS directory.")
    bold_bids = layout.get(subject=subject_list, suffix=[
                           'bold', 'cbv'], extension=['nii', 'nii.gz'])
    if len(bold_bids) == 0:
        raise ValueError(
            "No functional file with the suffix 'bold' or 'cbv' were found among the BIDS directory.")

    bold_dict = {}
    for bold in bold_bids:
        sub = bold.get_entities()['subject']
        try:
            ses = bold.get_entities()['session']
        except:
            ses = None

        try:
            run = bold.get_entities()['run']
        except:
            run = None

        if sub not in list(bold_dict.keys()):
            bold_dict[sub] = {}
        if ses not in list(bold_dict[sub].keys()):
            bold_dict[sub][ses] = {}

        bold_list = layout.get(subject=sub, session=ses, run=run, suffix=[
                               'bold', 'cbv'], extension=['nii', 'nii.gz'], return_type='filename')
        bold_dict[sub][ses][run] = bold_list

    # if not bold_only, then the bold_list and run_iter will be a dictionary with keys being the anat filename
    # otherwise, it will be a list of bold scans themselves
    for sub in list(bold_dict.keys()):
        for ses in list(bold_dict[sub].keys()):
            if not bold_only:
                anat_list = layout.get(subject=sub, session=ses, suffix=[
                                       'T2w', 'T1w'], extension=['nii', 'nii.gz'], return_type='filename')
                if len(anat_list) == 0:
                    raise ValueError(
                        f'Missing an anatomical image for sub {sub} and ses- {ses}')
                if len(anat_list) > 1:
                    raise ValueError(
                        f'Duplicate was found for the anatomical file sub- {sub}, ses- {ses}: {str(anat_list)}')
                file = anat_list[0]
                scan_list.append(file)
                filename_template = pathlib.Path(file).name.rsplit(".nii")[0]
                split_name.append(filename_template)
                scan_info.append({'subject_id': sub, 'session': ses})
                run_iter[filename_template] = []

            for run in list(bold_dict[sub][ses].keys()):
                bold_list = bold_dict[sub][ses][run]
                if len(bold_list) > 1:
                    raise ValueError(
                        f'Duplicate was found for bold files sub- {sub}, ses- {ses} and run {run}: {str(bold_list)}')
                file = bold_list[0]
                bold_scan_list.append(file)
                if bold_only:
                    scan_list.append(file)
                    filename_template = pathlib.Path(
                        file).name.rsplit(".nii")[0]
                    split_name.append(filename_template)
                    scan_info.append(
                        {'subject_id': sub, 'session': ses, 'run': run})
                else:
                    run_iter[filename_template].append(run)

    return split_name, scan_info, run_iter, scan_list, bold_scan_list


class BIDSDataGraberInputSpec(BaseInterfaceInputSpec):
    bids_dir = traits.Str(exists=True, mandatory=True,
                          desc="BIDS data directory")
    suffix = traits.List(exists=True, mandatory=True,
                         desc="Suffix to search for")
    scan_info = traits.Dict(exists=True, mandatory=True,
                            desc="Info required to find the scan")
    run = traits.Any(exists=True, desc="Run number")


class BIDSDataGraberOutputSpec(TraitedSpec):
    out_file = File(
        exists=True, desc="Selected file based on the provided parameters.")


class BIDSDataGraber(BaseInterface):
    """
    This interface will select a single scan from the BIDS directory based on the
    input specifications.
    """

    input_spec = BIDSDataGraberInputSpec
    output_spec = BIDSDataGraberOutputSpec

    def _run_interface(self, runtime):
        subject_id = self.inputs.scan_info['subject_id']
        session = self.inputs.scan_info['session']
        if 'run' in (self.inputs.scan_info.keys()):
            run = self.inputs.scan_info['run']
        else:
            run = self.inputs.run

        from bids.layout import BIDSLayout
        layout = BIDSLayout(self.inputs.bids_dir, validate=False)
        try:
            if run is None: # if there is no run spec to search, don't include it in the search
                file_list = layout.get(subject=subject_id, session=session, extension=[
                                  'nii', 'nii.gz'], suffix=self.inputs.suffix, return_type='filename')
            else:
                file_list = layout.get(subject=subject_id, session=session, run=run, extension=[
                                  'nii', 'nii.gz'], suffix=self.inputs.suffix, return_type='filename')
            if len(file_list) > 1:
                raise ValueError(f'Provided BIDS spec lead to duplicates: \
                    {str(self.inputs.suffix)} sub-{subject_id} ses-{session} run-{str(run)}')
            elif len(file_list)==0:
                raise ValueError(f'No file for found corresponding to the following BIDS spec: \
                    {str(self.inputs.suffix)} sub-{subject_id} ses-{session} run-{str(run)}')
        except:
            raise ValueError(f'Error with BIDS spec: \
                    {str(self.inputs.suffix)} sub-{subject_id} ses-{session} run-{str(run)}')

        setattr(self, 'out_file', file_list[0])

        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file')}


def init_bold_reference_wf(opts, name='gen_bold_ref'):

    from nipype.pipeline import engine as pe
    from nipype.interfaces import utility as niu

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

        import SimpleITK as sitk
        import numpy as np
        import logging
        log = logging.getLogger('root')

        in_nii = sitk.ReadImage(self.inputs.in_file,
                                self.inputs.rabies_data_type)
        data_slice = sitk.GetArrayFromImage(in_nii)[:50, :, :, :]

        n_volumes_to_discard = _get_vols_to_discard(in_nii)

        import pathlib  # Better path manipulation
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
        from rabies.preprocess_pkg.utils import run_command
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
        import os
        import SimpleITK as sitk
        from rabies.preprocess_pkg.utils import run_command

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
    import logging
    log = logging.getLogger('root')

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

        import os
        import SimpleITK as sitk
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

        import pathlib  # Better path manipulation
        split = pathlib.Path(self.inputs.name_source).name.rsplit(".nii")
        out_name = os.path.abspath(split[0]+'_slice_mc.nii.gz')
        sitk.WriteImage(resampled_timeseries, out_name)

        setattr(self, 'mc_corrected_bold', out_name)

        return runtime

    def _list_outputs(self):
        return {'mc_corrected_bold': getattr(self, 'mc_corrected_bold')}

def exec_applyTransforms(transforms, inverses, input_image, ref_image, output_image, mask=False):
    from rabies.preprocess_pkg.utils import run_command
    # tranforms is a list of transform files, set in order of call within antsApplyTransforms
    transform_string = ""
    for transform, inverse in zip(transforms, inverses):
        if transform=='NULL':
            continue
        elif bool(inverse):
            transform_string += f"-t [{transform},1] "
        else:
            transform_string += f"-t {transform} "

    if mask:
        interpolation = 'GenericLabel'
    else:
        interpolation = 'BSpline[5]'

    command = f'antsApplyTransforms -i {input_image} {transform_string}-n {interpolation} -r {ref_image} -o {output_image}'
    rc = run_command(command)
    if not os.path.isfile(output_image):
        raise ValueError(
            "Missing output image. Transform call failed: "+command)


class slice_applyTransformsInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="Input 4D EPI")
    ref_file = File(exists=True, mandatory=True,
                    desc="The reference 3D space to which the EPI will be warped.")
    transforms = traits.List(desc="List of transforms to apply to every volume.")
    inverses = traits.List(
        desc="Define whether some transforms must be inverse, with a boolean list where true defines inverse e.g.[0,1,0]")
    apply_motcorr = traits.Bool(
        default=True, desc="Whether to apply motion realignment.")
    motcorr_params = File(
        exists=True, desc="xforms from head motion estimation .csv file")
    resampling_dim = traits.Str(
        desc="Specification for the dimension of resampling.")
    rabies_data_type = traits.Int(mandatory=True,
                                  desc="Integer specifying SimpleITK data type.")


class slice_applyTransformsOutputSpec(TraitedSpec):
    out_files = traits.List(
        desc="warped images after the application of the transforms")


class slice_applyTransforms(BaseInterface):
    """
    This interface will apply a set of transforms to an input 4D EPI as well as motion realignment if specified.
    Susceptibility distortion correction can be applied through the provided transforms. A list of the corrected
    single volumes will be provided as outputs, and these volumes require to be merged to recover timeseries.
    """

    input_spec = slice_applyTransformsInputSpec
    output_spec = slice_applyTransformsOutputSpec

    def _run_interface(self, runtime):
        # resampling the reference image to the dimension of the EPI
        import SimpleITK as sitk
        import os
        from rabies.preprocess_pkg.utils import run_command

        img = sitk.ReadImage(self.inputs.in_file, self.inputs.rabies_data_type)

        if not self.inputs.resampling_dim == 'inputs_defined':
            shape = self.inputs.resampling_dim.split('x')
            spacing = (float(shape[0]), float(shape[1]), float(shape[2]))
        else:
            spacing = img.GetSpacing()[:3]
        resampled = resample_image_spacing(sitk.ReadImage(
            self.inputs.ref_file, self.inputs.rabies_data_type), spacing)
        sitk.WriteImage(resampled, 'resampled.nii.gz')

        # Splitting bold file into lists of single volumes
        [bold_volumes, num_volumes] = split_volumes(
            self.inputs.in_file, "bold_", self.inputs.rabies_data_type)

        if self.inputs.apply_motcorr:
            motcorr_params = self.inputs.motcorr_params
        ref_img = os.path.abspath('resampled.nii.gz')
        warped_volumes = []
        for x in range(0, num_volumes):
            warped_vol_fname = os.path.abspath(
                "deformed_volume" + str(x) + ".nii.gz")
            warped_volumes.append(warped_vol_fname)

            transforms = self.inputs.transforms
            inverses = self.inputs.inverses
            if self.inputs.apply_motcorr:
                command = f'antsMotionCorrStats -m {motcorr_params} -o motcorr_vol{x}.mat -t {x}'
                rc = run_command(command)

                transforms.append(f'motcorr_vol{x}.mat')
                inverses.append(0)

            exec_applyTransforms(transforms, inverses, bold_volumes[x], ref_img, warped_vol_fname, mask=False)
            # change image to specified data type
            sitk.WriteImage(sitk.ReadImage(warped_vol_fname,
                                           self.inputs.rabies_data_type), warped_vol_fname)

        setattr(self, 'out_files', warped_volumes)
        return runtime

    def _list_outputs(self):
        return {'out_files': getattr(self, 'out_files')}


def split_volumes(in_file, output_prefix, rabies_data_type):
    '''
    Takes as input a 4D .nii file and splits it into separate time series
    volumes by splitting on the 4th dimension
    '''
    import os
    in_nii = sitk.ReadImage(in_file, rabies_data_type)
    num_dimensions = len(in_nii.GetSize())
    num_volumes = in_nii.GetSize()[3]

    if num_dimensions != 4:
        raise ValueError("the input file must be of dimensions 4")

    volumes = []
    for x in range(0, num_volumes):
        data_slice = sitk.GetArrayFromImage(in_nii)[x, :, :, :]
        slice_fname = os.path.abspath(
            output_prefix + "vol" + str(x) + ".nii.gz")
        image_3d = copyInfo_3DImage(sitk.GetImageFromArray(
            data_slice, isVector=False), in_nii)
        sitk.WriteImage(image_3d, slice_fname)
        volumes.append(slice_fname)

    return [volumes, num_volumes]


class MergeInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiPath(File(exists=True), mandatory=True,
                              desc='input list of files to merge, listed in the order to merge')
    header_source = File(exists=True, mandatory=True,
                         desc='a Nifti file from which the header should be copied')
    clip_negative = traits.Bool(
        desc="Whether to clip out negative values.")
    rabies_data_type = traits.Int(mandatory=True,
                                  desc="Integer specifying SimpleITK data type.")


class MergeOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='output merged file')


class Merge(BaseInterface):
    """
    Takes a list of 3D Nifti files and merge them in the order listed.
    """

    input_spec = MergeInputSpec
    output_spec = MergeOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import SimpleITK as sitk

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.header_source).name.rsplit(".nii")

        sample_volume = sitk.ReadImage(
            self.inputs.in_files[0], self.inputs.rabies_data_type)
        length = len(self.inputs.in_files)
        shape = sitk.GetArrayFromImage(sample_volume).shape
        combined = np.zeros((length, shape[0], shape[1], shape[2]))

        i = 0
        for file in self.inputs.in_files:
            combined[i, :, :, :] = sitk.GetArrayFromImage(
                sitk.ReadImage(file, self.inputs.rabies_data_type))[:, :, :]
            i = i+1
        if (i != length):
            raise ValueError("Error occured with Merge.")
        combined_files = os.path.abspath(
            f"{filename_split[0]}_combined.nii.gz")

        if self.inputs.clip_negative:
            # clip potential negative values
            combined[(combined < 0).astype(bool)] = 0
        combined_image = sitk.GetImageFromArray(combined, isVector=False)

        # set metadata and affine for the newly constructed 4D image
        header_source = sitk.ReadImage(
            self.inputs.header_source, self.inputs.rabies_data_type)
        combined_image = copyInfo_4DImage(
            combined_image, sample_volume, header_source)
        sitk.WriteImage(combined_image, combined_files)

        setattr(self, 'out_file', combined_files)
        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file')}


def copyInfo_4DImage(image_4d, ref_3d, ref_4d):
    # function to establish metadata of an input 4d image. The ref_3d will provide
    # the information for the first 3 dimensions, and the ref_4d for the 4th.
    if ref_3d.GetMetaData('dim[0]') == '4':
        image_4d.SetSpacing(
            tuple(list(ref_3d.GetSpacing()[:3])+[ref_4d.GetSpacing()[3]]))
        image_4d.SetOrigin(
            tuple(list(ref_3d.GetOrigin()[:3])+[ref_4d.GetOrigin()[3]]))
        dim_3d = list(ref_3d.GetDirection())
        dim_4d = list(ref_4d.GetDirection())
        image_4d.SetDirection(
            tuple(dim_3d[:3]+[dim_4d[3]]+dim_3d[4:7]+[dim_4d[7]]+dim_3d[8:11]+dim_4d[11:]))
    elif ref_3d.GetMetaData('dim[0]') == '3':
        image_4d.SetSpacing(
            tuple(list(ref_3d.GetSpacing())+[ref_4d.GetSpacing()[3]]))
        image_4d.SetOrigin(
            tuple(list(ref_3d.GetOrigin())+[ref_4d.GetOrigin()[3]]))
        dim_3d = list(ref_3d.GetDirection())
        dim_4d = list(ref_4d.GetDirection())
        image_4d.SetDirection(
            tuple(dim_3d[:3]+[dim_4d[3]]+dim_3d[3:6]+[dim_4d[7]]+dim_3d[6:9]+dim_4d[11:]))
    else:
        raise ValueError('Unknown reference image dimensions.')
    return image_4d


def copyInfo_3DImage(image_3d, ref_3d):
    if ref_3d.GetMetaData('dim[0]') == '4':
        image_3d.SetSpacing(ref_3d.GetSpacing()[:3])
        image_3d.SetOrigin(ref_3d.GetOrigin()[:3])
        dim_3d = list(ref_3d.GetDirection())
        image_3d.SetDirection(tuple(dim_3d[:3]+dim_3d[4:7]+dim_3d[8:11]))
    elif ref_3d.GetMetaData('dim[0]') == '3':
        image_3d.SetSpacing(ref_3d.GetSpacing())
        image_3d.SetOrigin(ref_3d.GetOrigin())
        image_3d.SetDirection(ref_3d.GetDirection())
    else:
        raise ValueError('Unknown reference image dimensions.')
    return image_3d


def resample_image_spacing(image, output_spacing, resampling_interpolation='BSpline'):
    import SimpleITK as sitk
    import numpy as np
    dimension = 3
    identity = sitk.Transform(dimension, sitk.sitkIdentity)

    # Compute grid size based on the physical size and spacing.
    input_size = image.GetSize()
    sampling_ratio = np.asarray(image.GetSpacing())/np.asarray(output_spacing)
    output_size = [int(input_size[0]*sampling_ratio[0]), int(input_size[1]
                                                             * sampling_ratio[1]), int(input_size[2]*sampling_ratio[2])]

    # set default threader to platform to avoid freezing with MultiProc https://github.com/SimpleITK/SimpleITK/issues/1239
    sitk.ProcessObject_SetGlobalDefaultThreader('Platform')
    resampled_image = sitk.Resample(image, output_size, identity, sitk.sitkLinear,
                                    image.GetOrigin(), output_spacing, image.GetDirection())
    # clip potential negative values
    array = sitk.GetArrayFromImage(resampled_image)
    array[(array < 0).astype(bool)] = 0
    pos_resampled_image = sitk.GetImageFromArray(array, isVector=False)
    pos_resampled_image.CopyInformation(resampled_image)
    return pos_resampled_image


def convert_to_RAS(img_file, out_dir=None):
    # convert the input image to the RAS orientation convention
    import os
    import nibabel as nb
    img = nb.load(img_file)
    if nb.aff2axcodes(img.affine) == ('R', 'A', 'S'):
        return img_file
    else:
        import pathlib  # Better path manipulation
        split = pathlib.Path(img_file).name.rsplit(".nii")
        if out_dir is None:
            out_file = os.path.abspath(split[0]+'_RAS.nii.gz')
        else:
            out_file = out_dir+'/'+split[0]+'_RAS.nii.gz'
            if not os.path.isdir(out_dir):
                os.makedirs(out_dir)
        nb.as_closest_canonical(img).to_filename(out_file)
        return out_file


def resample_template(template_file, mask_file, file_list, spacing='inputs_defined', rabies_data_type=8):
    import os
    import SimpleITK as sitk
    import numpy as np
    from rabies.preprocess_pkg.utils import resample_image_spacing, run_command
    import logging
    log = logging.getLogger('root')

    if spacing == 'inputs_defined':
        file_list = list(np.asarray(file_list).flatten())
        img = sitk.ReadImage(file_list[0], rabies_data_type)
        low_dim = np.asarray(img.GetSpacing()[:3]).min()
        for file in file_list[1:]:
            img = sitk.ReadImage(file, rabies_data_type)
            new_low_dim = np.asarray(img.GetSpacing()[:3]).min()
            if new_low_dim < low_dim:
                low_dim = new_low_dim
        spacing = (low_dim, low_dim, low_dim)

        template_image = sitk.ReadImage(
            template_file, rabies_data_type)
        template_dim = template_image.GetSpacing()
        if np.asarray(template_dim[:3]).min() > low_dim:
            log.info("The template retains its original resolution.")
            return template_file, mask_file
    else:
        shape = spacing.split('x')
        spacing = (float(shape[0]), float(shape[1]), float(shape[2]))

    log.info(f"Resampling template to {spacing[0]}x{spacing[1]}x{spacing[2]}mm dimensions.")
    resampled_template = os.path.abspath("resampled_template.nii.gz")
    sitk.WriteImage(resample_image_spacing(sitk.ReadImage(
        template_file, rabies_data_type), spacing), resampled_template)

    # also resample the brain mask to ensure stable registrations further down
    resampled_mask = os.path.abspath("resampled_mask.nii.gz")
    command = f'antsApplyTransforms -d 3 -i {mask_file} -r {resampled_template} -o {resampled_mask} --verbose -n GenericLabel'
    rc = run_command(command)

    return resampled_template, resampled_mask


def run_command(command, verbose = False):
    # Run command and collect stdout
    # http://blog.endpoint.com/2015/01/getting-realtime-output-using-python.html # noqa
    import logging
    log = logging.getLogger('root')
    log.debug('Running: '+command)

    import subprocess
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
    if not out == '':
        if verbose:
            log.info(out)
        else:
            log.debug(out)
    if process.stderr is not None:
        if verbose:
            log.info(process.stderr)
        else:
            log.warning(process.stderr)
    rc = process.returncode
    return rc


def flatten_list(l):
    if type(l) == list:
        flattened = []
        for e in l:
            e = flatten_list(e)
            if type(e) == list:
                flattened += e
            else:
                flattened += [e]
        return flattened
    else:
        return l


def select_from_list(filename, filelist):
    filelist = flatten_list(filelist)
    import pathlib
    filename_template = pathlib.Path(filename).name.rsplit(".nii")[0]
    selected_file = None
    for file in filelist:
        if filename_template in file:
            if selected_file is None:
                selected_file = file
            else:
                raise ValueError(
                    f"Found duplicates for filename {filename_template}.")

    if selected_file is None:
        raise ValueError(f"No file associated with {filename_template} were found from list {filelist}.")
    else:
        return selected_file
