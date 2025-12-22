import os
import pathlib
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
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
    motion_estimation = pe.Node(sitkMotionCorr(level=opts.HMC_level,rabies_data_type=opts.data_type),
                         name='motion_correction', mem_gb=1.1*opts.scale_min_memory, n_procs=int(os.environ['RABIES_ITK_NUM_THREADS']))
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


class sitkMotionCorrInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='input BOLD time series')
    ref_file = File(exists=True, mandatory=True,
                    desc='ref file to realignment time series')
    level = traits.Int(desc="Select a level from 0 to 3, where each level is more stringent registration.")
    rabies_data_type = traits.Int(mandatory=True,
                                  desc="Integer specifying SimpleITK data type.")


class sitkMotionCorrOutputSpec(TraitedSpec):
    csv_params = File(
        exists=True, desc="csv files with the 6-parameters rigid body transformations")


class sitkMotionCorr(BaseInterface):
    """
    This interface performs motion realignment using simpleitk_timeseries_motion_correction. It takes a reference volume to which
    EPI volumes from the input 4D file are realigned based on a Rigid registration.
    """

    input_spec = sitkMotionCorrInputSpec
    output_spec = sitkMotionCorrOutputSpec

    def _run_interface(self, runtime):
        import SimpleITK as sitk
        from rabies.simpleitk_timeseries_motion_correction.motion import framewise_register_pair, write_transforms_to_csv
        moving = self.inputs.in_file
        ref_file = self.inputs.ref_file
        level = self.inputs.level
        filename_split = pathlib.Path(moving).name.rsplit(".nii")
        output_prefix = os.path.abspath(
            f'{filename_split[0]}_')

        transforms = framewise_register_pair(moving, ref_file, level=level, interpolation=sitk.sitkBSpline5)
        csv_param_file = output_prefix + f"moco.csv"
        write_transforms_to_csv(transforms, csv_param_file)
        setattr(self, 'csv_params', csv_param_file)

        return runtime

    def _list_outputs(self):
        return {'csv_params': getattr(self, 'csv_params')}


class EstimateMotionParamsInputSpec(BaseInterfaceInputSpec):
    motcorr_params = File(exists=True, mandatory=True,
                       desc="CSV file with the 6 rigid body parameters")
    boldspace_bold = File(exists=True, mandatory=True,
                      desc="BOLD image in its original space.")
    boldspace_brain_mask = File(exists=True, mandatory=True,
                      desc="Brain mask in the BOLD space.")


class EstimateMotionParamsOutputSpec(TraitedSpec):
    motion_params_csv = traits.File(desc="CSV file of motion parameters")
    FD_csv = traits.File(desc="CSV file with global framewise displacement.")


class EstimateMotionParams(BaseInterface):
    # motion_param_head_start
    """
    This interface generates estimations of framewise displacement, together with
    the expansion of the 6 motion parameters to include derivatives and squared parameters (Friston 24).
    Framewise displacement are computed within as follows:
        1. For each timepoint, the 3 Euler rotations and translations are converted to an affine matrix
        2. For each voxel within a brain mask representing the referential space post-motion realignment,
           the inverse transform is applied to generate a point pre-motion realignment.
        3. For framewise displacement, the distance is measured between the pre-correction points generated 
           from the current and the previous timeframe. Distance is measured in mm with the Euclidean distance.
    """
    # motion_param_head_end

    input_spec = EstimateMotionParamsInputSpec
    output_spec = EstimateMotionParamsOutputSpec

    def _run_interface(self, runtime):
        import numpy as np
        import os
        from rabies.utils import run_command
        import pathlib  # Better path manipulation
        from rabies.simpleitk_timeseries_motion_correction.framewise_displacement import calculate_framewise_displacement

        filename_split = pathlib.Path(self.inputs.boldspace_bold).name.rsplit(".nii")

        FD_csv = os.path.abspath(f"{filename_split[0]}_FD_file.csv")
        results = calculate_framewise_displacement(mask_file=self.inputs.boldspace_brain_mask, csv_file=self.inputs.motcorr_params, output_csv=FD_csv, verbose=False)

        motion_24,motion_24_header = motion_24_params(self.inputs.motcorr_params)
        # write into a .csv
        import pandas as pd
        df = pd.DataFrame(motion_24)
        df.columns = motion_24_header
        motion_params_csv = os.path.abspath(f"{filename_split[0]}_motion_params.csv")
        df.to_csv(motion_params_csv)

        setattr(self, 'FD_csv', FD_csv)
        setattr(self, 'motion_params_csv', motion_params_csv)

        return runtime

    def _list_outputs(self):
        return {'motion_params_csv': getattr(self, 'motion_params_csv'),
                'FD_csv': getattr(self, 'FD_csv')}


def motion_24_params(movpar_csv):
    '''
    motioncorr_24params: 6 head motion parameters, their temporal derivative, and the 12 corresponding squared items (Friston et al. 1996, Magn. Reson. Med.)
    '''
    import numpy as np
    import pandas as pd
    MOCO_df = pd.read_csv(movpar_csv)
    rigid_params = np.array([MOCO_df[param].values for param in ['EulerX','EulerY','EulerZ','TransX','TransY','TransZ']]).T

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
