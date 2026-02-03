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
            fields=['motcorr_params', 'hmc_qc_figure', 'hmc_qc_csv', 'hmc_qc_video', 'slice_corrected_bold']),
        name='outputnode')

    # Head motion correction (hmc)
    n_procs=int(os.environ['RABIES_ITK_NUM_THREADS'])
    motion_estimation = pe.Node(sitkMotionCorr(level=opts.HMC_level,rabies_data_type=opts.data_type),
                         name='motion_correction', mem_gb=1.1*opts.scale_min_memory, n_procs=n_procs)
    motion_estimation.plugin_args = {
        'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}
    
    workflow.connect([
        (inputnode, motion_estimation, [('ref_image', 'ref_file'),
                                        ('bold_file', 'in_file')]),
        (motion_estimation, outputnode, [
         ('csv_params', 'motcorr_params')]),
    ])

    if opts.gen_hmc_qc:
        hmc_qc_node = pe.Node(HMC_QC(figure_format=opts.figure_format),
                            name='hmc_qc_node', mem_gb=1.1*opts.scale_min_memory, n_procs=n_procs)

        workflow.connect([
            (inputnode, hmc_qc_node, [('ref_image', 'ref_file'),
                                            ('bold_file', 'in_file')]),
            (motion_estimation, hmc_qc_node, [
            ('csv_params', 'csv_params')]),
            (hmc_qc_node, outputnode, [
                ('out_figure', 'hmc_qc_figure'),
                ('out_csv', 'hmc_qc_csv'),
                ('video_file', 'hmc_qc_video'),
            ]),
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

        n_procs = int(os.environ['RABIES_ITK_NUM_THREADS']) if "RABIES_ITK_NUM_THREADS" in os.environ else os.cpu_count() # default to number of CPUs
        import SimpleITK as sitk
        from simpleitk_timeseries_motion_correction.motion import framewise_register_pair, write_transforms_to_csv
        moving = self.inputs.in_file
        ref_file = self.inputs.ref_file
        level = self.inputs.level
        filename_split = pathlib.Path(moving).name.rsplit(".nii")
        output_prefix = os.path.abspath(
            f'{filename_split[0]}_')

        transforms = framewise_register_pair(
            moving, 
            ref_file, 
            level=level, 
            interpolation=sitk.sitkBSpline5, 
            max_workers=n_procs)
        
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
        import pathlib  # Better path manipulation
        from simpleitk_timeseries_motion_correction.framewise_displacement import calculate_framewise_displacement

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


class HMC_QCInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='input BOLD time series')
    ref_file = File(exists=True, mandatory=True,
                    desc='ref file to realignment time series')
    csv_params = File(
        exists=True, desc="csv files with the 6-parameters rigid body transformations")
    figure_format = traits.Str(default='png',
        mandatory=True, exists=True, desc="png or svg")


class HMC_QCOutputSpec(TraitedSpec):
    out_figure = File(exists=True, desc="Output figure for QC visualization.")
    out_csv = File(exists=True, desc="Output CSV saving certain derivatives.")
    video_file = File(exists=True, desc="Output video showing timeseries pre/post correction.")


class HMC_QC(BaseInterface):
    """
    This interface pairs with the outputs from the head motion correction interface to generative a qualitative report
    that evaluate the performance of the correction.
    """

    input_spec = HMC_QCInputSpec
    output_spec = HMC_QCOutputSpec

    def _run_interface(self, runtime):
        import pandas as pd
        filename = pathlib.Path(self.inputs.in_file).name.rsplit(".nii")[0]
        figure_path = os.path.abspath(f'{filename}_HMC_QC.{self.inputs.figure_format}')
        csv_path = os.path.abspath(f'{filename}_derivatives.csv')

        n_procs = int(os.environ['RABIES_ITK_NUM_THREADS']) if "RABIES_ITK_NUM_THREADS" in os.environ else os.cpu_count() # default to number of CPUs

        derivatives_dict = HMC_derivatives(self.inputs.in_file, self.inputs.ref_file, self.inputs.csv_params, get_R2=False, n_procs=n_procs)

        # save some outputs to .csv
        key_l = ['D_Sc_preHMC', 'D_Sc_postHMC', 'mse_preHMC', 'mse_postHMC']
        pd.DataFrame(np.array([derivatives_dict[key].flatten() for key in key_l]).T, columns=key_l).to_csv(csv_path)

        fig = plot_motion_QC(derivatives_dict, self.inputs.ref_file, plot_R2=False)
        fig.savefig(figure_path, bbox_inches='tight')

        # write .webp file
        video_file = os.path.abspath(f'{filename}_HMC.webp')
        from simpleitk_timeseries_motion_correction.create_animation import main
        main(input_img=derivatives_dict['img_preHMC'], output_file=video_file, additional_input_imgs=[derivatives_dict['img_postHMC']], labels=['Before correction', 'After correction'], scale=2.0, fps=10)
        setattr(self, 'out_figure', figure_path)
        setattr(self, 'out_csv', csv_path)
        setattr(self, 'video_file', video_file)

        return runtime

    def _list_outputs(self):
        return {'out_figure': getattr(self, 'out_figure'),
                'out_csv': getattr(self, 'out_csv'),
                'video_file': getattr(self, 'video_file'),
                }


#############FUNCTIONS FOR HMC QC
import numpy as np
import SimpleITK as sitk
from rabies.utils import copyInfo_3DImage

def cosine_similarity(X,Y): 
    X_ = X.copy()
    X_ /= np.sqrt((X_ ** 2).sum(axis=0)) # first normalize from l2-norm
    Y_ = Y.copy()
    Y_ /= np.sqrt((Y_ ** 2).sum(axis=0)) # first normalize from l2-norm
    
    Sc = (X_.T.dot(Y_)) # then derive vector product across columns
    return Sc


def get_SD(img):
    array = sitk.GetArrayFromImage(img)
    std = array.std(axis=0)
    std_image = copyInfo_3DImage(
        sitk.GetImageFromArray(std, isVector=False), img)
    return std_image


def get_motion_R2(timeseries_img, translations,rotations):
    from rabies.analysis_pkg.analysis_math import closed_form
    timeseries_4d = sitk.GetArrayFromImage(timeseries_img)
    timeseries = timeseries_4d.reshape(timeseries_4d.shape[0],-1).astype(float)
    confounds_array = np.concatenate((translations,rotations), axis=1)
    
    X=confounds_array-confounds_array.mean(axis=0)
    Y=timeseries-timeseries.mean(axis=0)
    res = Y-X.dot(closed_form(X,Y))
    R2_spatial = 1-(res.var(axis=0)/Y.var(axis=0))

    R2_img = sitk.GetImageFromArray(R2_spatial.reshape(timeseries_4d.shape[1:]))
    R2_img = copyInfo_3DImage(R2_img, timeseries_img)
    return R2_img


def HMC_derivatives(in_img, in_ref, motcorr_params_file, n_procs=1, get_R2=False):
    import pandas as pd
    from simpleitk_timeseries_motion_correction.apply_transforms import read_transforms_from_csv, framewise_resample_volume

    # the input can be either a nifti file or an SITK image
    if isinstance(in_img, sitk.Image):
        img_preHMC = in_img
    elif os.path.isfile(in_img):
        img_preHMC = sitk.ReadImage(in_img)
    if isinstance(in_ref, sitk.Image):
        ref_img = in_ref
    elif os.path.isfile(in_ref):
        ref_img = sitk.ReadImage(in_ref)
    
    # prepare timeseries post-correction
    transforms = read_transforms_from_csv(motcorr_params_file)
    img_postHMC = framewise_resample_volume(
            img_preHMC, 
            ref_img, 
            transforms, 
            interpolation=sitk.sitkBSpline5, 
            clip_negative=True,
            extrapolator=False,
            max_workers=n_procs,
            )
    
    df = pd.read_csv(motcorr_params_file)
    translations = np.array([df[par] for par in ['TransX', 'TransY', 'TransZ']]).T
    rotations = np.array([df[par] for par in ['EulerX', 'EulerY', 'EulerZ']]).T

    num_volumes = img_preHMC.GetSize()[3]
    ref_img_array = sitk.GetArrayFromImage(ref_img).reshape(-1,1).astype(float)
    timeseries_preHMC = sitk.GetArrayFromImage(img_preHMC).reshape(num_volumes,-1).astype(float)
    timeseries_postHMC = sitk.GetArrayFromImage(img_postHMC).reshape(num_volumes,-1).astype(float)
    
    D_Sc_preHMC = 1-cosine_similarity(timeseries_preHMC.T,ref_img_array)
    D_Sc_postHMC = 1-cosine_similarity(timeseries_postHMC.T,ref_img_array)

    mse_preHMC = np.mean((timeseries_preHMC.T - ref_img_array)**2, axis=0) # taking mean square error
    mse_postHMC = np.mean((timeseries_postHMC.T - ref_img_array)**2, axis=0) # taking mean square error
    
    img_preHMC_SD = get_SD(img_preHMC)
    img_postHMC_SD = get_SD(img_postHMC)

    img_preHMC_R2 = get_motion_R2(img_preHMC, translations,rotations) if get_R2 else None
    img_postHMC_R2 = get_motion_R2(img_postHMC, translations,rotations) if get_R2 else None

    return {'img_preHMC':img_preHMC, 'img_postHMC':img_postHMC, 'translations':translations, 'rotations':rotations, 
            'D_Sc_preHMC':D_Sc_preHMC, 'D_Sc_postHMC':D_Sc_postHMC,'mse_preHMC':mse_preHMC, 'mse_postHMC':mse_postHMC, 
            'img_preHMC_SD':img_preHMC_SD, 'img_postHMC_SD':img_postHMC_SD, 'img_preHMC_R2':img_preHMC_R2, 'img_postHMC_R2':img_postHMC_R2}


def plot_motion_QC(derivatives_dict, ref_file, plot_R2=False):
    import matplotlib.pyplot as plt
    from rabies.visualization import plot_3d, otsu_scaling
    from matplotlib.gridspec import GridSpec
    nrows = 7
    if plot_R2:
        nrows +=3

    fig = plt.figure(figsize=(8*2,2*nrows), constrained_layout=True)
    gs = GridSpec(nrows, 3, figure=fig)
        
    ###PLOT timecourses
    ax_fused_l = [fig.add_subplot(gs[i, 0:3]) for i in range(4)]
    
    ax=ax_fused_l[0]
    ax.plot(derivatives_dict['translations'])
    ax.set_ylabel('mm', fontsize=15, color='white')
    ax.set_title('Translations', fontsize=20, color='white')
    ax.tick_params(labelsize=15)
    
    ax=ax_fused_l[1]
    ax.plot(derivatives_dict['rotations'])
    ax.set_ylabel('Euler angle (radians)', fontsize=15, color='white')
    ax.set_title('Rotations', fontsize=20, color='white')
    ax.tick_params(labelsize=15)
    
    ax=ax_fused_l[2]
    ax.plot(derivatives_dict['D_Sc_preHMC'])
    ax.plot(derivatives_dict['D_Sc_postHMC'])
    ax.set_ylabel('Cosine distance', fontsize=15, color='white')
    ax.set_title('Difference relative to reference volume', fontsize=20, color='white')
    ax.legend(['Before correction','After correction'], loc="upper right", fontsize=12)
    ax.tick_params(labelsize=15)
    
    ax=ax_fused_l[3]
    ax.plot(derivatives_dict['mse_preHMC'])
    ax.plot(derivatives_dict['mse_postHMC'])
    ax.set_ylabel('Mean square error', fontsize=15, color='white')
    ax.set_title('Difference relative to reference volume', fontsize=20, color='white')
    ax.legend(['Before correction','After correction'], loc="upper right", fontsize=12)
    ax.tick_params(labelsize=15)
    
    ###PLOT REFERENCE
    nrow=4
    ax_3d = [fig.add_subplot(gs[nrow+i, 0]) for i in range(3)]
    ax_3d[0].set_title('Reference image', fontsize=20, color='white')
    cbar_list = plot_3d(ax_3d,otsu_scaling(ref_file),fig=fig,vmin=0,vmax=1,cmap='gray', cbar=False)

    ###PLOT SD
    SD_array = sitk.GetArrayFromImage(derivatives_dict['img_preHMC_SD']).flatten()
    SD_array.sort()
    std_vmax=SD_array[int(len(SD_array)*0.99)]
    
    ax_3d = [fig.add_subplot(gs[nrow+i, 1]) for i in range(3)]
    ax_3d[0].set_title('Temporal standard deviation\n(before correction)', fontsize=20, color='white')
    cbar_list = plot_3d(ax_3d,derivatives_dict['img_preHMC_SD'],fig=fig,vmin=0,vmax=std_vmax,cmap='inferno', cbar=True)
    for cbar in cbar_list:
        cbar.ax.tick_params(labelsize=15)
    
    ax_3d = [fig.add_subplot(gs[nrow+i, 2]) for i in range(3)]
    ax_3d[0].set_title('Temporal standard deviation\n(after correction)', fontsize=20, color='white')
    cbar_list = plot_3d(ax_3d,derivatives_dict['img_postHMC_SD'],fig=fig,vmin=0,vmax=std_vmax,cmap='inferno', cbar=True)
    for cbar in cbar_list:
        cbar.ax.tick_params(labelsize=15)
    
    if plot_R2:
        nrow+=3

        ax_3d = [fig.add_subplot(gs[nrow+i, 0]) for i in range(3)]
        ax_3d[0].set_title('Reference image', fontsize=20, color='white')
        cbar_list = plot_3d(ax_3d,otsu_scaling(ref_file),fig=fig,vmin=0,vmax=1,cmap='gray', cbar=False)
        
        ax_3d = [fig.add_subplot(gs[nrow+i, 1]) for i in range(3)]
        ax_3d[0].set_title('R2 from 6 motion parameters\n(before correction)', fontsize=20, color='white')
        cbar_list = plot_3d(ax_3d,derivatives_dict['img_preHMC_R2'],fig=fig,vmin=0,vmax=1,cmap='inferno', cbar=True)
        for cbar in cbar_list:
            cbar.ax.tick_params(labelsize=15)
        
        ax_3d = [fig.add_subplot(gs[nrow+i, 2]) for i in range(3)]
        ax_3d[0].set_title('R2 from 6 motion parameters\n(after correction)', fontsize=20, color='white')
        cbar_list = plot_3d(ax_3d,derivatives_dict['img_postHMC_R2'],fig=fig,vmin=0,vmax=1,cmap='inferno', cbar=True)
        for cbar in cbar_list:
            cbar.ax.tick_params(labelsize=15)
        
    return fig