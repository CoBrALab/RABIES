from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function
from .utils import prep_CR
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

def init_confound_correction_wf(cr_opts, name="confound_correction_wf"):
    # confound_wf_head_start
    """
    This workflow applies the RABIES confound correction pipeline to preprocessed EPI timeseries. The correction steps are 
    orchestrated in line with recommendations from human litterature:   
    #1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)
    #2 - If --match_number_timepoints is selected, each scan is matched to the defined minimum_timepoint number of frames.
    #4 - Linear/Quadratic detrending of fMRI timeseries and nuisance regressors
    #4 - Apply ICA-AROMA.
    #5 - If frequency filtering and frame censoring are applied, simulate data in censored timepoints using the Lomb-Scargle periodogram, 
         as suggested in Power et al. (2014, Neuroimage), for both the fMRI timeseries and nuisance regressors prior to filtering.
    #6 - As recommended in Lindquist et al. (2019, Human brain mapping), make the nuisance regressors orthogonal
         to the temporal frequency filter.
    #7 - Apply highpass and/or lowpass filtering on the fMRI timeseries (with simulated timepoints).
    #8 - Re-apply the frame censoring mask onto filtered fMRI timeseries and nuisance regressors, taking out the
         simulated timepoints. Edge artefacts from frequency filtering can also be removed as recommended in Power et al. (2014, Neuroimage).
    #9 - Apply confound regression using the selected nuisance regressors.
    #10 - Scaling of timeseries variance.
    #11 - Apply Gaussian spatial smoothing.
    
    References:
        Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic 
            correlations in functional connectivity MRI networks arise from subject motion. Neuroimage, 59(3), 2142-2154.
        Power, J. D., Mitra, A., Laumann, T. O., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2014). Methods to detect, 
            characterize, and remove motion artifact in resting state fMRI. Neuroimage, 84, 320-341.
        Lindquist, M. A., Geuter, S., Wager, T. D., & Caffo, B. S. (2019). Modular preprocessing pipelines can reintroduce 
            artifacts into fMRI data. Human brain mapping, 40(8), 2358-2376.

    Workflow:
        parameters
            cr_opts: command line interface parameters from confound_correction

        inputs
            bold_file: preprocessed EPI timeseries
            brain_mask: brain mask overlapping with EPI timeseries
            csf_mask: CSF mask overlapping with EPI timeseries
            motion_params_csv: CSV file with motion regressors
            FD_file: CSV file with the framewise displacement

        outputs
            cleaned_path: the cleaned EPI timeseries
            aroma_out: folder with outputs from ICA-AROMA
            VE_file: variance explained (R^2) from confound regression at each voxel
            STD_file: standard deviation on the cleaned EPI timeseries
            CR_STD_file: standard deviation on the confound timeseries modelled during confound regression
            random_CR_STD_file_path: variance fitted by random regressors during confound regression
            corrected_CR_STD_file_path: CR_STD_file after substracting the variance fitted by random
                regressors.
            frame_mask_file: CSV file which records which frame were censored
            CR_data_dict: dictionary object storing extra data computed during confound correction
    """
    # confound_wf_head_end

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=[
                        'bold_file', 'brain_mask', 'WM_mask', 'CSF_mask', 'vascular_mask', 'motion_params_csv', 'FD_file', 'raw_input_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=[
                         'cleaned_path', 'aroma_out', 'VE_file', 'STD_file', 'CR_STD_file', 
                         'random_CR_STD_file_path', 'corrected_CR_STD_file_path', 'frame_mask_file', 'CR_data_dict']), name='outputnode')
    clean_image_node = pe.Node(CleanImage(cr_opts=cr_opts),
                           name='clean_image', mem_gb=5*cr_opts.scale_min_memory) # 5X memory as the timeseries is expanded into many arrays

    prep_CR_node = pe.Node(Function(input_names=['bold_file', 'motion_params_csv', 'FD_file', 'cr_opts'],
                                              output_names=['data_dict'],
                                              function=prep_CR),
                                     name='prep_CR', mem_gb=1)
    prep_CR_node.inputs.cr_opts = cr_opts

    workflow.connect([
        (inputnode, prep_CR_node, [
            ("bold_file", "bold_file"),
            ("motion_params_csv", "motion_params_csv"),
            ("FD_file", "FD_file"),
            ]),
        (inputnode, clean_image_node, [
            ("bold_file", "bold_file"),
            ("brain_mask", "brain_mask_file"),
            ("WM_mask", "WM_mask_file"),
            ("CSF_mask", "CSF_mask_file"),
            ("vascular_mask", "vascular_mask_file"),
            ("raw_input_file", "raw_input_file"),
            ]),
        (prep_CR_node, clean_image_node, [
            ("data_dict", "data_dict"),
            ]),
        (clean_image_node, outputnode, [
            ("cleaned_path", "cleaned_path"),
            ("VE_file_path", "VE_file"),
            ("STD_file_path", "STD_file"),
            ("CR_STD_file_path", "CR_STD_file"),
            ("random_CR_STD_file_path", "random_CR_STD_file_path"),
            ("corrected_CR_STD_file_path", "corrected_CR_STD_file_path"),
            ("frame_mask_file", "frame_mask_file"),
            ("data_dict", "CR_data_dict"),
            ("aroma_out", "aroma_out"),
            ]),
        ])

    return workflow


class CleanImageInputSpec(BaseInterfaceInputSpec):
    raw_input_file = File(exists=True, mandatory=True,
                      desc="The raw EPI scan before preprocessing.")
    bold_file = File(exists=True, mandatory=True,
                      desc="Timeseries to denoise.")
    data_dict = traits.Dict(
        exists=True, mandatory=True, desc="Dictionary with extra inputs.")
    brain_mask_file = File(exists=True, mandatory=True,
                      desc="Brain mask.")
    WM_mask_file = traits.Any(mandatory=True,
                      desc="WM mask.")
    CSF_mask_file = traits.Any(mandatory=True,
                      desc="CSF mask.")
    vascular_mask_file = traits.Any(mandatory=True,
                      desc="vascular mask.")
    cr_opts = traits.Any(
        exists=True, mandatory=True, desc="Processing specs.")

class CleanImageOutputSpec(TraitedSpec):
    cleaned_path = File(exists=True, mandatory=True,
                      desc="Cleaned timeseries.")
    VE_file_path = File(exists=True, mandatory=True,
                      desc="Variance explained map from confound regression.")
    STD_file_path = File(exists=True, mandatory=True,
                      desc="Temporal standard deviation map after confound correction, prior to standardization.")
    CR_STD_file_path = File(exists=True, mandatory=True,
                      desc="Temporal standard deviation on predicted confound timeseries.")
    random_CR_STD_file_path = File(exists=False, mandatory=False,
                      desc="Temporal standard deviation on predicted confound timeseries from random regressors.")
    corrected_CR_STD_file_path = File(exists=False, mandatory=False,
                      desc="Same as CR_STD_file_path, but the variance explained by random regressors was substracted.")
    frame_mask_file = File(exists=True, mandatory=True,
                      desc="Frame mask from temporal censoring.")
    data_dict = traits.Any(
        desc="A dictionary with key outputs.")
    aroma_out = traits.Any(
        desc="Output directory from ICA-AROMA.")

class CleanImage(BaseInterface):
    '''
    Apply a flexible confound correction pipeline in line with recommendations from
    human litterature. 
    
    #1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)
    #2 - If --match_number_timepoints is selected, each scan is matched to the defined minimum_timepoint number of frames.
    #4 - Linear/Quadratic detrending of fMRI timeseries and nuisance regressors
    #4 - Apply ICA-AROMA.
    #5 - If frequency filtering and frame censoring are applied, simulate data in censored timepoints using the Lomb-Scargle periodogram, 
         as suggested in Power et al. (2014, Neuroimage), for both the fMRI timeseries and nuisance regressors prior to filtering.
    #6 - As recommended in Lindquist et al. (2019, Human brain mapping), make the nuisance regressors orthogonal
         to the temporal frequency filter.
    #7 - Apply highpass and/or lowpass filtering on the fMRI timeseries (with simulated timepoints).
    #8 - Re-apply the frame censoring mask onto filtered fMRI timeseries and nuisance regressors, taking out the
         simulated timepoints. Edge artefacts from frequency filtering can also be removed as recommended in Power et al. (2014, Neuroimage).
    #9 - If selected, compute the WM/CSF/vascular signal or aCompCorr and add to list of regressors. This is computed post-AROMA and filtering to 
         minimize re-introduction of previously corrected signal fluctuations.
    #10 - Apply confound regression using the selected nuisance regressors.
    #11 - Scaling of timeseries.
    #12 - Apply Gaussian spatial smoothing.
    
    References:
        
        Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). 
        Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage, 59(3), 2142-2154.
        
        Power, J. D., Mitra, A., Laumann, T. O., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2014). 
        Methods to detect, characterize, and remove motion artifact in resting state fMRI. Neuroimage, 84, 320-341.
        
        Lindquist, M. A., Geuter, S., Wager, T. D., & Caffo, B. S. (2019). 
        Modular preprocessing pipelines can reintroduce artifacts into fMRI data. Human brain mapping, 40(8), 2358-2376.
    '''

    input_spec = CleanImageInputSpec
    output_spec = CleanImageOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import pandas as pd
        import SimpleITK as sitk
        from rabies.utils import recover_3D,recover_4D
        from rabies.confound_correction_pkg.utils import temporal_censoring,lombscargle_fill, exec_ICA_AROMA,butterworth, phase_randomized_regressors, smooth_image, remove_trend,compute_signal_regressors
        from rabies.analysis_pkg.analysis_functions import closed_form

        ### set null returns in case the workflow is interrupted
        empty_img = sitk.GetImageFromArray(np.empty([1,1]))
        empty_file = os.path.abspath('empty.nii.gz')
        sitk.WriteImage(empty_img, empty_file)

        setattr(self, 'cleaned_path', empty_file)
        setattr(self, 'VE_file_path', empty_file)
        setattr(self, 'STD_file_path', empty_file)
        setattr(self, 'CR_STD_file_path', empty_file)
        setattr(self, 'random_CR_STD_file_path', empty_file)
        setattr(self, 'corrected_CR_STD_file_path', empty_file)
        setattr(self, 'frame_mask_file', empty_file)
        setattr(self, 'data_dict', empty_file)
        setattr(self, 'aroma_out', empty_file)
        ###

        cr_opts = self.inputs.cr_opts
        data_dict = self.inputs.data_dict

        from nipype import logging
        nipype_log = logging.getLogger('nipype.workflow')

        cleaning_out_dict = clean_image(
            bold_file = self.inputs.bold_file,
            brain_mask_file = self.inputs.brain_mask_file,
            WM_mask_file = self.inputs.WM_mask_file,
            CSF_mask_file = self.inputs.CSF_mask_file,
            vascular_mask_file = self.inputs.vascular_mask_file,
            FD_trace=data_dict['FD_trace'],
            confounds_array=data_dict['confounds_array'],
            motion_params_csv=data_dict['motion_params_csv'],
            time_range=data_dict['time_range'],
            confounds_6rigid_array=data_dict['confounds_6rigid_array'],
            FD_censoring=cr_opts.frame_censoring['FD_censoring'], 
            FD_threshold=cr_opts.frame_censoring['FD_threshold'], 
            DVARS_censoring=cr_opts.frame_censoring['DVARS_censoring'], 
            minimum_timepoint=cr_opts.frame_censoring['minimum_timepoint'],
            TR=cr_opts.TR,
            detrending_order=cr_opts.detrending['order'], 
            detrending_time_interval=cr_opts.detrending['time_interval'], 
            apply_ica_aroma=cr_opts.ica_aroma['apply'], 
            ica_aroma_dim=cr_opts.ica_aroma['dim'], 
            ica_aroma_random_seed=cr_opts.ica_aroma['random_seed'],
            match_number_timepoints=cr_opts.match_number_timepoints, 
            highpass=cr_opts.highpass, 
            lowpass=cr_opts.lowpass, 
            edge_cutoff=cr_opts.edge_cutoff,
            nuisance_regressors = cr_opts.nuisance_regressors, 
            generate_CR_null=cr_opts.generate_CR_null,
            scale_variance_voxelwise=cr_opts.scale_variance_voxelwise,
            image_scaling=cr_opts.image_scaling,
            smoothing_filter=cr_opts.smoothing_filter,
            nipype_log=nipype_log,
            )

        if cleaning_out_dict is None:
            return runtime
        
        setattr(self, 'cleaned_path', cleaning_out_dict['cleaned_path'])
        setattr(self, 'VE_file_path', cleaning_out_dict['VE_file_path'])
        setattr(self, 'STD_file_path', cleaning_out_dict['STD_file_path'])
        setattr(self, 'CR_STD_file_path', cleaning_out_dict['CR_STD_file_path'])
        setattr(self, 'random_CR_STD_file_path', cleaning_out_dict['random_CR_STD_file_path'])
        setattr(self, 'corrected_CR_STD_file_path', cleaning_out_dict['corrected_CR_STD_file_path'])
        setattr(self, 'frame_mask_file', cleaning_out_dict['frame_mask_file'])
        setattr(self, 'data_dict', cleaning_out_dict)
        if cleaning_out_dict['aroma_out'] is not None:
            setattr(self, 'aroma_out', cleaning_out_dict['aroma_out'])

        return runtime

    def _list_outputs(self):
        return {'cleaned_path': getattr(self, 'cleaned_path'),
                'VE_file_path': getattr(self, 'VE_file_path'),
                'STD_file_path': getattr(self, 'STD_file_path'),
                'CR_STD_file_path': getattr(self, 'CR_STD_file_path'),
                'random_CR_STD_file_path': getattr(self, 'random_CR_STD_file_path'),
                'corrected_CR_STD_file_path': getattr(self, 'corrected_CR_STD_file_path'),
                'frame_mask_file': getattr(self, 'frame_mask_file'),
                'data_dict': getattr(self, 'data_dict'),
                'aroma_out': getattr(self, 'aroma_out'),
                }
    

def clean_image(bold_file, brain_mask_file, WM_mask_file, CSF_mask_file, vascular_mask_file,
                FD_trace, confounds_array, motion_params_csv, time_range, confounds_6rigid_array, # replacing data_dict
                FD_censoring=False, FD_threshold=0.05, DVARS_censoring=False, minimum_timepoint=3, TR='auto', # replacing cr_opts
                detrending_order=1, detrending_time_interval='all', 
                apply_ica_aroma=False, ica_aroma_dim=0, ica_aroma_random_seed=1,
                match_number_timepoints=False, highpass=None, lowpass=None, edge_cutoff=0,
                nuisance_regressors = [], generate_CR_null=False,
                scale_variance_voxelwise=False,image_scaling='grand_mean_scaling',
                smoothing_filter=None,
                nipype_log=None,
                ):
    import os
    import numpy as np
    import pandas as pd
    import SimpleITK as sitk
    from rabies.utils import recover_3D,recover_4D
    from rabies.confound_correction_pkg.utils import temporal_censoring,lombscargle_fill, exec_ICA_AROMA,butterworth, phase_randomized_regressors, smooth_image, remove_trend,compute_signal_regressors, get_DVARS
    from rabies.analysis_pkg.analysis_functions import closed_form

    cr_out = os.getcwd()
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

    volume_idx = sitk.GetArrayFromImage(sitk.ReadImage(brain_mask_file)).astype(bool)

    from rabies.utils import get_sitk_header
    header=get_sitk_header(bold_file)
    timeseries = sitk.GetArrayFromImage(sitk.ReadImage(bold_file, sitk.sitkFloat32))[:,volume_idx] # read directly as a 2D array
    timeseries = timeseries[time_range,:]

    if TR=='auto':
        TR = float(header.GetSpacing()[3])
    else:
        TR = float(TR)

    '''
    #1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)
    '''

    # compute the DVARS before denoising
    DVARS_trace = get_DVARS(timeseries)

    frame_mask = temporal_censoring(FD_trace, 
            FD_censoring, FD_threshold, DVARS_trace, DVARS_censoring, minimum_timepoint)
    if frame_mask is None:
        return None

    '''
    #2 - If --match_number_timepoints is selected, each scan is matched to the defined minimum_timepoint number of frames.
    '''
    if match_number_timepoints:
        if (not highpass is None) or (not lowpass is None):
            # if frequency filtering is applied, avoid selecting timepoints that would be removed with --edge_cutoff
            num_cut = int(edge_cutoff/TR)
            if not num_cut==0:
                frame_mask[:num_cut]=0
                frame_mask[-num_cut:]=0

                if frame_mask.sum()<int(minimum_timepoint):
                    if nipype_log:
                        nipype_log.warning(f"CONFOUND CORRECTION LEFT LESS THAN {str(minimum_timepoint)} VOLUMES. THIS SCAN WILL BE REMOVED FROM FURTHER PROCESSING.")
                    return None

        # randomly shuffle indices that haven't been censored, then remove an extra subset above --minimum_timepoint
        num_timepoints = len(frame_mask)
        time_idx=np.array(range(num_timepoints))
        perm = np.random.permutation(time_idx[frame_mask])
        # selecting the subset of extra timepoints, and censoring them
        subset_idx = perm[minimum_timepoint:]
        frame_mask[subset_idx]=0
        # keep track of the original number of timepoints for tDOF estimation, to evaluate latter if the correction was succesful
        number_extra_timepoints = len(subset_idx)
    else:
        number_extra_timepoints = 0

    timeseries = timeseries[frame_mask]
    confounds_array = confounds_array[frame_mask]

    '''
    #3 - Linear/Quadratic detrending of fMRI timeseries and nuisance regressors
    '''
    # apply detrending, after censoring
    # save grand mean prior to detrending
    timeseries_ = remove_trend(timeseries, frame_mask, order = detrending_order, time_interval = detrending_time_interval, keep_intercept=True)
    grand_mean = timeseries_.mean()
    voxelwise_mean = timeseries_.mean(axis=0)
    del timeseries_ # free space

    timeseries = remove_trend(timeseries, frame_mask, order = detrending_order, time_interval = detrending_time_interval, keep_intercept=False)
    confounds_array = remove_trend(confounds_array, frame_mask, order = detrending_order, time_interval = detrending_time_interval, keep_intercept=False)

    '''
    #4 - Apply ICA-AROMA.
    '''
    if apply_ica_aroma:
        # write intermediary output files for timeseries and 6 rigid body parameters
        inFile = f'{cr_out}/{filename_split[0]}_aroma_input.nii.gz'
        sitk.WriteImage(recover_4D(brain_mask_file, timeseries, bold_file), inFile)

        confounds_6rigid_array=confounds_6rigid_array[frame_mask,:]
        confounds_6rigid_array = remove_trend(confounds_6rigid_array, frame_mask, order = detrending_order, time_interval = detrending_time_interval, keep_intercept=False) # apply detrending to the confounds too
        df = pd.DataFrame(confounds_6rigid_array)
        df.columns = ['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3']
        mc_file = f'{cr_out}/{filename_split[0]}_aroma_input.csv'
        df.to_csv(mc_file)

        cleaned_file, aroma_out = exec_ICA_AROMA(inFile, mc_file, brain_mask_file, CSF_mask_file, TR, ica_aroma_dim, random_seed=ica_aroma_random_seed)
        timeseries = sitk.GetArrayFromImage(sitk.ReadImage(cleaned_file, sitk.sitkFloat32))[:,volume_idx] # read directly as a 2D array
    else:
        aroma_out = None

    if (not highpass is None) or (not lowpass is None):
        '''
        #5 - If frequency filtering and frame censoring are applied, simulate data in censored timepoints using the Lomb-Scargle periodogram, 
            as suggested in Power et al. (2014, Neuroimage), for both the fMRI timeseries and nuisance regressors prior to filtering.
        '''
        timeseries = lombscargle_fill(x=timeseries,time_step=TR,time_mask=frame_mask)
        confounds_array = lombscargle_fill(x=confounds_array,time_step=TR,time_mask=frame_mask)
        ### arrays are now interpolated

        '''
        #6 - As recommended in Lindquist et al. (2019, Human brain mapping), make the nuisance regressors orthogonal
            to the temporal filter.
        '''
        confounds_array = butterworth(confounds_array, TR=TR,
                                high_pass=highpass, low_pass=lowpass)

        '''
        #7 - Apply highpass and/or lowpass filtering on the fMRI timeseries (with simulated timepoints).
        '''

        timeseries = butterworth(timeseries, TR=TR,
                                high_pass=highpass, low_pass=lowpass)

        # correct for edge effects of the filters
        num_cut = int(edge_cutoff/TR)
        if len(frame_mask)<2*num_cut:
            raise ValueError(f"The timeseries are too short to remove {edge_cutoff}sec of data at each edge.")

        if not num_cut==0:
            frame_mask[:num_cut]=0
            frame_mask[-num_cut:]=0


        '''
        #8 - Re-apply the frame censoring mask onto filtered fMRI timeseries and nuisance regressors, taking out the
            simulated timepoints. Edge artefacts from frequency filtering can also be removed as recommended in Power et al. (2014, Neuroimage).
        '''
        # re-apply the masks to take out simulated data points, and take off the edges
        timeseries = timeseries[frame_mask]
        confounds_array = confounds_array[frame_mask]
    
    if frame_mask.sum()<int(minimum_timepoint):
        if nipype_log:
            nipype_log.warning(f"CONFOUND CORRECTION LEFT LESS THAN {str(minimum_timepoint)} VOLUMES. THIS SCAN WILL BE REMOVED FROM FURTHER PROCESSING.")
        return None

    '''
    #9 - If selected, compute the WM/CSF/vascular signal or aCompCorr and add to list of regressors. This is computed post-AROMA and filtering to 
        minimize re-introduction of previously corrected signal fluctuations.
    '''    
    # indices for each mask are first loaded in vector format that matches the timeseries array
    [brain_mask_idx, WM_mask_idx, CSF_mask_idx, vascular_mask_idx] = [sitk.GetArrayFromImage(sitk.ReadImage(mask_file)).astype(bool)[volume_idx] for mask_file in [brain_mask_file, WM_mask_file, CSF_mask_file, vascular_mask_file]]
    regressors_array = compute_signal_regressors(timeseries, nuisance_regressors, brain_mask_idx, WM_mask_idx, CSF_mask_idx, vascular_mask_idx)
    confounds_array = np.append(confounds_array,regressors_array,axis=1)

    '''
    #10 - Apply confound regression using the selected nuisance regressors.
    '''
    # voxels that have a NaN value are set to 0
    nan_voxels = np.isnan(timeseries).sum(axis=0)>1
    timeseries[:,nan_voxels] = 0

    # estimate the VE from the CR selection, or 6 rigid motion parameters if no CR is applied
    try:
        predicted = confounds_array.dot(closed_form(confounds_array,timeseries))
        residuals = timeseries-predicted
    except:
        if nipype_log:
            nipype_log.warning("SINGULAR MATRIX ERROR DURING CONFOUND REGRESSION. THIS SCAN WILL BE REMOVED FROM FURTHER PROCESSING.")
        return None

    # estimates of variance explained
    VE_total_ratio = 1-(residuals.var()/timeseries.var())
    VE_spatial = 1-(residuals.var(axis=0)/timeseries.var(axis=0))
    VE_temporal = 1-(residuals.var(axis=1)/timeseries.var(axis=1))

    if generate_CR_null:
        # estimate the fit from CR with randomized regressors as in BRIGHT AND MURPHY 2015
        randomized_confounds_array = phase_randomized_regressors(confounds_array, frame_mask, TR=TR)
        predicted_random = randomized_confounds_array.dot(closed_form(randomized_confounds_array,timeseries))
    else:
        predicted_random = np.empty([0,timeseries.shape[1]])

    if len(nuisance_regressors) > 0:
        # if confound regression is applied
        timeseries = residuals
    del residuals

    '''
    #11 - Scaling of timeseries.
    '''
    if scale_variance_voxelwise: # homogenize the variability distribution while preserving the same total variance
        if image_scaling=='voxelwise_standardization' or image_scaling=='voxelwise_mean':
            raise ValueError(f"Can't select --scale_variance_voxelwise with --image_scaling {image_scaling}.")
        else:
            temporal_std = timeseries.std(axis=0)
            total_std = timeseries.std()

            # homogenize variance across voxels
            timeseries /= temporal_std
            nan_voxels = np.isnan(timeseries).sum(axis=0)>1
            timeseries[:,nan_voxels] = 0
            scaled_total_std = timeseries.std()
            # rescale to preserve total variance
            timeseries = (timeseries/scaled_total_std)*total_std

            # the same scaling estimated from BOLD is applied to CR estimates to preserve accurate ratios of variance explained per voxel
            scaled_list = []
            for scaled_timeseries in [predicted,predicted_random]:
                scaled_timeseries /= temporal_std
                nan_voxels = np.isnan(scaled_timeseries).sum(axis=0)>1
                scaled_timeseries[:,nan_voxels] = 0
                scaled_timeseries = (scaled_timeseries/scaled_total_std)*total_std
                scaled_list.append(scaled_timeseries)
            [predicted,predicted_random] = scaled_list


    if image_scaling=='global_variance':
        scaling_factor = timeseries.std()
    elif image_scaling=='grand_mean_scaling':
        scaling_factor = grand_mean/100 # we scale BOLD in % fluctuations, hence dividing by 100
    elif image_scaling=='voxelwise_standardization':
        # each voxel is scaled according to its STD
        scaling_factor = timeseries.std(axis=0) 
    elif image_scaling=='voxelwise_mean':
        scaling_factor = voxelwise_mean/100 # we scale BOLD in % fluctuations, hence dividing by 100
    else:
        scaling_factor = 1

    # scale fMRI timeseries, as well as the CR prediction to keep consistent scaling
    scaled_list = []
    for scaled_timeseries in [timeseries,predicted,predicted_random]:
        scaled_timeseries /= scaling_factor
        nan_voxels = np.isnan(scaled_timeseries).sum(axis=0)>1
        scaled_timeseries[:,nan_voxels] = 0
        scaled_list.append(scaled_timeseries)
    [timeseries,predicted,predicted_random] = scaled_list

    # after variance scaling, compute the variability estimates
    temporal_std = timeseries.std(axis=0)
    predicted_std = predicted.std(axis=0)
    predicted_time = np.sqrt((predicted.T**2).mean(axis=0))
    predicted_global_std = predicted.std()

    if generate_CR_null:
        predicted_random_std = predicted_random.std(axis=0)

        # here we correct the previous STD estimates by substrating the variance explained by that of the overfitting with random regressors
        var_dif = predicted.var(axis=0) - predicted_random.var(axis=0)
        var_dif[var_dif<0] = 0 # when there's more variance explained in random regressors, set variance explained to 0
        corrected_predicted_std = np.sqrt(var_dif)
        del predicted_random
    del predicted

    # apply the frame mask to FD trace/DVARS to prepare outputs from the function
    DVARS_trace = DVARS_trace[frame_mask]
    FD_trace = FD_trace[frame_mask]

    # calculate temporal degrees of freedom left after confound correction
    num_timepoints = frame_mask.sum()
    if apply_ica_aroma:
        aroma_rm = (pd.read_csv(f'{aroma_out}/classification_overview.txt', sep='\t')['Motion/noise']).sum()
    else:
        aroma_rm = 0
    num_regressors = confounds_array.shape[1]
    tDOF = num_timepoints - (aroma_rm+num_regressors) + number_extra_timepoints

    # save output files
    VE_spatial_map = recover_3D(brain_mask_file, VE_spatial)
    STD_spatial_map = recover_3D(brain_mask_file, temporal_std)
    CR_STD_spatial_map = recover_3D(brain_mask_file, predicted_std)
    del VE_spatial, temporal_std, predicted_std
    timeseries_img = recover_4D(brain_mask_file, timeseries, bold_file)
    del timeseries
    if generate_CR_null:
        random_CR_STD_spatial_map = recover_3D(brain_mask_file, predicted_random_std)
        corrected_CR_STD_spatial_map = recover_3D(brain_mask_file, corrected_predicted_std)
        del predicted_random_std, corrected_predicted_std


    if smoothing_filter is not None:
        '''
        #12 - Apply Gaussian spatial smoothing.
        '''
        import nibabel as nb
        affine = nb.load(bold_file).affine[:3,:3] # still not sure how to match nibabel's affine reliably
        mask_img = sitk.ReadImage(brain_mask_file, sitk.sitkFloat32)
        timeseries_img = smooth_image(timeseries_img, affine, smoothing_filter, mask_img)

    cleaned_path = cr_out+'/'+filename_split[0]+'_cleaned.nii.gz'
    sitk.WriteImage(timeseries_img, cleaned_path)
    VE_file_path = cr_out+'/'+filename_split[0]+'_VE_map.nii.gz'
    sitk.WriteImage(VE_spatial_map, VE_file_path)
    STD_file_path = cr_out+'/'+filename_split[0]+'_STD_map.nii.gz'
    sitk.WriteImage(STD_spatial_map, STD_file_path)
    CR_STD_file_path = cr_out+'/'+filename_split[0]+'_CR_STD_map.nii.gz'
    sitk.WriteImage(CR_STD_spatial_map, CR_STD_file_path)
    frame_mask_file = cr_out+'/'+filename_split[0]+'_frame_censoring_mask.csv'
    pd.DataFrame(frame_mask).to_csv(frame_mask_file, index=False, header=['False = Masked Frames'])

    if generate_CR_null:
        random_CR_STD_file_path = cr_out+'/'+filename_split[0]+'_random_CR_STD_map.nii.gz'
        sitk.WriteImage(random_CR_STD_spatial_map, random_CR_STD_file_path)
        corrected_CR_STD_file_path = cr_out+'/'+filename_split[0]+'_corrected_CR_STD_map.nii.gz'
        sitk.WriteImage(corrected_CR_STD_spatial_map, corrected_CR_STD_file_path)
    else:
        random_CR_STD_file_path=''
        corrected_CR_STD_file_path=''

    return {
        'cleaned_path':cleaned_path, 'VE_file_path':VE_file_path, 'STD_file_path':STD_file_path, 'CR_STD_file_path':CR_STD_file_path, 'random_CR_STD_file_path':random_CR_STD_file_path, 'corrected_CR_STD_file_path':corrected_CR_STD_file_path, 'frame_mask_file':frame_mask_file, 'aroma_out':aroma_out,
        'TR':TR, 'FD_trace':FD_trace, 'DVARS':DVARS_trace, 'time_range':time_range, 'frame_mask':frame_mask, 'confounds_array':confounds_array, 'VE_temporal':VE_temporal, 'motion_params_csv':motion_params_csv, 'predicted_time':predicted_time, 'tDOF':tDOF, 'CR_global_std':predicted_global_std, 'VE_total_ratio':VE_total_ratio, 'voxelwise_mean':voxelwise_mean,
        }

