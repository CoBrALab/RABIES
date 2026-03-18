from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
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

    workflow.connect([
        (inputnode, clean_image_node, [
            ("bold_file", "bold_file"),
            ("brain_mask", "brain_mask_file"),
            ("WM_mask", "WM_mask_file"),
            ("CSF_mask", "CSF_mask_file"),
            ("vascular_mask", "vascular_mask_file"),
            ("raw_input_file", "raw_input_file"),
            ("motion_params_csv", "motion_params_csv"),
            ("FD_file", "FD_file"),
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
    brain_mask_file = File(exists=True, mandatory=True,
                      desc="Brain mask.")
    WM_mask_file = traits.Any(mandatory=True,
                      desc="WM mask.")
    CSF_mask_file = traits.Any(mandatory=True,
                      desc="CSF mask.")
    vascular_mask_file = traits.Any(mandatory=True,
                      desc="vascular mask.")
    FD_file = File(exists=True, mandatory=True,
                      desc="CSV with framewise displacement.")
    motion_params_csv = File(exists=True, mandatory=True,
                      desc="CSV with motion parameters.")
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

        from nipype import logging
        nipype_log = logging.getLogger('nipype.workflow')

        with np.errstate(invalid='ignore', divide='ignore'):
            cleaning_out = clean_image(
                input_bold = self.inputs.bold_file,
                brain_mask = self.inputs.brain_mask_file,
                WM_mask = self.inputs.WM_mask_file,
                CSF_mask = self.inputs.CSF_mask_file,
                vascular_mask = self.inputs.vascular_mask_file,
                FD_csv = self.inputs.FD_file,
                motion_params_csv=self.inputs.motion_params_csv,
                timeseries_interval = cr_opts.timeseries_interval,
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
                slicewise_correction_direction=cr_opts.slicewise_correction_direction,
                nipype_log=nipype_log,
                )

        if cleaning_out is None:
            return runtime
        else:
            [timeseries_img, 
            CR_data_dict, 
            VE_spatial_map, 
            STD_spatial_map, 
            CR_STD_spatial_map, 
            random_CR_STD_spatial_map, 
            corrected_CR_STD_spatial_map] = cleaning_out

        '''
        Generate all the output files
        '''
        cr_out = os.getcwd()
        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(self.inputs.raw_input_file).name.rsplit(".nii")

        cleaned_path = cr_out+'/'+filename_split[0]+'_cleaned.nii.gz'
        sitk.WriteImage(timeseries_img, cleaned_path)
        VE_file_path = cr_out+'/'+filename_split[0]+'_VE_map.nii.gz'
        sitk.WriteImage(VE_spatial_map, VE_file_path)
        STD_file_path = cr_out+'/'+filename_split[0]+'_STD_map.nii.gz'
        sitk.WriteImage(STD_spatial_map, STD_file_path)
        CR_STD_file_path = cr_out+'/'+filename_split[0]+'_CR_STD_map.nii.gz'
        sitk.WriteImage(CR_STD_spatial_map, CR_STD_file_path)
        frame_mask_file = cr_out+'/'+filename_split[0]+'_frame_censoring_mask.csv'
        pd.DataFrame(CR_data_dict['frame_mask']).to_csv(frame_mask_file, index=False, header=['False = Masked Frames'])

        if cr_opts.generate_CR_null:
            random_CR_STD_file_path = cr_out+'/'+filename_split[0]+'_random_CR_STD_map.nii.gz'
            sitk.WriteImage(random_CR_STD_spatial_map, random_CR_STD_file_path)
            corrected_CR_STD_file_path = cr_out+'/'+filename_split[0]+'_corrected_CR_STD_map.nii.gz'
            sitk.WriteImage(corrected_CR_STD_spatial_map, corrected_CR_STD_file_path)
        else:
            random_CR_STD_file_path=''
            corrected_CR_STD_file_path=''

        CR_data_dict['cleaned_path'] = cleaned_path
        CR_data_dict['VE_file_path'] = VE_file_path
        CR_data_dict['STD_file_path'] = STD_file_path
        CR_data_dict['CR_STD_file_path'] = CR_STD_file_path
        CR_data_dict['random_CR_STD_file_path'] = random_CR_STD_file_path
        CR_data_dict['corrected_CR_STD_file_path'] = corrected_CR_STD_file_path
        CR_data_dict['frame_mask_file'] = frame_mask_file
                
        setattr(self, 'cleaned_path', cleaned_path)
        setattr(self, 'VE_file_path', VE_file_path)
        setattr(self, 'STD_file_path', STD_file_path)
        setattr(self, 'CR_STD_file_path', CR_STD_file_path)
        setattr(self, 'random_CR_STD_file_path', random_CR_STD_file_path)
        setattr(self, 'corrected_CR_STD_file_path', corrected_CR_STD_file_path)
        setattr(self, 'frame_mask_file', frame_mask_file)
        setattr(self, 'data_dict', CR_data_dict)
        if CR_data_dict['aroma_out'] is not None:
            setattr(self, 'aroma_out', CR_data_dict['aroma_out'])

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
    

def clean_image(input_bold, brain_mask, FD_csv, motion_params_csv, # necessary input files
                WM_mask=None, CSF_mask=None, vascular_mask=None,
                timeseries_interval='0,end', FD_censoring=False, FD_threshold=0.05, DVARS_censoring=False, minimum_timepoint=3, TR='auto', # replacing cr_opts
                detrending_order=1, detrending_time_interval='all', 
                apply_ica_aroma=False, ica_aroma_dim=0, ica_aroma_random_seed=1,
                match_number_timepoints=False, highpass=None, lowpass=None, edge_cutoff=0,
                nuisance_regressors = [], generate_CR_null=False,
                scale_variance_voxelwise=False,image_scaling='grand_mean_scaling',
                smoothing_filter=None,
                slicewise_correction_direction = 'Off',
                nipype_log=None,
                ):
    '''
    input_bold: timeseries file or SITK image to clean

    slicewise_correction_direction: if applying slicewise correction, detrending, bandpass filtering, nuisance regression
    and smoothing are all applied on a per slice basis. In such case, the output array confound_array will be the mean of 
    regressors computed across slices, to output a representative average. Similarly VE_temporal and VE_total_ratio are
    averaged across slices. If there is not sufficient data point within the brain masks for nuisance regressors, those 
    regressors are excluded for that slice - meaning that certain corrections are not necessarily applied in each slice.
    '''
    import os
    import numpy as np
    import pandas as pd
    import SimpleITK as sitk
    from rabies.utils import recover_3D,recover_4D
    from rabies.analysis_pkg.analysis_functions import closed_form
    from . import utils as cr_utils

    '''
    The function can take as input either an adequately pre-loaded python object or a file path
    '''
    def read_input(input):
        if (input is None) or isinstance(input, sitk.Image) or isinstance(input, pd.DataFrame):
            return input
        elif os.path.isfile(input):
            if '.nii' in input:
                return sitk.ReadImage(input, sitk.sitkFloat32)
            elif '.csv' in input:
                return pd.read_csv(input)
            else:
                raise ValueError(f"The input {input} was neither a .nii or .csv file.")
        else:
            raise ValueError(f"The input {input} was neither a file nor a recognized python object.")

    [input_bold, 
     FD_csv, 
     motion_params_df, 
     brain_mask, 
     WM_mask, 
     CSF_mask, 
     vascular_mask] = [read_input(input) for input in [
         input_bold, 
         FD_csv, 
         motion_params_csv, 
         brain_mask, 
         WM_mask, 
         CSF_mask, 
         vascular_mask,
         ]]

    volume_idx = sitk.GetArrayFromImage(brain_mask).astype(bool)
     # save the header and original array size, since the image will be removed to space memory below
    from rabies.utils import get_geometry_header
    header_geo = get_geometry_header(input_bold)
    orig_4d_size = input_bold.GetSize() # copy the original array size

    if TR=='auto':
        TR = float(header_geo.GetSpacing()[3])
    else:
        TR = float(TR)

    # save specifically the 6 rigid parameters for AROMA
    motion6_regressors_array = cr_utils.select_motion_regressors(['mot_6'],motion_params_df)
    if len(nuisance_regressors)==0:
        motion_regressors_array = motion6_regressors_array
    else:
        motion_regressors_array = cr_utils.select_motion_regressors(nuisance_regressors,motion_params_df)

    FD_trace = FD_csv.get('MeanFD')

    time_range = cr_utils.prep_timeseries_interval(timeseries_interval, num_frames=orig_4d_size[3])

    timeseries = sitk.GetArrayFromImage(input_bold)[:,volume_idx] # read directly as a 2D array
    del input_bold # free memory
    motion_regressors_array = motion_regressors_array[time_range, :]
    motion6_regressors_array = motion6_regressors_array[time_range, :]
    FD_trace = FD_trace[time_range]
    timeseries = timeseries[time_range,:]

    '''
    #1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)
    '''

    # compute the DVARS before denoising
    DVARS_trace = cr_utils.get_DVARS(timeseries)

    frame_mask = cr_utils.temporal_censoring(FD_trace, 
            FD_censoring, FD_threshold, DVARS_trace, DVARS_censoring, minimum_timepoint)
    if frame_mask is None:
        return None

    '''
    #2 - If --match_number_timepoints is selected, each scan is matched to the defined minimum_timepoint number of frames.
    '''
    if match_number_timepoints:
        if (highpass is not None) or (lowpass is not None):
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
    motion_regressors_array = motion_regressors_array[frame_mask]

    '''
    #3 - Linear/Quadratic detrending of fMRI timeseries and nuisance regressors
    '''
    timeseries, voxelwise_intercept = cr_utils.remove_trend(timeseries, frame_mask, order = detrending_order, time_interval = detrending_time_interval)
    motion_regressors_array, _ = cr_utils.remove_trend(motion_regressors_array, frame_mask, order = detrending_order, time_interval = detrending_time_interval)

    '''
    #4 - Apply ICA-AROMA.
    '''
    if apply_ica_aroma:
        # write intermediary output files for timeseries and 6 rigid body parameters
        inFile = os.path.abspath('aroma_input_timeseries.nii.gz')
        sitk.WriteImage(recover_4D(brain_mask, timeseries, header_geo), inFile)
        brain_mask_file = os.path.abspath("aroma_brain_mask.nii.gz")
        sitk.WriteImage(brain_mask, brain_mask_file)
        CSF_mask_file = os.path.abspath("aroma_CSF_mask.nii.gz")
        sitk.WriteImage(CSF_mask, CSF_mask_file)

        motion6_regressors_array=motion6_regressors_array[frame_mask,:]
        motion6_regressors_array, _ = cr_utils.remove_trend(motion6_regressors_array, frame_mask, order = detrending_order, time_interval = detrending_time_interval) # apply detrending to the confounds too
        df = pd.DataFrame(motion6_regressors_array)
        df.columns = ['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3']
        mc_file = os.path.abspath('aroma_motion_params.csv')
        df.to_csv(mc_file)
        aroma_out = os.path.abspath('aroma_out/')
        cleaned_file, aroma_out = cr_utils.exec_ICA_AROMA(inFile, mc_file, brain_mask_file, CSF_mask_file, aroma_out, TR, ica_aroma_dim, random_seed=ica_aroma_random_seed)
        timeseries = sitk.GetArrayFromImage(sitk.ReadImage(cleaned_file, sitk.sitkFloat32))[:,volume_idx] # read directly as a 2D array
    else:
        aroma_out = None

    if (highpass is not None) or (lowpass is not None):
        '''
        #5 - If frequency filtering and frame censoring are applied, simulate data in censored timepoints using the Lomb-Scargle periodogram, 
            as suggested in Power et al. (2014, Neuroimage), for both the fMRI timeseries and nuisance regressors prior to filtering.
        '''
        timeseries = cr_utils.lombscargle_fill(x=timeseries,time_step=TR,time_mask=frame_mask)
        motion_regressors_array = cr_utils.lombscargle_fill(x=motion_regressors_array,time_step=TR,time_mask=frame_mask)
        ### arrays are now interpolated

        '''
        #6 - As recommended in Lindquist et al. (2019, Human brain mapping), make the nuisance regressors orthogonal
            to the temporal filter.
        '''
        motion_regressors_array = cr_utils.butterworth(motion_regressors_array, TR=TR,
                                high_pass=highpass, low_pass=lowpass)

        '''
        #7 - Apply highpass and/or lowpass filtering on the fMRI timeseries (with simulated timepoints).
        '''

        timeseries = cr_utils.butterworth(timeseries, TR=TR,
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
        motion_regressors_array = motion_regressors_array[frame_mask]
    
    if frame_mask.sum()<int(minimum_timepoint):
        if nipype_log:
            nipype_log.warning(f"CONFOUND CORRECTION LEFT LESS THAN {str(minimum_timepoint)} VOLUMES. THIS SCAN WILL BE REMOVED FROM FURTHER PROCESSING.")
        return None

    '''
    #9 - If selected, compute the WM/CSF/vascular signal or aCompCorr and add to list of regressors. This is computed post-AROMA and filtering to 
        minimize re-introduction of previously corrected signal fluctuations.
    #10 - Apply confound regression using the selected nuisance regressors.
    '''
    regress_out = cr_utils.nuisance_regression(
        timeseries, motion_regressors_array, TR, frame_mask, orig_4d_size, brain_mask, WM_mask, CSF_mask, vascular_mask, nuisance_regressors=nuisance_regressors, 
        slicewise_correction_direction=slicewise_correction_direction, generate_CR_null=generate_CR_null, nipype_log=nipype_log)
    if regress_out is None:
        return None
    timeseries, predicted, predicted_random, num_regressors, VE_temporal, VE_spatial, VE_total_ratio, cleaned_slice_l = regress_out

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
            del scaled_timeseries
            [predicted,predicted_random] = scaled_list


    if image_scaling=='global_variance':
        scaling_factor = timeseries.std()
    elif image_scaling=='grand_mean_scaling':
        scaling_factor = voxelwise_intercept.mean()/100 # we scale BOLD in % fluctuations, hence dividing by 100
    elif image_scaling=='voxelwise_standardization':
        # each voxel is scaled according to its STD
        scaling_factor = timeseries.std(axis=0) 
    elif image_scaling=='voxelwise_mean':
        scaling_factor = voxelwise_intercept/100 # we scale BOLD in % fluctuations, hence dividing by 100
    else:
        scaling_factor = 1

    # scale fMRI timeseries, as well as the CR prediction to keep consistent scaling
    scaled_list = []
    for scaled_timeseries in [timeseries,predicted,predicted_random]:
        scaled_timeseries /= scaling_factor
        nan_voxels = np.isnan(scaled_timeseries).sum(axis=0)>1
        scaled_timeseries[:,nan_voxels] = 0
        scaled_list.append(scaled_timeseries)
    del scaled_timeseries
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
    tDOF = num_timepoints - (aroma_rm+num_regressors) + number_extra_timepoints

    # smoothing takes an SITK image
    timeseries_img = recover_4D(brain_mask, timeseries, header_geo)
    del timeseries

    if smoothing_filter is not None:
        '''
        #12 - Apply Gaussian spatial smoothing.
        '''
        if not slicewise_correction_direction=='Off':
            # smooth 1 slice at a time
            for slice in cleaned_slice_l:
                # here the axis direction is in proper order since we are indexing a SITK image
                if slicewise_correction_direction=='RL':
                    timeseries_img[slice:slice+1,:,:,:] = cr_utils.smooth_image(
                        timeseries_img[slice:slice+1,:,:,:], smoothing_filter, brain_mask[slice:slice+1,:,:])
                elif slicewise_correction_direction=='AP':
                    timeseries_img[:,slice:slice+1,:,:] = cr_utils.smooth_image(
                        timeseries_img[:,slice:slice+1,:,:], smoothing_filter, brain_mask[:,slice:slice+1,:])
                elif slicewise_correction_direction=='SI':
                    timeseries_img[:,:,slice:slice+1,:] = cr_utils.smooth_image(
                        timeseries_img[:,:,slice:slice+1,:], smoothing_filter, brain_mask[:,:,slice:slice+1])
                else:
                    raise ValueError(f"slicewise_correction_direction must be one of 'RL', 'AP' or 'SI', got {slicewise_correction_direction} instead.")
        else:
            timeseries_img = cr_utils.smooth_image(
                timeseries_img, smoothing_filter, brain_mask)

    # save output files
    VE_spatial_map = recover_3D(brain_mask, VE_spatial)
    STD_spatial_map = recover_3D(brain_mask, temporal_std)
    CR_STD_spatial_map = recover_3D(brain_mask, predicted_std)
    del VE_spatial, temporal_std, predicted_std
    if generate_CR_null:
        random_CR_STD_spatial_map = recover_3D(brain_mask, predicted_random_std)
        corrected_CR_STD_spatial_map = recover_3D(brain_mask, corrected_predicted_std)
        del predicted_random_std, corrected_predicted_std
    else:
        random_CR_STD_spatial_map = None
        corrected_CR_STD_spatial_map = None

    CR_data_dict = {
        'TR':TR, 'FD_trace':FD_trace, 'DVARS':DVARS_trace, 'time_range':time_range, 'frame_mask':frame_mask, 'VE_temporal':VE_temporal, 
        'motion_params_df':motion_params_df, 'predicted_time':predicted_time, 'tDOF':tDOF, 'CR_global_std':predicted_global_std, 
        'VE_total_ratio':VE_total_ratio, 'voxelwise_mean':voxelwise_intercept, 'aroma_out':aroma_out,
        }
    return timeseries_img, CR_data_dict, VE_spatial_map, STD_spatial_map, CR_STD_spatial_map, random_CR_STD_spatial_map, corrected_CR_STD_spatial_map
