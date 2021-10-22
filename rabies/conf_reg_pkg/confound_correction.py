#!/usr/bin/env python3
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function
from .utils import prep_CR
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

def init_confound_correction_wf(cr_opts, name="confound_correction_wf"):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=[
                        'bold_file', 'brain_mask', 'csf_mask', 'confounds_file', 'FD_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=[
                         'cleaned_path', 'aroma_out', 'VE_file', 'frame_mask_file', 'CR_data_dict']), name='outputnode')

    regress_node = pe.Node(Regress(cr_opts=cr_opts),
                           name='regress', mem_gb=1)

    prep_CR_node = pe.Node(Function(input_names=['bold_file', 'confounds_file', 'FD_file', 'cr_opts'],
                                              output_names=['data_dict'],
                                              function=prep_CR),
                                     name='prep_CR', mem_gb=1)
    prep_CR_node.inputs.cr_opts = cr_opts

    workflow.connect([
        (inputnode, prep_CR_node, [
            ("bold_file", "bold_file"),
            ("confounds_file", "confounds_file"),
            ("FD_file", "FD_file"),
            ]),
        (inputnode, regress_node, [
            ("bold_file", "bold_file"),
            ("brain_mask", "brain_mask_file"),
            ("csf_mask", "CSF_mask_file"),
            ]),
        (prep_CR_node, regress_node, [
            ("data_dict", "data_dict"),
            ]),
        (regress_node, outputnode, [
            ("cleaned_path", "cleaned_path"),
            ("VE_file_path", "VE_file"),
            ("frame_mask_file", "frame_mask_file"),
            ("data_dict", "CR_data_dict"),
            ("aroma_out", "aroma_out"),
            ]),
        ])


    return workflow


class RegressInputSpec(BaseInterfaceInputSpec):
    bold_file = File(exists=True, mandatory=True,
                      desc="Timeseries to denoise.")
    data_dict = traits.Dict(
        exists=True, mandatory=True, desc="Dictionary with extra inputs.")
    brain_mask_file = File(exists=True, mandatory=True,
                      desc="Brain mask.")
    CSF_mask_file = File(exists=True, mandatory=True,
                      desc="CSF mask.")
    cr_opts = traits.Any(
        exists=True, mandatory=True, desc="Processing specs.")

class RegressOutputSpec(TraitedSpec):
    cleaned_path = File(exists=True, mandatory=True,
                      desc="Cleaned timeseries.")
    VE_file_path = File(exists=True, mandatory=True,
                      desc="Variance explained map from confound regression.")
    frame_mask_file = File(exists=True, mandatory=True,
                      desc="Frame mask from temporal censoring.")
    data_dict = traits.Any(
        desc="A dictionary with key outputs.")
    aroma_out = traits.Any(
        desc="Output directory from ICA-AROMA.")

class Regress(BaseInterface):
    '''
    Apply a flexible confound regression algorithm in line with recommendations from
    human litterature. 
    
    #1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)
    #2 - Detrend timeseries and confound regressors
    #3 - Apply ICA-AROMA.
    #4 - If filtering is applied, simulate censored timepoints as in Power et al. 2014
         for both the timeseries and confound regressors prior to filtering.
    #5 - As recommended in Lindquist et al. 2019, make the confound regressors orthogonal
         to the temporal filter.
    #6 - Apply bandpass filtering on the timeseries (with filled missing values), and 
         apply again the temporal mask onto output timeseries.
    #7 - Apply confound regression using the corrected regressors, while applying the 
         temporal masks to both the regressors and timeseries to remove simulated data
         points.
    #8 - Standardize timeseries
    #9 - Apply spatial smoothing.
    
    References:
        
        Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). 
        Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage, 59(3), 2142-2154.
        
        Power, J. D., Mitra, A., Laumann, T. O., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2014). 
        Methods to detect, characterize, and remove motion artifact in resting state fMRI. Neuroimage, 84, 320-341.
        
        Lindquist, M. A., Geuter, S., Wager, T. D., & Caffo, B. S. (2019). 
        Modular preprocessing pipelines can reintroduce artifacts into fMRI data. Human brain mapping, 40(8), 2358-2376.
    '''

    input_spec = RegressInputSpec
    output_spec = RegressOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import pandas as pd
        import SimpleITK as sitk
        import nilearn.image
        import nibabel as nb
        from nilearn.signal import butterworth
        from scipy.signal import detrend
        from rabies.conf_reg_pkg.utils import recover_3D,recover_4D,temporal_censoring,lombscargle_fill, exec_ICA_AROMA
        from rabies.analysis_pkg.analysis_functions import closed_form

        bold_file = self.inputs.bold_file
        brain_mask_file = self.inputs.brain_mask_file
        CSF_mask_file = self.inputs.CSF_mask_file
        data_dict = self.inputs.data_dict
        cr_opts = self.inputs.cr_opts

        FD_trace=data_dict['FD_trace']
        confounds_array=data_dict['confounds_array']
        confounds_file=data_dict['confounds_csv']
        time_range=data_dict['time_range']
        confounds_6rigid_array=data_dict['confounds_6rigid_array']

        cr_out = os.getcwd()
        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

        brain_mask = sitk.GetArrayFromImage(sitk.ReadImage(brain_mask_file, sitk.sitkFloat32))
        volume_indices = brain_mask.astype(bool)

        data_img = sitk.ReadImage(bold_file, sitk.sitkFloat32)
        data_array = sitk.GetArrayFromImage(data_img)
        num_volumes = data_array.shape[0]
        timeseries = np.zeros([num_volumes, volume_indices.sum()])
        for i in range(num_volumes):
            timeseries[i, :] = (data_array[i, :, :, :])[volume_indices]
        timeseries = timeseries[time_range,:]

        if cr_opts.TR=='auto':
            TR = float(data_img.GetSpacing()[3])
        else:
            TR = float(cr_opts.TR)

        '''
        #1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)
        '''
        data_dict = temporal_censoring(timeseries, data_dict,
                cr_opts.FD_censoring, cr_opts.FD_threshold, cr_opts.DVARS_censoring, cr_opts.minimum_timepoint)

        if data_dict is None:
            empty_img = sitk.GetImageFromArray(np.empty([1,1]))
            empty_file = os.path.abspath('empty.nii.gz')
            sitk.WriteImage(empty_img, empty_file)

            setattr(self, 'cleaned_path', empty_file)
            setattr(self, 'VE_file_path', empty_file)
            setattr(self, 'frame_mask_file', empty_file)
            setattr(self, 'data_dict', None)
            setattr(self, 'aroma_out', None)

            return runtime
        
        timeseries=data_dict['timeseries']
        FD_trace=data_dict['FD_trace']
        DVARS=data_dict['DVARS']
        frame_mask=data_dict['frame_mask']
        confounds_array=data_dict['confounds_array']

        '''
        #2 - Detrend timeseries and confound regressors
        '''
        # apply simple detrending, after censoring
        timeseries = detrend(timeseries,axis=0)
        confounds_array = detrend(confounds_array,axis=0) # apply detrending to the confounds too


        '''
        #3 - Apply ICA-AROMA.
        '''
        if cr_opts.run_aroma:
            # write intermediary output files for timeseries and 6 rigid body parameters
            timeseries_3d = recover_4D(brain_mask_file, timeseries, bold_file)
            inFile = f'{cr_out}/{filename_split[0]}_aroma_input.nii.gz'
            sitk.WriteImage(timeseries_3d, inFile)

            confounds_6rigid_array=confounds_6rigid_array[frame_mask,:]
            confounds_6rigid_array = detrend(confounds_6rigid_array,axis=0) # apply detrending to the confounds too
            df = pd.DataFrame(confounds_6rigid_array)
            df.columns = ['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3']
            mc_file = f'{cr_out}/{filename_split[0]}_aroma_input.csv'
            df.to_csv(mc_file)

            cleaned_file, aroma_out = exec_ICA_AROMA(inFile, mc_file, brain_mask_file, CSF_mask_file, TR, cr_opts.aroma_dim)
            setattr(self, 'aroma_out', aroma_out)

            data_img = sitk.ReadImage(cleaned_file, sitk.sitkFloat32)
            data_array = sitk.GetArrayFromImage(data_img)
            num_volumes = data_array.shape[0]
            timeseries = np.zeros([num_volumes, volume_indices.sum()])
            for i in range(num_volumes):
                timeseries[i, :] = (data_array[i, :, :, :])[volume_indices]
        else:
            setattr(self, 'aroma_out', None)

        if (not cr_opts.highpass is None) or (not cr_opts.lowpass is None):
            '''
            #4 - If filtering is applied, simulate censored timepoints as in Power et al. 2014
                for both the timeseries and confound regressors prior to filtering.
            '''
            timeseries_filled = lombscargle_fill(x=timeseries,time_step=TR,time_mask=frame_mask)
            confounds_filled = lombscargle_fill(x=confounds_array,time_step=TR,time_mask=frame_mask)

            '''
            #5 - As recommended in Lindquist et al. 2019, make the confound regressors orthogonal
                to the temporal filter.
            '''
            confounds_filtered = butterworth(confounds_filled, sampling_rate=1. / TR,
                                    low_pass=cr_opts.lowpass, high_pass=cr_opts.highpass)

            '''
            #6 - Apply bandpass filtering on the timeseries (with filled missing values), and 
                apply again the temporal mask onto output timeseries.
            '''

            timeseries_filtered = butterworth(timeseries_filled, sampling_rate=1. / TR,
                                    low_pass=cr_opts.lowpass, high_pass=cr_opts.highpass)

            # re-apply the masks to take out simulated data points    
            timeseries = timeseries_filtered[frame_mask]
            confounds_array = confounds_filtered[frame_mask]
        
        '''
        #7 - Apply confound regression using the corrected regressors, while applying the 
            temporal masks to both the regressors and timeseries to remove simulated data
            points.
        '''

        # estimate the VE from the CR selection, or 6 rigid motion parameters if no CR is applied
        X=confounds_array
        Y=timeseries
        try:
            res = Y-X.dot(closed_form(X,Y))
        except:
            import logging
            log = logging.getLogger('root')
            log.debug("SINGULAR MATRIX ERROR DURING CONFOUND REGRESSION. THIS SCAN WILL BE REMOVED FROM FURTHER PROCESSING.")
            empty_img = sitk.GetImageFromArray(np.empty([1,1]))
            empty_file = os.path.abspath('empty.nii.gz')
            sitk.WriteImage(empty_img, empty_file)

            setattr(self, 'cleaned_path', empty_file)
            setattr(self, 'VE_file_path', empty_file)
            setattr(self, 'frame_mask_file', empty_file)
            setattr(self, 'data_dict', None)

            return runtime


        VE_spatial = 1-(res.var(axis=0)/Y.var(axis=0))
        VE_temporal = 1-(res.var(axis=1)/Y.var(axis=1))

        if len(cr_opts.conf_list) > 0:
            # if confound regression is applied
            timeseries = res

        '''
        #8 - Standardize timeseries
        '''
        if cr_opts.standardize:
            timeseries = (timeseries-timeseries.mean(axis=0))/timeseries.std(axis=0)

        # save output files
        VE_spatial_map = recover_3D(brain_mask_file, VE_spatial)
        timeseries_3d = recover_4D(brain_mask_file, timeseries, bold_file)
        cleaned_path = cr_out+'/'+filename_split[0]+'_cleaned.nii.gz'
        sitk.WriteImage(timeseries_3d, cleaned_path)
        VE_file_path = cr_out+'/'+filename_split[0]+'_VE_map.nii.gz'
        sitk.WriteImage(VE_spatial_map, VE_file_path)
        frame_mask_file = cr_out+'/'+filename_split[0]+'_frame_censoring_mask.csv'
        pd.DataFrame(frame_mask).to_csv(frame_mask_file, index=False, header=['False = Masked Frames'])

        if cr_opts.smoothing_filter is not None:
            '''
            #9 - Apply spatial smoothing.
            '''
            timeseries_3d = nilearn.image.smooth_img(nb.load(cleaned_path), cr_opts.smoothing_filter)
            timeseries_3d.to_filename(cleaned_path)

        data_dict = {'FD_trace':FD_trace, 'DVARS':DVARS, 'time_range':time_range, 'frame_mask':frame_mask, 'confounds_array':confounds_array, 'VE_temporal':VE_temporal, 'confounds_csv':confounds_file}

        setattr(self, 'cleaned_path', cleaned_path)
        setattr(self, 'VE_file_path', VE_file_path)
        setattr(self, 'frame_mask_file', frame_mask_file)
        setattr(self, 'data_dict', data_dict)

        return runtime

    def _list_outputs(self):
        return {'cleaned_path': getattr(self, 'cleaned_path'),
                'VE_file_path': getattr(self, 'VE_file_path'),
                'frame_mask_file': getattr(self, 'frame_mask_file'),
                'data_dict': getattr(self, 'data_dict'),
                'aroma_out': getattr(self, 'aroma_out'),
                }
