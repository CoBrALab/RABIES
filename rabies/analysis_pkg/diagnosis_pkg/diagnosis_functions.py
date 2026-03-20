import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rabies.utils import recover_3D
from rabies.analysis_pkg import analysis_functions
import SimpleITK as sitk
import nilearn.plotting
from .analysis_QC import masked_plot, threshold_top_percent

def compute_spatiotemporal_features(CR_data_dict, sub_maps_data_dict, common_maps_data_dict, analysis_dict, 
                                    prior_bold_idx, prior_confound_idx,
                                    nativespace_analysis=False,resampling_specs=None,
                                    ):
    if resampling_specs is None:
        resampling_specs={}
    # sub_maps_data_dict is the maps_data_dict that overlaps with the subject, either in native or common space depending on confound correction
    # common_maps_data_dict contains files in commonspace
    temporal_info = {}
    spatial_info = {}
    temporal_info['name_source'] = CR_data_dict['name_source']
    spatial_info['name_source'] = CR_data_dict['name_source']
    # spatial maps will always be in commonspace for display
    spatial_info['mask_file'] = common_maps_data_dict['mask_file']

    # these timeseries can be in nativespace or commonspace depending on confound correction stage
    timeseries = CR_data_dict['timeseries']
    volume_indices = common_maps_data_dict['volume_indices']

    CR_CR_data_dict = CR_data_dict['CR_data_dict']

    temporal_info['FD_trace'] = CR_CR_data_dict['FD_trace']
    temporal_info['DVARS_pre_correction'] = CR_CR_data_dict['DVARS_pre_correction']
    temporal_info['mse_trace'] = CR_CR_data_dict['mse_trace']
    temporal_info['VE_temporal'] = CR_CR_data_dict['VE_temporal']


    ### SBC analysis
    if len(analysis_dict['seed_map_files'])>0:
        seed_list=[]
        for seed_map in analysis_dict['seed_map_files']:
            seed_list.append(np.asarray(
                sitk.GetArrayFromImage(sitk.ReadImage(seed_map)))[volume_indices])
        spatial_info['seed_map_list'] = seed_list
        time_list=[]
        for time_csv in analysis_dict['seed_timecourse_csv']:
            time_list.append(np.array(pd.read_csv(time_csv, header=None)).flatten())
        temporal_info['SBC_time'] = np.array(time_list).T # convert to time by network matrix
    else:
        spatial_info['seed_map_list'] = None
        temporal_info['SBC_time'] = None

    ### DR analysis
    DR_W = np.array(pd.read_csv(analysis_dict['dual_regression_timecourse_csv'], header=None))
    DR_array = sitk.GetArrayFromImage(
        sitk.ReadImage(analysis_dict['dual_regression_nii']))
    if len(DR_array.shape)==3: # if there was only one component, need to convert to 4D array
        DR_array = DR_array[np.newaxis,:,:,:]
    DR_C = np.zeros([DR_array.shape[0], volume_indices.sum()])
    for i in range(DR_array.shape[0]):
        DR_C[i, :] = (DR_array[i, :, :, :])[volume_indices]

    temporal_info['DR_all'] = DR_W
    temporal_info['DR_bold'] = DR_W[:, prior_bold_idx]
    temporal_info['DR_confound'] = DR_W[:, prior_confound_idx]

    '''Temporal Features'''
    global_signal = timeseries.mean(axis=1)
    temporal_info['global_signal'] = global_signal
    temporal_info['predicted_time'] = CR_CR_data_dict['predicted_time']
    # take regional timecourse by taking the RMS each 3D volume
    edge_idx = sub_maps_data_dict['edge_idx']
    #edge_trace = np.sqrt((timeseries.T[edge_idx]**2).mean(axis=0))
    edge_trace = timeseries[:,edge_idx].mean(axis=1)
    temporal_info['edge_trace'] = edge_trace
    if sub_maps_data_dict['WM_idx'] is not None:
        WM_idx = sub_maps_data_dict['WM_idx']
        #WM_trace = np.sqrt((timeseries.T[WM_idx]**2).mean(axis=0))
        WM_trace = timeseries[:,WM_idx].mean(axis=1)
        temporal_info['WM_trace'] = WM_trace
    else:
        temporal_info['WM_trace'] = None
    if sub_maps_data_dict['CSF_idx'] is not None:
        CSF_idx = sub_maps_data_dict['CSF_idx']
        #CSF_trace = np.sqrt((timeseries.T[CSF_idx]**2).mean(axis=0))
        CSF_trace = timeseries[:,CSF_idx].mean(axis=1)
        temporal_info['CSF_trace'] = CSF_trace
    else:
        temporal_info['CSF_trace'] = None

    '''Spatial Features'''
    GS_cov = (global_signal.reshape(-1,1)*timeseries).mean(axis=0) # calculate the covariance between global signal and each voxel
    GS_corr = analysis_functions.vcorrcoef(timeseries.T, global_signal)

    prior_fit_out = {'C': None, 'W': None}
    if not analysis_dict['NPR_prior_filename'] is None:
        prior_fit_out['W'] = np.array(pd.read_csv(analysis_dict['NPR_prior_timecourse_csv'], header=None))
        C_array = sitk.GetArrayFromImage(
            sitk.ReadImage(analysis_dict['NPR_prior_filename']))
        if len(C_array.shape)==3: # if there was only one component, need to convert to 4D array
            C_array = C_array[np.newaxis,:,:,:]

        C = np.zeros([C_array.shape[0], volume_indices.sum()])
        for i in range(C_array.shape[0]):
            C[i, :] = (C_array[i, :, :, :])[volume_indices]
        prior_fit_out['C'] = C

    if nativespace_analysis: # certain spatial maps were computed in nativespace and need resampling
        from rabies.utils import resample_volumes
        import tempfile
        tmppath = tempfile.mkdtemp()
        def resample_brain_map(brain_map):
            image_3d = recover_3D(sub_maps_data_dict['mask_file'], brain_map)
            resampled_img = resample_volumes(image_3d, common_maps_data_dict['anat_ref_file'], 
                                             transforms_3d_files=resampling_specs['transforms'], inverses_3d=resampling_specs['inverses'], 
                                             motcorr_params_file = None, interpolation=resampling_specs['interpolation'], 
                                             rabies_data_type=resampling_specs['rabies_data_type'], clip_negative=False)
            resampled_brain_map = sitk.GetArrayFromImage(resampled_img)[common_maps_data_dict['volume_indices']]
            return resampled_brain_map
            
        [GS_cov, GS_corr, VE_spatial, temporal_std, predicted_std] = [
            resample_brain_map(brain_map) for brain_map in [GS_cov, GS_corr, CR_data_dict['VE_spatial'], 
                                                            CR_data_dict['temporal_std'],CR_data_dict['predicted_std']]
            ]
        import shutil
        shutil.rmtree(tmppath, ignore_errors=True)
    else:
        VE_spatial, temporal_std, predicted_std = [CR_data_dict['VE_spatial'], CR_data_dict['temporal_std'],CR_data_dict['predicted_std']]


    spatial_info['prior_maps'] = common_maps_data_dict['prior_map_vectors'][prior_bold_idx]
    spatial_info['DR_BOLD'] = DR_C[prior_bold_idx]
    spatial_info['DR_all'] = DR_C

    spatial_info['NPR_maps'] = prior_fit_out['C']
    temporal_info['NPR_time'] = prior_fit_out['W']

    spatial_info['VE_spatial'] = VE_spatial
    spatial_info['temporal_std'] = temporal_std
    spatial_info['predicted_std'] = predicted_std
    spatial_info['GS_corr'] = GS_corr
    spatial_info['GS_cov'] = GS_cov

    return temporal_info, spatial_info


def temporal_external_formating(temporal_info):
    import os
    import pandas as pd
    import pathlib  # Better path manipulation

    filename_split = pathlib.Path(
        temporal_info['name_source']).name.rsplit(".nii")

    del temporal_info['DR_all'], temporal_info['DR_bold'],temporal_info['DR_confound'],temporal_info['NPR_time'],temporal_info['SBC_time']

    temporal_info_csv = os.path.abspath(filename_split[0]+'_temporal_info.csv')
    pd.DataFrame(temporal_info).to_csv(temporal_info_csv)
    return temporal_info_csv


def spatial_external_formating(spatial_info):
    import os
    import pathlib  # Better path manipulation
    import SimpleITK as sitk
    from rabies.utils import recover_3D

    mask_file = spatial_info['mask_file']
    filename_split = pathlib.Path(
        spatial_info['name_source']).name.rsplit(".nii")

    VE_filename = os.path.abspath(filename_split[0]+'_VE.nii.gz')
    sitk.WriteImage(recover_3D(
        mask_file, spatial_info['VE_spatial']), VE_filename)

    std_filename = os.path.abspath(filename_split[0]+'_tSTD.nii.gz')
    sitk.WriteImage(recover_3D(
        mask_file, spatial_info['temporal_std']), std_filename)

    predicted_std_filename = os.path.abspath(filename_split[0]+'_predicted_std.nii.gz')
    sitk.WriteImage(recover_3D(
        mask_file, spatial_info['predicted_std']), predicted_std_filename)

    GS_corr_filename = os.path.abspath(filename_split[0]+'_GS_corr.nii.gz')
    sitk.WriteImage(recover_3D(
        mask_file, spatial_info['GS_corr']), GS_corr_filename)

    GS_cov_filename = os.path.abspath(filename_split[0]+'_GS_cov.nii.gz')
    sitk.WriteImage(recover_3D(
        mask_file, spatial_info['GS_cov']), GS_cov_filename)

    return VE_filename, std_filename, predicted_std_filename, GS_corr_filename, GS_cov_filename


'''
Subject-level QC
'''


def carpetplot(timeseries, ax, frame_mask=None):
    vmax = timeseries.std()
    if frame_mask is not None:
        filled_timeseries = np.zeros([len(frame_mask), timeseries.shape[1]])
        filled_timeseries[frame_mask] = timeseries
        filled_timeseries[(frame_mask==False)] = None
        censoring_array = np.ones([len(frame_mask), timeseries.shape[1]])
        censoring_array[frame_mask] = None
        timeseries = filled_timeseries
    im = ax.imshow(timeseries.T, cmap='gray',
                   vmax=vmax, vmin=-vmax, aspect='auto')
    if frame_mask is not None:
        im = ax.imshow(censoring_array.T, cmap='Reds',
                       vmax=1, vmin=0, alpha=0.5, aspect='auto')
    return im


def plot_freqs(ax,timeseries, TR, frame_mask):
    if (frame_mask==0).sum()>0: # only simulate data points if some censoring was applied
        from rabies.confound_correction_pkg.utils import lombscargle_fill
        timeseries_ = lombscargle_fill(x=timeseries,time_step=TR,time_mask=frame_mask)
    else:
        timeseries_ = timeseries

    # voxels that have a NaN value are set to 0
    timeseries_[np.isnan(timeseries_)] = 0

    from scipy.signal import welch
    fs = 1/TR
    N = timeseries.shape[0]
    # N//4 : proxy to scale the resolution of the filter to the number of data points
    # minimum of 130 to be able to detect frequency spikes in the spectrum
    nperseg = max(N // 4, 130)
    freqs, psds = welch(timeseries_, fs=fs, nperseg=nperseg, axis=0)

    ps_mean = psds.mean(axis=1)
    ps_std = psds.std(axis=1)

    y_max = ps_mean[freqs>0.01].max()*1.5

    ax.plot(freqs, ps_mean)
    ax.fill_between(freqs, (ps_mean-ps_std),(ps_mean+ps_std), alpha=0.4)
    ax.set_ylim([0,y_max])
    #xticks = np.arange(0,freqs.max(),0.05)
    #ax.set_xticks(xticks)    

    x_max = freqs.max()
    xlim_edge_percent = 0.02
    x_lim = [-(x_max*xlim_edge_percent),x_max+(x_max*xlim_edge_percent)]
    ax.set_xlim(x_lim)


def scan_diagnosis(CR_data_dict, maps_data_dict, temporal_info, spatial_info, plot_seed_frequencies={}, brainmap_percent_threshold=10, display_censoring=True):
    TR = CR_data_dict['CR_data_dict']['TR']

    total_time_length = len(CR_data_dict['CR_data_dict']['frame_mask'])
    n_frames_per_fig = 1000 # max 1000 frames per figure
    num_figs = int(np.ceil(total_time_length/n_frames_per_fig))
    start_time = 0
    temporal_fig_list = []
    for n in range(num_figs):
        plot_time_subset = range(start_time, min(total_time_length, start_time+n_frames_per_fig))
        start_time+=n_frames_per_fig

        [
            timeseries, frame_mask, motion_params_df,
            FD_trace, DVARS, mse_trace,
            CR_VE_temporal, CR_var_temporal,
            global_signal, WM_trace, CSF_trace, edge_trace,
            DR_bold_time, DR_confound_time, SBC_time, NPR_time,
            ] = prep_temporal_subset_input(plot_time_subset, CR_data_dict, temporal_info)

        time0 = plot_time_subset.start
        fig,axes = temporal_diagnosis_plot(
            time0, timeseries, frame_mask, TR, motion_params_df,
            FD_trace, DVARS, mse_trace,
            CR_VE_temporal, CR_var_temporal,
            global_signal, WM_trace, CSF_trace, edge_trace,
            DR_bold_time, DR_confound_time, SBC_time, NPR_time,
            plot_seed_frequencies=plot_seed_frequencies, display_censoring=display_censoring)
        temporal_fig_list.append(fig)
        plt.close(fig)

    fig2 = spatial_diagnosis_plot(maps_data_dict, spatial_info, brainmap_percent_threshold=brainmap_percent_threshold)
    plt.close(fig2)
    return temporal_fig_list, fig2


def prep_temporal_subset_input(plot_time_subset, CR_data_dict, temporal_info):
    CR_CR_data_dict = CR_data_dict['CR_data_dict']

    frame_mask = CR_CR_data_dict['frame_mask'][plot_time_subset]
    # the motion parameters are the entire original length, so both the time_range from confound_correction and plot_time_subset must be applied
    time_range = CR_CR_data_dict['time_range']
    motion_params_df = CR_CR_data_dict['motion_params_df'][time_range.start:time_range.stop][plot_time_subset.start:plot_time_subset.stop]

    # the remaining timeseries were censored with frame_mask, and plot_time_subset is not aligned with the censored time axis
    first = plot_time_subset.start
    last = plot_time_subset.stop # .stop actually returns the max value for range(), which is excluded as an index
    # the idx range must be shifted by the number of censored frames before the first index to match the censored timeseries array
    frame_mask_full = CR_CR_data_dict['frame_mask']
    if first>0:
        first -= (frame_mask_full[:first]==False).sum()
    # similarly for the last index
    last -= (frame_mask_full[:last]==False).sum()
    plot_time_subset_censored = range(first,last) # we re-generated a corrected range

    timeseries = CR_data_dict['timeseries'][plot_time_subset_censored]

    FD_trace=temporal_info['FD_trace'].to_numpy()[plot_time_subset_censored]
    DVARS=temporal_info['DVARS_pre_correction'][plot_time_subset_censored]
    mse_trace=temporal_info['mse_trace'][plot_time_subset_censored]

    CR_VE_temporal=temporal_info['VE_temporal'][plot_time_subset_censored]
    CR_var_temporal=temporal_info['predicted_time'][plot_time_subset_censored]

    SBC_time=temporal_info['SBC_time']
    if SBC_time is not None: # might be None
        SBC_time=SBC_time[plot_time_subset_censored]
    DR_bold_time=temporal_info['DR_bold'][plot_time_subset_censored]
    DR_confound_time=temporal_info['DR_confound'][plot_time_subset_censored]
    NPR_time=temporal_info['NPR_time']
    if NPR_time is not None: # might be None
        NPR_time=NPR_time[plot_time_subset_censored]

    global_signal=temporal_info['global_signal'][plot_time_subset_censored]
    edge_trace=temporal_info['edge_trace'][plot_time_subset_censored]
    WM_trace=temporal_info['WM_trace'][plot_time_subset_censored]
    CSF_trace=temporal_info['CSF_trace'][plot_time_subset_censored]
    return [
        timeseries, frame_mask, motion_params_df,
        FD_trace, DVARS, mse_trace,
        CR_VE_temporal, CR_var_temporal,
        global_signal, WM_trace, CSF_trace, edge_trace,
        DR_bold_time, DR_confound_time, SBC_time, NPR_time,
    ]


def temporal_diagnosis_plot(
        time0, timeseries, frame_mask, TR, motion_params_df,
        FD_trace, DVARS, mse_trace,
        CR_VE_temporal, CR_var_temporal,
        global_signal, WM_trace, CSF_trace, edge_trace,
        DR_bold_time, DR_confound_time, SBC_time, NPR_time,
        plot_seed_frequencies={}, display_censoring=True):

    if display_censoring:
        num_timepoints = len(frame_mask)
    else:
        num_timepoints = frame_mask.sum()
    
    fig_width = max(1, int(np.ceil(num_timepoints / 40)))
    rows_height = [2,0.1,4,1,1,1,2,2,2,2]
    fig, axes = plt.subplots(
        nrows=len(rows_height), ncols=1,
        gridspec_kw={'height_ratios': rows_height},
        figsize=(fig_width, 24),
    )

    import matplotlib.transforms as mtransforms
    # set a consistent # of pixels to offset the axis legend, regardless of figure size
    offset = mtransforms.ScaledTranslation(70/72, 0, fig.dpi_scale_trans)

    a=0
    ax = axes[a]
    plot_freqs(ax,timeseries, TR, frame_mask)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_xticklabels(), fontsize=15)
    ax.set_xlabel('Frequency (Hz)', fontsize=20)
    ax.set_ylabel('Power (a.u.)', fontsize=20)

    freq_legend = ['Whole brain']
    # plot the spectrum for specific seeds
    seed_names = list(plot_seed_frequencies.keys())
    for seed_name in seed_names:
        seed_idx = plot_seed_frequencies[seed_name]
        plot_freqs(ax,SBC_time[:,[seed_idx]], TR, frame_mask)
        freq_legend.append(seed_name)
    ax.legend(
        freq_legend,
        loc='center left', 
        bbox_transform=ax.transAxes+offset,    
        fontsize=15, 
        bbox_to_anchor=(1, 0.5))
    a+=1

    ax = axes[a]
    ax.axis('off')
    a+=1 # skip the empty axis
    
    # set a common x time axis
    x = range(num_timepoints)
    xlim_edge_percent = 0.02
    x_lim = [-(num_timepoints*xlim_edge_percent),(num_timepoints-1)+(num_timepoints*xlim_edge_percent)]

    x_ticks = np.array([i*100 for i in range(int(num_timepoints/100)+1)])
    x_ticks_labels = x_ticks+time0
    x_ticks_labels = np.ceil(x_ticks_labels.astype(float)*TR).astype(int) # convert to real time in seconds
    for ax in axes[1:]:
        ax.set_xlim(x_lim)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks_labels)


    # carpet plot
    ax = axes[a]
    if display_censoring:
        im = carpetplot(timeseries, ax, frame_mask)
    else:
        im = carpetplot(timeseries, ax, frame_mask=None)
    ax.set_ylabel('Voxels', fontsize=20)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.axes.get_yaxis().set_ticks([])
    plt.setp(ax.get_xticklabels(), fontsize=15)
    plt.setp(ax.get_yticklabels(), fontsize=15)
    a+=1

    # plot the motion timecourses
    ax = axes[a]
    if display_censoring:
        ax.plot(x,motion_params_df['mov1'].to_numpy(), alpha=0.7)
        ax.plot(x,motion_params_df['mov2'].to_numpy(), alpha=0.7)
        ax.plot(x,motion_params_df['mov3'].to_numpy(), alpha=0.7)
    else:
        ax.plot(x,motion_params_df['mov1'].to_numpy()[frame_mask], alpha=0.7)
        ax.plot(x,motion_params_df['mov2'].to_numpy()[frame_mask], alpha=0.7)
        ax.plot(x,motion_params_df['mov3'].to_numpy()[frame_mask], alpha=0.7)

    ax.legend(
        ['translation 1', 'translation 2', 'translation 3'],
        loc='center left', 
        fontsize=15, 
        bbox_transform=ax.transAxes+offset,    
        bbox_to_anchor=(1, 0.5))
    ax.set_ylabel('mm', fontsize=20)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), fontsize=15)
    a+=1

    ax = axes[a]
    if display_censoring:
        ax.plot(x,motion_params_df['rot1'].to_numpy(), alpha=0.7)
        ax.plot(x,motion_params_df['rot2'].to_numpy(), alpha=0.7)
        ax.plot(x,motion_params_df['rot3'].to_numpy(), alpha=0.7)
    else:
        ax.plot(x,motion_params_df['rot1'].to_numpy()[frame_mask], alpha=0.7)
        ax.plot(x,motion_params_df['rot2'].to_numpy()[frame_mask], alpha=0.7)
        ax.plot(x,motion_params_df['rot3'].to_numpy()[frame_mask], alpha=0.7)
    ax.legend(
        ['Euler angle 1', 'Euler angle 2', 'Euler angle 3'],
        loc='center left', 
        fontsize=15, 
        bbox_transform=ax.transAxes+offset,    
        bbox_to_anchor=(1, 0.5))
    ax.set_ylabel('radians', fontsize=20)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), fontsize=15)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    a+=1

    def fill_censored_vector(vector,frame_mask):
        if not display_censoring:
            return vector
        filled = np.zeros(len(frame_mask))
        filled[frame_mask] = vector
        filled[frame_mask==False] = None
        return filled
    
    ax = axes[a]
    y = fill_censored_vector(mse_trace, frame_mask)
    ax.plot(x,y, 'r', alpha=0.7)
    ax.set_ylabel('MSE', fontsize=20)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), fontsize=15)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(
        ['Framewise\ndistance\nfrom mean'], 
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        bbox_transform=ax.transAxes+offset,    
        fontsize=15,
        )
    a+=1
    
    ax = axes[a]
    y = fill_censored_vector(FD_trace, frame_mask)
    ax.plot(x,y, 'r', alpha=0.7)
    ax.set_ylabel('FD (mm)', fontsize=20)
    ax_ = ax.twinx()

    y2 = fill_censored_vector(DVARS, frame_mask)
    ax_.plot(x,y2, 'b', alpha=0.7)
    ax_.set_ylabel('DVARS', fontsize=20)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax_.spines['top'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), fontsize=15)
    plt.setp(ax_.get_xticklabels(), visible=False)
    plt.setp(ax_.get_yticklabels(), fontsize=15)

    ax.legend(
        ['Framewise \nDisplacement (FD)'], 
        loc='center left',
        bbox_to_anchor=(1, 0.7),
        bbox_transform=ax.transAxes+offset,    
        fontsize=15,
        )
    ax_.legend(
        ['DVARS\n(before denoising)'], 
        loc='center left',
        bbox_to_anchor=(1, 0.3),
        bbox_transform=ax.transAxes+offset,    
        fontsize=15,
        )
    a+=1

    ax = axes[a]
    ax.plot(x,fill_censored_vector(CR_var_temporal, frame_mask), 'red', alpha=0.7)
    ax.set_ylabel('$CR_{var}$', fontsize=20)
    ax_ = ax.twinx()
    ax_.plot(x,fill_censored_vector(CR_VE_temporal, frame_mask), 'darkviolet', alpha=0.7)
    ax_.set_ylabel('CR $\mathregular{R^2}$', fontsize=20)
    ax_.spines['right'].set_visible(False)
    ax_.spines['top'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), fontsize=15)
    plt.setp(ax_.get_xticklabels(), visible=False)
    plt.setp(ax_.get_yticklabels(), fontsize=15)
    ax.legend(
        ['$CR_{var}$'], 
        loc='center left', 
        fontsize=15, 
        bbox_transform=ax.transAxes+offset,    
        bbox_to_anchor=(1, 0.7))
    ax_.legend(
        ['CR $\mathregular{R^2}$'], 
        loc='center left', 
        fontsize=15, 
        bbox_transform=ax.transAxes+offset,    
        bbox_to_anchor=(1, 0.3))
    ax_.set_ylim([0,1])
    a+=1
    
    ax = axes[a]
    ax.plot(x,fill_censored_vector(global_signal, frame_mask), alpha=0.7)
    ax.plot(x,fill_censored_vector(edge_trace, frame_mask), alpha=0.7)
    if WM_trace is not None:
        ax.plot(x,fill_censored_vector(WM_trace, frame_mask), alpha=0.7)
    else:
        ax.plot([],[])
    if CSF_trace is not None:
        ax.plot(x,fill_censored_vector(CSF_trace, frame_mask), alpha=0.7)
    else:
        ax.plot([],[])
    ax.set_ylabel('Mean signal', fontsize=20)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), fontsize=15)
    ax.legend(
        ['Whole-brain mask\n(global signal)', 'Edge Mask', 'WM Mask', 'CSF Mask'], 
        loc='center left', 
        fontsize=15, 
        bbox_transform=ax.transAxes+offset,    
        bbox_to_anchor=(1, 0.5))
    a+=1

    ax = axes[a]
    # we take the mean of the timecourse amplitude (absolute values) to summarized across all components
    import matplotlib.cm as cm
    Greens_colors = cm.Greens(np.linspace(0.5, 0.8, 2))
    Blues_colors = cm.Blues(np.linspace(0.5, 0.8, 2))
    Purples_colors = cm.Purples(np.linspace(0.5, 0.8, 3))
    YlOrRd_colors = cm.YlOrRd(np.linspace(0.3, 0.8, 4))
    legend=[]
    for network_time,name,color,scaler in zip(
            [DR_confound_time, DR_bold_time,SBC_time,NPR_time],
            ['DR confounds', 'DR networks','SBC networks','NPR networks'],
            [YlOrRd_colors[2], Blues_colors[1],Purples_colors[1],Greens_colors[1]],
            [0,-1,-1,-1]):
        if network_time is not None:
            # make sure the timecourses are normalized
            network_time = network_time.copy()
            network_time /= np.sqrt((network_time ** 2).mean(axis=0))
            network_time_avg = np.abs(network_time).mean(axis=1)
            # we introduce a scaler to prevent overlap of confound with network timecourses
            network_time_avg += scaler
            ax.plot(x,fill_censored_vector(network_time_avg, frame_mask), color=color, alpha=0.7)
            legend.append(name)

    ax.legend(
        legend, 
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        bbox_transform=ax.transAxes+offset,    
        fontsize=15,
        )

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Time (seconds)', fontsize=25)
    ax.set_ylabel('Mean amplitude \n(a.u.)', fontsize=20)
    ax.yaxis.set_ticklabels([])
    plt.setp(ax.get_xticklabels(), fontsize=15)
    plt.setp(ax.get_yticklabels(), fontsize=15)
    a+=1
    return fig,axes


def spatial_diagnosis_plot(maps_data_dict, spatial_info, brainmap_percent_threshold=10):
    template_file = maps_data_dict['anat_ref_file']
    dr_maps = spatial_info['DR_BOLD']
    SBC_maps = spatial_info['seed_map_list']
    NPR_maps = spatial_info['NPR_maps']
    mask_file = maps_data_dict['mask_file']

    # convert to empty list to read len() of 0
    if SBC_maps is None:
        SBC_maps=[]
    if NPR_maps is None:
        NPR_maps=[]

    nrows = 4+dr_maps.shape[0]+len(SBC_maps)+len(NPR_maps)

    fig2, axes2 = plt.subplots(nrows=nrows, ncols=3, figsize=(12*3, 2*nrows))
    plt.tight_layout()

    from rabies.visualization import otsu_scaling, plot_3d

    axes = axes2[0, :]
    scaled = otsu_scaling(template_file)
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    temporal_std = spatial_info['temporal_std']
    sitk_img = recover_3D(
        mask_file, temporal_std)

    # select vmax at 95th percentile value
    vector = temporal_std.flatten()
    vector.sort()
    vmax = vector[int(len(vector)*0.95)]
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=0, vmax=vmax,
            cmap='inferno', alpha=1, cbar=True, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 35
        cbar.set_label('Standard \n Deviation', fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)
    for ax in axes:
        ax.set_title('$\mathregular{BOLD_{SD}}$', fontsize=30, color='white')


    axes = axes2[1, :]
    scaled = otsu_scaling(template_file)
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    predicted_std = spatial_info['predicted_std']
    sitk_img = recover_3D(
        mask_file, predicted_std)

    # select vmax at 95th percentile value
    vector = predicted_std.flatten()
    vector.sort()
    vmax = vector[int(len(vector)*0.95)]
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=0, vmax=vmax,
            cmap='inferno', alpha=1, cbar=True, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 35
        cbar.set_label('Standard \n Deviation', fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)
    for ax in axes:
        ax.set_title('$\mathregular{CR_{SD}}$', fontsize=30, color='white')


    axes = axes2[2, :]
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    sitk_img = recover_3D(
        mask_file, spatial_info['VE_spatial'])
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=0, vmax=1, cmap='inferno',
            alpha=1, cbar=True, threshold=0.1, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label('$\mathregular{R^2}$', fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)
    for ax in axes:
        ax.set_title('CR $\mathregular{R^2}$', fontsize=30, color='white')

    axes = axes2[3, :]
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    sitk_img = recover_3D(
        mask_file, spatial_info['GS_cov'])
    # select vmax at 95th percentile value
    vector = spatial_info['GS_cov'].flatten()
    vector.sort()
    vmax = vector[int(len(vector)*0.95)]
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=-vmax, vmax=vmax, cmap='cold_hot',
            alpha=1, cbar=True, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Covariance", fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)
    for ax in axes:
        ax.set_title('Global Signal Covariance', fontsize=30, color='white')

    for i in range(dr_maps.shape[0]):
        axes = axes2[i+4, :]

        map = dr_maps[i, :]
        threshold = threshold_top_percent(map, top_percent=brainmap_percent_threshold)
        mask=np.abs(map)>=threshold # taking absolute values to include negative weights
        mask_img = recover_3D(mask_file,mask)        
        sitk_img = recover_3D(
            mask_file, map)
        cbar_list = masked_plot(fig2,axes, sitk_img, scaled, mask_img=mask_img, vmax=None)

        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 35
            cbar.set_label("Beta \nCoefficient", fontsize=17, rotation=270, color='white')
            cbar.ax.tick_params(labelsize=15)
        for ax in axes:
            ax.set_title(f'DR network {i}', fontsize=30, color='white')

    for i in range(len(SBC_maps)):
        axes = axes2[i+4+dr_maps.shape[0], :]

        map = SBC_maps[i]
        threshold = threshold_top_percent(map, top_percent=brainmap_percent_threshold)
        mask=np.abs(map)>=threshold # taking absolute values to include negative weights
        mask_img = recover_3D(mask_file,mask)        
        sitk_img = recover_3D(
            mask_file, map)
        cbar_list = masked_plot(fig2,axes, sitk_img, scaled, mask_img=mask_img, vmax=None)

        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 35
            cbar.set_label("Pearson r", fontsize=17, rotation=270, color='white')
            cbar.ax.tick_params(labelsize=15)
        for ax in axes:
            ax.set_title(f'SBC network {i}', fontsize=30, color='white')

    for i in range(len(NPR_maps)):
        axes = axes2[i+4+dr_maps.shape[0]+len(SBC_maps), :]

        map = NPR_maps[i, :]
        threshold = threshold_top_percent(map, top_percent=brainmap_percent_threshold)
        mask=np.abs(map)>=threshold # taking absolute values to include negative weights
        mask_img = recover_3D(mask_file,mask)        
        sitk_img = recover_3D(
            mask_file, map)
        cbar_list = masked_plot(fig2,axes, sitk_img, scaled, mask_img=mask_img, vmax=None)

        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 35
            cbar.set_label("Beta \nCoefficient", fontsize=17, rotation=270, color='white')
            cbar.ax.tick_params(labelsize=15)
        for ax in axes:
            ax.set_title(f'NPR network {i}', fontsize=30, color='white')

    return fig2
