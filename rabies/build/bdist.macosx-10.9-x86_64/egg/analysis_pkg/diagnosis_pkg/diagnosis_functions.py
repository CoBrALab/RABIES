import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rabies.utils import copyInfo_3DImage, recover_3D
from rabies.analysis_pkg import analysis_functions
import SimpleITK as sitk
import nilearn.plotting
from .analysis_QC import masked_plot, percent_threshold

def resample_mask(in_file, ref_file):
    transforms = []
    inverses = []
    # resampling the reference image to the dimension of the EPI
    from rabies.utils import run_command
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(
        in_file).name.rsplit(".nii")
    out_file = os.path.abspath(filename_split[0])+'_resampled.nii.gz'

    # tranforms is a list of transform files, set in order of call within antsApplyTransforms
    transform_string = ""
    for transform, inverse in zip(transforms, inverses):
        if transform == 'NULL':
            continue
        elif bool(inverse):
            transform_string += f"-t [{transform},1] "
        else:
            transform_string += f"-t {transform} "

    command = f'antsApplyTransforms -i {in_file} {transform_string}-n GenericLabel -r {ref_file} -o {out_file}'
    rc,c_out = run_command(command)
    return out_file


'''
Prepare the subject data
'''


def process_data(data_dict, analysis_dict, prior_bold_idx, prior_confound_idx):
    temporal_info = {}
    spatial_info = {}
    temporal_info['name_source'] = data_dict['name_source']
    spatial_info['name_source'] = data_dict['name_source']
    spatial_info['mask_file'] = data_dict['mask_file']

    timeseries = data_dict['timeseries']
    volume_indices = data_dict['volume_indices']
    edge_idx = data_dict['edge_idx']
    WM_idx = data_dict['WM_idx']
    CSF_idx = data_dict['CSF_idx']
    CR_data_dict = data_dict['CR_data_dict']

    not_edge_idx = (edge_idx == 0)*(WM_idx == 0)*(CSF_idx == 0)

    FD_trace = CR_data_dict['FD_trace']
    DVARS = CR_data_dict['DVARS']
    temporal_info['FD_trace'] = FD_trace
    temporal_info['DVARS'] = DVARS
    temporal_info['VE_temporal'] = CR_data_dict['VE_temporal']


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
        spatial_info['seed_map_list'] = []
        temporal_info['SBC_time'] = []

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
    # take regional timecourse from L2-norm
    WM_trace = np.sqrt((timeseries.T[WM_idx]**2).mean(axis=0))
    CSF_trace = np.sqrt((timeseries.T[CSF_idx]**2).mean(axis=0))
    edge_trace = np.sqrt((timeseries.T[edge_idx]**2).mean(axis=0))
    not_edge_trace = np.sqrt((timeseries.T[not_edge_idx]**2).mean(axis=0))
    temporal_info['WM_trace'] = WM_trace
    temporal_info['CSF_trace'] = CSF_trace
    temporal_info['edge_trace'] = edge_trace
    temporal_info['not_edge_trace'] = not_edge_trace
    temporal_info['predicted_time'] = CR_data_dict['predicted_time']

    '''Spatial Features'''
    global_signal = timeseries.mean(axis=1)
    GS_cov = (global_signal.reshape(-1,1)*timeseries).mean(axis=0) # calculate the covariance between global signal and each voxel
    GS_corr = analysis_functions.vcorrcoef(timeseries.T, global_signal)

    prior_fit_out = {'C': [], 'W': []}
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

    spatial_info['prior_maps'] = data_dict['prior_map_vectors'][prior_bold_idx]
    spatial_info['DR_BOLD'] = DR_C[prior_bold_idx]
    spatial_info['DR_all'] = DR_C

    spatial_info['NPR_maps'] = prior_fit_out['C']
    temporal_info['NPR_time'] = prior_fit_out['W']

    spatial_info['VE_spatial'] = data_dict['VE_spatial']
    spatial_info['temporal_std'] = data_dict['temporal_std']
    spatial_info['predicted_std'] = data_dict['predicted_std']
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


def grayplot_regional(timeseries_file, mask_file_dict, fig, ax):
    timeseries_4d = sitk.GetArrayFromImage(sitk.ReadImage(timeseries_file))

    WM_mask = sitk.GetArrayFromImage(
        sitk.ReadImage(mask_file_dict['WM_mask'])).astype(bool)
    CSF_mask = sitk.GetArrayFromImage(
        sitk.ReadImage(mask_file_dict['CSF_mask'])).astype(bool)
    right_hem_mask = sitk.GetArrayFromImage(
        sitk.ReadImage(mask_file_dict['right_hem_mask'])).astype(bool)
    left_hem_mask = sitk.GetArrayFromImage(
        sitk.ReadImage(mask_file_dict['left_hem_mask'])).astype(bool)

    grayplot_array = np.empty((0, timeseries_4d.shape[0]))
    slice_alt = np.array([])
    region_mask_label = np.array([])
    c = 0
    for mask_indices in [right_hem_mask, left_hem_mask, WM_mask, CSF_mask]:
        region_mask_label = np.append(
            region_mask_label, np.ones(mask_indices.sum())*c)
        c += 1
        token = False
        for i in range(mask_indices.shape[1]):
            grayplot_array = np.append(
                grayplot_array, timeseries_4d[:, :, i, :][mask_indices[:, i, :]], axis=0)
            slice_alt = np.append(slice_alt, np.ones(
                mask_indices[:, i, :].sum())*token)
            token = not token

    vmax = grayplot_array.std()
    im = ax.imshow(grayplot_array, cmap='gray',
                   vmax=vmax, vmin=-vmax, aspect='auto')
    return im, slice_alt, region_mask_label


def grayplot(timeseries, ax):
    grayplot_array = timeseries.T
    vmax = grayplot_array.std()
    im = ax.imshow(grayplot_array, cmap='gray',
                   vmax=vmax, vmin=-vmax, aspect='auto')
    return im


def plot_freqs(ax,timeseries, TR):
    freqs = np.fft.fftfreq(timeseries.shape[0], TR)
    idx = np.argsort(freqs)
    pos_idx = idx[freqs[idx]>0]

    ps = np.abs(np.fft.fft(timeseries.T))**2
    ps_mean = ps.mean(axis=0)
    ps_std = ps.std(axis=0)

    y_max = ps_mean[freqs>0.01].max()*1.5

    ax.plot(freqs[pos_idx], ps_mean[pos_idx])
    ax.fill_between(freqs[pos_idx], (ps_mean-ps_std)[pos_idx],(ps_mean+ps_std)[pos_idx], alpha=0.4)
    ax.set_ylim([0,y_max])
    ax.set_xlim([0,freqs.max()])
    xticks = np.arange(0,freqs.max(),0.05)
    ax.set_xticks(xticks)    


def scan_diagnosis(data_dict, temporal_info, spatial_info, regional_grayplot=False):
    timeseries = data_dict['timeseries']
    template_file = data_dict['template_file']
    CR_data_dict = data_dict['CR_data_dict']
    
    fig = plt.figure(figsize=(12, 24))

    ax0 = fig.add_subplot(4,2,3)
    ax0_f = fig.add_subplot(9,2,3) # plot frequencies

    ax1 = fig.add_subplot(16,2,17)
    ax1_ = fig.add_subplot(16,2,19)
    ax2 = fig.add_subplot(8,2,11)
    ax3 = fig.add_subplot(8,2,13)
    ax4 = fig.add_subplot(8,2,15)

    plot_freqs(ax0_f,timeseries, CR_data_dict['TR'])
    plt.setp(ax0_f.get_yticklabels(), visible=False)
    plt.setp(ax0_f.get_xticklabels(), fontsize=12)
    ax0_f.set_xlabel('Frequency (Hz)', fontsize=20)
    ax0_f.set_ylabel('Power (a.u.)', fontsize=20)

    # disable function
    regional_grayplot=False
    if regional_grayplot:
        
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax0)
        
        im, slice_alt, region_mask_label = grayplot_regional(
            bold_file, mask_file_dict, fig, ax0)
        ax0.yaxis.labelpad = 40
        ax_slice = divider.append_axes('left', size='5%', pad=0.0)
        ax_label = divider.append_axes('left', size='5%', pad=0.0)

        ax_slice.imshow(slice_alt.reshape(-1, 1), cmap='gray',
                        vmin=0, vmax=1.1, aspect='auto')
        ax_label.imshow(region_mask_label.reshape(-1, 1),
                        cmap='Spectral', aspect='auto')
        ax_slice.axis('off')
        ax_label.axis('off')

    else:
        im = grayplot(timeseries, ax0)

    ax0.set_ylabel('Voxels', fontsize=20)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_visible(False)
    ax0.axes.get_yaxis().set_ticks([])
    plt.setp(ax0.get_xticklabels(), fontsize=15)
    plt.setp(ax1.get_xticklabels(), visible=False)

    y = temporal_info['FD_trace'].to_numpy()
    x = range(len(y))
    ax0.set_xlim([0, len(y)-1])
    ax1.set_xlim([0, len(y)-1])
    ax1_.set_xlim([0, len(y)-1])
    ax2.set_xlim([0, len(y)-1])
    ax3.set_xlim([0, len(y)-1])
    ax4.set_xlim([0, len(y)-1])

    # plot the motion timecourses
    motion_params_csv = CR_data_dict['motion_params_csv']
    time_range = CR_data_dict['time_range']
    frame_mask = CR_data_dict['frame_mask']
    df = pd.read_csv(motion_params_csv)
    # take proper subset of timepoints
    ax1.plot(x,df['mov1'].to_numpy()[time_range][frame_mask])
    ax1.plot(x,df['mov2'].to_numpy()[time_range][frame_mask])
    ax1.plot(x,df['mov3'].to_numpy()[time_range][frame_mask])
    ax1.legend(['translation 1', 'translation 2', 'translation 3'],
               loc='center left', fontsize=15, bbox_to_anchor=(1.15, 0.5))
    ax1.set_ylabel('mm', fontsize=20)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax1_.plot(x,df['rot1'].to_numpy()[time_range][frame_mask])
    ax1_.plot(x,df['rot2'].to_numpy()[time_range][frame_mask])
    ax1_.plot(x,df['rot3'].to_numpy()[time_range][frame_mask])
    ax1_.legend(['Euler angle 1', 'Euler angle 2', 'Euler angle 3'],
                loc='center left', fontsize=15, bbox_to_anchor=(1.15, 0.5))
    ax1_.set_ylabel('radians', fontsize=20)
    plt.setp(ax1_.get_xticklabels(), visible=False)
    ax1_.spines['right'].set_visible(False)
    ax1_.spines['top'].set_visible(False)

    y = temporal_info['FD_trace'].to_numpy()
    ax2.plot(x,y, 'r')
    ax2.set_ylabel('FD (mm)', fontsize=20)
    DVARS = temporal_info['DVARS']
    DVARS[0] = None
    ax2_ = ax2.twinx()
    y2 = DVARS
    ax2_.plot(x,y2, 'b')
    ax2_.set_ylabel('DVARS', fontsize=20)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2_.spines['top'].set_visible(False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2_.get_xticklabels(), visible=False)
    ax2.legend(['Framewise \nDisplacement (FD)'
                ], loc='center left', fontsize=15, bbox_to_anchor=(1.15, 0.7))
    ax2_.legend(['DVARS'
                ], loc='center left', fontsize=15, bbox_to_anchor=(1.15, 0.3))

    ax3.plot(x,temporal_info['edge_trace'])
    ax3.plot(x,temporal_info['WM_trace'])
    ax3.plot(x,temporal_info['CSF_trace'])
    ax3.plot(x,temporal_info['predicted_time'])
    ax3.set_ylabel('Amplitude \n(L2-norm)', fontsize=20)
    ax3_ = ax3.twinx()
    ax3_.plot(x,temporal_info['VE_temporal'], 'darkviolet')
    ax3_.set_ylabel('CR $\mathregular{R^2}$', fontsize=20)
    ax3_.spines['right'].set_visible(False)
    ax3_.spines['top'].set_visible(False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3_.get_xticklabels(), visible=False)
    ax3.legend(['Edge Mask', 'WM Mask', 'CSF Mask', '$CR_{var}$'
                ], loc='center left', fontsize=15, bbox_to_anchor=(1.15, 0.7))
    ax3_.legend(['CR $\mathregular{R^2}$'
                ], loc='center left', fontsize=15, bbox_to_anchor=(1.15, 0.2))
    ax3_.set_ylim([0,1])

    # we take the mean of the timecourse amplitude (absolute values) to summarized across all components
    import matplotlib.cm as cm
    Greens_colors = cm.Greens(np.linspace(0.5, 0.8, 2))
    Blues_colors = cm.Blues(np.linspace(0.5, 0.8, 2))
    Purples_colors = cm.Purples(np.linspace(0.5, 0.8, 3))
    YlOrRd_colors = cm.YlOrRd(np.linspace(0.3, 0.8, 4))
    legend=[]
    for network_time,name,color,scaler in zip(
            [temporal_info['DR_confound'], temporal_info['DR_bold'],temporal_info['SBC_time'],temporal_info['NPR_time']],
            ['DR confounds', 'DR networks','SBC networks','NPR networks'],
            [YlOrRd_colors[2], Blues_colors[1],Purples_colors[1],Greens_colors[1]],
            [0,-1,-1,-1]):
        if len(network_time)>0:
            # make sure the timecourses are normalized
            network_time = network_time.copy()
            network_time /= np.sqrt((network_time ** 2).mean(axis=0))
            network_time_avg = np.abs(network_time).mean(axis=1)
            # we introduce a scaler to prevent overlap of confound with network timecourses
            network_time_avg += scaler
            ax4.plot(x,network_time_avg, color=color, alpha=0.6)
            legend.append(name)

    ax4.legend(legend, loc='center left', fontsize=15, bbox_to_anchor=(1.15, 0.5))
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.set_xlabel('Timepoint', fontsize=25)
    ax4.set_ylabel('Mean amplitude \n(a.u.)', fontsize=20)
    ax4.yaxis.set_ticklabels([])
    plt.setp(ax4.get_xticklabels(), fontsize=15)

    plt.setp(ax1.get_yticklabels(), fontsize=15)
    plt.setp(ax1_.get_yticklabels(), fontsize=15)
    plt.setp(ax2.get_yticklabels(), fontsize=15)
    plt.setp(ax2_.get_yticklabels(), fontsize=15)
    plt.setp(ax3.get_yticklabels(), fontsize=15)
    plt.setp(ax3_.get_yticklabels(), fontsize=15)
    plt.setp(ax4.get_yticklabels(), fontsize=15)


    dr_maps = spatial_info['DR_BOLD']
    SBC_maps = spatial_info['seed_map_list']
    NPR_maps = spatial_info['NPR_maps']
    mask_file = data_dict['mask_file']

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
        threshold = percent_threshold(map)
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
        threshold = percent_threshold(map)
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
        threshold = percent_threshold(map)
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

    return fig, fig2

