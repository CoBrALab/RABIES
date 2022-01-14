import os
import numpy as np
import nibabel as nb
import pandas as pd
import matplotlib.pyplot as plt
from rabies.analysis_pkg import analysis_functions
import SimpleITK as sitk
import nilearn.plotting
from .analysis_QC import masked_plot

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
    rc = run_command(command)
    return out_file


def compute_edge_mask(in_mask, out_file, num_edge_voxels=1):
    #custom function for computing edge mask from an input brain mask
    img = nb.load(in_mask)
    mask_array = np.asarray(img.dataobj)
    shape = mask_array.shape

    #iterate through all voxels from the three dimensions and look if it contains surrounding voxels
    edge_mask = np.zeros(shape, dtype=bool)
    num_voxel = 0
    while num_voxel < num_edge_voxels:
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    #only look if the voxel is part of the mask
                    if mask_array[x, y, z]:
                        if (mask_array[x-1:x+2, y-1:y+2, z-1:z+2] == 0).sum() > 0:
                            edge_mask[x, y, z] = 1
        mask_array = mask_array-edge_mask
        num_voxel += 1

    nb.Nifti1Image(edge_mask, img.affine, img.header).to_filename(out_file)


'''
Prepare the subject data
'''


def process_data(bold_file, data_dict, VE_file, STD_file, CR_STD_file, mask_file_dict, analysis_dict, prior_bold_idx, prior_confound_idx, dual_ICA=0):
    temporal_info = {}
    spatial_info = {}

    FD_trace = data_dict['FD_trace']
    DVARS = data_dict['DVARS']
    temporal_info['FD_trace'] = FD_trace
    temporal_info['DVARS'] = DVARS
    temporal_info['VE_temporal'] = data_dict['VE_temporal']

    WM_mask = np.asarray(
        nb.load(mask_file_dict['WM_mask']).dataobj).astype(bool)
    CSF_mask = np.asarray(
        nb.load(mask_file_dict['CSF_mask']).dataobj).astype(bool)

    edge_mask = np.asarray(
        nb.load(mask_file_dict['edge_mask']).dataobj).astype(bool)
    brain_mask = np.asarray(nb.load(mask_file_dict['brain_mask']).dataobj)
    volume_indices = brain_mask.astype(bool)
    edge_idx = edge_mask[volume_indices]
    WM_idx = WM_mask[volume_indices]
    CSF_idx = CSF_mask[volume_indices]
    not_edge_idx = (edge_idx == 0)*(WM_idx == 0)*(CSF_idx == 0)


    data_array = np.asarray(nb.load(bold_file).dataobj)
    timeseries = np.zeros([data_array.shape[3], volume_indices.sum()])
    for i in range(data_array.shape[3]):
        timeseries[i, :] = (data_array[:, :, :, i])[volume_indices]

    '''Temporal Features'''
    DR_W = np.array(pd.read_csv(analysis_dict['dual_regression_timecourse_csv'], header=None))
    DR_array = np.asarray(nb.load(analysis_dict['dual_regression_nii']).dataobj)
    DR_C = np.zeros([DR_array.shape[3], volume_indices.sum()])
    for i in range(DR_array.shape[3]):
        DR_C[i, :] = (DR_array[:, :, :, i])[volume_indices]

    temporal_info['DR_all'] = DR_W

    signal_trace = np.abs(DR_W[:, prior_bold_idx]).mean(axis=1)
    noise_trace = np.abs(DR_W[:, prior_confound_idx]).mean(axis=1)
    temporal_info['signal_trace'] = signal_trace
    temporal_info['noise_trace'] = noise_trace

    # take regional timecourse from L2-norm
    global_trace = np.sqrt((timeseries.T**2).mean(axis=0))
    WM_trace = np.sqrt((timeseries.T[WM_idx]**2).mean(axis=0))
    CSF_trace = np.sqrt((timeseries.T[CSF_idx]**2).mean(axis=0))
    edge_trace = np.sqrt((timeseries.T[edge_idx]**2).mean(axis=0))
    not_edge_trace = np.sqrt((timeseries.T[not_edge_idx]**2).mean(axis=0))
    temporal_info['global_trace'] = global_trace
    temporal_info['WM_trace'] = WM_trace
    temporal_info['CSF_trace'] = CSF_trace
    temporal_info['edge_trace'] = edge_trace
    temporal_info['not_edge_trace'] = not_edge_trace
    temporal_info['predicted_time'] = data_dict['predicted_time']

    '''Spatial Features'''
    global_signal = timeseries.mean(axis=1)
    GS_corr = analysis_functions.vcorrcoef(timeseries.T, global_signal)
    DVARS_corr = analysis_functions.vcorrcoef(timeseries.T[:, 1:], DVARS[1:])
    FD_corr = analysis_functions.vcorrcoef(timeseries.T, np.asarray(FD_trace))

    prior_fit_out = {'C': [], 'W': []}
    if dual_ICA > 0:
        prior_fit_out['W'] = np.array(pd.read_csv(analysis_dict['dual_ICA_timecourse_csv'], header=None))
        C_array = np.asarray(nb.load(analysis_dict['dual_ICA_filename']).dataobj)
        C = np.zeros([C_array.shape[3], volume_indices.sum()])
        for i in range(C_array.shape[3]):
            C[i, :] = (C_array[:, :, :, i])[volume_indices]
        prior_fit_out['C'] = C

    all_IC_array = np.asarray(nb.load(mask_file_dict['prior_maps']).dataobj)
    all_IC_vectors = np.zeros([all_IC_array.shape[3], volume_indices.sum()])
    for i in range(all_IC_array.shape[3]):
        all_IC_vectors[i, :] = (all_IC_array[:, :, :, i])[volume_indices]

    spatial_info['prior_maps'] = all_IC_vectors[prior_bold_idx]
    spatial_info['DR_BOLD'] = DR_C[prior_bold_idx]
    spatial_info['DR_all'] = DR_C

    spatial_info['dual_ICA_maps'] = prior_fit_out['C']
    temporal_info['dual_ICA_time'] = prior_fit_out['W']

    spatial_info['VE_spatial'] = np.asarray(
        nb.load(VE_file).dataobj)[volume_indices]
    spatial_info['temporal_std'] = np.asarray(
        nb.load(STD_file).dataobj)[volume_indices]
    spatial_info['predicted_std'] = np.asarray(
        nb.load(CR_STD_file).dataobj)[volume_indices]
    spatial_info['GS_corr'] = GS_corr
    spatial_info['DVARS_corr'] = DVARS_corr
    spatial_info['FD_corr'] = FD_corr

    return temporal_info, spatial_info


def temporal_external_formating(temporal_info, file_dict):
    import os
    import pandas as pd
    import pathlib  # Better path manipulation
    bold_file = file_dict['bold_file']
    filename_split = pathlib.Path(
        bold_file).name.rsplit(".nii")

    dual_regression_timecourse_csv = os.path.abspath(
        filename_split[0]+'_dual_regression_timecourse.csv')
    pd.DataFrame(temporal_info['DR_all']).to_csv(
        dual_regression_timecourse_csv, header=False, index=False)
    if len(temporal_info['dual_ICA_time']) > 0:
        dual_ICA_timecourse_csv = os.path.abspath(
            filename_split[0]+'_dual_ICA_timecourse.csv')
        pd.DataFrame(temporal_info['dual_ICA_time']).to_csv(
            dual_ICA_timecourse_csv, header=False, index=False)
    else:
        dual_ICA_timecourse_csv = None

    del temporal_info['DR_all'], temporal_info['dual_ICA_time']

    temporal_info_csv = os.path.abspath(filename_split[0]+'_temporal_info.csv')
    pd.DataFrame(temporal_info).to_csv(temporal_info_csv)
    return temporal_info_csv, dual_regression_timecourse_csv, dual_ICA_timecourse_csv


def spatial_external_formating(spatial_info, file_dict):
    import os
    import pathlib  # Better path manipulation
    from rabies.analysis_pkg import analysis_functions
    mask_file = file_dict['mask_file']
    bold_file = file_dict['bold_file']
    filename_split = pathlib.Path(
        bold_file).name.rsplit(".nii")

    VE_filename = os.path.abspath(filename_split[0]+'_VE.nii.gz')
    analysis_functions.recover_3D(
        mask_file, spatial_info['VE_spatial']).to_filename(VE_filename)

    std_filename = os.path.abspath(filename_split[0]+'_tSTD.nii.gz')
    analysis_functions.recover_3D(
        mask_file, spatial_info['temporal_std']).to_filename(std_filename)

    predicted_std_filename = os.path.abspath(filename_split[0]+'_predicted_std.nii.gz')
    analysis_functions.recover_3D(
        mask_file, spatial_info['predicted_std']).to_filename(predicted_std_filename)

    GS_corr_filename = os.path.abspath(filename_split[0]+'_GS_corr.nii.gz')
    analysis_functions.recover_3D(
        mask_file, spatial_info['GS_corr']).to_filename(GS_corr_filename)

    DVARS_corr_filename = os.path.abspath(
        filename_split[0]+'_DVARS_corr.nii.gz')
    analysis_functions.recover_3D(
        mask_file, spatial_info['DVARS_corr']).to_filename(DVARS_corr_filename)

    FD_corr_filename = os.path.abspath(filename_split[0]+'_FD_corr.nii.gz')
    analysis_functions.recover_3D(
        mask_file, spatial_info['FD_corr']).to_filename(FD_corr_filename)

    DR_maps_filename = os.path.abspath(filename_split[0]+'_DR_maps.nii.gz')
    analysis_functions.recover_3D_multiple(
        mask_file, spatial_info['DR_all']).to_filename(DR_maps_filename)

    if len(spatial_info['dual_ICA_maps']) > 0:
        import numpy as np
        dual_ICA_filename = os.path.abspath(
            filename_split[0]+'_dual_ICA.nii.gz')
        analysis_functions.recover_3D_multiple(mask_file, np.array(
            spatial_info['dual_ICA_maps'])).to_filename(dual_ICA_filename)
    else:
        dual_ICA_filename = None

    return VE_filename, std_filename, predicted_std_filename, GS_corr_filename, DVARS_corr_filename, FD_corr_filename, DR_maps_filename, dual_ICA_filename


'''
Subject-level QC
'''


def grayplot_regional(timeseries_file, mask_file_dict, fig, ax):
    timeseries_4d = np.asarray(nb.load(timeseries_file).dataobj)

    WM_mask = np.asarray(
        nb.load(mask_file_dict['WM_mask']).dataobj).astype(bool)
    CSF_mask = np.asarray(
        nb.load(mask_file_dict['CSF_mask']).dataobj).astype(bool)
    right_hem_mask = np.asarray(
        nb.load(mask_file_dict['right_hem_mask']).dataobj).astype(bool)
    left_hem_mask = np.asarray(
        nb.load(mask_file_dict['left_hem_mask']).dataobj).astype(bool)

    grayplot_array = np.empty((0, timeseries_4d.shape[3]))
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
                grayplot_array, timeseries_4d[:, i, :, :][mask_indices[:, i, :]], axis=0)
            slice_alt = np.append(slice_alt, np.ones(
                mask_indices[:, i, :].sum())*token)
            token = not token

    vmax = grayplot_array.std()
    im = ax.imshow(grayplot_array, cmap='gray',
                   vmax=vmax, vmin=-vmax, aspect='auto')
    return im, slice_alt, region_mask_label


def grayplot(timeseries_file, mask_file_dict, fig, ax):
    brain_mask = np.asarray(nb.load(mask_file_dict['brain_mask']).dataobj)
    volume_indices = brain_mask.astype(bool)

    data_array = np.asarray(nb.load(timeseries_file).dataobj)
    timeseries = np.zeros([data_array.shape[3], volume_indices.sum()])
    for i in range(data_array.shape[3]):
        timeseries[i, :] = (data_array[:, :, :, i])[volume_indices]

    grayplot_array = timeseries.T

    vmax = grayplot_array.std()
    im = ax.imshow(grayplot_array, cmap='gray',
                   vmax=vmax, vmin=-vmax, aspect='auto')
    return im


def scan_diagnosis(bold_file, mask_file_dict, temporal_info, spatial_info, CR_data_dict, regional_grayplot=False):
    template_file = mask_file_dict['template_file']
    
    fig = plt.figure(figsize=(6, 18))
    #fig.suptitle(name, fontsize=30, color='white')
    
    ax0 = fig.add_subplot(3,1,1)
    ax1 = fig.add_subplot(12,1,5)
    ax1_ = fig.add_subplot(12,1,6)
    ax2 = fig.add_subplot(6,1,4)
    ax3 = fig.add_subplot(6,1,5)
    ax4 = fig.add_subplot(6,1,6)

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
        im = grayplot(bold_file, mask_file_dict, fig, ax0)

    ax0.set_ylabel('Voxels', fontsize=20)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_visible(False)
    ax0.axes.get_yaxis().set_ticks([])
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
    confounds_csv = CR_data_dict['confounds_csv']
    time_range = CR_data_dict['time_range']
    frame_mask = CR_data_dict['frame_mask']
    df = pd.read_csv(confounds_csv)
    # take proper subset of timepoints
    ax1.plot(x,df['mov1'].to_numpy()[time_range][frame_mask])
    ax1.plot(x,df['mov2'].to_numpy()[time_range][frame_mask])
    ax1.plot(x,df['mov3'].to_numpy()[time_range][frame_mask])
    ax1.legend(['translation 1', 'translation 2', 'translation 3'],
               loc='center left', bbox_to_anchor=(1.15, 0.5))
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax1_.plot(x,df['rot1'].to_numpy()[time_range][frame_mask])
    ax1_.plot(x,df['rot2'].to_numpy()[time_range][frame_mask])
    ax1_.plot(x,df['rot3'].to_numpy()[time_range][frame_mask])
    ax1_.legend(['rotation 1', 'rotation 2', 'rotation 3'],
                loc='center left', bbox_to_anchor=(1.15, 0.5))
    plt.setp(ax1_.get_xticklabels(), visible=False)
    ax1_.spines['right'].set_visible(False)
    ax1_.spines['top'].set_visible(False)

    y = temporal_info['FD_trace'].to_numpy()
    ax2.plot(x,y, 'r')
    ax2.set_ylabel('FD in mm', fontsize=15)
    DVARS = temporal_info['DVARS']
    DVARS[0] = None
    ax2_ = ax2.twinx()
    y2 = DVARS
    ax2_.plot(x,y2, 'b')
    ax2_.set_ylabel('DVARS', fontsize=15)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2_.spines['top'].set_visible(False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2_.get_xticklabels(), visible=False)
    ax2.legend(['Framewise \nDisplacement (FD)'
                ], loc='center left', bbox_to_anchor=(1.15, 0.6))
    ax2_.legend(['DVARS'
                ], loc='center left', bbox_to_anchor=(1.15, 0.4))

    ax3.plot(x,temporal_info['edge_trace'])
    ax3.plot(x,temporal_info['WM_trace'])
    ax3.plot(x,temporal_info['CSF_trace'])
    ax3.plot(x,temporal_info['predicted_time'])
    ax3.set_ylabel('Mask L2-norm', fontsize=15)
    ax3_ = ax3.twinx()
    ax3_.plot(x,temporal_info['VE_temporal'], 'darkviolet')
    ax3_.set_ylabel('CR R^2', fontsize=15)
    ax3_.spines['right'].set_visible(False)
    ax3_.spines['top'].set_visible(False)
    plt.setp(ax3_.get_xticklabels(), visible=False)
    ax3.legend(['Edge Mask', 'WM Mask', 'CSF Mask', 'CR prediction'
                ], loc='center left', bbox_to_anchor=(1.15, 0.7))
    ax3_.legend(['CR R^2'
                ], loc='center left', bbox_to_anchor=(1.15, 0.3))
    ax3_.set_ylim([0,1])

    y = temporal_info['signal_trace']
    ax4.plot(x,y)
    ax4.plot(x,temporal_info['noise_trace'])
    ax4.legend(['BOLD components', 'Confound components'
                ], loc='center left', bbox_to_anchor=(1.15, 0.5))
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.set_xlabel('Timepoint', fontsize=15)
    ax4.set_ylabel('Abs. Beta \ncoefficients (Avg.)', fontsize=15)

    dr_maps = spatial_info['DR_BOLD']
    mask_file = mask_file_dict['brain_mask']

    nrows = 6+dr_maps.shape[0]

    fig2, axes2 = plt.subplots(nrows=nrows, ncols=3, figsize=(12*3, 2*nrows))
    plt.tight_layout()

    from rabies.visualization import otsu_scaling, plot_3d

    axes = axes2[0, :]
    scaled = otsu_scaling(template_file)
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    temporal_std = spatial_info['temporal_std']
    analysis_functions.recover_3D(
        mask_file, temporal_std).to_filename('temp_img.nii.gz')
    sitk_img = sitk.ReadImage('temp_img.nii.gz')

    # select vmax at 95th percentile value
    vector = temporal_std.flatten()
    vector.sort()
    vmax = vector[int(len(vector)*0.95)]
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=0, vmax=vmax,
            cmap='inferno', alpha=1, cbar=True, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 35
        cbar.set_label('Standard \n Deviation', fontsize=15, rotation=270, color='white')
    for ax in axes:
        ax.set_title('BOLD-Temporal s.d.', fontsize=30, color='white')


    axes = axes2[1, :]
    scaled = otsu_scaling(template_file)
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    predicted_std = spatial_info['predicted_std']
    analysis_functions.recover_3D(
        mask_file, predicted_std).to_filename('temp_img.nii.gz')
    sitk_img = sitk.ReadImage('temp_img.nii.gz')

    # select vmax at 95th percentile value
    vector = predicted_std.flatten()
    vector.sort()
    vmax = vector[int(len(vector)*0.95)]
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=0, vmax=vmax,
            cmap='inferno', alpha=1, cbar=True, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 35
        cbar.set_label('Standard \n Deviation', fontsize=15, rotation=270, color='white')
    for ax in axes:
        ax.set_title('CR-Temporal s.d.', fontsize=30, color='white')


    axes = axes2[2, :]
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    analysis_functions.recover_3D(
        mask_file, spatial_info['VE_spatial']).to_filename('temp_img.nii.gz')
    sitk_img = sitk.ReadImage('temp_img.nii.gz')
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=0, vmax=1, cmap='inferno',
            alpha=1, cbar=True, threshold=0.1, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label('CR R^2', fontsize=15, rotation=270, color='white')
    for ax in axes:
        ax.set_title('CR R^2', fontsize=30, color='white')

    axes = axes2[3, :]
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    analysis_functions.recover_3D(
        mask_file, spatial_info['GS_corr']).to_filename('temp_img.nii.gz')
    sitk_img = sitk.ReadImage('temp_img.nii.gz')
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=-1, vmax=1, cmap='cold_hot',
            alpha=1, cbar=True, threshold=0.1, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Pearson's' r", fontsize=15, rotation=270, color='white')
    for ax in axes:
        ax.set_title('Global Signal Correlation', fontsize=30, color='white')

    axes = axes2[4, :]
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    analysis_functions.recover_3D(
        mask_file, spatial_info['DVARS_corr']).to_filename('temp_img.nii.gz')
    sitk_img = sitk.ReadImage('temp_img.nii.gz')
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=-1, vmax=1, cmap='cold_hot',
            alpha=1, cbar=True, threshold=0.1, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Pearson's' r", fontsize=15, rotation=270, color='white')
    for ax in axes:
        ax.set_title('DVARS Correlation', fontsize=30, color='white')

    axes = axes2[5, :]
    plot_3d(axes, scaled, fig2, vmin=0, vmax=1,
            cmap='gray', alpha=1, cbar=False, num_slices=6)
    analysis_functions.recover_3D(
        mask_file, spatial_info['FD_corr']).to_filename('temp_img.nii.gz')
    sitk_img = sitk.ReadImage('temp_img.nii.gz')
    cbar_list = plot_3d(axes, sitk_img, fig2, vmin=-1, vmax=1, cmap='cold_hot',
            alpha=1, cbar=True, threshold=0.1, num_slices=6)
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Pearson's' r", fontsize=15, rotation=270, color='white')
    for ax in axes:
        ax.set_title('FD Correlation', fontsize=30, color='white')

    for i in range(dr_maps.shape[0]):
        axes = axes2[i+6, :]

        analysis_functions.recover_3D(
            mask_file, dr_maps[i, :]).to_filename('temp_img.nii.gz')
        sitk_img = sitk.ReadImage('temp_img.nii.gz')
        cbar_list = masked_plot(fig2,axes, sitk_img, scaled, percentile=0.015, vmax=None)

        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 35
            cbar.set_label("Beta \nCoefficient", fontsize=15, rotation=270, color='white')
        for ax in axes:
            ax.set_title(f'BOLD component {i}', fontsize=30, color='white')

    return fig, fig2

