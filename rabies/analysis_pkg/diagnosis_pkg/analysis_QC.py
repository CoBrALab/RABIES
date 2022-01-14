import SimpleITK as sitk
import pandas as pd
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
import nilearn
from rabies.visualization import otsu_scaling, plot_3d
from rabies.analysis_pkg.analysis_math import elementwise_spearman, dice_coefficient
from rabies.analysis_pkg.analysis_functions import recover_3D

import tempfile

def analysis_QC(FC_maps, consensus_network, mask_file, std_maps, CR_std_maps, VE_maps, tdof_list, template_file, fig_path):

    scaled = otsu_scaling(template_file)
        
    percentile=0.01
    smoothing=True
    
    maps = get_maps(consensus_network, FC_maps, std_maps, CR_std_maps, VE_maps, tdof_list, mask_file, smoothing)
    dataset_stats=eval_relationships(maps, mask_file, percentile=percentile)

    fig = plot_relationships(mask_file, scaled, maps, percentile=percentile)
    fig.savefig(fig_path, bbox_inches='tight')

    return dataset_stats

    
def get_maps(prior, prior_list, std_list, CR_std_maps, VE_list, tdof_list, mask_file, smoothing=False):
    volume_indices=np.asarray(nb.load(mask_file).dataobj).astype(bool)    

    Y=np.array(prior_list)
    average=Y.mean(axis=0)
    
    # compute MAD to be resistant to outliers
    mad = np.median(np.abs(Y-np.median(Y, axis=0)), axis=0)
    #network_var=Y.std(axis=0)
    network_var=mad
        
    X=np.array(std_list)
    corr_map_std = elementwise_spearman(X,Y)
    X=np.array(CR_std_maps)
    corr_map_CR_std = elementwise_spearman(X,Y)
    X=np.array(VE_list)
    corr_map_VE = elementwise_spearman(X,Y)
    
    # tdof effect; if there's no variability don't compute
    if np.array(tdof_list).std()==0:
        corr_map_tdof=None
    else:
        tdof = np.array(tdof_list).reshape(-1,1)
        corr_map_tdof = elementwise_spearman(tdof,Y)
    
    if smoothing:
        prior = np.array(nilearn.image.smooth_img(recover_3D(mask_file,prior), 0.3).dataobj)[volume_indices]
        average = np.array(nilearn.image.smooth_img(recover_3D(mask_file,average), 0.3).dataobj)[volume_indices]
        corr_map_std = np.array(nilearn.image.smooth_img(recover_3D(mask_file,corr_map_std), 0.3).dataobj)[volume_indices]
        corr_map_CR_std = np.array(nilearn.image.smooth_img(recover_3D(mask_file,corr_map_CR_std), 0.3).dataobj)[volume_indices]
        corr_map_VE = np.array(nilearn.image.smooth_img(recover_3D(mask_file,corr_map_VE), 0.3).dataobj)[volume_indices]
        if np.array(tdof_list).std()==0:
            corr_map_tdof=None
        else:
            corr_map_tdof = np.array(nilearn.image.smooth_img(recover_3D(mask_file,corr_map_tdof), 0.3).dataobj)[volume_indices]
    
    return prior, average, network_var, corr_map_std, corr_map_CR_std, corr_map_VE, corr_map_tdof
    

def eval_relationships(maps, mask_file, percentile=0.01):
    prior, average, network_var, corr_map_std, corr_map_CR_std, corr_map_VE, corr_map_tdof = maps

    # 1) the correlation between the prior and the fitted average
    #prior_cor = np.corrcoef(prior,average)[0,1]
    # 2) the weighted average of the temporal s.d. correlations
    #weights=prior/prior.sum()
    #weighted_std_corr=(weights*corr_map_std)[np.isnan(corr_map_std)==0].sum()
    #std_corr = np.corrcoef(prior[np.isnan(corr_map_std)==0], corr_map_std[np.isnan(corr_map_std)==0])[0,1]
    # 3) the correlation between the prior and the fitted average
    #VE_corr = np.corrcoef(prior[np.isnan(corr_map_VE)==0], corr_map_VE[np.isnan(corr_map_VE)==0])[0,1]
    #weighted_VE_corr=(weights*corr_map_VE)[np.isnan(corr_map_VE)==0].sum()
 
    tmppath = tempfile.mkdtemp()
    map_masks=[]
    for map in maps:
        if map is None:
            map_masks.append(None)
            tdof_avg_corr=None
            continue
        recover_3D(mask_file,map).to_filename(f'{tmppath}/temp_img.nii.gz')
        img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
        mask=percent_masking(img, percentile=percentile)
        map_masks.append(mask)
        
        # get also estimates of effect sizes for power estimation
        # defined by average correlation within the temporal s.d. map for temporal s.d./CR R^2
        if map is corr_map_std:
            std_avg_corr = sitk.GetArrayFromImage(img)[map_masks[0]].mean()
        if map is corr_map_CR_std:
            CR_std_avg_corr = sitk.GetArrayFromImage(img)[map_masks[0]].mean()
        if map is corr_map_VE:
            VE_avg_corr = sitk.GetArrayFromImage(img)[map_masks[0]].mean()
        if map is corr_map_tdof:
            tdof_avg_corr = sitk.GetArrayFromImage(img)[map_masks[0]].mean()
            
    prior_mask, average_mask, std_mask, corr_map_std_mask, corr_map_CR_std_mask, corr_map_VE_mask, tdof_mask = map_masks

    if tdof_mask is None:
        tdof_spec=None
    else:
        tdof_spec = dice_coefficient(prior_mask,tdof_mask)    

    return {'Overlap: Prior - Dataset avg.': dice_coefficient(prior_mask,average_mask), 
            'Overlap: Prior - Dataset MAD': dice_coefficient(prior_mask,std_mask),
            'Overlap: Prior - BOLD-Temporal s.d.': dice_coefficient(prior_mask,corr_map_std_mask), 
            'Overlap: Prior - CR-Temporal s.d.': dice_coefficient(prior_mask,corr_map_CR_std_mask), 
            'Overlap: Prior - tDOF': tdof_spec,
            'Avg.: BOLD-Temporal s.d.': std_avg_corr,
            'Avg.: CR-Temporal s.d.': CR_std_avg_corr,
            'Avg.: tDOF': tdof_avg_corr,
            } 


def plot_relationships(mask_file, scaled, maps, percentile=0.01):

    prior, average, network_var, corr_map_std, corr_map_CR_std, corr_map_VE, corr_map_tdof = maps

    fig,axes = plt.subplots(nrows=6, ncols=1,figsize=(12,2*6))
    
    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,prior).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[0]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=None, percentile=percentile)
    ax.set_title('Prior network', fontsize=25, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Prior measure", fontsize=12, rotation=270, color='white')

    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,average).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[1]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=None, percentile=percentile)
    ax.set_title('Dataset average', fontsize=25, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Mean", fontsize=12, rotation=270, color='white')

    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,network_var).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[2]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=None, percentile=percentile)
    ax.set_title('Dataset MAD', fontsize=25, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Median Absolute \nDeviation", fontsize=12, rotation=270, color='white')

    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,corr_map_std).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[3]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=1.0, percentile=percentile)
    ax.set_title('BOLD-Temporal s.d.', fontsize=25, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Spearman rho", fontsize=12, rotation=270, color='white')

    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,corr_map_CR_std).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[4]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=1.0, percentile=percentile)
    ax.set_title('CR-Temporal s.d.', fontsize=25, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Spearman rho", fontsize=12, rotation=270, color='white')

    ax=axes[5]
    if corr_map_tdof is None:
        ax.axis('off')
    else:
        tmppath = tempfile.mkdtemp()
        recover_3D(mask_file,corr_map_tdof).to_filename(f'{tmppath}/temp_img.nii.gz')
        img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
        cbar_list = masked_plot(fig,ax, img, scaled, vmax=1.0, percentile=percentile)
        ax.set_title('tDOF correlation', fontsize=25, color='white')
        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 20
            cbar.set_label("Spearman rho", fontsize=12, rotation=270, color='white')

    plt.tight_layout()

    return fig

def percent_masking(img, percentile):

    array=sitk.GetArrayFromImage(img)
    flat=array.flatten()
    flat.sort()
    idx=int((1-percentile)*len(flat))
    threshold = flat[idx]
    mask=array>=threshold

    return mask


def masked_plot(fig,axes, img, scaled, percentile=0.01, vmax=None):
    mask=percent_masking(img, percentile=percentile)
    
    masked=sitk.GetImageFromArray(sitk.GetArrayFromImage(img)*mask)
    masked.CopyInformation(img)
    
    data=sitk.GetArrayFromImage(img)
    if vmax is None:
        vmax=data.max()

    if not (type(axes) is np.ndarray or type(axes) is list):
        axes=[axes]
        planes = ('coronal')
    else:
        planes = ('sagittal', 'coronal', 'horizontal')
    plot_3d(axes,scaled,fig,vmin=0,vmax=1,cmap='gray', alpha=1, cbar=False, num_slices=6, planes=planes)
    # resample to match template
    sitk_img = sitk.Resample(masked, scaled)
    cbar_list = plot_3d(axes,sitk_img,fig,vmin=-vmax,vmax=vmax,cmap='cold_hot', alpha=1, cbar=True, threshold=vmax*0.001, num_slices=6, planes=planes)
    return cbar_list


def spatial_crosscorrelations(merged, scaled, mask_file, fig_path):

    dict_keys = ['temporal_std', 'predicted_std', 'VE_spatial', 'GS_corr',
                    'DVARS_corr', 'FD_corr', 'DR_BOLD', 'dual_ICA_maps']

    voxelwise_list = []
    for scan_data in merged:
        sub_list = [scan_data[key] for key in dict_keys]
        voxelwise_sub = np.array(sub_list[:6])
        if len(sub_list[7]) > 0:
            voxelwise_sub = np.concatenate(
                (voxelwise_sub, np.array(sub_list[6]), np.array(sub_list[7])), axis=0)
        else:
            voxelwise_sub = np.concatenate(
                (voxelwise_sub, np.array(sub_list[6])), axis=0)
        voxelwise_list.append(voxelwise_sub)
        num_prior_maps = len(sub_list[6])
    voxelwise_array = np.array(voxelwise_list)


    label_name = ['BOLD-Temporal s.d.', 'CR-Temporal s.d.', 'CR R^2',
                    'GS corr', 'DVARS corr', 'FD corr']
    label_name += [f'BOLD Dual Regression map {i}' for i in range(num_prior_maps)]
    label_name += [f'BOLD Dual ICA map {i}' for i in range(num_prior_maps)]

    ncols = 6
    fig, axes = plt.subplots(nrows=voxelwise_array.shape[1], ncols=ncols, figsize=(
        12*ncols, 2*voxelwise_array.shape[1]))
    for i, x_label in zip(range(voxelwise_array.shape[1]), label_name):
        for j, y_label in zip(range(ncols), label_name[:ncols]):
            ax = axes[i, j]
            if i <= j:
                ax.axis('off')
                continue

            X = voxelwise_array[:, i, :]
            Y = voxelwise_array[:, j, :]
            corr = elementwise_spearman(X, Y)

            plot_3d([ax], scaled, fig, vmin=0, vmax=1, cmap='gray',
                    alpha=1, cbar=False, num_slices=6, planes=('coronal'))
            recover_3D(
                mask_file, corr).to_filename('temp_img.nii.gz')
            sitk_img = sitk.ReadImage('temp_img.nii.gz')
            cbar_list = plot_3d([ax], sitk_img, fig, vmin=-1.0, vmax=1.0, cmap='cold_hot',
                    alpha=1, cbar=True, threshold=0.1, num_slices=6, planes=('coronal'))
            ax.set_title(f'Cross-correlation for \n{x_label} and {y_label}', fontsize=20, color='white')
            for cbar in cbar_list:
                cbar.ax.get_yaxis().labelpad = 20
                cbar.set_label("Spearman rho", fontsize=12, rotation=270, color='white')

    fig.savefig(fig_path,
                bbox_inches='tight')
