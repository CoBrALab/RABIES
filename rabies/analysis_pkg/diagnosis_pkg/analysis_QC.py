import os
import SimpleITK as sitk
import pandas as pd
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
import nilearn
from rabies.preprocess_pkg.preprocess_visual_QC import otsu_scaling, plot_3d
from rabies.analysis_pkg.analysis_math import elementwise_spearman, dice_coefficient
from rabies.analysis_pkg.analysis_functions import recover_3D

import tempfile
import shutil


def analysis_QC(FC_maps, consensus_network, mask_file, std_maps, VE_maps, template_file, fig_path):

    scaled = otsu_scaling(template_file)
        
    percentile=0.01
    smoothing=True
    threshold_spec=[[4,2],[4,2],[4,3],[4,3],[4,3]]
    
    maps = get_maps(consensus_network, FC_maps, std_maps, VE_maps, mask_file, smoothing)
    prior, average, network_var, corr_map_std, corr_map_VE = maps
    dataset_stats=eval_relationships(maps, mask_file, percentile=percentile, threshold_spec=threshold_spec)

    fig,axes = plt.subplots(nrows=5, ncols=1,figsize=(12,2*5))
    plot_relationships(fig,axes, mask_file, scaled, prior, average, network_var, corr_map_std, corr_map_VE, percentile=percentile, threshold_spec=threshold_spec)
    fig.savefig(fig_path, bbox_inches='tight')

    return dataset_stats

    
def get_maps(prior, prior_list, std_list, VE_list, mask_file, smoothing=False):
    
    volume_indices=np.asarray(nb.load(mask_file).dataobj).astype(bool)    
    
    Y=np.array(prior_list)
    average=Y.mean(axis=0)
    
    # compute MAD to be resistant to outliers
    mad = np.median(np.abs(Y-np.median(Y, axis=0)), axis=0)
    #network_var=Y.std(axis=0)
    network_var=mad
        
    X=np.array(std_list)
    corr_map_std = elementwise_spearman(X,Y)
    X=np.array(VE_list)
    corr_map_VE = elementwise_spearman(X,Y)
    
    if smoothing:
        prior = np.array(nilearn.image.smooth_img(recover_3D(mask_file,prior), 0.3).dataobj)[volume_indices]
        average = np.array(nilearn.image.smooth_img(recover_3D(mask_file,average), 0.3).dataobj)[volume_indices]
        corr_map_std = np.array(nilearn.image.smooth_img(recover_3D(mask_file,corr_map_std), 0.3).dataobj)[volume_indices]
        corr_map_VE = np.array(nilearn.image.smooth_img(recover_3D(mask_file,corr_map_VE), 0.3).dataobj)[volume_indices]
    
    return prior, average, network_var, corr_map_std, corr_map_VE
    

def eval_relationships(maps, mask_file, percentile=None, threshold_spec=[[4,2],[4,2],[4,2],[4,2],[4,2]], voxel_volume=0.2**3):
    prior, average, network_var, corr_map_std, corr_map_VE = maps

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
    for map,spec in zip(maps, threshold_spec):
        recover_3D(mask_file,map).to_filename(f'{tmppath}/temp_img.nii.gz')
        img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
        mask=masking(img, method='percent', percentile=percentile, threshold_spec=spec)
        map_masks.append(mask)
        
        # get also estimates of effect sizes for power estimation
        # defined by average correlation within the temporal s.d. map for temporal s.d./CR R^2
        if map is corr_map_std:
            std_avg_corr = sitk.GetArrayFromImage(img)[map_masks[0]].mean()
        if map is corr_map_VE:
            VE_avg_corr = sitk.GetArrayFromImage(img)[map_masks[0]].mean()
            
    prior_mask, average_mask, std_mask, corr_map_std_mask, corr_map_VE_mask = map_masks


    return {'criterion 1': dice_coefficient(prior_mask,average_mask), 
            'criterion 2': dice_coefficient(prior_mask,corr_map_std_mask), 
            'criterion 3': dice_coefficient(prior_mask,corr_map_VE_mask),
            'Network variability': dice_coefficient(prior_mask,std_mask),
            'Temporal s.d. avg. correlation': std_avg_corr,
            'CR R^2 avg. correlation': VE_avg_corr,
            'Temporal s.d. mask volume': voxel_volume*corr_map_std_mask.sum(),
            } 
    

def plot_relationships(fig,axes, mask_file, scaled, prior, average, network_var, corr_map_std, corr_map_VE, percentile=None, threshold_spec=[[4,2],[4,2],[4,2],[4,2],[4,2]]):
    
    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,prior).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[0]
    masked_plot(fig,[ax], img, scaled, hist=False, vmax=None, percentile=percentile, threshold_spec=threshold_spec[0])
    ax.set_title('Prior network', fontsize=25, color='white')

    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,average).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[1]
    masked_plot(fig,[ax], img, scaled, hist=False, vmax=None, percentile=percentile, threshold_spec=threshold_spec[1])
    ax.set_title('Dataset average', fontsize=25, color='white')

    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,network_var).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[2]
    masked_plot(fig,[ax], img, scaled, hist=False, vmax=None, percentile=percentile, threshold_spec=threshold_spec[2])
    ax.set_title('Network variability', fontsize=25, color='white')

    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,corr_map_std).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[3]
    masked_plot(fig,[ax], img, scaled, hist=False, vmax=0.7, percentile=percentile, threshold_spec=threshold_spec[3])
    ax.set_title('Temporal s.d. correlation', fontsize=25, color='white')

    tmppath = tempfile.mkdtemp()
    recover_3D(mask_file,corr_map_VE).to_filename(f'{tmppath}/temp_img.nii.gz')
    img = sitk.ReadImage(f'{tmppath}/temp_img.nii.gz')
    ax=axes[4]
    masked_plot(fig,[ax], img, scaled, hist=False, vmax=0.7, percentile=percentile, threshold_spec=threshold_spec[4])
    ax.set_title('CR R^2 correlation', fontsize=25, color='white')
    plt.tight_layout()


def masking(img, method='otsu', percentile=None, threshold_spec=[4,2]):
    if method=='otsu':
        mask=(otsu_mask(img, num_histograms=threshold_spec[0])>threshold_spec[1])
    elif method=='percent':
        if percentile is None:
            raise
        array=sitk.GetArrayFromImage(img)
        flat=array.flatten()
        flat.sort()
        idx=int((1-percentile)*len(flat))
        threshold = flat[idx]
        mask=array>=threshold
    else:
        raise
    return mask


def masked_plot(fig,axes, img, scaled, method='percent', percentile=None, hist=True, vmax=None, threshold_spec=[4,2]):
    mask=masking(img, method=method, percentile=percentile, threshold_spec=threshold_spec)
    
    masked=sitk.GetImageFromArray(sitk.GetArrayFromImage(img)*mask)
    masked.CopyInformation(img)
    
    data=sitk.GetArrayFromImage(img)
    if vmax is None:
        vmax=data.max()

    
    ax=axes[0]
    plot_3d([ax],scaled,fig,vmin=0,vmax=1,cmap='gray', alpha=1, cbar=False, num_slices=6, planes=('coronal'))
    # resample to match template
    sitk_img = sitk.Resample(masked, scaled)
    plot_3d([ax],sitk_img,fig,vmin=-vmax,vmax=vmax,cmap='cold_hot', alpha=1, cbar=True, threshold=vmax*0.001, num_slices=6, planes=('coronal'))
    
    if hist:
        ax=axes[1]
        ax.hist(np.abs(data).flatten(), density=True, bins=100)
        ax.set_ylim([0,1])
        ax.axvline(np.abs(data[mask]).min(), color='r')


def otsu_mask(img, num_histograms=1):
    tmppath = tempfile.mkdtemp()
    # we count both positive and negative values
    array=sitk.GetArrayFromImage(img)
    abs_array = np.abs(array)
    abs_img = sitk.GetImageFromArray(abs_array)
    abs_img.CopyInformation(img)
    sitk.WriteImage(abs_img, f'{tmppath}/temp_img.nii.gz')

    from rabies.preprocess_pkg.utils import run_command
    command = f'ThresholdImage 3 {tmppath}/temp_img.nii.gz {tmppath}/otsu_weight.nii.gz Otsu {num_histograms}'
    rc = run_command(command)
    
    mask=sitk.GetArrayFromImage(sitk.ReadImage(f'{tmppath}/otsu_weight.nii.gz'))
    shutil.rmtree(tmppath)
    return mask


def spatial_crosscorrelations(merged, scaled, mask_file, fig_path):

    dict_keys = ['temporal_std', 'VE_spatial', 'GS_corr',
                    'DVARS_corr', 'FD_corr', 'DR_BOLD', 'dual_ICA_maps']

    voxelwise_list = []
    for spatial_info in merged:
        sub_list = [spatial_info[key] for key in dict_keys]
        voxelwise_sub = np.array(sub_list[:5])
        if len(sub_list[6]) > 0:
            voxelwise_sub = np.concatenate(
                (voxelwise_sub, np.array(sub_list[5]), np.array(sub_list[6])), axis=0)
        else:
            voxelwise_sub = np.concatenate(
                (voxelwise_sub, np.array(sub_list[5])), axis=0)
        voxelwise_list.append(voxelwise_sub)
        num_prior_maps = len(sub_list[5])
    voxelwise_array = np.array(voxelwise_list)


    label_name = ['temporal_std', 'VE_spatial',
                    'GS_corr', 'DVARS_corr', 'FD_corr']
    label_name += [f'BOLD Dual Regression map {i}' for i in range(num_prior_maps)]
    label_name += [f'BOLD Dual ICA map {i}' for i in range(num_prior_maps)]

    ncols = 5
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
            plot_3d([ax], sitk_img, fig, vmin=-0.7, vmax=0.7, cmap='cold_hot',
                    alpha=1, cbar=True, threshold=0.1, num_slices=6, planes=('coronal'))
            ax.set_title(f'Cross-correlation for {x_label} and {y_label}', fontsize=15, color='white')

    fig.savefig(fig_path,
                bbox_inches='tight')
