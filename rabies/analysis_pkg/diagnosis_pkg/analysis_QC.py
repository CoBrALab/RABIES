import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
import nilearn
from rabies.visualization import otsu_scaling, plot_3d
from rabies.analysis_pkg.analysis_math import elementwise_spearman, elementwise_corrcoef, dice_coefficient
from rabies.utils import recover_3D
from rabies.confound_correction_pkg.utils import smooth_image
import tempfile


def analysis_QC(FC_maps, consensus_network, mask_file, corr_variable, variable_name, template_file, fig_path, non_parametric=False):

    scaled = otsu_scaling(template_file)
        
    percentile=0.01
    smoothing=True
    if non_parametric:
        name_list = ['Prior network', 'Dataset Median', 'Dataset MAD']+variable_name
        measure_list = ["Median", "Median Absolute \nDeviation", "Spearman rho"]
    else:
        name_list = ['Prior network', 'Dataset Mean', 'Dataset STD']+variable_name
        measure_list = ["Mean", "Standard \nDeviation", "Pearson r"]
    
    maps = get_maps(consensus_network, FC_maps, corr_variable, mask_file, smoothing, non_parametric=non_parametric)
    dataset_stats=eval_relationships(maps, name_list, mask_file, percentile=percentile)

    fig = plot_relationships(mask_file, scaled, maps, name_list, measure_list, percentile=percentile)
    fig.savefig(fig_path, bbox_inches='tight')

    return dataset_stats

    
def get_maps(prior, prior_list, corr_variable, mask_file, smoothing=False, non_parametric=False):

    maps = []
    maps.append(prior)
    volume_indices=sitk.GetArrayFromImage(sitk.ReadImage(mask_file)).astype(bool)    

    Y=np.array(prior_list)
    if non_parametric:
        average=np.median(Y, axis=0)
        # compute MAD to be resistant to outliers
        mad = np.median(np.abs(Y-np.median(Y, axis=0)), axis=0)
        network_var=mad
    else:
        average=Y.mean(axis=0)
        network_var=Y.std(axis=0)
    maps.append(average)
    maps.append(network_var)
        
    for variable in corr_variable:    
        X=np.array(variable)
        if non_parametric:
            corr_map = elementwise_spearman(X,Y)
        else:
            corr_map = elementwise_corrcoef(X,Y)
        maps.append(corr_map)
        
    if smoothing:
        import nibabel as nb
        affine = nb.load(mask_file).affine[:3,:3]
        for i in range(len(maps)):
            maps[i] = sitk.GetArrayFromImage(smooth_image(recover_3D(mask_file,maps[i]), affine, 0.3))[volume_indices]
    
    return maps


def eval_relationships(maps, name_list, mask_file, percentile=0.01):
    dataset_stats = {}
    map_masks=[]
    
    for i in range(len(maps)):
        img = recover_3D(mask_file,maps[i])
        mask=percent_masking(img, percentile=percentile)
        map_masks.append(mask)
        if i>0: #once we're past the prior, evaluate Dice relative to it
            dataset_stats[f'Overlap: Prior - {name_list[i]}'] = dice_coefficient(map_masks[0],mask)
        if i>2: # once we're past the MAD metric, for the corr maps, get an average correlation within the maks
            dataset_stats[f'Avg.: {name_list[i]}'] = sitk.GetArrayFromImage(img)[map_masks[0]].mean()    

    return dataset_stats


def plot_relationships(mask_file, scaled, maps, name_list, measure_list, percentile=0.01):

    nrows = len(name_list)
    fig,axes = plt.subplots(nrows=nrows, ncols=1,figsize=(12,2*nrows))
        
    img = recover_3D(mask_file,maps[0])
    ax=axes[0]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=None, percentile=percentile)
    ax.set_title('Prior network', fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Prior measure", fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    img = recover_3D(mask_file,maps[1])
    ax=axes[1]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=None, percentile=percentile)
    ax.set_title(name_list[1], fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label(measure_list[0], fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    img = recover_3D(mask_file,maps[2])
    ax=axes[2]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=None, percentile=percentile)
    ax.set_title(name_list[2], fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label(measure_list[1], fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    for i in range(3,len(maps)):

        img = recover_3D(mask_file,maps[i])
        ax=axes[i]
        cbar_list = masked_plot(fig,ax, img, scaled, vmax=1.0, percentile=percentile)
        ax.set_title(f'{name_list[i]} X network corr.', fontsize=30, color='white')
        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 20
            cbar.set_label(measure_list[2], fontsize=17, rotation=270, color='white')
            cbar.ax.tick_params(labelsize=15)

    plt.tight_layout()

    return fig

def percent_masking(img, percentile):

    array=np.abs(sitk.GetArrayFromImage(img)) # taking absolute values to include negative weights
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
