import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
import nilearn
from rabies.visualization import otsu_scaling, plot_3d
from rabies.analysis_pkg.analysis_math import elementwise_spearman, elementwise_corrcoef, dice_coefficient
from rabies.utils import recover_3D
from rabies.confound_correction_pkg.utils import smooth_image
import tempfile


def analysis_QC(FC_maps, consensus_network, mask_file, corr_variable, variable_name, template_file, non_parametric=False):

    scaled = otsu_scaling(template_file)
        
    smoothing=True
    if non_parametric:
        name_list = ['Prior network', 'Dataset Median', 'Dataset MAD']+variable_name
        measure_list = ["Median", "Median Absolute \nDeviation", "Spearman rho"]
    else:
        name_list = ['Prior network', 'Dataset Mean', 'Dataset STD']+variable_name
        measure_list = ["Mean", "Standard \nDeviation", "Pearson r"]
    
    maps = get_maps(consensus_network, FC_maps, corr_variable, mask_file, smoothing, non_parametric=non_parametric)
    dataset_stats, map_masks=eval_relationships(maps, name_list)

    fig = plot_relationships(mask_file, scaled, maps, map_masks, name_list, measure_list, thresholded=True)
    fig_unthresholded = plot_relationships(mask_file, scaled, maps, map_masks, name_list, measure_list, thresholded=False)

    return dataset_stats, fig, fig_unthresholded

    
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


def eval_relationships(maps, name_list):
    dataset_stats = {}
    map_masks=[]
    
    for i in range(len(maps)):
        map = maps[i]
        threshold = percent_threshold(map)
        mask=np.abs(map)>=threshold # taking absolute values to include negative weights
        map_masks.append(mask)
        if i>0: #once we're past the prior, evaluate Dice relative to it
            dataset_stats[f'Overlap: Prior - {name_list[i]}'] = dice_coefficient(map_masks[0],mask)
        if i>2: # once we're past the MAD metric, for the corr maps, get an average correlation within the maks
            dataset_stats[f'Avg.: {name_list[i]}'] = map[map_masks[0]].mean()    

    return dataset_stats, map_masks


def plot_relationships(mask_file, scaled, maps, map_masks, name_list, measure_list, thresholded=True):

    nrows = len(name_list)
    fig,axes = plt.subplots(nrows=nrows, ncols=1,figsize=(12,2*nrows))
        
    img = recover_3D(mask_file,maps[0])
    if thresholded:
        mask_img = recover_3D(mask_file,map_masks[0])
    else:
        mask_img = sitk.ReadImage(mask_file)
    ax=axes[0]
    cbar_list = masked_plot(fig,ax, img, scaled, mask_img=mask_img, vmax=None)
    ax.set_title('Prior network', fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Prior measure", fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    img = recover_3D(mask_file,maps[1])
    if thresholded:
        mask_img = recover_3D(mask_file,map_masks[1])
    else:
        mask_img = sitk.ReadImage(mask_file)
    ax=axes[1]
    cbar_list = masked_plot(fig,ax, img, scaled, mask_img=mask_img, vmax=None)
    ax.set_title(name_list[1], fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label(measure_list[0], fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    img = recover_3D(mask_file,maps[2])
    if thresholded:
        mask_img = recover_3D(mask_file,map_masks[2])
    else:
        mask_img = sitk.ReadImage(mask_file)
    ax=axes[2]
    cbar_list = masked_plot(fig,ax, img, scaled, mask_img=mask_img, vmax=None)
    ax.set_title(name_list[2], fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label(measure_list[1], fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    for i in range(3,len(maps)):

        img = recover_3D(mask_file,maps[i])
        if thresholded:
            mask=np.abs(maps[i])>=0.1 # we set the min correlation at 0.1
            mask_img = recover_3D(mask_file,mask)
        else:
            mask_img = sitk.ReadImage(mask_file)
        ax=axes[i]
        cbar_list = masked_plot(fig,ax, img, scaled, mask_img=mask_img, vmax=0.5)
        ax.set_title(f'{name_list[i]} X network corr.', fontsize=30, color='white')
        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 20
            cbar.set_label(measure_list[2], fontsize=17, rotation=270, color='white')
            cbar.ax.tick_params(labelsize=15)

    plt.tight_layout()

    return fig

'''
def threshold_distribution(array):
    x = array.copy()
    med = np.median(x)
    x -= med # center around median
    l2_std = np.sqrt(np.mean(x**2)) # take L2-norm after centering; similar to STD
    threshold = (l2_std*2)+med # set threshold as 2 STD away from the median
    return threshold
'''

def percent_threshold(array): # set threshold to be the top 4% of all voxels
    flat=array.flatten()
    flat.sort()
    idx=int((0.96)*len(flat))
    threshold = flat[idx]
    return threshold

def masked_plot(fig,axes, img, scaled, mask_img, vmax=None):
    masked = img*mask_img
    
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
