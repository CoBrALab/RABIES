import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
import nilearn
from rabies.visualization import otsu_scaling, plot_3d
from rabies.analysis_pkg.analysis_math import elementwise_spearman, dice_coefficient
from rabies.utils import recover_3D
from rabies.confound_correction_pkg.utils import smooth_image
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
    volume_indices=sitk.GetArrayFromImage(sitk.ReadImage(mask_file)).astype(bool)    

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
        import nibabel as nb
        affine = nb.load(mask_file).affine[:3,:3]
        prior = sitk.GetArrayFromImage(smooth_image(recover_3D(mask_file,prior), affine, 0.3))[volume_indices]
        average = sitk.GetArrayFromImage(smooth_image(recover_3D(mask_file,average), affine, 0.3))[volume_indices]
        corr_map_std = sitk.GetArrayFromImage(smooth_image(recover_3D(mask_file,corr_map_std), affine, 0.3))[volume_indices]
        corr_map_CR_std = sitk.GetArrayFromImage(smooth_image(recover_3D(mask_file,corr_map_CR_std), affine, 0.3))[volume_indices]
        corr_map_VE = sitk.GetArrayFromImage(smooth_image(recover_3D(mask_file,corr_map_VE), affine, 0.3))[volume_indices]
        if np.array(tdof_list).std()==0:
            corr_map_tdof=None
        else:
            corr_map_tdof = sitk.GetArrayFromImage(smooth_image(recover_3D(mask_file,corr_map_tdof), affine, 0.3))[volume_indices]
    
    return prior, average, network_var, corr_map_std, corr_map_CR_std, corr_map_VE, corr_map_tdof    


def eval_relationships(maps, mask_file, percentile=0.01):
    prior, average, network_var, corr_map_std, corr_map_CR_std, corr_map_VE, corr_map_tdof = maps
 
    tmppath = tempfile.mkdtemp()
    map_masks=[]
    for map in maps:
        if map is None:
            map_masks.append(None)
            tdof_avg_corr=None
            continue
        img = recover_3D(mask_file,map)
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
    
    img = recover_3D(mask_file,prior)
    ax=axes[0]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=None, percentile=percentile)
    ax.set_title('Prior network', fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Prior measure", fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    img = recover_3D(mask_file,average)
    ax=axes[1]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=None, percentile=percentile)
    ax.set_title('Dataset average', fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Mean", fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    img = recover_3D(mask_file,network_var)
    ax=axes[2]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=None, percentile=percentile)
    ax.set_title('Dataset MAD', fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Median Absolute \nDeviation", fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    img = recover_3D(mask_file,corr_map_std)
    ax=axes[3]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=1.0, percentile=percentile)
    ax.set_title('$\mathregular{BOLD_{SD}}$ X network corr.', fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Spearman rho", fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    img = recover_3D(mask_file,corr_map_CR_std)
    ax=axes[4]
    cbar_list = masked_plot(fig,ax, img, scaled, vmax=1.0, percentile=percentile)
    ax.set_title('$\mathregular{CR_{SD}}$ X network corr.', fontsize=30, color='white')
    for cbar in cbar_list:
        cbar.ax.get_yaxis().labelpad = 20
        cbar.set_label("Spearman rho", fontsize=17, rotation=270, color='white')
        cbar.ax.tick_params(labelsize=15)

    ax=axes[5]
    if corr_map_tdof is None:
        ax.axis('off')
    else:
        img = recover_3D(mask_file,corr_map_tdof)
        cbar_list = masked_plot(fig,ax, img, scaled, vmax=1.0, percentile=percentile)
        ax.set_title('tDOF correlation', fontsize=30, color='white')
        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 20
            cbar.set_label("Spearman rho", fontsize=17, rotation=270, color='white')
            cbar.ax.tick_params(labelsize=15)

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
