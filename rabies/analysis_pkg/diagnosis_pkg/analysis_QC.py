import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
import nilearn
from rabies.visualization import otsu_scaling, plot_3d
from rabies.analysis_pkg.analysis_math import elementwise_spearman, elementwise_corrcoef, dice_coefficient
from rabies.utils import recover_3D
from rabies.confound_correction_pkg.utils import smooth_image
import tempfile


def analysis_QC(FC_maps, consensus_network, mask_file, corr_variable, variable_name, template_file, non_parametric=False, top_percent=10, smoothing=False):

    scaled = otsu_scaling(template_file)
        
    if non_parametric:
        name_list = ['Prior network', 'Group average', 'Cross-scan variability']+variable_name
        measure_list = ["Median", "Median Absolute \nDeviation", "Spearman rho"]
    else:
        name_list = ['Prior network', 'Group average', 'Cross-scan variability']+variable_name
        measure_list = ["Mean", "Standard \nDeviation", "Pearson r"]
    
    maps = get_maps(consensus_network, FC_maps, corr_variable, mask_file, smoothing, non_parametric=non_parametric)
    dataset_stats, map_masks=eval_relationships(maps, name_list, top_percent=top_percent)

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
        mask_img = sitk.ReadImage(mask_file)
        for i in range(len(maps)):
            maps[i] = sitk.GetArrayFromImage(smooth_image(recover_3D(mask_file,maps[i]), 0.3, mask_img))[volume_indices]
    
    return maps


def eval_relationships(maps, name_list, top_percent=10):
    dataset_stats = {}
    map_masks=[]
    
    for i in range(len(maps)):
        map = maps[i]
        threshold = threshold_top_percent(map, top_percent=top_percent)
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


def threshold_top_percent(array, top_percent=10):
    # we take absolute values, since high negative values are also counted
    flat=np.abs(array).flatten()
    flat.sort()
    idx=int((1 - (top_percent/100))*len(flat))
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



###########################################
# PLOTTING DISTRIBUTIONS TO DETECT OUTLIERS
###########################################

from scipy import stats
import pandas as pd


def plot_density_2D(x,y, cm, ax, xlim, ylim):
    xmin, xmax = xlim
    ymin, ymax = ylim

    # Peform the kernel density estimate
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = stats.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    
    f[f<f.max()*0.2] = np.nan

    ax.contourf(xx, yy, f, cmap=cm)


def set_bounds(v,edge=1.5):
    half_range_=(v.max()-v.min())/2
    center=v.min()+half_range_
    min_ = center - (half_range_*edge)
    max_ = center + (half_range_*edge)
    return min_,max_


def detect_outliers(v, threshold=3.5):
    # taken from https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
    # v is a univariate vector for one variable
    v_ = v-np.median(v)
    mad = np.median(np.abs(v_))
    modified_z = (v_*0.6745)/mad
    return np.abs(modified_z)>threshold


def plot_density(v, bounds, outliers, ax, axis='x'):
    for outlier in [False,True]:
        idx = outliers==outlier
        if idx.sum()>0:
            v_=v[idx]
            if outlier:
                if axis=='x':
                    ax.bar(x=v_, height=ax.get_ylim()[1]*0.7, width=(bounds[1]-bounds[0])/70, color='orange')
                elif axis=='y':
                    ax.barh(y=v_, height=(bounds[1]-bounds[0])/70, width=ax.get_xlim()[1]*0.7, color='orange')
                else:
                    raise
            else:
                xs = np.linspace(bounds[0],bounds[1],200)
                density = stats.gaussian_kde(v_)
                if axis=='x':
                    ax.plot(xs,density(xs))
                    # here we append a 0 at each extreme to make sure the filling goes down to the axis at 0
                    ax.fill_between(np.append(np.append(xs.min(),xs),xs.max()), 
                                    np.append(np.append(0,density(xs)),0),
                                    alpha=0.3)
                    ax.set_ylim([0,set_bounds(density(xs),edge=1.2)[1]])
                elif axis=='y':
                    ax.plot(density(xs), xs)
                    # here we append a 0 at each extreme to make sure the filling goes down to the axis at 0
                    ax.fill_between(np.append(np.append(0,density(xs)),0),
                                    np.append(np.append(xs.min(),xs),xs.max()), 
                                    alpha=0.3)
                    ax.set_xlim([0,set_bounds(density(xs),edge=1.2)[1]])
                else:
                    raise


def plot_QC_distributions(net_metric_dict,qc_metric_dict, x_bounds_dict={}, y_bounds_dict={}, QC_inclusion=None, scan_QC_thresholds={}, outlier_threshold=3.5, write_corr=False):

    n_cols=len(net_metric_dict.keys())
    n_rows=len(qc_metric_dict.keys())

    if QC_inclusion is None:
        key, value = next(iter(net_metric_dict.items()))
        QC_inclusion = np.ones(len(value)).astype(bool)


    if n_cols==0 or n_rows==0:
        raise ValueError("Must have at minimum 1 row and 1 column.")

    height_ratios = [1]+[2]*n_rows # top row half height
    width_ratios = [2]*n_cols+[1] # last column half width
    fig, axes = plt.subplots(
        n_rows+1, n_cols+1,
        gridspec_kw={
            "height_ratios": height_ratios,
            "width_ratios": width_ratios
        },
        figsize=(np.array(width_ratios).sum()*3.5, np.array(height_ratios).sum()*3.5)
    )
    axes[0,-1].axis('off') # remove top right corner

    for ax in axes.flatten():
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.setp(ax.get_xticklabels(), fontsize=15)
        plt.setp(ax.get_yticklabels(), fontsize=15)


    for j, key_net in enumerate(net_metric_dict):
        axes[-1,j].set_xlabel(key_net, color='White', fontsize=25)
        if key_net in x_bounds_dict:
            x_bounds = x_bounds_dict[key_net]
        else:
            x_bounds = list(set_bounds(net_metric_dict[key_net]))
        for i in range(n_rows+1):
            axes[i,j].set_xlim(x_bounds)

    for i, key_qc in enumerate(qc_metric_dict):
        axes[i+1,0].set_ylabel(key_qc, color='White', fontsize=25)
        if key_qc in y_bounds_dict:
            y_bounds = y_bounds_dict[key_qc]
        else:
            y_bounds = list(set_bounds(qc_metric_dict[key_qc]))
        for j in range(n_cols+1):
            axes[i+1,j].set_ylim(y_bounds)

    # plot the thresholds for QC metrics
    for i, key_qc in enumerate(qc_metric_dict):
        if key_qc in scan_QC_thresholds: # if a threshold is set, plot the gray line
            for j, key_net in enumerate(net_metric_dict):
                ax = axes[i+1,j]
                x_bounds = list(ax.get_xlim())
                ax.plot(x_bounds,[scan_QC_thresholds[key_qc],scan_QC_thresholds[key_qc]], color='lightgray', linestyle='--')

    # plot the thresholds for network metrics
    for j, key_net in enumerate(net_metric_dict):
        if key_net in scan_QC_thresholds: # if a threshold is set, plot the gray line
            for i, key_qc in enumerate(qc_metric_dict):
                ax = axes[i+1,j]
                y_bounds = list(ax.get_ylim())
                ax.plot([scan_QC_thresholds[key_net],scan_QC_thresholds[key_net]], y_bounds, color='lightgray', linestyle='--')


    for j, key_net in enumerate(net_metric_dict):
        x = net_metric_dict[key_net]
        # don't include removed scans in side plots
        x_in = x[QC_inclusion]
        x_outliers = detect_outliers(x_in, threshold=outlier_threshold)

        # density plot for x axis along first row
        ax = axes[0,j]
        x_bounds = list(ax.get_xlim())
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        try:
            plot_density(v=x_in, bounds=x_bounds, outliers=x_outliers, ax=ax, axis='x')
        except np.linalg.LinAlgError:
            from nipype import logging
            log = logging.getLogger('nipype.workflow')
            log.warning("Singular matrix error when computing KDE. Density won't be shown.")

        for i, key_qc in enumerate(qc_metric_dict):
            y = qc_metric_dict[key_qc]
            # don't include removed scans in side plots
            y_in = y[QC_inclusion]
            y_outliers = detect_outliers(y_in, threshold=outlier_threshold)
            
            if j==0: # only need to plot the density during the first loop of rows
                ax = axes[i+1,-1]
                y_bounds = list(ax.get_ylim())
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                try:
                    plot_density(v=y_in, bounds=y_bounds, outliers=y_outliers, ax=ax, axis='y')
                except np.linalg.LinAlgError:
                    from nipype import logging
                    log = logging.getLogger('nipype.workflow')
                    log.warning("Singular matrix error when computing KDE. Density won't be shown.")

            # now manage the scatter/2D density plots
            ax = axes[i+1,j]
            union_outliers = x_outliers+y_outliers

            for outlier in [False,True]:
                idx = (union_outliers==outlier)
                x_=x_in[idx]
                y_=y_in[idx]
                if not outlier and idx.sum()>2: # only plot density if more than 2 data points
                    try:
                        plot_density_2D(x_,y_, cm='Blues', ax=ax, xlim=ax.get_xlim(), ylim=ax.get_ylim())
                    except np.linalg.LinAlgError:
                        from nipype import logging
                        log = logging.getLogger('nipype.workflow')
                        log.warning("Singular matrix error when computing KDE. Density won't be shown.")
                
                ax.scatter(x_,y_, s=30, edgecolor='black')
                if not outlier and write_corr:
                    corr = np.corrcoef(x_,y_)[0,1]
                    ax.text(
                        0.02, 0.98,
                        f"Scatter Correlation: {np.around(corr,3)}",
                        fontsize=20,
                        transform=ax.transAxes, # select coordinates of the axis
                        ha="left",
                        va="top"
                    )                    
                
            # show removed scans in gray
            ax.scatter(x[QC_inclusion==False],y[QC_inclusion==False], s=30, color='lightgray', edgecolor='black')

    return fig


def QC_distributions(prior_map,FC_maps,network_var,DR_conf_corr, total_CRsd, mean_FD_array, tdof_array, 
                     scan_name_list, scan_QC_thresholds, outlier_threshold=3.5, top_percent=10):

    threshold = threshold_top_percent(prior_map, top_percent=top_percent)
    prior_mask=np.abs(prior_map)>=threshold # taking absolute values to include negative weights
    dice_list=[]
    for i in range(len(scan_name_list)):
        threshold = threshold_top_percent(FC_maps[i,:], top_percent=top_percent)
        fit_mask=np.abs(FC_maps[i,:])>=threshold # taking absolute values to include negative weights
        dice = dice_coefficient(prior_mask,fit_mask)
        dice_list.append(dice)

    network_dice = np.array(dice_list)

    # apply QC thresholds
    QC_inclusion = np.ones(len(scan_name_list)).astype(bool)
    plot_QC_thresholds = {}
    if not scan_QC_thresholds['Dice'] is None:
        QC_inclusion *= (scan_QC_thresholds['Dice']<network_dice)
        plot_QC_thresholds['Network specificity (Dice)']=scan_QC_thresholds['Dice']
    if not scan_QC_thresholds['Conf'] is None:
        QC_inclusion *= (scan_QC_thresholds['Conf']>DR_conf_corr)
        plot_QC_thresholds['DR confound corr.\n(mean |pearson r|)']=scan_QC_thresholds['Conf']
    if scan_QC_thresholds['Amp']:
        while(True): # this is done iteratively, as new outliers may appear in the new distribution after removing a first set
            # we apply previous Dice/conf filters so removed scans don't contribute to estimating outliers
            Amp_QC = detect_outliers(network_var[QC_inclusion], threshold=outlier_threshold)==0
            if (Amp_QC==0).sum()==0: # if there are no more outliers, break from the loop
                break
            QC_inclusion[QC_inclusion] = Amp_QC
    
    net_metric_dict = {}
    if network_var is not None:
        net_metric_dict['Network amplitude']=np.array(network_var)
    net_metric_dict['Network specificity (Dice)']=np.array(network_dice)

    qc_metric_dict = {'DR confound corr.\n(mean |pearson r|)':np.array(DR_conf_corr), 
                      '$CR_{SD}$':np.array(total_CRsd), 
                      'Mean FD (mm)':np.array(mean_FD_array)}
    if tdof_array is not None: # exclude tdof_array if None
         qc_metric_dict['Degrees of freedom'] = np.array(tdof_array)
    x_bounds_dict = {'Network specificity (Dice)':[0,1.0]}
    y_bounds_dict = {'DR confound corr.\n(mean |pearson r|)':[0,1.0]}
    fig = plot_QC_distributions(
        net_metric_dict,qc_metric_dict, x_bounds_dict=x_bounds_dict, y_bounds_dict=y_bounds_dict, 
        QC_inclusion=QC_inclusion, scan_QC_thresholds=plot_QC_thresholds, outlier_threshold=outlier_threshold, 
        write_corr=True)

    df = pd.DataFrame(data=np.array([scan_name_list, QC_inclusion]).T, columns=['scan ID', 'QC inclusion?'])

    for feature,name in zip(
            [network_var,network_dice,DR_conf_corr, total_CRsd, mean_FD_array, tdof_array],
            ['Component variance','Dice overlap','DR confound corr.','$CR_{SD}$','mean FD','tDOF']):
        if feature is None:
            continue
        df[name] = feature

        outliers = np.empty(len(QC_inclusion))
        outliers[QC_inclusion] = detect_outliers(feature[QC_inclusion], threshold=outlier_threshold)
        
        outliers = outliers.astype(object) # convert to object to handle strings
        outliers[QC_inclusion==0]=''
        outliers[outliers==1]=True
        outliers[outliers==0]=False
        df[name+' - outlier?'] = outliers

    return fig,df,QC_inclusion
