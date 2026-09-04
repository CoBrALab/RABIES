import numpy as np

'''
SBC
'''

def compute_seed_FC(timeseries, seed_arr_dict):
    from .analysis_math import vcorrcoef
    SBC_dict = {}
    for seed_name, roi_mask_vec in seed_arr_dict.items():
        roi_mask = roi_mask_vec.astype(bool)
        seed_timeseries = timeseries[:, roi_mask].mean(axis=1)
        corrs = vcorrcoef(timeseries.T, seed_timeseries)
        corrs[np.isnan(corrs)] = 0
        SBC_dict[seed_name] = (corrs, seed_timeseries)
    return SBC_dict

'''
FC matrix
'''

def parcellated_FC_matrix(sub_timeseries, atlas_idx, roi_list):
    
    timeseries_dict = {}
    for i in roi_list:
        roi_mask = np.asarray(atlas_idx == i, dtype=bool)
        # extract the voxel timeseries within the mask, and take the mean ROI timeseries
        timeseries_dict[str(i)] = sub_timeseries[:,roi_mask].mean(axis=1)

    roi_labels = list(timeseries_dict.keys())
    sub_timeseries = []
    for roi in roi_labels:
        sub_timeseries.append(timeseries_dict[roi])
    corr_matrix = np.corrcoef(sub_timeseries)
    return corr_matrix,roi_labels


def plot_matrix(filename, corr_matrix):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    g = ax.imshow(corr_matrix, cmap='coolwarm', vmax=1, vmin=-1)
    ax.axis('off')
    cbar = plt.colorbar(g, ax=ax, shrink=0.5)
    cbar.set_label('R score', rotation=270, fontsize=10)
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight', dpi=150)


'''
ICA
'''

def run_group_ICA(bold_file_list, mask_file, dim, random_seed, background_image, disableMigp=False):
    import os
    import pandas as pd

    # create a filelist.txt
    file_path = os.path.abspath('filelist.txt')
    from rabies.utils import flatten_list
    merged = flatten_list(list(bold_file_list))
    df = pd.DataFrame(data=merged)
    df.to_csv(file_path, header=False, sep=',', index=False)

    from rabies.utils import run_command
    out_dir = os.path.abspath('group_melodic.ica')
    command = f'melodic -i {file_path} -m {mask_file} -o {out_dir} -d {dim} --report --seed={str(random_seed)} --bgimage={background_image}'
    if disableMigp:
        command+=' --disableMigp'
    rc,c_out = run_command(command)
    IC_file = out_dir+'/melodic_IC.nii.gz'
    return out_dir, IC_file


'''
NPR
'''
def compute_NPR(timeseries, prior_map_vectors, prior_bold_idx, optimize_NPR_dict,
                NPR_temporal_comp, NPR_spatial_comp,
                network_weighting, filename_split, figure_format):
    import os
    from .analysis_math import spatiotemporal_prior_fit

    C_prior = prior_map_vectors[prior_bold_idx,:].T

    if optimize_NPR_dict['apply']:
        modeling,_,_,_,optimize_report_fig = spatiotemporal_fit_converge(timeseries,C_prior,
                                    window_size=optimize_NPR_dict['window_size'],
                                    min_prior_corr=optimize_NPR_dict['min_prior_corr'],
                                    diff_thresh=optimize_NPR_dict['diff_thresh'],
                                    max_iter=optimize_NPR_dict['max_iter'], 
                                    compute_max=optimize_NPR_dict['compute_max'], 
                                    gen_report=True)
        optimize_report_file = os.path.abspath(f'{filename_split}_NPR_optimize.{figure_format}')
        optimize_report_fig.savefig(optimize_report_file, bbox_inches='tight')
    else:
        if NPR_temporal_comp<1: # make sure there is no negative number
            NPR_temporal_comp=0
        if NPR_spatial_comp<1: # make sure there is no negative number
            NPR_spatial_comp=0
        modeling = spatiotemporal_prior_fit(timeseries, C_prior, num_W=NPR_temporal_comp, num_C=NPR_spatial_comp)
        optimize_report_file=None

    # put together the spatial and temporal extra components
    W_extra = np.concatenate((modeling['W_spatial'],modeling['W_temporal']), axis=1)

    if network_weighting=='absolute':
        C_extra = np.concatenate((modeling['C_spatial']*modeling['S_spatial'],
            modeling['C_temporal']*modeling['S_temporal']), axis=1)
        C_fit = modeling['C_fitted_prior']*modeling['S_fitted_prior']
    elif network_weighting=='relative':
        C_extra = np.concatenate((modeling['C_spatial'],
            modeling['C_temporal']), axis=1)
        C_fit = modeling['C_fitted_prior']
    else:
        raise 
    W_fit = modeling['W_fitted_prior']
    return C_fit, W_fit, C_extra, W_extra, optimize_report_file


def eval_convergence(prior_corr_list,fit_diff_list,window_size=5,min_prior_corr=0.5,diff_thresh=0.04):
    if len(fit_diff_list)<window_size:
        return None
    for i in range(1,len(fit_diff_list)-window_size+1): # we start at 1 because idx 0 is nan
        if (np.array(prior_corr_list[i-1])<min_prior_corr).max():
            continue # if the preceding iteration does not pass corr threshold, search for another window
        window = np.array(fit_diff_list[i:i+window_size])
        if not (window>diff_thresh).max(): # if all diff are below threshold, return the first index
            return i-1 # we take -1 because the iteration preceding the window is ideal
    return None


def spatiotemporal_fit_converge(X,C_prior,window_size=5,min_prior_corr=0.5,diff_thresh=0.03,max_iter=20, compute_max=False, gen_report=False):
    '''
    NPR is conducted while incrementally increasing the number of components (alternating 
    between a spatial and a temporal component). 
    
    Convergence criterion 1: Iterations continue until the correlation between the fitted
    component and the prior does not reach the specified minimum.
    
    Convergence Criterion 2: At each iteration, the difference between the previous and new 
    output is evaluated (0=perfectly correlated; 1=uncorrelated). The forming set of 
    successive iterations (within a certain window length) is evaluated, and when a set 
    respects the convergence threshold for each iteration within the window, the iteration 
    preceding that window is selected as optimal output. We take the iteration preceding 
    the  window, as this corresponds to the last iteration which generated changes above 
    threshold. The sliding-window approach is employed to prevent falling within a local 
    minima, when further ameliorations may be possible with further iterations.
    
    When multiple priors are fitted, they are all simultaneously subjected to the evaluation
    of convergence, and as long as one prior fit does not meet the thresholds, iterations
    continue.
    
    X: timeseries of shape time by voxel
    C_prior: set of priors of shape voxel by num priors
    window_size: the size of the sliding window
    min_prior_corr: the minimum correlation with the prior
    diff_thresh: the convergence threshold
    max_iter: maximum number of components to derive
    '''
    
    from rabies.analysis_pkg.analysis_math import spatiotemporal_prior_fit

    num_prior = C_prior.shape[1]
    num_list=range(0,max_iter)

    fit_list=[]
    # first entry is NaN, since there is no comparison to make
    fit_diff_list=[np.repeat(np.nan,num_prior)] 
    NPR_dict_list=[]
    prior_corr_list=[]
    for num in num_list:
        num_C = int(num/2)
        if num%2:
            num_W = num_C+1
        else:
            num_W = num_C

        NPR_dict = spatiotemporal_prior_fit(X, C_prior, num_W=num_W, num_C=num_C)
        NPR_dict_list.append(NPR_dict)
        prior_corr_list.append(NPR_dict['corr_list'])
        C_fit = NPR_dict['C_fitted_prior']
        C_fit /= np.sqrt((C_fit ** 2).sum(axis=0))
        
        if len(fit_list)>0:
            ##### evaluate convergence
            C = C_fit
            C_prev = fit_list[-1]
            # each prior has a separate lim evaluation
            lim = np.abs(np.abs((C * C_prev).sum(axis=0)) - 1)
            fit_diff_list.append(lim)

        convergence_idx = eval_convergence(prior_corr_list,fit_diff_list,window_size=window_size,min_prior_corr=min_prior_corr,diff_thresh=diff_thresh)
        if not compute_max and not convergence_idx is None:
            if gen_report:
                fig = generate_convergence_report(convergence_idx,fit_diff_list,prior_corr_list,min_prior_corr=min_prior_corr,diff_thresh=diff_thresh)
            return NPR_dict_list[convergence_idx],convergence_idx,fit_diff_list,prior_corr_list,fig
        
        fit_list.append(C_fit)

    if gen_report:
        fig = generate_convergence_report(convergence_idx,fit_diff_list,prior_corr_list,min_prior_corr=min_prior_corr,diff_thresh=diff_thresh)
    if convergence_idx is None: # if there was no convergence, the last iteration is outputed
        convergence_idx = -1
    return NPR_dict_list[convergence_idx],None,fit_diff_list,prior_corr_list,fig


def generate_convergence_report(convergence_idx,fit_diff_list,prior_corr_list,min_prior_corr=0.5,diff_thresh=0.03):
    import matplotlib.pyplot as plt

    num_list=range(len(fit_diff_list))
    num_priors = len(prior_corr_list[0])
    
    fig,axes=plt.subplots(1,2,figsize=(10,5))
    
    axes[0].plot(num_list,prior_corr_list)
    axes[1].plot(num_list,fit_diff_list)
    if not convergence_idx is None:
        axes[0].scatter([num_list[convergence_idx]]*num_priors,prior_corr_list[convergence_idx], color='r', marker='*', s=80)
        axes[1].scatter([num_list[convergence_idx]]*num_priors,fit_diff_list[convergence_idx], color='r', marker='*', s=80)
        
    axes[0].set_ylim([0,1])
    axes[0].set_ylabel('Correlation with prior', color='white', fontsize=15)
    axes[1].set_ylim([0,min(diff_thresh*4,1)])
    axes[1].set_ylabel('Difference from previous', color='white', fontsize=15)
    
    axes[0].set_xlabel('Number of non-prior components', color='white', fontsize=15)
    axes[0].set_xlim([-1,max(num_list)+1])
    axes[1].set_xlabel('Number of non-prior components', color='white', fontsize=15)
    axes[1].set_xlim([-1,max(num_list)+1])
    
    axes[0].plot([-1,max(num_list)+1],[min_prior_corr,min_prior_corr], color='lightgray', linestyle='--', linewidth=2)
    axes[1].plot([-1,max(num_list)+1],[diff_thresh,diff_thresh], color='lightgray', linestyle='--', linewidth=2)
    
    legend_labels = [f'Prior {i+1}' for i in range(num_priors)]
    axes[0].legend(legend_labels, fontsize=12, loc='lower right')
    axes[1].legend(legend_labels, fontsize=12, loc='upper right')
    plt.tight_layout()
    return fig