'''
Constrained PCA (CPCA)
'''

import numpy as np
import matplotlib.pyplot as plt
from .analysis_math import closed_form
from sklearn.utils import check_random_state


def decorr(X,Y):
    # decorrelated Y relative to X
    return Y-X.dot(closed_form(X,Y))


def dual_OLS_single_fit(X, c_init=None, C_prior=None, tol=1e-6, max_iter=200, verbose=1):
    '''
    Derives one component orthogonal relative to C_prior, following dual_OLS convergence including C_prior in the regression
    '''

    # X: time by voxel matrix
    # c_init: can specify an voxel by component number matrix for initiating weights
    
    q=1
    # q defines the number of new sources to fit
    if c_init is None:
        random_state = check_random_state(None)
        c_init = random_state.normal(
            size=(X.shape[1], q))

    if C_prior is None:
        C_prior = np.zeros([X.shape[1], 0])
    C_prior /= np.sqrt((C_prior ** 2).sum(axis=0))
    
    C = c_init
    C /= np.sqrt((C ** 2).sum(axis=0))
    
    for i in range(max_iter):
        C_prev = C
        C_ = np.concatenate((C, C_prior), axis=1) # add in the prior to contribute to the fitting
        
        ##### first OLS convergence step
        W_ = closed_form(C_, X.T).T
        
        ##### second OLS convergence step
        C_ = closed_form(W_, X).T        
        C_ /= np.sqrt((C_ ** 2).sum(axis=0))

        C = C_[:,:q] # take out the fitted components
        W = W_[:,:q] # take out the fitted components
        
        ##### evaluate convergence
        lim = np.abs(np.abs((C * C_prev).sum(axis=0)) - 1).mean()
        if verbose > 2:
            print('lim:'+str(lim))
        if lim < tol:
            if verbose > 1:
                print(str(i)+' iterations to converge.')
            break
        if i == max_iter-1:
            if verbose > 0:
                print(
                    'Convergence failed. Consider increasing max_iter or decreasing tol.')
                
    
    return C, W, C_, W_


def complete_Wprior(X, w_init=None, C_prior=None, W_extra=None, tol=1e-6, max_iter=200, verbose=1):
    
    if C_prior is None:
        C_prior = np.zeros([X.shape[1], 0])
    C_prior /= np.sqrt((C_prior ** 2).sum(axis=0))
    
    num_priors = C_prior.shape[1]
    
    if w_init is None:
        random_state = check_random_state(None)
        w_init = random_state.normal(
            size=(X.shape[0], num_priors))
    
    
    if W_extra is None:
        W_extra = np.zeros([X.shape[1], 0])
    W_extra /= np.sqrt((W_extra ** 2).sum(axis=0))
    
    W = w_init
    W /= np.sqrt((W ** 2).sum(axis=0))
    
    for i in range(max_iter):
        W_prev = W
        
        W_ = np.concatenate((W, W_extra), axis=1)
                
        C_ = closed_form(W_, X).T
        C_ /= np.sqrt((C_ ** 2).sum(axis=0))
        
        C = C_[:,num_priors:]
        C_ = np.concatenate((C_prior, C), axis=1)
        
        W_ = closed_form(C_, X.T).T
        
        # impose W orthogonality
        W_[:,:num_priors] = decorr(W_extra,W[:,:num_priors])
        W_ /= np.sqrt((W_ ** 2).sum(axis=0))
       
        W = W_[:,:num_priors]        
 
        ##### evaluate convergence
        lim = np.abs(np.abs((W * W_prev).sum(axis=0)) - 1).mean()
        if verbose > 2:
            print('lim:'+str(lim))
        if lim < tol:
            if verbose > 1:
                print(str(i)+' iterations to converge.')
            break
        if i == max_iter-1:
            if verbose > 0:
                print(
                    'Convergence failed. Consider increasing max_iter or decreasing tol.')
    return W, C, W_, C_


def spatiotemporal_CPCA(X, C_prior, num_W, num_C, optim_dim=False, min_prior_corr=0.4, diff_thresh_t=0.05, diff_thresh_s=0.03):
    C_prior = C_prior.copy()

    Ct_extra = np.zeros([X.shape[1], 0])
    Cs_extra = np.zeros([X.shape[1], 0])

    C_prior_t = C_prior
    
    # derive num_W temporally orthogonal components
    for i in range(num_W):
        # deriving the orthogonal component W_extra
        C, W_extra, C_, W_ = dual_OLS_single_fit(X, c_init=None, C_prior=C_prior_t, tol=1e-6, max_iter=200, verbose=1)
        # uncovering the proper C_extra counterpart
        W, C_extra_t, W_, C_ = complete_Wprior(X, w_init=None, C_prior=C_prior_t, W_extra=W_extra, tol=1e-10, max_iter=200, verbose=1)
        # add new component to set of Ct_extra
        Ct_extra = np.concatenate((Ct_extra,C_extra_t), axis=1)
        # expend the C_prior set before next iteration
        C_prior_t = np.concatenate((C_prior_t,C_extra_t), axis=1)
        
    # the number of temporal component is optimized before fitting spatial components
    if optim_dim:
        optim_idx_t,prior_corr_list_t,CPCA_fit_diff_spacetime_t = find_optim_dim(X,C_prior,Ct_extra,Cs_extra,min_prior_corr,diff_thresh_t)
        # if convergence failed, then no optimal dimension is selected
        if not optim_idx_t is None:
            Ct_extra = Ct_extra[:,:optim_idx_t]

    Cs_extra = np.zeros([X.shape[1], 0])
    C_prior_s = np.concatenate((C_prior,Ct_extra), axis=1)
    # derive num_C temporally orthogonal components
    for j in range(num_C):
        # fit 1 new spatial component
        C_extra_s, W, C_, W_ = dual_OLS_single_fit(X, c_init=None, C_prior=C_prior_s, tol=1e-6, max_iter=200, verbose=1)
        # add new component to set of Cs_extra
        Cs_extra = np.concatenate((Cs_extra,C_extra_s), axis=1)
        C_prior_s = np.concatenate((C_prior_s,C_extra_s), axis=1)
        
    if optim_dim:
        optim_idx_s,prior_corr_list_s,CPCA_fit_diff_spacetime_s = find_optim_dim(X,C_prior,Cs_extra,Ct_extra,min_prior_corr,diff_thresh_s)
        # if convergence failed, then no optimal dimension is selected
        if not optim_idx_s is None:
            Cs_extra = Cs_extra[:,:optim_idx_s]
    
        fig,axes=plt.subplots(2,2,figsize=(10,10))
        if num_W>0:
            axes[0,0].set_title('Temporal CPCA optimization', color='white', fontsize=20)
            axes[0,1].set_title('Temporal CPCA optimization', color='white', fontsize=20)
            plot_convergence_report(axes[0,0],axes[0,1],optim_idx_t,CPCA_fit_diff_spacetime_t,prior_corr_list_t,min_prior_corr,diff_thresh_t)
        if num_C>0:
            axes[1,0].set_title('Spatial CPCA optimization', color='white', fontsize=20)
            axes[1,1].set_title('Spatial CPCA optimization', color='white', fontsize=20)
            plot_convergence_report(axes[1,0],axes[1,1],optim_idx_s,CPCA_fit_diff_spacetime_s,prior_corr_list_s,min_prior_corr,diff_thresh_s)
        plt.tight_layout()
    else:
        fig = None


    '''
    Compute dual regression to obtain a fit for the prior combined with the CPCA components
    '''

    # estimate a final set for C and W, offering a final linear fit to X = WtC
    C = np.concatenate((C_prior, Cs_extra, Ct_extra), axis=1)
    C /= np.sqrt((C ** 2).mean(axis=0))
    W = closed_form(C, X.T).T
    W /= np.sqrt((W ** 2).mean(axis=0)) # the temporal domain is variance-normalized so that the weights are contained in the spatial maps
    C = closed_form(W, X).T
    S = np.sqrt((C ** 2).mean(axis=0)) # the component variance/scaling is taken from the spatial maps
    C /= S # the spatial maps are variance normalized; the variance is stored in S

    num_priors = C_prior.shape[1]

    # Fitted priors are at the first indices
    C_fitted_prior = C[:,:num_priors]
    W_fitted_prior = W[:,:num_priors]
    S_fitted_prior = S[:num_priors]
    
    # temporal components are at the last indices
    C_temporal = C[:,num_priors+num_C:]
    W_temporal = W[:,num_priors+num_C:]
    S_temporal = S[num_priors+num_C:]

    # spatial components are in the middle
    C_spatial = C[:,num_priors:num_priors+num_C]
    W_spatial = W[:,num_priors:num_priors+num_C]
    S_spatial = S[num_priors:num_priors+num_C]
    
    # we thus output a model of the timeseries of the form X = W.dot((S*C).T)
    modeling = {'C_fitted_prior':C_fitted_prior, 'C_spatial':C_spatial, 'C_temporal':C_temporal, 
            'W_fitted_prior':W_fitted_prior, 'W_spatial':W_spatial, 'W_temporal':W_temporal, 
            'S_fitted_prior':S_fitted_prior, 'S_spatial':S_spatial, 'S_temporal':S_temporal, 
            }

    return modeling, Ct_extra, Cs_extra, fig


'''
Dimensionality optimization
'''

def find_optim_dim(X,C_prior,C_extra_optim,C_extra_fixed,min_prior_corr,diff_thresh):
    num_priors = C_prior.shape[1]
    
    C_fit_list=[]
    W_fit_list=[]
    for i in range(C_extra_optim.shape[1]+1): # up to +1, because first index is 0 components
        C = np.concatenate((C_prior,C_extra_optim[:,:i],C_extra_fixed), axis=1)
        W = closed_form(C, X.T).T
        W /= np.sqrt((W ** 2).mean(axis=0))
        C = closed_form(W, X).T
        C /= np.sqrt((C ** 2).mean(axis=0))
        C_fitted_prior = C[:,:num_priors]
        W_fitted_prior = W[:,:num_priors]

        C_fit_list.append(C_fitted_prior)
        W_fit_list.append(W_fitted_prior)

    FC_maps = np.array(C_fit_list).transpose(1,2,0)
    FC_times = np.array(W_fit_list).transpose(1,2,0)

    # estimate spatial correlation to prior
    from .analysis_math import elementwise_corrcoef
    prior_corr_list = np.array([elementwise_corrcoef(C_prior, FC_maps[:,:,i]) for i in range(FC_maps.shape[2])])

    # find last idx that did NOT pass threshold
    min_corr_idx_set = np.where(prior_corr_list<min_prior_corr)[0]
    if len(min_corr_idx_set)==0: # threshold passed everywhere
        min_corr_idx = 0
    else:
        min_corr_idx = min_corr_idx_set.max()+1 # we add 1 since we select the one after the last fail identified
        if min_corr_idx>=FC_maps.shape[2]:
            print("Minimal prior correlation was not met.")
            min_corr_idx = None # there was no iteration above threshold

    # estimate diff between subsequent iterations
    maps_dot = (FC_maps[:,:,1:] * FC_maps[:,:,:-1]).mean(axis=0)
    times_dot = (FC_times[:,:,1:] * FC_times[:,:,:-1]).mean(axis=0)
    # we concatenate an empty idx to have congruent dimension with prior_corr_list
    CPCA_fit_diff_spacetime = np.concatenate((np.array([np.repeat(np.nan,num_priors)]), np.abs(np.abs(maps_dot*times_dot) - 1).T),axis=0)

    # taking last idx passing threshold across all networks
    diff_idx = np.where(CPCA_fit_diff_spacetime>diff_thresh)[0]
    if len(diff_idx)==0:
        last_diff_idx = 0
    else:
        last_diff_idx = diff_idx.max()

    if min_corr_idx is None:
        optim_idx = None
    else:
        # select the iteration that pass min_corr and last_diff thresholds
        optim_idx = np.array([min_corr_idx,last_diff_idx]).max()
    
    return optim_idx,prior_corr_list,CPCA_fit_diff_spacetime
    

def plot_convergence_report(ax1,ax2,convergence_idx,fit_diff_list,prior_corr_list,min_prior_corr,diff_thresh):
    
    num_list=range(len(fit_diff_list))
    num_priors = len(prior_corr_list[0])
    ax1.plot(num_list,prior_corr_list)
    ax2.plot(num_list,fit_diff_list)
    if not convergence_idx is None:
        ax1.scatter([num_list[convergence_idx]]*num_priors,prior_corr_list[convergence_idx], color='r', marker='*', s=80)
        ax2.scatter([num_list[convergence_idx]]*num_priors,fit_diff_list[convergence_idx], color='r', marker='*', s=80)

    ax1.set_ylim([0,1])
    ax1.set_ylabel('Correlation with prior', color='white', fontsize=15)
    ax2.set_ylim([0,min(diff_thresh*4,1)])
    ax2.set_ylabel('Difference from previous', color='white', fontsize=15)

    ax1.set_xlabel('Number of CPCA components', color='white', fontsize=15)
    ax1.set_xlim([-1,max(num_list)+1])
    ax2.set_xlabel('Number of CPCA components', color='white', fontsize=15)
    ax2.set_xlim([-1,max(num_list)+1])

    ax1.plot([-1,max(num_list)+1],[min_prior_corr,min_prior_corr], color='lightgray', linestyle='--', linewidth=2)
    ax2.plot([-1,max(num_list)+1],[diff_thresh,diff_thresh], color='lightgray', linestyle='--', linewidth=2)

    legend_labels = [f'Prior {i+1}' for i in range(num_priors)]
    ax1.legend(legend_labels, fontsize=12, loc='lower right')
    ax2.legend(legend_labels, fontsize=12, loc='upper right')

