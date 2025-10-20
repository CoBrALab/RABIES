import numpy as np
import matplotlib.pyplot as plt
from ..analysis_math import closed_form


def cosine_similarity(X,Y=None): 
    # estimate cosine similarity across the columns of X and Y combined    
    if Y is None:
        arr = X.copy()
    else:
        arr = np.concatenate((X,Y),axis=1)
    arr /= np.sqrt((arr ** 2).sum(axis=0)) # first normalize from l2-norm
    Sc = (arr.T.dot(arr)) # then derive vector product across columns
    return Sc


def gen_report(X,C_prior, Cs, min_prior_sim=None, Dc_W_thresh=None, Dc_C_thresh=None):
    '''
    Generate CPCA fitting report, and derive an optimal value for n* using minimal 
    prior similarity threshold, and/or thresholds on Dc_Wnet and/or Dc_Cnet.
    If min_prior_sim=None, Dc_t_thresh=None, and Dc_s_thresh=None, then n_optim_idx 
    will be None (no optimal n* is derived).
    '''
    
    # generate the fitting report
    Snet_s_l,Snet_t_l,prior_sim_l,Dc_Cnet,Dc_Wnet,R2_s,R2_t = evaluate_fit(X,C_prior,Cs)

    prior_sim_l = np.array(prior_sim_l)
    
    n_optim_idx = optim_n(prior_sim_l,min_prior_sim,Dc_Cnet,Dc_Wnet,Dc_W_thresh,Dc_C_thresh)

    n_prior = C_prior.shape[1]
    N_max = Cs.shape[1]
    fig_list = plot_report(n_prior,N_max,n_optim_idx,min_prior_sim,Dc_W_thresh,Dc_C_thresh,prior_sim_l,
                    Snet_s_l,Snet_t_l,Dc_Cnet,Dc_Wnet,R2_s,R2_t)
    return n_optim_idx, fig_list
    

def evaluate_fit(X,C_prior,Cs):
    n_prior = C_prior.shape[1]
    N = Cs.shape[1]
    C_prior_ = np.concatenate((C_prior, Cs),axis=1) # expand prior
    W = closed_form(C_prior_, X.T).T

    # Wprior and Ws can be estimated with all Cs once, because adding new Cs do not impact their definition
    W_prior = W[:,:n_prior]
    W_prior /= np.sqrt((W_prior ** 2).mean(axis=0)) # normalization by root mean square so that the scale is estimated on C

    Ws = W[:,n_prior:]

    Cnet_l = []
    Wnet_l = []

    prior_sim_l = []

    SSe_l = []
    SSnet_s_l = []
    SSnet_t_l = []
    
    Snet_s_l = []
    Snet_t_l = []

    for i in range(N+1):
        Wt = Ws[:,:i]
        Ct = closed_form(Wt, X).T
        
        # compute residuals after spatial or temporal orthogonal filter
        X_s = X - Ws[:,:i].dot(Cs[:,:i].T)
        X_t = X - Wt.dot(Ct.T)
        
        # apply Wt correction
        Wnet = W_prior - Wt.dot(closed_form(Wt,W_prior))
        Wnet /= np.sqrt((Wnet ** 2).mean(axis=0)) # normalization by root mean square so that the scale is estimated on C
        Wnet_l.append(Wnet)
                        
        Cnet_s = closed_form(W_prior, X_s).T
        Cnet_l.append(Cnet_s)
        Cnet_t = closed_form(Wnet, X_t).T
        
        Snet_s_l.append(np.sqrt((Cnet_s ** 2).sum(axis=0)))
        Snet_t_l.append(np.sqrt((Cnet_t ** 2).sum(axis=0)))
        
        # residuals are the same between s and t versions
        X_net_s = W_prior.dot(Cnet_s.T)
        e = X_s - X_net_s
        SSe = (e ** 2).sum()
        SSe_l.append(SSe)

        prior_sim = [cosine_similarity(Cnet_s[:,[n]],C_prior[:,[n]])[0,1] for n in range(n_prior)]
        prior_sim_l.append(prior_sim)

        # we can calculate the entire variance of the subspace formed by the matrix product WC from the variance on each vector
        SSnet_s = (W_prior**2).sum(axis=0)*(Cnet_s**2).sum(axis=0)
        SSnet_s_l.append(SSnet_s)
        SSnet_t = (Wnet**2).sum(axis=0)*(Cnet_t**2).sum(axis=0)
        SSnet_t_l.append(SSnet_t)
        

    '''
    Cosine distance between iterations
    '''
    Cnet_arr = np.array(Cnet_l).transpose(1,0,2)
    Wnet_arr = np.array(Wnet_l).transpose(1,0,2)

    # l2 normalization on each component
    Cnet_arr = Cnet_arr/np.sqrt((Cnet_arr ** 2).sum(axis=0))
    Wnet_arr = Wnet_arr/np.sqrt((Wnet_arr ** 2).sum(axis=0))
    # then the vector product derives cosine similarity
    Sc_Cnet = (Cnet_arr[:,1:,:]*Cnet_arr[:,:-1,:]).sum(axis=0)
    Sc_Wnet = (Wnet_arr[:,1:,:]*Wnet_arr[:,:-1,:]).sum(axis=0)

    # we convert to cosine distance (positive valued)
    Dc_Cnet = np.concatenate((np.array([np.repeat(np.nan,n_prior)]), np.abs(np.abs(Sc_Cnet) - 1)),axis=0)
    Dc_Wnet = np.concatenate((np.array([np.repeat(np.nan,n_prior)]), np.abs(np.abs(Sc_Wnet) - 1)),axis=0)

    '''
    Compute variance losses
    '''
    SSe_loss = np.concatenate((np.array([np.nan]), np.array(SSe_l)[:-1]-np.array(SSe_l)[1:]),axis=0)
    SSnet_s_loss = np.concatenate((np.array([np.repeat(np.nan,n_prior)]), np.array(SSnet_s_l)[:-1]-np.array(SSnet_s_l)[1:]),axis=0)
    SSnet_t_loss = np.concatenate((np.array([np.repeat(np.nan,n_prior)]), np.array(SSnet_t_l)[:-1]-np.array(SSnet_t_l)[1:]),axis=0)
    
    R2_s = SSnet_s_loss/(SSnet_s_loss.T+SSe_loss).T
    R2_t = SSnet_t_loss/(SSnet_t_loss.T+SSe_loss).T
    
    return Snet_s_l,Snet_t_l,prior_sim_l,Dc_Cnet,Dc_Wnet,R2_s,R2_t


def optim_n(prior_sim_l,min_prior_sim,Dc_Cnet,Dc_Wnet,Dc_W_thresh,Dc_C_thresh):
    if (min_prior_sim is None) and (Dc_W_thresh is None) and (Dc_C_thresh is None):
        return None
    
    if not min_prior_sim is None:
        # find last idx that did NOT pass threshold
        min_corr_idx_set = np.where(prior_sim_l<min_prior_sim)[0]
        if len(min_corr_idx_set)==0: # threshold passed everywhere
            min_corr_idx = 0
        else:
            min_corr_idx = min_corr_idx_set.max()+1 # we add 1 since we select the one after the last fail identified
            if min_corr_idx>=prior_sim_l.shape[0]:
                print("Minimal prior correlation was not met.")
                return None # there was no iteration above threshold
    else:
        min_corr_idx = 0

    optim_idx_l = [min_corr_idx]
    # taking last idx passing threshold across all networks
    for Dc,Dc_thresh in zip([Dc_Cnet,Dc_Wnet],[Dc_C_thresh,Dc_W_thresh]):
        if Dc_thresh is None:
            continue
        Dc_idx = np.where(Dc>Dc_thresh)[0]
        if len(Dc_idx)==0:
            last_Dc_idx = 0
        else:
            last_Dc_idx = Dc_idx.max()
        optim_idx_l.append(last_Dc_idx)
        
    # select the iteration that pass min_corr and last_Dc_idx thresholds
    optim_idx = np.array(optim_idx_l).max()
    
    return optim_idx


def plot_report(n_prior,N_max,n_optim_idx,min_prior_sim,Dc_W_thresh,Dc_C_thresh,prior_sim_l,
                Snet_s_l,Snet_t_l,Dc_Cnet,Dc_Wnet,R2_s,R2_t):
    
    fig_list = []
    for prior_idx in range(n_prior):
        
        fig,axes = plt.subplots(2,2, figsize=(8,8), constrained_layout=True)

        ax = axes[0,0]
        ax.plot(np.array(Snet_s_l)[:,prior_idx])
        ax.plot(np.array(Snet_t_l)[:,prior_idx])
        ax.legend(['Orthogonal CPCA','Non-orthogonal CPCA'], fontsize=12)
        ax.set_title('Network amplitude', fontsize=15)
        ax.set_ylabel('L2-norm', fontsize=15)
        ax.set_ylim([0,max(max(np.array(Snet_s_l)[:,prior_idx]),max(np.array(Snet_t_l)[:,prior_idx]))*1.05])

        ax = axes[0,1]
        ax.plot(prior_sim_l[:,prior_idx])
        ax.set_ylim([0,1])
        ax.set_title('Prior similarity', fontsize=15)
        ax.set_ylabel('Cosine similarity', fontsize=15)
        if not min_prior_sim is None:
            ax.plot([-1,N_max+1],[min_prior_sim,min_prior_sim], color='lightgray', linestyle='--', linewidth=2)
        if not n_optim_idx is None:
            ax.scatter(n_optim_idx,prior_sim_l[n_optim_idx,prior_idx], color='r', marker='*', s=80)


        ax = axes[1,0]
        ax.plot(Dc_Cnet[:,prior_idx])
        ax.plot(Dc_Wnet[:,prior_idx])
        ax.legend(['Spatial map','Timecourse'], fontsize=12, loc='upper right')
        ax.set_ylim([0,0.4])
        ax.set_ylabel('Cosine distance', fontsize=15)
        ax.set_title('Network component change', fontsize=15)

        thresh_l = []
        if not Dc_W_thresh is None:
            thresh_l.append(Dc_W_thresh)
            ax.plot([-1,N_max+1],[Dc_W_thresh,Dc_W_thresh], color='bisque', linestyle='--', linewidth=2)
        if not Dc_C_thresh is None:
            thresh_l.append(Dc_C_thresh)
            ax.plot([-1,N_max+1],[Dc_C_thresh,Dc_C_thresh], color='lightsteelblue', linestyle='--', linewidth=2)
        if len(thresh_l)>0:
            ax.set_ylim([0,min(max(thresh_l)*4,1)])

        if not n_optim_idx is None:
            ax.scatter(n_optim_idx,Dc_Cnet[n_optim_idx,prior_idx], color='r', marker='*', s=80)

        ax = axes[1,1]
        ax.plot(R2_s[:,prior_idx])
        ax.plot(R2_t[:,prior_idx])
        ax.set_ylabel('Redundancy', fontsize=15)
        ax.set_title('CPCA and network redundancy', fontsize=15)
        ax.legend(['Orthogonal CPCA','Non-orthogonal CPCA'], fontsize=12, loc='upper right')
        ax.set_ylim([0,1])

        for i in range(2):
            axes[1,i].set_xlabel('Number of CPCA components', fontsize=15)
            plt.setp(axes[0,i].get_xticklabels(), visible=False)
        for ax_ in axes:
            for ax in ax_:
                ax.set_xlim([-1,N_max+1])
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
        
        fig_list.append(fig)
    return fig_list
