import numpy as np
from ..analysis_math import closed_form
from .report import gen_report
from sklearn.utils import check_random_state


def spatial_cpca(X, q=1, c_init=None, W_prior=None, tol=1e-6, max_iter=200, verbose=1):
    """
    Extracts spatially orthogonal complementary components from a data matrix using the first step of Complementary Principal Component Analysis (CPCA).
    
    Parameters:
        X (ndarray): Data matrix of shape (time points, voxels).
        q (int, optional): Number of spatial components to extract. Default is 1.
        c_init (ndarray, optional): Initial spatial weights matrix of shape (voxels, components). If None, initialized randomly.
        W_prior (ndarray, optional): Prior temporal weights matrix of shape (time points, components). If None, set to zeros.
        tol (float, optional): Convergence tolerance. Default is 1e-6.
        max_iter (int, optional): Maximum number of iterations. Default is 200.
        verbose (int, optional): Verbosity level for logging progress.
    
    Returns:
        Cs (ndarray): Extracted spatial components (voxels × q).
        Ws (ndarray): Corresponding temporal components (time points × q).
        C_ (ndarray): Combined components matrix after update.
    """
    # X: time by voxel matrix
    # c_init: can specify an voxel by component number matrix for initiating weights
    
    if c_init is None:
        random_state = check_random_state(None)
        c_init = random_state.normal(
            size=(X.shape[1], q))

    if W_prior is None:
        W_prior = np.zeros([X.shape[0], 0])
    W_prior /= np.sqrt((W_prior ** 2).mean(axis=0))
    
    Cs = c_init
    Cs /= np.sqrt((Cs ** 2).mean(axis=0))
    
    for i in range(max_iter):
        Cs_prev = Cs
        
        Ws = closed_form(Cs, X.T).T

        C_ = closed_form(np.concatenate((Ws, W_prior), axis=1), X).T
        Cs = C_[:,:q]

        Cs /= np.sqrt((Cs ** 2).mean(axis=0))
        
        ##### evaluate convergence
        lim = np.abs(np.abs((Cs * Cs_prev).mean(axis=0)) - 1).mean()
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
    
    return Cs,Ws,C_

def cpca(X, C_prior, sequence = ['s','s','s','s'], verbose=False):
    """
    Performs Complementary Principal Component Analysis (CPCA) to extract spatial and temporal components from a data matrix.
    
    This function iteratively derives spatial and temporal components according to a specified sequence, applies temporal correction to ensure orthogonality, normalizes component weights, and computes the final component matrices.
    
    Parameters:
        X (ndarray): Data matrix of shape (time points, voxels).
        C_prior (ndarray): Prior spatial components matrix.
        sequence (list of str, optional): List indicating the order and type of components to extract ('s' for spatial, 't' for temporal). Default is four spatial components.
        verbose (bool, optional): If True, prints progress information.
    
    Returns:
        Cnet (ndarray): Corrected prior spatial components.
        Wnet (ndarray): Corrected prior spatial weights.
        Cs (ndarray): Extracted spatial components.
        Ws (ndarray): Extracted spatial weights.
        Ct (ndarray): Extracted temporal components.
        Wt (ndarray): Extracted temporal weights.
        C (ndarray): Final combined component matrix.
        W (ndarray): Final combined weights matrix.
    """
    
    q=1 # one component is derived at a time to have deterministic outputs
    
    if C_prior is None:
        C_prior = np.zeros([X.shape[1], 0])
    C_prior_ = C_prior.copy() # C_prior_ will be including CPCA components
    n_prior = C_prior.shape[1]

    for i in range(len(sequence)): # derive all WstCst
        if verbose:
            print(i)
        W_prior_ = closed_form(C_prior_, X.T).T
        cs,ws,C_ = spatial_cpca(X, q=q, W_prior=W_prior_, tol=1e-6)
        C_prior_ = np.concatenate((C_prior_, cs),axis=1) # expand prior
        
    W_prior_ = closed_form(C_prior_, X.T).T
    W_prior = W_prior_[:,:n_prior]
    Wst = W_prior_[:,n_prior:]
    
    s_idx = np.where(np.array(sequence)=='s')[0]
    t_idx = np.where(np.array(sequence)=='t')[0]
    n_s = len(s_idx)
    n_t = len(t_idx)
    
    Ws = Wst[:,s_idx]
    Wt = Wst[:,t_idx]
    
    # apply Wt correction
    Wnet = W_prior - Wt.dot(closed_form(Wt,W_prior))
    
    # derive final model on X
    Wnet /= np.sqrt((Wnet ** 2).mean(axis=0)) # normalization, so that the scale is in C
    Ws /= np.sqrt((Ws ** 2).mean(axis=0)) # normalization, so that the scale is in C
    Wt /= np.sqrt((Wt ** 2).mean(axis=0)) # normalization, so that the scale is in C
    W = np.concatenate((Wnet, Ws, Wt),axis=1)
    
    C = closed_form(W, X).T
        
    Cnet = C[:,:n_prior]
    Cs = C[:,n_prior:n_prior+n_s]
    Ct = C[:,n_prior+n_s:]

        
    return Cnet,Wnet,Cs,Ws,Ct,Wt,C,W


def cpca_quick(X, C_prior, sequence = ['s','s','s','s'], tol=1e-10, verbose=False):
    """
    Performs CPCA with accelerated convergence by residualizing the data matrix after each spatial component extraction.
    
    This function iteratively extracts spatial components from the data matrix, subtracting each extracted component's contribution before the next iteration. After all components are derived according to the specified sequence, temporal correction and normalization are applied to the resulting weights and components.
    
    Parameters:
    	sequence (list of str): Specifies the order and type of components to extract ('s' for spatial, 't' for temporal).
    	tol (float): Convergence tolerance for spatial CPCA steps.
    
    Returns:
    	Cnet (ndarray): Corrected prior spatial components.
    	Wnet (ndarray): Corrected prior spatial weights.
    	Cs (ndarray): Extracted spatial components.
    	Ws (ndarray): Extracted spatial weights.
    	Ct (ndarray): Extracted temporal components.
    	Wt (ndarray): Extracted temporal weights.
    	C (ndarray): Final combined component matrix.
    	W (ndarray): Final combined weight matrix.
    """
    
    q=1 # one component is derived at a time to have deterministic outputs
    
    if C_prior is None:
        C_prior = np.zeros([X.shape[1], 0])
    n_prior = C_prior.shape[1]
    W_prior = closed_form(C_prior, X.T).T
    
    Cs = np.zeros([X.shape[1], 0])
    X_ = X.copy()
    for i in range(len(sequence)): # derive all WstCst
        if verbose:
            print(i)
        cs,ws,C_ = spatial_cpca(X_, q=q, W_prior=W_prior, tol=tol)
        ws = closed_form(cs, X_.T).T
        X_ = X_ - ws.dot(cs.T)
        Cs = np.concatenate((Cs, cs),axis=1) # expand prior
                
    C_prior_ = np.concatenate((C_prior, Cs),axis=1) # expand prior
    W_prior_ = closed_form(C_prior_, X.T).T
    W_prior = W_prior_[:,:n_prior]
    Wst = W_prior_[:,n_prior:]
    
    s_idx = np.where(np.array(sequence)=='s')[0]
    t_idx = np.where(np.array(sequence)=='t')[0]
    n_s = len(s_idx)
    n_t = len(t_idx)
    
    Ws = Wst[:,s_idx]
    Wt = Wst[:,t_idx]
    
    # apply Wt correction
    Wnet = W_prior - Wt.dot(closed_form(Wt,W_prior))
    
    # derive final model on X
    Wnet /= np.sqrt((Wnet ** 2).mean(axis=0)) # normalization, so that the scale is in C
    Ws /= np.sqrt((Ws ** 2).mean(axis=0)) # normalization, so that the scale is in C
    Wt /= np.sqrt((Wt ** 2).mean(axis=0)) # normalization, so that the scale is in C
    W = np.concatenate((Wnet, Ws, Wt),axis=1)
    
    C = closed_form(W, X).T
        
    Cnet = C[:,:n_prior]
    Cs = C[:,n_prior:n_prior+n_s]
    Ct = C[:,n_prior+n_s:]

        
    return Cnet,Wnet,Cs,Ws,Ct,Wt,C,W


def cpca_auto(X, C_prior, N_max, Wt_n='n', min_prior_sim=None, Dc_W_thresh=None, Dc_C_thresh=None):
    """
    Performs automated Complementary Principal Component Analysis (CPCA) with estimation of the optimal number of spatial components.
    
    This function runs CPCA with up to `N_max` spatial components, uses a fitting report to determine the optimal number of components, and constructs the final model accordingly. The number of temporal components used can be limited by `Wt_n`. Returns the corrected prior components and weights, spatial and temporal components and weights, final combined matrices, and diagnostic figures.
    
    Parameters:
        X (np.ndarray): Data matrix (time points × voxels).
        C_prior (np.ndarray): Prior spatial components matrix.
        N_max (int): Maximum number of spatial components to consider.
        Wt_n (int or 'n', optional): Maximum number of temporal components to use; if 'n', uses all available.
        min_prior_sim (float, optional): Minimum similarity threshold for prior components in the fitting report.
        Dc_W_thresh (float, optional): Threshold for temporal component diagnostics.
        Dc_C_thresh (float, optional): Threshold for spatial component diagnostics.
    
    Returns:
        Cnet (np.ndarray): Corrected prior spatial components.
        Wnet (np.ndarray): Corrected prior weights.
        Cs (np.ndarray): Derived spatial components.
        Ws (np.ndarray): Derived spatial weights.
        Ct (np.ndarray): Derived temporal components.
        Wt (np.ndarray): Derived temporal weights.
        C (np.ndarray): Final combined components matrix.
        W (np.ndarray): Final combined weights matrix.
        fig_list (list): List of diagnostic figures from the fitting report.
    """
    
    if Wt_n=='n':
        Wt_n = N_max
    elif not isinstance(Wt_n, int):
        raise ValueError("Wt_n must be an integer or 'n'.") 
    elif Wt_n<0:
        raise ValueError("Wt_n must be minimum 0.") 
        
    n_prior = C_prior.shape[1]
    sequence = ['s']*N_max
    Cnet,Wnet,Cs,Ws,Ct,Wt,C,W = cpca_quick(X, C_prior, sequence, verbose=False)
    
    # generate the fitting report
    n_optim_idx, fig_list = gen_report(X, C_prior, Cs, min_prior_sim=min_prior_sim, Dc_W_thresh=Dc_W_thresh, Dc_C_thresh=Dc_C_thresh)
    
    if n_optim_idx is None:
        # if no threshold was applied, or min_prior_sim was not met, we still generate a model using N_max
        n_optim_idx = N_max
    
    # derive model with optimized n*
    C_prior_ = np.concatenate((C_prior, Cs[:,:n_optim_idx]),axis=1)
    W_prior_ = closed_form(C_prior_, X.T).T
    W_prior = W_prior_[:,:n_prior]
    Wst = W_prior_[:,n_prior:]
    
    Wt = Wst[:,:Wt_n]
    Ws = Wst[:,Wt_n:]
    n_s = Ws.shape[1]
    
    # apply Wt correction
    Wnet = W_prior - Wt.dot(closed_form(Wt,W_prior))
    
    # derive final model on X
    Wnet /= np.sqrt((Wnet ** 2).mean(axis=0)) # normalization, so that the scale is in C
    Ws /= np.sqrt((Ws ** 2).mean(axis=0)) # normalization, so that the scale is in C
    Wt /= np.sqrt((Wt ** 2).mean(axis=0)) # normalization, so that the scale is in C
    W = np.concatenate((Wnet, Ws, Wt),axis=1)
    
    C = closed_form(W, X).T
        
    Cnet = C[:,:n_prior]
    Cs = C[:,n_prior:n_prior+n_s]
    Ct = C[:,n_prior+n_s:]

    return Cnet,Wnet,Cs,Ws,Ct,Wt,C,W,fig_list