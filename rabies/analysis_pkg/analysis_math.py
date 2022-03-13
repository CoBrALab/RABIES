import numpy as np

def vcorrcoef(X, y):  # return a correlation between each row of X with y
    Xm = np.reshape(np.mean(X, axis=1), (X.shape[0], 1))
    ym = np.mean(y)
    r_num = np.sum((X-Xm)*(y-ym), axis=1)
    r_den = np.sqrt(np.sum((X-Xm)**2, axis=1)*np.sum((y-ym)**2))
    r = r_num/r_den
    return r


def elementwise_corrcoef(X, Y):
    # X and Y are each of shape num_observations X num_element
    # computes the correlation between each element of X and Y
    Xm = X.mean(axis=0)
    Ym = Y.mean(axis=0)
    r_num = np.sum((X-Xm)*(Y-Ym), axis=0)

    r_den = np.sqrt(np.sum((X-Xm)**2, axis=0)*np.sum((Y-Ym)**2, axis=0))
    r = r_num/r_den
    return r


def elementwise_spearman(X,Y):    
    order = X.argsort(axis=0)
    X_ranks = order.argsort(axis=0)
    order = Y.argsort(axis=0)
    Y_ranks = order.argsort(axis=0)
    return elementwise_corrcoef(X_ranks, Y_ranks)
    

def dice_coefficient(mask1,mask2):
    dice = np.sum(mask1*mask2)*2.0 / (np.sum(mask1) + np.sum(mask2))
    return dice


'''
LINEAR REGRESSION --- CLOSED-FORM SOLUTION
'''


def closed_form(X, Y, intercept=False):  # functions that computes the Least Squares Estimates
    if intercept:
        X = np.concatenate((X, np.ones([X.shape[0], 1])), axis=1)
    return np.linalg.inv(X.transpose().dot(X)).dot(X.transpose()).dot(Y)


def mse(X, Y, w):  # function that computes the Mean Square Error (MSE)
    return np.mean((Y-np.matmul(X, w))**2)


'''
Convergence through alternating minimization using OLS
'''

from sklearn.utils import check_random_state

def dual_OLS_fit(X, q=1, c_init=None, C_prior=None, tol=1e-6, max_iter=200, verbose=1):
    # X: time by voxel matrix
    # q: number of new components to fit
    # c_init: can specify an voxel by component number matrix for initiating weights
    # C_prior: a voxel by component matrix of priors that are included in the fitting, but fixed as constant components
    
    if q<1:
        return np.zeros([X.shape[1], 0])
    
    # q defines the number of new sources to fit
    if c_init is None:
        random_state = check_random_state(None)
        c_init = random_state.normal(
            size=(X.shape[1], q))
        
    if C_prior is None:
        C_prior = np.zeros([X.shape[1], 0])
    C_prior /= np.sqrt((C_prior ** 2).sum(axis=0))
    num_prior = C_prior.shape[1]
    
    C = c_init
    C /= np.sqrt((C ** 2).sum(axis=0))
    
    for i in range(max_iter):
        C_prev = C
        C_ = np.concatenate((C, C_prior), axis=1) # add in the prior to contribute to the fitting
        
        ##### temporal convergence step
        W = closed_form(C_, X.T).T

        C_ = closed_form(W, X).T
            
        if num_prior>0:
            C = C_[:,:-num_prior] # take out the fitted components
        else:
            C = C_
        C /= np.sqrt((C ** 2).sum(axis=0))

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
    return C


'''

def closed_form_3d(X,Y):
    return np.matmul(np.matmul(np.linalg.inv(np.matmul(X.transpose(0,2,1),X)),X.transpose(0,2,1)),Y)

def lme_stats_3d(X,Y):
    #add an intercept
    X=np.concatenate((X,np.ones((X.shape[0],X.shape[1],1))),axis=2)
    [num_comparisons,num_observations,num_predictors] = X.shape
    [num_comparisons,num_observations,num_features] = Y.shape

    w=closed_form_3d(X,Y)

    residuals = Y-np.matmul(X, w)
    MSE = (((residuals)**2).sum(axis=1)/(num_observations-num_predictors))


    var_b = np.expand_dims(MSE, axis=1)*np.expand_dims(np.linalg.inv(np.matmul(X.transpose(0,2,1),X)).diagonal(axis1=1,axis2=2), axis=2)
    sd_b = np.sqrt(var_b) # standard error on the Betas
    ts_b = w/sd_b # calculate t-values for the Betas
    p_values =[2*(1-stats.t.cdf(np.abs(ts_b[:,i,:]),(num_observations-num_predictors))) for i in range(ts_b.shape[1])] # calculate a p-value map for each predictor

    return ts_b,p_values,w,residuals

non_nan_idx = (np.isnan(voxelwise_array).sum(axis=(0,1))==0)

# take out the voxels which have null values
X=voxelwise_array[:,:6,non_nan_idx].transpose(2,0,1)
Y=voxelwise_array[:,6:,non_nan_idx].transpose(2,0,1)


ts_b,p_values,w,residuals = lme_stats_3d(X,Y)

x_name=['Somatomotor','Dorsal Comp','DMN', 'Prior Modeling 1', 'Prior Modeling 2', 'Prior Modeling 3']
y_name=['group','temporal_std','VE_spatial','GS_corr','DVARS_corr','FD_corr']

fig,axes = plt.subplots(nrows=len(x_name), ncols=len(y_name),figsize=(12*len(y_name),3*len(x_name)))


for i,x_label in zip(range(len(x_name)),x_name):
    for j,y_label in zip(range(len(y_name)),y_name):
        ax=axes[i,j]

        stat_map=np.zeros(voxelwise_array.shape[2])
        stat_map[non_nan_idx]=ts_b[:,j,i]

        ax.set_title('T-value of {} on {}'.format(y_label,x_label), fontsize=15)
        plot_stat_map(analysis_functions.recover_3D(mask_file,stat_map),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y')
'''
