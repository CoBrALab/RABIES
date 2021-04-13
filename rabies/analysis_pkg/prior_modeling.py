from sklearn.utils import check_random_state
import numpy as np
import torch
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# functions that computes the Least Squares Estimates


def torch_closed_form(X, Y):
    return torch.matmul(torch.matmul(torch.inverse(torch.matmul(X.T, X)), X.T), Y)


def torch_dual_regression(C, X):
    # C is of shape num_voxelsxnum_ICs
    # X is of shape num_timepointsxnum_voxels

    # for one given volume, it's values can be expressed through a linear combination of the components
    # LR1: X = (CWt)t; or Xt = CWt
    W = torch_closed_form(C, X.T).T

    # for a given voxel timeseries, it's signal can be explained as a linear combination of the component timecourses
    # LR2: X = WCt
    C = torch_closed_form(W, X).T

    # L-2 norm normalization of the components
    std = torch.sqrt((C ** 2).sum(axis=0))
    C /= std
    # preserve the scaling of variance in W
    W *= std
    return C, W


'''
FastICA functions
'''


def pca_project(X, n_components):
    # X should be of shape num_samples X num_features
    # X will be project onto a low dimensional space through PCA
    # can consult this as reference https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
    # projection = XV = US, where X=USV⊤ from SVD
    U, s, V = torch.svd_lowrank(X, q=n_components, niter=20)
    # returns projection of shape n_samples X n_components, V corresponds to the principal components/directions
    return torch.matmul(X, V).float(), V.float(), s


def _sym_decorrelation(W):
    """ Symmetric decorrelation
    i.e. W <- (W * W.T) ^{-1/2} * W
    """
    s, u = torch.eig(torch.matmul(W, W.T), eigenvectors=True)
    # re-order the eigenvalues/vectors in increasing order of eigen value
    s, idx = torch.sort(s[:, 0])
    u = u[:, idx]
    # u (resp. s) contains the eigenvectors (resp. square roots of
    # the eigenvalues) of W * W.T
    return torch.matmul(torch.matmul(u * (1. / torch.sqrt(s)), u.T), W)


def _logcosh(x, fun_args=None):
    alpha = fun_args.get('alpha', 1.0)  # comment it out?

    x *= alpha
    x = torch.tanh(x)  # apply the tanh inplace
    gx = x

    g_x = torch.empty(x.shape[0], dtype=torch.float32).to(device)
    # XXX compute in chunks to avoid extra allocation
    for i, gx_i in enumerate(gx):  # please don't vectorize.
        g_x[i] = (alpha * (1 - gx_i ** 2)).mean()
    return gx, g_x


'''
Alternating minimization algorithm
'''


def parallel_fit(X, q=5, C_prior=None, W_prior=None, C_ortho=False, W_ortho=False, tol=1e-6, max_iter=200, verbose=1):
    # q defines the number of new sources to fit
    # the C_ortho and W_ortho options define whether to impose a spatial or temporal orthogonality constraint among the fitted sources
    random_state = check_random_state(None)
    c = torch.tensor(random_state.normal(
        size=(X.shape[1], q)), dtype=torch.float32).to(device)
    c /= torch.sqrt((c ** 2).sum(axis=0))
    w = torch_closed_form(c, X.T).T
    w /= torch.sqrt((w ** 2).sum(axis=0))

    # the C_prior and W_prior correspond to spatial and temporal priors respectively which will impose an orthogonality contraint
    # on the fitted sources in their respective dimension
    if C_prior == None:
        C_prior = torch.zeros(X.shape[1], 0).float().to(device)
    C_prior = _sym_decorrelation(C_prior.T).T  # impose orthogonality
    C_prior /= torch.sqrt((C_prior ** 2).sum(axis=0))

    if W_prior == None:
        W_prior = torch.zeros(X.shape[0], 0).float().to(device)
    W_prior = _sym_decorrelation(W_prior.T).T  # impose orthogonality
    W_prior /= torch.sqrt((W_prior ** 2).sum(axis=0))
    # take a seperate C equivalent for the temporal priors
    C_t = torch.pinverse(torch_closed_form(W_prior, X).T).T

    # only include fitted sources to be spatially constrained, add back spatial prior
    C = torch.cat((c, C_prior, C_t), axis=1)
    C /= torch.sqrt((C ** 2).sum(axis=0))

    for i in range(max_iter):
        c_prev = c
        w_prev = w

        ##### temporal convergence step
        W = torch_closed_form(C, X.T).T

        # take out the spatial priors which don't get temporally constrained
        W_s = W[:, q:q+C_prior.shape[1]].clone()
        # only include fitted sources to be temporally constrained, add back temporal prior
        W = torch.cat((W[:, :q], W_prior), axis=1)
        if W_ortho:  # impose across the entire temporal domain
            W = _sym_decorrelation(W.T).T
        else:  # introduce a partial decorrelation across time
            for j in range(q):
                W[:, j] -= torch.matmul(torch.matmul(W[:, j].T,
                                                     W_prior), W_prior.T).T
        # take out only the components that are being fitted
        w = W[:, :q].clone()
        w /= torch.sqrt((w ** 2).sum(axis=0))

        # append the temporal version of the spatial priors
        W = torch.cat((W, W_s), axis=1)

        ##### spatial convergence step
        # take the inverse because C is not orthogonal after OLS
        C = torch.pinverse(torch_closed_form(W, X).T).T

        # take out the temporal priors which don't get spatially constrained
        C_t = C[:, q:q+W_prior.shape[1]].clone()
        # only include fitted sources to be spatially constrained, add back spatial prior
        C = torch.cat((C[:, :q], C_prior), axis=1)
        if C_ortho:  # impose across the entire spatial domain
            C = _sym_decorrelation(C.T).T
        else:  # introduce a partial decorrelation across space
            for j in range(q):
                C[:, j] -= torch.matmul(torch.matmul(C[:, j].T,
                                                     C_prior), C_prior.T).T
        # append the spatial version of the temporal priors
        C = torch.cat((C, C_t), axis=1)
        # normalize over spatial domain to prevent unstability
        C /= torch.sqrt((C ** 2).sum(axis=0))
        c = C[:, :q].clone()

        ##### evaluate convergence
        w_lim = torch.abs(torch.abs((w * w_prev).sum(axis=0)) - 1).mean()
        c_lim = torch.abs(torch.abs((c * c_prev).sum(axis=0)) - 1).mean()
        lim = np.sqrt(c_lim.cpu().numpy()**2+w_lim.cpu().numpy()**2)
        if verbose > 2:
            print(' w:'+str(w_lim)+' c:'+str(c_lim)+' lim:'+str(lim))
        if lim < tol:
            if verbose > 1:
                print(str(i)+' iterations to converge.')
            break
        if i == max_iter-1:
            if verbose > 0:
                print(
                    'Convergence failed. Consider increasing max_iter or decreasing tol.')

    return C, W


def deflation_fit(X, q=1, c_init=None, C_convergence='OLS', C_prior=None, W_prior=None, W_ortho=False, tol=1e-6, max_iter=200, verbose=1):
    # q defines the number of new sources to fit
    if c_init is None:
        random_state = check_random_state(None)
        c_init = torch.tensor(random_state.normal(
            size=(X.shape[1], q)), dtype=torch.float32).to(device)

    # the C_prior and W_prior correspond to spatial and temporal priors respectively which will impose an orthogonality contraint
    # on the fitted sources in their respective dimension
    if C_prior is None:
        C_prior = torch.zeros(X.shape[1], 0).float().to(device)
    C_prior /= torch.sqrt((C_prior ** 2).sum(axis=0))

    if W_prior is None:
        W_prior = torch.zeros(X.shape[0], 0).float().to(device)
    W_prior /= torch.sqrt((W_prior ** 2).sum(axis=0))

    # initialize an empty C
    C = torch.zeros(X.shape[1], 0).float().to(device)
    for j in range(q):
        C_prev = torch.cat((C, C_prior), axis=1)
        c = c_init[:, j].clone().reshape(-1, 1)
        c /= torch.sqrt((c ** 2).sum(axis=0))

        # regress out the orthogonal dimensions already fitted
        X_ = X-torch.matmul(torch.matmul(X, C_prev), C_prev.T)

        for i in range(max_iter):
            c_prev = c

            w = torch.matmul(X_, c)
            if W_ortho and (W_prior.shape[1]>0):
                # impose complete temporal orthogonality, more similar to CR before but not quite the same
                w -= torch.matmul(W_prior, torch_closed_form(W_prior, w))
            W = torch.cat((w,W_prior),axis=1) # include the W priors in the convergence step

            if C_convergence == 'OLS':
                c = torch_closed_form(W, X_).T[:,0].clone().reshape(-1, 1) # take back only c
            elif C_convergence == 'ICA':
                gwtx, g_wtx = _logcosh(W[:,0].reshape(1,-1), {})
                c = ((X_.T * gwtx).mean(axis=1) - g_wtx.mean() * c.T).T
            else:
                raise

            # impose spatial orthogonality
            c -= torch.matmul(torch.matmul(c.T, C_prev), C_prev.T).T
            c /= torch.sqrt((c ** 2).sum(axis=0))

            ##### evaluate convergence
            lim = torch.abs(torch.abs((c * c_prev).sum(axis=0)) - 1).mean()
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
        C = torch.cat((C, c), axis=1)
    return C
