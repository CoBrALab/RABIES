from sklearn.utils import check_random_state
import numpy as np
from rabies.analysis_pkg.analysis_math import closed_form

# functions that computes the Least Squares Estimates


def _logcosh(x, fun_args=None):
    alpha = fun_args.get('alpha', 1.0)  # comment it out?

    x *= alpha
    x = np.tanh(x)  # apply the tanh inplace
    gx = x

    g_x = np.empty(x.shape[0])
    # XXX compute in chunks to avoid extra allocation
    for i, gx_i in enumerate(gx):  # please don't vectorize.
        g_x[i] = (alpha * (1 - gx_i ** 2)).mean()
    return gx, g_x


def deflation_fit(X, q=1, c_init=None, C_convergence='OLS', C_prior=None, W_prior=None, W_ortho=False, tol=1e-6, max_iter=200, verbose=1):
    # q defines the number of new sources to fit
    if c_init is None:
        random_state = check_random_state(None)
        c_init = random_state.normal(
            size=(X.shape[1], q))

    # the C_prior and W_prior correspond to spatial and temporal priors respectively which will impose an orthogonality contraint
    # on the fitted sources in their respective dimension
    if C_prior is None:
        C_prior = np.zeros([X.shape[1], 0])
    C_prior /= np.sqrt((C_prior ** 2).sum(axis=0))

    if W_prior is None:
        W_prior = np.zeros([X.shape[0], 0])
    W_prior /= np.sqrt((W_prior ** 2).sum(axis=0))

    # initialize an empty C
    C = np.zeros([X.shape[1], 0])
    for j in range(q):
        C_prev = np.concatenate((C, C_prior), axis=1)
        c = c_init[:, j].reshape(-1, 1)
        c /= np.sqrt((c ** 2).sum(axis=0))

        # regress out the orthogonal dimensions already fitted
        X_ = X-np.matmul(np.matmul(X, C_prev), C_prev.T)

        for i in range(max_iter):
            c_prev = c

            w = np.matmul(X_, c)
            if W_ortho and (W_prior.shape[1]>0):
                # impose complete temporal orthogonality, more similar to CR before but not quite the same
                w -= np.matmul(W_prior, closed_form(W_prior, w))
            W = np.concatenate((w,W_prior),axis=1) # include the W priors in the convergence step

            if C_convergence == 'OLS':
                c = closed_form(W, X_).T[:,0].reshape(-1, 1) # take back only c
            elif C_convergence == 'ICA':
                gwtx, g_wtx = _logcosh(W[:,0].reshape(1,-1), {})
                c = ((X_.T * gwtx).mean(axis=1) - g_wtx.mean() * c.T).T
            else:
                raise

            # impose spatial orthogonality
            c -= np.matmul(np.matmul(c.T, C_prev), C_prev.T).T
            c /= np.sqrt((c ** 2).sum(axis=0))

            ##### evaluate convergence
            lim = np.abs(np.abs((c * c_prev).sum(axis=0)) - 1).mean()
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
        C = np.concatenate((C, c), axis=1)
    return C

def dual_ICA_fit(timeseries, num_comp, all_IC_vectors, prior_bold_idx):
    prior_fit_out={'C':[],'W':[]}
    convergence_function = 'ICA'
    X=timeseries

    prior_networks = all_IC_vectors[prior_bold_idx,:].T

    C_prior=prior_networks
    C_conf = deflation_fit(X, q=num_comp, c_init=None, C_convergence=convergence_function,
                      C_prior=C_prior, W_prior=None, W_ortho=True, tol=1e-6, max_iter=200, verbose=1)
    for network in range(prior_networks.shape[1]):
        prior=prior_networks[:,network].reshape(-1,1)
        C_prior=np.concatenate((prior_networks[:,:network],prior_networks[:,network+1:],C_conf),axis=1)

        C_fit = deflation_fit(X, q=1, c_init=prior, C_convergence=convergence_function,
                              C_prior=C_prior, W_prior=None, W_ortho=True, tol=1e-6, max_iter=200, verbose=1)

        # make sure the sign of weights is the same as the prior
        corr = np.corrcoef(C_fit.flatten(), prior.flatten())[0, 1]
        if corr < 0:
            C_fit = C_fit*-1

        # the finalized C
        C = np.concatenate((C_fit, C_prior), axis=1)

        # L-2 norm normalization of the components
        C /= np.sqrt((C ** 2).sum(axis=0))
        W = closed_form(C,X.T, intercept=False).T
        # the components will contain the weighting/STD/singular value, and the timecourses are normalized
        C=C*W.std(axis=0)
        # normalize the component timecourses to unit variance
        W /= W.std(axis=0)

        prior_fit_out['C'].append(C[:,0])
        prior_fit_out['W'].append(W[:,0])
    return prior_fit_out
