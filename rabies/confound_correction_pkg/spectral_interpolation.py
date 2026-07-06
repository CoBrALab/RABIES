'''
Adaptation of the Lomb-Scargle periodogram for spectral interpolation of irregularly sampled time series,
as described in Mathias et al. 2004, and applied in Power et al. 2014

Mathias, A., Grond, F., Guardans, R., Seese, D., Canela, M., & Diebner, H. H. (2004). 
Algorithms for spectral analysis of irregularly sampled time series. Journal of Statistical Software, 11(1), 1-27.

Power, J. D., Mitra, A., Laumann, T. O., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2014).
Methods to detect, characterize, and remove motion artifact in resting state fMRI. Neuroimage, 84, 320-341.
'''

import numpy as np

def lombscargle_mathias(t, x, w):
    """Fit amplitude coefficients for a set of frequencies on irregularly
    sampled data using the Lomb-Scargle rotation.

    For each angular frequency in w, finds the time shift theta that makes
    cos(w*(t - theta)) and sin(w*(t - theta)) orthogonal over the sample
    times t, then computes the projection of x onto those two basis
    vectors. This is mathematically equivalent to the classical Lomb-Scargle
    periodogram but returns amplitude coefficients (cw, sw) rather than
    collapsing to a scalar power estimate — retaining the information needed
    to reconstruct the signal at arbitrary times.

    Parameters
    ----------
    t : array, shape (n,)     — irregular sample times in seconds
    x : array, shape (n,) or (n, n_voxels) — signal values at those times
    w : array, shape (n_freqs,) — angular frequencies

    Returns
    -------
    cw    : cosine amplitudes, shape (n_freqs,) or (n_freqs, n_voxels)
    sw    : sine amplitudes,   shape (n_freqs,) or (n_freqs, n_voxels)
    theta : phase shifts,      shape (n_freqs,)

    Docstring co-authored with Claude (Anthropic).
    """

    if (w == 0).sum()>0:
        raise ZeroDivisionError()

    # Check input sizes
    if t.shape[0] != x.shape[0]:
        raise ValueError("Input arrays do not have the same size.")

    # Create empty array for output periodogram
    cw = np.empty((len(w)))
    sw = np.empty((len(w)))
    theta = np.empty((len(w)))

    w_t = w[:,np.newaxis].dot(t[np.newaxis,:])
    theta = (1 / (2 * w)) * np.arctan2(
        np.sin(2 * w_t).sum(axis=1),
        np.cos(2 * w_t).sum(axis=1))

        
    wt = (t[:,np.newaxis]-theta[np.newaxis,:])*w
    c=np.cos(wt)
    s=np.sin(wt)
    
    cw = x.T.dot(c)/(c**2).sum(axis=0)
    sw = x.T.dot(s)/(s**2).sum(axis=0)
    
    return cw,sw,theta


def lombscargle_mathias_simulate(t, w, cw, sw,theta):
    # recover simulated timeseries for a given time vector t

    wt = (t[:,np.newaxis]-theta[np.newaxis,:])*w
    c=np.cos(wt)
    s=np.sin(wt)
    
    y = c.dot(cw.T)+s.dot(sw.T)
    return y


def spectral_interpolation(x, time_step, time_mask, rescale=True):
    """Fill censored timepoints by spectral interpolation using the Lomb-Scargle
    rotation method.

    Builds a frequency grid of all DFT harmonics up to the Nyquist frequency,
    fits cosine/sine amplitudes at each frequency using only the uncensored
    samples (via lombscargle_mathias), then evaluates the spectral model at
    every timepoint.  Censored positions in the output are replaced by the
    model prediction; uncensored positions retain their original values.

    Parameters
    ----------
    x : array_like, shape (n_good, n_voxels)
        Signal values at the uncensored timepoints only.  Must have the same
        number of rows as the number of True entries in time_mask.
    time_step : float
        Sampling interval in seconds (TR).
    time_mask : array_like of bool, shape (n_timepoints,)
        True = uncensored, False = censored / to be filled.
    rescale : bool, optional
        If True (default), rescale the spectral reconstruction so that its
        standard deviation matches that of x.
        This corrects for the mild amplitude inflation that arises when
        independent per-frequency fits are summed on gapped (non-uniformly
        sampled) data — a consequence of cross-frequency leakage breaking
        the orthogonality of the DFT basis.

    Returns
    -------
    x_fill : ndarray, shape (n_timepoints, n_voxels)
        Gap-filled timeseries; uncensored timepoints are identical to x.
    x_sim  : ndarray, shape (n_timepoints, n_voxels)
        Full spectral reconstruction evaluated at every timepoint (before
        gap-only replacement), useful for diagnosing fit quality.

    Docstring co-authored with Claude (Anthropic).
    """
    num_timepoints = len(time_mask)
    time = np.linspace(time_step, num_timepoints * time_step, num_timepoints)

    t_good = time[time_mask]
    t_span = t_good.max() - t_good.min() + time_step  # total duration: each sample covers one TR-length bin

    w = _build_frequency_grid_harmonic(t_span, time_step)
    cw, sw, theta = lombscargle_mathias(t_good, x, w)
    x_sim = lombscargle_mathias_simulate(time, w, cw, sw, theta)

    if rescale:
        # re-scale to match the original amplitude of the good points, 
        # since there is some oversampling in the frequency domain because 
        # of the time gaps
        x_sim -= x_sim[time_mask].mean(axis=0)
        x_sim *= x.std(axis=0) / x_sim[time_mask].std(axis=0)
        x_sim += x.mean(axis=0)


    x_fill = np.zeros(x_sim.shape)
    x_fill[time_mask, :] = x
    x_fill[~time_mask, :] = x_sim[~time_mask]
    return x_fill, x_sim


def _build_frequency_grid_harmonic(t_span, time_step):
    """
    Building a linear frequency grid of all DFT harmonics from lowest
    frequency of 2π/t_span (domega) to the Nyquist frequency of 
    π/time_step (omegamax)

    Every frequency is an exact integer multiple of domega = 2*pi/t_span,
    the Fourier resolution of the available data.  This is the classical
    orthogonal basis for spectral estimation: on complete uniformly-sampled
    data these sinusoids are exactly orthogonal; on gapped data they are as
    orthogonal as the data permits, minimising inter-frequency leakage when
    summing independent per-frequency fits.

    Parameters
    ----------
    t_span : float
        Time span of the good (uncensored) samples in seconds; defines
        the Fourier resolution domega = 2*pi/t_span.
    time_step : float
        Time step in seconds.

    Returns
    -------
    omegas : ndarray, all integer multiples of domega in [omegamin, omegamax].

    Docstring co-authored with Claude (Anthropic).
    """
    omegamax = np.pi / time_step # Nyquist frequency
    domega = 2.0 * np.pi / t_span # minimum frequency resolution
    k_max = int(omegamax / domega) # maximal number of harmonics
    return np.arange(1, int(k_max) + 1) * (domega)