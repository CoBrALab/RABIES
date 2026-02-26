import os
import numpy as np
import SimpleITK as sitk
from rabies.analysis_pkg.analysis_math import closed_form


def tree_list(dirName):
    # Get the list of all files in directory tree at given path
    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(dirName):
        listOfFiles += [os.path.join(dirpath, file) for file in filenames]
    return listOfFiles


def get_info_list(file_list):
    info_list = []
    for file in file_list:
        basename = os.path.basename(file)
        file_info = basename.split(
            '_run-')[0]+'_run-'+basename.split('_run-')[1][0]
        info_list.append(file_info)

    return info_list


def find_scans(scan_info, bold_files, brain_mask_files, confounds_files, csf_mask_files, FD_files):
    for file in bold_files:
        if scan_info in file:
            bold_file = file
            break
    for file in brain_mask_files:
        if scan_info in file:
            brain_mask_file = file
            break
    for file in confounds_files:
        if scan_info in file:
            confounds_file = file
            break
    for file in csf_mask_files:
        if scan_info in file:
            csf_mask = file
            break
    for file in FD_files:
        if scan_info in file:
            FD_file = file
            break
    return bold_file, brain_mask_file, confounds_file, csf_mask, FD_file


def exec_ICA_AROMA(inFile, mc_file, brain_mask, csf_mask, aroma_out=None, tr='auto', aroma_dim=0, random_seed=1):
    import os
    from rabies.confound_correction_pkg.utils import csv2par
    from rabies.confound_correction_pkg.mod_ICA_AROMA.ICA_AROMA_functions import run_ICA_AROMA
    import pathlib
    filename_split = pathlib.Path(inFile).name.rsplit(".nii")
    if not aroma_out:
        aroma_out = os.getcwd()+'/aroma_out'
    cleaned_file = aroma_out+f'/{filename_split[0]}_aroma.nii.gz'

    if tr=='auto':
        from rabies.utils import get_sitk_header
        header = get_sitk_header(os.path.abspath(inFile))
        tr = float(header.GetSpacing()[3])
    else:
        tr = float(tr)

    success, denoising_applied = run_ICA_AROMA(aroma_out, os.path.abspath(inFile), mc=csv2par(mc_file), TR=float(tr), mask=os.path.abspath(
        brain_mask), mask_csf=os.path.abspath(csf_mask), denType="nonaggr", melDir="", dim=str(aroma_dim), overwrite=True, random_seed=random_seed)
    if not success:
        return None, aroma_out

    # if no denoising was applied, return the input file
    if not denoising_applied:
        return inFile, aroma_out
    else:
        os.rename(aroma_out+'/denoised_func_data_nonaggr.nii.gz', cleaned_file)
        return cleaned_file, aroma_out


def csv2par(in_confounds):
    import pandas as pd
    df = pd.read_csv(in_confounds)
    new_df = pd.DataFrame(
        columns=['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3'])
    new_df['mov1'] = df['mov1']
    new_df['mov2'] = df['mov2']
    new_df['mov3'] = df['mov3']
    new_df['rot1'] = df['rot1']
    new_df['rot2'] = df['rot2']
    new_df['rot3'] = df['rot3']
    out_confounds = os.path.abspath(
        (os.path.basename(in_confounds).split('.')[0])+('.par'))
    new_df.to_csv(out_confounds, sep='\t', index=False, header=False)
    return out_confounds


def gen_FD_mask(FD_trace, scrubbing_threshold):
    '''
    Scrubbing based on FD: The frames that exceed the given threshold together with 1 back
    and 2 forward frames will be masked out from the data (as in Power et al. 2012)
    '''
    import numpy as np
    cutoff = np.asarray(FD_trace) >= scrubbing_threshold
    mask = np.ones(len(FD_trace)).astype(bool)
    for i in range(len(mask)):
        if cutoff[i]:
            mask[i-1:i+2] = 0
    return mask


def prep_timeseries_interval(timeseries_interval, num_frames):
    if timeseries_interval == 'all':
        raise ValueError(f"'all' is depreciated as an input for --timeseries_interval. Consult rabies confound_correction --help for new syntax. ")
    # select the subset of timeseries specified
    split = timeseries_interval.split(',')
    if not len(split)==2:
        raise ValueError(f"--timeseries_interval wasn't split into 2: {split}")
    begin=int(split[0]) # must be an integer
    end=split[1]
    if end=='end':
        end=num_frames # take the total number of time frames
    else:
        end=int(end)

    time_range = range(begin,end)
    return time_range


def get_DVARS(timeseries):
    # compute the DVARS before denoising
    # the first data point is set to 0
    derivative=np.concatenate((np.zeros((1, timeseries.shape[1])),timeseries[1:,:]-timeseries[:-1,:]))
    DVARS_trace=np.sqrt((derivative**2).mean(axis=1))
    return DVARS_trace


def temporal_censoring(FD_trace, 
        FD_censoring, FD_threshold, DVARS_trace, DVARS_censoring, minimum_timepoint):

    num_frames = len(FD_trace) # FD_trace length must correspond with the timeseries length
    frame_mask = np.ones(num_frames).astype(bool)
    if FD_censoring:
        FD_mask = gen_FD_mask(FD_trace, FD_threshold)
        frame_mask*=FD_mask
    if DVARS_censoring:
        # create a distribution where no timepoint falls more than 2.5 STD away from the mean
        mask1=np.zeros(len(DVARS_trace)).astype(bool)
        mask2=np.ones(len(DVARS_trace)).astype(bool)
        mask2[0]=False # remove the first timepoint, which is always 0
        while ((mask2!=mask1).sum()>0):
            mask1=mask2
            mean=DVARS_trace[mask1].mean()
            std=DVARS_trace[mask1].std()
            norm=(DVARS_trace-mean)/std
            mask2=np.abs(norm)<2.5
        DVARS_mask=mask2
        frame_mask*=DVARS_mask
    if frame_mask.sum()<int(minimum_timepoint):
        from nipype import logging
        log = logging.getLogger('nipype.workflow')
        log.warning(f"FD/DVARS CENSORING LEFT LESS THAN {str(minimum_timepoint)} VOLUMES. THIS SCAN WILL BE REMOVED FROM FURTHER PROCESSING.")
        return None

    return frame_mask


def select_motion_regressors(nuisance_regressors,motion_params_csv):
    import pandas as pd
    if ('mot_6' in nuisance_regressors) and ('mot_24' in nuisance_regressors):
        raise ValueError(
            "Can't select both the mot_6 and mot_24 options; must pick one.")

    if isinstance(motion_params_csv, pd.DataFrame):
        motion_params_df = motion_params_csv
    else:
        motion_params_df = pd.read_csv(motion_params_csv)
    keys = motion_params_df.keys()
    conf_keys = []
    # only inputing motion regressors at this stage
    for conf in nuisance_regressors:
        if conf == 'mot_6':
            conf_keys += ['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3']
        elif conf == 'mot_24':
            conf_keys += [s for s in keys if "rot" in s or "mov" in s]

    return np.asarray(motion_params_df[conf_keys])


def compute_signal_regressors(timeseries, nuisance_regressors, brain_mask_idx, WM_mask_idx, CSF_mask_idx, vascular_mask_idx):
    from nipype import logging
    log = logging.getLogger('nipype.workflow')

    # make sure there are no NaN voxels
    timeseries[np.isnan(timeseries)] = 0

    regressors_array = np.empty([timeseries.shape[0],0])
    for conf,mask_idx in zip(['WM_signal','CSF_signal','vascular_signal','global_signal'],
                                [WM_mask_idx,CSF_mask_idx,vascular_mask_idx,brain_mask_idx]):
        if conf in nuisance_regressors:
            if mask_idx is None:
                raise ValueError(f"Received a None value for mask_idx computing {conf}.")
            elif mask_idx.sum()==0:
                log.warning(f"0 voxels were found for the {conf} mask. No regressors are computed using this mask.")
            else:
                regressor_trace = timeseries.T[mask_idx].mean(axis=0)
                regressors_array = np.append(regressors_array,regressor_trace.reshape(-1,1),axis=1)
    
    if ('aCompCor_5' in nuisance_regressors) or ('aCompCor_percent' in nuisance_regressors):
        if ('aCompCor_5' in nuisance_regressors) and ('aCompCor_percent' in nuisance_regressors):
            raise ValueError(
                "Can't select both the aCompCor_5 and aCompCor_percent options; must pick one.")
        if 'aCompCor_5' in nuisance_regressors:
            method='aCompCor_5'
        elif 'aCompCor_percent' in nuisance_regressors:
            method='aCompCor_percent'
        else:
            raise
        if WM_mask_idx is None or CSF_mask_idx is None:
            raise ValueError(f"WM_mask or CSF_mask are None - cannot compute aCompCor.")

        combined_mask_idx = (WM_mask_idx+CSF_mask_idx) > 0
        if combined_mask_idx.sum()<5:
            log.warning(f"Less than 5 voxels were found within the combined WM and CSF masks - no aCompCor regressors will be computed.")
        else:
            masked_timeseries = timeseries[:,combined_mask_idx]

            from sklearn.decomposition import PCA
            if method == 'aCompCor_percent':
                pca = PCA()
                comp_timeseries = pca.fit_transform(masked_timeseries)
                explained_variance = pca.explained_variance_ratio_
                cum_var = 0
                num_comp = 0
                # evaluate the # of components to explain 50% of the variance
                while(cum_var <= 0.5):
                    cum_var += explained_variance[num_comp]
                    num_comp += 1
                log.info("Extracting "+str(num_comp)+" components for aCompCorr.")
                comp_timeseries = comp_timeseries[:,:num_comp]

            elif method == 'aCompCor_5':
                num_comp = 5
                pca = PCA(n_components=num_comp)
                comp_timeseries = pca.fit_transform(masked_timeseries)

            # double check that a singular matrix is not outputed
            A = comp_timeseries
            is_singular = np.linalg.matrix_rank(A) < min(A.shape)
            if not is_singular:
                regressors_array = np.append(regressors_array,comp_timeseries,axis=1)
            else:
                log.warning(f"aCompCor regressors yielded a singular matrix - these regressors are not included.")

    if regressors_array.shape[1]>0:
        # double check that a singular matrix is not outputed
        A = regressors_array
        is_singular = np.linalg.matrix_rank(A) < min(A.shape)
        if not is_singular:
            return regressors_array
        else:
            log.warning(f"Nuisance regressors yielded a singular matrix - an empty matrix is outputed instead.")
            return np.empty([timeseries.shape[0],0])
    else:
        return regressors_array


def remove_trend(timeseries, frame_mask, order=1 , time_interval='all'):
    '''The timeseries is already censored, we need to create a timeseries that is censored with the same time mask so that it matches'''
    #count number of non-censored timepoints
    num_timepoints = len(frame_mask)
    
    # we create a centered time axis, to fit an intercept value in the center - Centering improves numerical stability when fitting polynomials (because otherwise high powers of large time indices can explode)
    time=np.linspace(-num_timepoints/2,num_timepoints/2 , num_timepoints)
    time_postcensor = time[frame_mask]
    num_timepoints_postcensor = len(time_postcensor)
    
    # create the design matrix of regressors to remove, essentially just polynomials of the time array, up to the user-specified order (time, time**2, etc)
    X = []
    for k in range(1, order+1):
        X.append(time_postcensor**k) #append time array to list (list of arrays)
    X = np.column_stack(X) if X else np.empty((num_timepoints_postcensor, 0)) # turn list into 2D array; if X is an empty list, will evaluate as false, and the right side ('else') will run, otherwise, the left side, the stack will run

    # Always add intercept as last column (if X is empty, the output will only be a column of intercepts)
    X = np.column_stack([X, np.ones(num_timepoints_postcensor)])

    #fitting 
    Y=timeseries
    # in most cases, we want to fit the trend to the whole timeseries
    if time_interval == 'all':
        W = closed_form(X,Y) #solves the linear regression in closed form
    else:
        #but sometimes we want to compute the trend over only a portion of the timeseries
        lowcut = int(time_interval.split('-')[0])
        highcut = int(time_interval.split('-')[1])
        W = closed_form(X[lowcut:highcut, :],Y[lowcut:highcut, :])
    predicted = X.dot(W) #always subtract the trend over the whole timeseries
    
    residuals = (Y-predicted) # add back the intercept after
    fitted_intercept = W[-1,:] # predicted intercept

    return residuals, fitted_intercept
    

def lombscargle_mathias(t, x, w):

    '''
    Implementation of Lomb-Scargle periodogram as described in Mathias et al. 2004, 
    and applied in Power et al. 2014
    
    Mathias, A., Grond, F., Guardans, R., Seese, D., Canela, M., & Diebner, H. H. (2004). 
    Algorithms for spectral analysis of irregularly sampled time series. Journal of Statistical Software, 11(1), 1-27.
    
    Power, J. D., Mitra, A., Laumann, T. O., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2014). 
    Methods to detect, characterize, and remove motion artifact in resting state fMRI. Neuroimage, 84, 320-341.
    
    '''

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
    theta = (1/(2*w))*np.arctan( \
        np.sin(2*w_t).sum(axis=1)/ \
        np.cos(2*w_t).sum(axis=1))

        
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


def lombscargle_fill(x,time_step,time_mask):
    num_timepoints=len(time_mask)
    time=np.linspace(time_step,num_timepoints*time_step,num_timepoints)

    low_freq = 0.005
    high_freq = 1
    freqs=np.linspace(low_freq,high_freq,1000)
    w = freqs*2*np.pi

    t = time[time_mask]
    cw,sw,theta = lombscargle_mathias(t, x, w)
    
    # simulate over entire time
    t=time
    y = lombscargle_mathias_simulate(t, w, cw, sw,theta)
    
    # standardize according to masked data poaxis=0ints    
    y -= y[time_mask].mean(axis=0)
    y /= y[time_mask].std(axis=0)
    
    # re-scale according to original mean/std
    y *= x.std(axis=0)
    y += x.mean(axis=0)
    
    y_fill = np.zeros(y.shape)
    y_fill[time_mask,:] = x
    y_fill[time_mask==0,:] = y[(time_mask==0)]
    return y_fill


def butterworth(signals, TR, high_pass, low_pass):
    from scipy import signal
    
    critical_freq = []
    if high_pass is not None:
        btype = 'high'
        critical_freq.append(high_pass)
    
    if low_pass is not None:
        btype = 'low'
        critical_freq.append(low_pass)
    
    if len(critical_freq) == 2:
        btype = 'band'
    else:
        critical_freq = critical_freq[0]
    
    order=3
    sos = signal.butter(order, critical_freq, fs=1/TR, btype=btype, output='sos')
    return signal.sosfiltfilt(sos, signals, axis=0)


######### Taken from https://stackoverflow.com/questions/39543002/returning-a-real-valued-phase-scrambled-timeseries
def phaseScrambleTS(ts):
    from scipy.fftpack import fft, ifft
    """Returns a TS: original TS power is preserved; TS phase is shuffled."""
    fs = fft(ts, axis=0)
    pow_fs = np.abs(fs) ** 2.
    phase_fs = np.angle(fs)
    phase_fsr = phase_fs.copy()
    if len(ts) % 2 == 0:
        phase_fsr_lh = phase_fsr[1:int(len(phase_fsr)/2)]
    else:
        phase_fsr_lh = phase_fsr[1:int(len(phase_fsr)/2) + 1]
    np.random.shuffle(phase_fsr_lh)
    if len(ts) % 2 == 0:
        phase_fsr_rh = -phase_fsr_lh[::-1]
        phase_fsr = np.concatenate((np.array((phase_fsr[0],)), phase_fsr_lh,
                                    np.array((phase_fsr[int(len(phase_fsr)/2)],)),
                                    phase_fsr_rh))
    else:
        phase_fsr_rh = -phase_fsr_lh[::-1]
        phase_fsr = np.concatenate((np.array((phase_fsr[0],)), phase_fsr_lh, phase_fsr_rh))
    fsrp = np.sqrt(pow_fs) * (np.cos(phase_fsr) + 1j * np.sin(phase_fsr))
    tsrp = ifft(fsrp, axis=0)
    if not np.allclose(tsrp.imag, np.zeros(tsrp.shape)):
        max_imag = (np.abs(tsrp.imag)).max()
        imag_str = f'\nNOTE: a non-negligible imaginary component was discarded.\n\tMax: {max_imag}'
        print(imag_str)
    return tsrp.real


def phase_randomized_regressors(confounds_array, frame_mask, TR):
    """
    METHOD FROM BRIGHT AND MURPHY 2015
    "The true noise regressors were phase-randomised to create simulated noise regressors with similar frequency 
    distributions to the true noise regressors, but with no relation to the measured head motion or physiology. 
    To achieve this, the frequency spectra of the true noise regressors were obtained using a Fourier transform, 
    and the phase of each frequency in half of the spectrum was randomised and mirrored before performing the 
    inverse Fourier transform (this phase randomization was repeated until the temporal correlation with the true 
    regressor was r b 0.1). The entire set of resulting time-series was then orthogonalised to the complete set 
    of original regressors to make them independent from the true noise."

    """
    from rabies.confound_correction_pkg.utils import lombscargle_fill
    num_conf = confounds_array.shape[1]
    randomized_confounds_array = np.zeros(confounds_array.shape)
    for n in range(num_conf):
        corr=1
        iter=1
        while(corr>0.1):
            x=confounds_array[:,n:n+1]
            # fill missing datapoints to obtain a good reading of frequency spectrum
            y = lombscargle_fill(x,TR,frame_mask)
            # phase randomize
            y_r = phaseScrambleTS(y)
            # re-apply the time mask to the same number of timepoints
            y_m = y_r[frame_mask]
            corr = np.abs(np.corrcoef(x.T,y_m.T)[0,1])
            if iter>100: # set a maximum number of iterations
                from nipype import logging
                log = logging.getLogger('nipype.workflow')
                log.warning("Could not set uncorrelated random regressors!")
                break
            iter += 1
            
        #### impose orthogonality relative to every original regressor      
        y_m -= np.matmul(confounds_array, closed_form(confounds_array, y_m))
        randomized_confounds_array[:,n:n+1] = y_m
    return randomized_confounds_array


def smooth_image(img, fwhm, mask_img):
    # apply nilearn's Gaussian smoothing on a SITK image
    from nilearn.image.image import _smooth_array
    from rabies.utils import copyInfo_4DImage, copyInfo_3DImage

    dim = img.GetDimension()
    if not (dim==3 or dim==4):
        raise ValueError("Image must be 3D or 4D")

    # Build the affine matrix from SITK image
    affine = sitk_affine_lps(img)
    affine = affine[:3, :3]

    if dim==4:
        array_4d = sitk.GetArrayFromImage(img)
        arr = array_4d.transpose(3,2,1,0) # re-orient the array to match the affine
    elif dim==3:
        array_3d = sitk.GetArrayFromImage(img)
        arr = array_3d.transpose(2,1,0) # re-orient the array to match the affine
    smoothed_arr = _smooth_array(arr, affine, fwhm=fwhm, ensure_finite=True, copy=False)
    del arr
    
    # recover SITK orientation
    if dim==4:
        smoothed_arr = smoothed_arr.transpose(3,2,1,0)
    elif dim==3:
        smoothed_arr = smoothed_arr.transpose(2,1,0)
    
    # smoothing creates leakage around mask boundaries
    # correct for edge effects by dividing by the smoothed mask like FSL https://johnmuschelli.com/fslr/reference/fslsmooth.html
    mask_arr = sitk.GetArrayFromImage(mask_img).transpose(2,1,0) # re-orient the array to match the affine
    smoothed_mask = _smooth_array(mask_arr, affine, fwhm=fwhm, ensure_finite=True, copy=False)
    smoothed_mask = smoothed_mask.transpose(2,1,0) # recover SITK orientation
    mask_arr = mask_arr.transpose(2,1,0) # recover SITK orientation
    
    smoothed_arr /= smoothed_mask # correct for edge effect
    smoothed_arr *= mask_arr # re-apply mask to avoid leakage from smoothing
    smoothed_arr[np.isnan(smoothed_arr)] = 0
    
    # recover SITK img
    if dim==4:
        smoothed_img = copyInfo_4DImage(sitk.GetImageFromArray(
            smoothed_arr, isVector=False), img, img)
    elif dim==3:
        smoothed_img = copyInfo_3DImage(sitk.GetImageFromArray(
            smoothed_arr, isVector=False), img)

    return smoothed_img


def sitk_affine_lps(img):
    '''
    Generate the affine matrix from a SITK image.
    Note that SITK affines follow the LPS coordinate.
    '''
    dim = img.GetDimension()
    direction = np.array(img.GetDirection()).reshape(dim, dim)
    spacing = np.diag(img.GetSpacing())
    origin = np.array(img.GetOrigin())

    affine = np.eye(dim + 1)
    affine[:dim, :dim] = direction @ spacing
    affine[:dim, dim] = origin
    return affine