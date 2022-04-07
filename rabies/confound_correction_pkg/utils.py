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


def exec_ICA_AROMA(inFile, mc_file, brain_mask, csf_mask, tr, aroma_dim, random_seed=1):
    import os
    from rabies.confound_correction_pkg.utils import csv2par
    from rabies.confound_correction_pkg.mod_ICA_AROMA.ICA_AROMA_functions import run_ICA_AROMA
    import pathlib
    filename_split = pathlib.Path(inFile).name.rsplit(".nii")
    aroma_out = os.getcwd()+'/aroma_out'
    cleaned_file = aroma_out+f'/{filename_split[0]}_aroma.nii.gz'

    if tr=='auto':
        import SimpleITK as sitk
        img = sitk.ReadImage(os.path.abspath(inFile))
        tr = float(img.GetSpacing()[3])
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


def prep_CR(bold_file, confounds_file, FD_file, cr_opts):
    import pandas as pd
    import SimpleITK as sitk
    from rabies.confound_correction_pkg.utils import select_confound_timecourses

    # save specifically the 6 rigid parameters for AROMA
    confounds_6rigid_array = select_confound_timecourses(['mot_6'],confounds_file,FD_file)

    if len(cr_opts.conf_list)==0:
        confounds_array = confounds_6rigid_array
    else:
        confounds_array = select_confound_timecourses(cr_opts.conf_list,confounds_file,FD_file)

    FD_trace = pd.read_csv(FD_file).get('Mean')

    # select the subset of timeseries specified
    if not cr_opts.timeseries_interval == 'all':
        lowcut = int(cr_opts.timeseries_interval.split(',')[0])
        highcut = int(cr_opts.timeseries_interval.split(',')[1])
        confounds_array = confounds_array[lowcut:highcut, :]
        confounds_6rigid_array = confounds_6rigid_array[lowcut:highcut, :]
        FD_trace = FD_trace[lowcut:highcut]
        time_range = range(lowcut,highcut)
    else:
        time_range = range(sitk.ReadImage(bold_file).GetSize()[3])

    data_dict = {'FD_trace':FD_trace, 'confounds_array':confounds_array, 'confounds_6rigid_array':confounds_6rigid_array, 'confounds_csv':confounds_file, 'time_range':time_range}
    return data_dict

def temporal_censoring(timeseries, FD_trace, 
        FD_censoring, FD_threshold, DVARS_censoring, minimum_timepoint):

    # compute the DVARS before denoising
    derivative=np.concatenate((np.empty([1,timeseries.shape[1]]),timeseries[1:,:]-timeseries[:-1,:]))
    DVARS=np.sqrt((derivative**2).mean(axis=1))

    # apply the temporal censoring
    frame_mask = np.ones(timeseries.shape[0]).astype(bool)
    if FD_censoring:
        FD_mask = gen_FD_mask(FD_trace, FD_threshold)
        frame_mask*=FD_mask
    if DVARS_censoring:
        # create a distribution where no timepoint falls more than 2.5 STD away from the mean
        trace=DVARS
        mask1=np.zeros(len(trace)).astype(bool)
        mask2=np.ones(len(trace)).astype(bool)
        mask2[0]=False # remove the first timepoint, which is always 0
        while ((mask2!=mask1).sum()>0):
            mask1=mask2
            mean=trace[mask1].mean()
            std=trace[mask1].std()
            norm=(trace-mean)/std
            mask2=np.abs(norm)<2.5
        DVARS_mask=mask2
        frame_mask*=DVARS_mask
    if frame_mask.sum()<int(minimum_timepoint):
        from nipype import logging
        log = logging.getLogger('nipype.workflow')
        log.warning(f"FD/DVARS CENSORING LEFT LESS THAN {str(minimum_timepoint)} VOLUMES. THIS SCAN WILL BE REMOVED FROM FURTHER PROCESSING.")
        return None,None,None

    return frame_mask,FD_trace,DVARS


def select_confound_timecourses(conf_list,confounds_file,FD_file):
    import pandas as pd
    if ('mot_6' in conf_list) and ('mot_24' in conf_list):
        raise ValueError(
            "Can't select both the mot_6 and mot_24 options; must pick one.")

    confounds = pd.read_csv(confounds_file)
    keys = confounds.keys()
    conf_keys = []
    for conf in conf_list:
        if conf == 'mot_6':
            conf_keys += ['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3']
        elif conf == 'mot_24':
            conf_keys += [s for s in keys if "rot" in s or "mov" in s]
        elif conf == 'aCompCor':
            aCompCor_keys = [s for s in keys if "aCompCor" in s]
            from nipype import logging
            log = logging.getLogger('nipype.workflow')
            log.info('Applying aCompCor with '+len(aCompCor_keys)+' components.')
            conf_keys += aCompCor_keys
        elif conf == 'mean_FD':
            mean_FD = pd.read_csv(FD_file).get('Mean')
            confounds['mean_FD'] = mean_FD
            conf_keys += [conf]
        else:
            conf_keys += [conf]

    return np.asarray(confounds[conf_keys])


def remove_trend(timeseries, frame_mask, second_order=False, keep_intercept=False):
    num_timepoints = len(frame_mask)
    
    # we create a centered time axis, to fit an intercept value in the center
    time=np.linspace(-num_timepoints/2,num_timepoints/2 , num_timepoints)
    
    # evaluate the linear trend using an interpolated time axis based on previous censoring
    X = time[frame_mask].reshape(-1,1)
    if second_order:
        X = np.concatenate((X, X**2), axis=1) # second order polynomial
    X = np.concatenate((X, np.ones([X.shape[0], 1])), axis=1) # add an intercept at the end
    Y=timeseries
    W = closed_form(X,Y)
    
    predicted = X.dot(W)
    #plt.plot(time[frame_mask], timeseries.mean(axis=1))
    #plt.plot(time[frame_mask], predicted.mean(axis=1))
    
    res = (Y-predicted) # add back the intercept after
    if keep_intercept:
        fitted_intercept = X[:,-1:].dot(W[-1:,:]) # predicted intercept
        res += fitted_intercept
    #plt.plot(time[frame_mask], res.mean(axis=1))
    return res
    

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


def smooth_image(img, affine, fwhm):
    # apply nilearn's Gaussian smoothing on a SITK image
    from nilearn.image.image import _smooth_array
    from rabies.utils import copyInfo_4DImage, copyInfo_3DImage

    # the affine is a 3 by 3 matrix of the spacing*direction for the 3 spatial dimensions
    #spacing_3d = np.array(timeseries_img.GetSpacing())[:3]
    #direction_4d = np.array(timeseries_img.GetDirection())
    #direction_3d = np.array([list(direction_4d[:3]),list(direction_4d[4:7]),list(direction_4d[8:11])])
    #affine = direction_3d*spacing_3d

    dim = img.GetDimension()
    if dim==4:
        array_4d = sitk.GetArrayFromImage(img)
        arr = array_4d.transpose(3,2,1,0) # re-orient the array to match the affine
    elif dim==3:
        array_3d = sitk.GetArrayFromImage(img)
        arr = array_3d.transpose(2,1,0) # re-orient the array to match the affine
    smoothed_arr = _smooth_array(arr, affine, fwhm=fwhm, ensure_finite=True, copy=True)
    if dim==4:
        smoothed_arr = smoothed_arr.transpose(3,2,1,0)
        smoothed_img = copyInfo_4DImage(sitk.GetImageFromArray(
            smoothed_arr, isVector=False), img, img)
    elif dim==3:
        smoothed_arr = smoothed_arr.transpose(2,1,0)
        smoothed_img = copyInfo_3DImage(sitk.GetImageFromArray(
            smoothed_arr, isVector=False), img)

    return smoothed_img


def get_background_mask(bold_file, plotting=False):
    """
    Function that takes a 4D EPI timeseries and computes a background mask 
    excluding contributions from biological tissues.
    """
    import tempfile
    import matplotlib.pyplot as plt
    from rabies.utils import run_command
    from rabies.utils import copyInfo_3DImage
    tmppath = tempfile.mkdtemp()

    data_img = sitk.ReadImage(bold_file, sitk.sitkFloat32)
    data_array = sitk.GetArrayFromImage(data_img)
    median = np.median(data_array, axis=0)
    temporal_std = data_array.std(axis=0)
    
    median_img = copyInfo_3DImage(sitk.GetImageFromArray(
        median, isVector=False), data_img)
    sitk.WriteImage(median_img, f'{tmppath}/median.nii.gz')
    
    std_img = copyInfo_3DImage(sitk.GetImageFromArray(
        temporal_std, isVector=False), data_img)
    sitk.WriteImage(std_img, f'{tmppath}/std_map.nii.gz')
    
    ### first iteration of otsu thresholding, taking the background
    command = f'ThresholdImage 3 {tmppath}/median.nii.gz {tmppath}/weight.nii.gz Otsu 4'
    rc = run_command(command)
    otsu_img = sitk.ReadImage(f'{tmppath}/weight.nii.gz', sitk.sitkFloat32)
    otsu_array = sitk.GetArrayFromImage(otsu_img)
    background_mask = otsu_array<1
    
    not_zeros_indices = temporal_std!=0
    background_mask *= not_zeros_indices # make sure there are no 0s
    
    volume_img = copyInfo_3DImage(sitk.GetImageFromArray(
        background_mask.astype(float), isVector=False), otsu_img)
    sitk.WriteImage(volume_img, f'{tmppath}/background_mask.nii.gz')
    
    ### second iteration of otsu thresholding, taking the lower distribution
    # there were residual effects from the brain influencing the distribution
    command = f'ThresholdImage 3 {tmppath}/std_map.nii.gz {tmppath}/weight2.nii.gz Otsu 1 {tmppath}/background_mask.nii.gz'
    rc = run_command(command)
    otsu_img = sitk.ReadImage(f'{tmppath}/weight2.nii.gz', sitk.sitkFloat32)
    otsu_array = sitk.GetArrayFromImage(otsu_img)
    background_mask *= (otsu_array<2)
    mask_img = copyInfo_3DImage(sitk.GetImageFromArray(
        background_mask.astype(float), isVector=False), otsu_img)
    sitk.WriteImage(mask_img, f'{tmppath}/background_mask.nii.gz')

    from rabies.visualization import otsu_scaling, plot_3d
    
    fig_path = None
    if plotting:
        import pathlib
        filename_split = pathlib.Path(bold_file).name.rsplit(".nii")[0]

        fig = plt.figure(figsize=(24,6))
        #fig.suptitle(name, fontsize=30, color='white')
        ax1 = fig.add_subplot(3,2,1)
        ax2 = fig.add_subplot(3,2,3)
        ax3 = fig.add_subplot(3,2,5)
        ax4 = fig.add_subplot(1,4,3)

        ax4.hist(data_array[:,background_mask].flatten(), bins=100)
        ax4.set_title('Noise Distribution (should be Rician)', fontsize=20, color='white')

        
        from rabies.visualization import otsu_scaling, plot_3d
        planes = ('sagittal', 'coronal', 'horizontal')
        scaled = otsu_scaling(f'{tmppath}/std_map.nii.gz')
        
        plot_3d([ax1,ax2,ax3],scaled,fig,vmin=0,vmax=1,cmap='viridis', alpha=1, cbar=False, num_slices=6, planes=planes)
        plot_3d([ax1,ax2,ax3],mask_img,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.3, cbar=False, num_slices=6, planes=planes)
        ax1.set_title('Background masking (in red) over standard deviation map', fontsize=20, color='white')

        fig_path = os.path.abspath(f'{filename_split}_background_masking.png')
        fig.savefig(fig_path, bbox_inches='tight')

    return background_mask, data_array, fig_path
        
