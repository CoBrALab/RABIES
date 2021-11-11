import os
import numpy as np
import SimpleITK as sitk


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
    from rabies.conf_reg_pkg.utils import csv2par
    from rabies.conf_reg_pkg.mod_ICA_AROMA.ICA_AROMA_functions import run_ICA_AROMA
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

    success = run_ICA_AROMA(aroma_out, os.path.abspath(inFile), mc=csv2par(mc_file), TR=float(tr), mask=os.path.abspath(
        brain_mask), mask_csf=os.path.abspath(csf_mask), denType="nonaggr", melDir="", dim=str(aroma_dim), overwrite=True, random_seed=random_seed)
    if not success:
        return None, aroma_out
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
    from rabies.conf_reg_pkg.utils import select_confound_timecourses

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
        import logging
        log = logging.getLogger('root')
        log.info(f"FD/DVARS CENSORING LEFT LESS THAN {str(minimum_timepoint)} VOLUMES. THIS SCAN WILL BE REMOVED FROM FURTHER PROCESSING.")
        return None,None,None

    return frame_mask,FD_trace,DVARS


def recover_3D(mask_file, vector_map):
    from rabies.preprocess_pkg.utils import copyInfo_3DImage
    mask_img = sitk.ReadImage(mask_file, sitk.sitkFloat32)
    brain_mask = sitk.GetArrayFromImage(mask_img)
    volume_indices=brain_mask.astype(bool)
    volume=np.zeros(brain_mask.shape)
    volume[volume_indices]=vector_map
    volume_img = copyInfo_3DImage(sitk.GetImageFromArray(
        volume, isVector=False), mask_img)
    return volume_img

def recover_4D(mask_file, vector_maps, ref_4d):
    from rabies.preprocess_pkg.utils import copyInfo_4DImage
    #vector maps of shape num_volumeXnum_voxel
    mask_img = sitk.ReadImage(mask_file, sitk.sitkFloat32)
    brain_mask = sitk.GetArrayFromImage(mask_img)
    volume_indices=brain_mask.astype(bool)
    shape=(vector_maps.shape[0],brain_mask.shape[0],brain_mask.shape[1],brain_mask.shape[2])
    volumes=np.zeros(shape)
    for i in range(vector_maps.shape[0]):
        volume=volumes[i,:,:,:]
        volume[volume_indices]=vector_maps[i,:]
        volumes[i,:,:,:]=volume
    volume_img = copyInfo_4DImage(sitk.GetImageFromArray(
        volumes, isVector=False), mask_img, sitk.ReadImage(ref_4d))
    return volume_img

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
            import logging
            log = logging.getLogger('root')
            log.info('Applying aCompCor with '+len(aCompCor_keys)+' components.')
            conf_keys += aCompCor_keys
        elif conf == 'mean_FD':
            mean_FD = pd.read_csv(FD_file).get('Mean')
            confounds['mean_FD'] = mean_FD
            conf_keys += [conf]
        else:
            conf_keys += [conf]

    return np.asarray(confounds[conf_keys])


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
