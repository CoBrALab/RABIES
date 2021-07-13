import os
import numpy as np
import SimpleITK as sitk
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


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


def exec_ICA_AROMA(inFile, mc_file, brain_mask, csf_mask, tr, aroma_dim):
    import os
    from rabies.conf_reg_pkg.utils import csv2par
    from rabies.conf_reg_pkg.mod_ICA_AROMA.ICA_AROMA_functions import run_ICA_AROMA
    import pathlib
    filename_split = pathlib.Path(inFile).name.rsplit(".nii")
    aroma_out = os.getcwd()+'/aroma_out'
    cleaned_file = aroma_out+f'/{filename_split[0]}_aroma.nii.gz'

    run_ICA_AROMA(aroma_out, os.path.abspath(inFile), mc=csv2par(mc_file), TR=float(tr), mask=os.path.abspath(
        brain_mask), mask_csf=os.path.abspath(csf_mask), denType="nonaggr", melDir="", dim=str(aroma_dim), overwrite=True)
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
    and 4 forward frames will be masked out from the data (as in Power et al. 2012)
    '''
    import numpy as np
    cutoff = np.asarray(FD_trace) >= scrubbing_threshold
    mask = np.ones(len(FD_trace)).astype(bool)
    for i in range(len(mask)):
        if cutoff[i]:
            mask[i-1:i+4] = 0
    return mask


def prep_CR(bold_file, brain_mask_file, confounds_file, FD_file, cr_opts):
    import os
    import numpy as np
    import pandas as pd
    import SimpleITK as sitk
    from rabies.conf_reg_pkg.utils import select_confound_timecourses

    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

    if len(cr_opts.conf_list)==0:
        confounds_array = select_confound_timecourses(['mot_6'],confounds_file,FD_file)
    else:
        confounds_array = select_confound_timecourses(cr_opts.conf_list,confounds_file,FD_file)

    FD_trace = pd.read_csv(FD_file).get('Mean')

    # select the subset of timeseries specified
    if not cr_opts.timeseries_interval == 'all':
        lowcut = int(cr_opts.timeseries_interval.split(',')[0])
        highcut = int(cr_opts.timeseries_interval.split(',')[1])
        confounds_array = confounds_array[lowcut:highcut, :]
        FD_trace = FD_trace[lowcut:highcut]
        time_range = range(lowcut,highcut)
    else:
        time_range = range(sitk.ReadImage(bold_file).GetSize()[3])

    data_dict = {'FD_trace':FD_trace, 'confounds_array':confounds_array, 'confounds_csv':confounds_file, 'time_range':time_range}
    return data_dict

def temporal_filtering(timeseries, data_dict, TR, lowpass, highpass,
        FD_censoring, FD_threshold, DVARS_censoring, minimum_timepoint):
    FD_trace=data_dict['FD_trace']
    confounds_array=data_dict['confounds_array']

    # apply highpass/lowpass filtering
    import nilearn.signal
    timeseries = nilearn.signal.clean(timeseries, detrend=False, standardize=False, low_pass=lowpass,
                                      high_pass=highpass, confounds=None, t_r=TR)

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
        return None
    timeseries=timeseries[frame_mask,:]
    confounds_array=confounds_array[frame_mask,:]
    FD_trace=FD_trace[frame_mask]
    DVARS=DVARS[frame_mask]

    # apply simple detrending, after censoring
    from scipy.signal import detrend
    timeseries = detrend(timeseries,axis=0)
    confounds_array = detrend(confounds_array,axis=0) # apply detrending to the confounds too, like in nilearn's function

    data_dict = {'timeseries':timeseries,'FD_trace':FD_trace, 'DVARS':DVARS, 'frame_mask':frame_mask, 'confounds_array':confounds_array}
    return data_dict


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


def regress(bold_file, data_dict, brain_mask_file, cr_opts):
    import os
    import numpy as np
    import SimpleITK as sitk
    from rabies.conf_reg_pkg.utils import recover_3D,recover_4D,temporal_filtering
    from rabies.analysis_pkg.analysis_functions import closed_form

    FD_trace=data_dict['FD_trace']
    confounds_array=data_dict['confounds_array']
    confounds_file=data_dict['confounds_csv']
    time_range=data_dict['time_range']

    cr_out = os.getcwd()
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

    brain_mask = sitk.GetArrayFromImage(sitk.ReadImage(brain_mask_file, sitk.sitkFloat32))
    volume_indices = brain_mask.astype(bool)

    data_array = sitk.GetArrayFromImage(sitk.ReadImage(bold_file, sitk.sitkFloat32))
    num_volumes = data_array.shape[0]
    timeseries = np.zeros([num_volumes, volume_indices.sum()])
    for i in range(num_volumes):
        timeseries[i, :] = (data_array[i, :, :, :])[volume_indices]
    timeseries = timeseries[time_range,:]

    TR = float(cr_opts.TR.split('s')[0])
    data_dict = temporal_filtering(timeseries, data_dict, TR, cr_opts.lowpass, cr_opts.highpass,
            cr_opts.FD_censoring, cr_opts.FD_threshold, cr_opts.DVARS_censoring, cr_opts.minimum_timepoint)
    if data_dict is None:
        empty_img = sitk.GetImageFromArray(np.empty([1,1]))
        empty_file = os.path.abspath('empty.nii.gz')
        sitk.WriteImage(empty_img, empty_file)
        return empty_file,empty_file,None

    timeseries=data_dict['timeseries']
    FD_trace=data_dict['FD_trace']
    DVARS=data_dict['DVARS']
    frame_mask=data_dict['frame_mask']
    confounds_array=data_dict['confounds_array']

    # estimate the VE from the CR selection, or 6 rigid motion parameters if no CR is applied
    X=confounds_array
    Y=timeseries
    try:
        res = Y-X.dot(closed_form(X,Y))
    except:
        import logging
        log = logging.getLogger('root')
        log.debug("SINGULAR MATRIX ERROR DURING CONFOUND REGRESSION. THIS SCAN WILL BE REMOVED FROM FURTHER PROCESSING.")
        import SimpleITK as sitk
        empty_img = sitk.GetImageFromArray(np.empty([1,1]))
        empty_file = os.path.abspath('empty.nii.gz')
        sitk.WriteImage(empty_img, empty_file)
        return empty_file,empty_file,None

    VE_spatial = 1-(res.var(axis=0)/Y.var(axis=0))
    VE_temporal = 1-(res.var(axis=1)/Y.var(axis=1))

    if len(cr_opts.conf_list) > 0:
        timeseries = res
    if cr_opts.standardize:
        timeseries = (timeseries-timeseries.mean(axis=0))/timeseries.std(axis=0)

    VE_spatial_map = recover_3D(brain_mask_file, VE_spatial)
    timeseries_3d = recover_4D(brain_mask_file, timeseries, bold_file)
    cleaned_path = cr_out+'/'+filename_split[0]+'_cleaned.nii.gz'
    sitk.WriteImage(timeseries_3d, cleaned_path)
    VE_file_path = cr_out+'/'+filename_split[0]+'_VE_map.nii.gz'
    sitk.WriteImage(VE_spatial_map, VE_file_path)

    if cr_opts.smoothing_filter is not None:
        import nilearn.image
        import nibabel as nb
        timeseries_3d = nilearn.image.smooth_img(nb.load(cleaned_path), cr_opts.smoothing_filter)
        timeseries_3d.to_filename(cleaned_path)

    data_dict = {'FD_trace':FD_trace, 'DVARS':DVARS, 'time_range':time_range, 'frame_mask':frame_mask, 'confounds_array':confounds_array, 'VE_temporal':VE_temporal, 'confounds_csv':confounds_file}
    return cleaned_path, VE_file_path, data_dict
