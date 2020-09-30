import os
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
    cleaned_file = aroma_out+'/%s_aroma.nii.gz' % (filename_split[0])

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


def scrubbing(img, FD_file, scrubbing_threshold, timeseries_interval):
    '''
    Scrubbing based on FD: The frames that exceed the given threshold together with 1 back
    and 2 forward frames will be masked out from the data (as in Power et al. 2012)
    '''
    import numpy as np
    import nibabel as nb
    import pandas as pd
    mean_FD = pd.read_csv(FD_file).get('Mean')
    cutoff = np.asarray(mean_FD) >= scrubbing_threshold
    mask = np.ones(len(mean_FD))
    for i in range(len(mask)):
        if cutoff[i]:
            mask[i-1:i+2] = 0

    if not timeseries_interval == 'all':
        lowcut = int(timeseries_interval.split(',')[0])
        highcut = int(timeseries_interval.split(',')[1])
        mask = mask[lowcut:highcut]

    masked_img = np.asarray(img.dataobj)[:, :, :, mask.astype(bool)]
    return nb.Nifti1Image(masked_img, img.affine, img.header)


def select_timeseries(bold_file, timeseries_interval):
    if timeseries_interval == 'all':
        return bold_file
    else:
        import os
        import numpy as np
        import nibabel as nb
        img = nb.load(bold_file)
        lowcut = int(timeseries_interval.split(',')[0])
        highcut = int(timeseries_interval.split(',')[1])
        bold_file = os.path.abspath('selected_timeseries.nii.gz')
        nb.Nifti1Image(np.asarray(img.dataobj)[
                       :, :, :, lowcut:highcut], img.affine, img.header).to_filename(bold_file)
        return bold_file


def regress(bold_file, brain_mask_file, confounds_file, FD_file, conf_list, TR, lowpass, highpass, smoothing_filter,
            apply_scrubbing, scrubbing_threshold, timeseries_interval):
    import os
    import numpy as np
    import pandas as pd
    import nilearn.image
    from rabies.conf_reg_pkg.utils import scrubbing

    if ('mot_6' in conf_list) and ('mot_24' in conf_list):
        raise ValueError(
            "Can't select both the mot_6 and mot_24 options; must pick one.")

    cr_out = os.getcwd()
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

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
            print('Applying aCompCor with '+len(aCompCor_keys)+' components.')
            conf_keys += aCompCor_keys
        elif conf == 'mean_FD':
            mean_FD = pd.read_csv(FD_file).get('Mean')
            confounds['mean_FD'] = mean_FD
            conf_keys += [conf]
        else:
            conf_keys += [conf]

    confounds_array = np.asarray(confounds[conf_keys])

    '''
    Evaluate the variance explained (VE) by each regressor
    '''
    # functions that computes the Least Squares Estimates
    def closed_form(X, Y, intercept=False):
        if intercept:
            X = np.concatenate((X, np.ones([X.shape[0], 1])), axis=1)
        return np.linalg.inv(X.transpose().dot(X)).dot(X.transpose()).dot(Y)

    def linear_regression_var_explained(X, Y):
        X = (X-X.mean(axis=0))/X.std(axis=0)
        Y = (Y-Y.mean(axis=0))/Y.std(axis=0)
        # remove null values which may result from 0 in the array`
        X[np.isnan(X)] = 0
        Y[np.isnan(Y)] = 0
        # for each observation, it's values can be expressed through a linear combination of the predictors
        # take a bias into account in the model
        w = closed_form(X, Y, intercept=True)

        # add an intercept
        X_i = np.concatenate((X, np.ones([X.shape[0], 1])), axis=1)
        # the residuals after regressing out the predictors
        residuals = (Y-np.matmul(X_i, w))
        # mean square error after regression, for each observation independently
        MSE = np.mean((residuals**2), axis=0)
        # Original variance in each observation before regression
        TV = np.mean((Y-Y.mean(axis=0))**2, axis=0)

        VE = MSE/TV  # total variance explained in each observation
        VE[np.isnan(VE)] = 0

        # mean square error after regression, for each observation independently
        MSE = np.mean((residuals**2))
        # Original variance in each observation before regression
        TV = np.mean((Y-Y.mean(axis=0))**2)
        VE_tot = MSE/TV

        w = w[:-1, :]  # take out the intercept
        # now evaluate the portion of variance explained by each predictor,
        # scale all beta values to relative percentages of variance explained,
        # with the assumption that each beta value represents directly the relative contribution of each predictor
        # (guaranteed by the standardization of across feature)
        # Original variance in each observation before regression
        TV = np.sum((w-w.mean(axis=0))**2, axis=0)
        pred_VE = (w-w.mean(axis=0))**2/TV
        pred_VE[np.isnan(pred_VE)] = 0
        VE_observations = pred_VE*VE  # scale by the proportion of total variance explained

        # evaluate also variance explained from the entire data
        pred_variance = w.var(axis=1)  # variance along each dimension
        # evaluate the sum of variance present in this new dimensional space
        TV = pred_variance.sum()
        pred_VE = pred_variance/TV  # portion of variance explained from each dimension
        total_VE = pred_VE*VE_tot
        return total_VE, VE_observations, residuals

    import nibabel as nb
    brain_mask = np.asarray(nb.load(brain_mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)

    data_array = np.asarray(nb.load(bold_file).dataobj)
    timeseries = np.zeros([data_array.shape[3], volume_indices.sum()])
    for i in range(data_array.shape[3]):
        timeseries[i, :] = (data_array[:, :, :, i])[volume_indices]

    if not timeseries_interval == 'all':
        lowcut = int(timeseries_interval.split(',')[0])
        highcut = int(timeseries_interval.split(',')[1])
        confounds_array = confounds_array[lowcut:highcut, :]
        timeseries = timeseries[lowcut:highcut, :]

    total_VE, VE_observations, residuals = linear_regression_var_explained(
        confounds_array, timeseries)

    VE_dict = {}
    i = 0
    for VE, conf in zip(total_VE, conf_keys):
        VE_dict[conf] = VE_observations[i, :]
        print(conf+' explains '+str(round(VE, 3)*100)+'% of the variance.')
        i += 1

    import pickle
    VE_file = cr_out+'/'+filename_split[0]+'_VE_dict.pkl'
    with open(VE_file, 'wb') as handle:
        pickle.dump(VE_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # cleaning includes detrending, standardization
    if len(conf_list) > 0:
        cleaned = nilearn.image.clean_img(bold_file, detrend=True, standardize=True, low_pass=lowpass,
                                          high_pass=highpass, confounds=confounds_array, t_r=TR, mask_img=brain_mask_file)
    else:
        cleaned = nilearn.image.clean_img(bold_file, detrend=True, standardize=True,
                                          low_pass=lowpass, high_pass=highpass, confounds=None, t_r=TR, mask_img=brain_mask_file)

    if apply_scrubbing:
        cleaned = scrubbing(
            cleaned, FD_file, scrubbing_threshold, timeseries_interval)

    if smoothing_filter is not None:
        cleaned = nilearn.image.smooth_img(cleaned, smoothing_filter)

    cleaned_path = cr_out+'/'+filename_split[0]+'_cleaned.nii.gz'
    cleaned.to_filename(cleaned_path)
    return cleaned_path, bold_file, VE_file


class data_diagnosisInputSpec(BaseInterfaceInputSpec):
    bold_file = File(exists=True, mandatory=True,
                     desc='input BOLD time series')
    cleaned_path = File(exists=True, mandatory=True,
                        desc='ref file to realignment time series')
    brain_mask_file = File(exists=True, mandatory=True,
                           desc='ref file to realignment time series')


class data_diagnosisOutputSpec(TraitedSpec):
    mel_out = traits.Str(exists=True, desc="motion corrected time series")
    tSNR_file = File(exists=True, desc="motion corrected time series")


class data_diagnosis(BaseInterface):
    """

    """

    input_spec = data_diagnosisInputSpec
    output_spec = data_diagnosisOutputSpec

    def _run_interface(self, runtime):

        import os
        import nibabel as nb
        import numpy as np
        mel_out = os.path.abspath('melodic.ica/')
        os.mkdir(mel_out)
        command = 'melodic -i %s -o %s -m %s --report' % (
            self.inputs.cleaned_path, mel_out, self.inputs.brain_mask_file)
        os.system(command)
        img = nb.load(self.inputs.bold_file)
        array = np.asarray(img.dataobj)
        mean = array.mean(axis=3)
        std = array.std(axis=3)
        tSNR = np.divide(mean, std)
        tSNR_file = os.path.abspath('tSNR.nii.gz')
        nb.Nifti1Image(tSNR, img.affine, img.header).to_filename(tSNR_file)

        setattr(self, 'tSNR_file', tSNR_file)
        setattr(self, 'mel_out', mel_out)

        return runtime

    def _list_outputs(self):
        return {'tSNR_file': getattr(self, 'tSNR_file'),
                'mel_out': getattr(self, 'mel_out')}
