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


def exec_ICA_AROMA(inFile, outDir, mc_file, brain_mask, csf_mask, tr, aroma_dim):
    import os
    from rabies.conf_reg_pkg.mod_ICA_AROMA.ICA_AROMA_functions import run_ICA_AROMA
    run_ICA_AROMA(os.path.abspath(outDir), os.path.abspath(inFile), mc=os.path.abspath(mc_file), TR=float(tr), mask=os.path.abspath(
        brain_mask), mask_csf=os.path.abspath(csf_mask), denType="nonaggr", melDir="", dim=str(aroma_dim), overwrite=True)
    return os.path.abspath(outDir+'/denoised_func_data_nonaggr.nii.gz')


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


def regress(bold_file, brain_mask_file, confounds_file, csf_mask, FD_file, conf_list, TR, lowpass, highpass, smoothing_filter,
            run_aroma, aroma_dim, apply_scrubbing, scrubbing_threshold, timeseries_interval):
    import os
    import pandas as pd
    import numpy as np
    import nilearn.image
    from rabies.conf_reg_pkg.utils import scrubbing, exec_ICA_AROMA, csv2par

    cr_out = os.getcwd()

    confounds = pd.read_csv(confounds_file)
    keys = confounds.keys()
    confounds_list = []
    for conf in conf_list:
        if conf == 'mot_6':
            motion_keys = ['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3']
            for mov in motion_keys:
                confounds_list.append(np.asarray(confounds.get(mov)))
        elif conf == 'mot_24':
            motion_keys = [s for s in keys if "rot" in s or "mov" in s]
            for mov in motion_keys:
                confounds_list.append(np.asarray(confounds.get(mov)))
        elif conf == 'aCompCor':
            aCompCor_keys = [s for s in keys if "aCompCor" in s]
            print('Applying aCompCor with '+len(aCompCor_keys)+' components.')
            for aCompCor in aCompCor_keys:
                confounds_list.append(np.asarray(confounds.get(aCompCor)))
        elif conf == 'mean_FD':
            mean_FD = pd.read_csv(FD_file).get('Mean')
            confounds_list.append(np.asarray(mean_FD))
        else:
            confounds_list.append(np.asarray(confounds.get(conf)))

    '''
    what would be nice would be to have a print out of the variance explained for each regressor, to confirm it accounts for something
    '''

    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

    cleaning_input = bold_file

    # including detrending, standardization
    if run_aroma:
        aroma_out = cr_out+'/%s_aroma' % (filename_split[0])
        cleaning_input = exec_ICA_AROMA(cleaning_input, aroma_out, csv2par(
            confounds_file), brain_mask_file, csf_mask, TR, aroma_dim)
    if len(confounds_list) > 0:
        confounds_array = np.transpose(np.asarray(confounds_list))
        if not timeseries_interval == 'all':
            lowcut = int(timeseries_interval.split(',')[0])
            highcut = int(timeseries_interval.split(',')[1])
            confounds_array = confounds_array[lowcut:highcut, :]
        cleaned = nilearn.image.clean_img(cleaning_input, detrend=True, standardize=True, low_pass=lowpass,
                                          high_pass=highpass, confounds=confounds_array, t_r=TR, mask_img=brain_mask_file)
    else:
        cleaned = nilearn.image.clean_img(cleaning_input, detrend=True, standardize=True,
                                          low_pass=lowpass, high_pass=highpass, confounds=None, t_r=TR, mask_img=brain_mask_file)

    if smoothing_filter is not None:
        cleaned = nilearn.image.smooth_img(cleaned, smoothing_filter)

    if apply_scrubbing:
        cleaned = scrubbing(
            cleaned, FD_file, scrubbing_threshold, timeseries_interval)
    cleaned_path = cr_out+'/'+filename_split[0]+'_cleaned.nii.gz'
    cleaned.to_filename(cleaned_path)
    return cleaned_path, bold_file, cr_out


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
