import numpy as np
import SimpleITK as sitk
from .analysis_math import vcorrcoef,closed_form


'''
seed-based FC
'''

def seed_based_FC(dict_file, seed_dict, seed_name):
    import os
    import numpy as np
    import SimpleITK as sitk
    import pathlib
    from rabies.utils import run_command, recover_3D
    from rabies.analysis_pkg.analysis_math import vcorrcoef

    import pickle
    with open(dict_file, 'rb') as handle:
        data_dict = pickle.load(handle)
    bold_file = data_dict['bold_file']
    mask_file = data_dict['mask_file']
    timeseries = data_dict['timeseries']
    volume_indices = data_dict['volume_indices']

    seed_file = seed_dict[seed_name]
    resampled = os.path.abspath('resampled.nii.gz')
    command=f'antsApplyTransforms -i {seed_file} -r {mask_file} -o {resampled} -n GenericLabel'
    rc = run_command(command)
    roi_mask = sitk.GetArrayFromImage(sitk.ReadImage(resampled))[volume_indices].astype(bool)

    # extract the voxel timeseries within the mask, and take the mean ROI timeseries
    seed_timeseries = timeseries[:,roi_mask].mean(axis=1)

    corrs = vcorrcoef(timeseries.T, seed_timeseries)
    corrs[np.isnan(corrs)] = 0

    corr_map_img = recover_3D(mask_file, corrs)
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")[0]
    corr_map_file = os.path.abspath(
        filename_split+'_'+seed_name+'_corr_map.nii.gz')
    
    sitk.WriteImage(corr_map_img, corr_map_file)

    # also save the seed timecourse
    import pandas as pd
    seed_timecourse_csv = os.path.abspath(
        filename_split+'_'+seed_name+'_timecourse.csv')
    pd.DataFrame(seed_timeseries).to_csv(seed_timecourse_csv, header=False, index=False)

    return corr_map_file,seed_timecourse_csv


'''
FC matrix
'''


def run_FC_matrix(dict_file, figure_format, roi_type='parcellated'):
    import os
    import pandas as pd
    import SimpleITK as sitk
    import numpy as np
    import pathlib  # Better path manipulation
    from rabies.analysis_pkg.analysis_functions import parcellated_FC_matrix, plot_matrix

    import pickle
    with open(dict_file, 'rb') as handle:
        data_dict = pickle.load(handle)


    bold_file = data_dict['bold_file']
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")
    figname = os.path.abspath(filename_split[0]+f'_FC_matrix.{figure_format}')
    
    timeseries = data_dict['timeseries']
    atlas_idx = data_dict['atlas_idx']
    roi_list = data_dict['roi_list']

    if roi_type == 'parcellated':
        corr_matrix,roi_labels = parcellated_FC_matrix(timeseries, atlas_idx, roi_list)
        matrix_df = pd.DataFrame(corr_matrix, index=roi_labels, columns=roi_labels)
    elif roi_type == 'voxelwise':
        corr_matrix = np.corrcoef(timeseries.T)
        matrix_df = pd.DataFrame(corr_matrix)
    else:
        raise ValueError(
            f"Invalid --ROI_type provided: {roi_type}. Must be either 'parcellated' or 'voxelwise.'")
    plot_matrix(figname, corr_matrix)

    data_file = os.path.abspath(filename_split[0]+'_FC_matrix.csv')
    matrix_df.to_csv(data_file, sep=',')
    return data_file, figname


def parcellated_FC_matrix(sub_timeseries, atlas_idx, roi_list):
    
    timeseries_dict = {}
    for i in roi_list:
        roi_mask = np.asarray(atlas_idx == i, dtype=bool)
        # extract the voxel timeseries within the mask, and take the mean ROI timeseries
        timeseries_dict[str(i)] = sub_timeseries[:,roi_mask].mean(axis=1)

    roi_labels = list(timeseries_dict.keys())
    sub_timeseries = []
    for roi in roi_labels:
        sub_timeseries.append(timeseries_dict[roi])
    corr_matrix = np.corrcoef(sub_timeseries)
    return corr_matrix,roi_labels


def plot_matrix(filename, corr_matrix):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    g = ax.imshow(corr_matrix, cmap='coolwarm', vmax=1, vmin=-1)
    ax.axis('off')
    cbar = plt.colorbar(g, ax=ax, shrink=0.5)
    cbar.set_label('R score', rotation=270, fontsize=10)
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight', dpi=150)


'''
ICA
'''


def run_group_ICA(bold_file_list, mask_file, dim, random_seed):
    import os
    import pandas as pd

    # create a filelist.txt
    file_path = os.path.abspath('filelist.txt')
    from rabies.utils import flatten_list
    merged = flatten_list(list(bold_file_list))
    df = pd.DataFrame(data=merged)
    df.to_csv(file_path, header=False, sep=',', index=False)

    from rabies.utils import run_command
    out_dir = os.path.abspath('group_melodic.ica')
    command = f'melodic -i {file_path} -m {mask_file} -o {out_dir} -d {dim} --report --seed={str(random_seed)}'
    rc = run_command(command)
    IC_file = out_dir+'/melodic_IC.nii.gz'
    return out_dir, IC_file


def run_DR_ICA(dict_file,network_weighting):
    import os
    import pandas as pd
    import pathlib  # Better path manipulation
    import SimpleITK as sitk
    from rabies.utils import recover_4D
    from rabies.analysis_pkg.analysis_math import dual_regression

    import pickle
    with open(dict_file, 'rb') as handle:
        data_dict = pickle.load(handle)
    bold_file = data_dict['bold_file']
    mask_file = data_dict['mask_file']
    timeseries = data_dict['timeseries']
    prior_map_vectors = data_dict['prior_map_vectors']

    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

    DR = dual_regression(prior_map_vectors, timeseries)

    if network_weighting=='absolute':
        DR_C = DR['C']*DR['S']
    elif network_weighting=='relative':
        DR_C = DR['C']
    else:
        raise 

    dual_regression_timecourse_csv = os.path.abspath(filename_split[0]+'_dual_regression_timecourse.csv')
    pd.DataFrame(DR['W']).to_csv(dual_regression_timecourse_csv, header=False, index=False)

    # save the subjects' IC maps as .nii file
    DR_maps_filename = os.path.abspath(filename_split[0]+'_DR_maps.nii.gz')
    sitk.WriteImage(recover_4D(mask_file, DR_C.T, bold_file), DR_maps_filename)
    return DR_maps_filename, dual_regression_timecourse_csv


from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

class NeuralPriorRecoveryInputSpec(BaseInterfaceInputSpec):
    dict_file = File(exists=True, mandatory=True, desc="Dictionary with prepared analysis data.")
    prior_bold_idx = traits.List(desc="The index for the ICA components that correspond to bold sources.")
    NPR_temporal_comp = traits.Int(
        desc="number of data-driven temporal components to compute.")
    NPR_spatial_comp = traits.Int(
        desc="number of data-driven spatial components to compute.")
    network_weighting = traits.Str(
        desc="Whether to derive absolute or relative (variance-normalized) network maps.")

class NeuralPriorRecoveryOutputSpec(TraitedSpec):
    NPR_prior_timecourse_csv = File(
        exists=True, desc=".csv with timecourses from the fitted prior sources.")
    NPR_extra_timecourse_csv = File(
        exists=True, desc=".csv with timecourses from the scan-specific extra sources.")
    NPR_prior_filename = File(
        exists=True, desc=".nii file with spatial components from the fitted prior sources.")
    NPR_extra_filename = File(
        exists=True, desc=".nii file with spatial components from the scan-specific extra sources.")


class NeuralPriorRecovery(BaseInterface):
    """

    """

    input_spec = NeuralPriorRecoveryInputSpec
    output_spec = NeuralPriorRecoveryOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import pandas as pd
        import pathlib  # Better path manipulation
        from rabies.utils import recover_4D
        from rabies.analysis_pkg.analysis_math import spatiotemporal_prior_fit

        import pickle
        with open(self.inputs.dict_file, 'rb') as handle:
            data_dict = pickle.load(handle)
        bold_file = data_dict['bold_file']
        mask_file = data_dict['mask_file']
        timeseries = data_dict['timeseries']
        prior_map_vectors = data_dict['prior_map_vectors']
        network_weighting=self.inputs.network_weighting

        filename_split = pathlib.Path(
            bold_file).name.rsplit(".nii")

        C_prior = prior_map_vectors[self.inputs.prior_bold_idx,:].T
        NPR_temporal_comp = self.inputs.NPR_temporal_comp
        NPR_spatial_comp = self.inputs.NPR_spatial_comp
        if NPR_temporal_comp<1: # make sure there is no negative number
            NPR_temporal_comp=0
        if NPR_spatial_comp<1: # make sure there is no negative number
            NPR_spatial_comp=0
        modeling = spatiotemporal_prior_fit(timeseries, C_prior, num_W=NPR_temporal_comp, num_C=NPR_spatial_comp)

        # put together the spatial and temporal extra components
        W_extra = np.concatenate((modeling['W_spatial'],modeling['W_temporal']), axis=1)

        if network_weighting=='absolute':
            C_extra = np.concatenate((modeling['C_spatial']*modeling['S_spatial'],
                modeling['C_temporal']*modeling['S_temporal']), axis=1)
            C_fit = modeling['C_fitted_prior']*modeling['S_fitted_prior']
        elif network_weighting=='relative':
            C_extra = np.concatenate((modeling['C_spatial'],
                modeling['C_temporal']), axis=1)
            C_fit = modeling['C_fitted_prior']
        else:
            raise 

        NPR_prior_timecourse_csv = os.path.abspath(filename_split[0]+'_NPR_prior_timecourse.csv')
        pd.DataFrame(modeling['W_fitted_prior']).to_csv(NPR_prior_timecourse_csv, header=False, index=False)

        NPR_extra_timecourse_csv = os.path.abspath(filename_split[0]+'_NPR_extra_timecourse.csv')
        pd.DataFrame(W_extra).to_csv(NPR_extra_timecourse_csv, header=False, index=False)

        NPR_prior_filename = os.path.abspath(filename_split[0]+'_NPR_prior.nii.gz')
        sitk.WriteImage(recover_4D(mask_file,C_fit.T, bold_file), NPR_prior_filename)

        if (self.inputs.NPR_temporal_comp+self.inputs.NPR_spatial_comp)>0:
            NPR_extra_filename = os.path.abspath(filename_split[0]+'_NPR_extra.nii.gz')
            sitk.WriteImage(recover_4D(mask_file,C_extra.T, bold_file), NPR_extra_filename)
        else:
            empty_img = sitk.GetImageFromArray(np.empty([1,1]))
            empty_file = os.path.abspath('empty.nii.gz')
            sitk.WriteImage(empty_img, empty_file)
            NPR_extra_filename = empty_file

        setattr(self, 'NPR_prior_timecourse_csv', NPR_prior_timecourse_csv)
        setattr(self, 'NPR_extra_timecourse_csv', NPR_extra_timecourse_csv)
        setattr(self, 'NPR_prior_filename', NPR_prior_filename)
        setattr(self, 'NPR_extra_filename', NPR_extra_filename)

        return runtime

    def _list_outputs(self):
        return {'NPR_prior_timecourse_csv': getattr(self, 'NPR_prior_timecourse_csv'),
                'NPR_extra_timecourse_csv': getattr(self, 'NPR_extra_timecourse_csv'),
                'NPR_prior_filename': getattr(self, 'NPR_prior_filename'),
                'NPR_extra_filename': getattr(self, 'NPR_extra_filename'),
                }
