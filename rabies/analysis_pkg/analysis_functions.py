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
    rc,c_out = run_command(command)
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
    rc,c_out = run_command(command)
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

class ComplementaryPCAInputSpec(BaseInterfaceInputSpec):
    dict_file = File(exists=True, mandatory=True, desc="Dictionary with prepared analysis data.")
    prior_bold_idx = traits.List(desc="The index for the ICA components that correspond to bold sources.")
    CPCA_temporal_comp = traits.Int(
        desc="number of data-driven temporal components to compute.")
    CPCA_spatial_comp = traits.Int(
        desc="number of data-driven spatial components to compute.")
    optimize_CPCA_dict = traits.Dict(
        desc="Dictionary with options for optimizing CPCA convergence.")
    network_weighting = traits.Str(
        desc="Whether to derive absolute or relative (variance-normalized) network maps.")
    figure_format = traits.Str(
        desc="Select file format for figures.")

class ComplementaryPCAOutputSpec(TraitedSpec):
    CPCA_prior_timecourse_csv = File(
        exists=True, desc=".csv with timecourses from the fitted prior sources.")
    CPCA_extra_timecourse_csv = File(
        exists=True, desc=".csv with timecourses from the scan-specific extra sources.")
    CPCA_prior_filename = File(
        exists=True, desc=".nii file with spatial components from the fitted prior sources.")
    CPCA_extra_filename = File(
        exists=True, desc=".nii file with spatial components from the scan-specific extra sources.")
    optimize_report = traits.Any(
        exists=False, desc="The CPCA optimization report.")


class ComplementaryPCA(BaseInterface):
    """

    """

    input_spec = ComplementaryPCAInputSpec
    output_spec = ComplementaryPCAOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import pandas as pd
        import pathlib  # Better path manipulation
        from rabies.utils import recover_4D
        from rabies.analysis_pkg.CPCA import spatiotemporal_CPCA

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

        CPCA_temporal_comp = self.inputs.CPCA_temporal_comp
        CPCA_spatial_comp = self.inputs.CPCA_spatial_comp
        if CPCA_temporal_comp<1: # make sure there is no negative number
            CPCA_temporal_comp=0
        if CPCA_spatial_comp<1: # make sure there is no negative number
            CPCA_spatial_comp=0
        optimize_CPCA_dict = self.inputs.optimize_CPCA_dict
        modeling, Ct_extra, Cs_extra, optimize_report_fig = spatiotemporal_CPCA(timeseries, C_prior, 
                                                                num_W=CPCA_temporal_comp, 
                                                                num_C=CPCA_spatial_comp, 
                                                                optim_dim=optimize_CPCA_dict['apply'],
                                                                min_prior_corr=optimize_CPCA_dict['min_prior_corr'], 
                                                                diff_thresh_t=optimize_CPCA_dict['diff_thresh_t'], 
                                                                diff_thresh_s=optimize_CPCA_dict['diff_thresh_s']) 
        
        if optimize_CPCA_dict['apply']:
            optimize_report_file = os.path.abspath(f'{filename_split[0]}_CPCA_optimize.{self.inputs.figure_format}')
            optimize_report_fig.savefig(optimize_report_file, bbox_inches='tight')
        else:
            optimize_report_file=None

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

        CPCA_prior_timecourse_csv = os.path.abspath(filename_split[0]+'_CPCA_prior_timecourse.csv')
        pd.DataFrame(modeling['W_fitted_prior']).to_csv(CPCA_prior_timecourse_csv, header=False, index=False)

        CPCA_extra_timecourse_csv = os.path.abspath(filename_split[0]+'_CPCA_extra_timecourse.csv')
        pd.DataFrame(W_extra).to_csv(CPCA_extra_timecourse_csv, header=False, index=False)

        CPCA_prior_filename = os.path.abspath(filename_split[0]+'_CPCA_prior.nii.gz')
        sitk.WriteImage(recover_4D(mask_file,C_fit.T, bold_file), CPCA_prior_filename)

        if (self.inputs.CPCA_temporal_comp+self.inputs.CPCA_spatial_comp)>0:
            CPCA_extra_filename = os.path.abspath(filename_split[0]+'_CPCA_extra.nii.gz')
            sitk.WriteImage(recover_4D(mask_file,C_extra.T, bold_file), CPCA_extra_filename)
        else:
            empty_img = sitk.GetImageFromArray(np.empty([1,1]))
            empty_file = os.path.abspath('empty.nii.gz')
            sitk.WriteImage(empty_img, empty_file)
            CPCA_extra_filename = empty_file

        setattr(self, 'CPCA_prior_timecourse_csv', CPCA_prior_timecourse_csv)
        setattr(self, 'CPCA_extra_timecourse_csv', CPCA_extra_timecourse_csv)
        setattr(self, 'CPCA_prior_filename', CPCA_prior_filename)
        setattr(self, 'CPCA_extra_filename', CPCA_extra_filename)
        setattr(self, 'optimize_report', optimize_report_file)

        return runtime

    def _list_outputs(self):
        return {'CPCA_prior_timecourse_csv': getattr(self, 'CPCA_prior_timecourse_csv'),
                'CPCA_extra_timecourse_csv': getattr(self, 'CPCA_extra_timecourse_csv'),
                'CPCA_prior_filename': getattr(self, 'CPCA_prior_filename'),
                'CPCA_extra_filename': getattr(self, 'CPCA_extra_filename'),
                'optimize_report': getattr(self, 'optimize_report'),
                }
