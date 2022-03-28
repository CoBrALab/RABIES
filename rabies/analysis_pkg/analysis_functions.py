import numpy as np
import SimpleITK as sitk
from .analysis_math import vcorrcoef,closed_form


def resample_IC_file(in_file, ref_file, transforms = [], inverses = []):
    # resampling the reference image to the dimension of the EPI
    import SimpleITK as sitk
    import os
    from rabies.utils import split_volumes, copyInfo_4DImage, exec_applyTransforms
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(
        in_file).name.rsplit(".nii")
    out_file = os.path.abspath(filename_split[0])+'_resampled.nii.gz'

    # Splitting bold file into lists of single volumes
    [volumes_list, num_volumes] = split_volumes(
        in_file, "bold_", sitk.sitkFloat32)

    warped_volumes = []
    for x in range(0, num_volumes):
        warped_vol_fname = os.path.abspath(
            "deformed_volume" + str(x) + ".nii.gz")
        warped_volumes.append(warped_vol_fname)
        exec_applyTransforms(transforms=transforms, inverses=inverses, input_image=volumes_list[x], ref_image=ref_file, output_image=warped_vol_fname, mask=False)

    sample_volume = sitk.ReadImage(
        warped_volumes[0])
    shape = sitk.GetArrayFromImage(sample_volume).shape
    combined = np.zeros((num_volumes, shape[0], shape[1], shape[2]))

    i = 0
    for file in warped_volumes:
        combined[i, :, :, :] = sitk.GetArrayFromImage(
            sitk.ReadImage(file))[:, :, :]
        i = i+1
    if (i != num_volumes):
        raise ValueError("Error occured with Merge.")

    combined_image = sitk.GetImageFromArray(combined, isVector=False)

    # set metadata and affine for the newly constructed 4D image
    header_source = sitk.ReadImage(
        in_file)
    combined_image = copyInfo_4DImage(
        combined_image, sample_volume, header_source)
    sitk.WriteImage(combined_image, out_file)
    return out_file


'''
seed-based FC
'''

def seed_based_FC(bold_file, brain_mask, seed_dict, seed_name):
    import os
    import SimpleITK as sitk
    import pathlib
    from rabies.analysis_pkg.analysis_functions import seed_corr

    seed_file = seed_dict[seed_name]
    corr_map_img = seed_corr(bold_file, brain_mask, seed_file)

    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")[0]
    corr_map_file = os.path.abspath(
        filename_split+'_'+seed_name+'_corr_map.nii.gz')
    
    sitk.WriteImage(corr_map_img, corr_map_file)
    return corr_map_file


def seed_corr(bold_file, mask_file, seed):
    import os
    from rabies.utils import run_command, recover_3D
    
    mask_img = sitk.ReadImage(mask_file)
    brain_mask = sitk.GetArrayFromImage(mask_img)
    volume_indices = brain_mask.astype(bool)

    data_img = sitk.ReadImage(bold_file)
    data_array = sitk.GetArrayFromImage(data_img)
    num_volumes = data_array.shape[0]
    timeseries = np.zeros([num_volumes, volume_indices.sum()])
    for i in range(num_volumes):
        timeseries[i, :] = (data_array[i, :, :, :])[volume_indices]
        

    resampled = os.path.abspath('resampled.nii.gz')
    command=f'antsApplyTransforms -i {seed} -r {mask_file} -o {resampled} -n GenericLabel'
    rc = run_command(command)

    roi_mask = sitk.GetArrayFromImage(sitk.ReadImage(resampled))[volume_indices].astype(bool)

    # extract the voxel timeseries within the mask, and take the mean ROI timeseries
    seed_timeseries = timeseries[:,roi_mask].mean(axis=1)

    corrs = vcorrcoef(timeseries.T, seed_timeseries)
    corrs[np.isnan(corrs)] = 0

    corr_map_img = recover_3D(mask_file, corrs)
    return corr_map_img


'''
FC matrix
'''


def run_FC_matrix(bold_file, mask_file, atlas, roi_type='parcellated'):
    import os
    import pandas as pd
    import SimpleITK as sitk
    import numpy as np
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")
    figname = os.path.abspath(filename_split[0]+'_FC_matrix.png')

    from rabies.analysis_pkg.analysis_functions import parcellated_FC_matrix, plot_matrix
    
    mask_img = sitk.ReadImage(mask_file)
    brain_mask = sitk.GetArrayFromImage(mask_img)
    volume_indices = brain_mask.astype(bool)

    data_img = sitk.ReadImage(bold_file)
    data_array = sitk.GetArrayFromImage(data_img)
    num_volumes = data_array.shape[0]
    timeseries = np.zeros([num_volumes, volume_indices.sum()])
    for i in range(num_volumes):
        timeseries[i, :] = (data_array[i, :, :, :])[volume_indices]
    
    if roi_type == 'parcellated':
        corr_matrix = parcellated_FC_matrix(timeseries, volume_indices, atlas)
    elif roi_type == 'voxelwise':
        corr_matrix = np.corrcoef(timeseries.T)
    else:
        raise ValueError(
            f"Invalid --ROI_type provided: {roi_type}. Must be either 'parcellated' or 'voxelwise.'")
    plot_matrix(figname, corr_matrix)

    data_file = os.path.abspath(filename_split[0]+'_FC_matrix.csv')
    df = pd.DataFrame(corr_matrix)
    df.to_csv(data_file, sep=',')
    return data_file, figname


def parcellated_FC_matrix(sub_timeseries, volume_indices, atlas):

    atlas_data = sitk.GetArrayFromImage(sitk.ReadImage(atlas))[volume_indices]
    max_int = atlas_data.max()
    
    timeseries_dict = {}
    for i in range(1, max_int+1):
        if np.max(i == atlas_data):  # taking a ROI only if it has labeled voxels
            roi_mask = np.asarray(atlas_data == i, dtype=bool)
            # extract the voxel timeseries within the mask, and take the mean ROI timeseries
            timeseries_dict[str(i)] = sub_timeseries[:,roi_mask].mean(axis=1)
    
    roi_labels = timeseries_dict.keys()
    sub_timeseries = []
    for roi in roi_labels:
        sub_timeseries.append(timeseries_dict[roi])
    corr_matrix = np.corrcoef(sub_timeseries)
    return corr_matrix


def plot_matrix(filename, corr_matrix):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    vmax = np.abs(corr_matrix).max()
    vmin = -vmax
    g = ax.imshow(corr_matrix, cmap='coolwarm', vmax=vmax, vmin=vmin)
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


def run_DR_ICA(bold_file, mask_file, IC_file):
    import os
    import pandas as pd
    import pathlib  # Better path manipulation
    import numpy as np
    import SimpleITK as sitk
    from rabies.utils import recover_4D
    from rabies.analysis_pkg.analysis_functions import resample_IC_file, dual_regression
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

    # check if the IC_file has the same dimensions as bold_file
    if not sitk.ReadImage(bold_file).GetSize()[:-1]==sitk.ReadImage(IC_file).GetSize()[:-1]:
        from nipype import logging
        log = logging.getLogger('nipype.workflow')
        log.info('Resampling file with IC components to match the scan dimensionality.')
        IC_file = resample_IC_file(IC_file, mask_file, transforms = [], inverses = [])

    mask_img = sitk.ReadImage(mask_file)
    brain_mask = sitk.GetArrayFromImage(mask_img)
    volume_indices = brain_mask.astype(bool)

    data_img = sitk.ReadImage(bold_file)
    data_array = sitk.GetArrayFromImage(data_img)
    num_volumes = data_array.shape[0]
    timeseries = np.zeros([num_volumes, volume_indices.sum()])
    for i in range(num_volumes):
        timeseries[i, :] = (data_array[i, :, :, :])[volume_indices]

    data_img = sitk.ReadImage(IC_file)
    data_array = sitk.GetArrayFromImage(data_img)
    num_ICs = data_array.shape[0]
    all_IC_vectors = np.zeros([num_ICs, volume_indices.sum()])
    for i in range(num_ICs):
        all_IC_vectors[i, :] = (data_array[i, :, :, :])[volume_indices]

    DR = dual_regression(all_IC_vectors, timeseries)

    dual_regression_timecourse_csv = os.path.abspath(filename_split[0]+'_dual_regression_timecourse.csv')
    pd.DataFrame(DR['W']).to_csv(dual_regression_timecourse_csv, header=False, index=False)

    # save the subjects' IC maps as .nii file
    DR_maps_filename = os.path.abspath(filename_split[0]+'_DR_maps.nii.gz')
    sitk.WriteImage(recover_4D(mask_file, DR['C'], bold_file), DR_maps_filename)
    return DR_maps_filename, dual_regression_timecourse_csv


def dual_regression(all_IC_vectors, timeseries):
    ### compute dual regression
    ### Here, we adopt an approach where the algorithm should explain the data
    ### as a linear combination of spatial maps. The data itself, is only temporally
    ### detrended, and not spatially centered, which could cause inconsistencies during
    ### linear regression according to https://mandymejia.com/2018/03/29/the-role-of-centering-in-dual-regression/#:~:text=Dual%20regression%20requires%20centering%20across%20time%20and%20space&text=time%20points.,each%20time%20course%20at%20zero
    ### The fMRI timeseries aren't assumed theoretically to be spatially centered, and
    ### this measure would be removing global signal variations which we are interested in.
    ### Thus we prefer to avoid this step here, despite modelling limitations.
    X = all_IC_vectors.T
    Y = timeseries.T
    # for one given volume, it's values can be expressed through a linear combination of the components
    W = closed_form(X, Y, intercept=False).T
    W /= np.sqrt((W ** 2).sum(axis=0)) # the temporal domain is variance-normalized so that the weights are contained in the spatial maps
    # for a given voxel timeseries, it's signal can be explained a linear combination of the component timecourses
    C = closed_form(W, Y.T, intercept=False)
    DR = {'C':C, 'W':W}
    return DR


from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

class NeuralPriorRecoveryInputSpec(BaseInterfaceInputSpec):
    bold_file = File(exists=True, mandatory=True, desc="Timeseries to fit.")
    mask_file = File(exists=True, mandatory=True, desc="Brain mask.")
    prior_maps = File(exists=True, mandatory=True, desc="MELODIC ICA components to use.")
    prior_bold_idx = traits.List(desc="The index for the ICA components that correspond to bold sources.")
    NPR_temporal_comp = traits.Int(
        desc="number of data-driven temporal components to compute.")
    NPR_spatial_comp = traits.Int(
        desc="number of data-driven spatial components to compute.")

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
        filename_split = pathlib.Path(
            self.inputs.bold_file).name.rsplit(".nii")

        resampled_priors = resample_IC_file(self.inputs.prior_maps, self.inputs.mask_file)


        mask_img = sitk.ReadImage(self.inputs.mask_file)
        brain_mask = sitk.GetArrayFromImage(mask_img)
        volume_indices = brain_mask.astype(bool)

        data_img = sitk.ReadImage(self.inputs.bold_file)
        data_array = sitk.GetArrayFromImage(data_img)
        num_volumes = data_array.shape[0]
        timeseries = np.zeros([num_volumes, volume_indices.sum()])
        for i in range(num_volumes):
            timeseries[i, :] = (data_array[i, :, :, :])[volume_indices]

        data_img = sitk.ReadImage(resampled_priors)
        data_array = sitk.GetArrayFromImage(data_img)
        num_ICs = data_array.shape[0]
        all_IC_vectors = np.zeros([num_ICs, volume_indices.sum()])
        for i in range(num_ICs):
            all_IC_vectors[i, :] = (data_array[i, :, :, :])[volume_indices]

        C_prior = all_IC_vectors[self.inputs.prior_bold_idx,:].T
        NPR_temporal_comp = self.inputs.NPR_temporal_comp
        NPR_spatial_comp = self.inputs.NPR_spatial_comp
        if NPR_temporal_comp<1: # make sure there is no negative number
            NPR_temporal_comp=0
        if NPR_spatial_comp<1: # make sure there is no negative number
            NPR_spatial_comp=0
        modeling = spatiotemporal_prior_fit(timeseries, C_prior, num_W=NPR_temporal_comp, num_C=NPR_spatial_comp)

        # put together the spatial and temporal extra componentsZ
        W_extra = np.concatenate((modeling['W_spatial'],modeling['W_temporal']), axis=1)
        C_extra = np.concatenate((modeling['C_spatial'],modeling['C_temporal']), axis=1)

        NPR_prior_timecourse_csv = os.path.abspath(filename_split[0]+'_NPR_prior_timecourse.csv')
        pd.DataFrame(modeling['W_fitted_prior']).to_csv(NPR_prior_timecourse_csv, header=False, index=False)

        NPR_extra_timecourse_csv = os.path.abspath(filename_split[0]+'_NPR_extra_timecourse.csv')
        pd.DataFrame(W_extra).to_csv(NPR_extra_timecourse_csv, header=False, index=False)

        NPR_prior_filename = os.path.abspath(filename_split[0]+'_NPR_prior.nii.gz')
        sitk.WriteImage(recover_4D(self.inputs.mask_file,modeling['C_fitted_prior'].T, self.inputs.bold_file), NPR_prior_filename)

        if (self.inputs.NPR_temporal_comp+self.inputs.NPR_spatial_comp)>0:
            NPR_extra_filename = os.path.abspath(filename_split[0]+'_NPR_extra.nii.gz')
            sitk.WriteImage(recover_4D(self.inputs.mask_file,C_extra.T, self.inputs.bold_file), NPR_extra_filename)
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
