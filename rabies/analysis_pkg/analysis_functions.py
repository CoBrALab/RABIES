import numpy as np
import nibabel as nb
from .analysis_math import vcorrcoef,closed_form


def recover_3D(mask_file, vector_map):
    brain_mask = np.asarray(nb.load(mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)
    volume = np.zeros(brain_mask.shape)
    volume[volume_indices] = vector_map
    volume_img = nb.Nifti1Image(volume, nb.load(
        mask_file).affine, nb.load(mask_file).header)
    return volume_img


def recover_3D_multiple(mask_file, vector_maps):
    # vector maps of shape num_volumeXnum_voxel
    brain_mask = np.asarray(nb.load(mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)
    shape = (brain_mask.shape[0], brain_mask.shape[1],
             brain_mask.shape[2], vector_maps.shape[0])
    volumes = np.zeros(shape)
    for i in range(vector_maps.shape[0]):
        volume = volumes[:, :, :, i]
        volume[volume_indices] = vector_maps[i, :]
        volumes[:, :, :, i] = volume
    volume_img = nb.Nifti1Image(volumes, nb.load(
        mask_file).affine, nb.load(mask_file).header)
    return volume_img


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
    import nibabel as nb
    import numpy as np
    import pathlib
    from rabies.analysis_pkg.analysis_functions import seed_corr

    seed_file = seed_dict[seed_name]

    mask_array = np.asarray(nb.load(brain_mask).dataobj)
    mask_vector = mask_array.reshape(-1)
    mask_indices = (mask_vector==True)
    mask_vector = mask_vector.astype(float)

    mask_vector[mask_indices] = seed_corr(bold_file, brain_mask, seed_file)
    corr_map = mask_vector.reshape(mask_array.shape)

    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")[0]
    corr_map_file = os.path.abspath(
        filename_split+'_'+seed_name+'_corr_map.nii.gz')
    nb.Nifti1Image(corr_map, nb.load(brain_mask).affine, nb.load(
        brain_mask).header).to_filename(corr_map_file)
    return corr_map_file


def seed_corr(bold_file, mask_file, seed):
    import os
    from rabies.utils import run_command
    
    brain_mask = np.asarray(nb.load(mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)

    timeseries_array = np.asarray(nb.load(bold_file).dataobj)
    sub_timeseries = np.zeros(
        [timeseries_array.shape[3], volume_indices.sum()])
    for t in range(timeseries_array.shape[3]):
        sub_timeseries[t, :] = timeseries_array[:, :, :, t][volume_indices]
    
        
    resampled = os.path.abspath('resampled.nii.gz')
    command=f'antsApplyTransforms -i {seed} -r {mask_file} -o {resampled} -n GenericLabel'
    rc = run_command(command)

    roi_mask = np.asarray(nb.load(resampled).dataobj)[volume_indices].astype(bool)

    # extract the voxel timeseries within the mask, and take the mean ROI timeseries
    seed_timeseries = sub_timeseries[:,roi_mask].mean(axis=1)

    corrs = vcorrcoef(sub_timeseries.T, seed_timeseries)
    corrs[np.isnan(corrs)] = 0
    return corrs


'''
FC matrix
'''


def run_FC_matrix(bold_file, mask_file, atlas, roi_type='parcellated'):
    import os
    import pandas as pd
    import numpy as np
    import nibabel as nb
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")
    figname = os.path.abspath(filename_split[0]+'_FC_matrix.png')

    from rabies.analysis_pkg.analysis_functions import parcellated_FC_matrix, plot_matrix
    
    brain_mask = np.asarray(nb.load(mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)

    timeseries_array = np.asarray(nb.load(bold_file).dataobj)
    sub_timeseries = np.zeros(
        [timeseries_array.shape[3], volume_indices.sum()])
    for t in range(timeseries_array.shape[3]):
        sub_timeseries[t, :] = timeseries_array[:, :, :, t][volume_indices]
    
    
    if roi_type == 'parcellated':
        corr_matrix = parcellated_FC_matrix(sub_timeseries, volume_indices, atlas)
    elif roi_type == 'voxelwise':
        corr_matrix = np.corrcoef(sub_timeseries.T)
    else:
        raise ValueError(
            f"Invalid --ROI_type provided: {roi_type}. Must be either 'parcellated' or 'voxelwise.'")
    plot_matrix(figname, corr_matrix)

    data_file = os.path.abspath(filename_split[0]+'_FC_matrix.csv')
    df = pd.DataFrame(corr_matrix)
    df.to_csv(data_file, sep=',')
    return data_file, figname


def parcellated_FC_matrix(sub_timeseries, volume_indices, atlas):

    atlas_data = np.asarray(nb.load(atlas).dataobj)[volume_indices]
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


def resample_4D(input_4d, ref_file):
    import os
    import pathlib  # Better path manipulation
    import SimpleITK as sitk
    from rabies.utils import run_command, split_volumes, Merge, copyInfo_3DImage
    rabies_data_type=sitk.sitkFloat32


    # check if the IC_file has the same dimensions as bold_file
    img_array = sitk.GetArrayFromImage(
        sitk.ReadImage(ref_file))[0, :, :, :]
    image_3d = copyInfo_3DImage(sitk.GetImageFromArray(
        img_array, isVector=False), sitk.ReadImage(ref_file))
    new_ref = 'temp_ref.nii.gz'
    sitk.WriteImage(image_3d, 'temp_ref.nii.gz')

    filename_split = pathlib.Path(
        input_4d).name.rsplit(".nii")

    # Splitting into list of single volumes
    [split_volumes_files, num_volumes] = split_volumes(
        input_4d, "split_", rabies_data_type)

    resampled_volumes = []
    for x in range(0, num_volumes):
        resampled_vol_fname = os.path.abspath(
            "resampled_volume" + str(x) + ".nii.gz")
        resampled_volumes.append(resampled_vol_fname)

        command = f'antsApplyTransforms -i {split_volumes_files[x]} -n BSpline[5] -r {new_ref} -o {resampled_vol_fname}'
        rc = run_command(command)
        # change image to specified data type
        sitk.WriteImage(sitk.ReadImage(resampled_vol_fname, rabies_data_type), resampled_vol_fname)

    out=Merge(in_files=resampled_volumes, header_source=input_4d, rabies_data_type=rabies_data_type, clip_negative=False).run()
    return out.outputs.out_file

def run_DR_ICA(bold_file, mask_file, IC_file):
    import os
    import pandas as pd
    import pathlib  # Better path manipulation
    import numpy as np
    import nibabel as nb
    import SimpleITK as sitk
    from rabies.analysis_pkg.analysis_functions import recover_3D_multiple, resample_4D, dual_regression
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

    # check if the IC_file has the same dimensions as bold_file
    if not sitk.ReadImage(bold_file).GetSize()[:-1]==sitk.ReadImage(IC_file).GetSize()[:-1]:
        from nipype import logging
        log = logging.getLogger('nipype.workflow')
        log.info('Resampling file with IC components to match the scan dimensionality.')
        IC_file = resample_4D(IC_file, bold_file)

    brain_mask = np.asarray(nb.load(mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)

    timeseries_array = np.asarray(nb.load(bold_file).dataobj)
    sub_timeseries = np.zeros(
        [timeseries_array.shape[3], volume_indices.sum()])
    for t in range(timeseries_array.shape[3]):
        sub_timeseries[t, :] = timeseries_array[:, :, :, t][volume_indices]

    all_IC_array = np.asarray(nb.load(IC_file).dataobj)
    all_IC_vectors = np.zeros([all_IC_array.shape[3], volume_indices.sum()])
    for i in range(all_IC_array.shape[3]):
        all_IC_vectors[i, :] = (all_IC_array[:, :, :, i])[volume_indices]

    DR = dual_regression(all_IC_vectors, sub_timeseries)

    dual_regression_timecourse_csv = os.path.abspath(filename_split[0]+'_dual_regression_timecourse.csv')
    pd.DataFrame(DR['W']).to_csv(dual_regression_timecourse_csv, header=False, index=False)

    # save the subjects' IC maps as .nii file
    DR_maps_filename = os.path.abspath(filename_split[0]+'_DR_maps.nii.gz')
    recover_3D_multiple(mask_file,DR['C']).to_filename(DR_maps_filename)
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
    # normalize the component timecourses to unit variance
    W /= W.std(axis=0)
    # for a given voxel timeseries, it's signal can be explained a linear combination of the component timecourses
    C = closed_form(W, Y.T, intercept=False)
    DR = {'C':C, 'W':W}
    return DR


from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

class dual_ICA_wrapperInputSpec(BaseInterfaceInputSpec):
    bold_file = File(exists=True, mandatory=True, desc="Timeseries to fit.")
    mask_file = File(exists=True, mandatory=True, desc="Brain mask.")
    prior_maps = File(exists=True, mandatory=True, desc="MELODIC ICA components to use.")
    prior_bold_idx = traits.List(desc="The index for the ICA components that correspond to bold sources.")
    num_comp = traits.Int(
        desc="number of components to compute from dual ICA.")

class dual_ICA_wrapperOutputSpec(TraitedSpec):
    dual_ICA_timecourse_csv = File(
        exists=True, desc=".csv with timecourses from dual ICA")
    dual_ICA_filename = File(
        exists=True, desc=".nii file of the Dual ICA maps")

class dual_ICA_wrapper(BaseInterface):
    """

    """

    input_spec = dual_ICA_wrapperInputSpec
    output_spec = dual_ICA_wrapperOutputSpec

    def _run_interface(self, runtime):
        import os
        import numpy as np
        import nibabel as nb
        import pandas as pd
        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.bold_file).name.rsplit(".nii")

        resampled_priors = resample_IC_file(self.inputs.prior_maps, self.inputs.mask_file)


        brain_mask=np.asarray(nb.load(self.inputs.mask_file).dataobj)
        volume_indices=brain_mask.astype(bool)

        data_array=np.asarray(nb.load(self.inputs.bold_file).dataobj)
        timeseries=np.zeros([data_array.shape[3],volume_indices.sum()])
        for i in range(data_array.shape[3]):
            timeseries[i,:]=(data_array[:,:,:,i])[volume_indices]

        all_IC_array=np.asarray(nb.load(resampled_priors).dataobj)
        all_IC_vectors=np.zeros([all_IC_array.shape[3],volume_indices.sum()])
        for i in range(all_IC_array.shape[3]):
            all_IC_vectors[i,:]=(all_IC_array[:,:,:,i])[volume_indices]

        from rabies.analysis_pkg.prior_modeling import dual_ICA_fit
        prior_fit_out = dual_ICA_fit(timeseries, self.inputs.num_comp, all_IC_vectors, self.inputs.prior_bold_idx)

        dual_ICA_timecourse_csv = os.path.abspath(filename_split[0]+'_dual_ICA_timecourse.csv')
        pd.DataFrame(prior_fit_out['W']).to_csv(dual_ICA_timecourse_csv, header=False, index=False)

        dual_ICA_filename = os.path.abspath(filename_split[0]+'_dual_ICA.nii.gz')
        recover_3D_multiple(self.inputs.mask_file,np.array(prior_fit_out['C'])).to_filename(dual_ICA_filename)

        setattr(self, 'dual_ICA_timecourse_csv', dual_ICA_timecourse_csv)
        setattr(self, 'dual_ICA_filename', dual_ICA_filename)

        return runtime

    def _list_outputs(self):
        return {'dual_ICA_timecourse_csv': getattr(self, 'dual_ICA_timecourse_csv'),
                'dual_ICA_filename': getattr(self, 'dual_ICA_filename')}
