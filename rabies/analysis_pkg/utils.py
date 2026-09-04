def compute_edge_mask(mask_array, num_edge_voxels=1):
    import numpy as np
    #custom function for computing edge mask from an input brain mask
    shape = mask_array.shape

    #iterate through all voxels from the three dimensions and look if it contains surrounding voxels
    edge_mask = np.zeros(shape, dtype=bool)
    num_voxel = 0
    while num_voxel < num_edge_voxels:
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    #only look if the voxel is part of the mask
                    if mask_array[x, y, z]:
                        if (mask_array[x-1:x+2, y-1:y+2, z-1:z+2] == 0).sum() > 0:
                            edge_mask[x, y, z] = 1
        mask_array = mask_array-edge_mask
        num_voxel += 1

    return edge_mask


def load_resample_analysis_maps(mask_file, anat_ref_file,
        transform_list=[], inverse_list=[], seed_dict={}, prior_maps=None, atlas_file=None,
        cleaned_bold_file=None, CR_data_dict=None, remove_EPI_avg=False,
        WM_mask_file=None, CSF_mask_file=None, compute_edge_idx=False,
        interpolation=None, rabies_data_type=None, output_prefix=''):
    '''
    Loads, exactly once, every NIfTI input a scan's FC analyses and/or diagnosis may need,
    as vectorized arrays masked to the brain mask. seed_dict/prior_maps/atlas_file are
    resampled onto the analysis space (defined by mask_file/anat_ref_file) if transform_list
    is non-empty; each is skipped (left None/empty) if not provided, so the same function
    serves callers that only need a subset -- e.g. FC_analysis never passes WM_mask_file/
    CSF_mask_file/compute_edge_idx, which are diagnosis-only.

    Returns a dict with keys: volume_indices, timeseries, seed_arr_dict, prior_map_vectors,
    atlas_idx, roi_list, WM_idx, CSF_idx, edge_idx.
    '''
    import os
    import numpy as np
    import SimpleITK as sitk
    from rabies.utils import resample_volumes, antsApplyTransforms

    if interpolation is None:
        interpolation = sitk.sitkLinear
    if rabies_data_type is None:
        rabies_data_type = sitk.sitkFloat32

    mask_array = sitk.GetArrayFromImage(sitk.ReadImage(mask_file))
    volume_indices = mask_array.astype(bool)

    timeseries = None
    if cleaned_bold_file is not None:
        timeseries = sitk.GetArrayFromImage(sitk.ReadImage(cleaned_bold_file))[:, volume_indices]
        if remove_EPI_avg:
            timeseries = timeseries - CR_data_dict['voxelwise_mean']

    seed_arr_dict = {}
    for seed_name, seed_file in seed_dict.items():
        resampled_seed_file = os.path.abspath(f'{output_prefix}_{seed_name}_resampled.nii.gz')
        antsApplyTransforms(transforms=transform_list, inverses=inverse_list,
            input_image=seed_file, ref_image=anat_ref_file, output_filename=resampled_seed_file,
            interpolation='GenericLabel', rabies_data_type=sitk.sitkInt16, clip_negative=False)
        seed_arr_dict[seed_name] = sitk.GetArrayFromImage(sitk.ReadImage(resampled_seed_file))[volume_indices]

    prior_map_vectors = None
    if prior_maps is not None:
        resampled_4D_img = resample_volumes(in_img=prior_maps, in_ref=anat_ref_file,
            transforms_3d_files=transform_list, inverses_3d=inverse_list, motcorr_params_file=None,
            interpolation=interpolation, rabies_data_type=rabies_data_type, clip_negative=False)
        prior_map_vectors = sitk.GetArrayFromImage(resampled_4D_img)[:, volume_indices]

    atlas_idx = None
    roi_list = None
    if atlas_file is not None:
        resampled_atlas_file = os.path.abspath(f'{output_prefix}_atlas_resampled.nii.gz')
        antsApplyTransforms(transforms=transform_list, inverses=inverse_list,
            input_image=atlas_file, ref_image=anat_ref_file, output_filename=resampled_atlas_file,
            interpolation='GenericLabel', rabies_data_type=sitk.sitkInt16, clip_negative=False)
        atlas_idx = sitk.GetArrayFromImage(sitk.ReadImage(resampled_atlas_file))[volume_indices]

        # original (unresampled) atlas defines which ROI integers are present
        atlas_data = sitk.GetArrayFromImage(sitk.ReadImage(str(atlas_file))).astype(int)
        roi_list = [i for i in range(1, atlas_data.max() + 1) if np.max(i == atlas_data)]

    WM_idx = None
    if WM_mask_file is not None:
        WM_idx = sitk.GetArrayFromImage(sitk.ReadImage(WM_mask_file))[volume_indices].astype(bool)

    CSF_idx = None
    if CSF_mask_file is not None:
        CSF_idx = sitk.GetArrayFromImage(sitk.ReadImage(CSF_mask_file))[volume_indices].astype(bool)

    edge_idx = None
    if compute_edge_idx:
        edge_idx = compute_edge_mask(mask_array, num_edge_voxels=1)[volume_indices]

    return {
        'volume_indices': volume_indices,
        'timeseries': timeseries,
        'seed_arr_dict': seed_arr_dict,
        'prior_map_vectors': prior_map_vectors,
        'atlas_idx': atlas_idx,
        'roi_list': roi_list,
        'WM_idx': WM_idx,
        'CSF_idx': CSF_idx,
        'edge_idx': edge_idx,
    }
