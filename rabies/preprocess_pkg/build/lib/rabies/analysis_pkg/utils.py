import numpy as np
import SimpleITK as sitk


def resample_prior_maps(in_file, ref_file, transforms = [], inverses = []):
    # resampling the reference image to the dimension of the EPI
    import SimpleITK as sitk
    import os
    from rabies.utils import split_volumes, exec_applyTransforms

    # Splitting bold file into lists of single volumes
    [volumes_list, num_volumes] = split_volumes(
        in_file, "bold_", sitk.sitkFloat32)

    warped_volumes = []
    for x in range(0, num_volumes):
        warped_vol_fname = os.path.abspath(
            "deformed_volume" + str(x) + ".nii.gz")
        warped_volumes.append(warped_vol_fname)
        exec_applyTransforms(transforms=transforms, inverses=inverses, input_image=volumes_list[x], ref_image=ref_file, output_image=warped_vol_fname, interpolation='Linear')

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

    return combined


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