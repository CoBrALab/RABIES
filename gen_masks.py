import sys
import os

atlas_file = os.path.abspath(sys.argv[1])
csv_labels = os.path.abspath(sys.argv[2])
prefix = str(sys.argv[3])


def compute_masks(atlas, csv_labels, prefix):
    import numpy as np
    import pandas as pd
    import os
    import SimpleITK as sitk

    '''extract rows from the csv file'''
    label_info = pd.read_csv(csv_labels)
    right_id = label_info['right label']
    left_id = label_info['left label']
    tissue_type = label_info['tissue type']

    '''create lists with the WM/CSF label numbers'''
    WM_labels = []
    CSF_labels = []

    for x in range(1, len(tissue_type)):
        if tissue_type[x] == 'WM':
            WM_labels.append(right_id[x])
            WM_labels.append(left_id[x])
        elif tissue_type[x] == 'CSF':
            CSF_labels.append(right_id[x])
            CSF_labels.append(left_id[x])

    # generate a list of the atlas labels which correspond to WM or CSF regions
    WM_labels_list = np.zeros(len(WM_labels))
    CSF_labels_list = np.zeros(len(CSF_labels))

    for x in range(0, len(WM_labels)):
        WM_labels_list[x] = int(WM_labels[x])

    for x in range(0, len(CSF_labels)):
        CSF_labels_list[x] = int(CSF_labels[x])

    '''extract the voxels which fall into the labels'''
    atlas_img = sitk.ReadImage(os.path.abspath(atlas), sitk.sitkInt32)
    atlas_data = sitk.GetArrayFromImage(atlas_img)
    shape = atlas_data.shape
    vector_atlas = np.transpose(np.reshape(
        atlas_data, [shape[0]*shape[1]*shape[2]]))
    WM_voxels = np.zeros(vector_atlas.shape)
    CSF_voxels = np.zeros(vector_atlas.shape)

    # for each WM/CSF label, establish which voxels are part of these labels
    for x in range(0, len(WM_labels_list)):
        WM_voxels = WM_voxels+(vector_atlas == WM_labels_list[x])

    WM_voxels = (WM_voxels >= 1).astype(int)  # binarize the mask

    for x in range(0, len(CSF_labels_list)):
        CSF_voxels = CSF_voxels+(vector_atlas == CSF_labels_list[x])

    CSF_voxels = (CSF_voxels >= 1).astype(int)  # binarize the mask

    WM_mask = np.reshape(WM_voxels, shape)
    CSF_mask = np.reshape(CSF_voxels, shape)

    WM_mask_file = '%s_WM_mask.nii.gz' % (prefix,)
    CSF_mask_file = '%s_CSF_mask.nii.gz' % (prefix,)

    WM_mask_img = sitk.GetImageFromArray(
        WM_mask.astype('int16'), isVector=False)
    WM_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(WM_mask_img, WM_mask_file)
    CSF_mask_img = sitk.GetImageFromArray(
        CSF_mask.astype('int16'), isVector=False)
    CSF_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(CSF_mask_img, CSF_mask_file)

    '''Erode the masks'''
    from scipy.ndimage.morphology import binary_erosion
    eroded_WM_mask = binary_erosion(WM_mask, iterations=1)
    eroded_CSF_mask = binary_erosion(CSF_mask, iterations=1)

    eroded_WM_mask_file = '%s_eroded_WM_mask.nii.gz' % (prefix,)
    eroded_CSF_mask_file = '%s_eroded_CSF_mask.nii.gz' % (prefix,)

    eroded_WM_mask_img = sitk.GetImageFromArray(
        eroded_WM_mask.astype('int16'), isVector=False)
    eroded_WM_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(eroded_WM_mask_img, eroded_WM_mask_file)
    eroded_CSF_mask_img = sitk.GetImageFromArray(
        eroded_CSF_mask.astype('int16'), isVector=False)
    eroded_CSF_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(eroded_CSF_mask_img, eroded_CSF_mask_file)

    return WM_mask_file, CSF_mask_file, eroded_WM_mask_file, eroded_CSF_mask_file


out = compute_masks(atlas_file, csv_labels, prefix)
