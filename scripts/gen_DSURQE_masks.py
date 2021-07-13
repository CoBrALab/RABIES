#! /usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk

atlas_file = os.path.abspath(sys.argv[1])
csv_labels = os.path.abspath(sys.argv[2])
prefix = str(sys.argv[3])


def compute_masks(atlas, csv_labels, prefix):
    df = pd.read_csv(csv_labels)

    '''create lists with the region label numbers'''

    GM_right_labels = df['right label'][np.array(df['tissue type'] == 'GM')]
    GM_left_labels = df['left label'][np.array(df['tissue type'] == 'GM')]

    # take only labels that are specific to one axis to avoid overlap
    right_hem_labels = list(GM_right_labels[GM_right_labels != GM_left_labels])
    left_hem_labels = list(GM_left_labels[GM_right_labels != GM_left_labels])

    GM_labels = \
        list(df['right label'][np.array(df['tissue type'] == 'GM')]) + \
        list(df['left label'][np.array(df['tissue type'] == 'GM')])

    WM_labels = \
        list(df['right label'][np.array(df['tissue type'] == 'WM')]) + \
        list(df['left label'][np.array(df['tissue type'] == 'WM')])

    CSF_labels = \
        list(df['right label'][np.array(df['tissue type'] == 'CSF')]) + \
        list(df['left label'][np.array(df['tissue type'] == 'CSF')])

    '''extract the voxels which fall into the labels'''
    atlas_img = sitk.ReadImage(os.path.abspath(atlas), sitk.sitkInt32)
    atlas_data = sitk.GetArrayFromImage(atlas_img)
    shape = atlas_data.shape

    right_hem_mask = np.zeros(shape).astype(bool)
    left_hem_mask = np.zeros(shape).astype(bool)
    GM_mask = np.zeros(shape).astype(bool)
    WM_mask = np.zeros(shape).astype(bool)
    CSF_mask = np.zeros(shape).astype(bool)

    for i in range(atlas_data.max()+1):
        roi_mask = atlas_data == i
        if i in right_hem_labels:
            right_hem_mask += roi_mask
        if i in left_hem_labels:
            left_hem_mask += roi_mask
        if i in GM_labels:
            GM_mask += roi_mask
        if i in WM_labels:
            WM_mask += roi_mask
        if i in CSF_labels:
            CSF_mask += roi_mask

    GM_mask_file = f'{prefix}_GM_mask.nii.gz'
    WM_mask_file = f'{prefix}_WM_mask.nii.gz'
    CSF_mask_file = f'{prefix}_CSF_mask.nii.gz'
    right_hem_mask_file = f'{prefix}_right_hem_mask.nii.gz'
    left_hem_mask_file = f'{prefix}_left_hem_mask.nii.gz'

    GM_mask_img = sitk.GetImageFromArray(
        GM_mask.astype('int16'), isVector=False)
    GM_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(GM_mask_img, GM_mask_file)

    WM_mask_img = sitk.GetImageFromArray(
        WM_mask.astype('int16'), isVector=False)
    WM_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(WM_mask_img, WM_mask_file)

    CSF_mask_img = sitk.GetImageFromArray(
        CSF_mask.astype('int16'), isVector=False)
    CSF_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(CSF_mask_img, CSF_mask_file)

    right_hem_mask_img = sitk.GetImageFromArray(
        right_hem_mask.astype('int16'), isVector=False)
    right_hem_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(right_hem_mask_img, right_hem_mask_file)

    left_hem_mask_img = sitk.GetImageFromArray(
        left_hem_mask.astype('int16'), isVector=False)
    left_hem_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(left_hem_mask_img, left_hem_mask_file)

    '''Erode the masks'''
    from scipy.ndimage.morphology import binary_erosion
    eroded_WM_mask = binary_erosion(WM_mask, iterations=1)
    eroded_CSF_mask = binary_erosion(CSF_mask, iterations=1)

    eroded_WM_mask_file = f'{prefix}_eroded_WM_mask.nii.gz'
    eroded_CSF_mask_file = f'{prefix}_eroded_CSF_mask.nii.gz'

    eroded_WM_mask_img = sitk.GetImageFromArray(
        eroded_WM_mask.astype('int16'), isVector=False)
    eroded_WM_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(eroded_WM_mask_img, eroded_WM_mask_file)
    eroded_CSF_mask_img = sitk.GetImageFromArray(
        eroded_CSF_mask.astype('int16'), isVector=False)
    eroded_CSF_mask_img.CopyInformation(atlas_img)
    sitk.WriteImage(eroded_CSF_mask_img, eroded_CSF_mask_file)

    return [GM_mask_file, WM_mask_file, CSF_mask_file, right_hem_mask_file,
            left_hem_mask_file, eroded_WM_mask_file, eroded_CSF_mask_file]


out = compute_masks(atlas_file, csv_labels, prefix)
