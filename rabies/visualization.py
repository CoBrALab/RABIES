import os
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
# set a dark background
plt.rcParams.update({
    "lines.color": "white",
    "patch.edgecolor": "white",
    "text.color": "black",
    "axes.facecolor": "white",
    "axes.edgecolor": "lightgray",
    "axes.labelcolor": "white",
    "xtick.color": "white",
    "ytick.color": "white",
    "grid.color": "lightgray",
    "figure.facecolor": "black",
    "figure.edgecolor": "black",
    "savefig.facecolor": "black",
    "savefig.edgecolor": "black"})

from rabies.utils import run_command


def otsu_scaling(image_file):
    img = sitk.ReadImage(image_file)
    array = sitk.GetArrayFromImage(img)

    # select a smart vmax for the image display to enhance contrast
    command = f'ThresholdImage 3 {image_file} otsu_weight.nii.gz Otsu 4'
    rc = run_command(command)

    # clip off the background
    mask = sitk.GetArrayFromImage(sitk.ReadImage('otsu_weight.nii.gz'))
    voxel_subset=array[mask>1.0]

    # select a maximal value which encompasses 90% of the voxels in the mask
    voxel_subset.sort()
    vmax=voxel_subset[int(len(voxel_subset)*0.9)]

    scaled = array/vmax
    scaled_img=sitk.GetImageFromArray(scaled, isVector=False)
    scaled_img.CopyInformation(img)
    return scaled_img


def plot_3d(axes,sitk_img,fig,vmin=0,vmax=1,cmap='gray', alpha=1, cbar=False, threshold=None, planes=('sagittal', 'coronal', 'horizontal'), num_slices=4, slice_spacing=0.1):
    physical_dimensions = (np.array(sitk_img.GetSpacing())*np.array(sitk_img.GetSize()))[::-1] # invert because the array is inverted indices
    array=sitk.GetArrayFromImage(sitk_img)

    array[array==0]=None # set 0 values to be empty

    if not threshold is None:
        array[np.abs(array)<threshold]=None

    slice_0 = (1.0-((num_slices-1)*slice_spacing))/2
    slice_fractions=[slice_0]
    for i in range(1,num_slices):
        slice_fractions.append(slice_0+(i*slice_spacing))

    cbar_list = []
    
    ax_number=0
    if 'sagittal' in planes:
        ax=axes[ax_number]
        ax_number+=1
        empty_slice = np.array([np.nan]).repeat(array.shape[0])[:,np.newaxis]
        slices=empty_slice
        for s in slice_fractions:
            slice=array[::-1,:,int(array.shape[2]*s)]
            slices=np.concatenate((slices,slice,empty_slice),axis=1)
        pos = ax.imshow(slices, extent=[0,physical_dimensions[1]*num_slices,0,physical_dimensions[0]], vmin=vmin, vmax=vmax,cmap=cmap, alpha=alpha, interpolation='none')
        ax.axis('off')
        if cbar:
            cbar_list.append(fig.colorbar(pos, ax=ax))

    if 'coronal' in planes:
        ax=axes[ax_number]
        ax_number+=1
        empty_slice = np.array([np.nan]).repeat(array.shape[0])[:,np.newaxis]
        slices=empty_slice
        for s in slice_fractions:
            slice=array[::-1,int(array.shape[1]*s),:]
            slices=np.concatenate((slices,slice,empty_slice),axis=1)
        pos = ax.imshow(slices, extent=[0,physical_dimensions[2]*num_slices,0,physical_dimensions[0]], vmin=vmin, vmax=vmax,cmap=cmap, alpha=alpha, interpolation='none')
        ax.axis('off')
        if cbar:
            cbar_list.append(fig.colorbar(pos, ax=ax))

    if 'horizontal' in planes:
        ax=axes[ax_number]
        ax_number+=1
        empty_slice = np.array([np.nan]).repeat(array.shape[1])[:,np.newaxis]
        slices=empty_slice
        for s in slice_fractions:
            slice=array[int(array.shape[0]*s),::-1,:]
            slices=np.concatenate((slices,slice,empty_slice),axis=1)
        pos = ax.imshow(slices, extent=[0,physical_dimensions[2]*num_slices,0,physical_dimensions[1]], vmin=vmin, vmax=vmax,cmap=cmap, alpha=alpha, interpolation='none')
        ax.axis('off')
        if cbar:
            cbar_list.append(fig.colorbar(pos, ax=ax))
    return cbar_list
