from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
import os
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt

class PlotOverlapInputSpec(BaseInterfaceInputSpec):
    moving = File(exists=True, mandatory=True,
                  desc="Moving image from registration.")
    fixed = File(exists=True, mandatory=True,
                 desc="Fixed image from registration.")
    out_dir = traits.Str(mandatory=True, desc="Directory for QC outputs.")
    name_source = traits.Str(mandatory=True, desc="Input file template for naming outputs.")


class PlotOverlapOutputSpec(TraitedSpec):
    out_png = File(exists=True, desc="Output png.")


class PlotOverlap(BaseInterface):

    input_spec = PlotOverlapInputSpec
    output_spec = PlotOverlapOutputSpec

    def _run_interface(self, runtime):
        import os
        import pathlib
        filename_template = pathlib.Path(self.inputs.name_source).name.rsplit(".nii")[0]

        script_path = 'plot_overlap.sh'
        os.makedirs(self.inputs.out_dir, exist_ok=True)
        out_name = self.inputs.out_dir+'/' + \
            filename_template+'_registration.png'

        from rabies.preprocess_pkg.utils import run_command
        command = 'bash %s %s %s %s' % (
            script_path, self.inputs.moving, self.inputs.fixed, out_name)
        rc = run_command(command)

        setattr(self, 'out_png', out_name)
        return runtime

    def _list_outputs(self):
        return {'out_png': getattr(self, 'out_png')}


def otsu_scaling(image):
    import numpy as np
    import SimpleITK as sitk
    img = sitk.ReadImage(image)
    array = sitk.GetArrayFromImage(img)

    # select a smart vmax for the image display to enhance contrast
    from rabies.preprocess_pkg.utils import run_command
    command = 'ThresholdImage 3 %s otsu_weight.nii.gz Otsu 4' % (image)
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

    ax_number=0
    if 'sagittal' in planes:
        ax=axes[ax_number]
        ax_number+=1
        slices=np.empty([array.shape[0],1])
        for s in slice_fractions:
            slice=array[::-1,:,int(array.shape[2]*s)]
            slices=np.concatenate((slices,slice,np.empty([array.shape[0],1])),axis=1)
        pos = ax.imshow(slices, extent=[0,physical_dimensions[1]*num_slices,0,physical_dimensions[0]], vmin=vmin, vmax=vmax,cmap=cmap, alpha=alpha, interpolation='none')
        ax.axis('off')
        if cbar:
            fig.colorbar(pos, ax=ax)

    if 'coronal' in planes:
        ax=axes[ax_number]
        ax_number+=1
        slices=np.empty([array.shape[0],1])
        for s in slice_fractions:
            slice=array[::-1,int(array.shape[1]*s),:]
            slices=np.concatenate((slices,slice,np.empty([array.shape[0],1])),axis=1)
        pos = ax.imshow(slices, extent=[0,physical_dimensions[2]*num_slices,0,physical_dimensions[0]], vmin=vmin, vmax=vmax,cmap=cmap, alpha=alpha, interpolation='none')
        ax.axis('off')
        if cbar:
            fig.colorbar(pos, ax=ax)

    if 'horizontal' in planes:
        ax=axes[ax_number]
        ax_number+=1
        slices=np.empty([array.shape[1],1])
        for s in slice_fractions:
            slice=array[int(array.shape[0]*s),::-1,:]
            slices=np.concatenate((slices,slice,np.empty([array.shape[1],1])),axis=1)
        pos = ax.imshow(slices, extent=[0,physical_dimensions[2]*num_slices,0,physical_dimensions[1]], vmin=vmin, vmax=vmax,cmap=cmap, alpha=alpha, interpolation='none')
        ax.axis('off')
        if cbar:
            fig.colorbar(pos, ax=ax)

def plot_reg(image1,image2, name_source, out_dir):
    import os
    import pathlib
    filename_template = pathlib.Path(name_source).name.rsplit(".nii")[0]
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_dir+'/'+ \
        filename_template

    import matplotlib.pyplot as plt
    from rabies.preprocess_pkg.preprocess_visual_QC import plot_3d,otsu_scaling
    fig,axes = plt.subplots(nrows=2, ncols=3, figsize=(12*3,2*2))
    plt.tight_layout()

    scaled = otsu_scaling(image1)
    display1,display2,display3 = plot_3d(scaled,axes[0,:], cmap='gray')
    display1.add_edges(image2)
    display2.add_edges(image2)
    display3.add_edges(image2)

    scaled = otsu_scaling(image2)
    display1,display2,display3 = plot_3d(scaled,axes[1,:], cmap='gray')
    display1.add_edges(image1)
    display2.add_edges(image1)
    display3.add_edges(image1)
    fig.savefig('%s_registration.png' % (prefix), bbox_inches='tight')


def template_info(anat_template, opts, out_dir):
    import os
    from nilearn import plotting
    import matplotlib.pyplot as plt
    from rabies.preprocess_pkg.preprocess_visual_QC import plot_3d,otsu_scaling
    brain_mask = str(opts.brain_mask)
    WM_mask = str(opts.WM_mask)
    CSF_mask = str(opts.CSF_mask)
    vascular_mask = str(opts.vascular_mask)
    labels = str(opts.labels)
    os.makedirs(out_dir, exist_ok=True)

    import SimpleITK as sitk
    # make sure that masks are binary
    for mask in [brain_mask,WM_mask,vascular_mask]:
        img = sitk.ReadImage(mask)
        array = sitk.GetArrayFromImage(img)
        if ((array!=1)*(array!=0)).sum()>0:
            raise ValueError("The file %s is not a binary mask. Non-binary masks cannot be processed." % (mask))

    scaled = otsu_scaling(anat_template)

    fig,axes = plt.subplots(nrows=3, ncols=6, figsize=(4*6,2*2))

    axes[0,0].set_title('Anatomical Template', fontsize=20)
    plot_3d(axes[:,0],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    # plot brain mask
    mask = brain_mask
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    axes[0,1].set_title('Brain Mask', fontsize=20)
    plot_3d(axes[:,1],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,1],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.3, cbar=False)
    # plot WM mask
    mask = WM_mask
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    axes[0,2].set_title('WM Mask', fontsize=20)
    plot_3d(axes[:,2],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,2],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.5, cbar=False)
    # plot CSF mask
    mask = CSF_mask
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    axes[0,3].set_title('CSF Mask', fontsize=20)
    plot_3d(axes[:,3],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,3],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.5, cbar=False)
    # plot VASC mask
    mask = vascular_mask
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    axes[0,4].set_title('Vascular Mask', fontsize=20)
    plot_3d(axes[:,4],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,4],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.5, cbar=False)

    # plot labels
    mask = labels
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    axes[0,5].set_title('Atlas Labels', fontsize=20)
    plot_3d(axes[:,5],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,5],sitk_mask,fig=fig,vmin=1,vmax=sitk.GetArrayFromImage(sitk_mask).max(),cmap='rainbow', alpha=0.5, cbar=False)
    plt.tight_layout()

    fig.savefig(out_dir+'/template_files.png', bbox_inches='tight')


def temporal_features(bold_file, confounds_csv, FD_csv, rabies_data_type, name_source, out_dir):
    import os
    import pathlib
    filename_template = pathlib.Path(name_source).name.rsplit(".nii")[0]
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_dir+'/'+ \
        filename_template

    import numpy as np
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    from rabies.preprocess_pkg.preprocess_visual_QC import plot_3d
    from rabies.preprocess_pkg.utils import copyInfo_3DImage
    fig,axes = plt.subplots(nrows=3, ncols=3, figsize=(20,5))
    # plot the motion timecourses
    import pandas as pd
    df = pd.read_csv(confounds_csv)
    ax = axes[0,0]
    ax.plot(df['mov1'])
    ax.plot(df['mov2'])
    ax.plot(df['mov3'])
    ax.legend(['mov1','mov2','mov3'])
    ax.set_title('Translation parameters', fontsize=20)
    ax = axes[1,0]
    ax.plot(df['rot1'])
    ax.plot(df['rot2'])
    ax.plot(df['rot3'])
    ax.legend(['rot1','rot2','rot3'])
    ax.set_title('Rotation parameters', fontsize=20)

    df = pd.read_csv(FD_csv)
    ax=axes[2,0]
    ax.plot(df['Mean'], color='r')
    ax.set_title('Framewise Displacement', fontsize=20)

    plt.tight_layout()

    # calculate STD and tSNR map on preprocessed timeseries
    img = sitk.ReadImage(bold_file, rabies_data_type)
    array = sitk.GetArrayFromImage(img)
    mean = array.mean(axis=0)
    std = array.std(axis=0)
    std_filename = os.path.abspath('tSTD.nii.gz')
    std_image = copyInfo_3DImage(
        sitk.GetImageFromArray(std, isVector=False), img)
    sitk.WriteImage(std_image, std_filename)

    tSNR = np.divide(mean, std)
    tSNR[np.isnan(tSNR)]=0
    tSNR_filename = os.path.abspath('tSNR.nii.gz')
    tSNR_image = copyInfo_3DImage(
        sitk.GetImageFromArray(tSNR, isVector=False), img)
    sitk.WriteImage(tSNR_image, tSNR_filename)

    axes[0,1].set_title('Temporal STD', fontsize=20)
    std=std.flatten()
    std.sort()
    std_vmax = std[int(len(std)*0.95)]
    plot_3d(axes[:,1],std_image,fig=fig,vmin=0,vmax=std_vmax,cmap='inferno', cbar=True)
    axes[0,2].set_title('Temporal SNR', fontsize=20)
    plot_3d(axes[:,2],tSNR_image,fig=fig,vmin=0,vmax=tSNR.max(),cmap='Spectral', cbar=True)

    fig.savefig('%s_temporal_features.png' % (prefix), bbox_inches='tight')

    return std_filename, tSNR_filename


def denoising_diagnosis(raw_img,init_denoise,warped_mask,final_denoise, name_source, out_dir):
    import os
    import pathlib
    import SimpleITK as sitk
    filename_template = pathlib.Path(name_source).name.rsplit(".nii")[0]
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_dir+'/'+ \
        filename_template

    import matplotlib.pyplot as plt
    from rabies.preprocess_pkg.preprocess_visual_QC import plot_3d,otsu_scaling
    fig,axes = plt.subplots(nrows=3, ncols=4, figsize=(12*4,2*3))
    plt.tight_layout()

    scaled = otsu_scaling(raw_img)
    axes[0,0].set_title('Raw EPI', fontsize=20)
    plot_3d(axes[:,0],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')

    scaled = otsu_scaling(init_denoise)
    axes[0,1].set_title('Initial Denoising', fontsize=20)
    plot_3d(axes[:,1],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')

    axes[0,2].set_title('Resampled Mask', fontsize=20)
    plot_3d(axes[:,2],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')
    plot_3d(axes[:,2],sitk.ReadImage(warped_mask,sitk.sitkFloat32),fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.3, cbar=False)

    scaled = otsu_scaling(final_denoise)
    axes[0,3].set_title('Final Denoising', fontsize=20)
    plot_3d(axes[:,3],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')

    fig.savefig('%s_denoising.png' % (prefix), bbox_inches='tight')
