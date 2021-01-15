from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


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

        import rabies
        dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
        script_path = dir_path+'/shell_scripts/plot_overlap.sh'
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
    import nibabel as nb
    img=nb.load(image)
    array=np.asarray(img.dataobj)

    # select a smart vmax for the image display to enhance contrast
    from rabies.preprocess_pkg.utils import run_command
    command = 'ThresholdImage 3 %s otsu_weight.nii.gz Otsu 4' % (image)
    rc = run_command(command)

    # clip off the background
    mask = np.asarray(nb.load('otsu_weight.nii.gz').dataobj)
    voxel_subset=array[mask>1.0]

    # select a maximal value which encompasses 90% of the voxels in the mask
    voxel_subset.sort()
    vmax=voxel_subset[int(len(voxel_subset)*0.9)]

    # re-scale the values to be within a range of -1 and 1
    scaled = ((array/vmax)*2)-1
    return nb.Nifti1Image(scaled, img.affine, img.header)

def plot_3d(image,axes,vmax=1,cmap='gray', cbar=False):
    from nilearn import plotting
    ax=axes[0]
    display1 = plotting.plot_stat_map(image,bg_img=image, axes=ax, cmap=cmap, cut_coords=4, display_mode='x', vmax=vmax, threshold=None, draw_cross=False, colorbar=cbar)
    ax=axes[1]
    display2 = plotting.plot_stat_map(image,bg_img=image, axes=ax, cmap=cmap, cut_coords=4, display_mode='y', vmax=vmax, threshold=None, draw_cross=False, colorbar=cbar)
    ax=axes[2]
    display3 = plotting.plot_stat_map(image,bg_img=image, axes=ax, cmap=cmap, cut_coords=4, display_mode='z', vmax=vmax, threshold=None, draw_cross=False, colorbar=cbar)
    return display1,display2,display3

def plot_reg(image1,image2, name_source, out_dir):
    import os
    import pathlib
    filename_template = pathlib.Path(name_source).name.rsplit(".nii")[0]
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_dir+'/'+ \
        filename_template

    import matplotlib.pyplot as plt
    from rabies.preprocess_pkg.visual_diagnosis import plot_3d,otsu_scaling
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

def template_diagnosis(anat_template, opts, out_dir):
    import os
    from nilearn import plotting
    import matplotlib.pyplot as plt
    from rabies.preprocess_pkg.visual_diagnosis import plot_3d,otsu_scaling
    brain_mask = str(opts.brain_mask)
    WM_mask = str(opts.WM_mask)
    CSF_mask = str(opts.CSF_mask)
    vascular_mask = str(opts.vascular_mask)
    labels = str(opts.labels)
    os.makedirs(out_dir, exist_ok=True)

    scaled = otsu_scaling(anat_template)

    fig,axes = plt.subplots(nrows=6, ncols=3, figsize=(12*3,2*6))
    plt.tight_layout()

    display1,display2,display3 = plot_3d(scaled,axes[0,:], cmap='gray')
    # plot brain mask
    mask = brain_mask
    display1,display2,display3 = plot_3d(scaled,axes[1,:], cmap='gray')
    display1.add_overlay(mask, cmap=plotting.cm.red_transparent)
    display2.add_overlay(mask, cmap=plotting.cm.red_transparent)
    display3.add_overlay(mask, cmap=plotting.cm.red_transparent)
    # plot WM mask
    mask = WM_mask
    display1,display2,display3 = plot_3d(scaled,axes[2,:], cmap='gray')
    display1.add_overlay(mask, cmap=plotting.cm.red_transparent)
    display2.add_overlay(mask, cmap=plotting.cm.red_transparent)
    display3.add_overlay(mask, cmap=plotting.cm.red_transparent)
    # plot CSF mask
    mask = CSF_mask
    display1,display2,display3 = plot_3d(scaled,axes[3,:], cmap='gray')
    display1.add_overlay(mask, cmap=plotting.cm.red_transparent)
    display2.add_overlay(mask, cmap=plotting.cm.red_transparent)
    display3.add_overlay(mask, cmap=plotting.cm.red_transparent)
    # plot VASC mask
    mask = vascular_mask
    display1,display2,display3 = plot_3d(scaled,axes[4,:], cmap='gray')
    display1.add_overlay(mask, cmap=plotting.cm.red_transparent)
    display2.add_overlay(mask, cmap=plotting.cm.red_transparent)
    display3.add_overlay(mask, cmap=plotting.cm.red_transparent)

    # plot labels
    mask = labels
    display1,display2,display3 = plot_3d(scaled,axes[5,:], cmap='gray')
    display1.add_overlay(mask, cmap='rainbow')
    display2.add_overlay(mask, cmap='rainbow')
    display3.add_overlay(mask, cmap='rainbow')
    fig.savefig(out_dir+'/template_diagnosis.png', bbox_inches='tight')

def temporal_diagnosis(bold_file, confounds_csv, FD_csv, rabies_data_type, name_source, out_dir):
    import os
    import pathlib
    filename_template = pathlib.Path(name_source).name.rsplit(".nii")[0]
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_dir+'/'+ \
        filename_template

    import numpy as np
    import SimpleITK as sitk
    from nilearn import plotting
    import matplotlib.pyplot as plt
    from rabies.preprocess_pkg.visual_diagnosis import plot_3d,otsu_scaling
    from rabies.preprocess_pkg.utils import copyInfo_3DImage
    fig,axes = plt.subplots(nrows=3, ncols=3, figsize=(12*3,3*3))
    # plot the motion timecourses
    import pandas as pd
    df = pd.read_csv(confounds_csv)
    ax = axes[0,0]
    ax.plot(df['mov1'])
    ax.plot(df['mov2'])
    ax.plot(df['mov3'])
    ax.legend(['mov1','mov2','mov3'])
    ax.set_title('Translation parameters', fontsize=20)
    ax = axes[0,1]
    ax.plot(df['rot1'])
    ax.plot(df['rot2'])
    ax.plot(df['rot3'])
    ax.legend(['rot1','rot2','rot3'])
    ax.set_title('Rotation parameters', fontsize=20)

    df = pd.read_csv(FD_csv)
    ax=axes[0,2]
    ax.plot(df['Mean'], color='r')
    ax.set_title('Framewise Displacement', fontsize=20)

    plt.tight_layout()

    # calculate STD and tSNR map on preprocessed timeseries
    img = sitk.ReadImage(bold_file, rabies_data_type)
    array = sitk.GetArrayFromImage(img)
    mean = array.mean(axis=0)
    std = array.std(axis=0)
    std_filename = os.path.abspath('tSTD.nii.gz')
    image_3d = copyInfo_3DImage(
        sitk.GetImageFromArray(std, isVector=False), img)
    sitk.WriteImage(image_3d, std_filename)

    tSNR = np.divide(mean, std)
    tSNR_filename = os.path.abspath('tSNR.nii.gz')
    image_3d = copyInfo_3DImage(
        sitk.GetImageFromArray(tSNR, isVector=False), img)
    sitk.WriteImage(image_3d, tSNR_filename)

    plot_3d(std_filename,axes[1,:],vmax=std.max(),cmap=plotting.cm.cold_hot, cbar=True)
    plot_3d(tSNR_filename,axes[2,:],vmax=tSNR.max(),cmap='Spectral', cbar=True)

    fig.savefig('%s_temporal_diagnosis.png' % (prefix), bbox_inches='tight')

    return std_filename, tSNR_filename


def denoising_diagnosis(raw_img,init_denoise,warped_mask,final_denoise, name_source, out_dir):
    import os
    import pathlib
    filename_template = pathlib.Path(name_source).name.rsplit(".nii")[0]
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_dir+'/'+ \
        filename_template

    from nilearn import plotting
    import matplotlib.pyplot as plt
    from rabies.preprocess_pkg.visual_diagnosis import plot_3d,otsu_scaling
    fig,axes = plt.subplots(nrows=4, ncols=3, figsize=(12*3,2*4))
    plt.tight_layout()

    scaled = otsu_scaling(raw_img)
    display1,display2,display3 = plot_3d(scaled,axes[0,:], cmap='viridis')
    display1,display2,display3 = plot_3d(scaled,axes[2,:], cmap='viridis')
    display1.add_overlay(warped_mask, cmap=plotting.cm.red_transparent)
    display2.add_overlay(warped_mask, cmap=plotting.cm.red_transparent)
    display3.add_overlay(warped_mask, cmap=plotting.cm.red_transparent)

    scaled = otsu_scaling(init_denoise)
    display1,display2,display3 = plot_3d(scaled,axes[1,:], cmap='viridis')

    scaled = otsu_scaling(final_denoise)
    display1,display2,display3 = plot_3d(scaled,axes[3,:], cmap='viridis')

    fig.savefig('%s_denoising.png' % (prefix), bbox_inches='tight')
