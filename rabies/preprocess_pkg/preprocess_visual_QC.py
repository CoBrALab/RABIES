from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
import os
import pathlib
import matplotlib.pyplot as plt
from rabies.utils import run_command

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
        filename_template = pathlib.Path(self.inputs.name_source).name.rsplit(".nii")[0]

        script_path = 'plot_overlap.sh'
        os.makedirs(self.inputs.out_dir, exist_ok=True)
        out_name = self.inputs.out_dir+'/' + \
            filename_template+'_registration.png'

        command = f'{script_path} {self.inputs.moving} {self.inputs.fixed} {out_name}'
        rc = run_command(command)

        setattr(self, 'out_png', out_name)
        return runtime

    def _list_outputs(self):
        return {'out_png': getattr(self, 'out_png')}


def plot_reg(image1,image2, name_source, out_dir):
    import os
    import pathlib
    filename_template = pathlib.Path(name_source).name.rsplit(".nii")[0]
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_dir+'/'+ \
        filename_template

    import matplotlib.pyplot as plt
    from rabies.visualization import plot_3d, otsu_scaling
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
    fig.savefig(f'{prefix}_registration.png', bbox_inches='tight')


def template_info(anat_template, opts, out_dir):
    import os
    import SimpleITK as sitk
    # set default threader to platform to avoid freezing with MultiProc https://github.com/SimpleITK/SimpleITK/issues/1239
    sitk.ProcessObject_SetGlobalDefaultThreader('Platform')
    from nilearn import plotting
    import matplotlib.pyplot as plt
    from rabies.visualization import plot_3d, otsu_scaling
    brain_mask = str(opts.brain_mask)
    WM_mask = str(opts.WM_mask)
    CSF_mask = str(opts.CSF_mask)
    vascular_mask = str(opts.vascular_mask)
    labels = str(opts.labels)
    os.makedirs(out_dir, exist_ok=True)

    scaled = otsu_scaling(anat_template)

    fig,axes = plt.subplots(nrows=3, ncols=6, figsize=(4*6,2*2))

    axes[0,0].set_title('Anatomical Template', fontsize=30, color='white')
    plot_3d(axes[:,0],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    # plot brain mask
    mask = brain_mask
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    # resample mask to match template
    sitk_mask = sitk.Resample(sitk_mask, scaled)
    axes[0,1].set_title('Brain Mask', fontsize=30, color='white')
    plot_3d(axes[:,1],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,1],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.3, cbar=False)
    # plot WM mask
    mask = WM_mask
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    # resample mask to match template
    sitk_mask = sitk.Resample(sitk_mask, scaled)
    axes[0,2].set_title('WM Mask', fontsize=30, color='white')
    plot_3d(axes[:,2],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,2],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.5, cbar=False)
    # plot CSF mask
    mask = CSF_mask
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    # resample mask to match template
    sitk_mask = sitk.Resample(sitk_mask, scaled)
    axes[0,3].set_title('CSF Mask', fontsize=30, color='white')
    plot_3d(axes[:,3],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,3],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.5, cbar=False)
    # plot VASC mask
    mask = vascular_mask
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    # resample mask to match template
    sitk_mask = sitk.Resample(sitk_mask, scaled)
    axes[0,4].set_title('Vascular Mask', fontsize=30, color='white')
    plot_3d(axes[:,4],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,4],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.5, cbar=False)

    # plot labels
    mask = labels
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    # resample mask to match template
    sitk_mask = sitk.Resample(sitk_mask, scaled)
    axes[0,5].set_title('Atlas Labels', fontsize=30, color='white')
    plot_3d(axes[:,5],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:,5],sitk_mask,fig=fig,vmin=1,vmax=sitk.GetArrayFromImage(sitk_mask).max(),cmap='rainbow', alpha=0.5, cbar=False)
    plt.tight_layout()

    fig.savefig(out_dir+'/template_files.png', bbox_inches='tight')


def template_masking(template, mask, out_dir):
    import os
    import SimpleITK as sitk
    # set default threader to platform to avoid freezing with MultiProc https://github.com/SimpleITK/SimpleITK/issues/1239
    sitk.ProcessObject_SetGlobalDefaultThreader('Platform')
    from nilearn import plotting
    import matplotlib.pyplot as plt
    from rabies.visualization import plot_3d, otsu_scaling

    os.makedirs(out_dir, exist_ok=True)

    scaled = otsu_scaling(template)

    fig,axes = plt.subplots(nrows=3, ncols=1, figsize=(4,2*2))

    # plot brain mask
    sitk_mask = sitk.ReadImage(
        mask, sitk.sitkFloat32)
    # resample mask to match template
    sitk_mask = sitk.Resample(sitk_mask, scaled)
    plot_3d(axes[:],scaled,fig=fig,vmin=0,vmax=1,cmap='gray')
    plot_3d(axes[:],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.3, cbar=False)
    plt.tight_layout()
    fig.savefig(out_dir+'/template_masking.png', bbox_inches='tight')

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
    from rabies.visualization import plot_3d
    from rabies.utils import copyInfo_3DImage
    fig,axes = plt.subplots(nrows=3, ncols=3, figsize=(20,5))
    # plot the motion timecourses
    import pandas as pd
    df = pd.read_csv(confounds_csv)
    ax = axes[0,0]
    ax.plot(df['mov1'])
    ax.plot(df['mov2'])
    ax.plot(df['mov3'])
    ax.legend(['mov1','mov2','mov3'])
    ax.set_title('Translation parameters', fontsize=30, color='white')
    ax = axes[1,0]
    ax.plot(df['rot1'])
    ax.plot(df['rot2'])
    ax.plot(df['rot3'])
    ax.legend(['rot1','rot2','rot3'])
    ax.set_title('Rotation parameters', fontsize=30, color='white')

    df = pd.read_csv(FD_csv)
    ax=axes[2,0]
    ax.plot(df['Mean'], color='r')
    ax.set_title('Framewise Displacement', fontsize=30, color='white')

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

    axes[0,1].set_title('Temporal STD', fontsize=30, color='white')
    std=std.flatten()
    std.sort()
    std_vmax = std[int(len(std)*0.95)]
    plot_3d(axes[:,1],std_image,fig=fig,vmin=0,vmax=std_vmax,cmap='inferno', cbar=True)
    axes[0,2].set_title('Temporal SNR', fontsize=30, color='white')
    plot_3d(axes[:,2],tSNR_image,fig=fig,vmin=0,vmax=tSNR.max(),cmap='Spectral', cbar=True)

    fig.savefig(f'{prefix}_temporal_features.png', bbox_inches='tight')

    return std_filename, tSNR_filename


def inho_cor_diagnosis(raw_img,init_denoise,warped_mask,final_denoise, name_source, out_dir):
    import os
    import pathlib
    import SimpleITK as sitk
    # set default threader to platform to avoid freezing with MultiProc https://github.com/SimpleITK/SimpleITK/issues/1239
    sitk.ProcessObject_SetGlobalDefaultThreader('Platform')
    filename_template = pathlib.Path(name_source).name.rsplit(".nii")[0]
    os.makedirs(out_dir, exist_ok=True)
    prefix = out_dir+'/'+ \
        filename_template

    import matplotlib.pyplot as plt
    from rabies.visualization import plot_3d, otsu_scaling
    fig,axes = plt.subplots(nrows=3, ncols=4, figsize=(12*4,2*3))

    scaled = otsu_scaling(raw_img)
    axes[0,0].set_title('Raw Image', fontsize=30, color='white')
    #add_filenames(axes[-1,0], {'File':raw_img})
    plot_3d(axes[:,0],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')

    axes[0,2].set_title('Resampled Mask', fontsize=30, color='white')
    #add_filenames(axes[-1,2], {'Mask File':warped_mask,'EPI File':raw_img})
    plot_3d(axes[:,2],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')
    sitk_mask = sitk.ReadImage(warped_mask,sitk.sitkFloat32)
    # resample mask to match template
    sitk_mask = sitk.Resample(sitk_mask, scaled)
    plot_3d(axes[:,2],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.3, cbar=False)

    scaled = otsu_scaling(init_denoise)
    axes[0,1].set_title('Initial Correction', fontsize=30, color='white')
    #add_filenames(axes[-1,1], {'File':init_denoise})
    plot_3d(axes[:,1],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')

    scaled = otsu_scaling(final_denoise)
    axes[0,3].set_title('Final Correction', fontsize=30, color='white')
    #add_filenames(axes[-1,3], {'File':final_denoise})
    plot_3d(axes[:,3],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')

    plt.tight_layout()
    fig.savefig(f'{prefix}_inho_cor.png', bbox_inches='tight')

def add_filenames(ax, file_dict, line_length=40):
    txt=""
    for key in list(file_dict.keys()):
        txt+=key+": "
        file=file_dict[key]
        i=0
        while(i<len(file)):
            txt+=file[i:i+line_length]+"\n"
            i+=line_length

    ax.text(0.5, -0.5,txt[:10], color='white', fontsize=15,
         horizontalalignment='center',
         verticalalignment='bottom',
         transform = ax.transAxes)
