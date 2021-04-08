import os
import numpy as np
import nibabel as nb
import pickle
import matplotlib.pyplot as plt
from rabies.analysis_pkg import analysis_functions
from nilearn.plotting import plot_stat_map

#import torch
#import prior_modeling
#device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


'''
Workflow structure:
Step 1: resample the WM/CSF/GM/right/left hem masks and the group-ICA onto subject space, or if commonspace only do it once onto the template
*make it as an option to select the DSURQE-generated regional masks or use the WM/CSF from the inputs
Step 2: prep the subject-level data
Step 3: generate teh subject-level figures
Step 4: run the group level
'''

from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


class ScanDiagnosisInputSpec(BaseInterfaceInputSpec):
    bold_file = File(exists=True, mandatory=True, desc="4D EPI file")
    brain_mask_file = File(exists=True, mandatory=True, desc="Brain mask.")
    WM_mask_file = File(exists=True, mandatory=True, desc="WM mask.")
    CSF_mask_file = File(exists=True, mandatory=True, desc="CSF mask.")
    IC_file = File(exists=True, mandatory=True, desc="MELODIC ICA components to use.")
    IC_bold_idx = traits.List(desc="The index for the ICA components that correspond to bold sources.")
    IC_confound_idx = traits.List(desc="The index for the ICA components that correspond to confounding sources.")
    CR_data_dict = traits.Dict(
        desc="specify if should detect and remove dummy scans, and use these volumes as reference image.")
    DSURQE_regions = traits.Bool(
        desc="Whether to use the regional masks generated from the DSURQE atlas for the grayplots outputs. Requires using the DSURQE template for preprocessing.")
    transforms = traits.List(default = [], desc="List of transforms to apply to the IC file and other commonspace masks.")
    inverses = traits.List(default = [],
        desc="Define whether some transforms must be inverse, with a boolean list where true defines inverse e.g.[0,1,0]")


class ScanDiagnosisOutputSpec(TraitedSpec):
    figure_path = File(
        exists=True, desc="Output figure from the scan diagnosis")
    temporal_info = traits.Dict(desc="A dictionary regrouping the temporal features.")
    spatial_info = traits.Dict(desc="A dictionary regrouping the spatial features.")

class ScanDiagnosis(BaseInterface):
    """

    """

    input_spec = ScanDiagnosisInputSpec
    output_spec = ScanDiagnosisOutputSpec

    def _run_interface(self, runtime):

        brain_mask_file = self.inputs.brain_mask_file
        WM_mask_file = self.inputs.WM_mask_file
        CSF_mask_file = self.inputs.CSF_mask_file

        if self.inputs.DSURQE_regions:
            right_hem_mask_file = resample_mask(os.environ['RABIES']+'/template_files/DSURQE_100micron_right_hem_mask.nii.gz',
                    brain_mask_file, self.inputs.transforms, self.inputs.inverses)
            left_hem_mask_file = resample_mask(os.environ['RABIES']+'/template_files/DSURQE_100micron_left_hem_mask.nii.gz',
                    brain_mask_file, self.inputs.transforms, self.inputs.inverses)
        else:
            right_hem_mask_file=''
            left_hem_mask_file=''

        IC_file = resample_IC_file(self.inputs.IC_file, brain_mask_file, self.inputs.transforms, self.inputs.inverses)

        edge_mask_file = os.path.abspath('edge_mask.nii.gz')
        compute_edge_mask(brain_mask_file,edge_mask_file, num_edge_voxels=1)
        mask_file_dict = {'brain_mask':brain_mask_file, 'WM_mask':WM_mask_file, 'CSF_mask':CSF_mask_file, 'edge_mask':edge_mask_file, 'right_hem_mask':right_hem_mask_file, 'left_hem_mask':left_hem_mask_file, 'IC_file':IC_file}

        temporal_info,spatial_info = process_data(self.inputs.bold_file, self.inputs.CR_data_dict, mask_file_dict, self.inputs.IC_bold_idx, self.inputs.IC_confound_idx, prior_fit=False,prior_fit_options=[])

        import pathlib
        filename_template = pathlib.Path(self.inputs.bold_file).name.rsplit(".nii")[0]+'_scan_diagnosis.png'
        figure_path = os.path.abspath(filename_template)

        scan_diagnosis(self.inputs.bold_file,mask_file_dict,temporal_info,spatial_info, figure_path=figure_path, regional_grayplot=self.inputs.DSURQE_regions)

        setattr(self, 'figure_path', figure_path)
        setattr(self, 'temporal_info', temporal_info)
        setattr(self, 'spatial_info', spatial_info)

        return runtime

    def _list_outputs(self):
        return {'figure_path': getattr(self, 'figure_path'),
                'temporal_info': getattr(self, 'temporal_info'),
                'spatial_info': getattr(self, 'spatial_info'), }


def resample_mask(in_file, ref_file, transforms, inverses):
    # resampling the reference image to the dimension of the EPI
    from rabies.preprocess_pkg.utils import run_command
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(
        in_file).name.rsplit(".nii")
    out_file = os.path.abspath(filename_split[0])+'_resampled.nii.gz'

    # tranforms is a list of transform files, set in order of call within antsApplyTransforms
    transform_string = ""
    for transform, inverse in zip(transforms, inverses):
        if transform=='NULL':
            continue
        elif bool(inverse):
            transform_string += "-t [%s,1] " % (transform,)
        else:
            transform_string += "-t %s " % (transform,)

    command = 'antsApplyTransforms -i %s %s-n GenericLabel -r %s -o %s' % (
        in_file, transform_string, ref_file, out_file)
    rc = run_command(command)
    return out_file


def resample_IC_file(in_file, ref_file, transforms, inverses):
    # resampling the reference image to the dimension of the EPI
    import SimpleITK as sitk
    import os
    from rabies.preprocess_pkg.utils import run_command,split_volumes,copyInfo_4DImage
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(
        in_file).name.rsplit(".nii")
    out_file = os.path.abspath(filename_split[0])+'_resampled.nii.gz'

    # tranforms is a list of transform files, set in order of call within antsApplyTransforms
    transform_string = ""
    for transform, inverse in zip(transforms, inverses):
        if transform=='NULL':
            continue
        elif bool(inverse):
            transform_string += "-t [%s,1] " % (transform,)
        else:
            transform_string += "-t %s " % (transform,)

    # Splitting bold file into lists of single volumes
    [volumes_list, num_volumes] = split_volumes(
        in_file, "bold_", sitk.sitkFloat32)

    warped_volumes = []
    for x in range(0, num_volumes):
        warped_vol_fname = os.path.abspath(
            "deformed_volume" + str(x) + ".nii.gz")
        warped_volumes.append(warped_vol_fname)

        command = 'antsApplyTransforms -i %s %s-n BSpline[5] -r %s -o %s' % (
            volumes_list[x], transform_string, ref_file, warped_vol_fname)
        rc = run_command(command)


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


def compute_edge_mask(in_mask,out_file, num_edge_voxels=1):
    #custom function for computing edge mask from an input brain mask
    img=nb.load(in_mask)
    mask_array=np.asarray(img.dataobj)
    shape=mask_array.shape

    #iterate through all voxels from the three dimensions and look if it contains surrounding voxels
    edge_mask=np.zeros(shape, dtype=bool)
    num_voxel=0
    while num_voxel<num_edge_voxels:
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    #only look if the voxel is part of the mask
                    if mask_array[x,y,z]:
                        if (mask_array[x-1:x+2,y-1:y+2,z-1:z+2]==0).sum()>0:
                            edge_mask[x,y,z]=1
        mask_array=mask_array-edge_mask
        num_voxel+=1

    nb.Nifti1Image(edge_mask, img.affine, img.header).to_filename(out_file)

'''
Prepare the subject data
'''

def process_data(bold_file, data_dict, mask_file_dict, IC_bold_idx, IC_confound_idx, prior_fit=False,prior_fit_options=[]):
    temporal_info={}
    spatial_info={}

    FD_trace = data_dict['FD_trace']
    DVARS = data_dict['DVARS']
    temporal_info['FD_trace'] = FD_trace
    temporal_info['DVARS'] = DVARS
    temporal_info['VE_temporal'] = data_dict['VE_temporal']
    spatial_info['VE_spatial'] = data_dict['VE_spatial']

    WM_mask=np.asarray(nb.load(mask_file_dict['WM_mask']).dataobj).astype(bool)
    CSF_mask=np.asarray(nb.load(mask_file_dict['CSF_mask']).dataobj).astype(bool)

    edge_mask=np.asarray(nb.load(mask_file_dict['edge_mask']).dataobj).astype(bool)
    brain_mask=np.asarray(nb.load(mask_file_dict['brain_mask']).dataobj)
    volume_indices=brain_mask.astype(bool)
    edge_idx=edge_mask[volume_indices]
    WM_idx = WM_mask[volume_indices]
    CSF_idx = CSF_mask[volume_indices]
    not_edge_idx=(edge_idx==0)*(WM_idx==0)*(CSF_idx==0)

    data_array=np.asarray(nb.load(bold_file).dataobj)
    timeseries=np.zeros([data_array.shape[3],volume_indices.sum()])
    for i in range(data_array.shape[3]):
        timeseries[i,:]=(data_array[:,:,:,i])[volume_indices]

    all_IC_array=np.asarray(nb.load(mask_file_dict['IC_file']).dataobj)
    all_IC_vectors=np.zeros([all_IC_array.shape[3],volume_indices.sum()])
    for i in range(all_IC_array.shape[3]):
        all_IC_vectors[i,:]=(all_IC_array[:,:,:,i])[volume_indices]


    '''Temporal Features'''
    X = all_IC_vectors.T
    Y = timeseries.T
    # for one given volume, it's values can be expressed through a linear combination of the components
    w = analysis_functions.closed_form(X, Y, intercept=True) # take a bias into account in the model
    w = w[:-1,:] # take out the intercept

    signal_trace = np.abs(w[IC_bold_idx,:]).mean(axis=0)
    noise_trace = np.abs(w[IC_confound_idx,:]).mean(axis=0)
    temporal_info['signal_trace']=signal_trace
    temporal_info['noise_trace']=noise_trace

    # take regional timecourse from L2-norm
    global_trace=np.sqrt((timeseries.T**2).mean(axis=0))
    WM_trace=np.sqrt((timeseries.T[WM_idx]**2).mean(axis=0))
    CSF_trace=np.sqrt((timeseries.T[CSF_idx]**2).mean(axis=0))
    edge_trace=np.sqrt((timeseries.T[edge_idx]**2).mean(axis=0))
    not_edge_trace=np.sqrt((timeseries.T[not_edge_idx]**2).mean(axis=0))
    temporal_info['global_trace']=global_trace
    temporal_info['WM_trace']=WM_trace
    temporal_info['CSF_trace']=CSF_trace
    temporal_info['edge_trace']=edge_trace
    temporal_info['not_edge_trace']=not_edge_trace

    '''Spatial Features'''
    temporal_std=timeseries.std(axis=0)
    global_signal = timeseries.mean(axis=1)
    GS_corr = analysis_functions.vcorrcoef(timeseries.T, global_signal)
    DVARS_corr = analysis_functions.vcorrcoef(timeseries.T, DVARS)
    FD_corr = analysis_functions.vcorrcoef(timeseries.T, np.asarray(FD_trace))

    dr_maps = analysis_functions.dual_regression(all_IC_vectors, timeseries)
    signal_maps=dr_maps[IC_bold_idx]

    prior_out=[]
    if prior_fit:
        [prior,num_comp,convergence_function] = prior_fit_options
        X=torch.tensor(timeseries).float().to(device)
        #Wcr=torch.tensor(scan_data[scan]['CR_time']).float().to(device)
        #X_cr = X-torch.matmul(Wcr,prior_modeling.torch_closed_form(Wcr,X))

        prior_networks = torch.tensor(all_IC_vectors[IC_bold_idx,:].T).float().to(device)

        C_prior=prior_networks
        C_conf = prior_modeling.deflation_fit(X, q=num_comp, c_init=None, C_convergence=convergence_function,
                          C_prior=C_prior, W_prior=None, W_ortho=True, tol=1e-6, max_iter=200, verbose=1)

        for network in range(prior_networks.shape[1]):
            prior=prior_networks[:,network].reshape(-1,1)
            C_prior=torch.cat((prior_networks[:,:network],prior_networks[:,network+1:],C_conf),axis=1)

            C_fit = prior_modeling.deflation_fit(X, q=1, c_init=prior, C_convergence=convergence_function,
                                  C_prior=C_prior, W_prior=None, W_ortho=True, tol=1e-6, max_iter=200, verbose=1)

            # make sure the sign of weights is the same as the prior
            corr = np.corrcoef(C_fit.cpu().flatten(), prior.cpu().flatten())[0, 1]
            if corr < 0:
                C_fit = C_fit*-1

            # the finalized C
            C = torch.cat((C_fit, C_prior), axis=1)

            # L-2 norm normalization of the components
            C /= torch.sqrt((C ** 2).sum(axis=0))
            W = prior_modeling.torch_closed_form(C,X.T).T
            C_norm=C*W.std(axis=0)

            prior_out.append([C_norm[:,0],corr])

    spatial_info['temporal_std'] = temporal_std
    spatial_info['GS_corr'] = GS_corr
    spatial_info['DVARS_corr'] = DVARS_corr
    spatial_info['FD_corr'] = FD_corr
    spatial_info['DR_maps'] = signal_maps
    spatial_info['prior_modeling_maps'] = prior_out

    return temporal_info,spatial_info

'''
Subject-level QC
'''

def grayplot_regional(timeseries_file,mask_file_dict,fig,ax,fontsize=20):
    timeseries_4d=np.asarray(nb.load(timeseries_file).dataobj)

    WM_mask=np.asarray(nb.load(mask_file_dict['WM_mask']).dataobj).astype(bool)
    CSF_mask=np.asarray(nb.load(mask_file_dict['CSF_mask']).dataobj).astype(bool)
    right_hem_mask=np.asarray(nb.load(mask_file_dict['right_hem_mask']).dataobj).astype(bool)
    left_hem_mask=np.asarray(nb.load(mask_file_dict['left_hem_mask']).dataobj).astype(bool)

    grayplot_array=np.empty((0,timeseries_4d.shape[3]))
    slice_alt=np.array([])
    region_mask_label=np.array([])
    c=0
    for mask_indices in [right_hem_mask,left_hem_mask,WM_mask,CSF_mask]:
        region_mask_label=np.append(region_mask_label,np.ones(mask_indices.sum())*c)
        c+=1
        token=False
        for i in range(mask_indices.shape[1]):
            grayplot_array=np.append(grayplot_array,timeseries_4d[:,i,:,:][mask_indices[:,i,:]],axis=0)
            slice_alt=np.append(slice_alt,np.ones(mask_indices[:,i,:].sum())*token)
            token= not token

    vmax=grayplot_array.std()
    im = ax.imshow(grayplot_array, cmap='gray', vmax=vmax, vmin=-vmax)
    ax.set_aspect(grayplot_array.shape[1]/grayplot_array.shape[0])
    ax.set_xlabel('Timepoint', fontsize=fontsize)
    # increase tick size on x axis
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize/2)

    ax.set_ylabel('Voxel', fontsize=fontsize)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axes.get_yaxis().set_ticks([])
    ax.yaxis.labelpad = 40
    cbar = fig.colorbar(im, ax=ax,pad=0.02, fraction=0.046)
    cbar.set_label('Voxel Intensity', fontsize=fontsize)


    ax2 = fig.add_axes([0.018, 0.140, 0.2, 0.725])
    ax2.imshow(slice_alt.reshape(-1,1), cmap='gray', vmin=0, vmax=1.1)
    ax2.set_aspect(0.0015)
    ax2.axis('off')

    ax3 = fig.add_axes([0.004, 0.140, 0.2, 0.725])
    ax3.imshow(region_mask_label.reshape(-1,1), cmap='Spectral')
    ax3.set_aspect(0.0015)
    ax3.axis('off')

def grayplot(timeseries_file,mask_file_dict,fig,ax,fontsize=20):
    brain_mask=np.asarray(nb.load(mask_file_dict['brain_mask']).dataobj)
    volume_indices=brain_mask.astype(bool)

    data_array=np.asarray(nb.load(timeseries_file).dataobj)
    timeseries=np.zeros([data_array.shape[3],volume_indices.sum()])
    for i in range(data_array.shape[3]):
        timeseries[i,:]=(data_array[:,:,:,i])[volume_indices]

    grayplot_array=timeseries.T

    vmax=grayplot_array.std()
    im = ax.imshow(grayplot_array, cmap='gray', vmax=vmax, vmin=-vmax)
    ax.set_aspect(grayplot_array.shape[1]/grayplot_array.shape[0])
    ax.set_xlabel('Timepoint', fontsize=fontsize)
    # increase tick size on x axis
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize/2)

    ax.set_ylabel('Voxel', fontsize=fontsize)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axes.get_yaxis().set_ticks([])
    ax.yaxis.labelpad = 40
    cbar = fig.colorbar(im, ax=ax,pad=0.02, fraction=0.046)
    cbar.set_label('Voxel Intensity', fontsize=fontsize)

def scan_diagnosis(bold_file,mask_file_dict,temporal_info,spatial_info, figure_path=None, regional_grayplot=False):

    fig1,axes1 = plt.subplots(nrows=4, ncols=1, gridspec_kw = {'height_ratios':[2,1,1,1]}, figsize=(5,10), sharex=True)
    #fig3,axes3 = plt.subplots(nrows=8, ncols=n_sub,figsize=(12*n_sub,3*8))

    if regional_grayplot:
        grayplot_regional(bold_file,mask_file_dict,fig1,axes1[0],fontsize=20)
    else:
        grayplot(bold_file,mask_file_dict,fig1,axes1[0],fontsize=20)

    ax1 = axes1[1]
    #ax1.set_title(name, fontsize=15)
    ax2 = ax1.twinx()
    ax1.plot(temporal_info['FD_trace'], 'r')
    ax2.plot(temporal_info['DVARS'], 'b')
    ax2.plot(temporal_info['global_trace'], 'g')
    ax1.set_ylabel('Framewise Displacement (FD)', color='r', fontsize=12)
    ax2.set_ylabel('DVARS', color='b', fontsize=12)
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)

    ax1.set_ylim([0.0,0.2])
    ax2.set_ylim([0.0,2.0])

    ax=axes1[2]
    ax.plot(temporal_info['edge_trace'])
    ax.plot(temporal_info['WM_trace'])
    ax.plot(temporal_info['CSF_trace'])
    ax.plot(temporal_info['not_edge_trace'])
    ax.plot(temporal_info['VE_temporal'])
    ax.set_ylabel('Mask L2-norm', fontsize=12)
    ax.legend(['Edge', 'WM', 'CSF', 'Not Edge', 'CR R^2'])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim([0.0,1])

    ax=axes1[2]
    ax.plot(temporal_info['signal_trace'])
    ax.plot(temporal_info['noise_trace'])
    ax.set_ylabel('Mean Component Weight', fontsize=12)
    ax.legend(['Signal components', 'Noise components'])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim([0.0,0.09])

    fig1.tight_layout()
    if figure_path is not None:
        fig1.savefig(figure_path, bbox_inches='tight')

    '''
    ax=axes3[0,n]
    ax.set_title(name+' - Temporal STD', fontsize=15)
    plot_stat_map(analysis_functions.recover_3D(mask_file,temporal_std),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y', vmax=1, colorbar=False)
    ax=axes3[1,n]
    ax.set_title(name+' - CR R^2', fontsize=15)
    plot_stat_map(analysis_functions.recover_3D(mask_file,VE_spatial),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y', vmax=1, colorbar=False)
    ax=axes3[2,n]
    ax.set_title(name+' - Global Signal Correlation', fontsize=15)
    plot_stat_map(analysis_functions.recover_3D(mask_file,GS_corr),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y', vmax=0.7)
    ax=axes3[3,n]
    ax.set_title(name+' - DVARS Correlation', fontsize=15)
    plot_stat_map(analysis_functions.recover_3D(mask_file,DVARS_corr),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y', vmax=0.7)
    ax=axes3[4,n]
    ax.set_title(name+' - FD Correlation', fontsize=15)
    plot_stat_map(analysis_functions.recover_3D(mask_file,FD_corr),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y', vmax=0.7)

    ax=axes3[5,n]
    ax.set_title(name+' - Somatomotor', fontsize=15)
    plot_stat_map(analysis_functions.recover_3D(mask_file,signal_map1),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y', vmax=1)
    ax=axes3[6,n]
    ax.set_title(name+' - Dorsal Comp', fontsize=15)
    plot_stat_map(analysis_functions.recover_3D(mask_file,signal_map2),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y', vmax=1)
    ax=axes3[7,n]
    ax.set_title(name+' - DMN', fontsize=15)
    plot_stat_map(analysis_functions.recover_3D(mask_file,signal_map3),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y', vmax=1)

    n+=1

    fig1.tight_layout()
    if figure_path is not None:
        fig1.savefig(figure_path+'_temporal_analysis.png', bbox_inches='tight')
        fig3.savefig(figure_path+'_spatial_analysis.png', bbox_inches='tight')
    '''


'''
Group-level QC
'''

'''
def elementwise_corrcoef(X, Y):
    # X and Y are each of shape num_observations X num_element
    # computes the correlation between each element of X and Y
    Xm = X.mean(axis=0)
    Ym = Y.mean(axis=0)
    r_num = np.sum((X-Xm)*(Y-Ym), axis=0)

    r_den = np.sqrt(np.sum((X-Xm)**2, axis=0)*np.sum((Y-Ym)**2, axis=0))
    r = r_num/r_den
    return r

voxelwise_list=[]
for voxelwise_sub in voxelwise_sub_list:
    C_list=voxelwise_sub[-1]
    voxelwise_sub=np.array(voxelwise_sub[:-1])
    for C,corr in C_list:
        #if corr>0.4:
        voxelwise_sub=np.concatenate((voxelwise_sub,C.cpu().reshape(1,-1)),axis=0)
    voxelwise_list.append(voxelwise_sub)
voxelwise_array=np.array(voxelwise_list)

x_name=['temporal_std','VE_spatial','GS_corr','DVARS_corr','FD_corr','Somatomotor','Dorsal Comp','DMN', 'Prior Modeling 1', 'Prior Modeling 2', 'Prior Modeling 3']
y_name=['temporal_std','VE_spatial','GS_corr','DVARS_corr','FD_corr','Somatomotor','Dorsal Comp','DMN', 'Prior Modeling 1', 'Prior Modeling 2', 'Prior Modeling 3']

ncols=voxelwise_array.shape[1]-6
fig,axes = plt.subplots(nrows=voxelwise_array.shape[1], ncols=ncols,figsize=(12*ncols,3*voxelwise_array.shape[1]))

for i,x_label in zip(range(voxelwise_array.shape[1]),x_name):
    for j,y_label in zip(range(ncols),y_name[:-6]):
        ax=axes[i,j]
        if i<=j:
            ax.axis('off')
            continue

        X=voxelwise_array[:,i,:]
        Y=voxelwise_array[:,j,:]
        corr=elementwise_corrcoef(X, Y)

        ax.set_title('Cross-correlation for %s and %s' % (x_label,y_label), fontsize=15)
        plot_stat_map(analysis_functions.recover_3D(mask_file,corr),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y', vmax=0.7)


def closed_form_3d(X,Y):
    return np.matmul(np.matmul(np.linalg.inv(np.matmul(X.transpose(0,2,1),X)),X.transpose(0,2,1)),Y)

def lme_stats_3d(X,Y):
    #add an intercept
    X=np.concatenate((X,np.ones((X.shape[0],X.shape[1],1))),axis=2)
    [num_comparisons,num_observations,num_predictors] = X.shape
    [num_comparisons,num_observations,num_features] = Y.shape

    w=closed_form_3d(X,Y)

    residuals = Y-np.matmul(X, w)
    MSE = (((residuals)**2).sum(axis=1)/(num_observations-num_predictors))


    var_b = np.expand_dims(MSE, axis=1)*np.expand_dims(np.linalg.inv(np.matmul(X.transpose(0,2,1),X)).diagonal(axis1=1,axis2=2), axis=2)
    sd_b = np.sqrt(var_b) # standard error on the Betas
    ts_b = w/sd_b # calculate t-values for the Betas
    p_values =[2*(1-stats.t.cdf(np.abs(ts_b[:,i,:]),(num_observations-num_predictors))) for i in range(ts_b.shape[1])] # calculate a p-value map for each predictor

    return ts_b,p_values,w,residuals

non_nan_idx = (np.isnan(voxelwise_array).sum(axis=(0,1))==0)

# take out the voxels which have null values
X=voxelwise_array[:,:6,non_nan_idx].transpose(2,0,1)
Y=voxelwise_array[:,6:,non_nan_idx].transpose(2,0,1)


ts_b,p_values,w,residuals = lme_stats_3d(X,Y)

x_name=['Somatomotor','Dorsal Comp','DMN', 'Prior Modeling 1', 'Prior Modeling 2', 'Prior Modeling 3']
y_name=['group','temporal_std','VE_spatial','GS_corr','DVARS_corr','FD_corr']

fig,axes = plt.subplots(nrows=len(x_name), ncols=len(y_name),figsize=(12*len(y_name),3*len(x_name)))


for i,x_label in zip(range(len(x_name)),x_name):
    for j,y_label in zip(range(len(y_name)),y_name):
        ax=axes[i,j]

        stat_map=np.zeros(voxelwise_array.shape[2])
        stat_map[non_nan_idx]=ts_b[:,j,i]

        ax.set_title('T-value of %s on %s' % (y_label,x_label), fontsize=15)
        plot_stat_map(analysis_functions.recover_3D(mask_file,stat_map),bg_img='DSURQE.nii.gz', axes=ax, cut_coords=(0,1,2,3,4,5), display_mode='y')
'''
