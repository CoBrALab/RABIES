import numpy as np
import nibabel as nb
import os

def get_CAPs(data,volumes,n_clusters):
    from sklearn.cluster import KMeans
    kmeans=KMeans(n_clusters=n_clusters, n_init=10, max_iter=300)
    kmeans.fit(data)
    cluster_labels=kmeans.labels_
    CAPs=[]
    for cluster in range(n_clusters):
        timepoints=(cluster_labels==cluster)
        CAPs.append(volumes[timepoints,:].mean(axis=0))
    return CAPs, cluster_labels

def recover_3D(mask_file, vector_map):
    brain_mask=np.asarray(nb.load(mask_file).dataobj)
    volume_indices=brain_mask.astype(bool)
    volume=np.zeros(brain_mask.shape)
    volume[volume_indices]=vector_map
    volume_img=nb.Nifti1Image(volume, nb.load(mask_file).affine, nb.load(mask_file).header)
    return volume_img

def recover_3D_mutiple(mask_file, vector_maps):
    #vector maps of shape num_volumeXnum_voxel
    brain_mask=np.asarray(nb.load(mask_file).dataobj)
    volume_indices=brain_mask.astype(bool)
    shape=(brain_mask.shape[0],brain_mask.shape[1],brain_mask.shape[2],vector_maps.shape[0])
    volumes=np.zeros(shape)
    for i in range(vector_maps.shape[0]):
        volume=volumes[:,:,:,i]
        volume[volume_indices]=vector_maps[i,:]
        volumes[:,:,:,i]=volume
    volume_img=nb.Nifti1Image(volumes, nb.load(mask_file).affine, nb.load(mask_file).header)
    return volume_img

def threshold_maps(vector_maps, fraction):
    #vector_maps of shape map_numxvoxels
    num_voxels = int(vector_maps.shape[1]*fraction)
    thresholded_maps = np.zeros(vector_maps.shape)
    binary_maps = np.zeros(vector_maps.shape)
    for i in range(vector_maps.shape[0]):
        vector_map=vector_maps[i,:]
        idx=vector_map.argsort()[-num_voxels:]
        thresholded_maps[i,idx]=vector_map[idx]
        binary_maps[i,idx]=1
    return [thresholded_maps,binary_maps.astype(bool)]

'''
FC matrix
'''

def run_FC_matrix(bold_file,mask_file,atlas, dim_type='parcellated'):
    import os
    import pickle
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")
    figname=os.path.abspath(filename_split[0]+'_FC_matrix.png')

    if dim_type=='parcellated':
        corr_matrix=parcellated_FC_matrix(bold_file, atlas)
    elif dim_type=='voxelwise':
        corr_matrix=voxelwise_FC_matrix(bold_file, mask_file)
    plot_matrix(figname,corr_matrix)

    data_file=os.path.abspath(filename_split[0]+'_FC_matrix.pkl')
    with open(data_file, 'wb') as handle:
        pickle.dump(corr_matrix, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return data_file,figname

def voxelwise_FC_matrix(bold_file, mask_file):
    brain_mask=np.asarray(nb.load(mask_file).dataobj)
    volume_indices=brain_mask.astype(bool)

    timeseries_array=np.asarray(nb.load(bold_file).dataobj)
    sub_timeseries=np.zeros([timeseries_array.shape[3],volume_indices.sum()])
    for t in range(timeseries_array.shape[3]):
        sub_timeseries[t,:]=timeseries_array[:,:,:,t][volume_indices]

    corr_matrix=np.corrcoef(sub_timeseries.T)
    return corr_matrix

def extract_timeseries(bold_file,atlas):
    from nilearn.input_data import NiftiMasker
    atlas_img=nb.load(atlas)
    atlas_data=np.asarray(atlas_img.dataobj)
    max_int=atlas_data.max()

    timeseries_dict={}
    for i in range(1,max_int+1):
        if np.max(i==atlas_data): #taking a ROI only if it has labeled voxels
            roi_mask=np.asarray(atlas_data==i, dtype=int)
            roi_mask_nb=nb.Nifti1Image(roi_mask, nb.load(atlas).affine, nb.load(atlas).header)
            masker = NiftiMasker(mask_img=roi_mask_nb, standardize=False, verbose=0)
            voxel_timeseries = masker.fit_transform(bold_file) #extract the voxel timeseries within the mask
            roi_timeseries=np.mean(voxel_timeseries, axis=1) #take the mean ROI timeseries
            timeseries_dict[str(i)]=roi_timeseries
    return timeseries_dict


def parcellated_FC_matrix(bold_file, atlas):
    timeseries_dict=extract_timeseries(bold_file,atlas)
    roi_labels=timeseries_dict.keys()
    sub_timeseries=[]
    for roi in roi_labels:
        sub_timeseries.append(timeseries_dict[roi])
    corr_matrix=np.corrcoef(sub_timeseries)
    return corr_matrix

def plot_matrix(filename,corr_matrix):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,1, figsize=(4,4))
    vmax=np.abs(corr_matrix).max()
    vmin=-vmax
    g=ax.imshow(corr_matrix, cmap='coolwarm', vmax=vmax, vmin=vmin)
    ax.axis('off')
    cbar = plt.colorbar(g, ax=ax, shrink=0.5)
    cbar.set_label('R score', rotation=270, fontsize=10)
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight',dpi=150)


'''
ICA
'''
def run_group_ICA(bold_file_list, mask_file, dim, tr):
    import os
    import pandas as pd

    #create a filelist.txt
    file_path=os.path.abspath('filelist.txt')
    import itertools
    merged = list(itertools.chain.from_iterable(bold_file_list))
    df = pd.DataFrame(data=merged)
    df.to_csv(file_path, header=False, sep=',',index=False)

    from rabies.preprocess_bold_pkg.utils import run_command
    out_dir=os.path.abspath('group_melodic.ica')
    command='melodic -i %s -m %s -o %s --tr %s -d %s --report' % (file_path,mask_file,out_dir, tr, dim)
    rc = run_command(command)
    IC_file=out_dir+'/melodic_IC.nii.gz'
    return out_dir, IC_file


def run_DR_ICA(bold_file, mask_file, IC_file):
    import os
    import pickle
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

    sub_ICs=sub_DR_ICA(bold_file, mask_file, IC_file)

    data_file=os.path.abspath(filename_split[0]+'_DR_ICA.pkl')
    with open(data_file, 'wb') as handle:
        pickle.dump(sub_ICs, handle, protocol=pickle.HIGHEST_PROTOCOL)

    #save the subjects' IC maps as .nii file
    nii_file=os.path.abspath(filename_split[0]+'_DR_ICA.nii.gz')
    recover_3D_mutiple(mask_file, sub_ICs).to_filename(nii_file)
    return data_file,nii_file

def sub_DR_ICA(bold_file, mask_file, IC_file):
    brain_mask=np.asarray(nb.load(mask_file).dataobj)
    volume_indices=brain_mask.astype(bool)

    timeseries_array=np.asarray(nb.load(bold_file).dataobj)
    sub_timeseries=np.zeros([timeseries_array.shape[3],volume_indices.sum()])
    for t in range(timeseries_array.shape[3]):
        sub_timeseries[t,:]=timeseries_array[:,:,:,t][volume_indices]

    all_IC_array=np.asarray(nb.load(IC_file).dataobj)
    all_IC_vectors=np.zeros([all_IC_array.shape[3],volume_indices.sum()])
    for i in range(all_IC_array.shape[3]):
        all_IC_vectors[i,:]=(all_IC_array[:,:,:,i])[volume_indices]

    sub_ICs=dual_regression(all_IC_vectors, sub_timeseries)
    return sub_ICs

'''
LINEAR REGRESSION --- CLOSED-FORM SOLUTION
'''

#functions that computes the Least Squares Estimates
def closed_form(X,Y, intercept=False):
    if intercept:
        X=np.concatenate((X,np.ones([X.shape[0],1])),axis=1)
    return np.linalg.inv(X.transpose().dot(X)).dot(X.transpose()).dot(Y)

#functions that computes the Mean Square Error (MSE)
def mse(X, Y, w):
    return np.mean((Y-np.matmul(X, w))**2)

def dual_regression(IC_vectors, timeseries):
    #IC_vectors is of shape num_ICxnum_voxels
    #timeseries is of shape num_timepointsxnum_voxels
    X=IC_vectors.transpose()
    Y=timeseries.transpose()

    #for one given volume, it's values can be expressed through a linear combination of the components ()
    w=closed_form(X,Y, intercept=False)

    #normalize the component timecourses to unit variance
    std=w.std(axis=1)
    for i in range(w.shape[0]):
        w[i,:]=w[i,:]/std[i]

    #for a given voxel timeseries, it's signal can be explained a linear combination of the component timecourses
    X=w.transpose()
    Y=timeseries
    return closed_form(X,Y, intercept=False) #return recovered components of dim num_ICsxnum_voxels
