import numpy as np
import nibabel as nb


def seed_based_FC(bold_file, brain_mask, seed_dict, seed_name):
    import os
    import nibabel as nb
    import numpy as np
    import pathlib
    from rabies.analysis_pkg.analysis_functions import seed_corr

    seed_file = seed_dict[seed_name]

    mask_array = np.asarray(nb.load(brain_mask).dataobj)
    mask_vector = mask_array.reshape(-1)
    mask_indices = (mask_vector==True)
    mask_vector = mask_vector.astype(float)

    mask_vector[mask_indices] = seed_corr(bold_file, brain_mask, seed_file)
    corr_map = mask_vector.reshape(mask_array.shape)

    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")[0]
    corr_map_file = os.path.abspath(
        filename_split+'_'+seed_name+'_corr_map.nii.gz')
    nb.Nifti1Image(corr_map, nb.load(brain_mask).affine, nb.load(
        brain_mask).header).to_filename(corr_map_file)
    return corr_map_file


def seed_corr(bold_file, brain_mask, seed):
    import os
    from nilearn.input_data import NiftiMasker

    resampled = os.path.abspath('resampled.nii.gz')
    os.system('antsApplyTransforms -i %s -r %s -o %s -n GenericLabel' %
              (seed, brain_mask, resampled))

    masker = NiftiMasker(mask_img=nb.load(resampled), standardize=False, verbose=0)
    # extract the voxel timeseries within the mask
    voxel_seed_timeseries = masker.fit_transform(bold_file)
    # take the mean ROI timeseries
    seed_timeseries = np.mean(voxel_seed_timeseries, axis=1)

    mask_array = np.asarray(nb.load(brain_mask).dataobj)
    mask_vector = mask_array.reshape(-1)
    mask_indices = (mask_vector==True)

    timeseries_array = np.asarray(nb.load(bold_file).dataobj)
    sub_timeseries = np.zeros([mask_indices.sum(), timeseries_array.shape[3]])
    for t in range(timeseries_array.shape[3]):
        sub_timeseries[:, t] = (
            timeseries_array[:, :, :, t].reshape(-1))[mask_indices]

    corrs = vcorrcoef(sub_timeseries, seed_timeseries)
    corrs[np.isnan(corrs)] = 0
    return corrs


def vcorrcoef(X, y):  # return a correlation between each row of X with y
    Xm = np.reshape(np.mean(X, axis=1), (X.shape[0], 1))
    ym = np.mean(y)
    r_num = np.sum((X-Xm)*(y-ym), axis=1)
    r_den = np.sqrt(np.sum((X-Xm)**2, axis=1)*np.sum((y-ym)**2))
    r = r_num/r_den
    return r


def get_CAPs(data, volumes, n_clusters):
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=n_clusters, n_init=10, max_iter=300)
    kmeans.fit(data)
    cluster_labels = kmeans.labels_
    CAPs = []
    for cluster in range(n_clusters):
        timepoints = (cluster_labels == cluster)
        CAPs.append(volumes[timepoints, :].mean(axis=0))
    return CAPs, cluster_labels


def recover_3D(mask_file, vector_map):
    brain_mask = np.asarray(nb.load(mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)
    volume = np.zeros(brain_mask.shape)
    volume[volume_indices] = vector_map
    volume_img = nb.Nifti1Image(volume, nb.load(
        mask_file).affine, nb.load(mask_file).header)
    return volume_img


def recover_3D_multiple(mask_file, vector_maps):
    # vector maps of shape num_volumeXnum_voxel
    brain_mask = np.asarray(nb.load(mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)
    shape = (brain_mask.shape[0], brain_mask.shape[1],
             brain_mask.shape[2], vector_maps.shape[0])
    volumes = np.zeros(shape)
    for i in range(vector_maps.shape[0]):
        volume = volumes[:, :, :, i]
        volume[volume_indices] = vector_maps[i, :]
        volumes[:, :, :, i] = volume
    volume_img = nb.Nifti1Image(volumes, nb.load(
        mask_file).affine, nb.load(mask_file).header)
    return volume_img


def threshold_maps(vector_maps, fraction):
    # vector_maps of shape map_numxvoxels
    num_voxels = int(vector_maps.shape[1]*fraction)
    thresholded_maps = np.zeros(vector_maps.shape)
    binary_maps = np.zeros(vector_maps.shape)
    for i in range(vector_maps.shape[0]):
        vector_map = vector_maps[i, :]
        idx = vector_map.argsort()[-num_voxels:]
        thresholded_maps[i, idx] = vector_map[idx]
        binary_maps[i, idx] = 1
    return [thresholded_maps, binary_maps.astype(bool)]


'''
FC matrix
'''


def run_FC_matrix(bold_file, mask_file, atlas, roi_type='parcellated'):
    import os
    import pandas as pd
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")
    figname = os.path.abspath(filename_split[0]+'_FC_matrix.png')

    from rabies.analysis_pkg.analysis_functions import parcellated_FC_matrix, voxelwise_FC_matrix, plot_matrix
    if roi_type == 'parcellated':
        corr_matrix = parcellated_FC_matrix(bold_file, atlas)
    elif roi_type == 'voxelwise':
        corr_matrix = voxelwise_FC_matrix(bold_file, mask_file)
    else:
        raise ValueError(
            "Invalid --ROI_type provided: %s. Must be either 'parcellated' or 'voxelwise.'" % (roi_type))
    plot_matrix(figname, corr_matrix)

    data_file = os.path.abspath(filename_split[0]+'_FC_matrix.csv')
    df = pd.DataFrame(corr_matrix)
    df.to_csv(data_file, sep=',')
    return data_file, figname


def voxelwise_FC_matrix(bold_file, mask_file):
    brain_mask = np.asarray(nb.load(mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)

    timeseries_array = np.asarray(nb.load(bold_file).dataobj)
    sub_timeseries = np.zeros(
        [timeseries_array.shape[3], volume_indices.sum()])
    for t in range(timeseries_array.shape[3]):
        sub_timeseries[t, :] = timeseries_array[:, :, :, t][volume_indices]

    corr_matrix = np.corrcoef(sub_timeseries.T)
    return corr_matrix


def extract_timeseries(bold_file, atlas):
    from nilearn.input_data import NiftiMasker
    atlas_img = nb.load(atlas)
    atlas_data = np.asarray(atlas_img.dataobj)
    max_int = atlas_data.max()

    timeseries_dict = {}
    for i in range(1, max_int+1):
        if np.max(i == atlas_data):  # taking a ROI only if it has labeled voxels
            roi_mask = np.asarray(atlas_data == i, dtype=int)
            roi_mask_nb = nb.Nifti1Image(roi_mask, nb.load(
                atlas).affine, nb.load(atlas).header)
            masker = NiftiMasker(mask_img=roi_mask_nb,
                                 standardize=False, verbose=0)
            # extract the voxel timeseries within the mask
            voxel_timeseries = masker.fit_transform(bold_file)
            # take the mean ROI timeseries
            roi_timeseries = np.mean(voxel_timeseries, axis=1)
            timeseries_dict[str(i)] = roi_timeseries
    return timeseries_dict


def parcellated_FC_matrix(bold_file, atlas):
    timeseries_dict = extract_timeseries(bold_file, atlas)
    roi_labels = timeseries_dict.keys()
    sub_timeseries = []
    for roi in roi_labels:
        sub_timeseries.append(timeseries_dict[roi])
    corr_matrix = np.corrcoef(sub_timeseries)
    return corr_matrix


def plot_matrix(filename, corr_matrix):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    vmax = np.abs(corr_matrix).max()
    vmin = -vmax
    g = ax.imshow(corr_matrix, cmap='coolwarm', vmax=vmax, vmin=vmin)
    ax.axis('off')
    cbar = plt.colorbar(g, ax=ax, shrink=0.5)
    cbar.set_label('R score', rotation=270, fontsize=10)
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight', dpi=150)


'''
ICA
'''


def run_group_ICA(bold_file_list, mask_file, dim, tr):
    import os
    import pandas as pd

    # create a filelist.txt
    file_path = os.path.abspath('filelist.txt')
    from rabies.preprocess_pkg.utils import flatten_list
    merged = flatten_list(list(bold_file_list))
    df = pd.DataFrame(data=merged)
    df.to_csv(file_path, header=False, sep=',', index=False)

    from rabies.preprocess_pkg.utils import run_command
    out_dir = os.path.abspath('group_melodic.ica')
    command = 'melodic -i %s -m %s -o %s --tr=%s -d %s --report' % (
        file_path, mask_file, out_dir, tr, dim)
    rc = run_command(command)
    IC_file = out_dir+'/melodic_IC.nii.gz'
    return out_dir, IC_file


def resample_4D(input_4d, ref_file):
    import os
    import pathlib  # Better path manipulation
    import SimpleITK as sitk
    from rabies.preprocess_pkg.utils import run_command, split_volumes, Merge, copyInfo_3DImage
    rabies_data_type=sitk.sitkFloat32


    # check if the IC_file has the same dimensions as bold_file
    img_array = sitk.GetArrayFromImage(
        sitk.ReadImage(ref_file))[0, :, :, :]
    image_3d = copyInfo_3DImage(sitk.GetImageFromArray(
        img_array, isVector=False), sitk.ReadImage(ref_file))
    new_ref = 'temp_ref.nii.gz'
    sitk.WriteImage(image_3d, 'temp_ref.nii.gz')

    filename_split = pathlib.Path(
        input_4d).name.rsplit(".nii")

    # Splitting into list of single volumes
    [split_volumes_files, num_volumes] = split_volumes(
        input_4d, "split_", rabies_data_type)

    resampled_volumes = []
    for x in range(0, num_volumes):
        resampled_vol_fname = os.path.abspath(
            "resampled_volume" + str(x) + ".nii.gz")
        resampled_volumes.append(resampled_vol_fname)

        command = 'antsApplyTransforms -i %s -n BSpline[5] -r %s -o %s' % (
            split_volumes_files[x], new_ref, resampled_vol_fname)
        rc = run_command(command)
        # change image to specified data type
        sitk.WriteImage(sitk.ReadImage(resampled_vol_fname, rabies_data_type), resampled_vol_fname)

    out=Merge(in_files=resampled_volumes, header_source=input_4d, rabies_data_type=rabies_data_type, clip_negative=False).run()
    return out.outputs.out_file

def run_DR_ICA(bold_file, mask_file, IC_file):
    import os
    import pandas as pd
    import pathlib  # Better path manipulation
    import SimpleITK as sitk
    from rabies.analysis_pkg.analysis_functions import sub_DR_ICA, recover_3D_multiple, resample_4D
    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")

    # check if the IC_file has the same dimensions as bold_file
    if not sitk.ReadImage(bold_file).GetSize()[:-1]==sitk.ReadImage(IC_file).GetSize()[:-1]:
        import logging
        log = logging.getLogger('root')
        log.info('Resampling file with IC components to match the scan dimensionality.')
        IC_file = resample_4D(IC_file, bold_file)

    sub_ICs = sub_DR_ICA(bold_file, mask_file, IC_file)

    data_file = os.path.abspath(filename_split[0]+'_DR_ICA.csv')
    df=pd.DataFrame(sub_ICs.T)
    df.to_csv(data_file,sep=',')

    # save the subjects' IC maps as .nii file
    nii_file = os.path.abspath(filename_split[0]+'_DR_ICA.nii.gz')
    recover_3D_multiple(mask_file, sub_ICs).to_filename(nii_file)
    return data_file, nii_file


def sub_DR_ICA(bold_file, mask_file, IC_file):
    brain_mask = np.asarray(nb.load(mask_file).dataobj)
    volume_indices = brain_mask.astype(bool)

    timeseries_array = np.asarray(nb.load(bold_file).dataobj)
    sub_timeseries = np.zeros(
        [timeseries_array.shape[3], volume_indices.sum()])
    for t in range(timeseries_array.shape[3]):
        sub_timeseries[t, :] = timeseries_array[:, :, :, t][volume_indices]

    all_IC_array = np.asarray(nb.load(IC_file).dataobj)
    all_IC_vectors = np.zeros([all_IC_array.shape[3], volume_indices.sum()])
    for i in range(all_IC_array.shape[3]):
        all_IC_vectors[i, :] = (all_IC_array[:, :, :, i])[volume_indices]

    sub_ICs = dual_regression(all_IC_vectors, sub_timeseries)
    return sub_ICs


'''
LINEAR REGRESSION --- CLOSED-FORM SOLUTION
'''


def closed_form(X, Y, intercept=False):  # functions that computes the Least Squares Estimates
    if intercept:
        X = np.concatenate((X, np.ones([X.shape[0], 1])), axis=1)
    return np.linalg.inv(X.transpose().dot(X)).dot(X.transpose()).dot(Y)


def mse(X, Y, w):  # function that computes the Mean Square Error (MSE)
    return np.mean((Y-np.matmul(X, w))**2)


def dual_regression(IC_vectors, timeseries):
    # IC_vectors is of shape num_ICxnum_voxels
    # timeseries is of shape num_timepointsxnum_voxels
    X = IC_vectors.transpose()
    Y = timeseries.transpose()

    # spatial and temporal centering of the matrices as suggested here https://mandymejia.com/2018/03/29/the-role-of-centering-in-dual-regression/#:~:text=Dual%20regression%20requires%20centering%20across%20time%20and%20space&text=time%20points.,each%20time%20course%20at%20zero).
    # spatial centering of the group ICs
    X = X-X.mean(axis=0)
    # spatial and temporal centering of the timeseries
    Y = Y-Y.mean(axis=0)
    Y = (Y.T-Y.mean(axis=1)).T

    # for one given volume, it's values can be expressed through a linear combination of the components ()
    w = closed_form(X, Y, intercept=False)

    # normalize the component timecourses to unit variance
    std = w.std(axis=1)
    for i in range(w.shape[0]):
        w[i, :] = w[i, :]/std[i]

    # for a given voxel timeseries, it's signal can be explained a linear combination of the component timecourses
    X = w.transpose()
    Y = timeseries
    # return recovered components of dim num_ICsxnum_voxels
    return closed_form(X, Y, intercept=False)
