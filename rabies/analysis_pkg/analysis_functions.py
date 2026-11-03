import numpy as np

'''
SBC
'''

def compute_seed_FC(timeseries, seed_arr_dict):
    from .analysis_math import vcorrcoef
    SBC_dict = {}
    for seed_name, roi_mask_vec in seed_arr_dict.items():
        roi_mask = roi_mask_vec.astype(bool)
        seed_timeseries = timeseries[:, roi_mask].mean(axis=1)
        corrs = vcorrcoef(timeseries.T, seed_timeseries)
        corrs[np.isnan(corrs)] = 0
        SBC_dict[seed_name] = (corrs, seed_timeseries)
    return SBC_dict

'''
FC matrix
'''

def parcellated_FC_matrix(sub_timeseries, atlas_idx, roi_list):
    
    timeseries_dict = {}
    for i in roi_list:
        roi_mask = np.asarray(atlas_idx == i, dtype=bool)
        # extract the voxel timeseries within the mask, and take the mean ROI timeseries
        timeseries_dict[str(i)] = sub_timeseries[:,roi_mask].mean(axis=1)

    roi_labels = list(timeseries_dict.keys())
    sub_timeseries = []
    for roi in roi_labels:
        sub_timeseries.append(timeseries_dict[roi])
    corr_matrix = np.corrcoef(sub_timeseries)
    return corr_matrix,roi_labels


def plot_matrix(filename, corr_matrix):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    g = ax.imshow(corr_matrix, cmap='coolwarm', vmax=1, vmin=-1)
    ax.axis('off')
    cbar = plt.colorbar(g, ax=ax, shrink=0.5)
    cbar.set_label('R score', rotation=270, fontsize=10)
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight', dpi=150)


'''
ICA
'''

def run_group_ICA(bold_file_list, mask_file, dim, random_seed, background_image, disableMigp=False):
    import os
    import pandas as pd

    # create a filelist.txt
    file_path = os.path.abspath('filelist.txt')
    from rabies.utils import flatten_list
    merged = flatten_list(list(bold_file_list))
    df = pd.DataFrame(data=merged)
    df.to_csv(file_path, header=False, sep=',', index=False)

    from rabies.utils import run_command
    out_dir = os.path.abspath('group_melodic.ica')
    command = f'melodic -i {file_path} -m {mask_file} -o {out_dir} -d {dim} --report --seed={str(random_seed)} --bgimage={background_image}'
    if disableMigp:
        command+=' --disableMigp'
    rc,c_out = run_command(command)
    IC_file = out_dir+'/melodic_IC.nii.gz'
    return out_dir, IC_file
