import sys
import numpy as np
import csv
import os
import nibabel as nb


csv_labels=os.path.abspath('/data/chamal/projects/Gabriel_DG/data/atlases/labels/DSURQE_40micron_R_mapping.csv')
atlas=nb.load(os.path.abspath(sys.argv[1]))
prefix=os.path.abspath(sys.argv[2])

'''extract rows from the csv file'''
temp = []
with open(csv_labels) as csvfile:
    file_data = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in file_data:
        temp.append(row)

'''create lists with the WM/CSF label numbers'''
WM_labels=[]
CSF_labels=[]

for x in range(1, len(temp)):
    if temp[x][4]=='WM':
        WM_labels.append(temp[x][1])
        WM_labels.append(temp[x][2])
    elif temp[x][4]=='CSF':
        CSF_labels.append(temp[x][1])
        CSF_labels.append(temp[x][2])

array_WM_labels=np.zeros(len(WM_labels))
array_CSF_labels=np.zeros(len(CSF_labels))

for x in range(0, len(WM_labels)):
    array_WM_labels[x]=int(WM_labels[x])


for x in range(0, len(CSF_labels)):
    array_CSF_labels[x]=int(CSF_labels[x])

'''extract the voxels which falls into the labels'''
data=atlas.dataobj
shape=data.shape
resampled_data=np.transpose(np.reshape(data, [shape[0]*shape[1]*shape[2]]))
WM_voxels=np.zeros(resampled_data.shape)
CSF_voxels=np.zeros(resampled_data.shape)

for x in range(0,len(array_WM_labels)):
    WM_voxels=WM_voxels+(resampled_data==array_WM_labels[x])

WM_voxels=(WM_voxels>=1).astype(int)

for x in range(0,len(array_CSF_labels)):
    CSF_voxels=CSF_voxels+(resampled_data==array_CSF_labels[x])

CSF_voxels=(CSF_voxels>=1).astype(int)

WM_mask=np.reshape(WM_voxels, shape)
CSF_mask=np.reshape(CSF_voxels, shape)

nb.Nifti1Image(WM_mask, atlas.affine, atlas.header).to_filename(prefix+'WM_mask.nii.gz')
nb.Nifti1Image(CSF_mask, atlas.affine, atlas.header).to_filename(prefix+'CSF_mask.nii.gz')
