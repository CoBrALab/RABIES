from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.utility import Function


def init_anat_mask_prep_wf(csv_labels, name='anat_prep_mask_wf'):
    '''
    This workflow will take the output masks and labels from pydpiper for each
    subject, the transform of each subject,
    '''

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=["subject_id", "session", 'labels']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['WM_mask', 'CSF_mask', 'eroded_WM_mask', 'eroded_CSF_mask']), name='outputnode')

    compute_anat_masks = pe.Node(Function(input_names=['atlas', 'csv_labels', "subject_id", "session"],
                              output_names=['WM_mask_file', 'CSF_mask_file', 'eroded_WM_mask_file', 'eroded_CSF_mask_file'],
                              function=compute_masks),
                     name='compute_anat_masks')
    compute_anat_masks.inputs.csv_labels = csv_labels

    workflow.connect([
        (inputnode, compute_anat_masks, [
            ("labels", "atlas"),
            ("session", "session"),
            ("subject_id", "subject_id"),
            ]),
        (compute_anat_masks, outputnode, [
            ("WM_mask_file", "WM_mask"),
            ("eroded_WM_mask_file", "eroded_WM_mask"),
            ("CSF_mask_file", "CSF_mask"),
            ("eroded_CSF_mask_file", "eroded_CSF_mask"),
            ]),
    ])

    return workflow

def apply_transform(subject_id, session, reference_image,transforms,input_image, file_spec):
    import os
    cwd = os.getcwd()
    output_image='%s/%s_ses-%s_%s' % (cwd, subject_id, session, file_spec)
    os.system('antsApplyTransforms -d 3 -i %s -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,transforms,reference_image,output_image,))

    return output_image

def compute_masks(atlas, csv_labels, subject_id, session):
    import numpy as np
    import pandas as pd
    import os
    import nibabel as nb
    cwd = os.getcwd()

    '''extract rows from the csv file'''
    label_info=pd.read_csv(csv_labels)
    right_id=label_info['right label']
    left_id=label_info['left label']
    tissue_type=label_info['tissue type']

    '''create lists with the WM/CSF label numbers'''
    WM_labels=[]
    CSF_labels=[]

    for x in range(1, len(tissue_type)):
        if tissue_type[x]=='WM':
            WM_labels.append(right_id[x])
            WM_labels.append(left_id[x])
        elif tissue_type[x]=='CSF':
            CSF_labels.append(right_id[x])
            CSF_labels.append(left_id[x])

    #generate a list of the atlas labels which correspond to WM or CSF regions
    WM_labels_list=np.zeros(len(WM_labels))
    CSF_labels_list=np.zeros(len(CSF_labels))

    for x in range(0, len(WM_labels)):
        WM_labels_list[x]=int(WM_labels[x])

    for x in range(0, len(CSF_labels)):
        CSF_labels_list[x]=int(CSF_labels[x])

    '''extract the voxels which falls into the labels'''
    atlas_img=nb.load(os.path.abspath(atlas))
    atlas_data=atlas_img.dataobj
    shape=atlas_data.shape
    vector_atlas=np.transpose(np.reshape(atlas_data, [shape[0]*shape[1]*shape[2]]))
    WM_voxels=np.zeros(vector_atlas.shape)
    CSF_voxels=np.zeros(vector_atlas.shape)

    #for each WM/CSF label, establish which voxels are part of these labels
    for x in range(0,len(WM_labels_list)):
        WM_voxels=WM_voxels+(vector_atlas==WM_labels_list[x])
    WM_voxels=(WM_voxels>=1).astype(int) #binarize the mask

    for x in range(0,len(CSF_labels_list)):
        CSF_voxels=CSF_voxels+(vector_atlas==CSF_labels_list[x])
    CSF_voxels=(CSF_voxels>=1).astype(int) #binarize the mask

    WM_mask=np.reshape(WM_voxels, shape)
    CSF_mask=np.reshape(CSF_voxels, shape)

    WM_mask_file='%s/%s_ses-%s_full_WM_mask.nii.gz' % (cwd, subject_id, session,)
    CSF_mask_file='%s/%s_ses-%s_full_CSF_mask.nii.gz' % (cwd, subject_id, session,)

    nb.Nifti1Image(WM_mask, atlas_img.affine, atlas_img.header).to_filename(WM_mask_file)
    nb.Nifti1Image(CSF_mask, atlas_img.affine, atlas_img.header).to_filename(CSF_mask_file)

    '''Erode the masks'''
    from scipy.ndimage.morphology import binary_erosion
    eroded_WM_mask=binary_erosion(WM_mask, iterations=1)
    eroded_CSF_mask=binary_erosion(CSF_mask, iterations=1)

    eroded_WM_mask_file='%s/%s_ses-%s_eroded_WM_mask.nii.gz' % (cwd, subject_id, session,)
    eroded_CSF_mask_file='%s/%s_ses-%s_eroded_CSF_mask.nii.gz' % (cwd, subject_id, session,)

    nb.Nifti1Image(eroded_WM_mask, atlas_img.affine, atlas_img.header).to_filename(eroded_WM_mask_file)
    nb.Nifti1Image(eroded_CSF_mask, atlas_img.affine, atlas_img.header).to_filename(eroded_CSF_mask_file)

    return WM_mask_file, CSF_mask_file, eroded_WM_mask_file, eroded_CSF_mask_file


def mnc2nii(mnc_file):
    import os
    cwd = os.getcwd()
    basename=os.path.basename(mnc_file).split('.')[0]
    os.system('mnc2nii %s %s/%s.nii' % (mnc_file,cwd,basename))
    os.system('gzip *.nii')
    return '%s/%s.nii.gz' % (cwd,basename)
