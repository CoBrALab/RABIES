from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.utility import Function


def init_anat_mask_prep_wf(csv_labels, name='anat_prep_mask_wf'):
    '''
    This workflow will take the output masks and labels from pydpiper for each
    subject, the transform of each subject,
    '''

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=["subject_id", "session", 'anat_preproc', 'nlin_transform', 'nlin_mask', 'labels']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['resampled_mask', 'resampled_labels', 'WM_mask', 'CSF_mask', 'eroded_WM_mask', 'eroded_CSF_mask']), name='outputnode')

    apply_transform_mask = pe.Node(Function(input_names=["subject_id", "session", "reference_image",'transforms','input_image', 'file_spec'],
                              output_names=["output_image"],
                              function=apply_transform),
                     name='apply_transform_mask')
    apply_transform_mask.inputs.file_spec='anat_mask.mnc'

    apply_transform_labels = pe.Node(Function(input_names=["subject_id", "session", "reference_image",'transforms','input_image', 'file_spec'],
                              output_names=["output_image"],
                              function=apply_transform),
                     name='apply_transform_labels')
    apply_transform_labels.inputs.file_spec='anat_labels.mnc'

    labels2nii = pe.Node(Function(input_names=['mnc_file'],
                              output_names=['nii_file'],
                              function=mnc2nii),
                     name='labels2nii')

    mask2nii = pe.Node(Function(input_names=['mnc_file'],
                              output_names=['nii_file'],
                              function=mnc2nii),
                     name='mask2nii')

    compute_anat_masks = pe.Node(Function(input_names=['atlas', 'csv_labels', "subject_id", "session"],
                              output_names=['WM_mask_file', 'CSF_mask_file', 'eroded_WM_mask_file', 'eroded_CSF_mask_file'],
                              function=compute_masks),
                     name='compute_anat_masks')
    compute_anat_masks.inputs.csv_labels = csv_labels

    workflow.connect([
        (inputnode, apply_transform_mask, [
            ("subject_id", "subject_id"),
            ("session", "session"),
            ("anat_preproc", "reference_image"),
            ('nlin_transform', 'transforms'),
            ('nlin_mask', 'input_image'),
            ]),
        (inputnode, apply_transform_labels, [
            ("subject_id", "subject_id"),
            ("session", "session"),
            ("anat_preproc", "reference_image"),
            ('nlin_transform', 'transforms'),
            ('labels', 'input_image'),
            ]),
        (apply_transform_mask, mask2nii, [("output_image", "mnc_file")]),
        (mask2nii, outputnode, [("nii_file", "resampled_mask")]),
        (apply_transform_labels, labels2nii, [("output_image", "mnc_file")]),
        (labels2nii, outputnode, [("nii_file", "resampled_labels")]),
        (labels2nii, compute_anat_masks, [("nii_file", "atlas")]),
        (inputnode, compute_anat_masks, [
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
    import csv
    import os
    import nibabel as nb
    cwd = os.getcwd()

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
    eroded_WM_mask=binary_erosion(CSF_mask, iterations=1)
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
