from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.utility import Function


def init_anat_mask_prep_wf(name='anat_prep_mask_wf'):
    '''
    This workflow will take the output masks and labels from pydpyper for each
    subject, the transform of each subject,
    '''

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=["subject_id", 'anat_preproc', 'nlin_transform', 'nlin_mask', 'labels']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['resampled_mask', 'resampled_labels']), name='outputnode')

    apply_transform_mask = pe.Node(Function(input_names=["subject_id", "reference_image",'transforms','input_image', 'file_spec'],
                              output_names=["output_image"],
                              function=apply_transform),
                     name='apply_transform_mask')
    apply_transform_mask.inputs.file_spec='_anat_mask.mnc'

    apply_transform_labels = pe.Node(Function(input_names=["subject_id", "reference_image",'transforms','input_image', 'file_spec'],
                              output_names=["output_image"],
                              function=apply_transform),
                     name='apply_transform_labels')
    apply_transform_labels.inputs.file_spec='_anat_labels.mnc'


    workflow.connect([
        (inputnode, apply_transform_mask, [
            ("anat_preproc", "reference_image"),
            ('nlin_transform', 'transforms'),
            ('nlin_mask', 'input_image'),
            ]),
        (inputnode, apply_transform_labels, [
            ("anat_preproc", "reference_image"),
            ('nlin_transform', 'transforms'),
            ('labels', 'input_image'),
            ]),
        (apply_transform_mask, outputnode, [("output_image", "resampled_mask")]),
        (apply_transform_labels, outputnode, [("output_image", "resampled_labels")]),
    ])

    return workflow

def apply_transform(subject_id, reference_image,transforms,input_image, file_spec):
    import os
    cwd = os.getcwd()
    output_image=cwd+'/'+subject_id+file_spec
    os.system('antsApplyTransforms -d 3 -i %s -t %s -r %s -o %s --verbose -n GenericLabel' % (input_image,transforms,reference_image,output_image,))

    return output_image
