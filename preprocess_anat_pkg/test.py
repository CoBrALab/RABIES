from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.utility import Function

def init_test_wf(name='test_wf'):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['subject_id', 'anat']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_anat']), name='outputnode')


    test_node = pe.Node(Function(input_names=['file','id'],
                              output_names=['mnc_file'],
                              function=test_func),
                     name='test_node')

    workflow.connect([
        (inputnode, test_node, [("anat", "file"), ('subject_id', 'id')]),
        (test_node, outputnode, [("mnc_file", "out_anat")]),
    ])
    return workflow


def test_func(file, id):
    import os
    cwd = os.getcwd()
    os.system('mkdir -p %s/mnc_anat/' % (cwd,))
    os.system('nii2mnc %s %s/mnc_anat/%s.mnc' % (file, cwd,id))
    mnc_file='%s/mnc_anat/%s.mnc' % (cwd,id)

    return mnc_file
