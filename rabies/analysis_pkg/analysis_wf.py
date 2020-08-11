import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function

def init_analysis_wf(name="analysis_wf"):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['bold_file_list', 'brain_mask', 'atlas_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['cleaned_path', 'aroma_out', 'mel_out','tSNR_file','corr_map_list']), name='outputnode')

    print_file_node = pe.Node(Function(input_names=['file_list'],
                              output_names=['out'],
                              function=print_file),
                     name='print_file', mem_gb=1)

    workflow.connect([
        (inputnode, print_file_node, [
            ("bold_file_list", "file_list"),
            ]),
        ])

    return workflow

def print_file(file_list):
    print(file_list)
