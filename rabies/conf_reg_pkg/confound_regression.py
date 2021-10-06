#!/usr/bin/env python3
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function
from .utils import regress, prep_CR, exec_ICA_AROMA

def init_confound_regression_wf(cr_opts, name="confound_regression_wf"):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=[
                        'bold_file', 'brain_mask', 'csf_mask', 'confounds_file', 'FD_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=[
                         'cleaned_path', 'aroma_out', 'VE_file', 'frame_mask_file', 'CR_data_dict']), name='outputnode')

    regress_node = pe.Node(Function(input_names=['bold_file', 'data_dict', 'brain_mask_file', 'cr_opts'],
                                    output_names=['cleaned_path', 'VE_file_path', 'frame_mask_file', 'data_dict'],
                                    function=regress),
                           name='regress', mem_gb=1)
    regress_node.inputs.cr_opts = cr_opts

    prep_CR_node = pe.Node(Function(input_names=['bold_file', 'confounds_file', 'FD_file', 'cr_opts'],
                                              output_names=['data_dict'],
                                              function=prep_CR),
                                     name='prep_CR', mem_gb=1)
    prep_CR_node.inputs.cr_opts = cr_opts

    workflow.connect([
        (inputnode, prep_CR_node, [
            ("bold_file", "bold_file"),
            ("confounds_file", "confounds_file"),
            ("FD_file", "FD_file"),
            ]),
        (inputnode, regress_node, [
            ("brain_mask", "brain_mask_file"),
            ]),
        (prep_CR_node, regress_node, [
            ("data_dict", "data_dict"),
            ]),
        (regress_node, outputnode, [
            ("cleaned_path", "cleaned_path"),
            ("VE_file_path", "VE_file"),
            ("frame_mask_file", "frame_mask_file"),
            ("data_dict", "CR_data_dict"),
            ]),
        ])

    if cr_opts.run_aroma:
        ica_aroma_node = pe.Node(Function(input_names=['inFile', 'mc_file', 'brain_mask', 'csf_mask', 'tr', 'aroma_dim'],
                                          output_names=['cleaned_file', 'aroma_out'],
                                          function=exec_ICA_AROMA),
                                 name='ica_aroma', mem_gb=1)
        ica_aroma_node.inputs.tr = cr_opts.TR
        ica_aroma_node.inputs.aroma_dim = cr_opts.aroma_dim

        workflow.connect([
            (inputnode, ica_aroma_node, [
                ("bold_file", "inFile"),
                ("brain_mask", "brain_mask"),
                ("confounds_file", "mc_file"),
                ("csf_mask", "csf_mask"),
                ]),
            (ica_aroma_node, regress_node, [
                ("cleaned_file", "bold_file"),
                ]),
            (ica_aroma_node, outputnode, [
                ("aroma_out", "aroma_out"),
                ]),
            ])
    else:
        workflow.connect([
            (inputnode, regress_node, [
                ("bold_file", "bold_file"),
                ]),
            ])

    return workflow
