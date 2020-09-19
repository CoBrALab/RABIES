#!/usr/bin/env python3
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function
from .utils import regress, data_diagnosis, select_timeseries


def init_confound_regression_wf(lowpass=None, highpass=None, smoothing_filter=0.3, run_aroma=False, aroma_dim=0, conf_list=[],
                                TR='1.0s', apply_scrubbing=False, scrubbing_threshold=0.1, timeseries_interval='all', diagnosis_output=False, name="confound_regression_wf"):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=[
                        'bold_file', 'brain_mask', 'csf_mask', 'confounds_file', 'FD_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=[
                         'cleaned_path', 'cr_out', 'mel_out', 'tSNR_file']), name='outputnode')

    regress_node = pe.Node(Function(input_names=['bold_file', 'brain_mask_file', 'confounds_file', 'csf_mask', 'FD_file', 'conf_list',
                                                 'TR', 'lowpass', 'highpass', 'smoothing_filter', 'run_aroma', 'aroma_dim', 'apply_scrubbing', 'scrubbing_threshold', 'timeseries_interval'],
                                    output_names=['cleaned_path',
                                                  'bold_file', 'cr_out'],
                                    function=regress),
                           name='regress', mem_gb=1)
    regress_node.inputs.conf_list = conf_list
    regress_node.inputs.TR = float(TR.split('s')[0])
    regress_node.inputs.lowpass = lowpass
    regress_node.inputs.highpass = highpass
    regress_node.inputs.smoothing_filter = smoothing_filter
    regress_node.inputs.run_aroma = run_aroma
    regress_node.inputs.aroma_dim = aroma_dim
    regress_node.inputs.apply_scrubbing = apply_scrubbing
    regress_node.inputs.scrubbing_threshold = scrubbing_threshold
    regress_node.inputs.timeseries_interval = timeseries_interval

    workflow.connect([
        (inputnode, regress_node, [
            ("brain_mask", "brain_mask_file"),
            ("confounds_file", "confounds_file"),
            ("csf_mask", "csf_mask"),
            ("FD_file", "FD_file"),
            ]),
        (regress_node, outputnode, [
            ("cleaned_path", "cleaned_path"),
            ("cr_out", "cr_out"),
            ]),
        ])

    if not timeseries_interval == 'all':
        select_timeseries_node = pe.Node(Function(input_names=['bold_file', 'timeseries_interval'],
                                                  output_names=['bold_file'],
                                                  function=select_timeseries),
                                         name='select_timeseries', mem_gb=1)
        select_timeseries_node.inputs.timeseries_interval = timeseries_interval

        workflow.connect([
            (inputnode, select_timeseries_node, [
                ("bold_file", "bold_file"),
                ]),
            (select_timeseries_node, regress_node, [
                ("bold_file", "bold_file"),
                ]),
            ])
    else:
        workflow.connect([
            (inputnode, regress_node, [
                ("bold_file", "bold_file"),
                ]),
            ])

    if diagnosis_output:
        data_diagnosis_node = pe.Node(data_diagnosis(),
                                      name='data_diagnosis', mem_gb=1)

        workflow.connect([
            (inputnode, data_diagnosis_node, [
                ("brain_mask", "brain_mask_file"),
                ]),
            (regress_node, data_diagnosis_node, [
                ("cleaned_path", "cleaned_path"),
                ("bold_file", "bold_file"),
                ]),
            (data_diagnosis_node, outputnode, [
                ("mel_out", "mel_out"),
                ("tSNR_file", "tSNR_file"),
                ]),
            ])

    return workflow
