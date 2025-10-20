import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.utility import Function
from .utils import BIDSDataGraber, convert_to_RAS, correct_oblique_affine, convert_3dWarp, apply_autobox, log_transform_nii
from .commonspace_reg import join_iterables

def init_prep_input_wf(data_dir_path, split_name, scan_info, run_iter, opts, name='prep_input_wf'):
    '''
    This workflow organizes the entire processing.

    **Parameters**

        data_dir_path
            Path to the input data directory with proper BIDS folder structure.
        split_name
            input inherited from the prep_bids_iter function
        scan_info
            input inherited from the prep_bids_iter function
        run_iter
            input inherited from the prep_bids_iter function
        opts
            parser options for preprocess

    **Outputs**

        input_bold_list
            A list of input bold files
        input_anat_list
            A list of input anat files
        prep_bold_list
            A list of bold files formatted for preprocessing
        prep_anat_list
            A list of anat files formatted for preprocessing
        joined_scan_info
            A list of scan_info parameters matching bold files
        joined_run
            A list of run parameters matching bold files
        anat_joined_scan_info
            A list of scan_info parameters matching anat files
    '''

    workflow = pe.Workflow(name=name)

    # set output node
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['input_bold_list', 'input_anat_list', 'prep_bold_list', 'prep_anat_list', 'joined_scan_info', 'joined_run', 'anat_joined_scan_info']),
        name='outputnode')

    # setting up all iterables
    input_main_split = pe.Node(niu.IdentityInterface(fields=['split_name', 'scan_info']),
                         name="input_main_split")
    input_main_split.iterables = [('split_name', split_name),
                            ('scan_info', scan_info)]
    input_main_split.synchronize = True

    bold_selectfiles = pe.Node(BIDSDataGraber(bids_dir=data_dir_path, bids_filter=opts.bids_filter['func']),
                               name='bold_selectfiles')

    # node to conver input image to consistent RAS orientation
    bold_convert_to_RAS_node = pe.Node(Function(input_names=['img_file'],
                                                output_names=['RAS_file'],
                                                function=convert_to_RAS),
                                       name='bold_convert_to_RAS')

    if not opts.oblique2card=='none':
        if opts.oblique2card=='3dWarp':
            correct_oblique=convert_3dWarp
        elif opts.oblique2card=='affine':
            correct_oblique=correct_oblique_affine
        bold_oblique2card_node = pe.Node(Function(input_names=['input'],
                                                    output_names=['output'],
                                                    function=correct_oblique),
                                        name='bold_oblique2card')
        workflow.connect([
            (bold_selectfiles, bold_oblique2card_node, [
                ('out_file', 'input'),
                ]),
            (bold_oblique2card_node, bold_convert_to_RAS_node, [
                ('output', 'img_file'),
                ]),
            ])
    else:
        workflow.connect([
            (bold_selectfiles, bold_convert_to_RAS_node, [
                ('out_file', 'img_file'),
                ]),
            ])

    format_bold_buffer = pe.Node(niu.IdentityInterface(fields=['formatted_bold']),
                                        name="format_bold_buffer")

    if opts.bold_autobox: # apply AFNI's 3dAutobox
        bold_autobox = pe.Node(Function(input_names=['in_file'],
                                                    output_names=['out_file'],
                                                    function=apply_autobox),
                                        name='bold_autobox')
        workflow.connect([
            (bold_convert_to_RAS_node, bold_autobox, [
                ('RAS_file', 'in_file'),
                ]),
            (bold_autobox, format_bold_buffer, [
                ('out_file', 'formatted_bold'),
                ]),
            ])
    else:
        workflow.connect([
            (bold_convert_to_RAS_node, format_bold_buffer, [
                ("RAS_file", "formatted_bold"),
                ]),
            ])

    workflow.connect([
        (input_main_split, bold_selectfiles, [
            ("scan_info", "scan_info"),
            ]),
        ])

    if not opts.bold_only:
        input_run_split = pe.Node(niu.IdentityInterface(fields=['run', 'split_name']),
                            name="input_run_split")
        input_run_split.itersource = ('input_main_split', 'split_name')
        input_run_split.iterables = [('run', run_iter)]

        anat_selectfiles = pe.Node(BIDSDataGraber(bids_dir=data_dir_path, bids_filter=opts.bids_filter['anat']),
                                   name='anat_selectfiles')
        anat_selectfiles.inputs.run = None

        anat_convert_to_RAS_node = pe.Node(Function(input_names=['img_file'],
                                                    output_names=['RAS_file'],
                                                    function=convert_to_RAS),
                                           name='anat_convert_to_RAS')

        if not opts.oblique2card=='none':
            if opts.oblique2card=='3dWarp':
                correct_oblique=convert_3dWarp
            elif opts.oblique2card=='affine':
                correct_oblique=correct_oblique_affine
            anat_oblique2card_node = pe.Node(Function(input_names=['input'],
                                                        output_names=['output'],
                                                        function=correct_oblique),
                                            name='anat_oblique2card')
            workflow.connect([
                (anat_selectfiles, anat_oblique2card_node, [
                    ('out_file', 'input'),
                    ]),
                (anat_oblique2card_node, anat_convert_to_RAS_node, [
                    ('output', 'img_file'),
                    ]),
                ])
        else:
            workflow.connect([
                (anat_selectfiles, anat_convert_to_RAS_node, [
                    ('out_file', 'img_file'),
                    ]),
                ])

        format_anat_buffer = pe.Node(niu.IdentityInterface(fields=['formatted_anat']),
                                            name="format_anat_buffer")


        if opts.log_transform:
            log_anat_node = pe.Node(Function(input_names=['in_nii'],
                                                            output_names=[
                                                                'log_nii'],
                                                            function=log_transform_nii),
                                                    name='log_anat_node')

        if opts.anat_autobox: # apply AFNI's 3dAutobox
            anat_autobox = pe.Node(Function(input_names=['in_file'],
                                                        output_names=['out_file'],
                                                        function=apply_autobox),
                                            name='anat_autobox')
            workflow.connect([
                (anat_convert_to_RAS_node, anat_autobox, [
                    ('RAS_file', 'in_file'),
                    ]),
                ])
            if opts.log_transform:
                workflow.connect([
                    (anat_autobox, log_anat_node, [
                        ("out_file", "in_nii"),
                        ]),
                    (log_anat_node, format_anat_buffer, [
                        ("log_nii", "formatted_anat"),
                        ]),
                    ])
            else:
                workflow.connect([
                    (anat_autobox, format_anat_buffer, [
                        ('out_file', 'formatted_anat'),
                        ]),
                    ])
        else:
            if opts.log_transform:
                workflow.connect([
                    (anat_convert_to_RAS_node, log_anat_node, [
                        ("RAS_file", "in_nii"),
                        ]),
                    (log_anat_node, format_anat_buffer, [
                        ("log_nii", "formatted_anat"),
                        ]),
                    ])
            else:
                workflow.connect([
                    (anat_convert_to_RAS_node, format_anat_buffer, [
                        ("RAS_file", "formatted_anat"),
                        ]),
                    ])


        workflow.connect([
            (input_main_split, input_run_split, [
                ("split_name", "split_name"),
                ]),
            (input_main_split, anat_selectfiles,
             [("scan_info", "scan_info")]),
            (input_run_split, bold_selectfiles, [
                ("run", "run"),
                ]),
            ])

    if not opts.bold_only:
        workflow, source_join_anat_list, merged_join_anat_list = join_iterables(workflow=workflow, joinsource_list=['input_main_split'], node_prefix='anat_list', num_inputs=3)
        workflow.connect([
            (format_anat_buffer, source_join_anat_list, [
                ("formatted_anat", "file_list0"),
                ]),
            (anat_selectfiles, source_join_anat_list, [
                ("out_file", "file_list1"),
                ]),
            (input_main_split, source_join_anat_list, [
                ("scan_info", "file_list2"),
                ]),
            (merged_join_anat_list, outputnode, [
                ("file_list0", "prep_anat_list"),
                ("file_list1", "input_anat_list"),
                ("file_list2", "anat_joined_scan_info"),
                ]),
            ])
        joinsource_list=['input_run_split','input_main_split']
        bold_join_inputs = 4
    else:
        joinsource_list=['input_main_split']
        bold_join_inputs = 3
    workflow, source_join_bold_list, merged_join_bold_list = join_iterables(workflow=workflow, joinsource_list=joinsource_list, node_prefix='bold_list', num_inputs=bold_join_inputs)

    workflow.connect([
        (format_bold_buffer, source_join_bold_list, [
            ("formatted_bold", "file_list0"),
            ]),
        (bold_selectfiles, source_join_bold_list, [
            ("out_file", "file_list1"),
            ]),
        (input_main_split, source_join_bold_list, [
            ("scan_info", "file_list2"),
            ]),
        (merged_join_bold_list, outputnode, [
            ("file_list0", "prep_bold_list"),
            ("file_list1", "input_bold_list"),
            ("file_list2", "joined_scan_info"),
            ]),
        ])
    if opts.bold_only:
        workflow.connect([
            (merged_join_bold_list, outputnode, [
                ("file_list0", "prep_anat_list"),
                ("file_list1", "input_anat_list"),
                ("file_list2", "anat_joined_scan_info"),
                ]),
            ])
    else:
        workflow.connect([
            (input_run_split, source_join_bold_list, [
                ("run", "file_list3"),
                ]),
            (merged_join_bold_list, outputnode, [
                ("file_list3", "joined_run"),
                ]),
            ])
 
    return workflow


def match_iterables(scan_info, run, joined_scan_info, joined_run, joined_file_list):
    from rabies.utils import flatten_list
    joined_file_list = flatten_list(joined_file_list)
    joined_scan_info = flatten_list(joined_scan_info)
    if run is not None:
        joined_run = flatten_list(joined_run)
    else: # make run redundant with scan_info so that the code below runs 
        joined_run = joined_scan_info
        run = scan_info

    if not len(joined_file_list)==len(joined_scan_info):
        raise ValueError(f'The file list {joined_file_list} is not the same length as the scan info list {joined_scan_info}.')
    if not len(joined_file_list)==len(joined_run):
        raise ValueError(f'The file list {joined_file_list} is not the same length as the run list {joined_run}.')

    selected_file = None
    for file_i,scan_info_i,run_i in zip(joined_file_list,joined_scan_info, joined_run):
        if (scan_info_i == scan_info) and (run_i == run):
            if selected_file is None:
                selected_file = file_i
            else:
                raise ValueError(f'Duplicate files found for {scan_info} and run {run}.')
    if selected_file is None:
        raise ValueError(f'No file found for {scan_info} and run {run}.')
    return selected_file
