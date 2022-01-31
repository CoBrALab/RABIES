import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
from nipype.interfaces.afni import Autobox
from .inho_correction import init_inho_correction_wf
from .commonspace_reg import init_commonspace_reg_wf
from .bold_main_wf import init_bold_main_wf
from .utils import BIDSDataGraber, prep_bids_iter, convert_to_RAS, resample_template
from . import preprocess_visual_QC

def init_main_wf(data_dir_path, output_folder, opts, name='main_wf'):
    '''
    This workflow organizes the entire processing.

    **Parameters**

        data_dir_path
            Path to the input data directory with proper BIDS folder structure.
        output_folder
            path to output folder for the workflow and datasink
        opts
            parser options for preprocess
        cr_opts
            parser options for confound_correction
        analysis_opts
            parser options for analysis

    **Outputs**


        input_bold
            Input EPIs to the preprocessing
        commonspace_resampled_template
            the anatomical commonspace template after initial resampling
        anat_preproc
            Preprocessed anatomical image after bias field correction and denoising
        anat_mask
            Brain mask inherited from the common space registration
        anat_labels
            Anatomical labels inherited from the common space registration
        WM_mask
            Eroded WM mask inherited from the common space registration
        CSF_mask
            Eroded CSF mask inherited from the common space registration
        initial_bold_ref
            Initial EPI median volume subsequently used as 3D reference EPI volume
        inho_cor_bold
            3D reference EPI volume after bias field correction
        bold_to_anat_affine
            affine transform from the EPI space to the anatomical space
        bold_to_anat_warp
            non-linear transform from the EPI space to the anatomical space
        bold_to_anat_inverse_warp
            inverse non-linear transform from the EPI space to the anatomical space
        inho_cor_bold_warped2anat
            Bias field corrected 3D EPI volume warped to the anatomical space
        native_corrected_bold
            Preprocessed EPI resampled to match the anatomical space for
            susceptibility distortion correction
        corrected_bold_ref
            3D ref EPI volume from the native EPI timeseries
        confounds_csv
            .csv file with measured confound timecourses, including global signal,
            WM signal, CSF signal, 6 rigid body motion parameters + their first
            temporal derivate + the 12 parameters squared (24 motion parameters),
            and aCompCorr timecourses
        FD_voxelwise
            Voxelwise framewise displacement (FD) measures that can be integrated
            to future confound regression.
            These measures are computed from antsMotionCorrStats.
        pos_voxelwise
            Voxel distancing across time based on rigid body movement parameters,
            which can be integrated for a voxelwise motion regression
            These measures are computed from antsMotionCorrStats.
        FD_csv
            .csv file with global framewise displacement (FD) measures
        bold_brain_mask
            EPI brain mask for native corrected bold
        bold_WM_mask
            EPI WM mask for native corrected bold
        bold_CSF_mask
            EPI CSF mask for native corrected bold
        bold_labels
            EPI anatomical labels for native corrected bold
        commonspace_bold
            Motion and SDC-corrected EPI timeseries resampled into common space
            by applying transforms from the anatomical common space registration
        commonspace_mask
            EPI brain mask for commonspace bold
        commonspace_WM_mask
            EPI WM mask for commonspace bold
        commonspace_CSF_mask
            EPI CSF mask for commonspace bold
        commonspace_vascular_mask
            EPI vascular mask for commonspace bold
        commonspace_labels
            EPI anatomical labels for commonspace bold
        std_filename
            temporal STD map of the preprocessed timeseries
        tSNR_filename
            temporal SNR map of the preprocessed timeseries
    '''

    workflow = pe.Workflow(name=name)

    # set output node
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['input_bold', 'commonspace_resampled_template', 'anat_preproc', 'initial_bold_ref', 'inho_cor_bold', 'bold_to_anat_affine',
                'bold_to_anat_warp', 'bold_to_anat_inverse_warp', 'inho_cor_bold_warped2anat', 'native_bold', 'native_bold_ref', 'confounds_csv',
                'FD_voxelwise', 'pos_voxelwise', 'FD_csv', 'native_brain_mask', 'native_WM_mask', 'native_CSF_mask', 'native_labels',
                'commonspace_bold', 'commonspace_mask', 'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_vascular_mask',
                'commonspace_labels', 'std_filename', 'tSNR_filename']),
        name='outputnode')

    # Datasink - creates output folder for important outputs
    bold_datasink = pe.Node(DataSink(base_directory=output_folder,
                                     container="bold_datasink"),
                            name="bold_datasink")

    anat_datasink = pe.Node(DataSink(base_directory=output_folder,
                                     container="anat_datasink"),
                            name="anat_datasink")

    transforms_datasink = pe.Node(DataSink(base_directory=output_folder,
                                           container="transforms_datasink"),
                                  name="transforms_datasink")

    confounds_datasink = pe.Node(DataSink(base_directory=output_folder,
                                          container="confounds_datasink"),
                                 name="confounds_datasink")

    import bids
    bids.config.set_option('extension_initial_dot', True)
    layout = bids.layout.BIDSLayout(data_dir_path, validate=False)
    split_name, scan_info, run_iter, scan_list, bold_scan_list = prep_bids_iter(
        layout, opts.bold_only)

    # setting up all iterables
    main_split = pe.Node(niu.IdentityInterface(fields=['split_name', 'scan_info']),
                         name="main_split")
    main_split.iterables = [('split_name', split_name),
                            ('scan_info', scan_info)]
    main_split.synchronize = True

    bold_selectfiles = pe.Node(BIDSDataGraber(bids_dir=data_dir_path, suffix=[
                               'bold', 'cbv']), name='bold_selectfiles')

    # node to conver input image to consistent RAS orientation
    bold_convert_to_RAS_node = pe.Node(Function(input_names=['img_file'],
                                                output_names=['RAS_file'],
                                                function=convert_to_RAS),
                                       name='bold_convert_to_RAS')

    format_bold_buffer = pe.Node(niu.IdentityInterface(fields=['formatted_bold']),
                                        name="format_bold_buffer")

    if opts.bold_autobox: # apply AFNI's 3dAutobox
        bold_autobox = pe.Node(Autobox(padding=1, outputtype='NIFTI_GZ'),
                                name="bold_autobox")
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

    # Resample the anatomical template according to the resolution of the provided input data
    resample_template_node = pe.Node(Function(input_names=['template_file', 'mask_file', 'file_list', 'spacing', 'rabies_data_type'],
                                              output_names=[
                                                  'resampled_template', 'resampled_mask'],
                                              function=resample_template),
                                     name='resample_template', mem_gb=1*opts.scale_min_memory)
    resample_template_node.inputs.template_file = str(opts.anat_template)
    resample_template_node.inputs.mask_file = str(opts.brain_mask)
    resample_template_node.inputs.spacing = opts.anatomical_resampling
    resample_template_node.inputs.file_list = scan_list
    resample_template_node.inputs.rabies_data_type = opts.data_type

    # calculate the number of scans that will be registered
    num_scan = len(scan_list)
    num_procs = min(opts.local_threads, num_scan)

    EPI_target_buffer = pe.Node(niu.IdentityInterface(fields=['EPI_template', 'EPI_mask']),
                                        name="EPI_target_buffer")

    if opts.fast_commonspace:
        # if fast commonspace, then the inputs iterables are not merged
        source_join_common_reg = pe.Node(niu.IdentityInterface(fields=['file_list0', 'file_list1']),
                                            name="fast_commonreg_buffer")
        merged_join_common_reg = source_join_common_reg
    else:
        workflow, source_join_common_reg, merged_join_common_reg = join_iterables(workflow=workflow, joinsource_list=['main_split'], node_prefix='commonspace_reg', num_inputs=2)

    commonspace_reg_wf = init_commonspace_reg_wf(opts=opts, output_folder=output_folder, transforms_datasink=transforms_datasink, num_procs=num_procs, name='commonspace_reg_wf')

    bold_main_wf = init_bold_main_wf(opts=opts)

    # organizing visual QC outputs
    template_diagnosis = pe.Node(Function(input_names=['anat_template', 'opts', 'out_dir'],
                                       function=preprocess_visual_QC.template_info),
                              name='template_info')
    template_diagnosis.inputs.opts = opts
    template_diagnosis.inputs.out_dir = output_folder+'/preprocess_QC_report/template_files/'

    bold_inho_cor_diagnosis = pe.Node(Function(input_names=['raw_img','init_denoise','warped_mask','final_denoise', 'name_source', 'out_dir'],
                                       function=preprocess_visual_QC.inho_cor_diagnosis),
                              name='bold_inho_cor_diagnosis')
    bold_inho_cor_diagnosis.inputs.out_dir = output_folder+'/preprocess_QC_report/bold_inho_cor/'

    temporal_diagnosis = pe.Node(Function(input_names=['bold_file', 'confounds_csv', 'FD_csv', 'rabies_data_type', 'name_source', 'out_dir'],
                                          output_names=[
                                            'std_filename', 'tSNR_filename'],
                                       function=preprocess_visual_QC.temporal_features),
                              name='temporal_features')
    temporal_diagnosis.inputs.out_dir = output_folder+'/preprocess_QC_report/temporal_features/'
    temporal_diagnosis.inputs.rabies_data_type = opts.data_type

    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (main_split, bold_selectfiles, [
            ("scan_info", "scan_info"),
            ]),
        (bold_selectfiles, bold_convert_to_RAS_node, [
            ('out_file', 'img_file'),
            ]),
        (bold_selectfiles, outputnode, [
            ('out_file', 'input_bold'),
            ]),
        (format_bold_buffer, bold_main_wf, [
            ("formatted_bold", "inputnode.bold"),
            ]),
        (resample_template_node, template_diagnosis, [
            ("resampled_template", "anat_template"),
            ]),
        (resample_template_node, commonspace_reg_wf, [
            ("resampled_template", "inputnode.atlas_anat"),
            ("resampled_mask", "inputnode.atlas_mask"),
            ]),
        (resample_template_node, bold_main_wf, [
            ("resampled_template", "inputnode.commonspace_ref"),
            ]),
        (resample_template_node, outputnode, [
            ("resampled_template", "commonspace_resampled_template"),
            ]),
        (merged_join_common_reg, commonspace_reg_wf, [
            ("file_list0", "inputnode.moving_image_list"),
            ("file_list1", "inputnode.moving_mask_list"),
            ]),
        (commonspace_reg_wf, bold_main_wf, [
            ("outputnode.native_to_commonspace_transform_list", "inputnode.native_to_commonspace_transform_list"),
            ("outputnode.native_to_commonspace_inverse_list", "inputnode.native_to_commonspace_inverse_list"),
            ("outputnode.commonspace_to_native_transform_list", "inputnode.commonspace_to_native_transform_list"),
            ("outputnode.commonspace_to_native_inverse_list", "inputnode.commonspace_to_native_inverse_list"),
            ]),
        (bold_main_wf, outputnode, [
            ("outputnode.bold_ref", "initial_bold_ref"),
            ("outputnode.corrected_EPI", "inho_cor_bold"),
            ("outputnode.native_brain_mask", "native_brain_mask"),
            ("outputnode.native_WM_mask", "native_WM_mask"),
            ("outputnode.native_CSF_mask", "native_CSF_mask"),
            ("outputnode.native_labels", "native_labels"),
            ("outputnode.confounds_csv", "confounds_csv"),
            ("outputnode.FD_voxelwise", "FD_voxelwise"),
            ("outputnode.pos_voxelwise", "pos_voxelwise"),
            ("outputnode.FD_csv", "FD_csv"),
            ('outputnode.bold_to_anat_affine', 'bold_to_anat_affine'),
            ('outputnode.bold_to_anat_warp', 'bold_to_anat_warp'),
            ('outputnode.bold_to_anat_inverse_warp', 'bold_to_anat_inverse_warp'),
            ("outputnode.output_warped_bold", "inho_cor_bold_warped2anat"),
            ("outputnode.native_bold", "native_bold"),
            ("outputnode.native_bold_ref", "native_bold_ref"),
            ("outputnode.commonspace_bold", "commonspace_bold"),
            ("outputnode.commonspace_mask", "commonspace_mask"),
            ("outputnode.commonspace_WM_mask", "commonspace_WM_mask"),
            ("outputnode.commonspace_CSF_mask", "commonspace_CSF_mask"),
            ("outputnode.commonspace_vascular_mask", "commonspace_vascular_mask"),
            ("outputnode.commonspace_labels", "commonspace_labels"),
            ]),
        (bold_main_wf, bold_inho_cor_diagnosis, [
            ("outputnode.bold_ref", "raw_img"),
            ("outputnode.init_denoise", "init_denoise"),
            ("outputnode.corrected_EPI", "final_denoise"),
            ("outputnode.denoise_mask", "warped_mask"),
            ]),
        (bold_selectfiles, bold_inho_cor_diagnosis,
         [("out_file", "name_source")]),
        (bold_main_wf, temporal_diagnosis, [
            ("outputnode.commonspace_bold", "bold_file"),
            ("outputnode.confounds_csv", "confounds_csv"),
            ("outputnode.FD_csv", "FD_csv"),
            ]),
        (bold_selectfiles, temporal_diagnosis,
         [("out_file", "name_source")]),
        (temporal_diagnosis, outputnode, [
            ("tSNR_filename", "tSNR_filename"),
            ("std_filename", "std_filename"),
            ]),
        ])

    if not opts.bold_only:
        run_split = pe.Node(niu.IdentityInterface(fields=['run', 'split_name']),
                            name="run_split")
        run_split.itersource = ('main_split', 'split_name')
        run_split.iterables = [('run', run_iter)]

        anat_selectfiles = pe.Node(BIDSDataGraber(bids_dir=data_dir_path, suffix=[
                                   'T2w', 'T1w']), name='anat_selectfiles')
        anat_selectfiles.inputs.run = None

        anat_convert_to_RAS_node = pe.Node(Function(input_names=['img_file'],
                                                    output_names=['RAS_file'],
                                                    function=convert_to_RAS),
                                           name='anat_convert_to_RAS')

        format_anat_buffer = pe.Node(niu.IdentityInterface(fields=['formatted_anat']),
                                            name="format_anat_buffer")

        if opts.anat_autobox: # apply AFNI's 3dAutobox
            anat_autobox = pe.Node(Autobox(padding=1, outputtype='NIFTI_GZ'),
                                    name="anat_autobox")
            workflow.connect([
                (anat_convert_to_RAS_node, anat_autobox, [
                    ('RAS_file', 'in_file'),
                    ]),
                (anat_autobox, format_anat_buffer, [
                    ('out_file', 'formatted_anat'),
                    ]),
                ])
        else:
            workflow.connect([
                (anat_convert_to_RAS_node, format_anat_buffer, [
                    ("RAS_file", "formatted_anat"),
                    ]),
                ])

        # setting anat preprocessing nodes
        anat_inho_cor_wf = init_inho_correction_wf(opts=opts, image_type='structural', name="anat_inho_cor_wf")

        workflow.connect([
            (main_split, run_split, [
                ("split_name", "split_name"),
                ]),
            (main_split, anat_selectfiles,
             [("scan_info", "scan_info")]),
            (run_split, bold_selectfiles, [
                ("run", "run"),
                ]),
            (anat_selectfiles, anat_convert_to_RAS_node,
             [("out_file", "img_file")]),
            (format_anat_buffer, anat_inho_cor_wf, [
                ("formatted_anat", "inputnode.target_img"),
                ("formatted_anat", "inputnode.name_source"),
                ]),
            (resample_template_node, anat_inho_cor_wf, [
                ("resampled_template", "inputnode.anat_ref"),
                ("resampled_mask", "inputnode.anat_mask"),
                ]),
            (anat_inho_cor_wf, bold_main_wf, [
                ("outputnode.corrected", "inputnode.coreg_anat"),
                ]),
            (commonspace_reg_wf, bold_main_wf, [
                ("outputnode.native_mask", "inputnode.coreg_mask"),
                ]),
            (EPI_target_buffer, bold_main_wf, [
                ("EPI_template", "inputnode.inho_cor_anat"),
                ("EPI_mask", "inputnode.inho_cor_mask"),
                ]),
            (anat_inho_cor_wf, source_join_common_reg, [
                ("outputnode.corrected", "file_list0"),
                ("outputnode.denoise_mask", "file_list1"),
                ]),
            (anat_inho_cor_wf, commonspace_reg_wf, [
                ("outputnode.corrected", "inputnode_iterable.iter_name"),
                ]),
            (anat_inho_cor_wf, EPI_target_buffer, [
                ("outputnode.corrected", "EPI_template"),
                ]),
            (commonspace_reg_wf, EPI_target_buffer, [
                ("outputnode.native_mask", 'EPI_mask'),
                ]),
            ])

        if not opts.anat_inho_cor_method=='disable':
            anat_inho_cor_diagnosis = pe.Node(Function(input_names=['raw_img','init_denoise','warped_mask','final_denoise', 'name_source', 'out_dir'],
                                               function=preprocess_visual_QC.inho_cor_diagnosis),
                                      name='anat_inho_cor_diagnosis')
            anat_inho_cor_diagnosis.inputs.out_dir = output_folder+'/preprocess_QC_report/anat_inho_cor/'

            workflow.connect([
                (format_anat_buffer, anat_inho_cor_diagnosis, [
                    ("formatted_anat", "raw_img"),
                    ]),
                (anat_selectfiles, anat_inho_cor_diagnosis, [
                    ("out_file", "name_source"),
                    ]),
                (anat_inho_cor_wf, anat_inho_cor_diagnosis, [
                    ("outputnode.init_denoise", "init_denoise"),
                    ("outputnode.corrected", "final_denoise"),
                    ("outputnode.denoise_mask", "warped_mask"),
                    ]),
                ])

    else:
        inho_cor_bold_main_wf = init_bold_main_wf(
            inho_cor_only=True, name='inho_cor_bold_main_wf', opts=opts)

        workflow.connect([
            (format_bold_buffer, inho_cor_bold_main_wf, [
                ("formatted_bold", "inputnode.bold"),
                ]),
            (EPI_target_buffer, inho_cor_bold_main_wf, [
                ("EPI_template", "inputnode.inho_cor_anat"),
                ("EPI_mask", "inputnode.inho_cor_mask"),
                ]),
            (inho_cor_bold_main_wf, bold_main_wf, [
                ("transitionnode.bold_file", "transitionnode.bold_file"),
                ("transitionnode.bold_ref", "transitionnode.bold_ref"),
                ("transitionnode.init_denoise", "transitionnode.init_denoise"),
                ("transitionnode.denoise_mask", "transitionnode.denoise_mask"),
                ("transitionnode.corrected_EPI", "transitionnode.corrected_EPI"),
                ]),
            (inho_cor_bold_main_wf, source_join_common_reg, [
                ("transitionnode.corrected_EPI", "file_list0"),
                ("transitionnode.denoise_mask", "file_list1"),
                ]),
            (inho_cor_bold_main_wf, commonspace_reg_wf, [
                ("transitionnode.corrected_EPI", "inputnode_iterable.iter_name"),
                ]),
            (resample_template_node, EPI_target_buffer, [
                ("resampled_template", "EPI_template"),
                ("resampled_mask", "EPI_mask"),
                ]),
            ])

    if not opts.bold_only:
        PlotOverlap_EPI2Anat_node = pe.Node(
            preprocess_visual_QC.PlotOverlap(), name='PlotOverlap_EPI2Anat')
        PlotOverlap_EPI2Anat_node.inputs.out_dir = output_folder+'/preprocess_QC_report/EPI2Anat'
        workflow.connect([
            (bold_selectfiles, PlotOverlap_EPI2Anat_node,
             [("out_file", "name_source")]),
            (anat_inho_cor_wf, PlotOverlap_EPI2Anat_node,
             [("outputnode.corrected", "fixed")]),
            (outputnode, PlotOverlap_EPI2Anat_node, [
                ("inho_cor_bold_warped2anat", "moving"),  # warped EPI to anat
                ]),
            ])

    # fill the datasinks
    workflow.connect([
        (bold_selectfiles, bold_datasink, [
            ("out_file", "input_bold"),
            ]),
        (outputnode, confounds_datasink, [
            ("confounds_csv", "confounds_csv"),  # confounds file
            ("FD_voxelwise", "FD_voxelwise"),
            ("pos_voxelwise", "pos_voxelwise"),
            ("FD_csv", "FD_csv"),
            ]),
        (outputnode, bold_datasink, [
            ("initial_bold_ref", "initial_bold_ref"),  # inspect initial bold ref
            ("inho_cor_bold", "inho_cor_bold"),  # inspect bias correction
            ("native_brain_mask", "native_brain_mask"),  # get the EPI labels
            ("native_WM_mask", "native_WM_mask"),  # get the EPI labels
            ("native_CSF_mask", "native_CSF_mask"),  # get the EPI labels
            ("native_labels", "native_labels"),  # get the EPI labels
            # warped EPI to anat
            ("inho_cor_bold_warped2anat", "inho_cor_bold_warped2anat"),
            # resampled EPI after motion realignment and SDC
            ("native_bold", "native_bold"),
            # resampled EPI after motion realignment and SDC
            ("native_bold_ref", "native_bold_ref"),
            # resampled EPI after motion realignment and SDC
            ("commonspace_bold", "commonspace_bold"),
            ("commonspace_mask", "commonspace_mask"),
            ("commonspace_WM_mask", "commonspace_WM_mask"),
            ("commonspace_CSF_mask", "commonspace_CSF_mask"),
            ("commonspace_vascular_mask", "commonspace_vascular_mask"),
            ("commonspace_labels", "commonspace_labels"),
            ("tSNR_filename", "tSNR_map_preprocess"),
            ("std_filename", "std_map_preprocess"),
            ("commonspace_resampled_template", "commonspace_resampled_template"),
            ]),
        ])

    if not opts.bold_only:
        workflow.connect([
            (anat_inho_cor_wf, anat_datasink, [
                ("outputnode.corrected", "anat_preproc"),
                ]),
            (outputnode, transforms_datasink, [
                ('bold_to_anat_affine', 'bold_to_anat_affine'),
                ('bold_to_anat_warp', 'bold_to_anat_warp'),
                ('bold_to_anat_inverse_warp', 'bold_to_anat_inverse_warp'),
                ]),
            ])

    return workflow


def join_iterables(workflow, joinsource_list, node_prefix, num_inputs=1):

    field_list=[]
    for j in range(num_inputs):
        field_list.append(f'file_list{j}')

    i=0
    for joinsource in joinsource_list:
        joinnode = pe.JoinNode(niu.IdentityInterface(fields=field_list),
                                            name=f"{node_prefix}_{joinsource}_joinnode",
                                            joinsource=joinsource,
                                            joinfield=field_list)
        if i==0:
            source_join = joinnode
        else:
            for field in field_list:
                workflow.connect([
                    (joinnode_prev, joinnode, [
                        (field, field),
                        ]),
                    ])

        joinnode_prev = joinnode
        i+=1

    merged_join = joinnode

    return workflow, source_join, merged_join
