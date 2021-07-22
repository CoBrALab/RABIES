import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
from nipype.interfaces.afni import Autobox
from .preprocess_pkg.inho_correction import init_inho_correction_wf
from .preprocess_pkg.commonspace_reg import init_commonspace_reg_wf,GenerateTemplate
from .preprocess_pkg.bold_main_wf import init_bold_main_wf
from .preprocess_pkg.utils import BIDSDataGraber, prep_bids_iter, convert_to_RAS
from .preprocess_pkg import preprocess_visual_QC


def init_main_wf(data_dir_path, output_folder, opts, cr_opts=None, analysis_opts=None, name='main_wf'):
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
            parser options for confound_regression
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
        affine_bold2anat
            affine transform from the EPI space to the anatomical space
        warp_bold2anat
            non-linear transform from the EPI space to the anatomical space
        inverse_warp_bold2anat
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
        fields=['input_bold', 'commonspace_resampled_template', 'anat_preproc', 'initial_bold_ref', 'inho_cor_bold', 'affine_bold2anat',
                'warp_bold2anat', 'inverse_warp_bold2anat', 'inho_cor_bold_warped2anat', 'native_bold', 'native_bold_ref', 'confounds_csv',
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

    from bids.layout import BIDSLayout
    layout = BIDSLayout(data_dir_path, validate=False)
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
    from rabies.preprocess_pkg.utils import resample_template
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
    if opts.local_threads < num_scan:
        num_scan = opts.local_threads

    EPI_target_buffer = pe.Node(niu.IdentityInterface(fields=['EPI_template', 'EPI_mask']),
                                        name="EPI_target_buffer")

    if opts.fast_commonspace:
        # if fast commonspace, then the inputs iterables are not merged
        source_join_common_reg = pe.Node(niu.IdentityInterface(fields=['file_list0', 'file_list1']),
                                            name="fast_commonreg_buffer")
        merged_join_common_reg = source_join_common_reg
    else:
        workflow, source_join_common_reg, merged_join_common_reg = join_iterables(workflow=workflow, joinsource_list=['main_split'], node_prefix='commonspace_reg', num_inputs=2)

    commonspace_reg_wf = init_commonspace_reg_wf(opts=opts, output_folder=output_folder, transforms_datasink=transforms_datasink, num_scan=num_scan, name='commonspace_reg_wf')

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
            ("file_list0", "inputnode.moving_image"),
            ("file_list1", "inputnode.moving_mask"),
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
            ('outputnode.affine_bold2anat', 'affine_bold2anat'),
            ('outputnode.warp_bold2anat', 'warp_bold2anat'),
            ('outputnode.inverse_warp_bold2anat', 'inverse_warp_bold2anat'),
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
        anat_inho_cor_wf = init_inho_correction_wf(opts=opts, image_type='structural', inho_cor_method=opts.anat_inho_cor_method, name="anat_inho_cor_wf")

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
            ])

        if not opts.robust_bold_inho_cor:
            workflow.connect([
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
            ])

        if not opts.robust_bold_inho_cor:
            workflow.connect([
                (resample_template_node, EPI_target_buffer, [
                    ("resampled_template", "EPI_template"),
                    ("resampled_mask", "EPI_mask"),
                    ]),
                ])

    if opts.robust_bold_inho_cor:
        if opts.fast_commonspace:
            raise ValueError("--fast_commonspace prevents to gain any benefit from --robust_bold_inho_cor")

        inho_cor_robust_wf = init_bold_main_wf(
            inho_cor_only=True, name='inho_cor_robust_wf', opts=opts)

        generate_EPI_template_wf = init_commonspace_reg_wf(opts=opts, output_folder=output_folder, transforms_datasink=transforms_datasink,
                                                           num_scan=num_scan, name='generate_EPI_template_wf')

        EPI_template_masking = pe.Node(Function(input_names=['fixed_mask', 'moving_image', 'inverse_warp', 'affine'],
                                        output_names=['new_mask'],
                                        function=transform_mask),
                               name='EPI_template_masking')
        EPI_template_masking.inputs.fixed_mask = str(opts.brain_mask)

        if not opts.bold_only:
            workflow, source_join_robust_cor, merged_join_robust_cor = join_iterables(workflow=workflow, joinsource_list=['run_split','main_split'], node_prefix='robust_bold_inho_cor', num_inputs=2)

            workflow.connect([
                (anat_inho_cor_wf, inho_cor_robust_wf, [
                    ("outputnode.corrected", "inputnode.inho_cor_anat"),
                    ]),
                (commonspace_reg_wf, inho_cor_robust_wf, [
                    ("outputnode.native_mask", 'inputnode.inho_cor_mask'),
                    ]),
                ])
        else:
            workflow, source_join_robust_cor, merged_join_robust_cor = join_iterables(workflow=workflow, joinsource_list=['main_split'], node_prefix='robust_bold_inho_cor', num_inputs=2)
            workflow.connect([
                (resample_template_node, inho_cor_robust_wf, [
                    ("resampled_template", "inputnode.inho_cor_anat"),
                    ("resampled_mask", "inputnode.inho_cor_mask"),
                    ]),
                ])

        workflow.connect([
            (format_bold_buffer, inho_cor_robust_wf, [
                ("formatted_bold", "inputnode.bold"),
                ]),
            (inho_cor_robust_wf, source_join_robust_cor, [
                ("transitionnode.corrected_EPI", "file_list0"),
                ("transitionnode.denoise_mask", "file_list1"),
                ]),
            (inho_cor_robust_wf, generate_EPI_template_wf, [
                ("transitionnode.corrected_EPI", "inputnode_iterable.iter_name"),
                ]),
            (resample_template_node, generate_EPI_template_wf, [
                ("resampled_template", "inputnode.atlas_anat"),
                ]),
            (merged_join_robust_cor, generate_EPI_template_wf, [
                ("file_list0", "inputnode.moving_image"),
                ("file_list1", "inputnode.moving_mask"),
                ]),
            (generate_EPI_template_wf, EPI_template_masking, [
                ("outputnode.unbiased_template", "moving_image"),
                ]),
            (generate_EPI_template_wf, EPI_template_masking, [
                ("outputnode.to_atlas_affine", "affine"),
                ("outputnode.to_atlas_inverse_warp", "inverse_warp"),
                ]),
            (generate_EPI_template_wf, EPI_target_buffer, [
                ("outputnode.unbiased_template", "EPI_template"),
                ]),
            (EPI_template_masking, EPI_target_buffer, [
                ("new_mask", "EPI_mask"),
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


    # Integrate confound regression
    if cr_opts is not None:
        workflow, confound_regression_wf = integrate_confound_regression(
            workflow, outputnode, cr_opts, bold_only=opts.bold_only)

        # Integrate analysis
        if analysis_opts is not None:
            workflow = integrate_analysis(
                workflow, outputnode, confound_regression_wf, analysis_opts, opts.bold_only, cr_opts.commonspace_analysis, bold_scan_list, opts)

    elif opts.rabies_step == 'preprocess':
        # only fill datasinks if the related workflow is running

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
                ]),
            ])

        if not opts.bold_only:
            workflow.connect([
                (anat_inho_cor_wf, anat_datasink, [
                    ("outputnode.corrected", "anat_preproc"),
                 ]),
                (outputnode, transforms_datasink, [
                    ('affine_bold2anat', 'affine_bold2anat'),
                    ('warp_bold2anat', 'warp_bold2anat'),
                    ('inverse_warp_bold2anat', 'inverse_warp_bold2anat'),
                    ]),
                ])

    return workflow


def integrate_confound_regression(workflow, outputnode, cr_opts, bold_only):
    cr_output = os.path.abspath(str(cr_opts.output_dir))

    from rabies.conf_reg_pkg.confound_regression import init_confound_regression_wf
    confound_regression_wf = init_confound_regression_wf(cr_opts=cr_opts, name=cr_opts.output_name)

    workflow.connect([
        (outputnode, confound_regression_wf, [
            ("confounds_csv", "inputnode.confounds_file"),  # confounds file
            ("FD_csv", "inputnode.FD_file"),
            ]),
        ])

    if bold_only and not cr_opts.commonspace_analysis:
        raise ValueError(
            'Must select --commonspace option for running confound regression on outputs from --bold_only.')

    if cr_opts.commonspace_analysis:
        workflow.connect([
            (outputnode, confound_regression_wf, [
                ("commonspace_bold", "inputnode.bold_file"),
                ("commonspace_mask", "inputnode.brain_mask"),
                ("commonspace_CSF_mask", "inputnode.csf_mask"),
                ]),
            ])
    else:
        workflow.connect([
            (outputnode, confound_regression_wf, [
                ("native_bold", "inputnode.bold_file"),
                ("native_brain_mask", "inputnode.brain_mask"),
                ("native_CSF_mask", "inputnode.csf_mask"),
                ]),
            ])

    if cr_opts.rabies_step == 'confound_regression':
        confound_regression_datasink = pe.Node(DataSink(base_directory=cr_output,
                                                        container=cr_opts.output_name+"_datasink"),
                                               name=cr_opts.output_name+"_datasink")
        workflow.connect([
            (confound_regression_wf, confound_regression_datasink, [
                ("outputnode.cleaned_path", "cleaned_timeseries"),
                ]),
            ])
        if cr_opts.run_aroma:
            workflow.connect([
                (confound_regression_wf, confound_regression_datasink, [
                    ("outputnode.aroma_out", "aroma_out"),
                    ]),
                ])

    return workflow, confound_regression_wf


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


def transit_iterables(workflow, prep_dict_node, scan_list, bold_only, bold_scan_list, node_prefix=''):
    # this function adds nodes to the workflow which first joins the iterables from preprocessing
    # and then creates a new iterable list than can be customized for analysis

    joinnode_main = pe.JoinNode(niu.IdentityInterface(fields=['file_list']),
                                         name=node_prefix+'_joinnode_main',
                                         joinsource='main_split',
                                         joinfield=['file_list'])

    scan_split_name = get_iterable_scan_list(scan_list, bold_scan_list)

    analysis_split = pe.Node(niu.IdentityInterface(fields=['scan_split_name']),
                         name=node_prefix+"_split")
    analysis_split.iterables = [('scan_split_name', scan_split_name)]

    find_iterable_node = pe.Node(Function(input_names=['file_list', 'scan_split_name'],
                                           output_names=[
                                               'file'],
                                       function=find_iterable),
                              name=node_prefix+"_find_iterable")

    workflow.connect([
        (analysis_split, find_iterable_node, [
            ("scan_split_name", "scan_split_name"),
            ]),
        (joinnode_main, find_iterable_node, [
            ("file_list", "file_list"),
            ]),
        ])

    if not bold_only:
        joinnode_run = pe.JoinNode(niu.IdentityInterface(fields=['file_list']),
                                            name=node_prefix+"_joinnode_run",
                                            joinsource='run_split',
                                            joinfield=['file_list'])

        workflow.connect([
            (joinnode_run, joinnode_main, [
                ("file_list", "file_list"),
                ]),
            (prep_dict_node, joinnode_run, [
                ("prep_dict", "file_list"),
                ]),
            ])
    else:
        workflow.connect([
            (prep_dict_node, joinnode_main, [
                ("prep_dict", "file_list"),
                ]),
            ])

    return workflow,find_iterable_node, joinnode_main,analysis_split


def integrate_analysis(workflow, outputnode, confound_regression_wf, analysis_opts, bold_only, commonspace_bold, bold_scan_list, opts):
    analysis_output = os.path.abspath(str(analysis_opts.output_dir))

    from rabies.analysis_pkg.analysis_wf import init_analysis_wf
    analysis_wf = init_analysis_wf(
        opts=analysis_opts, commonspace_cr=commonspace_bold, seed_list=analysis_opts.seed_list, name=analysis_opts.output_name)

    analysis_datasink = pe.Node(DataSink(base_directory=analysis_output,
                                         container=analysis_opts.output_name+"_datasink"),
                                name=analysis_opts.output_name+"_datasink")

    def prep_dict(bold_file, CR_data_dict, VE_file, mask_file, WM_mask_file, CSF_mask_file, atlas_file, preprocess_anat_template, name_source):
        return {'bold_file':bold_file, 'CR_data_dict':CR_data_dict, 'VE_file':VE_file, 'mask_file':mask_file, 'WM_mask_file':WM_mask_file, 'CSF_mask_file':CSF_mask_file, 'atlas_file':atlas_file, 'preprocess_anat_template':preprocess_anat_template, 'name_source':name_source}
    prep_dict_node = pe.Node(Function(input_names=['bold_file', 'CR_data_dict', 'VE_file', 'mask_file', 'WM_mask_file', 'CSF_mask_file', 'atlas_file', 'preprocess_anat_template', 'name_source'],
                                           output_names=[
                                               'prep_dict'],
                                       function=prep_dict),
                              name=analysis_opts.output_name+'_prep_dict')

    def read_dict(prep_dict):
        return prep_dict['bold_file'],prep_dict['mask_file'],prep_dict['atlas_file']
    read_dict_node = pe.Node(Function(input_names=['prep_dict'],
                                           output_names=['bold_file', 'mask_file', 'atlas_file'],
                                       function=read_dict),
                              name=analysis_opts.output_name+'_read_dict')

    workflow,find_iterable_node, joinnode_main,analysis_split = transit_iterables(workflow, prep_dict_node, analysis_opts.scan_list, bold_only, bold_scan_list, node_prefix=analysis_opts.output_name)

    analysis_split_joinnode = pe.JoinNode(niu.IdentityInterface(fields=['file_list', 'mask_file']),
                                         name=analysis_opts.output_name+'_split_joinnode',
                                         joinsource=analysis_split.name,
                                         joinfield=['file_list'])

    if commonspace_bold or bold_only:
        workflow.connect([
            (outputnode, prep_dict_node, [
                ("commonspace_mask", "mask_file"),
                ("commonspace_labels", "atlas_file"),
                ("commonspace_WM_mask", "WM_mask_file"),
                ("commonspace_CSF_mask", "CSF_mask_file"),
                ("commonspace_resampled_template", "preprocess_anat_template"),
                ]),
            ])
    else:
        workflow.connect([
            (outputnode, prep_dict_node, [
                ("native_brain_mask", "subject_inputnode.mask_file"),
                ("native_labels", "subject_inputnode.atlas_file"),
                ]),
            ])

    workflow.connect([
        (outputnode, prep_dict_node, [
            ("input_bold", "name_source"),
            ]),
        (confound_regression_wf, prep_dict_node, [
            ("outputnode.cleaned_path", "bold_file"),
            ("outputnode.CR_data_dict", "CR_data_dict"),
            ("outputnode.VE_file", "VE_file"),
            ]),
        (find_iterable_node, read_dict_node, [
            ("file", "prep_dict"),
            ]),
        (analysis_split_joinnode, analysis_wf, [
            ("file_list", "group_inputnode.bold_file_list"),
            ("mask_file", "group_inputnode.commonspace_mask"),
            ]),
        (read_dict_node, analysis_split_joinnode, [
            ("mask_file", "mask_file"),
            ]),
        (read_dict_node, analysis_split_joinnode, [
            ("bold_file", "file_list"),
            ]),
        (read_dict_node, analysis_wf, [
            ("bold_file", "subject_inputnode.bold_file"),
            ("mask_file", "subject_inputnode.mask_file"),
            ("atlas_file", "subject_inputnode.atlas_file"),
            ]),
        (analysis_wf, analysis_datasink, [
            ("outputnode.group_ICA_dir", "group_ICA_dir"),
            ("outputnode.matrix_data_file", "matrix_data_file"),
            ("outputnode.matrix_fig", "matrix_fig"),
            ("outputnode.corr_map_file", "seed_correlation_maps"),
            ]),
        ])

    if analysis_opts.data_diagnosis:
        if not (commonspace_bold or bold_only):
            raise ValueError("--commonspace_analysis outputs are currently required for running data_diagnosis")

        if os.path.basename(opts.anat_template)=='DSURQE_40micron_average.nii.gz':
            DSURQE_regions=True
        else:
            DSURQE_regions=False

        from rabies.analysis_pkg.data_diagnosis import ScanDiagnosis, PrepMasks, DatasetDiagnosis, temporal_external_formating, spatial_external_formating
        ScanDiagnosis_node = pe.Node(ScanDiagnosis(prior_bold_idx=analysis_opts.prior_bold_idx,
            prior_confound_idx=analysis_opts.prior_confound_idx,
                dual_ICA = analysis_opts.dual_ICA, DSURQE_regions=DSURQE_regions),
            name=analysis_opts.output_name+'_ScanDiagnosis')

        PrepMasks_node = pe.Node(PrepMasks(prior_maps=os.path.abspath(str(analysis_opts.prior_maps)), DSURQE_regions=DSURQE_regions),
            name=analysis_opts.output_name+'_PrepMasks')

        DatasetDiagnosis_node = pe.Node(DatasetDiagnosis(),
            name=analysis_opts.output_name+'_DatasetDiagnosis')

        temporal_external_formating_node = pe.Node(Function(input_names=['temporal_info', 'file_dict'],
                                               output_names=[
                                                   'temporal_info_csv', 'dual_regression_timecourse_csv', 'dual_ICA_timecourse_csv'],
                                           function=temporal_external_formating),
                                  name=analysis_opts.output_name+'_temporal_external_formating')

        spatial_external_formating_node = pe.Node(Function(input_names=['spatial_info', 'file_dict'],
                                               output_names=[
                                                   'std_filename', 'GS_corr_filename', 'DVARS_corr_filename', 'FD_corr_filename', 'DR_maps_filename', 'dual_ICA_filename'],
                                           function=spatial_external_formating),
                                  name=analysis_opts.output_name+'_spatial_external_formating')

        data_diagnosis_split_joinnode = pe.JoinNode(niu.IdentityInterface(fields=['spatial_info_list']),
                                             name=analysis_opts.output_name+'_diagnosis_split_joinnode',
                                             joinsource=analysis_split.name,
                                             joinfield=['spatial_info_list'])

        workflow.connect([
            (joinnode_main, PrepMasks_node, [
                ("file_list", "mask_dict_list"),
                ]),
            (PrepMasks_node, ScanDiagnosis_node, [
                ("mask_file_dict", "mask_file_dict"),
                ]),
            (find_iterable_node, ScanDiagnosis_node, [
                ("file", "file_dict"),
                ]),
            (ScanDiagnosis_node, data_diagnosis_split_joinnode, [
                ("spatial_info", "spatial_info_list"),
                ]),
            (data_diagnosis_split_joinnode, DatasetDiagnosis_node, [
                ("spatial_info_list", "spatial_info_list"),
                ]),
            (PrepMasks_node, DatasetDiagnosis_node, [
                ("mask_file_dict", "mask_file_dict"),
                ]),
            (find_iterable_node, temporal_external_formating_node, [
                ("file", "file_dict"),
                ]),
            (ScanDiagnosis_node, temporal_external_formating_node, [
                ("temporal_info", "temporal_info"),
                ]),
            (find_iterable_node, spatial_external_formating_node, [
                ("file", "file_dict"),
                ]),
            (ScanDiagnosis_node, spatial_external_formating_node, [
                ("spatial_info", "spatial_info"),
                ]),
            ])

        workflow.connect([
            (ScanDiagnosis_node, analysis_datasink, [
                ("figure_temporal_diagnosis", "figure_temporal_diagnosis"),
                ("figure_spatial_diagnosis", "figure_spatial_diagnosis"),
                ("VE_file", "VE_file"),
                ]),
            (DatasetDiagnosis_node, analysis_datasink, [
                ("figure_dataset_diagnosis", "figure_dataset_diagnosis"),
                ]),
            (temporal_external_formating_node, analysis_datasink, [
                ("temporal_info_csv", "temporal_info_csv"),
                ("dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
                ]),
            (spatial_external_formating_node, analysis_datasink, [
                ("std_filename", "temporal_std_nii"),
                ("GS_corr_filename", "GS_corr_nii"),
                ("DVARS_corr_filename", "DVARS_corr_nii"),
                ("FD_corr_filename", "FD_corr_nii"),
                ("DR_maps_filename", "dual_regression_nii"),
                ]),
            ])

        if analysis_opts.dual_ICA>0:
            workflow.connect([
                (spatial_external_formating_node, analysis_datasink, [
                    ("dual_ICA_filename", "dual_ICA_nii"),
                    ]),
                (temporal_external_formating_node, analysis_datasink, [
                    ("dual_ICA_timecourse_csv", "dual_ICA_timecourse_csv"),
                    ]),
                ])
    else:
        workflow.connect([
            (analysis_wf, analysis_datasink, [
                ("outputnode.DR_nii_file", "dual_regression_nii"),
                ("outputnode.dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
                ("outputnode.dual_ICA_timecourse_csv", "dual_ICA_timecourse_csv"),
                ("outputnode.dual_ICA_filename", "dual_ICA_filename"),
                ]),
            ])


    return workflow


def integrate_data_diagnosis(workflow, outputnode, confound_regression_wf, data_diagnosis_opts, bold_only, commonspace_bold, bold_scan_list):

    data_diagnosis_output = os.path.abspath(str(data_diagnosis_opts.output_dir))

    from rabies.analysis_pkg.data_diagnosis import ScanDiagnosis, PrepMasks, DatasetDiagnosis, temporal_external_formating, spatial_external_formating
    ScanDiagnosis_node = pe.Node(ScanDiagnosis(prior_bold_idx=data_diagnosis_opts.prior_bold_idx,
        prior_confound_idx=data_diagnosis_opts.prior_confound_idx,
            dual_ICA = data_diagnosis_opts.dual_ICA, DSURQE_regions=data_diagnosis_opts.DSURQE_regions),
        name=data_diagnosis_opts.output_name+'_ScanDiagnosis')

    PrepMasks_node = pe.Node(PrepMasks(prior_maps=os.path.abspath(str(data_diagnosis_opts.prior_maps)), DSURQE_regions=data_diagnosis_opts.DSURQE_regions),
        name=data_diagnosis_opts.output_name+'_PrepMasks')

    DatasetDiagnosis_node = pe.Node(DatasetDiagnosis(),
        name=data_diagnosis_opts.output_name+'_DatasetDiagnosis')

    temporal_external_formating_node = pe.Node(Function(input_names=['temporal_info', 'file_dict'],
                                           output_names=[
                                               'temporal_info_csv', 'dual_regression_timecourse_csv', 'dual_ICA_timecourse_csv'],
                                       function=temporal_external_formating),
                              name=data_diagnosis_opts.output_name+'_temporal_external_formating')

    spatial_external_formating_node = pe.Node(Function(input_names=['spatial_info', 'file_dict'],
                                           output_names=[
                                               'std_filename', 'GS_corr_filename', 'DVARS_corr_filename', 'FD_corr_filename', 'DR_maps_filename', 'dual_ICA_filename'],
                                       function=spatial_external_formating),
                              name=data_diagnosis_opts.output_name+'_spatial_external_formating')

    def prep_dict(bold_file, CR_data_dict, VE_file, brain_mask_file, WM_mask_file, CSF_mask_file, preprocess_anat_template, name_source):
        return {'bold_file':bold_file, 'CR_data_dict':CR_data_dict, 'VE_file':VE_file, 'brain_mask_file':brain_mask_file, 'WM_mask_file':WM_mask_file, 'CSF_mask_file':CSF_mask_file, 'preprocess_anat_template':preprocess_anat_template, 'name_source':name_source}
    prep_dict_node = pe.Node(Function(input_names=['bold_file', 'CR_data_dict', 'VE_file', 'brain_mask_file', 'WM_mask_file', 'CSF_mask_file', 'preprocess_anat_template', 'name_source'],
                                           output_names=[
                                               'prep_dict'],
                                       function=prep_dict),
                              name=data_diagnosis_opts.output_name+'_prep_dict')

    workflow, find_iterable_node, joinnode_main, analysis_split = transit_iterables(workflow, prep_dict_node, data_diagnosis_opts.scan_list, bold_only, bold_scan_list, node_prefix=data_diagnosis_opts.output_name)

    data_diagnosis_split_joinnode = pe.JoinNode(niu.IdentityInterface(fields=['spatial_info_list']),
                                         name=data_diagnosis_opts.output_name+'_split_joinnode',
                                         joinsource=analysis_split.name,
                                         joinfield=['spatial_info_list'])

    workflow.connect([
        (outputnode, prep_dict_node, [
            ("commonspace_mask", "brain_mask_file"),
            ("commonspace_WM_mask", "WM_mask_file"),
            ("commonspace_CSF_mask", "CSF_mask_file"),
            ("commonspace_resampled_template", "preprocess_anat_template"),
            ("input_bold", "name_source"),
            ]),
        (confound_regression_wf, prep_dict_node, [
            ("outputnode.cleaned_path", "bold_file"),
            ("outputnode.CR_data_dict", "CR_data_dict"),
            ("outputnode.VE_file", "VE_file"),
            ]),
        (joinnode_main, PrepMasks_node, [
            ("file_list", "mask_dict_list"),
            ]),
        (PrepMasks_node, ScanDiagnosis_node, [
            ("mask_file_dict", "mask_file_dict"),
            ]),
        (find_iterable_node, ScanDiagnosis_node, [
            ("file", "file_dict"),
            ]),
        (ScanDiagnosis_node, data_diagnosis_split_joinnode, [
            ("spatial_info", "spatial_info_list"),
            ]),
        (data_diagnosis_split_joinnode, DatasetDiagnosis_node, [
            ("spatial_info_list", "spatial_info_list"),
            ]),
        (PrepMasks_node, DatasetDiagnosis_node, [
            ("mask_file_dict", "mask_file_dict"),
            ]),
        (find_iterable_node, temporal_external_formating_node, [
            ("file", "file_dict"),
            ]),
        (ScanDiagnosis_node, temporal_external_formating_node, [
            ("temporal_info", "temporal_info"),
            ]),
        (find_iterable_node, spatial_external_formating_node, [
            ("file", "file_dict"),
            ]),
        (ScanDiagnosis_node, spatial_external_formating_node, [
            ("spatial_info", "spatial_info"),
            ]),
        ])

    data_diagnosis_datasink = pe.Node(DataSink(base_directory=data_diagnosis_output,
                                     container=data_diagnosis_opts.output_name+"_datasink"),
                            name=data_diagnosis_opts.output_name+"_datasink")
    workflow.connect([
        (ScanDiagnosis_node, data_diagnosis_datasink, [
            ("figure_temporal_diagnosis", "figure_temporal_diagnosis"),
            ("figure_spatial_diagnosis", "figure_spatial_diagnosis"),
            ("VE_file", "VE_file"),
            ]),
        (DatasetDiagnosis_node, data_diagnosis_datasink, [
            ("figure_dataset_diagnosis", "figure_dataset_diagnosis"),
            ]),
        (temporal_external_formating_node, data_diagnosis_datasink, [
            ("temporal_info_csv", "temporal_info_csv"),
            ("dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
            ]),
        (spatial_external_formating_node, data_diagnosis_datasink, [
            ("std_filename", "temporal_std_nii"),
            ("GS_corr_filename", "GS_corr_nii"),
            ("DVARS_corr_filename", "DVARS_corr_nii"),
            ("FD_corr_filename", "FD_corr_nii"),
            ("DR_maps_filename", "dual_regression_nii"),
            ]),
        ])

    if data_diagnosis_opts.dual_ICA>0:
        workflow.connect([
            (spatial_external_formating_node, data_diagnosis_datasink, [
                ("dual_ICA_filename", "dual_ICA_nii"),
                ]),
            (temporal_external_formating_node, data_diagnosis_datasink, [
                ("dual_ICA_timecourse_csv", "dual_ICA_timecourse_csv"),
                ]),
            ])

    return workflow


def transform_mask(fixed_mask, moving_image, inverse_warp, affine):
    # function to transform atlas masks to individual anatomical scans
    import os
    from rabies.preprocess_pkg.utils import run_command
    cwd = os.getcwd()

    import pathlib  # Better path manipulation
    filename_template = pathlib.Path(moving_image).name.rsplit(".nii")[0]

    new_mask = f'{cwd}/{filename_template}_mask.nii.gz'
    if inverse_warp=='NULL':
        transforms = [affine]
        inverses = [1]
    else:
        transforms = [affine,inverse_warp]
        inverses = [1,0]

    from rabies.preprocess_pkg.utils import exec_applyTransforms
    exec_applyTransforms(transforms, inverses, fixed_mask, moving_image, new_mask, mask=True)
    return new_mask


def get_iterable_scan_list(scan_list, bold_scan_list):
    # prep the subset of scans on which the analysis will be run
    import numpy as np
    import pandas as pd
    scan_split_name=[]
    if os.path.isfile(os.path.abspath(scan_list[0])):
        if '.nii' in pathlib.Path(scan_list[0]).name:
            for scan in scan_list:
                scan_split_name.append(pathlib.Path(scan).name.rsplit(".nii")[0])
        else:
            # read the file as a .txt
            scan_list = np.array(pd.read_csv(os.path.abspath(scan_list[0]), header=None)).flatten()
            for scan in scan_list:
                scan_split_name.append(pathlib.Path(scan).name.rsplit(".nii")[0])

    elif scan_list[0]=='all':
        for scan in bold_scan_list:
            scan_split_name.append(pathlib.Path(scan).name.rsplit(".nii")[0])
    else:
        raise ValueError("The scan_list input had improper format.")
    return scan_split_name


def find_iterable(file_list, scan_split_name):
    # find the proper iterable index on the file_list based on the
    # correspondence between the input scan names and the iterable split name
    from rabies.preprocess_pkg.utils import flatten_list
    file_list = flatten_list(list(file_list))
    for file in file_list:
        if scan_split_name in file['name_source']:
            import os
            if os.path.basename(file['bold_file']) == 'empty.nii.gz':
                raise ValueError(f"FD_censoring and/or DVARS_censoring during confound regression resulted in an empty file for scan {file['name_source']}. \
                                You will have to specify a list of scans that are not empty with --scan_list option.")
            return file
    raise ValueError(f"No matching file was found for {scan_split_name}")
