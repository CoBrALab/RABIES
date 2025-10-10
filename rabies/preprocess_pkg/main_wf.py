import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
from .inho_correction import init_inho_correction_wf
from .commonspace_reg import init_commonspace_reg_wf,inherit_unbiased_files
from .bold_main_wf import init_bold_main_wf
from .utils import BIDSDataGraber, prep_bids_iter, convert_to_RAS, correct_oblique_affine, convert_3dWarp, apply_autobox, resample_template,log_transform_nii
from . import preprocess_visual_QC
from .prep_input import init_prep_input_wf,match_iterables

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
        motion_params_csv
            .csv file with measured motion timecourses, used as regressors for confound
            correction: 6 rigid body motion parameters + their first temporal derivate 
            + the 12 parameters squared (24 motion parameters)
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
                'bold_to_anat_warp', 'bold_to_anat_inverse_warp', 'inho_cor_bold_warped2anat', 'native_bold', 'native_bold_ref', 'motion_params_csv',
                'FD_voxelwise', 'pos_voxelwise', 'FD_csv', 'native_brain_mask', 'native_WM_mask', 'native_CSF_mask', 'native_vascular_mask', 'native_labels',
                'commonspace_bold', 'commonspace_mask', 'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_vascular_mask',
                'commonspace_labels', 'std_filename', 'tSNR_filename', 'raw_brain_mask']),
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

    motion_datasink = pe.Node(DataSink(base_directory=output_folder,
                                          container="motion_datasink"),
                                 name="motion_datasink")

    import bids
    bids.config.set_option('extension_initial_dot', True)
    try:
        layout = bids.layout.BIDSLayout(data_dir_path, validate=True)
    except Exception as e:
        from nipype import logging
        log = logging.getLogger('nipype.workflow')
        log.warning(f"The BIDS compliance failed: {e} \n\nRABIES will run anyway; double-check that the right files were picked up for processing.\n")
        layout = bids.layout.BIDSLayout(data_dir_path, validate=False)

    split_name, scan_info, run_iter, structural_scan_list, bold_scan_list = prep_bids_iter(
        layout, opts.bids_filter, opts.bold_only, inclusion_list=opts.inclusion_ids, exclusion_list=opts.exclusion_ids)
    '''***details on outputs from prep_bids_iter:
    split_name: a list of strings, providing a sensible name to distinguish each iterable, 
        and also necessary to link up the run iterables with a specific session later.
    scan_info: a list of dictionary including the subject ID and session # for a given 
        iterable from split_name
    run_iter: a list of dictionary, where the keys correspond to a session split from 
        split_name, and the value is a list of runs for that split. This manages iterables
        for runs.
    structural_scan_list: list of structural file names used for commonspace registration; used for resample_template and managing # threads
    bold_scan_list: list of functional file names; used for resample_template and managing # threads
    '''
    number_structural_scans = len(structural_scan_list)
    number_functional_scans = len(bold_scan_list)

    # setting up all iterables
    main_split = pe.Node(niu.IdentityInterface(fields=['split_name', 'scan_info']),
                         name="main_split")
    main_split.iterables = [('split_name', split_name),
                            ('scan_info', scan_info)]
    main_split.synchronize = True

    '''
    First execute a workflow that finds and formats inputs. 
    This workflow includes a first split and merge of iterables, which is necessary to conduct separately to avoid
    downstream errors with iterables.
    '''
    prep_input_wf = init_prep_input_wf(data_dir_path, split_name, scan_info, run_iter, opts, name='prep_input_wf')

    '''
    MATCHING INPUTS TO MAIN ITERABLES
    '''
    format_bold_node = pe.Node(Function(input_names=['scan_info', 'run', 'joined_scan_info', 'joined_run', 'joined_file_list'],
                                                output_names=['selected_file'],
                                                function=match_iterables),
                                        name='format_bold')
    
    input_bold_node = pe.Node(Function(input_names=['scan_info', 'run', 'joined_scan_info', 'joined_run', 'joined_file_list'],
                                                output_names=['selected_file'],
                                                function=match_iterables),
                                        name='input_bold')

    workflow.connect([
        (prep_input_wf, format_bold_node, [
            ("outputnode.prep_bold_list", "joined_file_list"),
            ("outputnode.joined_scan_info", "joined_scan_info"),
            ]),
        (main_split, format_bold_node, [
            ("scan_info", "scan_info"),
            ]),
        (prep_input_wf, input_bold_node, [
            ("outputnode.input_bold_list", "joined_file_list"),
            ("outputnode.joined_scan_info", "joined_scan_info"),
            ]),
        (main_split, input_bold_node, [
            ("scan_info", "scan_info"),
            ]),
        ])
    
    if opts.bold_only:
        format_bold_node.inputs.run = None
        input_bold_node.inputs.run = None
        format_bold_node.inputs.joined_run = None
        input_bold_node.inputs.joined_run = None
    else:
        run_split = pe.Node(niu.IdentityInterface(fields=['run', 'split_name']),
                            name="run_split")
        run_split.itersource = ('main_split', 'split_name')
        run_split.iterables = [('run', run_iter)]

        workflow.connect([
            (main_split, run_split, [
                ("split_name", "split_name"),
                ]),
            (prep_input_wf, format_bold_node, [
                ("outputnode.joined_run", "joined_run"),
                ]),
            (run_split, format_bold_node, [
                ("run", "run"),
                ]),
            (prep_input_wf, input_bold_node, [
                ("outputnode.joined_run", "joined_run"),
                ]),
            (run_split, input_bold_node, [
                ("run", "run"),
                ]),
            ])

        format_anat_node = pe.Node(Function(input_names=['scan_info', 'run', 'joined_scan_info', 'joined_run', 'joined_file_list'],
                                                    output_names=['selected_file'],
                                                    function=match_iterables),
                                            name='format_anat')
        format_anat_node.inputs.run = None
        format_anat_node.inputs.joined_run = None
        
        input_anat_node = pe.Node(Function(input_names=['scan_info', 'run', 'joined_scan_info', 'joined_run', 'joined_file_list'],
                                                    output_names=['selected_file'],
                                                    function=match_iterables),
                                            name='input_anat')
        input_anat_node.inputs.run = None
        input_anat_node.inputs.joined_run = None
        
        workflow.connect([
            (prep_input_wf, format_anat_node, [
                ("outputnode.prep_anat_list", "joined_file_list"),
                ("outputnode.anat_joined_scan_info", "joined_scan_info"),
                ]),
            (main_split, format_anat_node, [
                ("scan_info", "scan_info"),
                ]),
            (prep_input_wf, input_anat_node, [
                ("outputnode.input_anat_list", "joined_file_list"),
                ("outputnode.anat_joined_scan_info", "joined_scan_info"),
                ]),
            (main_split, input_anat_node, [
                ("scan_info", "scan_info"),
                ]),
            ])


    if opts.inherit_unbiased_template=='none':
        inherit_unbiased=False
        # Resample the anatomical template according to the resolution of the provided input data
        resample_template_node = pe.Node(Function(input_names=['opts', 'structural_scan_list', 'bold_scan_list'],
                                                output_names=[
                                                    'registration_template', 'registration_mask', 'commonspace_template'],
                                                function=resample_template),
                                        name='resample_template', mem_gb=1*opts.scale_min_memory)
        resample_template_node.inputs.opts = opts

        workflow.connect([
            (prep_input_wf, resample_template_node, [
                ("outputnode.prep_anat_list", "structural_scan_list"),
                ("outputnode.prep_bold_list", "bold_scan_list"),
                ]),
            ])

    else: # inherit the atlas files from previous run
        inherit_unbiased=True
        opts, inherit_dict = inherit_unbiased_files(opts.inherit_unbiased_template, opts)
        resample_template_node = pe.Node(niu.IdentityInterface(fields=['registration_template', 'registration_mask', 'commonspace_template']),
                                            name="resample_template")
        resample_template_node.inputs.registration_template = inherit_dict['registration_template']
        resample_template_node.inputs.registration_mask = inherit_dict['registration_mask']
        resample_template_node.inputs.commonspace_template = inherit_dict['commonspace_template']

    # calculate the number of scans that will be registered
    num_procs = min(opts.local_threads, number_structural_scans)

    EPI_target_buffer = pe.Node(niu.IdentityInterface(fields=['EPI_template', 'EPI_mask']),
                                        name="EPI_target_buffer")

    commonspace_reg_wf = init_commonspace_reg_wf(opts=opts, commonspace_reg_opts=opts.commonspace_reg, inherit_unbiased=inherit_unbiased,
                                                 output_folder=output_folder, transforms_datasink=transforms_datasink, num_procs=num_procs, output_datasinks=True, 
                                                 joinsource_list=['main_split'], name='commonspace_reg_wf')
    if inherit_unbiased:
        commonspace_reg_wf.inputs.inherit_unbiased_inputnode.unbiased_template = inherit_dict['unbiased_template']
        commonspace_reg_wf.inputs.inherit_unbiased_inputnode.unbiased_mask = inherit_dict['unbiased_mask']
        commonspace_reg_wf.inputs.inherit_unbiased_inputnode.unbiased_to_atlas_affine = inherit_dict['unbiased_to_atlas_affine']
        commonspace_reg_wf.inputs.inherit_unbiased_inputnode.unbiased_to_atlas_warp = inherit_dict['unbiased_to_atlas_warp']
        commonspace_reg_wf.inputs.inherit_unbiased_inputnode.unbiased_to_atlas_inverse_warp = inherit_dict['unbiased_to_atlas_inverse_warp']
        commonspace_reg_wf.inputs.inherit_unbiased_inputnode.warped_unbiased = inherit_dict['warped_unbiased']

    bold_main_wf = init_bold_main_wf(opts=opts, output_folder=output_folder, number_functional_scans=number_functional_scans)

    # organizing visual QC outputs
    template_diagnosis = pe.Node(Function(input_names=['anat_template', 'opts', 'out_dir', 'figure_format'],
                                       function=preprocess_visual_QC.template_info),
                              name='template_info')
    template_diagnosis.inputs.opts = opts
    template_diagnosis.inputs.out_dir = output_folder+'/preprocess_QC_report/template_files/'
    template_diagnosis.inputs.figure_format = opts.figure_format

    bold_inho_cor_diagnosis = pe.Node(Function(input_names=['raw_img','init_denoise','warped_mask','final_denoise', 'name_source', 'out_dir', 'figure_format'],
                                       function=preprocess_visual_QC.inho_cor_diagnosis),
                              name='bold_inho_cor_diagnosis')
    bold_inho_cor_diagnosis.inputs.out_dir = output_folder+'/preprocess_QC_report/bold_inho_cor/'
    bold_inho_cor_diagnosis.inputs.figure_format = opts.figure_format

    temporal_diagnosis = pe.Node(Function(input_names=['bold_file', 'motion_params_csv', 'FD_csv', 'rabies_data_type', 'name_source', 'out_dir', 'figure_format'],
                                          output_names=[
                                            'std_filename', 'tSNR_filename'],
                                       function=preprocess_visual_QC.temporal_features),
                              name='temporal_features')
    temporal_diagnosis.inputs.out_dir = output_folder+'/preprocess_QC_report/temporal_features/'
    temporal_diagnosis.inputs.rabies_data_type = opts.data_type
    temporal_diagnosis.inputs.figure_format = opts.figure_format

    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (format_bold_node, bold_main_wf, [
            ("selected_file", "inputnode.bold"),
            ]),
        (resample_template_node, template_diagnosis, [
            ("registration_template", "anat_template"),
            ]),
        (resample_template_node, commonspace_reg_wf, [
            ("registration_template", "template_inputnode.template_anat"),
            ("registration_mask", "template_inputnode.template_mask"),
            ]),
        (resample_template_node, bold_main_wf, [
            ("commonspace_template", "inputnode.commonspace_ref"),
            ]),
        (resample_template_node, outputnode, [
            ("commonspace_template", "commonspace_resampled_template"),
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
            ("outputnode.native_vascular_mask", "native_vascular_mask"),
            ("outputnode.native_labels", "native_labels"),
            ("outputnode.motion_params_csv", "motion_params_csv"),
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
            ("outputnode.raw_brain_mask", "raw_brain_mask"),
            ]),
        (bold_main_wf, bold_inho_cor_diagnosis, [
            ("outputnode.bold_ref", "raw_img"),
            ("outputnode.init_denoise", "init_denoise"),
            ("outputnode.corrected_EPI", "final_denoise"),
            ("outputnode.denoise_mask", "warped_mask"),
            ]),
        (input_bold_node, bold_inho_cor_diagnosis,
         [("selected_file", "name_source")]),
        (bold_main_wf, temporal_diagnosis, [
            ("outputnode.commonspace_bold", "bold_file"),
            ("outputnode.motion_params_csv", "motion_params_csv"),
            ("outputnode.FD_csv", "FD_csv"),
            ]),
        (input_bold_node, temporal_diagnosis,
         [("selected_file", "name_source")]),
        (temporal_diagnosis, outputnode, [
            ("tSNR_filename", "tSNR_filename"),
            ("std_filename", "std_filename"),
            ]),
        ])

    if not opts.bold_only:
        # setting anat preprocessing nodes
        anat_inho_cor_wf = init_inho_correction_wf(opts=opts, image_type='structural', output_folder=output_folder, num_procs=num_procs, name="anat_inho_cor_wf")

        workflow.connect([
            (format_anat_node, anat_inho_cor_wf, [
                ("selected_file", "inputnode.target_img"),
                ("selected_file", "inputnode.name_source"),
                ]),
            (resample_template_node, anat_inho_cor_wf, [
                ("registration_template", "inputnode.anat_ref"),
                ("registration_mask", "inputnode.anat_mask"),
                ]),
            (resample_template_node, anat_inho_cor_wf, [
                ("registration_template", "template_inputnode.template_anat"),
                ("registration_mask", "template_inputnode.template_mask"),
                ]),
            (anat_inho_cor_wf, bold_main_wf, [
                ("outputnode.corrected", "inputnode.coreg_anat"),
                ]),
            (commonspace_reg_wf, bold_main_wf, [
                ("outputnode.native_mask", "inputnode.coreg_mask"),
                ("outputnode.unbiased_template", "template_inputnode.template_anat"),
                ("outputnode.unbiased_mask", "template_inputnode.template_mask"),
                ]),
            (anat_inho_cor_wf, commonspace_reg_wf, [
                ("outputnode.corrected", "inputnode.moving_image"),
                ("outputnode.denoise_mask", "inputnode.moving_mask"),
                ]),
            (anat_inho_cor_wf, EPI_target_buffer, [
                ("outputnode.corrected", "EPI_template"),
                ]),
            (commonspace_reg_wf, EPI_target_buffer, [
                ("outputnode.native_mask", 'EPI_mask'),
                ]),
            (EPI_target_buffer, bold_main_wf, [
                ("EPI_template", "inputnode.inho_cor_anat"),
                ("EPI_mask", "inputnode.inho_cor_mask"),
                ]),
            ])

        if not opts.anat_inho_cor['method']=='disable':
            anat_inho_cor_diagnosis = pe.Node(Function(input_names=['raw_img','init_denoise','warped_mask','final_denoise', 'name_source', 'out_dir', 'figure_format'],
                                               function=preprocess_visual_QC.inho_cor_diagnosis),
                                      name='anat_inho_cor_diagnosis')
            anat_inho_cor_diagnosis.inputs.out_dir = output_folder+'/preprocess_QC_report/anat_inho_cor/'
            anat_inho_cor_diagnosis.inputs.figure_format = opts.figure_format

            workflow.connect([
                (format_anat_node, anat_inho_cor_diagnosis, [
                    ("selected_file", "raw_img"),
                    ]),
                (input_anat_node, anat_inho_cor_diagnosis, [
                    ("selected_file", "name_source"),
                    ]),
                (anat_inho_cor_wf, anat_inho_cor_diagnosis, [
                    ("outputnode.init_denoise", "init_denoise"),
                    ("outputnode.corrected", "final_denoise"),
                    ("outputnode.denoise_mask", "warped_mask"),
                    ]),
                ])

    else:
        inho_cor_bold_main_wf = init_bold_main_wf(
            output_folder=output_folder, number_functional_scans=number_functional_scans, inho_cor_only=True, name='inho_cor_bold_main_wf', opts=opts)

        workflow.connect([
            (resample_template_node, inho_cor_bold_main_wf, [
                ("registration_template", "template_inputnode.template_anat"),
                ("registration_mask", "template_inputnode.template_mask"),
                ]),
            (format_bold_node, inho_cor_bold_main_wf, [
                ("selected_file", "inputnode.bold"),
                ]),
            (EPI_target_buffer, inho_cor_bold_main_wf, [
                ("EPI_template", "inputnode.inho_cor_anat"),
                ("EPI_mask", "inputnode.inho_cor_mask"),
                ]),
            (inho_cor_bold_main_wf, bold_main_wf, [
                ("transitionnode.bold_file", "transitionnode.bold_file"),
                ("transitionnode.isotropic_bold_file", "transitionnode.isotropic_bold_file"),
                ("transitionnode.bold_ref", "transitionnode.bold_ref"),
                ("transitionnode.init_denoise", "transitionnode.init_denoise"),
                ("transitionnode.denoise_mask", "transitionnode.denoise_mask"),
                ("transitionnode.corrected_EPI", "transitionnode.corrected_EPI"),
                ("transitionnode.log_bold", "transitionnode.log_bold"),
                ]),
            (inho_cor_bold_main_wf, commonspace_reg_wf, [
                ("transitionnode.corrected_EPI", "inputnode.moving_image"),
                ("transitionnode.denoise_mask", "inputnode.moving_mask"),
                ]),
            (resample_template_node, EPI_target_buffer, [
                ("registration_template", "EPI_template"),
                ("registration_mask", "EPI_mask"),
                ]),
            ])

    if not opts.bold_only:
        PlotOverlap_EPI2Anat_node = pe.Node(
            preprocess_visual_QC.PlotOverlap(), name='PlotOverlap_EPI2Anat')
        PlotOverlap_EPI2Anat_node.inputs.out_dir = output_folder+'/preprocess_QC_report/EPI2Anat'
        workflow.connect([
            (input_bold_node, PlotOverlap_EPI2Anat_node,
             [("selected_file", "name_source")]),
            (anat_inho_cor_wf, PlotOverlap_EPI2Anat_node,
             [("outputnode.corrected", "fixed")]),
            (outputnode, PlotOverlap_EPI2Anat_node, [
                ("inho_cor_bold_warped2anat", "moving"),  # warped EPI to anat
                ]),
            ])

    # fill the datasinks
    workflow.connect([
        (input_bold_node, bold_datasink, [
            ("selected_file", "input_bold"),
            ]),
        (outputnode, motion_datasink, [
            ("motion_params_csv", "motion_params_csv"),  # confounds file
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
            ("native_vascular_mask", "native_vascular_mask"),  # get the EPI labels
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
            ("raw_brain_mask", "raw_brain_mask"),
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
