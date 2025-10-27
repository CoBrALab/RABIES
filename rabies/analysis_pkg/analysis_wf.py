import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function

from rabies.utils import ResampleVolumes
from .analysis_functions import run_group_ICA, run_DR_ICA, run_FC_matrix, seed_based_FC


def init_analysis_wf(opts, nativespace_analysis=False, name="analysis_wf"):

    workflow = pe.Workflow(name=name)
    subject_inputnode = pe.Node(niu.IdentityInterface(
        fields=['CR_dict_file', 'maps_dict_file', 'token', 'native_to_commonspace_transform_list', 'native_to_commonspace_inverse_list']), name='subject_inputnode')
    group_inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file_list', 'commonspace_mask', 'commonspace_template', 'token']), name='group_inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['group_ICA_dir', 'IC_file', 'dual_regression_timecourse_csv',
                                                       'DR_nii_file', 'DR_nii_file_resampled', 'matrix_data_file', 'matrix_fig', 'corr_map_file', 'corr_map_file_resampled', 'seed_timecourse_csv', 'joined_corr_map_file', 'joined_seed_timecourse_csv',
                                                       'sub_token', 'group_token','NPR_prior_timecourse_csv', 'NPR_extra_timecourse_csv',
                                                       'NPR_prior_filename', 'NPR_extra_filename', 'NPR_optimize_report']), name='outputnode')

    # connect the nodes so that they exist even without running analysis
    workflow.connect([
        (subject_inputnode, outputnode, [
            ("token", "sub_token"),
            ]),
        (group_inputnode, outputnode, [
            ("token", "group_token"),
            ]),
        ])

    if len(opts.seed_list) > 0:
        seed_based_FC_node = pe.Node(Function(input_names=['CR_dict_file', 'maps_dict_file', 'seed_name'],
                                              output_names=['corr_map_file', 'seed_timecourse_csv'],
                                              function=seed_based_FC),
                                     name='seed_based_FC', mem_gb=1*opts.scale_min_memory)
        seed_based_FC_node.iterables = ('seed_name', opts.seed_name_list)

        # create a joinnode to provide also combined seed maps
        seed_FC_joinnode = pe.JoinNode(niu.IdentityInterface(fields=['joined_corr_map_file', 'joined_seed_timecourse_csv']),
                                                name='seed_FC_joinnode',
                                                joinsource='seed_based_FC',
                                                joinfield=['joined_corr_map_file', 'joined_seed_timecourse_csv'])


        workflow.connect([
            (subject_inputnode, seed_based_FC_node, [
                ("CR_dict_file", "CR_dict_file"),
                ("maps_dict_file", "maps_dict_file"),
                ]),
            (seed_based_FC_node, outputnode, [
                ("corr_map_file", "corr_map_file"),
                ("seed_timecourse_csv", "seed_timecourse_csv"),
                ]),
            (seed_based_FC_node, seed_FC_joinnode, [
                ("seed_timecourse_csv", "joined_seed_timecourse_csv"),
                ]),
            (seed_FC_joinnode, outputnode, [
                ("joined_corr_map_file", "joined_corr_map_file"),
                ("joined_seed_timecourse_csv", "joined_seed_timecourse_csv"),
                ]),
            ])


        if opts.resample_to_commonspace:
            SBC_transform_node = pe.Node(ResampleVolumes(
                resampling_dim='ref_file', interpolation=opts.interpolation,
                rabies_data_type=opts.data_type, apply_motcorr=False, clip_negative=False), 
                name='SBC_to_commonspace')
            
            workflow.connect([
                (subject_inputnode, SBC_transform_node, [
                    ("native_to_commonspace_transform_list", "transforms"),
                    ("native_to_commonspace_inverse_list", "inverses"),
                    ]),
                (group_inputnode, SBC_transform_node, [
                    ("commonspace_template", "ref_file"),
                    ]),
                (seed_based_FC_node, SBC_transform_node, [
                    ("corr_map_file", "in_file"),
                    ("corr_map_file", "name_source"),
                    ]),
                (SBC_transform_node, outputnode, [
                    ("resampled_file", "corr_map_file_resampled"),
                    ]),
                (SBC_transform_node, seed_FC_joinnode, [
                    ("resampled_file", "joined_corr_map_file"),
                    ]),
                ])
        else:
            workflow.connect([
                (seed_based_FC_node, seed_FC_joinnode, [
                    ("corr_map_file", "joined_corr_map_file"),
                    ]),
                ])



    if opts.DR_ICA or opts.data_diagnosis:

        DR_ICA = pe.Node(Function(input_names=['CR_dict_file', 'maps_dict_file', 'network_weighting'],
                                  output_names=['DR_maps_filename', 'dual_regression_timecourse_csv'],
                                  function=run_DR_ICA),
                         name='DR_ICA', mem_gb=1*opts.scale_min_memory)
        DR_ICA.inputs.network_weighting = opts.network_weighting

        workflow.connect([
            (subject_inputnode, DR_ICA, [
                ("CR_dict_file", "CR_dict_file"),
                ("maps_dict_file", "maps_dict_file"),
                ]),
            (DR_ICA, outputnode, [
                ("dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
                ("DR_maps_filename", "DR_nii_file"),
                ]),
            ])

        if opts.resample_to_commonspace:
            DR_transform_node = pe.Node(ResampleVolumes(
                resampling_dim='ref_file', interpolation=opts.interpolation,
                rabies_data_type=opts.data_type, apply_motcorr=False, clip_negative=False), 
                name='DR_to_commonspace')
            
            workflow.connect([
                (subject_inputnode, DR_transform_node, [
                    ("native_to_commonspace_transform_list", "transforms"),
                    ("native_to_commonspace_inverse_list", "inverses"),
                    ]),
                (group_inputnode, DR_transform_node, [
                    ("commonspace_template", "ref_file"),
                    ]),
                (DR_ICA, DR_transform_node, [
                    ("DR_maps_filename", "in_file"),
                    ("DR_maps_filename", "name_source"),
                    ]),
                (DR_transform_node, outputnode, [
                    ("resampled_file", "DR_nii_file_resampled"),
                    ]),
                ])

    if opts.group_ica['apply']:
        group_ICA = pe.Node(Function(input_names=['bold_file_list', 'mask_file', 'dim', 'random_seed', 'background_image', 'disableMigp'],
                                     output_names=['out_dir', 'IC_file'],
                                     function=run_group_ICA),
                            name='group_ICA', mem_gb=1*opts.scale_min_memory)
        group_ICA.inputs.dim = opts.group_ica['dim']
        group_ICA.inputs.random_seed = opts.group_ica['random_seed']
        group_ICA.inputs.disableMigp = opts.group_ica['disableMigp']

        workflow.connect([
            (group_inputnode, group_ICA, [
                ("bold_file_list", "bold_file_list"),
                ("commonspace_mask", "mask_file"),
                ("commonspace_template", "background_image"),
                ]),
            (group_ICA, outputnode, [
                ("IC_file", "IC_file"),
                ("out_dir", "group_ICA_dir"),
                ]),
            ])

    if (opts.NPR_temporal_comp>-1) or (opts.NPR_spatial_comp>-1) or opts.optimize_NPR['apply']:
        from .analysis_functions import NeuralPriorRecovery

        NPR_node = pe.Node(NeuralPriorRecovery(
                            NPR_temporal_comp=opts.NPR_temporal_comp, 
                            NPR_spatial_comp=opts.NPR_spatial_comp, 
                            optimize_NPR_dict=opts.optimize_NPR, 
                            prior_bold_idx = opts.prior_bold_idx,
                            network_weighting=opts.network_weighting,
                            figure_format=opts.figure_format),
                            name='NPR_node', mem_gb=1*opts.scale_min_memory)

        workflow.connect([
            (subject_inputnode, NPR_node, [
                ("CR_dict_file", "CR_dict_file"),
                ("maps_dict_file", "maps_dict_file"),
                ]),
            (NPR_node, outputnode, [
                ("NPR_prior_timecourse_csv", "NPR_prior_timecourse_csv"),
                ("NPR_extra_timecourse_csv", "NPR_extra_timecourse_csv"),
                ("NPR_prior_filename", "NPR_prior_filename"),
                ("NPR_extra_filename", "NPR_extra_filename"),
                ]),
            ])
        
    if opts.optimize_NPR['apply']:
        workflow.connect([
            (NPR_node, outputnode, [
                ("optimize_report", "NPR_optimize_report"),
                ]),
            ])


    if opts.FC_matrix:
        FC_matrix = pe.Node(Function(input_names=['CR_dict_file', 'maps_dict_file', 'figure_format', 'roi_type'],
                                     output_names=['data_file', 'figname'],
                                     function=run_FC_matrix),
                            name='FC_matrix', mem_gb=1*opts.scale_min_memory)
        FC_matrix.inputs.roi_type = opts.ROI_type
        FC_matrix.inputs.figure_format = opts.figure_format

        workflow.connect([
            (subject_inputnode, FC_matrix, [
                ("CR_dict_file", "CR_dict_file"),
                ("maps_dict_file", "maps_dict_file"),
                ]),
            (FC_matrix, outputnode, [
                ("data_file", "matrix_data_file"),
                ("figname", "matrix_fig"),
                ]),
            ])

    return workflow
