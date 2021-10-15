
def integrate_data_diagnosis(workflow, outputnode, confound_correction_wf, data_diagnosis_opts, bold_only, commonspace_bold, bold_scan_list):

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
        (confound_correction_wf, prep_dict_node, [
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
