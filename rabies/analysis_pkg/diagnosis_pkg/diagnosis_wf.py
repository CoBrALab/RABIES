import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function

from rabies.analysis_pkg.diagnosis_pkg.interfaces import ScanDiagnosis, PrepMasks, DatasetDiagnosis
from rabies.analysis_pkg.diagnosis_pkg.diagnosis_functions import temporal_external_formating, spatial_external_formating


def init_diagnosis_wf(analysis_opts, commonspace_bold, preprocess_opts, scan_split_name, name="diagnosis_wf"):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['mask_dict_list', 'file_dict', 'analysis_dict']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['figure_temporal_diagnosis', 'figure_spatial_diagnosis', 
                                                       'dataset_diagnosis', 'temporal_info_csv', 'spatial_VE_nii', 'temporal_std_nii', 'GS_corr_nii',
                                                       'temporal_std_nii', 'GS_corr_nii', 'DVARS_corr_nii', 'FD_corr_nii']), name='outputnode')

    if not (commonspace_bold or preprocess_opts.bold_only):
        raise ValueError("Cannot currently select --nativespace_analysis for running data_diagnosis")

    if os.path.basename(preprocess_opts.anat_template)=='DSURQE_40micron_average.nii.gz':
        DSURQE_regions=True
    else:
        DSURQE_regions=False

    ScanDiagnosis_node = pe.Node(ScanDiagnosis(prior_bold_idx=analysis_opts.prior_bold_idx,
        prior_confound_idx=analysis_opts.prior_confound_idx,
            dual_ICA = analysis_opts.dual_ICA, DSURQE_regions=DSURQE_regions),
        name='ScanDiagnosis')

    PrepMasks_node = pe.Node(PrepMasks(prior_maps=os.path.abspath(str(analysis_opts.prior_maps)), DSURQE_regions=DSURQE_regions),
        name='PrepMasks')

    temporal_external_formating_node = pe.Node(Function(input_names=['temporal_info', 'file_dict'],
                                            output_names=[
                                                'temporal_info_csv', 'dual_regression_timecourse_csv', 'dual_ICA_timecourse_csv'],
                                        function=temporal_external_formating),
                                name='temporal_external_formating')

    spatial_external_formating_node = pe.Node(Function(input_names=['spatial_info', 'file_dict'],
                                            output_names=[
                                                'VE_filename', 'std_filename', 'predicted_std_filename', 'GS_corr_filename', 'DVARS_corr_filename', 'FD_corr_filename', 'DR_maps_filename', 'dual_ICA_filename'],
                                        function=spatial_external_formating),
                                name='spatial_external_formating')


    workflow.connect([
        (inputnode, PrepMasks_node, [
            ("mask_dict_list", "mask_dict_list"),
            ]),
        (PrepMasks_node, ScanDiagnosis_node, [
            ("mask_file_dict", "mask_file_dict"),
            ]),
        (inputnode, ScanDiagnosis_node, [
            ("file_dict", "file_dict"),
            ("analysis_dict", "analysis_dict"),
            ]),
        (inputnode, temporal_external_formating_node, [
            ("file_dict", "file_dict"),
            ]),
        (ScanDiagnosis_node, temporal_external_formating_node, [
            ("temporal_info", "temporal_info"),
            ]),
        (inputnode, spatial_external_formating_node, [
            ("file_dict", "file_dict"),
            ]),
        (ScanDiagnosis_node, spatial_external_formating_node, [
            ("spatial_info", "spatial_info"),
            ]),
        ])

    workflow.connect([
        (ScanDiagnosis_node, outputnode, [
            ("figure_temporal_diagnosis", "figure_temporal_diagnosis"),
            ("figure_spatial_diagnosis", "figure_spatial_diagnosis"),
            ]),
        (temporal_external_formating_node, outputnode, [
            ("temporal_info_csv", "temporal_info_csv"),
            ("dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
            ]),
        (spatial_external_formating_node, outputnode, [
            ("VE_filename", "spatial_VE_nii"),
            ("std_filename", "temporal_std_nii"),
            ("predicted_std_filename", "CR_prediction_std_nii"),
            ("GS_corr_filename", "GS_corr_nii"),
            ("DVARS_corr_filename", "DVARS_corr_nii"),
            ("FD_corr_filename", "FD_corr_nii"),
            ("DR_maps_filename", "dual_regression_nii"),
            ]),
        ])

    if not len(scan_split_name)<3:

        def prep_scan_data(spatial_info, analysis_dict, file_dict, mask_file_dict):
            scan_data={}

            dict_keys = ['temporal_std', 'predicted_std', 'VE_spatial', 'GS_corr',
                            'DVARS_corr', 'FD_corr', 'DR_BOLD', 'dual_ICA_maps', 'prior_maps' ]
            for key in dict_keys:
                scan_data[key] = spatial_info[key]

            scan_data['tDOF'] = file_dict['CR_data_dict']['tDOF']

            import numpy as np
            import nibabel as nb
            mask_file = mask_file_dict['brain_mask']
            brain_mask = np.asarray(nb.load(mask_file).dataobj)
            volume_indices = brain_mask.astype(bool)

            seed_list=[]
            for seed_map in analysis_dict['seed_map_files']:
                seed_list.append(np.asarray(
                    nb.load(seed_map).dataobj)[volume_indices])
            scan_data['seed_list'] = seed_list
            return scan_data

        prep_scan_data_node = pe.Node(Function(input_names=['spatial_info', 'analysis_dict', 'file_dict', 'mask_file_dict'],
                                            output_names=['scan_data'],
                                        function=prep_scan_data),
                                name='prep_scan_data_node')

        # calculate the number of scans combined in diagnosis
        num_scan = len(scan_split_name)
        num_procs = min(analysis_opts.local_threads, num_scan)

        data_diagnosis_split_joinnode = pe.JoinNode(niu.IdentityInterface(fields=['scan_data_list']),
                                                name='diagnosis_split_joinnode',
                                                joinsource='main_split',
                                                joinfield=['scan_data_list'],
                                                n_procs=num_procs, mem_gb=1*num_scan*analysis_opts.scale_min_memory)

        DatasetDiagnosis_node = pe.Node(DatasetDiagnosis(),
            name='DatasetDiagnosis',
            n_procs=num_procs, mem_gb=1*num_scan*analysis_opts.scale_min_memory)
        DatasetDiagnosis_node.inputs.seed_prior_maps = analysis_opts.seed_prior_list

        workflow.connect([
            (inputnode, prep_scan_data_node, [
                ("file_dict", "file_dict"),
                ("analysis_dict", "analysis_dict"),
                ]),
            (ScanDiagnosis_node, prep_scan_data_node, [
                ("spatial_info", "spatial_info"),
                ]),
            (PrepMasks_node, prep_scan_data_node, [
                ("mask_file_dict", "mask_file_dict"),
                ]),
            (prep_scan_data_node, data_diagnosis_split_joinnode, [
                ("scan_data", "scan_data_list"),
                ]),
            (data_diagnosis_split_joinnode, DatasetDiagnosis_node, [
                ("scan_data_list", "scan_data_list"),
                ]),
            (PrepMasks_node, DatasetDiagnosis_node, [
                ("mask_file_dict", "mask_file_dict"),
                ]),
            (DatasetDiagnosis_node, outputnode, [
                ("dataset_diagnosis", "dataset_diagnosis"),
                ]),
            ])
    else:
        from nipype import logging
        log = logging.getLogger('nipype.workflow')
        log.warning(
            "Cannot run statistics on a sample size smaller than 3, so dataset diagnosis is not run.")


    return workflow
