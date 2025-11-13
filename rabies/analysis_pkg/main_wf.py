import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
from rabies.analysis_pkg.diagnosis_pkg.diagnosis_wf import init_diagnosis_wf
from rabies.analysis_pkg.analysis_wf import init_analysis_wf
from rabies.utils import fill_split_dict, get_workflow_dict


def init_main_analysis_wf(preprocess_opts, cr_opts, analysis_opts):

    if preprocess_opts.labels is None and analysis_opts.FC_matrix and analysis_opts.ROI_type=='parcellated':
        raise ValueError(
            'No labels were inherited from preprocessing to derive parcellated connectivity matrices.')

    workflow = pe.Workflow(name='analysis_main_wf')

    conf_output = os.path.abspath(str(analysis_opts.confound_correction_out))

    split_dict, split_name_list, target_list = read_confound_workflow(conf_output, cr_opts)

    if len(split_name_list)==0:
        raise ValueError(f"""
            No outputs were founds from the confound correction stage. 
            All scans may have been removed for not meeting the minimum_timepoint threshold 
            when applying --frame_censoring. Outputs will be named empty.nii.gz if this is 
            the case.
            """)

    if analysis_opts.group_ica['apply'] and cr_opts.nativespace_analysis and not cr_opts.resample_to_commonspace:
        raise ValueError(
            'No commonspace timeseries are found for running group-ICA. Try applying --resample_to_commonspace to generate these timeseries when using --nativespace_analysis.')

    if analysis_opts.data_diagnosis and cr_opts.nativespace_analysis:
        analysis_opts.resample_to_commonspace = True # resampling to common is necessary for --data_diagnosis report

    # filter inclusion/exclusion lists
    from rabies.utils import filter_scan_inclusion, filter_scan_exclusion
    split_name_list = filter_scan_inclusion(analysis_opts.inclusion_ids, split_name_list)
    split_name_list = filter_scan_exclusion(analysis_opts.exclusion_ids, split_name_list)

    # setting up iterables from the BOLD scan splits
    main_split = pe.Node(niu.IdentityInterface(fields=['split_name']),
                         name="main_split")
    main_split.iterables = [('split_name', split_name_list)]

    # function to read each elements from the dictionary of confound correction outputs
    def read_dict(split_dict, split_name, target_list):
        return [split_dict[split_name][target] for target in target_list]

    # set outputnode from confound correction
    conf_outputnode = pe.Node(Function(input_names=['split_dict', 'split_name', 'target_list'],
                                           output_names=target_list,
                                       function=read_dict),
                              name='conf_outputnode')
    conf_outputnode.inputs.split_dict = split_dict
    conf_outputnode.inputs.target_list = target_list

    # prepare dictionary for seed files
    seed_dict = {}
    seed_name_list = []
    for file in analysis_opts.seed_list:
        file = os.path.abspath(file)
        if not os.path.isfile(file):
            raise ValueError(
                f"Provide seed file path {file} doesn't exists.")
        seed_name = pathlib.Path(file).name.rsplit(".nii")[0]
        seed_name_list.append(seed_name)
        seed_dict[seed_name] = file
    analysis_opts.seed_name_list = seed_name_list # store for developing iterables later

    if not os.path.isfile(str(analysis_opts.prior_maps)):
        raise ValueError("--prior_maps doesn't exists.")

    load_CR_dict_node = pe.Node(Function(input_names=['maps_dict_file', 'cleaned_bold_file', 'CR_data_dict', 'VE_file', 'STD_file', 'CR_STD_file', 'name_source'],
                                           output_names=[
                                               'CR_dict_file'],
                                       function=load_CR_input_dict),
                              name='load_CR_dict_node')

    # prepare analysis workflow
    analysis_output = os.path.abspath(str(analysis_opts.output_dir))
    analysis_wf = init_analysis_wf(
        opts=analysis_opts)

    '''
    PREPARE DATASINKS
    '''
    if cr_opts.nativespace_analysis:
        # prepare analysis datasink for nativespace computations
        nativespace_analysis_datasink = pe.Node(DataSink(base_directory=analysis_output,
                                            container="nativespace_analysis_datasink"),
                                    name="nativespace_analysis_datasink")
        main_analysis_datasink = nativespace_analysis_datasink

    if analysis_opts.resample_to_commonspace or not cr_opts.nativespace_analysis or analysis_opts.group_ica['apply']:
        # prepare analysis datasink for commonspace computations OR nativespace computations resampled to commonspace
        commonspace_analysis_datasink = pe.Node(DataSink(base_directory=analysis_output,
                                            container="commonspace_analysis_datasink"),
                                    name="commonspace_analysis_datasink")
        if not cr_opts.nativespace_analysis:
            main_analysis_datasink = commonspace_analysis_datasink
        if analysis_opts.group_ica['apply']:
            workflow.connect([
                (analysis_wf, commonspace_analysis_datasink, [
                    ("outputnode.group_ICA_dir", "group_ICA_dir"),
                    ]),
                ])

    data_diagnosis_datasink = pe.Node(DataSink(base_directory=analysis_output,
                                         container="data_diagnosis_datasink"),
                                name="data_diagnosis_datasink")
    
    '''
    CORE WORKFLOW CONNECTIONS
    '''

    workflow.connect([
        (main_split, conf_outputnode, [
            ("split_name", "split_name"),
            ]),
        (conf_outputnode, load_CR_dict_node, [
            ("cleaned_path", "cleaned_bold_file"),
            ("data_dict", "CR_data_dict"),
            ("VE_file_path", "VE_file"),
            ("STD_file_path", "STD_file"),
            ("CR_STD_file_path", "CR_STD_file"),
            ("input_bold", "name_source"),
            ]),
        (load_CR_dict_node, analysis_wf, [
            ("CR_dict_file", "subject_inputnode.CR_dict_file"),
            ]),
        (analysis_wf, main_analysis_datasink, [
            ("outputnode.matrix_data_file", "matrix_data_file"),
            ("outputnode.matrix_fig", "matrix_fig"),
            ("outputnode.corr_map_file", "seed_correlation_maps"),
            ("outputnode.seed_timecourse_csv", "seed_timecourse_csv"),
            ("outputnode.DR_nii_file", "dual_regression_nii"),
            ("outputnode.dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
            ("outputnode.NPR_prior_timecourse_csv", "NPR_prior_timecourse_csv"),
            ("outputnode.NPR_extra_timecourse_csv", "NPR_extra_timecourse_csv"),
            ("outputnode.NPR_prior_filename", "NPR_prior_filename"),
            ("outputnode.NPR_extra_filename", "NPR_extra_filename"),
            ("outputnode.NPR_optimize_report", "NPR_optimize_report"),
            ]),
        ])
    
    analysis_wf.inputs.group_inputnode.commonspace_template = split_dict[split_name_list[0]]["commonspace_resampled_template"]

    '''
    CONDITIONAL NODES AND CONNECTIONS
    '''

    # if inputs are in commmonspace, then analysis is computed with commonspace files
    # if --data_diagnosis, then commonspace files are needed for plotting results
    if not cr_opts.nativespace_analysis or analysis_opts.data_diagnosis:
        load_maps_dict_common_node = pe.Node(Function(input_names=['mask_file', 'WM_mask_file', 'CSF_mask_file', 'atlas_file', 'atlas_ref', 'anat_ref_file',
                                                            'seed_dict', 'prior_maps', 'transform_list','inverse_list', 'name_source', 'interpolation', 'rabies_data_type'],
                                            output_names=[
                                                'maps_dict_file'],
                                        function=load_maps_dict),
                                name='load_maps_dict_common_node')
        load_maps_dict_common_node.inputs.atlas_ref = preprocess_opts.labels
        load_maps_dict_common_node.inputs.seed_dict = seed_dict
        load_maps_dict_common_node.inputs.prior_maps = os.path.abspath(analysis_opts.prior_maps)
        load_maps_dict_common_node.inputs.interpolation = analysis_opts.interpolation
        load_maps_dict_common_node.inputs.rabies_data_type = analysis_opts.data_type

        split_name = split_name_list[0] # can take the commonspace files from any subject, they are all identical
        load_maps_dict_common_node.inputs.mask_file = split_dict[split_name]["commonspace_mask"]
        load_maps_dict_common_node.inputs.WM_mask_file = split_dict[split_name]["commonspace_WM_mask"]
        load_maps_dict_common_node.inputs.CSF_mask_file = split_dict[split_name]["commonspace_CSF_mask"]
        load_maps_dict_common_node.inputs.atlas_file = split_dict[split_name]["commonspace_labels"]
        load_maps_dict_common_node.inputs.anat_ref_file = split_dict[split_name]["commonspace_resampled_template"]
        load_maps_dict_common_node.inputs.transform_list = []
        load_maps_dict_common_node.inputs.inverse_list = []
        load_maps_dict_common_node.inputs.name_source = 'commonspace.nii'

        # only if inputs are in commonspace then analysis computations are also in commonspace
        if not cr_opts.nativespace_analysis:
            workflow.connect([
                (load_maps_dict_common_node, analysis_wf, [
                    ("maps_dict_file", "subject_inputnode.maps_dict_file"),
                    ]),
                (load_maps_dict_common_node, load_CR_dict_node, [
                    ("maps_dict_file", "maps_dict_file"),
                    ]),
                ])

    # if inputs are in native space, then subject specific maps are loaded
    if cr_opts.nativespace_analysis:
        load_maps_dict_native_node = pe.Node(Function(input_names=['mask_file', 'WM_mask_file', 'CSF_mask_file', 'atlas_file', 'atlas_ref', 'anat_ref_file', 
                                                            'seed_dict', 'prior_maps', 'transform_list','inverse_list', 'name_source', 'interpolation', 'rabies_data_type'],
                                            output_names=[
                                                'maps_dict_file'],
                                        function=load_maps_dict),
                                name='load_maps_dict_native_node')
        load_maps_dict_native_node.inputs.atlas_ref = preprocess_opts.labels
        load_maps_dict_native_node.inputs.seed_dict = seed_dict
        load_maps_dict_native_node.inputs.prior_maps = os.path.abspath(analysis_opts.prior_maps)
        load_maps_dict_native_node.inputs.interpolation = analysis_opts.interpolation
        load_maps_dict_native_node.inputs.rabies_data_type = analysis_opts.data_type

        workflow.connect([
            (conf_outputnode, load_maps_dict_native_node, [
                ("native_brain_mask", "mask_file"),
                ("native_WM_mask", "WM_mask_file"),
                ("native_CSF_mask", "CSF_mask_file"),
                ("native_labels", "atlas_file"),
                ("native_bold_ref", "anat_ref_file"),
                ("commonspace_to_native_transform_list", "transform_list"),
                ("commonspace_to_native_inverse_list", "inverse_list"),
                ("input_bold", "name_source"),
                ]),
            (load_maps_dict_native_node, analysis_wf, [
                ("maps_dict_file", "subject_inputnode.maps_dict_file"),
                ]),
            (load_maps_dict_native_node, load_CR_dict_node, [
                ("maps_dict_file", "maps_dict_file"),
                ]),
            (conf_outputnode, analysis_wf, [
                ("native_to_commonspace_transform_list", "subject_inputnode.native_to_commonspace_transform_list"),
                ("native_to_commonspace_inverse_list", "subject_inputnode.native_to_commonspace_inverse_list"),
                ]),
            ])

    if not cr_opts.nativespace_analysis or cr_opts.resample_to_commonspace:
        # this node can only exist if commonspace timeseries were previously generated
        analysis_split_joinnode = pe.JoinNode(niu.IdentityInterface(fields=['file_list', 'mask_file']),
                                            name='analysis_split_joinnode',
                                            joinsource='main_split',
                                            joinfield=['file_list'])

        workflow.connect([
            (conf_outputnode, analysis_split_joinnode, [
                ("commonspace_mask", "mask_file"),
                ]),
            (analysis_split_joinnode, analysis_wf, [
                ("file_list", "group_inputnode.bold_file_list"),
                ("mask_file", "group_inputnode.commonspace_mask"),
                ]),
            ])

        if not cr_opts.nativespace_analysis:
            workflow.connect([
                (conf_outputnode, analysis_split_joinnode, [
                    ("cleaned_path", "file_list"),
                    ]),
                ])
        elif cr_opts.resample_to_commonspace:
            workflow.connect([
                (conf_outputnode, analysis_split_joinnode, [
                    ("cleaned_timeseries_commonspace", "file_list"),
                    ]),
                ])

    # handling the resampling of nativespace computations into commonspace outputs        
    if cr_opts.nativespace_analysis and analysis_opts.resample_to_commonspace:
        workflow.connect([
            (analysis_wf, commonspace_analysis_datasink, [
                ("outputnode.corr_map_file_resampled", "seed_correlation_maps_resampled"),
                ]),
            (analysis_wf, commonspace_analysis_datasink, [
                ("outputnode.DR_nii_file_resampled", "dual_regression_nii_resampled"),
                ]),
            ])

    if analysis_opts.data_diagnosis:

        def prep_analysis_dict(seed_map_files, seed_timecourse_csv, dual_regression_nii, dual_regression_timecourse_csv, NPR_prior_timecourse_csv, NPR_extra_timecourse_csv, NPR_prior_filename, NPR_extra_filename):
            return {'seed_map_files':seed_map_files, 'seed_timecourse_csv':seed_timecourse_csv, 'dual_regression_nii':dual_regression_nii, 'dual_regression_timecourse_csv':dual_regression_timecourse_csv, 
                    'NPR_prior_timecourse_csv':NPR_prior_timecourse_csv, 'NPR_extra_timecourse_csv':NPR_extra_timecourse_csv, 
                    'NPR_prior_filename':NPR_prior_filename, 'NPR_extra_filename':NPR_extra_filename}
        prep_analysis_dict_node = pe.Node(Function(input_names=['seed_map_files', 'seed_timecourse_csv', 'dual_regression_nii', 'dual_regression_timecourse_csv', 'NPR_prior_timecourse_csv', 'NPR_extra_timecourse_csv', 'NPR_prior_filename', 'NPR_extra_filename'],
                                            output_names=[
                                                'analysis_dict'],
                                        function=prep_analysis_dict),
                                name='analysis_prep_analysis_dict')

        diagnosis_wf = init_diagnosis_wf(analysis_opts, cr_opts.nativespace_analysis, preprocess_opts, split_name_list, name="diagnosis_wf")

        workflow.connect([
            (load_CR_dict_node, diagnosis_wf, [
                ("CR_dict_file", "inputnode.CR_dict_file"),
                ]),
            (load_maps_dict_common_node, diagnosis_wf, [
                ("maps_dict_file", "inputnode.common_maps_dict_file"),
                ]),
            (analysis_wf, prep_analysis_dict_node, [
                ("outputnode.dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
                ]),
            (prep_analysis_dict_node, diagnosis_wf, [
                ("analysis_dict", "inputnode.analysis_dict"),
                ]),
            (diagnosis_wf, data_diagnosis_datasink, [
                ("outputnode.figure_temporal_diagnosis", "figure_temporal_diagnosis"),
                ("outputnode.figure_spatial_diagnosis", "figure_spatial_diagnosis"),
                ("outputnode.analysis_QC", "analysis_QC"),
                ("outputnode.temporal_info_csv", "temporal_info_csv"),
                ("outputnode.spatial_VE_nii", "spatial_VE_nii"),
                ("outputnode.temporal_std_nii", "temporal_std_nii"),
                ("outputnode.CR_prediction_std_nii", "CR_prediction_std_nii"),
                ("outputnode.GS_cov_nii", "GS_cov_nii"),
                ]),
            ])
        if cr_opts.nativespace_analysis:
            workflow.connect([
                (load_maps_dict_native_node, diagnosis_wf, [
                    ("maps_dict_file", "inputnode.sub_maps_dict_file"),
                    ]),
                (analysis_wf, prep_analysis_dict_node, [
                    ("outputnode.DR_nii_file_resampled", "dual_regression_nii"),
                    ]),
                (conf_outputnode, diagnosis_wf, [
                    ("native_to_commonspace_transform_list", "inputnode.native_to_commonspace_transform_list"),
                    ("native_to_commonspace_inverse_list", "inputnode.native_to_commonspace_inverse_list"),
                    ]),
                ])
        else:
            workflow.connect([
                (load_maps_dict_common_node, diagnosis_wf, [
                    ("maps_dict_file", "inputnode.sub_maps_dict_file"),
                    ]),
                (analysis_wf, prep_analysis_dict_node, [
                    ("outputnode.DR_nii_file", "dual_regression_nii"),
                    ]),
                ])

        if (analysis_opts.NPR_temporal_comp>-1) or (analysis_opts.NPR_spatial_comp>-1) or analysis_opts.optimize_NPR['apply']:
            workflow.connect([
                (analysis_wf, prep_analysis_dict_node, [
                    ("outputnode.NPR_prior_timecourse_csv", "NPR_prior_timecourse_csv"),
                    ("outputnode.NPR_extra_timecourse_csv", "NPR_extra_timecourse_csv"),
                    ("outputnode.NPR_prior_filename", "NPR_prior_filename"),
                    ("outputnode.NPR_extra_filename", "NPR_extra_filename"),
                    ]),
                ])
        else:
            prep_analysis_dict_node.inputs.NPR_prior_timecourse_csv = None
            prep_analysis_dict_node.inputs.NPR_extra_timecourse_csv = None
            prep_analysis_dict_node.inputs.NPR_prior_filename = None
            prep_analysis_dict_node.inputs.NPR_extra_filename = None

        if len(analysis_opts.seed_list) > 0:
            workflow.connect([
                (analysis_wf, prep_analysis_dict_node, [
                    ("outputnode.joined_corr_map_file", "seed_map_files"),
                    ("outputnode.joined_seed_timecourse_csv", "seed_timecourse_csv"),
                    ]),
                ])
        else:
            prep_analysis_dict_node.inputs.seed_map_files = []
            prep_analysis_dict_node.inputs.seed_timecourse_csv = None

    return workflow


# this function handles masks/maps that can be either common across subjects in commonspace, or resampled into individual spaces for nativespace analyses
# in the case of commonspace, this node only runs once, hence saving computations
# anat_ref_file is an anatomical 3D image that defines the resampling space inputted to antsApplyTransform
def load_maps_dict(mask_file, WM_mask_file, CSF_mask_file, atlas_file, atlas_ref, anat_ref_file,
                   seed_dict, prior_maps, transform_list, inverse_list, name_source, interpolation, rabies_data_type):
    import pickle
    import numpy as np
    import SimpleITK as sitk
    import os
    import pathlib  # Better path manipulation
    from rabies.utils import applyTransforms_4D, applyTransforms_3D
    from rabies.analysis_pkg.utils import compute_edge_mask
    mask_img = sitk.ReadImage(mask_file)
    mask_array = sitk.GetArrayFromImage(mask_img)
    volume_indices = mask_array.astype(bool)
    
    if WM_mask_file is None:
        WM_idx = None
    else:
        WM_idx = sitk.GetArrayFromImage(
            sitk.ReadImage(WM_mask_file))[volume_indices].astype(bool)
    if CSF_mask_file is None:
        CSF_idx = None
    else:
        CSF_idx = sitk.GetArrayFromImage(
            sitk.ReadImage(CSF_mask_file))[volume_indices].astype(bool)
    if atlas_file is None:
        atlas_idx = None
    else:
        atlas_idx = sitk.GetArrayFromImage(sitk.ReadImage(atlas_file))[volume_indices]

    edge_idx = compute_edge_mask(mask_array, num_edge_voxels=1)[volume_indices]

    resampled_4D_img = applyTransforms_4D(in_img=prior_maps, ref_file=anat_ref_file, transforms_3D=transform_list, inverses_3D=inverse_list, 
                                          motcorr_affine_list=None, interpolation=interpolation, rabies_data_type=rabies_data_type, clip_negative=False)
    resampled_maps = sitk.GetArrayFromImage(resampled_4D_img)
    prior_map_vectors = resampled_maps[:,volume_indices] # we return the 2D format of map number by voxels

    # resample each seed into the right space and then load the seed as an array
    seed_arr_dict = {}
    seed_name_list = list(seed_dict.keys())
    for seed_name in seed_name_list:
        seed_file = seed_dict[seed_name]
        resampled_seed_file = os.path.abspath(f'{seed_name}_resampled.nii.gz')
        applyTransforms_3D(transforms = transform_list, inverses = inverse_list, 
                        input_image = seed_file, ref_image = anat_ref_file, output_filename = resampled_seed_file, interpolation='GenericLabel', rabies_data_type=sitk.sitkInt16, clip_negative=False)
        seed_arr_dict[seed_name] = sitk.GetArrayFromImage(sitk.ReadImage(resampled_seed_file))[volume_indices]

    # prepare the list ROI numbers from the atlas for FC matrices
    # get the complete set of original ROI integers 
    if atlas_ref is None:
        roi_list=None
    else:
        atlas_data = sitk.GetArrayFromImage(sitk.ReadImage(str(atlas_ref))).astype(int)
        roi_list = []
        for i in range(1, atlas_data.max()+1):
            if np.max(i == atlas_data):  # include integers that do have labelled voxels
                roi_list.append(i)

    maps_dict = {'mask_file':mask_file, 'volume_indices':volume_indices, 'WM_idx':WM_idx, 'CSF_idx':CSF_idx,
                 'atlas_idx':atlas_idx, 'edge_idx':edge_idx, 'anat_ref_file':anat_ref_file,
                 'seed_arr_dict':seed_arr_dict, 'prior_map_vectors':prior_map_vectors, 'roi_list':roi_list}
    
    filename_split = pathlib.Path(name_source).name.rsplit(".nii")
    maps_dict_file = os.path.abspath(f'{filename_split[0]}_maps_dict.pkl')
    with open(maps_dict_file, 'wb') as handle:
        pickle.dump(maps_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return maps_dict_file


# this function loads subject-specific data from confound correction
def load_CR_input_dict(maps_dict_file, cleaned_bold_file, CR_data_dict, VE_file, STD_file, CR_STD_file, name_source):
    import pickle
    import pathlib
    import os
    import numpy as np
    import SimpleITK as sitk

    with open(maps_dict_file, 'rb') as handle:
        maps_dict = pickle.load(handle)

    volume_indices = maps_dict['volume_indices']

    data_img = sitk.ReadImage(cleaned_bold_file)
    data_array = sitk.GetArrayFromImage(data_img)
    num_volumes = data_array.shape[0]
    timeseries = np.zeros([num_volumes, volume_indices.sum()])
    for i in range(num_volumes):
        timeseries[i, :] = (data_array[i, :, :, :])[volume_indices]

    VE_spatial = sitk.GetArrayFromImage(
        sitk.ReadImage(VE_file))[volume_indices]
    temporal_std = sitk.GetArrayFromImage(
        sitk.ReadImage(STD_file))[volume_indices]
    predicted_std = sitk.GetArrayFromImage(
        sitk.ReadImage(CR_STD_file))[volume_indices]

    sub_dict = {'bold_file':cleaned_bold_file, 'name_source':name_source, 'CR_data_dict':CR_data_dict, 
            'timeseries':timeseries, 'VE_spatial':VE_spatial, 'temporal_std':temporal_std,
            'predicted_std':predicted_std}
    
    filename_split = pathlib.Path(name_source).name.rsplit(".nii")
    CR_dict_file = os.path.abspath(f'{filename_split[0]}_data_dict.pkl')
    with open(CR_dict_file, 'wb') as handle:
        pickle.dump(sub_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return CR_dict_file


def read_confound_workflow(conf_output, cr_opts):
    from nipype import logging
    log = logging.getLogger('nipype.workflow')

    conf_workflow_file = f'{conf_output}/rabies_confound_correction_workflow.pkl'

    node_dict = get_workflow_dict(conf_workflow_file)

    preproc_outputnode_name = 'confound_correction_main_wf.preproc_outputnode'
    match_targets = {'input_bold':[preproc_outputnode_name, 'input_bold'],
                    'commonspace_bold':[preproc_outputnode_name, 'commonspace_bold'],
                    'commonspace_mask':[preproc_outputnode_name, 'commonspace_mask'],
                    'commonspace_WM_mask':[preproc_outputnode_name, 'commonspace_WM_mask'],
                    'commonspace_CSF_mask':[preproc_outputnode_name, 'commonspace_CSF_mask'],
                    'commonspace_vascular_mask':[preproc_outputnode_name, 'commonspace_vascular_mask'],
                    'commonspace_labels':[preproc_outputnode_name, 'commonspace_labels'],
                    'motion_params_csv':[preproc_outputnode_name, 'motion_params_csv'],
                    'FD_csv':[preproc_outputnode_name, 'FD_csv'],
                    'commonspace_resampled_template':[preproc_outputnode_name, 'commonspace_resampled_template'],
                    'cleaned_path':['confound_correction_main_wf.confound_correction_wf.regress', 'cleaned_path'],
                    'data_dict':['confound_correction_main_wf.confound_correction_wf.regress', 'data_dict'],
                    'VE_file_path':['confound_correction_main_wf.confound_correction_wf.regress', 'VE_file_path'],
                    'STD_file_path':['confound_correction_main_wf.confound_correction_wf.regress', 'STD_file_path'],
                    'CR_STD_file_path':['confound_correction_main_wf.confound_correction_wf.regress', 'CR_STD_file_path'],
                    'random_CR_STD_file_path':['confound_correction_main_wf.confound_correction_wf.regress', 'random_CR_STD_file_path'],
                    'corrected_CR_STD_file_path':['confound_correction_main_wf.confound_correction_wf.regress', 'corrected_CR_STD_file_path'],
                    }

    if cr_opts.nativespace_analysis:
        match_targets.update({'native_bold':[preproc_outputnode_name, 'native_bold'],
                        'native_bold_ref':[preproc_outputnode_name, 'native_bold_ref'],
                        'native_brain_mask':[preproc_outputnode_name, 'native_brain_mask'],
                        'native_WM_mask':[preproc_outputnode_name, 'native_WM_mask'],
                        'native_CSF_mask':[preproc_outputnode_name, 'native_CSF_mask'],
                        'native_labels':[preproc_outputnode_name, 'native_labels'],
                        'commonspace_to_native_transform_list':[preproc_outputnode_name, 'commonspace_to_native_transform_list'],
                        'commonspace_to_native_inverse_list':[preproc_outputnode_name, 'commonspace_to_native_inverse_list'],
                        'native_to_commonspace_transform_list':[preproc_outputnode_name, 'native_to_commonspace_transform_list'],
                        'native_to_commonspace_inverse_list':[preproc_outputnode_name, 'native_to_commonspace_inverse_list'],
                        })
        if cr_opts.resample_to_commonspace:
            match_targets.update({'cleaned_timeseries_commonspace':['confound_correction_main_wf.cleaned_bold_to_commonspace', 'resampled_file'],
                        })

    split_dict = {}
    split_name = []
    # preparing a new iterative node where each BOLD scan is a different split
    [unit_bold, output_bold] = match_targets['input_bold']
    bold_dict = node_dict[unit_bold]
    # fill each BOLD scan split with proper affiliated outputs from preprocessing
    fill_split_dict(bold_dict, output_bold, split_name, split_dict, [], node_dict, match_targets)

    target_list = list(match_targets.keys())

    # don't include scans that were removed during confound correction
    corrected_split_name=[]
    remove_list = []
    import pathlib
    for name in split_name:
        filename = pathlib.Path(split_dict[name]['cleaned_path']).name
        if 'empty' in filename:
            remove_list.append(name)
            del split_dict[name]
        else:
            corrected_split_name.append(name)
    split_name = corrected_split_name

    if len(remove_list)>0:
        scan_list_str = ''
        for name in remove_list:
            scan_list_str += f'\n        - {name}'
        log.warning(f"""
        The following scans were not included for analysis as the file was empty: {scan_list_str}
        This is likely due to not meeting the minimum_timepoints threshold from --frame_censoring.
                    """)

    return split_dict, split_name, target_list


def find_split(scan, split_name):
    for split in split_name:
        if split in scan:
            return split
    raise ValueError(f"No previous file name is matching {scan}")