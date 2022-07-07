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

    workflow = pe.Workflow(name='analysis_main_wf')

    conf_output = os.path.abspath(str(analysis_opts.confound_correction_out))

    split_dict, split_name_list, target_list = read_confound_workflow(conf_output, nativespace=cr_opts.nativespace_analysis)

    # update split_name according to the --scan_list option
    split_name_list = get_iterable_scan_list(analysis_opts.scan_list, split_name_list)

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

    workflow.connect([
        (main_split, conf_outputnode, [
            ("split_name", "split_name"),
            ]),
        ])


    # prepare analysis workflow
    analysis_output = os.path.abspath(str(analysis_opts.output_dir))
    commonspace_bold = not cr_opts.nativespace_analysis
    analysis_wf = init_analysis_wf(
        opts=analysis_opts, commonspace_cr=commonspace_bold)

    # prepare analysis datasink
    analysis_datasink = pe.Node(DataSink(base_directory=analysis_output,
                                         container="analysis_datasink"),
                                name="analysis_datasink")
    if analysis_opts.FC_matrix:
        if os.path.isfile(str(analysis_opts.ROI_csv)):
            analysis_datasink.inputs.matrix_ROI_csv = str(analysis_opts.ROI_csv)

    data_diagnosis_datasink = pe.Node(DataSink(base_directory=analysis_output,
                                         container="data_diagnosis_datasink"),
                                name="data_diagnosis_datasink")

    load_maps_dict_node = pe.Node(Function(input_names=['mask_file', 'WM_mask_file', 'CSF_mask_file', 'atlas_file', 'atlas_ref', 'preprocess_anat_template', 'prior_maps', 'transform_list','inverse_list'],
                                           output_names=[
                                               'maps_dict'],
                                       function=load_maps_dict),
                              name='load_maps_dict_node')
    load_maps_dict_node.inputs.atlas_ref = str(preprocess_opts.labels)
    if not os.path.isfile(str(analysis_opts.prior_maps)):
        raise ValueError("--prior_maps doesn't exists.")
    else:
        load_maps_dict_node.inputs.prior_maps = os.path.abspath(analysis_opts.prior_maps)


    if commonspace_bold or preprocess_opts.bold_only:
        split_name = split_name_list[0] # can take the commonspace files from any subject, they are all identical
        load_maps_dict_node.inputs.mask_file = split_dict[split_name]["commonspace_mask"]
        load_maps_dict_node.inputs.WM_mask_file = split_dict[split_name]["commonspace_WM_mask"]
        load_maps_dict_node.inputs.CSF_mask_file = split_dict[split_name]["commonspace_CSF_mask"]
        load_maps_dict_node.inputs.atlas_file = split_dict[split_name]["commonspace_labels"]
        load_maps_dict_node.inputs.preprocess_anat_template = split_dict[split_name]["commonspace_resampled_template"]
        load_maps_dict_node.inputs.transform_list = []
        load_maps_dict_node.inputs.inverse_list = []
    else:
        workflow.connect([
            (conf_outputnode, load_maps_dict_node, [
                ("native_brain_mask", "mask_file"),
                ("native_WM_mask", "WM_mask_file"),
                ("native_CSF_mask", "CSF_mask_file"),
                ("native_labels", "atlas_file"),
                ("anat_preproc", "preprocess_anat_template"),
                ("commonspace_to_native_transform_list", "transform_list"),
                ("commonspace_to_native_inverse_list", "inverse_list"),
                ]),
            ])


    load_sub_dict_node = pe.Node(Function(input_names=['maps_dict', 'bold_file', 'CR_data_dict', 'VE_file', 'STD_file', 'CR_STD_file', 'random_CR_STD_file', 'corrected_CR_STD_file', 'name_source'],
                                           output_names=[
                                               'dict_file'],
                                       function=load_sub_input_dict),
                              name='load_sub_dict_node')

    workflow.connect([
        (conf_outputnode, load_sub_dict_node, [
            ("cleaned_path", "bold_file"),
            ("data_dict", "CR_data_dict"),
            ("VE_file_path", "VE_file"),
            ("STD_file_path", "STD_file"),
            ("CR_STD_file_path", "CR_STD_file"),
            ("random_CR_STD_file_path", "random_CR_STD_file"),
            ("corrected_CR_STD_file_path", "corrected_CR_STD_file"),
            ("input_bold", "name_source"),
            ]),
        (load_maps_dict_node, load_sub_dict_node, [
            ("maps_dict", "maps_dict"),
            ]),
        ])


    analysis_split_joinnode = pe.JoinNode(niu.IdentityInterface(fields=['file_list', 'mask_file']),
                                         name='analysis_split_joinnode',
                                         joinsource='main_split',
                                         joinfield=['file_list'])

    if commonspace_bold or preprocess_opts.bold_only:
        workflow.connect([
            (conf_outputnode, analysis_split_joinnode, [
                ("commonspace_mask", "mask_file"),
                ]),
            ])
    else:
        workflow.connect([
            (conf_outputnode, analysis_split_joinnode, [
                ("native_brain_mask", "mask_file"),
                ]),
            ])

    workflow.connect([
        (conf_outputnode, analysis_split_joinnode, [
            ("cleaned_path", "file_list"),
            ]),
        (load_sub_dict_node, analysis_wf, [
            ("dict_file", "subject_inputnode.dict_file"),
            ]),
        (analysis_split_joinnode, analysis_wf, [
            ("file_list", "group_inputnode.bold_file_list"),
            ("mask_file", "group_inputnode.commonspace_mask"),
            ]),
        (analysis_wf, analysis_datasink, [
            ("outputnode.group_ICA_dir", "group_ICA_dir"),
            ("outputnode.matrix_data_file", "matrix_data_file"),
            ("outputnode.matrix_fig", "matrix_fig"),
            ("outputnode.corr_map_file", "seed_correlation_maps"),
            ("outputnode.DR_nii_file", "dual_regression_nii"),
            ("outputnode.dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
            ("outputnode.NPR_prior_timecourse_csv", "NPR_prior_timecourse_csv"),
            ("outputnode.NPR_extra_timecourse_csv", "NPR_extra_timecourse_csv"),
            ("outputnode.NPR_prior_filename", "NPR_prior_filename"),
            ("outputnode.NPR_extra_filename", "NPR_extra_filename"),
            ]),
        ])

    if analysis_opts.data_diagnosis:

        def prep_analysis_dict(seed_map_files, dual_regression_nii, dual_regression_timecourse_csv, NPR_prior_timecourse_csv, NPR_extra_timecourse_csv, NPR_prior_filename, NPR_extra_filename):
            return {'seed_map_files':seed_map_files, 'dual_regression_nii':dual_regression_nii, 'dual_regression_timecourse_csv':dual_regression_timecourse_csv, 
                    'NPR_prior_timecourse_csv':NPR_prior_timecourse_csv, 'NPR_extra_timecourse_csv':NPR_extra_timecourse_csv, 
                    'NPR_prior_filename':NPR_prior_filename, 'NPR_extra_filename':NPR_extra_filename}
        prep_analysis_dict_node = pe.Node(Function(input_names=['seed_map_files', 'dual_regression_nii', 'dual_regression_timecourse_csv', 'NPR_prior_timecourse_csv', 'NPR_extra_timecourse_csv', 'NPR_prior_filename', 'NPR_extra_filename'],
                                            output_names=[
                                                'analysis_dict'],
                                        function=prep_analysis_dict),
                                name='analysis_prep_analysis_dict')

        diagnosis_wf = init_diagnosis_wf(analysis_opts, commonspace_bold, preprocess_opts, split_name_list, name="diagnosis_wf")

        workflow.connect([
            (load_sub_dict_node, diagnosis_wf, [
                ("dict_file", "inputnode.dict_file"),
                ]),
            (analysis_wf, prep_analysis_dict_node, [
                ("outputnode.DR_nii_file", "dual_regression_nii"),
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
                ("outputnode.random_CR_std_nii", "random_CR_std_nii"),
                ("outputnode.corrected_CR_std_nii", "corrected_CR_std_nii"),
                ("outputnode.GS_cov_nii", "GS_cov_nii"),
                ]),
            ])
        if (analysis_opts.NPR_temporal_comp>-1) or (analysis_opts.NPR_spatial_comp>-1):
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
                    ]),
                ])
        else:
            prep_analysis_dict_node.inputs.seed_map_files = []

    return workflow


# this function handles masks/maps that can be either common across subjects in commonspace, or resampled into individual spaces
def load_maps_dict(mask_file, WM_mask_file, CSF_mask_file, atlas_file, atlas_ref, preprocess_anat_template, prior_maps, transform_list, inverse_list):
    import numpy as np
    import SimpleITK as sitk
    import os
    import pathlib  # Better path manipulation
    from rabies.utils import resample_image_spacing
    from rabies.analysis_pkg.utils import resample_prior_maps, compute_edge_mask
    mask_img = sitk.ReadImage(mask_file)
    mask_array = sitk.GetArrayFromImage(mask_img)
    volume_indices = mask_array.astype(bool)

    WM_idx = sitk.GetArrayFromImage(
        sitk.ReadImage(WM_mask_file))[volume_indices].astype(bool)
    CSF_idx = sitk.GetArrayFromImage(
        sitk.ReadImage(CSF_mask_file))[volume_indices].astype(bool)
    atlas_idx = sitk.GetArrayFromImage(sitk.ReadImage(atlas_file))[volume_indices]

    edge_idx = compute_edge_mask(mask_array, num_edge_voxels=1)[volume_indices]

    # the reference anatomical image (either the native space anat scan or commonspace template) is resampled to match the EPI resolution for plotting during --data_diagnosis
    resampled = resample_image_spacing(sitk.ReadImage(preprocess_anat_template), mask_img.GetSpacing())
    filename_split = pathlib.Path(preprocess_anat_template).name.rsplit(".nii")
    template_file = os.path.abspath(f'{filename_split[0]}_display_template.nii.gz')
    sitk.WriteImage(resampled, template_file)

    resampled_maps = resample_prior_maps(prior_maps, mask_file, transforms = transform_list, inverses = inverse_list)
    prior_map_vectors = resampled_maps[:,volume_indices] # we return the 2D format of map number by voxels

    # prepare the list ROI numbers from the atlas for FC matrices
    # get the complete set of original ROI integers 
    atlas_data = sitk.GetArrayFromImage(sitk.ReadImage(atlas_ref)).astype(int)
    roi_list = []
    for i in range(1, atlas_data.max()+1):
        if np.max(i == atlas_data):  # include integers that do have labelled voxels
            roi_list.append(i)

    return {'mask_file':mask_file, 'volume_indices':volume_indices, 'WM_idx':WM_idx, 'CSF_idx':CSF_idx, 
            'atlas_idx':atlas_idx, 'edge_idx':edge_idx, 'template_file':template_file, 'prior_map_vectors':prior_map_vectors, 'roi_list':roi_list}


# this function loads subject-specific data
def load_sub_input_dict(maps_dict, bold_file, CR_data_dict, VE_file, STD_file, CR_STD_file, random_CR_STD_file, corrected_CR_STD_file, name_source):
    import pickle
    import pathlib
    import os
    import numpy as np
    import SimpleITK as sitk

    volume_indices = maps_dict['volume_indices']

    data_img = sitk.ReadImage(bold_file)
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
    random_CR_std = sitk.GetArrayFromImage(
        sitk.ReadImage(random_CR_STD_file))[volume_indices]
    corrected_CR_std = sitk.GetArrayFromImage(
        sitk.ReadImage(corrected_CR_STD_file))[volume_indices]

    sub_dict = {'bold_file':bold_file, 'name_source':name_source, 'CR_data_dict':CR_data_dict, 
            'timeseries':timeseries, 'VE_spatial':VE_spatial, 'temporal_std':temporal_std,
            'predicted_std':predicted_std, 'random_CR_std':random_CR_std, 
            'corrected_CR_std':corrected_CR_std}
    
    # add all the maps_dict into the sub_dict
    for k in maps_dict.keys():
        sub_dict[k] = maps_dict[k]

    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")
    dict_file = os.path.abspath(f'{filename_split[0]}_data_dict.pkl')
    with open(dict_file, 'wb') as handle:
        pickle.dump(sub_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return dict_file


def read_confound_workflow(conf_output, nativespace=False):

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
                    'confounds_csv':[preproc_outputnode_name, 'confounds_csv'],
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

    if nativespace:
        match_targets.update({'native_bold':[preproc_outputnode_name, 'native_bold'],
                        'native_brain_mask':[preproc_outputnode_name, 'native_brain_mask'],
                        'native_WM_mask':[preproc_outputnode_name, 'native_WM_mask'],
                        'native_CSF_mask':[preproc_outputnode_name, 'native_CSF_mask'],
                        'native_labels':[preproc_outputnode_name, 'native_labels'],
                        'anat_preproc':[preproc_outputnode_name, 'anat_preproc'],
                        'commonspace_to_native_transform_list':[preproc_outputnode_name, 'commonspace_to_native_transform_list'],
                        'commonspace_to_native_inverse_list':[preproc_outputnode_name, 'commonspace_to_native_inverse_list'],
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
    import pathlib
    for name in split_name:
        filename = pathlib.Path(split_dict[name]['cleaned_path']).name
        if 'empty' in filename:
            del split_dict[name]
        else:
            corrected_split_name.append(name)
    split_name = corrected_split_name

    return split_dict, split_name, target_list


def get_iterable_scan_list(scan_list, split_name):
    # prep the subset of scans on which the analysis will be run
    import numpy as np
    import pandas as pd
    if os.path.isfile(os.path.abspath(scan_list[0])):
        updated_split_name=[]
        if '.nii' in pathlib.Path(scan_list[0]).name:
            for scan in scan_list:
                updated_split_name.append(find_split(scan, split_name))
        else:
            # read the file as a .txt
            scan_list = np.array(pd.read_csv(os.path.abspath(scan_list[0]), header=None)).flatten()
            for scan in scan_list:
                updated_split_name.append(find_split(scan, split_name))
    elif scan_list[0]=='all':
        updated_split_name = split_name
    else:
        raise ValueError("The scan_list input had improper format.")
    return updated_split_name


def find_split(scan, split_name):
    for split in split_name:
        if split in scan:
            return split
    raise ValueError(f"No previous file name is matching {scan}")