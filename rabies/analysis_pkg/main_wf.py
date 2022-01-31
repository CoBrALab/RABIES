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

    split_dict, split_name, target_list = read_confound_workflow(conf_output, cr_opts, nativespace=cr_opts.nativespace_analysis)

    # setting up iterables from the BOLD scan splits
    main_split = pe.Node(niu.IdentityInterface(fields=['split_name']),
                         name="main_split")
    main_split.iterables = [('split_name', split_name)]

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


    analysis_output = os.path.abspath(str(analysis_opts.output_dir))
    commonspace_bold = not cr_opts.nativespace_analysis
    # prepare analysis workflow
    analysis_wf = init_analysis_wf(
        opts=analysis_opts, commonspace_cr=commonspace_bold)

    # prepare analysis datasink
    analysis_datasink = pe.Node(DataSink(base_directory=analysis_output,
                                         container="analysis_datasink"),
                                name="analysis_datasink")

    data_diagnosis_datasink = pe.Node(DataSink(base_directory=analysis_output,
                                         container="data_diagnosis_datasink"),
                                name="data_diagnosis_datasink")


    input_buffer_node = pe.Node(niu.IdentityInterface(fields=['bold_file', 'mask_file','atlas_file','WM_mask_file',
                         'CSF_mask_file', 'preprocess_anat_template']),
                         name="input_buffer")

    workflow.connect([
        (conf_outputnode, input_buffer_node, [
            ("cleaned_path", "bold_file"),
            ]),
        ])

    if commonspace_bold or preprocess_opts.bold_only:
        workflow.connect([
            (conf_outputnode, input_buffer_node, [
                ("commonspace_mask", "mask_file"),
                ("commonspace_labels", "atlas_file"),
                ("commonspace_WM_mask", "WM_mask_file"),
                ("commonspace_CSF_mask", "CSF_mask_file"),
                ("commonspace_resampled_template", "preprocess_anat_template"),
                ]),
            ])
    else:
        workflow.connect([
            (conf_outputnode, input_buffer_node, [
                ("native_brain_mask", "mask_file"),
                ("native_labels", "atlas_file"),
                ]),
            ])


    # concatenate the different inputs into a dictionary
    def prep_dict(bold_file, CR_data_dict, VE_file, STD_file, CR_STD_file, mask_file, WM_mask_file, CSF_mask_file, atlas_file, preprocess_anat_template, name_source):
        return {'bold_file':bold_file, 'CR_data_dict':CR_data_dict, 'VE_file':VE_file, 'STD_file':STD_file, 'CR_STD_file':CR_STD_file, 'mask_file':mask_file, 'WM_mask_file':WM_mask_file, 'CSF_mask_file':CSF_mask_file, 'atlas_file':atlas_file, 'preprocess_anat_template':preprocess_anat_template, 'name_source':name_source}
    prep_dict_node = pe.Node(Function(input_names=['bold_file', 'CR_data_dict', 'VE_file', 'STD_file', 'CR_STD_file', 'mask_file', 'WM_mask_file', 'CSF_mask_file', 'atlas_file', 'preprocess_anat_template', 'name_source'],
                                           output_names=[
                                               'prep_dict'],
                                       function=prep_dict),
                              name='prep_dict')
    workflow.connect([
        (conf_outputnode, prep_dict_node, [
            ("data_dict", "CR_data_dict"),
            ("VE_file_path", "VE_file"),
            ("STD_file_path", "STD_file"),
            ("CR_STD_file_path", "CR_STD_file"),
            ("input_bold", "name_source"),
            ]),
        (input_buffer_node, prep_dict_node, [
            ("bold_file", "bold_file"),
            ("mask_file", "mask_file"),
            ("WM_mask_file", "WM_mask_file"),
            ("CSF_mask_file", "CSF_mask_file"),
            ("atlas_file", "atlas_file"),
            ("preprocess_anat_template", "preprocess_anat_template"),
            ]),
        ])


    analysis_split_joinnode = pe.JoinNode(niu.IdentityInterface(fields=['file_list', 'mask_file', 'prep_dict_list']),
                                         name='analysis_split_joinnode',
                                         joinsource='main_split',
                                         joinfield=['file_list', 'prep_dict_list'])


    workflow.connect([
        (input_buffer_node, analysis_split_joinnode, [
            ("mask_file", "mask_file"),
            ("bold_file", "file_list"),
            ]),
        (input_buffer_node, analysis_wf, [
            ("bold_file", "subject_inputnode.bold_file"),
            ("mask_file", "subject_inputnode.mask_file"),
            ("atlas_file", "subject_inputnode.atlas_file"),
            ]),
        (prep_dict_node, analysis_split_joinnode, [
            ("prep_dict", "prep_dict_list"),
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
            ("outputnode.dual_ICA_timecourse_csv", "dual_ICA_timecourse_csv"),
            ("outputnode.dual_ICA_filename", "dual_ICA_filename"),
            ]),
        ])

    if analysis_opts.data_diagnosis:

        def prep_analysis_dict(seed_map_files, dual_regression_nii, dual_regression_timecourse_csv, dual_ICA_timecourse_csv, dual_ICA_filename):
            return {'seed_map_files':seed_map_files, 'dual_regression_nii':dual_regression_nii, 'dual_regression_timecourse_csv':dual_regression_timecourse_csv, 'dual_ICA_timecourse_csv':dual_ICA_timecourse_csv, 'dual_ICA_filename':dual_ICA_filename}
        prep_analysis_dict_node = pe.Node(Function(input_names=['seed_map_files', 'dual_regression_nii', 'dual_regression_timecourse_csv', 'dual_ICA_timecourse_csv', 'dual_ICA_filename'],
                                            output_names=[
                                                'analysis_dict'],
                                        function=prep_analysis_dict),
                                name='analysis_prep_analysis_dict')

        diagnosis_wf = init_diagnosis_wf(analysis_opts, commonspace_bold, preprocess_opts, split_name, name="diagnosis_wf")

        workflow.connect([
            (analysis_split_joinnode, diagnosis_wf, [
                ("prep_dict_list", "inputnode.mask_dict_list"),
                ]),
            (analysis_wf, prep_analysis_dict_node, [
                ("outputnode.DR_nii_file", "dual_regression_nii"),
                ("outputnode.dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
                ]),
            (prep_dict_node, diagnosis_wf, [
                ("prep_dict", "inputnode.file_dict"),
                ]),
            (prep_analysis_dict_node, diagnosis_wf, [
                ("analysis_dict", "inputnode.analysis_dict"),
                ]),
            (diagnosis_wf, data_diagnosis_datasink, [
                ("outputnode.figure_temporal_diagnosis", "figure_temporal_diagnosis"),
                ("outputnode.figure_spatial_diagnosis", "figure_spatial_diagnosis"),
                ("outputnode.dataset_diagnosis", "dataset_diagnosis"),
                ("outputnode.temporal_info_csv", "temporal_info_csv"),
                ("outputnode.spatial_VE_nii", "spatial_VE_nii"),
                ("outputnode.temporal_std_nii", "temporal_std_nii"),
                ("outputnode.GS_corr_nii", "GS_corr_nii"),
                ("outputnode.DVARS_corr_nii", "DVARS_corr_nii"),
                ("outputnode.FD_corr_nii", "FD_corr_nii"),
                ]),
            ])
        if analysis_opts.dual_ICA>0:
            workflow.connect([
                (analysis_wf, prep_analysis_dict_node, [
                    ("outputnode.dual_ICA_timecourse_csv", "dual_ICA_timecourse_csv"),
                    ("outputnode.dual_ICA_filename", "dual_ICA_filename"),
                    ]),
                ])
        else:
            prep_analysis_dict_node.inputs.dual_ICA_timecourse_csv = None
            prep_analysis_dict_node.inputs.dual_ICA_filename = None

        if len(analysis_opts.seed_list) > 0:
            workflow.connect([
                (analysis_wf, prep_analysis_dict_node, [
                    ("outputnode.joined_corr_map_file", "seed_map_files"),
                    ]),
                ])
        else:
            prep_analysis_dict_node.inputs.seed_map_files = []

    return workflow



def read_confound_workflow(conf_output, cr_opts, nativespace=False):

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
                    }

    if nativespace:
        match_targets.update({'native_bold':[preproc_outputnode_name, 'native_bold'],
                        'native_brain_mask':[preproc_outputnode_name, 'native_brain_mask'],
                        'native_WM_mask':[preproc_outputnode_name, 'native_WM_mask'],
                        'native_CSF_mask':[preproc_outputnode_name, 'native_CSF_mask'],
                        'native_labels':[preproc_outputnode_name, 'native_labels'],
                        })

    split_dict = {}
    split_name = []
    # preparing a new iterative node where each BOLD scan is a different split
    [unit_bold, output_bold] = match_targets['input_bold']
    bold_dict = node_dict[unit_bold]
    # fill each BOLD scan split with proper affiliated outputs from preprocessing
    fill_split_dict(bold_dict, output_bold, split_name, split_dict, [], node_dict, match_targets)

    target_list = list(match_targets.keys())

    return split_dict, split_name, target_list
