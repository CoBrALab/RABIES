import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
from rabies.analysis_pkg.diagnosis_pkg.diagnosis_wf import init_diagnosis_wf
from rabies.analysis_pkg.analysis_wf import init_analysis_wf


def detached_confound_correction_wf(preprocess_opts, cr_opts):
    from rabies.conf_reg_pkg.confound_correction import init_confound_correction_wf

    workflow = pe.Workflow(name=cr_opts.output_name+'_main_post_wf')

    preproc_output = os.path.abspath(str(cr_opts.preprocess_out))

    if cr_opts.read_datasink:
        split_dict, split_name, target_list = read_preproc_datasinks(preproc_output, nativespace=cr_opts.nativespace_analysis)
    else:
        split_dict, split_name, target_list = read_preproc_workflow(preproc_output, nativespace=cr_opts.nativespace_analysis)

    # setting up iterables from the BOLD scan splits
    main_split = pe.Node(niu.IdentityInterface(fields=['split_name']),
                         name="main_split")
    main_split.iterables = [('split_name', split_name)]

    # set output node from preprocessing
    def read_dict(split_dict, split_name, target_list):
        return [split_dict[split_name][target] for target in target_list]
    preproc_outputnode = pe.Node(Function(input_names=['split_dict', 'split_name', 'target_list'],
                                           output_names=target_list,
                                       function=read_dict),
                              name='preproc_outputnode')
    preproc_outputnode.inputs.split_dict = split_dict
    preproc_outputnode.inputs.target_list = target_list

    # need to set a buffer function which will be holding the preproc_outputnode outputs, 
    # so that it is saved in the workflow graph and can be read later during analysis
    def buffer_outputnode(input_bold=None, commonspace_bold=None, commonspace_mask=None, commonspace_WM_mask=None,
        commonspace_CSF_mask=None, commonspace_vascular_mask=None, commonspace_labels=None, confounds_csv=None,
        FD_csv=None, FD_voxelwise=None, pos_voxelwise=None, commonspace_resampled_template=None, native_bold=None, 
        native_brain_mask=None, native_WM_mask=None, native_CSF_mask=None, native_labels=None):
        return
    buffer_outputnode_node = pe.Node(Function(input_names=target_list,
                                           output_names=[],
                                       function=buffer_outputnode),
                              name='buffer_outputnode')
    for target in target_list:
        workflow.connect([
            (preproc_outputnode, buffer_outputnode_node, [
                (target, target),
                ]),
            ])


    confound_correction_wf = init_confound_correction_wf(cr_opts=cr_opts)

    workflow.connect([
        (main_split, preproc_outputnode, [
            ("split_name", "split_name"),
            ]),
        (preproc_outputnode, confound_correction_wf, [
            ("confounds_csv", "inputnode.confounds_file"),  # confounds file
            ("FD_csv", "inputnode.FD_file"),
            ]),
        ])

    if preprocess_opts.bold_only and cr_opts.nativespace_analysis:
        raise ValueError(
            'Must not select --nativespace_analysis option for running confound regression on outputs from --bold_only.')

    if cr_opts.nativespace_analysis:
        workflow.connect([
            (preproc_outputnode, confound_correction_wf, [
                ("native_bold", "inputnode.bold_file"),
                ("native_brain_mask", "inputnode.brain_mask"),
                ("native_CSF_mask", "inputnode.csf_mask"),
                ]),
            ])
    else:
        workflow.connect([
            (preproc_outputnode, confound_correction_wf, [
                ("commonspace_bold", "inputnode.bold_file"),
                ("commonspace_mask", "inputnode.brain_mask"),
                ("commonspace_CSF_mask", "inputnode.csf_mask"),
                ]),
            ])

    cr_output = os.path.abspath(str(cr_opts.output_dir))

    confound_correction_datasink = pe.Node(DataSink(base_directory=cr_output,
                                                    container=cr_opts.output_name+"_datasink"),
                                            name=cr_opts.output_name+"_datasink")
    workflow.connect([
        (confound_correction_wf, confound_correction_datasink, [
            ("outputnode.cleaned_path", "cleaned_timeseries"),
            ]),
        ])
    if cr_opts.run_aroma:
        workflow.connect([
            (confound_correction_wf, confound_correction_datasink, [
                ("outputnode.aroma_out", "aroma_out"),
                ]),
            ])
    if cr_opts.DVARS_censoring or cr_opts.FD_censoring:
        workflow.connect([
            (confound_correction_wf, confound_correction_datasink, [
                ("outputnode.frame_mask_file", "frame_censoring_mask"),
                ]),
            ])

    return workflow


def detached_analysis_wf(preprocess_opts, cr_opts, analysis_opts):

    workflow = pe.Workflow(name=analysis_opts.output_name+'_main_post_wf')

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
                                         container=analysis_opts.output_name+"_datasink"),
                                name=analysis_opts.output_name+"_datasink")


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
                                name=analysis_opts.output_name+'_prep_analysis_dict')

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
            (diagnosis_wf, analysis_datasink, [
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


def read_preproc_datasinks(preproc_output, nativespace=False):
    import pathlib
    import glob

    template_file = glob.glob(f'{preproc_output}/bold_datasink/commonspace_resampled_template/*')
    if len(template_file)==1:
        template_file = template_file[0]
    else:
        raise ValueError(f"Multiple files were found in {preproc_output}/bold_datasink/commonspace_resampled_template/"
                        "but there should only be one template file.")

    split_dict = {}
    bold_scan_list = get_files_from_tree(f'{preproc_output}/bold_datasink/input_bold')
    split_name = []
    for f in bold_scan_list:
        name = pathlib.Path(f).name.rsplit(".nii")[0]
        split_name.append(name)
        split_dict[name]={}
        split_dict[name]['commonspace_resampled_template']=template_file

    directory_list = [['bold_datasink','input_bold'],
        ['bold_datasink','commonspace_bold'], ['bold_datasink','commonspace_mask'], ['bold_datasink','commonspace_WM_mask'],
        ['bold_datasink','commonspace_CSF_mask'], ['bold_datasink','commonspace_vascular_mask'], ['bold_datasink','commonspace_labels'],
        ['confounds_datasink','confounds_csv'], ['confounds_datasink','FD_voxelwise'], ['confounds_datasink','pos_voxelwise'], ['confounds_datasink','FD_csv']]

    if nativespace:
        directory_list+=[['bold_datasink','native_bold'], ['bold_datasink','native_brain_mask'],
            ['bold_datasink','native_WM_mask'], ['bold_datasink','native_CSF_mask'], ['bold_datasink','native_labels']]
        

    target_list=['commonspace_resampled_template']
    for datasink,target in directory_list:

        if not os.path.isdir(f'{preproc_output}/{datasink}/{target}'):
            raise ValueError(f"The directory {preproc_output}/{datasink}/{target} does not exist. Make sure that all required "
                "datasink outputs are available. If --bold_only was selected, there are no native space outputs available.")
        target_list.append(target)
        file_list = get_files_from_tree(f'{preproc_output}/{datasink}/{target}')
        for f in file_list:
            for split in split_name:
                if split in f:
                    split_dict[split][target]=f
                    break

    return split_dict, split_name, target_list


def get_files_from_tree(startpath):
    file_list=[]
    for root, dirs, files in os.walk(startpath):
        for f in files:
            file_list.append(f'{root}/{f}')
    return file_list


def read_preproc_workflow(preproc_output, nativespace=False):

    preproc_workflow_file = f'{preproc_output}/rabies_preprocess_workflow.pkl'

    node_dict = get_workflow_dict(preproc_workflow_file)

    match_targets = {'input_bold':['main_wf.bold_selectfiles', 'out_file'],
                    'commonspace_bold':['main_wf.bold_main_wf.bold_commonspace_trans_wf.merge', 'out_file'],
                    'commonspace_mask':['main_wf.bold_main_wf.bold_commonspace_trans_wf.Brain_mask_EPI', 'EPI_mask'],
                    'commonspace_WM_mask':['main_wf.bold_main_wf.bold_commonspace_trans_wf.WM_mask_EPI', 'EPI_mask'],
                    'commonspace_CSF_mask':['main_wf.bold_main_wf.bold_commonspace_trans_wf.CSF_mask_EPI', 'EPI_mask'],
                    'commonspace_vascular_mask':['main_wf.bold_main_wf.bold_commonspace_trans_wf.vascular_mask_EPI', 'EPI_mask'],
                    'commonspace_labels':['main_wf.bold_main_wf.bold_commonspace_trans_wf.prop_labels_EPI', 'EPI_mask'],
                    'confounds_csv':['main_wf.bold_main_wf.bold_confs_wf.estimate_confounds', 'confounds_csv'],
                    'FD_voxelwise':['main_wf.bold_main_wf.bold_confs_wf.estimate_confounds', 'FD_voxelwise'],
                    'pos_voxelwise':['main_wf.bold_main_wf.bold_confs_wf.estimate_confounds', 'pos_voxelwise'],
                    'FD_csv':['main_wf.bold_main_wf.bold_confs_wf.estimate_confounds', 'FD_csv'],
                    'commonspace_resampled_template':['main_wf.resample_template', 'resampled_template'],
                    }
    if nativespace:
        match_targets.update({'native_bold':['main_wf.bold_main_wf.bold_commonspace_trans_wf.merge', 'out_file'],
                        'native_brain_mask':['main_wf.bold_main_wf.bold_commonspace_trans_wf.Brain_mask_EPI', 'EPI_mask'],
                        'native_WM_mask':['main_wf.bold_main_wf.bold_commonspace_trans_wf.WM_mask_EPI', 'EPI_mask'],
                        'native_CSF_mask':['main_wf.bold_main_wf.bold_commonspace_trans_wf.CSF_mask_EPI', 'EPI_mask'],
                        'native_labels':['main_wf.bold_main_wf.bold_commonspace_trans_wf.prop_labels_EPI', 'EPI_mask'],
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


def read_confound_workflow(conf_output, cr_opts, nativespace=False):

    conf_workflow_file = f'{conf_output}/rabies_confound_correction_workflow.pkl'

    node_dict = get_workflow_dict(conf_workflow_file)

    preproc_outputnode_name = f'{cr_opts.output_name}_main_post_wf.preproc_outputnode'
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
                    'cleaned_path':[f'{cr_opts.output_name}_main_post_wf.confound_correction_wf.regress', 'cleaned_path'],
                    'data_dict':[f'{cr_opts.output_name}_main_post_wf.confound_correction_wf.regress', 'data_dict'],
                    'VE_file_path':[f'{cr_opts.output_name}_main_post_wf.confound_correction_wf.regress', 'VE_file_path'],
                    'STD_file_path':[f'{cr_opts.output_name}_main_post_wf.confound_correction_wf.regress', 'STD_file_path'],
                    'CR_STD_file_path':[f'{cr_opts.output_name}_main_post_wf.confound_correction_wf.regress', 'CR_STD_file_path'],
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


def fill_split_dict(d, output_bold, split_name, split_dict, keys, node_dict, match_targets):
    if isinstance(d, dict):
        for key in list(d.keys()):
            fill_split_dict(d[key], output_bold, split_name, split_dict, keys+[key], node_dict, match_targets)
    else:
        f = d.result.outputs.get()[output_bold]
        split = pathlib.Path(f).name.rsplit(".nii")[0]
        split_name.append(split)
        split_dict[split]={}
        target_list = list(match_targets.keys())
        for target in target_list:
            [unit, output] = match_targets[target]
            node = retrieve_node(node_dict[unit], keys)
            split_dict[split][target] = node.result.outputs.get()[output]
        
def retrieve_node(d, keys):
    if isinstance(d, dict):
        return retrieve_node(d[keys[0]], keys[1:])
    else:
        return d


def get_workflow_dict(workflow_file):
    import pickle
    with open(workflow_file, 'rb') as handle:
        graph = pickle.load(handle)
    
    node_list = list(graph.nodes)
    node_dict = {}
    for node in node_list:
        key_l = [node.fullname]+node.parameterization
        fill_node_dict(node_dict, key_l, node)
    return node_dict


def fill_node_dict(d, key_l, e):
    if len(key_l)>0:
        key = key_l[0]
        if not (key in list(d.keys())):
            d[key] = {}
        d[key] = fill_node_dict(d[key], key_l[1:], e)
        return d
    else:
        return e
