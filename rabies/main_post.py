import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
from rabies.analysis_pkg.diagnosis_pkg.diagnosis_wf import init_diagnosis_wf
from rabies.analysis_pkg.analysis_wf import init_analysis_wf


def detached_confound_correction_wf(preprocess_opts, cr_opts, analysis_opts):

    workflow = pe.Workflow(name=cr_opts.output_name+'_main_post_wf')

    preproc_output = os.path.abspath(str(cr_opts.preprocess_out))

    split_dict, split_name, target_list, bold_scan_list = read_preproc_datasinks(preproc_output, nativespace=cr_opts.nativespace_analysis)

    # setting up iterables
    main_split = pe.Node(niu.IdentityInterface(fields=['split_name']),
                         name="main_split")
    main_split.iterables = [('split_name', split_name)]

    # set output node from preprocessing
    def read_dict(split_dict, split_name, target_list):
        return [split_dict[split_name][target] for target in target_list]
    outputnode = pe.Node(Function(input_names=['split_dict', 'split_name', 'target_list'],
                                           output_names=target_list,
                                       function=read_dict),
                              name=cr_opts.output_name+'_outputnode')
    outputnode.inputs.split_dict = split_dict
    outputnode.inputs.target_list = target_list

    workflow.connect([
        (main_split, outputnode, [
            ("split_name", "split_name"),
            ]),
        ])

    workflow, confound_correction_wf = integrate_confound_correction(workflow, outputnode, cr_opts, preprocess_opts.bold_only)

    # Integrate analysis
    if analysis_opts is not None:
        workflow = integrate_analysis(
            workflow, outputnode, confound_correction_wf, analysis_opts, True, not cr_opts.nativespace_analysis, bold_scan_list, preprocess_opts)

    return workflow


def integrate_confound_correction(workflow, outputnode, cr_opts, bold_only):
    cr_output = os.path.abspath(str(cr_opts.output_dir))

    from rabies.conf_reg_pkg.confound_correction import init_confound_correction_wf
    confound_correction_wf = init_confound_correction_wf(cr_opts=cr_opts, name=cr_opts.output_name)

    workflow.connect([
        (outputnode, confound_correction_wf, [
            ("confounds_csv", "inputnode.confounds_file"),  # confounds file
            ("FD_csv", "inputnode.FD_file"),
            ]),
        ])

    if bold_only and cr_opts.nativespace_analysis:
        raise ValueError(
            'Must not select --nativespace_analysis option for running confound regression on outputs from --bold_only.')

    if cr_opts.nativespace_analysis:
        workflow.connect([
            (outputnode, confound_correction_wf, [
                ("native_bold", "inputnode.bold_file"),
                ("native_brain_mask", "inputnode.brain_mask"),
                ("native_CSF_mask", "inputnode.csf_mask"),
                ]),
            ])
    else:
        workflow.connect([
            (outputnode, confound_correction_wf, [
                ("commonspace_bold", "inputnode.bold_file"),
                ("commonspace_mask", "inputnode.brain_mask"),
                ("commonspace_CSF_mask", "inputnode.csf_mask"),
                ]),
            ])

    if cr_opts.rabies_step == 'confound_correction':
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

    return workflow, confound_correction_wf


def integrate_analysis(workflow, outputnode, confound_correction_wf, analysis_opts, bold_only, commonspace_bold, bold_scan_list, opts):
    analysis_output = os.path.abspath(str(analysis_opts.output_dir))

    analysis_wf = init_analysis_wf(
        opts=analysis_opts, commonspace_cr=commonspace_bold, name=analysis_opts.output_name)

    analysis_datasink = pe.Node(DataSink(base_directory=analysis_output,
                                         container=analysis_opts.output_name+"_datasink"),
                                name=analysis_opts.output_name+"_datasink")

    def prep_dict(bold_file, CR_data_dict, VE_file, STD_file, mask_file, WM_mask_file, CSF_mask_file, atlas_file, preprocess_anat_template, name_source):
        return {'bold_file':bold_file, 'CR_data_dict':CR_data_dict, 'VE_file':VE_file, 'STD_file':STD_file, 'mask_file':mask_file, 'WM_mask_file':WM_mask_file, 'CSF_mask_file':CSF_mask_file, 'atlas_file':atlas_file, 'preprocess_anat_template':preprocess_anat_template, 'name_source':name_source}
    prep_dict_node = pe.Node(Function(input_names=['bold_file', 'CR_data_dict', 'VE_file', 'STD_file', 'mask_file', 'WM_mask_file', 'CSF_mask_file', 'atlas_file', 'preprocess_anat_template', 'name_source'],
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

    workflow,find_iterable_node, joinnode_main,analysis_split,scan_split_name = transit_iterables(workflow, prep_dict_node, analysis_opts.scan_list, bold_only, bold_scan_list, node_prefix=analysis_opts.output_name)

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
        (confound_correction_wf, prep_dict_node, [
            ("outputnode.cleaned_path", "bold_file"),
            ("outputnode.CR_data_dict", "CR_data_dict"),
            ("outputnode.VE_file", "VE_file"),
            ("outputnode.STD_file", "STD_file"),
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

        diagnosis_wf = init_diagnosis_wf(analysis_opts, commonspace_bold, opts, analysis_split, scan_split_name, name=analysis_opts.output_name+"diagnosis_wf")

        workflow.connect([
            (analysis_wf, prep_analysis_dict_node, [
                ("outputnode.DR_nii_file", "dual_regression_nii"),
                ("outputnode.dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
                ]),
            (joinnode_main, diagnosis_wf, [
                ("file_list", "inputnode.mask_dict_list"),
                ]),
            (find_iterable_node, diagnosis_wf, [
                ("file", "inputnode.file_dict"),
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

    return workflow,find_iterable_node, joinnode_main,analysis_split, scan_split_name

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

    return split_dict, split_name, target_list, bold_scan_list


def get_files_from_tree(startpath):
    file_list=[]
    for root, dirs, files in os.walk(startpath):
        for f in files:
            file_list.append(f'{root}/{f}')
    return file_list
