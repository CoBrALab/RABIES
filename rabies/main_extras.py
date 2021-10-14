import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function


def detached_confound_regression_wf(preprocess_opts, cr_opts, analysis_opts):

    workflow = pe.Workflow(name=cr_opts.output_name)

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

    workflow, confound_regression_wf = integrate_confound_regression(workflow, outputnode, cr_opts, preprocess_opts.bold_only)

    # Integrate analysis
    if analysis_opts is not None:
        workflow = integrate_analysis(
            workflow, outputnode, confound_regression_wf, analysis_opts, True, not cr_opts.nativespace_analysis, bold_scan_list, preprocess_opts)

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

    if bold_only and cr_opts.nativespace_analysis:
        raise ValueError(
            'Must not select --nativespace_analysis option for running confound regression on outputs from --bold_only.')

    if cr_opts.nativespace_analysis:
        workflow.connect([
            (outputnode, confound_regression_wf, [
                ("native_bold", "inputnode.bold_file"),
                ("native_brain_mask", "inputnode.brain_mask"),
                ("native_CSF_mask", "inputnode.csf_mask"),
                ]),
            ])
    else:
        workflow.connect([
            (outputnode, confound_regression_wf, [
                ("commonspace_bold", "inputnode.bold_file"),
                ("commonspace_mask", "inputnode.brain_mask"),
                ("commonspace_CSF_mask", "inputnode.csf_mask"),
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
        if cr_opts.DVARS_censoring or cr_opts.FD_censoring:
            workflow.connect([
                (confound_regression_wf, confound_regression_datasink, [
                    ("outputnode.frame_mask_file", "frame_censoring_mask"),
                    ]),
                ])

    return workflow, confound_regression_wf


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
            raise ValueError("Cannot currently select --nativespace_analysis for running data_diagnosis")

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

    split_dict = {}
    bold_scan_list = get_files_from_tree(f'{preproc_output}/bold_datasink/input_bold')
    split_name = []
    for f in bold_scan_list:
        name = pathlib.Path(f).name.rsplit(".nii")[0]
        split_name.append(name)
        split_dict[name]={}
        split_dict[name]['commonspace_resampled_template']=f'{preproc_output}/bold_datasink/commonspace_resampled_template/resampled_template.nii.gz'

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
