import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
from rabies.utils import fill_split_dict, get_workflow_dict


def init_main_confound_correction_wf(preprocess_opts, cr_opts):
    from rabies.confound_correction_pkg.confound_correction import init_confound_correction_wf

    workflow = pe.Workflow(name='confound_correction_main_wf')

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
                                                    container="confound_correction_datasink"),
                                            name="confound_correction_datasink")
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
