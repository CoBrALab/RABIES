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

    if preprocess_opts.bold_only and cr_opts.nativespace_analysis:
        raise ValueError(
            'Must not select --nativespace_analysis option for running confound regression on outputs from --bold_only.')

    if cr_opts.read_datasink:
        split_dict, split_name, target_list = read_preproc_datasinks(preproc_output, nativespace=cr_opts.nativespace_analysis, fast_commonspace=preprocess_opts.commonspace_reg['fast_commonspace'], atlas_reg_script=preprocess_opts.commonspace_reg['template_registration'], voxelwise_motion=preprocess_opts.voxelwise_motion)
    else:
        split_dict, split_name, target_list = read_preproc_workflow(preproc_output, nativespace=cr_opts.nativespace_analysis)

    # filter inclusion/exclusion lists
    from rabies.utils import filter_scan_inclusion, filter_scan_exclusion
    split_name = filter_scan_inclusion(cr_opts.inclusion_ids, split_name)
    split_name = filter_scan_exclusion(cr_opts.exclusion_ids, split_name)

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
        commonspace_CSF_mask=None, commonspace_vascular_mask=None, commonspace_labels=None, motion_params_csv=None,
        FD_csv=None, FD_voxelwise=None, pos_voxelwise=None, commonspace_resampled_template=None, native_bold=None, 
        native_brain_mask=None, native_WM_mask=None, native_CSF_mask=None, native_vascular_mask=None, native_labels=None,
        anat_preproc=None, commonspace_to_native_transform_list=None, commonspace_to_native_inverse_list=None):
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
            ("motion_params_csv", "inputnode.motion_params_csv"),  # confounds file
            ("FD_csv", "inputnode.FD_file"),
            ("input_bold", "inputnode.raw_input_file"),
            ]),
        ])

    if cr_opts.nativespace_analysis:
        workflow.connect([
            (preproc_outputnode, confound_correction_wf, [
                ("native_bold", "inputnode.bold_file"),
                ("native_brain_mask", "inputnode.brain_mask"),
                ("native_WM_mask", "inputnode.WM_mask"),
                ("native_CSF_mask", "inputnode.CSF_mask"),
                ("native_vascular_mask", "inputnode.vascular_mask"),
                ]),
            ])
    else:
        workflow.connect([
            (preproc_outputnode, confound_correction_wf, [
                ("commonspace_bold", "inputnode.bold_file"),
                ("commonspace_mask", "inputnode.brain_mask"),
                ("commonspace_WM_mask", "inputnode.WM_mask"),
                ("commonspace_CSF_mask", "inputnode.CSF_mask"),
                ("commonspace_vascular_mask", "inputnode.vascular_mask"),
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
    if cr_opts.ica_aroma['apply']:
        workflow.connect([
            (confound_correction_wf, confound_correction_datasink, [
                ("outputnode.aroma_out", "aroma_out"),
                ]),
            ])
    if cr_opts.frame_censoring['DVARS_censoring'] or cr_opts.frame_censoring['FD_censoring']:
        workflow.connect([
            (confound_correction_wf, confound_correction_datasink, [
                ("outputnode.frame_mask_file", "frame_censoring_mask"),
                ]),
            ])
        
    if cr_opts.generate_CR_null:
        plot_CR_overfit_node = pe.Node(Function(input_names=['mask_file', 'STD_file_path', 'CR_STD_file_path', 'random_CR_STD_file_path', 'corrected_CR_STD_file_path', 'figure_format'],
                                            output_names=['figure_path'],
                                        function=plot_CR_overfit),
                                name='plot_CR_overfit_node')
        plot_CR_overfit_node.inputs.figure_format = cr_opts.figure_format

        workflow.connect([
            (confound_correction_wf, plot_CR_overfit_node, [
                ("outputnode.STD_file", "STD_file_path"),
                ("outputnode.CR_STD_file", "CR_STD_file_path"),
                ("outputnode.random_CR_STD_file_path", "random_CR_STD_file_path"),
                ("outputnode.corrected_CR_STD_file_path", "corrected_CR_STD_file_path"),
                ]),
            (plot_CR_overfit_node, confound_correction_datasink, [
                ("figure_path", "plot_CR_overfit"),
                ]),
            ])
        
        if cr_opts.nativespace_analysis:
            workflow.connect([
                (preproc_outputnode, plot_CR_overfit_node, [
                    ("native_brain_mask", "mask_file"),
                    ]),
                ])
        else:
            workflow.connect([
                (preproc_outputnode, plot_CR_overfit_node, [
                    ("commonspace_mask", "mask_file"),
                    ]),
                ])



    return workflow



def read_preproc_datasinks(preproc_output, nativespace=False, fast_commonspace=False, atlas_reg_script='SyN', voxelwise_motion=False):
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
        ['motion_datasink','motion_params_csv'], ['motion_datasink','FD_csv']]
    
    if voxelwise_motion:
        directory_list+=[['motion_datasink','FD_voxelwise'], ['motion_datasink','pos_voxelwise']]

    if nativespace:
        directory_list+=[['bold_datasink','native_bold'], ['bold_datasink','native_brain_mask'],
            ['bold_datasink','native_WM_mask'], ['bold_datasink','native_CSF_mask'], ['bold_datasink','native_vascular_mask'], ['bold_datasink','native_labels']]

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

    if nativespace:
        ###
        # For the anat_preproc and transforms, there needs to be a different file matching, where files may be named based on the anat
        # scan, so here we match the BIDS specs. The transforms to native space are put together into a prepared list of transforms.
        ###
        directory_list=[['anat_datasink','anat_preproc']]
        if fast_commonspace:
            directory_list+=[['transforms_datasink','native_to_atlas_affine']]
            if atlas_reg_script=='SyN':
                directory_list+=[['transforms_datasink','native_to_atlas_inverse_warp']]
        else:
            directory_list+=[['transforms_datasink','unbiased_to_atlas_affine'], 
                ['transforms_datasink','native_to_unbiased_affine'], ['transforms_datasink','native_to_unbiased_inverse_warp']]
            if atlas_reg_script=='SyN':
                directory_list+=[['transforms_datasink','unbiased_to_atlas_inverse_warp']]

        from bids.layout import parse_file_entities
        for datasink,target in directory_list:

            if not os.path.isdir(f'{preproc_output}/{datasink}/{target}'):
                raise ValueError(f"The directory {preproc_output}/{datasink}/{target} does not exist. Make sure that all required "
                    "datasink outputs are available. If --bold_only was selected, there are no native space outputs available.")
            target_list.append(target)
            file_list = get_files_from_tree(f'{preproc_output}/{datasink}/{target}')
            for split in split_name:
                # for the unbiased to atlas transforms it is not subject-specific
                if target in ['unbiased_to_atlas_affine', 'unbiased_to_atlas_inverse_warp']:
                    if not len(file_list)==1:
                        raise ValueError(f"There should be only a single transform from unbiased to atlas space. Instead there is {file_list}.")
                    split_dict[split][target]=file_list[0]
                else:
                    parser_split = parse_file_entities('/'+split)
                    for f in file_list:
                        parser_f = parse_file_entities(f)
                        # subject info is mandatory
                        if not 'subject' in list(parser_f.keys()):
                            raise ValueError(f"The file {f} is missing the 'subject' BIDS specs.")
                        if not 'subject' in list(parser_split.keys()):
                            raise ValueError(f"The file {split} is missing the 'subject' BIDS specs.")
                        if parser_f['subject']==parser_split['subject']:
                            # if session was specified in the BOLD file, it should also be in the target file
                            if 'session' in list(parser_split.keys()):
                                if not 'session' in list(parser_f.keys()):
                                    raise ValueError(f"The file {f} is missing the 'session' BIDS specs.")
                                if parser_f['session']==parser_split['session']:
                                    split_dict[split][target]=f
                                    break
                            else:
                                split_dict[split][target]=f
                                break

        for split in split_name:
            if fast_commonspace:
                to_atlas_affine = split_dict[split]['native_to_atlas_affine']
                if atlas_reg_script=='SyN':
                    to_atlas_inverse_warp = split_dict[split]['native_to_atlas_inverse_warp']
                    commonspace_to_native_transform_list=[to_atlas_affine,to_atlas_inverse_warp]
                    commonspace_to_native_inverse_list=[1,0]
                else:
                    commonspace_to_native_transform_list=[to_atlas_affine]
                    commonspace_to_native_inverse_list=[1]
            else:
                native_to_unbiased_inverse_warp = split_dict[split]['native_to_unbiased_inverse_warp']
                native_to_unbiased_affine = split_dict[split]['native_to_unbiased_affine']
                to_atlas_affine = split_dict[split]['unbiased_to_atlas_affine']
                if atlas_reg_script=='SyN':
                    to_atlas_inverse_warp = split_dict[split]['unbiased_to_atlas_inverse_warp']
                    commonspace_to_native_transform_list=[native_to_unbiased_affine,native_to_unbiased_inverse_warp,to_atlas_affine,to_atlas_inverse_warp]
                    commonspace_to_native_inverse_list=[1,0,1,0]
                else:
                    commonspace_to_native_transform_list=[native_to_unbiased_affine,native_to_unbiased_inverse_warp,to_atlas_affine]
                    commonspace_to_native_inverse_list=[1,0,1]

            split_dict[split]['commonspace_to_native_transform_list'] = commonspace_to_native_transform_list
            split_dict[split]['commonspace_to_native_inverse_list'] = commonspace_to_native_inverse_list
        target_list += ['commonspace_to_native_transform_list', 'commonspace_to_native_inverse_list']

        if fast_commonspace:
            if atlas_reg_script=='SyN':
                target_list.remove('native_to_atlas_inverse_warp')
            target_list.remove('native_to_atlas_affine')
        else:
            target_list.remove('native_to_unbiased_inverse_warp')
            target_list.remove('native_to_unbiased_affine')
            if atlas_reg_script=='SyN':
                target_list.remove('unbiased_to_atlas_inverse_warp')
            target_list.remove('unbiased_to_atlas_affine')

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
                    'motion_params_csv':['main_wf.bold_main_wf.estimate_motion_node', 'motion_params_csv'],
                    'FD_voxelwise':['main_wf.bold_main_wf.estimate_motion_node', 'FD_voxelwise'],
                    'pos_voxelwise':['main_wf.bold_main_wf.estimate_motion_node', 'pos_voxelwise'],
                    'FD_csv':['main_wf.bold_main_wf.estimate_motion_node', 'FD_csv'],
                    'commonspace_resampled_template':['main_wf.resample_template', 'commonspace_template'],
                    }
    if nativespace:
        match_targets.update({'native_bold':['main_wf.bold_main_wf.bold_native_trans_wf.merge', 'out_file'],
                        'native_brain_mask':['main_wf.bold_main_wf.bold_native_trans_wf.Brain_mask_EPI', 'EPI_mask'],
                        'native_WM_mask':['main_wf.bold_main_wf.bold_native_trans_wf.WM_mask_EPI', 'EPI_mask'],
                        'native_CSF_mask':['main_wf.bold_main_wf.bold_native_trans_wf.CSF_mask_EPI', 'EPI_mask'],
                        'native_vascular_mask':['main_wf.bold_main_wf.bold_native_trans_wf.vascular_mask_EPI', 'EPI_mask'],
                        'native_labels':['main_wf.bold_main_wf.bold_native_trans_wf.prop_labels_EPI', 'EPI_mask'],
                        'anat_preproc':['main_wf.anat_inho_cor_wf.InhoCorrection', 'corrected'],
                        'commonspace_to_native_transform_list':['main_wf.commonspace_reg_wf.prep_commonspace_transform', 'commonspace_to_native_transform_list'],
                        'commonspace_to_native_inverse_list':['main_wf.commonspace_reg_wf.prep_commonspace_transform', 'commonspace_to_native_inverse_list'],
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


def plot_CR_overfit(mask_file, STD_file_path, CR_STD_file_path, random_CR_STD_file_path, corrected_CR_STD_file_path, figure_format):

    for file in STD_file_path, CR_STD_file_path, random_CR_STD_file_path, corrected_CR_STD_file_path:
        if 'empty' in file:
            return file

    import os
    import matplotlib.pyplot as plt
    import SimpleITK as sitk
    nrows = 4
    fig, axes = plt.subplots(nrows=nrows, ncols=3, figsize=(12*3, 2*nrows))
    plt.tight_layout()

    from rabies.visualization import plot_3d
    volume_indices = sitk.GetArrayFromImage(sitk.ReadImage(mask_file, sitk.sitkFloat32)).astype(bool)

    def get_vmax(sitk_img, volume_indices):
        # select vmax at 95th percentile value
        vector = sitk.GetArrayFromImage(sitk_img)[volume_indices].flatten()
        vector.sort()
        vmax = vector[int(len(vector)*0.95)]
        return vmax
    
    title_list = ['$\mathregular{BOLD_{SD}}$', '$\mathregular{CR_{SD}}$', 
                  '$\mathregular{CR_{SD}}$ random', '$\mathregular{CR_{SD}}$ corrected']
    row=0
    for file, title in zip([STD_file_path, CR_STD_file_path, random_CR_STD_file_path, corrected_CR_STD_file_path], title_list):
        sitk_img = sitk.ReadImage(file)
        cbar_list = plot_3d(axes[row, :], sitk_img, fig, vmin=0, vmax=get_vmax(sitk_img, volume_indices),
                    cmap='inferno', alpha=1, cbar=True, num_slices=6)
        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 35
            cbar.set_label('Standard \n Deviation', fontsize=17, rotation=270, color='white')
            cbar.ax.tick_params(labelsize=15)
        for ax in axes[row, :]:
            ax.set_title(title, fontsize=30, color='white')
        row +=1

    import pathlib
    filename_template = pathlib.Path(STD_file_path).name.rsplit("_STD_map.nii.gz")[0]
    figure_path = os.path.abspath(filename_template)+f'_CR_overfit.{figure_format}'
    fig.savefig(figure_path, bbox_inches='tight')

    return figure_path

