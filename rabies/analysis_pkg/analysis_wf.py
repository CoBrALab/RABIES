import os
import pathlib
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function

from .analysis_functions import run_group_ICA, run_DR_ICA, run_FC_matrix, seed_based_FC


def init_analysis_wf(opts, commonspace_cr=False, seed_list=[], name="analysis_wf"):

    workflow = pe.Workflow(name=name)
    subject_inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file', 'mask_file', 'atlas_file', 'token']), name='subject_inputnode')
    group_inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file_list', 'commonspace_mask', 'token']), name='group_inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['group_ICA_dir', 'IC_file', 'DR_data_file',
                                                       'DR_nii_file', 'matrix_data_file', 'matrix_fig', 'corr_map_file', 'sub_token', 'group_token']), name='outputnode')

    # connect the nodes so that they exist even without running analysis
    workflow.connect([
        (subject_inputnode, outputnode, [
            ("token", "sub_token"),
            ]),
        (group_inputnode, outputnode, [
            ("token", "group_token"),
            ]),
        ])

    if len(seed_list) > 0:
        if not commonspace_cr:
            raise ValueError(
                'Outputs from confound regression must be in commonspace to run seed-based analysis. Try running confound regression again with --commonspace_bold.')
        seed_based_FC_node = pe.Node(Function(input_names=['bold_file', 'brain_mask', 'seed_dict', 'seed_name'],
                                              output_names=['corr_map_file'],
                                              function=seed_based_FC),
                                     name='seed_based_FC', mem_gb=1)
        seed_dict = {}
        name_list = []
        for file in seed_list:
            file = os.path.abspath(file)
            if not os.path.isfile(file):
                raise ValueError(
                    "Provide seed file path %s doesn't exists." % (file))
            seed_name = pathlib.Path(file).name.rsplit(".nii")[0]
            name_list.append(seed_name)
            seed_dict[seed_name] = file
        seed_based_FC_node.iterables = ('seed_name', name_list)
        seed_based_FC_node.inputs.seed_dict = seed_dict

        workflow.connect([
            (subject_inputnode, seed_based_FC_node, [
                ("bold_file", "bold_file"),
                ("mask_file", "brain_mask"),
                ]),
            (seed_based_FC_node, outputnode, [
                ("corr_map_file", "corr_map_file"),
                ]),
            ])

    include_group_ICA = opts.group_ICA

    if opts.DR_ICA:
        if not commonspace_cr:
            raise ValueError(
                'Outputs from confound regression must be in commonspace to run dual regression. Try running confound regression again with --commonspace_bold.')

        DR_ICA = pe.Node(Function(input_names=['bold_file', 'mask_file', 'IC_file'],
                                  output_names=['data_file', 'nii_file'],
                                  function=run_DR_ICA),
                         name='DR_ICA', mem_gb=1)

        workflow.connect([
            (subject_inputnode, DR_ICA, [
                ("bold_file", "bold_file"),
                ("mask_file", "mask_file"),
                ]),
            (DR_ICA, outputnode, [
                ("data_file", "DR_data_file"),
                ("nii_file", "DR_nii_file"),
                ]),
            ])

        if opts.IC_file is None:
            import logging
            log = logging.getLogger('root')
            log.info(
                'Group-ICA will be run on the processed dataset since no previous group-ICA file was provided.')
            include_group_ICA = True
        elif not os.path.isfile(str(opts.IC_file)):
            raise ValueError("--IC_file doesn't exists.")
        else:
            DR_ICA.inputs.IC_file = os.path.abspath(opts.IC_file)

    if include_group_ICA:
        if not commonspace_cr:
            raise ValueError(
                'Outputs from confound regression must be in commonspace to run group-ICA. Try running confound regression again with --commonspace_bold.')
        group_ICA = pe.Node(Function(input_names=['bold_file_list', 'mask_file', 'dim', 'tr'],
                                     output_names=['out_dir', 'IC_file'],
                                     function=run_group_ICA),
                            name='group_ICA', mem_gb=1)
        group_ICA.inputs.tr = float(opts.TR.split('s')[0])
        group_ICA.inputs.dim = opts.dim

        workflow.connect([
            (group_inputnode, group_ICA, [
                ("bold_file_list", "bold_file_list"),
                ("commonspace_mask", "mask_file"),
                ]),
            (group_ICA, outputnode, [
                ("IC_file", "IC_file"),
                ("out_dir", "group_ICA_dir"),
                ]),
            ])

        if (opts.IC_file is None) and opts.DR_ICA:
            workflow.connect([
                (group_ICA, DR_ICA, [
                    ("IC_file", "IC_file"),
                    ]),
                ])

    if opts.FC_matrix:
        FC_matrix = pe.Node(Function(input_names=['bold_file', 'mask_file', 'atlas', 'roi_type'],
                                     output_names=['data_file', 'figname'],
                                     function=run_FC_matrix),
                            name='FC_matrix', mem_gb=1)
        FC_matrix.inputs.roi_type = opts.ROI_type

        workflow.connect([
            (subject_inputnode, FC_matrix, [
                ("bold_file", "bold_file"),
                ("mask_file", "mask_file"),
                ("atlas_file", "atlas"),
                ]),
            (FC_matrix, outputnode, [
                ("data_file", "matrix_data_file"),
                ("figname", "matrix_fig"),
                ]),
            ])

    return workflow
