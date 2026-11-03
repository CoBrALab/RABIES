import os
import pathlib
import numpy as np
import pandas as pd
import SimpleITK as sitk
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

from rabies.utils import ResampleVolumes
from .analysis_functions import run_group_ICA


def init_analysis_wf(analysis_opts, cr_opts, name="analysis_wf"):

    workflow = pe.Workflow(name=name)
    subject_inputnode = pe.Node(niu.IdentityInterface(
        fields=['cleaned_bold_file', 'name_source', 'CR_data_dict', 
                'mask_file', 'anat_ref_file', 'to_analysis_space_transform_list', 'to_analysis_space_inverse_list',
                'token', 'native_to_commonspace_transform_list', 'native_to_commonspace_inverse_list']), name='subject_inputnode')
    group_inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file_list', 'commonspace_mask', 'commonspace_template', 'token']), name='group_inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['group_ICA_dir', 'IC_file', 'dual_regression_timecourse_csv',
                                                       'DR_nii_file', 'DR_nii_file_resampled', 'matrix_data_file', 'matrix_fig', 'corr_map_file_list', 'corr_map_file_resampled', 'seed_timecourse_csv_list',
                                                       'sub_token', 'group_token']), name='outputnode')

    # connect the nodes so that they exist even without running analysis
    workflow.connect([
        (subject_inputnode, outputnode, [
            ("token", "sub_token"),
            ]),
        (group_inputnode, outputnode, [
            ("token", "group_token"),
            ]),
        ])

    run_seed_FC = len(analysis_opts.seed_list) > 0

    # a single node loads every input this scan's FC analyses need exactly once, and
    # computes whichever of seed-based FC/dual regression/FC matrix were requested
    if run_seed_FC or analysis_opts.DR_ICA or analysis_opts.FC_matrix:

        FC_analysis_node = pe.Node(FCAnalysis(
            analysis_opts=analysis_opts, cr_opts=cr_opts,
            ), name='FC_analysis', mem_gb=1*analysis_opts.scale_min_memory)

        workflow.connect([
            (subject_inputnode, FC_analysis_node, [
                ("cleaned_bold_file", "cleaned_bold_file"),
                ("name_source", "name_source"),
                ("CR_data_dict", "CR_data_dict"),
                ("mask_file", "mask_file"),
                ("anat_ref_file", "anat_ref_file"),
                ("to_analysis_space_transform_list", "to_analysis_space_transform_list"),
                ("to_analysis_space_inverse_list", "to_analysis_space_inverse_list"),
                ]),
            ])

        if run_seed_FC:
            workflow.connect([
                (FC_analysis_node, outputnode, [
                    ("corr_map_file_list", "corr_map_file_list"),
                    ("seed_timecourse_csv_list", "seed_timecourse_csv_list"),
                    ]),
                ])

            if analysis_opts.resample_to_commonspace:
                SBC_transform_node = pe.MapNode(ResampleVolumes(
                    resampling_dim='ref_file', interpolation=analysis_opts.interpolation_sitk,
                    rabies_data_type=analysis_opts.data_type, apply_motcorr=False, clip_negative=False),
                    iterfield=['in_file', 'name_source'],
                    name='SBC_to_commonspace')

                workflow.connect([
                    (subject_inputnode, SBC_transform_node, [
                        ("native_to_commonspace_transform_list", "transforms_3d_files"),
                        ("native_to_commonspace_inverse_list", "inverses_3d"),
                        ]),
                    (group_inputnode, SBC_transform_node, [
                        ("commonspace_template", "ref_file"),
                        ]),
                    (FC_analysis_node, SBC_transform_node, [
                        ("corr_map_file_list", "in_file"),
                        ("corr_map_file_list", "name_source"),
                        ]),
                    (SBC_transform_node, outputnode, [
                        ("resampled_file", "corr_map_file_resampled"),
                        ]),
                    ])

        if analysis_opts.DR_ICA:
            workflow.connect([
                (FC_analysis_node, outputnode, [
                    ("dual_regression_timecourse_csv", "dual_regression_timecourse_csv"),
                    ("DR_maps_filename", "DR_nii_file"),
                    ]),
                ])

            if analysis_opts.resample_to_commonspace:
                DR_transform_node = pe.Node(ResampleVolumes(
                    resampling_dim='ref_file', interpolation=analysis_opts.interpolation_sitk,
                    rabies_data_type=analysis_opts.data_type, apply_motcorr=False, clip_negative=False),
                    n_procs=int(os.environ['RABIES_ITK_NUM_THREADS']),
                    name='DR_to_commonspace')

                workflow.connect([
                    (subject_inputnode, DR_transform_node, [
                        ("native_to_commonspace_transform_list", "transforms_3d_files"),
                        ("native_to_commonspace_inverse_list", "inverses_3d"),
                        ]),
                    (group_inputnode, DR_transform_node, [
                        ("commonspace_template", "ref_file"),
                        ]),
                    (FC_analysis_node, DR_transform_node, [
                        ("DR_maps_filename", "in_file"),
                        ("DR_maps_filename", "name_source"),
                        ]),
                    (DR_transform_node, outputnode, [
                        ("resampled_file", "DR_nii_file_resampled"),
                        ]),
                    ])

        if analysis_opts.FC_matrix:
            workflow.connect([
                (FC_analysis_node, outputnode, [
                    ("matrix_data_file", "matrix_data_file"),
                    ("matrix_fig", "matrix_fig"),
                    ]),
                ])

    if analysis_opts.group_ica['apply']:
        group_ICA = pe.Node(Function(input_names=['bold_file_list', 'mask_file', 'dim', 'random_seed', 'background_image', 'disableMigp'],
                                     output_names=['out_dir', 'IC_file'],
                                     function=run_group_ICA),
                            name='group_ICA', mem_gb=1*analysis_opts.scale_min_memory)
        group_ICA.inputs.dim = analysis_opts.group_ica['dim']
        group_ICA.inputs.random_seed = analysis_opts.group_ica['random_seed']
        group_ICA.inputs.disableMigp = analysis_opts.group_ica['disableMigp']

        workflow.connect([
            (group_inputnode, group_ICA, [
                ("bold_file_list", "bold_file_list"),
                ("commonspace_mask", "mask_file"),
                ("commonspace_template", "background_image"),
                ]),
            (group_ICA, outputnode, [
                ("IC_file", "IC_file"),
                ("out_dir", "group_ICA_dir"),
                ]),
            ])

    return workflow


class FCAnalysisInputSpec(BaseInterfaceInputSpec):
    # confound correction outputs for this scan
    cleaned_bold_file = File(exists=True, mandatory=True,
        desc="Cleaned/confound-corrected 4D EPI timeseries for this scan, in the analysis space (native or commonspace).")
    name_source = File(exists=True, mandatory=True,
        desc="The original raw EPI file, used only to name outputs consistently with the rest of RABIES.")
    CR_data_dict = traits.Dict(
        desc="Confound correction info dict for this scan; only 'voxelwise_mean' is read from it, to undo the intercept re-added by --keep_EPI_average.")

    # definition of the analysis space (native or commonspace, resolved by the caller)
    mask_file = File(exists=True, mandatory=True, desc="Brain mask in the analysis space.")
    anat_ref_file = File(exists=True, mandatory=True,
        desc="Anatomical reference defining the analysis space, onto which priors/seeds/atlas are resampled.")
    to_analysis_space_transform_list = traits.List([], usedefault=True,
        desc="Transforms mapping priors/seeds/atlas (in their original commonspace) onto the analysis space. Empty if already in commonspace.")
    to_analysis_space_inverse_list = traits.List([], usedefault=True,
        desc="Booleans indicating which of to_analysis_space_transform_list must be applied as an inverse.")
    analysis_opts = traits.Any(
        exists=True, mandatory=True, desc="RABIES analysis parameters specs.")
    cr_opts = traits.Any(
        exists=True, mandatory=True, desc="RABIES confound correction parameters specs.")


class FCAnalysisOutputSpec(TraitedSpec):
    corr_map_file_list = traits.List(desc="Seed-based FC map for each seed in seed_dict.")
    seed_timecourse_csv_list = traits.List(desc="Seed timecourse .csv for each seed in seed_dict.")
    DR_maps_filename = traits.Any(desc=".nii file with the dual regression spatial maps.")
    dual_regression_timecourse_csv = traits.Any(desc=".csv with the dual regression timecourses.")
    matrix_data_file = traits.Any(desc=".csv with the FC matrix.")
    matrix_fig = traits.Any(desc="Figure of the FC matrix.")


class FCAnalysis(BaseInterface):
    """
    Central per-scan FC analysis. Reads the raw RABIES outputs for one scan once,
    then computes every requested scan-level analysis -- seed-based FC, dual 
    regression, FC matrix.
    """

    input_spec = FCAnalysisInputSpec
    output_spec = FCAnalysisOutputSpec

    def _run_interface(self, runtime):
        from rabies.utils import recover_3D, recover_4D

        analysis_opts = self.inputs.analysis_opts
        cr_opts = self.inputs.cr_opts

        seed_dict = {seed_name: seed_file for seed_file, seed_name in zip(analysis_opts.seed_list, analysis_opts.seed_name_list)}

        SBC_out, DR_out, FC_matrix_df = fc_analysis(
            self.inputs.cleaned_bold_file, self.inputs.name_source, self.inputs.mask_file, self.inputs.anat_ref_file, # minimal inputs
            seed_dict=seed_dict, # SBC options
            FC_matrix=analysis_opts.FC_matrix, atlas_file=analysis_opts.ROI_labels_file, ROI_type=analysis_opts.ROI_type, # FC matrix options
            prior_maps=analysis_opts.prior_maps, DR_ICA=analysis_opts.DR_ICA, network_weighting=analysis_opts.network_weighting, # prior/dual regression options
            CR_data_dict=self.inputs.CR_data_dict, remove_EPI_avg=cr_opts.keep_EPI_average, # whether to remove the EPI average from the cleaned timeseries before analysis
            to_analysis_space_transform_list=self.inputs.to_analysis_space_transform_list, to_analysis_space_inverse_list=self.inputs.to_analysis_space_inverse_list, interpolation=analysis_opts.interpolation_sitk, # resampling to a different space
            rabies_data_type=analysis_opts.data_type,
            )


        '''
        Prepare output files
        '''
        filename_split = pathlib.Path(self.inputs.name_source).name.rsplit(".nii")[0]

        corr_map_file_l = []
        seed_timecourse_csv_l = []
        if SBC_out is not None:
            for seed_name, (corr_map, seed_timeseries) in SBC_out.items():
                corr_map_img = recover_3D(self.inputs.mask_file, corr_map)
                corr_map_file = os.path.abspath(f'{filename_split}_{seed_name}_corr_map.nii.gz')
                sitk.WriteImage(corr_map_img, corr_map_file)

                seed_timecourse_csv = os.path.abspath(f'{filename_split}_{seed_name}_timecourse.csv')
                pd.DataFrame(seed_timeseries).to_csv(seed_timecourse_csv, header=False, index=False)

                corr_map_file_l.append(corr_map_file)
                seed_timecourse_csv_l.append(seed_timecourse_csv)

        if DR_out is None:
            DR_maps_filename = None
            dual_regression_timecourse_csv = None
        else:
            DR_C, DR_W = DR_out
            dual_regression_timecourse_csv = os.path.abspath(f'{filename_split}_dual_regression_timecourse.csv')
            pd.DataFrame(DR_W).to_csv(dual_regression_timecourse_csv, header=False, index=False)

            DR_maps_filename = os.path.abspath(f'{filename_split}_DR_maps.nii.gz')
            sitk.WriteImage(recover_4D(self.inputs.mask_file, DR_C.T, self.inputs.cleaned_bold_file), DR_maps_filename)

        if FC_matrix_df is not None:
            from .analysis_functions import plot_matrix
            matrix_fig = os.path.abspath(f'{filename_split}_FC_matrix.{analysis_opts.figure_format}')
            plot_matrix(matrix_fig, FC_matrix_df.values)

            matrix_data_file = os.path.abspath(f'{filename_split}_FC_matrix.csv')
            FC_matrix_df.to_csv(matrix_data_file, sep=',')
        else:
            matrix_data_file = None
            matrix_fig = None


        setattr(self, 'corr_map_file_list', corr_map_file_l)
        setattr(self, 'seed_timecourse_csv_list', seed_timecourse_csv_l)
        setattr(self, 'DR_maps_filename', DR_maps_filename)
        setattr(self, 'dual_regression_timecourse_csv', dual_regression_timecourse_csv)
        setattr(self, 'matrix_data_file', matrix_data_file)
        setattr(self, 'matrix_fig', matrix_fig)

        return runtime


    def _list_outputs(self):
        return {
            'corr_map_file_list': getattr(self, 'corr_map_file_list', []),
            'seed_timecourse_csv_list': getattr(self, 'seed_timecourse_csv_list', []),
            'DR_maps_filename': getattr(self, 'DR_maps_filename', None),
            'dual_regression_timecourse_csv': getattr(self, 'dual_regression_timecourse_csv', None),
            'matrix_data_file': getattr(self, 'matrix_data_file', None),
            'matrix_fig': getattr(self, 'matrix_fig', None),
        }



def fc_analysis(
    cleaned_bold_file, name_source, mask_file, anat_ref_file, # minimal inputs
    seed_dict={}, # SBC options
    FC_matrix=False, atlas_file=None, ROI_type='parcellated', # FC matrix options
    prior_maps=None, DR_ICA=False, network_weighting='absolute', # prior/dual regression options
    CR_data_dict={}, remove_EPI_avg=False, # whether to remove the EPI average from the cleaned timeseries before analysis
    to_analysis_space_transform_list=[], to_analysis_space_inverse_list=[], interpolation=sitk.sitkLinear, # resampling to a different space
    rabies_data_type=sitk.sitkFloat32,
    ):
    """
    Core function for computing scan-level FC analyses -- seed-based connectivity (SBC),
    dual regression (DR), and/or an FC matrix -- from a single
    load of a scan's cleaned timeseries and reference maps. Each analysis is only run if
    its trigger parameters are set; otherwise its corresponding output is None.

    Parameters
    ----------

    cleaned_bold_file : filepath
        Cleaned/confound-corrected 4D EPI timeseries for this scan, in the analysis space
        (native or commonspace).

    name_source : filepath
        The original raw EPI file, used only to name intermediate outputs consistently
        with the rest of RABIES.

    mask_file : filepath
        Brain mask in the analysis space.

    anat_ref_file : filepath
        Anatomical reference defining the analysis space, onto which seed_dict/prior_maps/
        atlas_file are resampled.

    seed_dict : dict, default={}
        Dictionary of seed_name:seed_file pairs to compute SBC for. If empty, SBC is
        skipped.

    FC_matrix : bool, default=False
        Whether to compute a whole-brain FC matrix.

    atlas_file : filepath, default=None
        Parcellation file, required for a parcellated (ROI_type='parcellated') FC matrix.

    ROI_type : str, default='parcellated'
        Either 'parcellated' or 'voxelwise' FC matrix.

    prior_maps : filepath, default=None
        4D file of network priors, required for DR_ICA.

    DR_ICA : bool, default=False
        Whether to compute dual regression against prior_maps.

    network_weighting : str, default='absolute'
        Whether DR network maps are 'absolute' or 'relative' (variance-normalized).

    CR_data_dict : dict, default={}
        Confound correction info dict for this scan; only 'voxelwise_mean' is read from
        it, and only if remove_EPI_avg is True.

    remove_EPI_avg : bool, default=False
        Whether the voxelwise average was kept in the cleaned timeseries
        (--keep_EPI_average) and must be subtracted back out, using
        CR_data_dict['voxelwise_mean'], prior to analysis.

    to_analysis_space_transform_list : list, default=[]
        Transforms mapping seed_dict/prior_maps/atlas_file, in their original commonspace,
        onto the analysis space. Empty if already in commonspace.

    to_analysis_space_inverse_list : list, default=[]
        Booleans indicating which of to_analysis_space_transform_list must be applied as an inverse.

    interpolation : SimpleITK interpolator, default=sitk.sitkLinear
        Interpolator used to resample prior_maps.

    rabies_data_type : SimpleITK pixel type, default=sitk.sitkFloat32
        Data type used for resampling.

    Returns
    -------

    SBC_out : dict or None
        {seed_name: (corr_map, seed_timeseries)} for each seed in seed_dict, or None if
        seed_dict is empty.

    DR_out : tuple or None
        (DR_C, DR_W), the dual regression spatial maps and timecourses, or None if DR_ICA
        is False.

    FC_matrix_df : pd.DataFrame or None
        The FC matrix, or None if FC_matrix is False.
    """
    from .analysis_math import dual_regression
    from .analysis_functions import compute_seed_FC, parcellated_FC_matrix
    from .utils import load_resample_analysis_maps

    filename_split = pathlib.Path(name_source).name.rsplit(".nii")[0]

    # prior_maps and atlas_file are replaced by None to skip the resampling step inside 
    # load_resample_analysis_maps() if they are not used by their respective analysis functions
    prior_maps = prior_maps if DR_ICA else None
    atlas_file = atlas_file if FC_matrix else None

    loaded = load_resample_analysis_maps(
        mask_file, anat_ref_file,
        transform_list=to_analysis_space_transform_list, inverse_list=to_analysis_space_inverse_list,
        seed_dict=seed_dict, prior_maps=prior_maps, atlas_file=atlas_file,
        cleaned_bold_file=cleaned_bold_file, CR_data_dict=CR_data_dict, remove_EPI_avg=remove_EPI_avg,
        interpolation=interpolation, rabies_data_type=rabies_data_type, output_prefix=filename_split,
        )
    timeseries = loaded['timeseries']
    seed_arr_dict = loaded['seed_arr_dict']
    prior_map_vectors = loaded['prior_map_vectors']
    atlas_idx = loaded['atlas_idx']
    roi_list = loaded['roi_list']

    '''
    FC computations
    '''
    if len(seed_arr_dict.keys())>0:
        SBC_out = compute_seed_FC(timeseries, seed_arr_dict)
    else:
        SBC_out = None

    if DR_ICA:
        if prior_map_vectors is None:
            raise ValueError("DR_ICA=True requires prior_maps to be provided.")
        DR = dual_regression(prior_map_vectors, timeseries)
        if network_weighting == 'absolute':
            DR_C = DR['C'] * DR['S']
        elif network_weighting == 'relative':
            DR_C = DR['C']
        else:
            raise ValueError(f"Invalid network_weighting: {network_weighting}")
        DR_out = (DR_C, DR['W'])
    else:
        DR_out = None


    if not FC_matrix:
        FC_matrix_df = None
    else:
        if ROI_type == 'parcellated':
            if atlas_idx is None or roi_list is None:
                raise ValueError("ROI_type='parcellated' requires atlas_file to be provided.")
            corr_matrix, roi_labels = parcellated_FC_matrix(timeseries, atlas_idx, roi_list)
            FC_matrix_df = pd.DataFrame(corr_matrix, index=roi_labels, columns=roi_labels)
        elif ROI_type == 'voxelwise':
            corr_matrix = np.corrcoef(timeseries.T)
            FC_matrix_df = pd.DataFrame(corr_matrix)
        else:
            raise ValueError(f"Invalid --ROI_type provided: {ROI_type}. Must be either 'parcellated' or 'voxelwise.'")

    return SBC_out, DR_out, FC_matrix_df

