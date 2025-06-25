from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from .hmc import init_bold_hmc_wf,EstimateMotionParams
from .bold_ref import init_bold_reference_wf
from .resampling import init_bold_preproc_trans_wf
from .stc import init_bold_stc_wf
from .inho_correction import init_inho_correction_wf
from .registration import init_cross_modal_reg_wf
from nipype.interfaces.utility import Function
from .utils import apply_despike

def init_bold_main_wf(opts, output_folder, number_functional_scans, inho_cor_only=False, name='bold_main_wf'):
    """
    This workflow controls the functional preprocessing stages of the pipeline when both
    functional and anatomical images are provided.

    **Parameters**

        opts
            parser options for preprocess
        inho_cor_only
            whether to run the bias correction steps, or further processing steps.

    **Inputs**

        bold
            Input BOLD series NIfTI file
        coreg_anat
            Anatomical reference for BOLD alignment
        coreg_mask
            Brain mask for anatomical reference
        WM_mask
            WM mask inherited from the common space registration
        CSF_mask
            CSF mask inherited from the common space registration
        vascular_mask
            vascular mask inherited from the common space registration
        labels
            Anatomical labels inherited from the common space registration
        unbiased_to_atlas_affine
            affine transform from the dataset template space to the commonspace space
        unbiased_to_atlas_warp
            non-linear transform from the dataset template space to the commonspace space
        native_to_unbiased_affine
            affine transform from the subject anatomical to the dataset template space
        native_to_unbiased_warp
            non-linear transform from the subject anatomical to the dataset template space
        commonspace_ref
            commonspace anatomical template

    **Outputs**

        input_bold
            The provided input BOLD file
        bold_ref
            Initial EPI median volume subsequently used as 3D reference EPI volume
        motcorr_params
            motion parameters file provided from antsMotionCorr
        init_denoise
            Corrected 3D ref EPI after initial correction step
        denoise_mask
            resampled mask used for final denoising
        corrected_EPI
            3D reference EPI volume after bias field correction
        output_warped_bold
            Bias field corrected 3D EPI volume warped to the anatomical space
        bold_to_anat_affine
            affine transform from the EPI space to the anatomical space
        bold_to_anat_warp
            non-linear transform from the EPI space to the anatomical space
        bold_to_anat_inverse_warp
            inverse non-linear transform from the EPI space to the anatomical space
        resampled_bold
            Original BOLD timeseries resampled through motion realignment and
            susceptibility distortion correction based on registration to the
            anatomical image
        resampled_ref_bold
            3D median EPI volume from the resampled native BOLD timeseries
        motion_params_csv
            .csv file with measured motion timecourses, used as regressors for confound
            correction: 6 rigid body motion parameters + their first temporal derivate 
            + the 12 parameters squared (24 motion parameters)
        FD_voxelwise
            Voxelwise framewise displacement (FD) measures that can be integrated
            to future confound regression.
            These measures are computed from antsMotionCorrStats.
        pos_voxelwise
            Voxel distancing across time based on rigid body movement parameters,
            which can be integrated for a voxelwise motion regression
            These measures are computed from antsMotionCorrStats.
        FD_csv
            .csv file with global framewise displacement (FD) measures
        EPI_brain_mask
            EPI brain mask for resampled bold
        EPI_WM_mask
            EPI WM mask for resampled bold
        EPI_CSF_mask
            EPI CSF mask for resampled bold
        EPI_labels
            EPI anatomical labels for resampled bold
        commonspace_bold
            Motion and SDC-corrected EPI timeseries resampled into common space
            by applying transforms from the anatomical common space registration
        commonspace_mask
            EPI brain mask for commonspace bold
        commonspace_WM_mask
            EPI WM mask for commonspace bold
        commonspace_CSF_mask
            EPI CSF mask for commonspace bold
        commonspace_vascular_mask
            EPI vascular mask for commonspace bold
        commonspace_labels
            EPI anatomical labels for commonspace bold
    """

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
                fields=['bold', 'inho_cor_anat', 'inho_cor_mask', 'coreg_anat', 'coreg_mask',
                        'native_to_commonspace_transform_list','native_to_commonspace_inverse_list',
                        'commonspace_to_native_transform_list','commonspace_to_native_inverse_list',
                        'commonspace_ref']),
                        name="inputnode")

    outputnode = pe.Node(niu.IdentityInterface(
                fields=['input_bold', 'bold_ref', 'motcorr_params', 'init_denoise', 'denoise_mask', 'corrected_EPI',
                        'output_warped_bold', 'bold_to_anat_affine', 'bold_to_anat_warp', 'bold_to_anat_inverse_warp',
                        'native_bold', 'native_bold_ref', 'native_brain_mask', 'native_WM_mask', 'native_CSF_mask', 'native_vascular_mask', 'native_labels',
                        'motion_params_csv', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv', 'commonspace_bold', 'commonspace_mask',
                        'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_vascular_mask', 'commonspace_labels',
                        'raw_brain_mask']),
                name='outputnode')

    boldbuffer = pe.Node(niu.IdentityInterface(fields=['bold_file']),
                         name="boldbuffer")

    # this node will serve as a relay of outputs from the inho_cor main_wf to the inputs for the rest of the main_wf for bold_only
    transitionnode = pe.Node(niu.IdentityInterface(fields=['bold_file', 'isotropic_bold_file', 'bold_ref', 'init_denoise', 'denoise_mask', 'corrected_EPI']),
                             name="transitionnode")

    if inho_cor_only or (not opts.bold_only):
        template_inputnode = pe.Node(niu.IdentityInterface(fields=['template_anat', 'template_mask']),
                                            name="template_inputnode")


        bold_reference_wf = init_bold_reference_wf(opts=opts)

        num_procs = min(opts.local_threads, number_functional_scans)
        inho_cor_wf = init_inho_correction_wf(opts=opts, image_type='EPI', output_folder=output_folder, num_procs=num_procs, name="bold_inho_cor_wf")

        if opts.isotropic_HMC:
            def resample_isotropic(bold_file, rabies_data_type):
                import numpy as np
                import SimpleITK as sitk
                import os
                import pathlib
                from rabies.utils import resample_image_spacing_4d
                image_4d = sitk.ReadImage(bold_file, rabies_data_type)
                # the image gets resample to the dimension of the axis with highest resolution
                min_dim = np.array(image_4d.GetSpacing()[:3]).min()
                output_spacing = (min_dim,min_dim,min_dim)
                resampled = resample_image_spacing_4d(image_4d, output_spacing, clip_negative=True)
                filename_split = pathlib.Path(
                    bold_file).name.rsplit(".nii")
                isotropic_bold_file = os.path.abspath(filename_split[0]+'_isotropic.nii.gz')
                sitk.WriteImage(resampled, isotropic_bold_file)
                return isotropic_bold_file

            isotropic_resampling_node = pe.Node(Function(input_names=['bold_file', 'rabies_data_type'],
                                                            output_names=[
                                                                'isotropic_bold_file'],
                                                            function=resample_isotropic),
                                                    name='resample_isotropic')
            isotropic_resampling_node.inputs.rabies_data_type = opts.data_type

            workflow.connect([
                (boldbuffer, isotropic_resampling_node, [
                    ('bold_file', 'bold_file'),
                    ]),
                (isotropic_resampling_node, bold_reference_wf, [
                    ('isotropic_bold_file', 'inputnode.bold_file'),
                    ]),
                (isotropic_resampling_node, transitionnode, [
                    ('isotropic_bold_file', 'isotropic_bold_file'),
                    ]),
                ])
        else:
            workflow.connect([
                (boldbuffer, bold_reference_wf, [
                    ('bold_file', 'inputnode.bold_file'),
                    ]),
                ])

        if opts.apply_despiking:
            despike = pe.Node(Function(input_names=['in_file'],
                                                        output_names=['out_file'],
                                                        function=apply_despike),
                                            name='despike')
            workflow.connect([
                (inputnode, despike, [('bold', 'in_file')]),
                (despike, boldbuffer, [('out_file', 'bold_file')]),
                ])
        else:
            workflow.connect([
                (inputnode, boldbuffer, [('bold', 'bold_file')]),
                ])

        if opts.detect_dummy:
            def remove_dummy(bold_file):
                import SimpleITK as sitk
                import pathlib
                import os
                from rabies.utils import copyInfo_4DImage
                from rabies.preprocess_pkg.bold_ref import _get_vols_to_discard

                in_nii = sitk.ReadImage(bold_file)
                data_array = sitk.GetArrayFromImage(in_nii)
                n_volumes_to_discard = _get_vols_to_discard(in_nii)
                if (not n_volumes_to_discard == 0):
                    filename_split = pathlib.Path(bold_file).name.rsplit(".nii")
                    out_bold_file = os.path.abspath(
                        f'{filename_split[0]}_cropped_dummy.nii.gz')
                    img_array = data_array[n_volumes_to_discard:, :, :, :]

                    image_4d = copyInfo_4DImage(sitk.GetImageFromArray(
                        img_array, isVector=False), in_nii, in_nii)
                    sitk.WriteImage(image_4d, out_bold_file)
                else:
                    out_bold_file = bold_file
                return out_bold_file

            remove_dummy_node = pe.Node(Function(input_names=['bold_file'],
                                                            output_names=[
                                                                'out_bold_file'],
                                                            function=remove_dummy),
                                                    name='remove_dummy')
            workflow.connect([
                (boldbuffer, remove_dummy_node, [
                    ('bold_file', 'bold_file'),
                    ]),
                (remove_dummy_node, transitionnode, [
                    ('out_bold_file', 'bold_file'),
                    ]),
                ])
        else:
            workflow.connect([
                (boldbuffer, transitionnode, [
                    ('bold_file', 'bold_file'),
                    ]),
                ])

        workflow.connect([
            (inputnode, inho_cor_wf, [
                ('inho_cor_anat', 'inputnode.anat_ref'),
                ('inho_cor_mask', 'inputnode.anat_mask'),
                ('bold', 'inputnode.name_source'),
                ]),
            (template_inputnode, inho_cor_wf, [
                ("template_anat", "template_inputnode.template_anat"),
                ("template_mask", "template_inputnode.template_mask"),
                ]),
            (bold_reference_wf, inho_cor_wf, [
                ('outputnode.ref_image', 'inputnode.target_img'),
                ]),
            (bold_reference_wf, transitionnode, [
                ('outputnode.ref_image', 'bold_ref'),
                ]),
            (inho_cor_wf, transitionnode, [
                ('outputnode.init_denoise', 'init_denoise'),
                ('outputnode.denoise_mask', 'denoise_mask'),
                ('outputnode.corrected', 'corrected_EPI'),
                ]),
            ])

    if inho_cor_only:
        return workflow

    bold_stc_wf = init_bold_stc_wf(opts=opts)

    # HMC on the BOLD
    bold_hmc_wf = init_bold_hmc_wf(opts=opts)

    bold_commonspace_trans_wf = init_bold_preproc_trans_wf(opts=opts, resampling_dim='ref_file', name='bold_commonspace_trans_wf')
    bold_commonspace_trans_wf.inputs.inputnode.mask_transforms_list = []
    bold_commonspace_trans_wf.inputs.inputnode.mask_inverses = []

    estimate_motion_node = pe.Node(EstimateMotionParams(),
                                name='estimate_motion_node', mem_gb=2.3*opts.scale_min_memory)
    estimate_motion_node.plugin_args = {
        'qsub_args': f'-pe smp {str(2*opts.min_proc)}', 'overwrite': True}


    def prep_resampling_transforms(native_to_commonspace_transform_list,native_to_commonspace_inverse_list, bold_to_anat_warp, bold_to_anat_inverse_warp, bold_to_anat_affine, 
                                   commonspace_to_native_transform_list, commonspace_to_native_inverse_list, bold_only=False):
        if bold_only:
            to_commonspace_transform_list = native_to_commonspace_transform_list
            to_commonspace_inverse_list = native_to_commonspace_inverse_list
            commonspace_to_raw_transform_list = commonspace_to_native_transform_list
            commonspace_to_raw_inverse_list = commonspace_to_native_inverse_list
            raw_to_native_transform_list = None
            raw_to_native_inverse_list = None
        else:
            to_commonspace_transform_list = native_to_commonspace_transform_list+[bold_to_anat_warp, bold_to_anat_affine]
            to_commonspace_inverse_list = native_to_commonspace_inverse_list+[0,0]
            commonspace_to_raw_transform_list = [bold_to_anat_affine, bold_to_anat_inverse_warp]+commonspace_to_native_transform_list
            commonspace_to_raw_inverse_list = [1,0]+commonspace_to_native_inverse_list
            raw_to_native_transform_list = [bold_to_anat_warp, bold_to_anat_affine]
            raw_to_native_inverse_list = [0, 0]

        return to_commonspace_transform_list, to_commonspace_inverse_list, raw_to_native_transform_list, raw_to_native_inverse_list, commonspace_to_raw_transform_list, commonspace_to_raw_inverse_list

    prep_resampling_transforms_node = pe.Node(Function(input_names=['native_to_commonspace_transform_list','native_to_commonspace_inverse_list', 'bold_to_anat_warp', 'bold_to_anat_inverse_warp', 'bold_to_anat_affine', 
                                                                    'commonspace_to_native_transform_list', 'commonspace_to_native_inverse_list', 'bold_only'],
                                                    output_names=[
                                                        'to_commonspace_transform_list','to_commonspace_inverse_list', 'raw_to_native_transform_list', 'raw_to_native_inverse_list',
                                                        'commonspace_to_raw_transform_list', 'commonspace_to_raw_inverse_list'],
                                                    function=prep_resampling_transforms),
                                            name='prep_resampling_transforms')
    prep_resampling_transforms_node.inputs.bold_only = opts.bold_only

    if not opts.bold_only:
        cross_modal_reg_wf = init_cross_modal_reg_wf(opts=opts)

        bold_native_trans_wf = init_bold_preproc_trans_wf(opts=opts, resampling_dim=opts.nativespace_resampling, name='bold_native_trans_wf')

        workflow.connect([
            (inputnode, cross_modal_reg_wf, [
                ('coreg_anat', 'inputnode.anat_ref'),
                ('coreg_mask', 'inputnode.anat_mask')]),
            (inputnode, bold_native_trans_wf, [
                ('commonspace_to_native_transform_list', 'inputnode.mask_transforms_list'),
                ('commonspace_to_native_inverse_list', 'inputnode.mask_inverses'),
                ('bold', 'inputnode.name_source'),
                ]),
            (transitionnode, cross_modal_reg_wf, [
                ('corrected_EPI', 'inputnode.ref_bold_brain'),
                ('denoise_mask', 'inputnode.moving_mask'),
                ]),
            (cross_modal_reg_wf, outputnode, [
                ('outputnode.bold_to_anat_affine', 'bold_to_anat_affine'),
                ('outputnode.bold_to_anat_warp', 'bold_to_anat_warp'),
                ('outputnode.bold_to_anat_inverse_warp', 'bold_to_anat_inverse_warp'),
                ('outputnode.output_warped_bold', 'output_warped_bold'),
                ]),
            (cross_modal_reg_wf, prep_resampling_transforms_node, [
                ('outputnode.bold_to_anat_affine', 'bold_to_anat_affine'),
                ('outputnode.bold_to_anat_warp', 'bold_to_anat_warp'),
                ('outputnode.bold_to_anat_inverse_warp', 'bold_to_anat_inverse_warp'),
                ]),
            (prep_resampling_transforms_node, bold_native_trans_wf, [
                ('raw_to_native_transform_list', 'inputnode.transforms_list'),
                ('raw_to_native_inverse_list', 'inputnode.inverses'),
                ('commonspace_to_raw_transform_list', 'inputnode.commonspace_to_raw_transform_list'),
                ('commonspace_to_raw_inverse_list', 'inputnode.commonspace_to_raw_inverse_list'),
                ]),
            (transitionnode, bold_native_trans_wf, [
                ('bold_ref', 'inputnode.raw_bold_ref'),
                ]),
            (cross_modal_reg_wf, bold_native_trans_wf, [
                ('outputnode.output_warped_bold', 'inputnode.ref_file')]),
            (bold_hmc_wf, bold_native_trans_wf, [
             ('outputnode.motcorr_params', 'inputnode.motcorr_params')]),
            (bold_native_trans_wf, outputnode, [
                ('outputnode.bold', 'native_bold'),
                ('outputnode.bold_ref','native_bold_ref'),
                ('outputnode.brain_mask', 'native_brain_mask'),
                ('outputnode.WM_mask', 'native_WM_mask'),
                ('outputnode.CSF_mask', 'native_CSF_mask'),
                ('outputnode.vascular_mask', 'native_vascular_mask'),
                ('outputnode.labels', 'native_labels'),
                ]),
            ])

    else:
        prep_resampling_transforms_node.inputs.bold_to_anat_warp = None
        prep_resampling_transforms_node.inputs.bold_to_anat_inverse_warp = None
        prep_resampling_transforms_node.inputs.bold_to_anat_affine = None

    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (inputnode, prep_resampling_transforms_node, [
            ('native_to_commonspace_transform_list', 'native_to_commonspace_transform_list'),
            ('native_to_commonspace_inverse_list', 'native_to_commonspace_inverse_list'),
            ('commonspace_to_native_transform_list', 'commonspace_to_native_transform_list'),
            ('commonspace_to_native_inverse_list', 'commonspace_to_native_inverse_list'),
            ]),
        (transitionnode, bold_stc_wf, [
            ('bold_file', 'inputnode.bold_file'),
            ]),
        (transitionnode, bold_hmc_wf, [
            ('bold_ref', 'inputnode.ref_image'),
            ]),
        (bold_hmc_wf, outputnode, [
            ('outputnode.motcorr_params', 'motcorr_params')]),
        (transitionnode, outputnode, [
            ('bold_ref', 'bold_ref'),
            ('init_denoise', 'init_denoise'),
            ('denoise_mask', 'denoise_mask'),
            ('corrected_EPI', 'corrected_EPI'),
            ]),
        (bold_hmc_wf, estimate_motion_node, [
            ('outputnode.motcorr_params', 'motcorr_params'),
            ]),
        (estimate_motion_node, outputnode, [
            ('motion_params_csv', 'motion_params_csv'),
            ('FD_csv', 'FD_csv'),
            ]),
        (prep_resampling_transforms_node, bold_commonspace_trans_wf, [
            ('to_commonspace_transform_list', 'inputnode.transforms_list'),
            ('to_commonspace_inverse_list', 'inputnode.inverses'),
            ('commonspace_to_raw_transform_list', 'inputnode.commonspace_to_raw_transform_list'),
            ('commonspace_to_raw_inverse_list', 'inputnode.commonspace_to_raw_inverse_list'),
            ]),
        (transitionnode, bold_commonspace_trans_wf, [
            ('bold_ref', 'inputnode.raw_bold_ref'),
            ]),
        (bold_commonspace_trans_wf, estimate_motion_node, [
            ('outputnode.raw_brain_mask', 'raw_brain_mask'),
            ]),
        (bold_hmc_wf, bold_commonspace_trans_wf, [
         ('outputnode.motcorr_params', 'inputnode.motcorr_params')]),
        (inputnode, bold_commonspace_trans_wf, [
            ('bold', 'inputnode.name_source'),
            ('commonspace_ref', 'inputnode.ref_file'),
            ]),
        (bold_commonspace_trans_wf, outputnode, [
            ('outputnode.bold', 'commonspace_bold'),
            ('outputnode.brain_mask', 'commonspace_mask'),
            ('outputnode.WM_mask', 'commonspace_WM_mask'),
            ('outputnode.CSF_mask', 'commonspace_CSF_mask'),
            ('outputnode.vascular_mask', 'commonspace_vascular_mask'),
            ('outputnode.labels', 'commonspace_labels'),
            ('outputnode.raw_brain_mask', 'raw_brain_mask'),
            ]),
        ])

    if opts.apply_slice_mc:
        workflow.connect([
            (bold_stc_wf, bold_hmc_wf, [
             ('outputnode.stc_file', 'inputnode.bold_file')]),
            (bold_hmc_wf, bold_commonspace_trans_wf, [
             ('outputnode.slice_corrected_bold', 'inputnode.bold_file')]),
        ])
        if not opts.bold_only:
            workflow.connect([
                (bold_hmc_wf, bold_native_trans_wf, [
                 ('outputnode.slice_corrected_bold', 'inputnode.bold_file')]),
            ])
    else:
        workflow.connect([
            (bold_stc_wf, bold_commonspace_trans_wf, [
             ('outputnode.stc_file', 'inputnode.bold_file')]),
        ])
        if not opts.bold_only:
            workflow.connect([
                (bold_stc_wf, bold_native_trans_wf, [
                 ('outputnode.stc_file', 'inputnode.bold_file')]),
            ])
        if opts.isotropic_HMC:
            workflow.connect([
                (transitionnode, bold_hmc_wf, [
                    ('isotropic_bold_file', 'inputnode.bold_file')]),
                ])
        else:
            workflow.connect([
                (transitionnode, bold_hmc_wf, [
                    ('bold_file', 'inputnode.bold_file')]),
                ])


    if opts.isotropic_HMC:
        workflow.connect([
            (transitionnode, estimate_motion_node, [
                ('isotropic_bold_file', 'raw_bold'),
                ]),
            ])
    else:
        workflow.connect([
            (inputnode, estimate_motion_node, [
                ('bold', 'raw_bold'),
                ]),
            ])

    if opts.voxelwise_motion:
        workflow.connect([
            (estimate_motion_node, outputnode, [
                ('FD_voxelwise', 'FD_voxelwise'),
                ('pos_voxelwise', 'pos_voxelwise'),
                ]),
            ])

    return workflow
