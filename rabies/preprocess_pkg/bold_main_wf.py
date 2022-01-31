from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, afni

from .hmc import init_bold_hmc_wf
from .bold_ref import init_bold_reference_wf
from .resampling import init_bold_preproc_trans_wf
from .stc import init_bold_stc_wf
from .inho_correction import init_inho_correction_wf
from .registration import init_cross_modal_reg_wf
from .confounds import init_bold_confs_wf
from nipype.interfaces.utility import Function


def init_bold_main_wf(opts, inho_cor_only=False, name='bold_main_wf'):
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
        confounds_csv
            .csv file with measured confound timecourses, including global signal,
            WM signal, CSF signal, 6 rigid body motion parameters + their first
            temporal derivate + the 12 parameters squared (24 motion parameters),
            and aCompCorr timecourses
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
                        'native_bold', 'native_bold_ref', 'native_brain_mask', 'native_WM_mask', 'native_CSF_mask', 'native_labels',
                        'confounds_csv', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv', 'commonspace_bold', 'commonspace_mask',
                        'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_vascular_mask', 'commonspace_labels']),
                name='outputnode')

    boldbuffer = pe.Node(niu.IdentityInterface(fields=['bold_file']),
                         name="boldbuffer")

    # this node will serve as a relay of outputs from the inho_cor main_wf to the inputs for the rest of the main_wf for bold_only
    transitionnode = pe.Node(niu.IdentityInterface(fields=['bold_file', 'bold_ref', 'init_denoise', 'denoise_mask', 'corrected_EPI']),
                             name="transitionnode")

    if inho_cor_only or (not opts.bold_only):
        bold_reference_wf = init_bold_reference_wf(opts=opts)
        inho_cor_wf = init_inho_correction_wf(opts=opts, image_type='EPI', name="bold_inho_cor_wf")

        if opts.apply_despiking:
            despike = pe.Node(
                afni.Despike(outputtype='NIFTI_GZ'),
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
            workflow.connect([
                (bold_reference_wf, transitionnode, [
                    ('outputnode.bold_file', 'bold_file'),
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
            (boldbuffer, bold_reference_wf, [
                ('bold_file', 'inputnode.bold_file'),
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

    bold_commonspace_trans_wf = init_bold_preproc_trans_wf(opts=opts, resampling_dim=opts.commonspace_resampling, name='bold_commonspace_trans_wf')
    bold_commonspace_trans_wf.inputs.inputnode.mask_transforms_list = []
    bold_commonspace_trans_wf.inputs.inputnode.mask_inverses = []

    bold_confs_wf = init_bold_confs_wf(opts=opts, name="bold_confs_wf")

    if not opts.bold_only:
        def commonspace_transforms(to_commonspace_transform_list,to_commonspace_inverse_list, bold_to_anat_warp, bold_to_anat_affine):
            # transforms_list,inverses
            return to_commonspace_transform_list+[bold_to_anat_warp, bold_to_anat_affine], to_commonspace_inverse_list+[0,0]
        bold_to_commonspace_transforms = pe.Node(Function(input_names=['to_commonspace_transform_list','to_commonspace_inverse_list', 'bold_to_anat_warp', 'bold_to_anat_affine'],
                                                       output_names=[
                                                           'to_commonspace_transform_list','to_commonspace_inverse_list'],
                                                       function=commonspace_transforms),
                                              name='bold_to_commonspace_transforms')

        cross_modal_reg_wf = init_cross_modal_reg_wf(opts=opts)

        def SyN_coreg_transforms_prep(bold_to_anat_warp, bold_to_anat_affine):
            # transforms_list,inverses
            return [bold_to_anat_warp, bold_to_anat_affine], [0, 0]
        transforms_prep = pe.Node(Function(input_names=['bold_to_anat_warp', 'bold_to_anat_affine'],
                                           output_names=[
                                               'transforms_list', 'inverses'],
                                           function=SyN_coreg_transforms_prep),
                                  name='transforms_prep')

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
            (cross_modal_reg_wf, transforms_prep, [
                ('outputnode.bold_to_anat_affine', 'bold_to_anat_affine'),
                ('outputnode.bold_to_anat_warp', 'bold_to_anat_warp'),
                ]),
            (transforms_prep, bold_native_trans_wf, [
                ('transforms_list', 'inputnode.transforms_list'),
                ('inverses', 'inputnode.inverses'),
                ]),
            (cross_modal_reg_wf, bold_native_trans_wf, [
                ('outputnode.output_warped_bold', 'inputnode.ref_file')]),
            (cross_modal_reg_wf, bold_to_commonspace_transforms, [
                ('outputnode.bold_to_anat_affine', 'bold_to_anat_affine'),
                ('outputnode.bold_to_anat_warp', 'bold_to_anat_warp'),
                ]),
            (bold_hmc_wf, bold_native_trans_wf, [
             ('outputnode.motcorr_params', 'inputnode.motcorr_params')]),
            (bold_native_trans_wf, bold_confs_wf, [
                ('outputnode.bold', 'inputnode.bold'),
                ('outputnode.bold_ref','inputnode.ref_bold'),
                ('outputnode.brain_mask', 'inputnode.brain_mask'),
                ('outputnode.WM_mask', 'inputnode.WM_mask'),
                ('outputnode.CSF_mask', 'inputnode.CSF_mask'),
                ('outputnode.vascular_mask', 'inputnode.vascular_mask'),
                ]),
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
        bold_to_commonspace_transforms = pe.Node(niu.IdentityInterface(fields=['to_commonspace_transform_list','to_commonspace_inverse_list']),
                                 name="bold_to_commonspace_transforms")

        workflow.connect([
            (bold_commonspace_trans_wf, bold_confs_wf, [
                ('outputnode.bold', 'inputnode.bold'),
                ('outputnode.bold_ref','inputnode.ref_bold'),
                ('outputnode.brain_mask', 'inputnode.brain_mask'),
                ('outputnode.WM_mask', 'inputnode.WM_mask'),
                ('outputnode.CSF_mask', 'inputnode.CSF_mask'),
                ('outputnode.vascular_mask', 'inputnode.vascular_mask'),
                ]),
            ])


    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (inputnode, bold_to_commonspace_transforms, [
            ('native_to_commonspace_transform_list', 'to_commonspace_transform_list'),
            ('native_to_commonspace_inverse_list', 'to_commonspace_inverse_list'),
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
        (bold_hmc_wf, bold_confs_wf, [
            ('outputnode.motcorr_params', 'inputnode.movpar_file'),
            ]),
        (bold_confs_wf, outputnode, [
            ('outputnode.confounds_csv', 'confounds_csv'),
            ('outputnode.FD_csv', 'FD_csv'),
            ('outputnode.FD_voxelwise', 'FD_voxelwise'),
            ('outputnode.pos_voxelwise', 'pos_voxelwise'),
            ]),
        (bold_to_commonspace_transforms, bold_commonspace_trans_wf, [
            ('to_commonspace_transform_list', 'inputnode.transforms_list'),
            ('to_commonspace_inverse_list', 'inputnode.inverses'),
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
            (transitionnode, bold_hmc_wf, [
                ('bold_file', 'inputnode.bold_file')]),
            (bold_stc_wf, bold_commonspace_trans_wf, [
             ('outputnode.stc_file', 'inputnode.bold_file')]),
        ])
        if not opts.bold_only:
            workflow.connect([
                (bold_stc_wf, bold_native_trans_wf, [
                 ('outputnode.stc_file', 'inputnode.bold_file')]),
            ])

    return workflow
