from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, afni

from .hmc import init_bold_hmc_wf
from .utils import init_bold_reference_wf
from .resampling import init_bold_preproc_trans_wf, init_bold_commonspace_trans_wf
from .stc import init_bold_stc_wf
from .bias_correction import bias_correction_wf
from .registration import init_bold_reg_wf
from .confounds import init_bold_confs_wf
from nipype.interfaces.utility import Function


def init_bold_main_wf(opts, bias_cor_only=False, name='bold_main_wf'):
    """
    This workflow controls the functional preprocessing stages of the pipeline when both
    functional and anatomical images are provided.

    **Parameters**

        apply_despiking
            whether to apply despiking using AFNI's 3dDespike https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html.
        tr
            repetition time for the EPI
        tpattern
            specification for the within TR slice acquisition method. The input is fed to AFNI's 3dTshift
        no_STC
            whether to apply slice timing correction (STC) or not
        detect_dummy
            whether to detect and remove dummy volumes at the beginning of the EPI Sequences
        slice_mc
            whether to apply slice-specific motion correction through 2D registration of each slice, which can improve the correction
            of within-TR motion
        coreg_script
            path to registration script for EPI to anat coregistraion. The script must
            follow the template structure of registration scripts in shell_scripts/.
            Default is set to 'SyN' registration.
        nativespace_resampling
            Specified dimensions for the resampling of the corrected EPI in native space.
        commonspace_resampling
            Specified dimensions for the resampling of the corrected EPI in common space.

    **Inputs**

        bold
            Input BOLD series NIfTI file
        anat_ref
            Preprocessed anatomical image after bias field correction and denoising
        anat_mask
            Brain mask inherited from the common space registration
        WM_mask
            Eroded WM mask inherited from the common space registration
        CSF_mask
            Eroded CSF mask inherited from the common space registration
        labels
            Anatomical labels inherited from the common space registration
        commonspace_transforms_list
            list of transforms to be applied to resample to commonspace
        commonspace_inverses
            Specification for the application of inverse affine transform for
            the provided commonspace transforms

    **Outputs**

        input_bold
            The provided input BOLD file
        bold_ref
            Initial EPI median volume subsequently used as 3D reference EPI volume
        motcorr_params
            motion parameters file provided from antsMotionCorr
        corrected_EPI
            3D reference EPI volume after bias field correction
        itk_bold_to_anat
            Composite transforms from the EPI space to the anatomical space
        itk_anat_to_bold
            Composite transforms from the anatomical space to the EPI space
        output_warped_bold
            Bias field corrected 3D EPI volume warped to the anatomical space
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
        commonspace_labels
            EPI anatomical labels for commonspace bold
    """

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['subject_id', 'bold', 'anat_ref', 'anat_mask', 'WM_mask', 'CSF_mask', 'vascular_mask', 'labels', 'template_to_common_affine', 'template_to_common_warp', 'anat_to_template_affine', 'anat_to_template_warp', 'template_anat']),
                        name="inputnode")

    outputnode = pe.Node(niu.IdentityInterface(
                fields=['input_bold', 'bold_ref', 'motcorr_params', 'init_denoise', 'denoise_mask', 'corrected_EPI', 'output_warped_bold', 'affine_bold2anat', 'warp_bold2anat', 'inverse_warp_bold2anat', 'resampled_bold', 'resampled_ref_bold', 'EPI_brain_mask', 'EPI_WM_mask', 'EPI_CSF_mask', 'EPI_labels',
                        'confounds_csv', 'FD_voxelwise', 'pos_voxelwise', 'FD_csv', 'commonspace_bold', 'commonspace_mask', 'commonspace_WM_mask', 'commonspace_CSF_mask', 'commonspace_vascular_mask', 'commonspace_labels']),
                name='outputnode')

    boldbuffer = pe.Node(niu.IdentityInterface(fields=['bold_file']),
                         name="boldbuffer")

    # this node will serve as a relay of outputs from the bias_cor main_wf to the inputs for the rest of the main_wf for bold_only
    transitionnode = pe.Node(niu.IdentityInterface(fields=['bold_file', 'bold_ref', 'init_denoise', 'denoise_mask', 'corrected_EPI']),
                             name="transitionnode")

    if bias_cor_only or (not opts.bold_only):
        bold_reference_wf = init_bold_reference_wf(opts=opts)
        bias_cor_wf = bias_correction_wf(opts=opts)

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
            (inputnode, bias_cor_wf, [
                ('anat_ref', 'inputnode.anat'),
                ('anat_mask', 'inputnode.anat_mask'),
                ('bold', 'inputnode.name_source'),
                ]),
            (boldbuffer, bold_reference_wf, [
                ('bold_file', 'inputnode.bold_file'),
                ]),
            (bold_reference_wf, bias_cor_wf, [
                ('outputnode.ref_image', 'inputnode.ref_EPI'),
                ]),
            (bold_reference_wf, transitionnode, [
                ('outputnode.ref_image', 'bold_ref'),
                ]),
            (bias_cor_wf, transitionnode, [
                ('outputnode.init_denoise', 'init_denoise'),
                ('outputnode.denoise_mask', 'denoise_mask'),
                ('outputnode.corrected_EPI', 'corrected_EPI'),
                ]),
            ])

    if opts.bold_only and bias_cor_only:
        return workflow

    bold_stc_wf = init_bold_stc_wf(opts=opts)

    # HMC on the BOLD
    bold_hmc_wf = init_bold_hmc_wf(opts=opts)

    if not opts.bold_only:
        def commonspace_transforms(template_to_common_warp, template_to_common_affine, anat_to_template_warp, anat_to_template_affine, warp_bold2anat, affine_bold2anat):
            # transforms_list,inverses
            return [template_to_common_warp, template_to_common_affine, anat_to_template_warp, anat_to_template_affine, warp_bold2anat, affine_bold2anat], [0, 0, 0, 0, 0, 0]
        commonspace_transforms_prep = pe.Node(Function(input_names=['template_to_common_warp', 'template_to_common_affine', 'anat_to_template_warp', 'anat_to_template_affine', 'warp_bold2anat', 'affine_bold2anat'],
                                                       output_names=[
                                                           'transforms_list', 'inverses'],
                                                       function=commonspace_transforms),
                                              name='commonspace_transforms_prep')
    else:
        def commonspace_transforms(template_to_common_warp, template_to_common_affine, anat_to_template_warp, anat_to_template_affine):
            # transforms_list,inverses
            return [template_to_common_warp, template_to_common_affine, anat_to_template_warp, anat_to_template_affine], [0, 0, 0, 0]
        commonspace_transforms_prep = pe.Node(Function(input_names=['template_to_common_warp', 'template_to_common_affine', 'anat_to_template_warp', 'anat_to_template_affine', ],
                                                       output_names=[
                                                           'transforms_list', 'inverses'],
                                                       function=commonspace_transforms),
                                              name='commonspace_transforms_prep')

    bold_commonspace_trans_wf = init_bold_commonspace_trans_wf(opts=opts)

    bold_confs_wf = init_bold_confs_wf(opts=opts, name="bold_confs_wf")

    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (inputnode, commonspace_transforms_prep, [
            ("template_to_common_affine", "template_to_common_affine"),
            ("template_to_common_warp", "template_to_common_warp"),
            ("anat_to_template_affine", "anat_to_template_affine"),
            ("anat_to_template_warp", "anat_to_template_warp"),
            ]),
        (inputnode, bold_confs_wf, [('anat_mask', 'inputnode.t1_mask'),
                                    ('WM_mask', 'inputnode.WM_mask'),
                                    ('CSF_mask', 'inputnode.CSF_mask'),
                                    ('vascular_mask', 'inputnode.vascular_mask'),
                                    ('labels', 'inputnode.t1_labels'),
                                    ('bold', 'inputnode.name_source'),
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
            ('outputnode.brain_mask', 'EPI_brain_mask'),
            ('outputnode.WM_mask', 'EPI_WM_mask'),
            ('outputnode.CSF_mask', 'EPI_CSF_mask'),
            ('outputnode.EPI_labels', 'EPI_labels'),
            ('outputnode.confounds_csv', 'confounds_csv'),
            ('outputnode.FD_csv', 'FD_csv'),
            ('outputnode.FD_voxelwise', 'FD_voxelwise'),
            ('outputnode.pos_voxelwise', 'pos_voxelwise'),
            ]),
        (commonspace_transforms_prep, bold_commonspace_trans_wf, [
            ('transforms_list', 'inputnode.transforms_list'),
            ('inverses', 'inputnode.inverses'),
            ]),
        (bold_hmc_wf, bold_commonspace_trans_wf, [
         ('outputnode.motcorr_params', 'inputnode.motcorr_params')]),
        (inputnode, bold_commonspace_trans_wf, [
            ('bold', 'inputnode.name_source'),
            ('template_anat', 'inputnode.ref_file'),
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

    if not opts.bold_only:
        bold_reg_wf = init_bold_reg_wf(opts=opts)

        def SyN_coreg_transforms_prep(warp_bold2anat, affine_bold2anat):
            # transforms_list,inverses
            return [warp_bold2anat, affine_bold2anat], [0, 0]
        transforms_prep = pe.Node(Function(input_names=['warp_bold2anat', 'affine_bold2anat'],
                                           output_names=[
                                               'transforms_list', 'inverses'],
                                           function=SyN_coreg_transforms_prep),
                                  name='transforms_prep')

        # Apply transforms in 1 shot
        bold_bold_trans_wf = init_bold_preproc_trans_wf(opts=opts)

        workflow.connect([
            (inputnode, bold_reg_wf, [
                ('anat_ref', 'inputnode.anat_ref'),
                ('anat_mask', 'inputnode.anat_mask')]),
            (inputnode, bold_bold_trans_wf, [
                ('bold', 'inputnode.name_source')]),
            (transitionnode, bold_reg_wf, [
                ('corrected_EPI', 'inputnode.ref_bold_brain')]),
            (bold_reg_wf, outputnode, [
                ('outputnode.affine_bold2anat', 'affine_bold2anat'),
                ('outputnode.warp_bold2anat', 'warp_bold2anat'),
                ('outputnode.inverse_warp_bold2anat', 'inverse_warp_bold2anat'),
                ('outputnode.output_warped_bold', 'output_warped_bold'),
                ]),
            (bold_reg_wf, transforms_prep, [
                ('outputnode.affine_bold2anat', 'affine_bold2anat'),
                ('outputnode.warp_bold2anat', 'warp_bold2anat'),
                ]),
            (transforms_prep, bold_bold_trans_wf, [
                ('transforms_list', 'inputnode.transforms_list'),
                ('inverses', 'inputnode.inverses'),
                ]),
            (bold_reg_wf, bold_bold_trans_wf, [
                ('outputnode.output_warped_bold', 'inputnode.ref_file')]),
            (bold_reg_wf, commonspace_transforms_prep, [
                ('outputnode.affine_bold2anat', 'affine_bold2anat'),
                ('outputnode.warp_bold2anat', 'warp_bold2anat'),
                ]),
            (bold_hmc_wf, bold_bold_trans_wf, [
             ('outputnode.motcorr_params', 'inputnode.motcorr_params')]),
            (bold_bold_trans_wf, outputnode, [
                ('outputnode.bold_ref', 'resampled_ref_bold'),
                ('outputnode.bold', 'resampled_bold'),
                ]),
            (bold_bold_trans_wf, bold_confs_wf, [('outputnode.bold', 'inputnode.bold'),
                                                 ('outputnode.bold_ref',
                                                  'inputnode.ref_bold'),
                                                 ]),
            ])
    else:
        workflow.connect([
            (bold_commonspace_trans_wf, bold_confs_wf, [('outputnode.bold', 'inputnode.bold'),
                                                        ('outputnode.bold_ref',
                                                         'inputnode.ref_bold'),
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
                (bold_hmc_wf, bold_bold_trans_wf, [
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
                (bold_stc_wf, bold_bold_trans_wf, [
                 ('outputnode.stc_file', 'inputnode.bold_file')]),
            ])

    return workflow
