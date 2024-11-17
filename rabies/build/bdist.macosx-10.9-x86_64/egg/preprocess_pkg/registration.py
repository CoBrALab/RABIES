from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function


def init_cross_modal_reg_wf(opts, name='cross_modal_reg_wf'):
    # cross_modal_reg_head_start
    """
    The input volumetric EPI image is registered non-linearly to an associated structural MRI image.
    The non-linear transform estimates the correction for EPI susceptibility distortions (Wang et al., 2017).

    References:
        Wang, S., Peterson, D. J., Gatenby, J. C., Li, W., Grabowski, T. J., & Madhyastha, T. M. (2017). 
            Evaluation of Field Map and Nonlinear Registration Methods for Correction of Susceptibility Artifacts 
            in Diffusion MRI. Frontiers in Neuroinformatics, 11, 17.

    Command line interface parameters:                  
        --bold2anat_coreg BOLD2ANAT_COREG
                            Specify the registration script for cross-modal alignment between the EPI and structural
                            images. This operation is responsible for correcting EPI susceptibility distortions.
                            * masking: With this option, the brain masks obtained from the EPI inhomogeneity correction 
                            step are used to support registration.
                            *** Specify 'true' or 'false'. 
                            * brain_extraction: conducts brain extraction prior to registration using the EPI masks from 
                            inhomogeneity correction. This will enhance brain edge-matching, but requires good quality 
                            masks. This should be selected along the 'masking' option.
                            *** Specify 'true' or 'false'. 
                            * keep_mask_after_extract: If using brain_extraction, use the mask to compute the registration metric
                            within the mask only. Choose to prevent stretching of the images beyond the limit of the brain mask
                            (e.g. if the moving and target images don't have the same brain coverage).
                            *** Specify 'true' or 'false'.
                            * registration: Specify a registration script.
                            *** Rigid: conducts only rigid registration.
                            *** Affine: conducts Rigid then Affine registration.
                            *** SyN: conducts Rigid, Affine then non-linear registration.
                            *** no_reg: skip registration.
                            (default: masking=false,brain_extraction=false,keep_mask_after_extract=false,registration=SyN)

    Workflow:
        parameters
            opts: command line interface parameters

        inputs
            ref_bold_brain: volumetric EPI image to register
            anat_ref: the target structural image
            anat_mask: the brain mask of the structural image
            moving_mask: a EPI mask inherited from inhomogeneity correction

        outputs
            bold_to_anat_affine: affine transform from the EPI to the anatomical image
            bold_to_anat_warp: non-linear transform from the EPI to the anatomical image
            bold_to_anat_inverse_warp: inverse non-linear transform from the EPI to the anatomical image
            output_warped_bold: the EPI image warped onto the structural image
    """
    # cross_modal_reg_head_end

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['ref_bold_brain', 'anat_ref', 'anat_mask', 'moving_mask']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=[
            'bold_to_anat_affine', 'bold_to_anat_warp', 'bold_to_anat_inverse_warp', 'output_warped_bold']),
        name='outputnode'
    )

    run_reg = pe.Node(Function(input_names=["reg_method", "brain_extraction", "keep_mask_after_extract", "moving_image", "moving_mask", "fixed_image",
                                            "fixed_mask", "rabies_data_type"],
                               output_names=['bold_to_anat_affine', 'bold_to_anat_warp',
                                             'bold_to_anat_inverse_warp', 'output_warped_bold'],
                               function=run_antsRegistration), name='EPI_Coregistration', mem_gb=3*opts.scale_min_memory)

    # don't use brain extraction without a moving mask
    run_reg.inputs.reg_method = opts.bold2anat_coreg['registration']
    run_reg.inputs.brain_extraction = opts.bold2anat_coreg['brain_extraction']
    run_reg.inputs.keep_mask_after_extract = opts.bold2anat_coreg['keep_mask_after_extract']
    run_reg.inputs.rabies_data_type = opts.data_type
    run_reg.plugin_args = {
        'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}

    if opts.bold2anat_coreg['masking']:
        workflow.connect([
            (inputnode, run_reg, [
                ('moving_mask', 'moving_mask')]),
            ])


    workflow.connect([
        (inputnode, run_reg, [
            ('ref_bold_brain', 'moving_image'),
            ('anat_ref', 'fixed_image'),
            ('anat_mask', 'fixed_mask')]),
        (run_reg, outputnode, [
            ('bold_to_anat_affine', 'bold_to_anat_affine'),
            ('bold_to_anat_warp', 'bold_to_anat_warp'),
            ('bold_to_anat_inverse_warp', 'bold_to_anat_inverse_warp'),
            ('output_warped_bold', 'output_warped_bold'),
            ]),
        ])

    return workflow


def run_antsRegistration(reg_method, brain_extraction=False, keep_mask_after_extract=False, moving_image='NULL', moving_mask='NULL', fixed_image='NULL', fixed_mask='NULL', rabies_data_type=8):
    import os
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(moving_image).name.rsplit(".nii")

    from rabies.preprocess_pkg.registration import define_reg_script
    reg_call = define_reg_script(reg_method)

    if reg_method == 'Rigid' or reg_method == 'Affine' or reg_method == 'SyN':
        if brain_extraction:
            reg_call+=" --mask-extract"
            if keep_mask_after_extract:
                reg_call+=" --keep-mask-after-extract"
        command = f"{reg_call} --moving-mask {moving_mask} --fixed-mask {fixed_mask} --resampled-output {filename_split[0]}_output_warped_image.nii.gz {moving_image} {fixed_image} {filename_split[0]}_output_"
    else:
        command = f'{reg_call} {moving_image} {moving_mask} {fixed_image} {fixed_mask} {filename_split[0]}'
    from rabies.utils import run_command
    rc,c_out = run_command(command)

    cwd = os.getcwd()
    warped_image = f'{cwd}/{filename_split[0]}_output_warped_image.nii.gz'
    affine = f'{cwd}/{filename_split[0]}_output_0GenericAffine.mat'
    warp = f'{cwd}/{filename_split[0]}_output_1Warp.nii.gz'
    inverse_warp = f'{cwd}/{filename_split[0]}_output_1InverseWarp.nii.gz'
    if not os.path.isfile(warped_image) or not os.path.isfile(affine):
        raise ValueError(
            'REGISTRATION ERROR: OUTPUT FILES MISSING. Make sure the provided registration script runs properly.')
    if not os.path.isfile(warp) or not os.path.isfile(inverse_warp):
        from nipype import logging
        log = logging.getLogger('nipype.workflow')
        log.debug('No non-linear warp files as output. Assumes linear registration.')
        warp = 'NULL'
        inverse_warp = 'NULL'

    import SimpleITK as sitk
    sitk.WriteImage(sitk.ReadImage(warped_image, rabies_data_type), warped_image)

    return [affine, warp, inverse_warp, warped_image]


def define_reg_script(reg_option):
    import os
    if reg_option == 'Rigid':
        reg_call = "antsRegistration_affine_SyN.sh --linear-type rigid --skip-nonlinear"
    elif reg_option == 'Affine':
        reg_call = "antsRegistration_affine_SyN.sh --linear-type affine --skip-nonlinear"
    elif reg_option == 'SyN':
        reg_call = "antsRegistration_affine_SyN.sh --linear-type affine"
    elif reg_option == 'no_reg':
        reg_call = 'null_nonlin.sh'
    else:
        raise ValueError(
            'The registration option must be among Rigid,Affine,SyN or NULL.')
    return reg_call
