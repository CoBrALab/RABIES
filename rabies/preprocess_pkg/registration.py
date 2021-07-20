from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function


def init_cross_modal_reg_wf(opts, name='cross_modal_reg_wf'):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['ref_bold_brain', 'anat_ref', 'anat_mask']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=[
            'affine_bold2anat', 'warp_bold2anat', 'inverse_warp_bold2anat', 'output_warped_bold']),
        name='outputnode'
    )

    run_reg = pe.Node(Function(input_names=["reg_method", "moving_image", "moving_mask", "fixed_image",
                                            "fixed_mask", "rabies_data_type"],
                               output_names=['affine_bold2anat', 'warp_bold2anat',
                                             'inverse_warp_bold2anat', 'output_warped_bold'],
                               function=run_antsRegistration), name='EPI_Coregistration', mem_gb=3*opts.scale_min_memory)
    run_reg.inputs.reg_method = opts.coreg_script
    run_reg.inputs.rabies_data_type = opts.data_type
    run_reg.plugin_args = {
        'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}

    workflow.connect([
        (inputnode, run_reg, [
            ('ref_bold_brain', 'moving_image'),
            ('anat_ref', 'fixed_image'),
            ('anat_mask', 'fixed_mask')]),
        (run_reg, outputnode, [
            ('affine_bold2anat', 'affine_bold2anat'),
            ('warp_bold2anat', 'warp_bold2anat'),
            ('inverse_warp_bold2anat', 'inverse_warp_bold2anat'),
            ('output_warped_bold', 'output_warped_bold'),
            ]),
        ])

    return workflow


def run_antsRegistration(reg_method, moving_image='NULL', moving_mask='NULL', fixed_image='NULL', fixed_mask='NULL', rabies_data_type=8):
    import os
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(moving_image).name.rsplit(".nii")

    from rabies.preprocess_pkg.registration import define_reg_script
    reg_call = define_reg_script(reg_method)

    print(moving_mask)

    if reg_method == 'Rigid' or reg_method == 'Affine' or reg_method == 'SyN':
        command = f"{reg_call} --moving-mask {moving_mask} --fixed-mask {fixed_mask} --resampled-output {filename_split[0]}_output_warped_image.nii.gz {moving_image} {fixed_image} {filename_split[0]}_output_"
    else:
        command = f'{reg_call} {moving_image} {moving_mask} {fixed_image} {fixed_mask} {filename_split[0]}'
    from rabies.preprocess_pkg.utils import run_command
    rc = run_command(command)

    cwd = os.getcwd()
    warped_image = f'{cwd}/{filename_split[0]}_output_warped_image.nii.gz'
    affine = f'{cwd}/{filename_split[0]}_output_0GenericAffine.mat'
    warp = f'{cwd}/{filename_split[0]}_output_1Warp.nii.gz'
    inverse_warp = f'{cwd}/{filename_split[0]}_output_1InverseWarp.nii.gz'
    if not os.path.isfile(warped_image) or not os.path.isfile(affine):
        raise ValueError(
            'REGISTRATION ERROR: OUTPUT FILES MISSING. Make sure the provided registration script runs properly.')
    if not os.path.isfile(warp) or not os.path.isfile(inverse_warp):
        import logging
        log = logging.getLogger('root')
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
    elif reg_option == 'multiRAT':
        reg_call = 'multiRAT_registration.sh'
    elif reg_option == 'null_nonlin':
        reg_call = 'null_nonlin.sh'
    else:
        '''
        For user-provided antsRegistration command.
        '''
        if os.path.isfile(reg_option):
            reg_call = reg_option
        else:
            raise ValueError(
                'REGISTRATION ERROR: THE REG SCRIPT FILE DOES NOT EXISTS')
    return reg_call
