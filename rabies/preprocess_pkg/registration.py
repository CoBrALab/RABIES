from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function


def init_bold_reg_wf(coreg_script='SyN', rabies_data_type=8, rabies_mem_scale=1.0, min_proc=1, name='bold_reg_wf'):
    """
    This workflow registers the reference BOLD image to anat-space, using
    antsRegistration, either applying Affine registration only, or the
    combination of Affine followed by non-linear registration using the SyN
    algorithm, which may apply distortion correction through the registration
    to the structural image.

    **Parameters**

        SyN_reg : bool
            Determine whether SyN registration will used or not, and uses the
            transform from the registration as SDC transforms to transform from
            bold to anat

    **Inputs**

        ref_bold_brain
            Reference image to which BOLD series is aligned
        anat_preproc
            Bias-corrected structural template image
        anat_mask
            Mask of the skull-stripped template image

    **Outputs**

        itk_bold_to_anat
            Transform from ``ref_bold_brain`` to anat space
        itk_anat_to_bold
            Transform from anat space to BOLD space (ITK format)
        output_warped_bold
            output warped image from antsRegistration

    """

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['ref_bold_brain', 'anat_preproc', 'anat_mask']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=[
            'affine_bold2anat', 'warp_bold2anat', 'inverse_warp_bold2anat', 'output_warped_bold']),
        name='outputnode'
    )

    run_reg = pe.Node(Function(input_names=["reg_method", "moving_image", "fixed_image",
                                            "anat_mask", "rabies_data_type"],
                               output_names=['affine_bold2anat', 'warp_bold2anat',
                                             'inverse_warp_bold2anat', 'output_warped_bold'],
                               function=run_antsRegistration), name='EPI_Coregistration', mem_gb=3*rabies_mem_scale)
    run_reg.inputs.reg_method = coreg_script
    run_reg.inputs.rabies_data_type = rabies_data_type
    run_reg.plugin_args = {
        'qsub_args': '-pe smp %s' % (str(3*min_proc)), 'overwrite': True}

    workflow.connect([
        (inputnode, run_reg, [
            ('ref_bold_brain', 'moving_image'),
            ('anat_preproc', 'fixed_image'),
            ('anat_mask', 'anat_mask')]),
        (run_reg, outputnode, [
            ('affine_bold2anat', 'affine_bold2anat'),
            ('warp_bold2anat', 'warp_bold2anat'),
            ('inverse_warp_bold2anat', 'inverse_warp_bold2anat'),
            ('output_warped_bold', 'output_warped_bold'),
            ]),
        ])

    return workflow


def run_antsRegistration(reg_method, moving_image='NULL', fixed_image='NULL', anat_mask='NULL', rabies_data_type=8):
    import os
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(moving_image).name.rsplit(".nii")

    from rabies.preprocess_pkg.registration import define_reg_script
    reg_script = define_reg_script(reg_method)

    command = 'bash %s %s %s %s %s %s' % (
        reg_script, moving_image, fixed_image, anat_mask, filename_split[0], reg_method)
    from rabies.preprocess_pkg.utils import run_command
    rc = run_command(command)

    cwd = os.getcwd()
    warped_image = '%s/%s_output_warped_image.nii.gz' % (
        cwd, filename_split[0],)
    affine = '%s/%s_output_0GenericAffine.mat' % (cwd, filename_split[0],)
    warp = '%s/%s_output_1Warp.nii.gz' % (cwd, filename_split[0],)
    inverse_warp = '%s/%s_output_1InverseWarp.nii.gz' % (
        cwd, filename_split[0],)
    if not os.path.isfile(warped_image) or not os.path.isfile(affine):
        raise ValueError(
            'REGISTRATION ERROR: OUTPUT FILES MISSING. Make sure the provided registration script runs properly.')
    if not os.path.isfile(warp) or not os.path.isfile(inverse_warp):
        print('No non-linear warp files as output. Assumes linear registration.')
        warp = 'NULL'
        inverse_warp = 'NULL'

    import SimpleITK as sitk
    sitk.WriteImage(sitk.ReadImage(warped_image, rabies_data_type), warped_image)

    return [affine, warp, inverse_warp, warped_image]


def define_reg_script(reg_option):
    import os
    import rabies
    dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
    if reg_option == 'Rigid' or reg_option == 'Affine' or reg_option == 'SyN':
        reg_script = dir_path+'/shell_scripts/generic_registration.sh'
    elif reg_option == 'light_SyN':
        reg_script = dir_path+'/shell_scripts/light_SyN_registration.sh'
    elif reg_option == 'heavy_SyN':
        reg_script = dir_path+'/shell_scripts/heavy_SyN_registration.sh'
    elif reg_option == 'multiRAT':
        reg_script = dir_path+'/shell_scripts/multiRAT_registration.sh'
    else:
        '''
        For user-provided antsRegistration command.
        '''
        if os.path.isfile(reg_option):
            reg_script = reg_option
        else:
            raise ValueError(
                'REGISTRATION ERROR: THE REG SCRIPT FILE DOES NOT EXISTS')
    return reg_script
