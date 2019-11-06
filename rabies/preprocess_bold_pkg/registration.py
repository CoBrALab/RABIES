import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function

def init_bold_reg_wf(coreg_script='SyN', name='bold_reg_wf'):

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
            'itk_bold_to_anat', 'itk_anat_to_bold', 'output_warped_bold']),
        name='outputnode'
    )


    run_reg = pe.Node(Function(input_names=["reg_script", "moving_image", "fixed_image",
                                            "anat_mask"],
                   output_names=['itk_bold_to_anat', 'itk_anat_to_bold', 'output_warped_bold'],
                   function=run_antsRegistration), name='EPI_Coregistration')
    run_reg.inputs.reg_script=coreg_script
    run_reg.plugin_args = {'qsub_args': '-pe smp 4', 'overwrite': True}

    workflow.connect([
        (inputnode, run_reg, [
            ('ref_bold_brain', 'moving_image'),
            ('anat_preproc', 'fixed_image'),
            ('anat_mask', 'anat_mask')]),
        (run_reg, outputnode, [
            ('itk_bold_to_anat', 'itk_bold_to_anat'),
            ('itk_anat_to_bold', 'itk_anat_to_bold'),
            ('output_warped_bold', 'output_warped_bold'),
            ]),
        ])

    return workflow


def run_antsRegistration(reg_script, moving_image='NULL', fixed_image='NULL', anat_mask='NULL'):
    import os
    filename_template=os.path.basename(moving_image).split('.')[0]

    if os.path.isfile(reg_script):
        reg_script_path=reg_script
    else:
        raise ValueError('REGISTRATION ERROR: THE REG SCRIPT FILE DOES NOT EXISTS')
    os.system('bash %s %s %s %s %s' % (reg_script_path,moving_image, fixed_image, anat_mask, filename_template))

    cwd=os.getcwd()
    warped_image='%s/%s_output_warped_image.nii.gz' % (cwd, filename_template)
    inverse_composite_transform='%s/%s_output_InverseComposite.h5' % (cwd, filename_template)
    composite_transform='%s/%s_output_Composite.h5' % (cwd, filename_template)
    if not os.path.isfile(warped_image) or not os.path.isfile(inverse_composite_transform) or not os.path.isfile(composite_transform):
        raise ValueError('REGISTRATION ERROR: OUTPUT FILES MISSING.')

    return [composite_transform, inverse_composite_transform, warped_image]
