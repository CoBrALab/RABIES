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

def run_registration_interface(interface):
    from nipype.interfaces.base import CommandLine
    res=interface.run()
    return [res.outputs.composite_transform, res.outputs.inverse_composite_transform, res.outputs.warped_image]


def run_antsRegistration(reg_script='Affine', moving_image='NULL', fixed_image='NULL', anat_mask='NULL'):
    import os
    subject_id=os.path.basename(moving_image).split('_ses-')[0]
    session=os.path.basename(moving_image).split('_ses-')[1][0]
    run=os.path.basename(moving_image).split('_run-')[1][0]
    filename_template = '%s_ses-%s_run-%s' % (subject_id, session, run)

    if reg_script=='SyN':
        import mfp
        dir_path = os.path.dirname(os.path.realpath(mfp.__file__))
        reg_script_path=dir_path+'/shell_scripts/SyN_registration.sh'
        '''
        EPI=$1
        anat_file=$2
        mask=$3
        filename_template=$4

        antsRegistration -d 3 \
        --verbose -o [${filename_template}_output_,${filename_template}_output_warped_image.nii.gz] \
        -t Rigid[0.1] -m Mattes[$anat_file,$EPI,1,64,None] \
        -c 1000x500x250x100x50x25 -s 8x4x2x1x0.5x0 -f 6x5x4x3x2x1 --masks [NULL,NULL] \
        -t Similarity[0.1] -m Mattes[$anat_file,$EPI,1,64,None] \
        -c 100x50x25 -s 1x0.5x0 -f 3x2x1 --masks [$mask,NULL] \
        -t Affine[0.1] -m Mattes[$anat_file,$EPI,1,64,None] \
        -c 100x50x25 -s 1x0.5x0 -f 3x2x1 --masks [$mask,$mask] \
        -t SyN[ 0.2, 3.0, 0.0 ] -m Mattes[$anat_file,$EPI, 1, 64 ] \
        -c [ 40x20x0, 1e-06, 6 ] -s 2x1x0 -f 4x2x1 --masks [$mask,$mask] \
        --interpolation BSpline[5] -z 1 -u 0 -a 1
        '''

    elif reg_script=='Affine':
        import mfp
        dir_path = os.path.dirname(os.path.realpath(mfp.__file__))
        reg_script_path=dir_path+'/shell_scripts/Affine_registration.sh'
        '''
        EPI=$1
        anat_file=$2
        mask=$3
        filename_template=$4

        antsRegistration -d 3 \
        --verbose -o [${filename_template}_output_,${filename_template}_output_warped_image.nii.gz] \
        -t Rigid[0.1] -m Mattes[$anat_file,$EPI,1,64,None] \
        -c 1000x500x250x100x50x25 -s 8x4x2x1x0.5x0 -f 6x5x4x3x2x1 --masks [NULL,NULL] \
        -t Similarity[0.1] -m Mattes[$anat_file,$EPI,1,64,None] \
        -c 100x50x25 -s 1x0.5x0 -f 3x2x1 --masks [$mask,NULL] \
        -t Affine[0.1] -m Mattes[$anat_file,$EPI,1,64,None] \
        -c 100x50x25 -s 1x0.5x0 -f 3x2x1 --masks [$mask,$mask] \
        --interpolation BSpline[5] -z 1 -u 0 -a 1
        '''

    elif reg_script=='Rigid':
        import mfp
        dir_path = os.path.dirname(os.path.realpath(mfp.__file__))
        reg_script_path=dir_path+'/shell_scripts/Rigid_registration.sh'
        '''
        EPI=$1
        anat_file=$2
        mask=$3
        filename_template=$4

        antsRegistration -d 3 \
        --verbose -o [${filename_template}_output_,${filename_template}_output_warped_image.nii.gz] \
        -t Rigid[0.1] -m Mattes[$anat_file,$EPI,1,64,None] \
        -c 1000x500x250x100x50x25 -s 8x4x2x1x0.5x0 -f 6x5x4x3x2x1 --masks [NULL,NULL] \
        --interpolation BSpline[5] -z 1 -u 0 -a 1
        '''

    else:
        '''
        For user-provided antsRegistration command.
        '''
        if os.path.isfile(reg_script):
            reg_script_path=reg_script
        else:
            raise ValueError('REGISTRATION ERROR: THE REG SCRIPT FILE DOES NOT EXISTS')
    os.system('bash %s %s %s %s %s' % (reg_script_path,moving_image, fixed_image, anat_mask, filename_template))

    cwd=os.getcwd()
    warped_image='%s/%s_output_warped_image.nii.gz' % (cwd, filename_template)
    inverse_composite_transform='%s/%s_output_InverseComposite.h5' % (cwd, filename_template)
    composite_transform='%s/%s_output_Composite.h5' % (cwd, filename_template)

    return [composite_transform, inverse_composite_transform, warped_image]
