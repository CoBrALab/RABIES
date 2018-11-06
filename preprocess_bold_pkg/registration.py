import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function

def init_bold_reg_wf(SyN_reg=False, name='bold_reg_wf'):

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

    prep_reg = pe.Node(Function(input_names=["reg_type", "moving_image", "fixed_image",
                                            "anat_mask"],
                   output_names=["interface"],
                   function=setup_antsCoRegistration), name='prep_EPI_Coregistration')
    run_reg = pe.Node(Function(input_names=["interface"],
                   output_names=['itk_bold_to_anat', 'itk_anat_to_bold', 'output_warped_bold'],
                   function=run_registration_interface), name='run_EPI_Coregistration')

    if SyN_reg:
        prep_reg.inputs.reg_type='SyN'
    else:
        prep_reg.inputs.reg_type='Affine'

    workflow.connect([
        (inputnode, prep_reg, [
            ('ref_bold_brain', 'moving_image'),
            ('anat_preproc', 'fixed_image'),
            ('anat_mask', 'anat_mask')]),
        (prep_reg, run_reg, [('interface', 'interface')]),
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


def setup_antsCoRegistration(reg_type='Affine', moving_image='NULL', fixed_image='NULL', anat_mask='NULL'):
    from nipype.interfaces.ants import Registration

    reg = Registration()


    if reg_type=='SyN':
        '''
        antsRegistration -d 3 \
        --verbose -o [output_,output_warped_image.nii.gz] \
        -t Rigid[0.1] -m Mattes[$template,$bias_cor,1,64,None] \
        -c 1000x500x250x100x50x25 -s 8x4x2x1x0.5x0 -f 6x5x4x3x2x1 --masks [NULL,NULL] \
        -t Similarity[0.1] -m Mattes[$template,$bias_cor,1,64,None] \
        -c 100x50x25 -s 1x0.5x0 -f 3x2x1 --masks [$mask,NULL] \
        -t Affine[0.1] -m Mattes[$template,warped.nii.gz,1,64,None] \
        -c 100x50x25 -s 1x0.5x0 -f 3x2x1 --masks [$mask,$mask] \
        -t SyN[ 0.2, 3.0, 0.0 ] -m Mattes[$template,aff_warped.nii.gz, 1, 64 ] \
        -c [ 40x20x0, 1e-06, 6 ] -s 2x1x0 -f 4x2x1 --masks [$mask,$mask] \
        --interpolation BSpline[5]

        '''

        #-d
        reg.inputs.dimension = 3
        #-m
        reg.inputs.metric = ['Mattes']*4
        reg.inputs.metric_weight = [1]*4
        reg.inputs.radius_or_number_of_bins = [64]*4
        #-t
        reg.inputs.transforms = ['Rigid', 'Similarity', 'Affine', 'SyN']
        reg.inputs.transform_parameters = [(0.1,), (0.1,), (0.1,), (0.2, 3.0, 0.0)]
        reg.inputs.sampling_strategy = ['None','None', 'None', 'None']
        #-c
        reg.inputs.number_of_iterations = [[1000, 500, 250, 100, 50, 25], [100, 50, 25], [100, 50, 25], [40, 20, 0]]
        reg.inputs.convergence_window_size = [10, 10, 10, 6]
        #-s
        reg.inputs.smoothing_sigmas = [[8,4,2,1,0.5,0], [1,0.5,0], [1,0.5,0], [2,1,0]]
        #-f
        reg.inputs.shrink_factors = [[6,5,4,3,2,1], [3,2,1], [3,2,1], [4,2,1]]
        #-i 0
        reg.inputs.initialize_transforms_per_stage = False
        #-o
        reg.inputs.output_transform_prefix = "output_"
        #-z 1
        reg.inputs.collapse_output_transforms = True #collapses transforms into a single File
        #interpolation
        reg.inputs.interpolation = 'BSpline'
        reg.inputs.interpolation_parameters = (5,)
        #-u 0
        reg.inputs.use_histogram_matching = False
        #--masks
        reg.inputs.moving_image_masks = ['NULL', 'NULL', anat_mask, anat_mask]
        reg.inputs.fixed_image_masks = ['NULL', anat_mask, anat_mask, anat_mask]

        reg.inputs.moving_image = moving_image
        reg.inputs.fixed_image = fixed_image
        reg.inputs.write_composite_transform=True
        reg.inputs.output_warped_image = 'output_warped_image.nii.gz'

    elif reg_type=='Affine':
        '''
        antsRegistration -d 3 \
        --verbose -o [output_,output_warped_image.nii.gz] \
        -t Rigid[0.1] -m Mattes[$template,$bias_cor,1,64,None] \
        -c 1000x500x250x100x50x25 -s 8x4x2x1x0.5x0 -f 6x5x4x3x2x1 --masks [NULL,NULL] \
        -t Affine[0.1] -m Mattes[$template,warped.nii.gz,1,64,None] \
        -c 100x50x25 -s 1x0.5x0 -f 3x2x1 --masks [$mask,$mask] \
        --interpolation BSpline[5]

        '''

        #-d
        reg.inputs.dimension = 3
        #-m
        reg.inputs.metric = ['Mattes']*2
        reg.inputs.metric_weight = [1]*2
        reg.inputs.radius_or_number_of_bins = [64]*2
        #-t
        reg.inputs.transforms = ['Rigid', 'Affine']
        reg.inputs.transform_parameters = [(0.1,), (0.1,)]
        reg.inputs.sampling_strategy = ['None','None']
        #-c
        reg.inputs.number_of_iterations = [[1000, 500, 250, 100, 50, 25], [100, 50, 25]]
        #-s
        reg.inputs.smoothing_sigmas = [[8,4,2,1,0.5,0], [1,0.5,0]]
        #-f
        reg.inputs.shrink_factors = [[6,5,4,3,2,1,0], [3,2,1]]
        #-i 0
        reg.inputs.initialize_transforms_per_stage = False
        #-o
        reg.inputs.output_transform_prefix = "output_"
        #-z 1
        reg.inputs.collapse_output_transforms = True #collapses transforms into a single File
        #interpolation Linear
        reg.inputs.interpolation = 'BSpline'
        reg.inputs.interpolation_parameters = (5,)
        #-u 0
        reg.inputs.use_histogram_matching = False
        #--masks
        reg.inputs.moving_image_masks = ['NULL', anat_mask]
        reg.inputs.fixed_image_masks = ['NULL', anat_mask]

        reg.inputs.moving_image = moving_image
        reg.inputs.fixed_image = fixed_image
        reg.inputs.write_composite_transform=True
        reg.inputs.output_warped_image = 'output_warped_image.nii.gz'

    elif reg_type=='Rigid':
        '''
        antsRegistration -d 3 \
        --verbose -o [output_,output_warped_image.nii.gz] \
        -t Rigid[0.1] -m Mattes[$template,$bias_cor,1,64,None] \
        -c 1000x500x250x100x50x25 -s 8x4x2x1x0.5x0 -f 6x5x4x3x2x1 --masks [NULL,NULL] \
        --interpolation BSpline[5]

        '''
        reg = Registration(dimension=3, metric = ['Mattes'], metric_weight = [1], radius_or_number_of_bins = [64],
                        transforms = ['Rigid'], transform_parameters = [(0.1,)], sampling_strategy = ['None'],
                        number_of_iterations = [[1000, 500, 250, 100, 50, 25]], smoothing_sigmas = [[8,4,2,1,0.5,0]],
                        shrink_factors = [[6,5,4,3,2,1]], initialize_transforms_per_stage = False, output_transform_prefix = "output_",
                        collapse_output_transforms = True, interpolation = 'BSpline', interpolation_parameters = (5,), use_histogram_matching = [False],
                        moving_image_masks = 'NULL', fixed_image_masks = 'NULL', moving_image = moving_image, fixed_image = fixed_image,
                        write_composite_transform=True, output_warped_image = 'output_warped_image.nii.gz')

    return reg
