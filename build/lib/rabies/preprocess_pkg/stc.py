from nipype.pipeline import engine as pe
from nipype.interfaces.utility import Function
from nipype.interfaces import utility as niu


def init_bold_stc_wf(opts, name='bold_stc_wf'):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['stc_file']), name='outputnode')

    if opts.apply_STC:
        slice_timing_correction_node = pe.Node(Function(input_names=['in_file', 'tr', 'tpattern', 'stc_axis', 'interp_method', 'rabies_data_type'],
                                                        output_names=[
                                                            'out_file'],
                                                        function=slice_timing_correction),
                                               name='slice_timing_correction', mem_gb=1.5*opts.scale_min_memory)
        slice_timing_correction_node.inputs.tr = opts.TR
        slice_timing_correction_node.inputs.tpattern = opts.tpattern
        slice_timing_correction_node.inputs.stc_axis = opts.stc_axis
        slice_timing_correction_node.inputs.interp_method = opts.interp_method
        slice_timing_correction_node.inputs.rabies_data_type = opts.data_type
        slice_timing_correction_node.plugin_args = {
            'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}

        workflow.connect([
            (inputnode, slice_timing_correction_node, [('bold_file', 'in_file')]),
            (slice_timing_correction_node,
             outputnode, [('out_file', 'stc_file')]),
        ])
    else:
        workflow.connect([
            (inputnode, outputnode, [('bold_file', 'stc_file')]),
        ])

    return workflow


def slice_timing_correction(in_file, tr='auto', tpattern='alt-z', stc_axis='Y', interp_method = 'fourier', rabies_data_type=8):
    '''
    This functions applies slice-timing correction on the anterior-posterior
    slice acquisition direction. The input image, assumed to be in RAS orientation
    (accoring to nibabel; note that the nibabel reading of RAS corresponds to
     LPI for AFNI). The A and S dimensions will be swapped to apply AFNI's
    3dTshift STC with a quintic fitting function, which can only be applied to
    the Z dimension of the data matrix. The corrected image is then re-created with
    proper axes and the corrected timeseries.

    **Inputs**

        in_file
            BOLD series NIfTI file in RAS orientation.
        tr
            TR of the BOLD image.
        tpattern
            Input to AFNI's 3dTshift -tpattern option, which specifies the
            directionality of slice acquisition, or whether it is sequential or
            interleaved.
        stc_axis
            Can specify over which axis between X,Y and Z slices were acquired
        interp_method
            Input to AFNI's 3dTshift that specifies the interpolation method to be used.

    **Outputs**

        out_file
            Slice-timing corrected BOLD series NIfTI file

    '''

    import os
    import SimpleITK as sitk
    import numpy as np

    img = sitk.ReadImage(in_file, rabies_data_type)

    if tr=='auto':
        tr = str(img.GetSpacing()[3])+'s'
    else:
        tr = str(tr)+'s'

    # get image data
    img_array = sitk.GetArrayFromImage(img)
    shape = img_array.shape
    # 3dTshift applies STC on the Z axis, so swap the desired axis in position before applying STC
    if stc_axis=='Z':
        # no need to swap axes if it is Z
        target_file = in_file
    else:
        target_file = os.path.abspath('STC_swap.nii.gz')
        if stc_axis=='Y':
            new_array = img_array.transpose(0,2,1,3)
        elif stc_axis=='X':
            new_array = img_array.transpose(0,3,2,1)
        else:
            raise ValueError('Wrong axis name.')

        image_out = sitk.GetImageFromArray(new_array, isVector=False)
        sitk.WriteImage(image_out, target_file)

    command = f'3dTshift -{interp_method} -prefix temp_tshift.nii.gz -tpattern {tpattern} -TR {tr} {target_file}'
    from rabies.utils import run_command
    rc,c_out = run_command(command)

    tshift_img = sitk.ReadImage(
        'temp_tshift.nii.gz', rabies_data_type)
    tshift_array = sitk.GetArrayFromImage(tshift_img)

    new_array = np.zeros(shape)
    if stc_axis=='Z':
        new_array = tshift_array
    elif stc_axis=='Y':
        new_array = tshift_array.transpose(0,2,1,3)
    elif stc_axis=='X':
        new_array = tshift_array.transpose(0,3,2,1)
    image_out = sitk.GetImageFromArray(new_array, isVector=False)

    from rabies.utils import copyInfo_4DImage
    image_out = copyInfo_4DImage(image_out, img, img)

    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(in_file).name.rsplit(".nii")
    out_file = os.path.abspath(filename_split[0]+'_tshift.nii.gz')
    sitk.WriteImage(image_out, out_file)
    return out_file
