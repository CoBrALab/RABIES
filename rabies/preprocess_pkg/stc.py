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
        slice_timing_correction_node = pe.Node(Function(input_names=['in_file', 'args_3dTshift', 'stc_axis', 'rabies_data_type'],
                                                        output_names=[
                                                            'out_file'],
                                                        function=slice_timing_correction),
                                               name='slice_timing_correction', mem_gb=1.5*opts.scale_min_memory)
        slice_timing_correction_node.inputs.args_3dTshift = opts.args_3dTshift
        slice_timing_correction_node.inputs.stc_axis = opts.stc_axis
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


def slice_timing_correction(in_file, args_3dTshift = '-Fourier -tpattern alt-z', stc_axis='Y', rabies_data_type=8):
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
        args_3dTshift
            Set of custom parameters to input to 3dTshift
        stc_axis
            Can specify over which axis between X,Y and Z slices were acquired

    **Outputs**

        out_file
            Slice-timing corrected BOLD series NIfTI file

    '''

    import os
    import SimpleITK as sitk
    import numpy as np

    from rabies.utils import get_sitk_header
    img_header = get_sitk_header(in_file)

    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(in_file).name.rsplit(".nii")
    out_file = os.path.abspath(filename_split[0]+'_tshift.nii.gz')

    '''
    3dTshift applies STC on the Z axis, so swap the desired axis in position before applying STC
    '''
    if stc_axis=='Z':
        # no need to swap axes if it is Z
        target_file = in_file
    else:
        img_array = sitk.GetArrayFromImage(sitk.ReadImage(in_file, rabies_data_type))
        target_file = os.path.abspath('STC_swap.nii.gz')
        if stc_axis=='Y':
            new_array = img_array.transpose(0,2,1,3)
        elif stc_axis=='X':
            new_array = img_array.transpose(0,3,2,1)
        else:
            raise ValueError('Wrong axis name.')
        image_out = sitk.GetImageFromArray(new_array, isVector=False)
        sitk.WriteImage(image_out, target_file)
        del img_array, image_out, new_array

    '''
    Run 3dTshift
    '''
    if not '-TR' in args_3dTshift: # read the header if not provided explicitely
        TR = img_header.GetSpacing()[3]
        args_3dTshift += f" -TR {TR}s"
    if '-prefix' in args_3dTshift:
        raise ValueError("Do not input -prefix to 3dTshift - RABIES uses it internally.")
    
    command = f'3dTshift -prefix {out_file} {args_3dTshift} {target_file}'
    from rabies.utils import run_command
    rc,c_out = run_command(command)

    '''
    Re-transpose back the axes if they were shifted
    '''
    if not stc_axis=='Z':
        tshift_array = sitk.GetArrayFromImage(sitk.ReadImage(out_file, rabies_data_type))
        if stc_axis=='Y':
            new_array = tshift_array.transpose(0,2,1,3)
        elif stc_axis=='X':
            new_array = tshift_array.transpose(0,3,2,1)
        from rabies.utils import copyInfo_4DImage
        image_out = copyInfo_4DImage(sitk.GetImageFromArray(new_array, isVector=False), img_header, img_header)
        sitk.WriteImage(image_out, out_file)

    return out_file
