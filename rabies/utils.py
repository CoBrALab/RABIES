import os
import pathlib  # Better path manipulation
import SimpleITK as sitk
import numpy as np
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, InputMultiPath, BaseInterface
)

######################
#IMAGE MANIPULATION
######################


def recover_3D(mask_file, vector_map):
    mask_img = sitk.ReadImage(mask_file)
    brain_mask = sitk.GetArrayFromImage(mask_img)
    volume_indices=brain_mask.astype(bool)
    volume=np.zeros(brain_mask.shape)
    volume[volume_indices]=vector_map
    volume_img = copyInfo_3DImage(sitk.GetImageFromArray(
        volume, isVector=False), mask_img)
    return volume_img


def recover_4D(mask_file, vector_maps, ref_4d):
    #vector maps of shape num_volumeXnum_voxel
    mask_img = sitk.ReadImage(mask_file)
    brain_mask = sitk.GetArrayFromImage(mask_img)
    volume_indices=brain_mask.astype(bool)
    shape=(vector_maps.shape[0],brain_mask.shape[0],brain_mask.shape[1],brain_mask.shape[2])
    volumes=np.zeros(shape)
    for i in range(vector_maps.shape[0]):
        volume=volumes[i,:,:,:]
        volume[volume_indices]=vector_maps[i,:]
        volumes[i,:,:,:]=volume
    volume_img = copyInfo_4DImage(sitk.GetImageFromArray(
        volumes, isVector=False), mask_img, sitk.ReadImage(ref_4d))
    return volume_img


def resample_image_spacing(image, output_spacing):
    dimension = 3
    identity = sitk.Transform(dimension, sitk.sitkIdentity)

    # Compute grid size based on the physical size and spacing.
    input_size = image.GetSize()
    sampling_ratio = np.asarray(image.GetSpacing())/np.asarray(output_spacing)
    output_size = [int(input_size[0]*sampling_ratio[0]), int(input_size[1]
                                                             * sampling_ratio[1]), int(input_size[2]*sampling_ratio[2])]

    # set default threader to platform to avoid freezing with MultiProc https://github.com/SimpleITK/SimpleITK/issues/1239
    sitk.ProcessObject_SetGlobalDefaultThreader('Platform')
    # the resampling is done with BSpline, the interpolator order is 3
    # other BSpline options are blurry, see https://discourse.itk.org/t/resample-produces-blurry-results-when-just-cropping/4473
    resampled_image = sitk.Resample(image, output_size, identity, sitk.sitkBSpline,
                                    image.GetOrigin(), output_spacing, image.GetDirection())
    # clip potential negative values
    array = sitk.GetArrayFromImage(resampled_image)
    array[(array < 0).astype(bool)] = 0
    pos_resampled_image = sitk.GetImageFromArray(array, isVector=False)
    pos_resampled_image.CopyInformation(resampled_image)
    return pos_resampled_image


def resample_image_spacing_4d(image_4d, output_spacing, clip_negative=True): 
    # output_spacing should be specifying the resolution of the 3 spatial dimensions in mm
    if not image_4d.GetDimension()==4:
        raise ValueError("The image must be 4 dimensional.")
    
    # Compute grid size based on the physical size and spacing.
    input_size = image_4d.GetSize()
    sampling_ratio = np.asarray(image_4d.GetSpacing())[:3]/np.asarray(output_spacing)
    output_size = [int(input_size[0]*sampling_ratio[0]), int(input_size[1]
                                                                * sampling_ratio[1]), int(input_size[2]*sampling_ratio[2])]

    origin = image_4d.GetOrigin()[:3]

    direction_4d = image_4d.GetDirection()
    direction_3d = tuple(direction_4d[:3]+direction_4d[4:7]+direction_4d[8:11])

    dimension = 3
    identity = sitk.Transform(dimension, sitk.sitkIdentity)
    # set default threader to platform to avoid freezing with MultiProc https://github.com/SimpleITK/SimpleITK/issues/1239
    sitk.ProcessObject_SetGlobalDefaultThreader('Platform')
    resampled_list=[]
    for i in range(input_size[3]):
        resampled_image = sitk.Resample(image_4d[:,:,:,i], output_size, identity, sitk.sitkLinear,
                                        origin, output_spacing, direction_3d)
        resampled_list.append(resampled_image)
    combined = sitk.JoinSeries(resampled_list) 
    resampled_4d = copyInfo_4DImage(combined, resampled_image, image_4d)

    if clip_negative:
        # clip potential negative values
        array = sitk.GetArrayFromImage(resampled_4d)
        array[(array < 0).astype(bool)] = 0
        pos_resampled_image = sitk.GetImageFromArray(array, isVector=False)
        pos_resampled_image.CopyInformation(resampled_4d)
        resampled_4d = pos_resampled_image
    return resampled_4d


def copyInfo_4DImage(image_4d, ref_3d, ref_4d):
    # function to establish metadata of an input 4d image. The ref_3d will provide
    # the information for the first 3 dimensions, and the ref_4d for the 4th.
    if ref_3d.GetDimension() == 4:
        image_4d.SetSpacing(
            tuple(list(ref_3d.GetSpacing()[:3])+[ref_4d.GetSpacing()[3]]))
        image_4d.SetOrigin(
            tuple(list(ref_3d.GetOrigin()[:3])+[ref_4d.GetOrigin()[3]]))
        dim_3d = list(ref_3d.GetDirection())
        dim_4d = list(ref_4d.GetDirection())
        image_4d.SetDirection(
            tuple(dim_3d[:3]+[dim_4d[3]]+dim_3d[4:7]+[dim_4d[7]]+dim_3d[8:11]+dim_4d[11:]))
    elif ref_3d.GetDimension() == 3:
        image_4d.SetSpacing(
            tuple(list(ref_3d.GetSpacing())+[ref_4d.GetSpacing()[3]]))
        image_4d.SetOrigin(
            tuple(list(ref_3d.GetOrigin())+[ref_4d.GetOrigin()[3]]))
        dim_3d = list(ref_3d.GetDirection())
        dim_4d = list(ref_4d.GetDirection())
        image_4d.SetDirection(
            tuple(dim_3d[:3]+[dim_4d[3]]+dim_3d[3:6]+[dim_4d[7]]+dim_3d[6:9]+dim_4d[11:]))
    else:
        raise ValueError('Unknown reference image dimensions.')
    return image_4d


def copyInfo_3DImage(image_3d, ref_3d):
    if ref_3d.GetDimension() == 4:
        image_3d.SetSpacing(ref_3d.GetSpacing()[:3])
        image_3d.SetOrigin(ref_3d.GetOrigin()[:3])
        dim_3d = list(ref_3d.GetDirection())
        image_3d.SetDirection(tuple(dim_3d[:3]+dim_3d[4:7]+dim_3d[8:11]))
    elif ref_3d.GetDimension() == 3:
        image_3d.SetSpacing(ref_3d.GetSpacing())
        image_3d.SetOrigin(ref_3d.GetOrigin())
        image_3d.SetDirection(ref_3d.GetDirection())
    else:
        raise ValueError('Unknown reference image dimensions.')
    return image_3d


class slice_applyTransformsInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="Input 4D EPI")
    ref_file = File(exists=True, mandatory=True,
                    desc="The reference 3D space to which the EPI will be warped.")
    transforms = traits.List(desc="List of transforms to apply to every volume.")
    inverses = traits.List(
        desc="Define whether some transforms must be inverse, with a boolean list where true defines inverse e.g.[0,1,0]")
    apply_motcorr = traits.Bool(
        default=True, desc="Whether to apply motion realignment.")
    motcorr_params = File(
        exists=True, desc="xforms from head motion estimation .csv file")
    resampling_dim = traits.Str(
        desc="Specify the image dimension of post-resampling.")
    interpolation = traits.Str(
        desc="Select the interpolator for antsApplyTransform.")
    rabies_data_type = traits.Int(mandatory=True,
                                  desc="Integer specifying SimpleITK data type.")


class slice_applyTransformsOutputSpec(TraitedSpec):
    out_files = traits.List(
        desc="warped images after the application of the transforms")


class slice_applyTransforms(BaseInterface):
    """
    This interface will apply a set of transforms to an input 4D EPI as well as motion realignment if specified.
    Susceptibility distortion correction can be applied through the provided transforms. A list of the corrected
    single volumes will be provided as outputs, and these volumes require to be merged to recover timeseries.
    """

    input_spec = slice_applyTransformsInputSpec
    output_spec = slice_applyTransformsOutputSpec

    def _run_interface(self, runtime):
        # resampling the reference image to the dimension of the EPI

        img = sitk.ReadImage(self.inputs.in_file, self.inputs.rabies_data_type)

        if self.inputs.resampling_dim == 'ref_file': # with 'ref_file' the reference file is unaltered, and thus directly defines the commonspace resolution
            ref_file = self.inputs.ref_file
        else:
            if not self.inputs.resampling_dim == 'inputs_defined':
                shape = self.inputs.resampling_dim.split('x')
                spacing = (float(shape[0]), float(shape[1]), float(shape[2]))
            else:
                spacing = img.GetSpacing()[:3]
            resampled = resample_image_spacing(sitk.ReadImage(
                self.inputs.ref_file, self.inputs.rabies_data_type), spacing)
            ref_file = os.path.abspath('resampled.nii.gz')
            sitk.WriteImage(resampled, ref_file)

        # Splitting bold file into lists of single volumes
        [bold_volumes, num_volumes] = split_volumes(
            self.inputs.in_file, "bold_", self.inputs.rabies_data_type)

        if self.inputs.apply_motcorr:
            motcorr_params = self.inputs.motcorr_params

        warped_volumes = []

        orig_transforms = self.inputs.transforms
        orig_inverses = self.inputs.inverses
        for x in range(0, num_volumes):
            warped_vol_fname = os.path.abspath(
                "deformed_volume" + str(x) + ".nii.gz")
            warped_volumes.append(warped_vol_fname)

            if self.inputs.apply_motcorr:
                command = f'antsMotionCorrStats -m {motcorr_params} -o motcorr_vol{x}.mat -t {x}'
                rc,c_out = run_command(command)

                transforms = orig_transforms+[f'motcorr_vol{x}.mat']
                inverses = orig_inverses+[0]
            else:
                transforms = orig_transforms
                inverses = orig_inverses

            exec_applyTransforms(transforms, inverses, bold_volumes[x], ref_file, warped_vol_fname, interpolation=self.inputs.interpolation)
            # change image to specified data type
            sitk.WriteImage(sitk.ReadImage(warped_vol_fname,
                                           self.inputs.rabies_data_type), warped_vol_fname)

        setattr(self, 'out_files', warped_volumes)
        return runtime

    def _list_outputs(self):
        return {'out_files': getattr(self, 'out_files')}


def exec_applyTransforms(transforms, inverses, input_image, ref_image, output_image, interpolation):
    # tranforms is a list of transform files, set in order of call within antsApplyTransforms
    transform_string = ""
    for transform, inverse in zip(transforms, inverses):
        if transform=='NULL':
            continue
        elif bool(inverse):
            transform_string += f"-t [{transform},1] "
        else:
            transform_string += f"-t {transform} "

    command = f'antsApplyTransforms -i {input_image} {transform_string}-n {interpolation} -r {ref_image} -o {output_image}'
    rc,c_out = run_command(command)
    if not os.path.isfile(output_image):
        raise ValueError(
            "Missing output image. Transform call failed: "+command)


def split_volumes(in_file, output_prefix, rabies_data_type):
    '''
    Takes as input a 4D .nii file and splits it into separate time series
    volumes by splitting on the 4th dimension
    '''
    in_nii = sitk.ReadImage(in_file, rabies_data_type)
    num_dimensions = len(in_nii.GetSize())
    num_volumes = in_nii.GetSize()[3]

    if num_dimensions != 4:
        raise ValueError("the input file must be of dimensions 4")

    volumes = []
    for x in range(0, num_volumes):
        data_slice = sitk.GetArrayFromImage(in_nii)[x, :, :, :]
        slice_fname = os.path.abspath(
            output_prefix + "vol" + str(x) + ".nii.gz")
        image_3d = copyInfo_3DImage(sitk.GetImageFromArray(
            data_slice, isVector=False), in_nii)
        sitk.WriteImage(image_3d, slice_fname)
        volumes.append(slice_fname)

    return [volumes, num_volumes]


class MergeInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiPath(File(exists=True), mandatory=True,
                              desc='input list of files to merge, listed in the order to merge')
    header_source = File(exists=True, mandatory=True,
                         desc='a Nifti file from which the header should be copied')
    clip_negative = traits.Bool(
        desc="Whether to clip out negative values.")
    rabies_data_type = traits.Int(mandatory=True,
                                  desc="Integer specifying SimpleITK data type.")


class MergeOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='output merged file')


class Merge(BaseInterface):
    """
    Takes a list of 3D Nifti files and merge them in the order listed.
    """

    input_spec = MergeInputSpec
    output_spec = MergeOutputSpec

    def _run_interface(self, runtime):
        filename_split = pathlib.Path(
            self.inputs.header_source).name.rsplit(".nii")

        sample_volume = sitk.ReadImage(
            self.inputs.in_files[0], self.inputs.rabies_data_type)
        length = len(self.inputs.in_files)
        shape = sitk.GetArrayFromImage(sample_volume).shape
        combined = np.zeros((length, shape[0], shape[1], shape[2]))

        i = 0
        for file in self.inputs.in_files:
            combined[i, :, :, :] = sitk.GetArrayFromImage(
                sitk.ReadImage(file, self.inputs.rabies_data_type))[:, :, :]
            i = i+1
        if (i != length):
            raise ValueError("Error occured with Merge.")
        combined_files = os.path.abspath(
            f"{filename_split[0]}_combined.nii.gz")

        if self.inputs.clip_negative:
            # clip potential negative values
            combined[(combined < 0).astype(bool)] = 0
        combined_image = sitk.GetImageFromArray(combined, isVector=False)

        # set metadata and affine for the newly constructed 4D image
        header_source = sitk.ReadImage(
            self.inputs.header_source, self.inputs.rabies_data_type)
        combined_image = copyInfo_4DImage(
            combined_image, sample_volume, header_source)
        sitk.WriteImage(combined_image, combined_files)

        setattr(self, 'out_file', combined_files)
        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file')}


######################
#GENERAL
######################

def run_command(command, verbose = False):
    # Run command and collect stdout
    # http://blog.endpoint.com/2015/01/getting-realtime-output-using-python.html # noqa
    from nipype import logging
    log = logging.getLogger('nipype.workflow')
    log.debug('Running: '+command)

    import subprocess
    try:
        process = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=True,
            shell=True,
            )
    except Exception as e:
        log.warning(e.output.decode("utf-8"))
        raise

    c_out = process.stdout.decode("utf-8")
    if not c_out == '':
        if verbose:
            log.info(c_out)
        else:
            log.debug(c_out)
    if process.stderr is not None:
        if verbose:
            log.info(process.stderr)
        else:
            log.warning(process.stderr)
    rc = process.returncode
    if len(c_out)>0 and c_out[-1]=='\n': # remove the extra break point, can affect string construction
        c_out = c_out[:-1] 
    return rc,c_out


def flatten_list(l):
    if type(l) == list:
        flattened = []
        for e in l:
            e = flatten_list(e)
            if type(e) == list:
                flattened += e
            else:
                flattened += [e]
        return flattened
    else:
        return l


def filter_scan_exclusion(exclusion_list, split_name):
    # the function removes a list of scan IDs from split_name
    
    # exclusion_list: the input provided by the user
    # split_name: a list of all scan IDs that were found
    
    import numpy as np
    import pandas as pd
    if os.path.isfile(os.path.abspath(exclusion_list[0])):
        updated_split_name=[]
        if not '.nii' in pathlib.Path(exclusion_list[0]).name:
            # read the file as a .txt
            exclusion_list = np.array(pd.read_csv(os.path.abspath(exclusion_list[0]), header=None)).flatten()
        for split in split_name:
            exclude = False
            for scan in exclusion_list:
                if split in scan:
                    exclude = True
            if not exclude:
                updated_split_name.append(split)
    elif exclusion_list[0]=='none':
        updated_split_name = split_name
    else:
        raise ValueError(f"The --exclusion_ids {exclusion_list} input had improper format. It must the full path to a .txt or .nii files.")
    
    if len(updated_split_name)==0:
        raise ValueError(f"""
            No scans are left after scan exclusion!
            """)

    return updated_split_name


def filter_scan_inclusion(inclusion_list, split_name):
    # the function will update the list of scan IDs in split_name to correspond to inclusion/exclusion list
    
    # inclusion_list: the input provided by the user
    # split_name: a list of all scan IDs that were found
    
    import numpy as np
    import pandas as pd
    if os.path.isfile(os.path.abspath(inclusion_list[0])):
        updated_split_name=[]
        if '.nii' in pathlib.Path(inclusion_list[0]).name:
            for scan in inclusion_list:
                updated_split_name.append(find_split(scan, split_name))
        else:
            # read the file as a .txt
            inclusion_list = np.array(pd.read_csv(os.path.abspath(inclusion_list[0]), header=None)).flatten()
            for scan in inclusion_list:
                updated_split_name.append(find_split(scan, split_name))
    elif inclusion_list[0]=='all':
        updated_split_name = split_name
    else:
        raise ValueError(f"The --inclusion_ids {inclusion_list} input had improper format. It must the full path to a .txt or .nii files, or 'all' to keep all scans.")
    return updated_split_name


def find_split(scan, split_name):
    for split in split_name:
        if split in scan:
            return split
    raise ValueError(f"No previous file name is matching {scan}")


######################
#FUNCTIONS TO READ WORKFLOW GRAPH
######################

def fill_split_dict(d, output_bold, split_name, split_dict, keys, node_dict, match_targets):
    if isinstance(d, dict):
        for key in list(d.keys()):
            fill_split_dict(d[key], output_bold, split_name, split_dict, keys+[key], node_dict, match_targets)
    else:
        f = d.result.outputs.get()[output_bold]
        split = pathlib.Path(f).name.rsplit(".nii")[0]
        split_name.append(split)
        split_dict[split]={}
        target_list = list(match_targets.keys())
        for target in target_list:
            [unit, output] = match_targets[target]
            node = retrieve_node(node_dict[unit], keys)
            split_dict[split][target] = node.result.outputs.get()[output]
        
def retrieve_node(d, keys):
    if isinstance(d, dict):
        return retrieve_node(d[keys[0]], keys[1:])
    else:
        return d


def get_workflow_dict(workflow_file):
    import pickle
    with open(workflow_file, 'rb') as handle:
        graph = pickle.load(handle)
    
    node_list = list(graph.nodes)
    node_dict = {}
    for node in node_list:
        key_l = [node.fullname]+node.parameterization
        fill_node_dict(node_dict, key_l, node)
    return node_dict


def fill_node_dict(d, key_l, e):
    if len(key_l)>0:
        key = key_l[0]
        if not (key in list(d.keys())):
            d[key] = {}
        d[key] = fill_node_dict(d[key], key_l[1:], e)
        return d
    else:
        return e


######################
#DEBUGGING
######################

def generate_token_data(tmppath, number_scans):
    # this function generates fake scans at low resolution for quick testing and debugging

    os.makedirs(tmppath+'/inputs', exist_ok=True)

    if 'XDG_DATA_HOME' in os.environ.keys():
        rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
    else:
        rabies_path = os.environ['HOME']+'/.local/share/rabies'

    template = f"{rabies_path}/DSURQE_40micron_average.nii.gz"
    mask = f"{rabies_path}/DSURQE_40micron_mask.nii.gz"
    melodic_file = f"{rabies_path}/melodic_IC.nii.gz"

    spacing = (float(1), float(1), float(1))  # resample to 1mmx1mmx1mm
    resampled_template = resample_image_spacing(sitk.ReadImage(template), spacing)
    # generate template masks
    resampled_mask = resample_image_spacing(sitk.ReadImage(mask), spacing)
    array = sitk.GetArrayFromImage(resampled_mask)
    array[array < 1] = 0
    array[array > 1] = 1
    binarized = sitk.GetImageFromArray(array, isVector=False)
    binarized.CopyInformation(resampled_mask)
    sitk.WriteImage(binarized, tmppath+'/inputs/token_mask.nii.gz')
    array[:, :, :6] = 0
    binarized = sitk.GetImageFromArray(array, isVector=False)
    binarized.CopyInformation(resampled_mask)
    sitk.WriteImage(binarized, tmppath+'/inputs/token_mask_half.nii.gz')

    # generate fake scans from the template
    array = sitk.GetArrayFromImage(resampled_template)
    array_4d = np.repeat(array[np.newaxis, :, :, :], 15, axis=0)
    
    melodic_img = sitk.ReadImage(melodic_file)
    network1_map = sitk.GetArrayFromImage(sitk.Resample(melodic_img[:,:,:,5], resampled_template))
    network2_map = sitk.GetArrayFromImage(sitk.Resample(melodic_img[:,:,:,19], resampled_template))
    time1 = np.random.normal(0, array_4d.mean()/100, array_4d.shape[0]) # network timecourse; scale is 1% of image intensity
    time2 = np.random.normal(0, array_4d.mean()/100, array_4d.shape[0]) # network timecourse; scale is 1% of image intensity
    # creating fake network timeseries
    network1_time = (np.repeat(network1_map[np.newaxis, :, :, :], 15, axis=0).T*time1).T
    network2_time = (np.repeat(network2_map[np.newaxis, :, :, :], 15, axis=0).T*time2).T

    for i in range(number_scans):
        # generate anatomical scan
        sitk.WriteImage(resampled_template, tmppath+f'/inputs/sub-token{i+1}_T1w.nii.gz')
        # generate functional scan
        array_4d_ = array_4d + network1_time + network2_time + np.random.normal(0, array_4d.mean()/ 100, array_4d.shape)  # add gaussian noise; scale is 1% of the mean intensity of the template
        sitk.WriteImage(sitk.GetImageFromArray(array_4d_, isVector=False),
                        tmppath+f'/inputs/sub-token{i+1}_bold.nii.gz')

        # necessary to read matrix orientation properly at the analysis stage
        sitk.WriteImage(copyInfo_4DImage(sitk.ReadImage(tmppath+f'/inputs/sub-token{i+1}_bold.nii.gz'), sitk.ReadImage(tmppath
                        + f'/inputs/sub-token{i+1}_T1w.nii.gz'), sitk.ReadImage(tmppath+f'/inputs/sub-token{i+1}_bold.nii.gz')), tmppath+f'/inputs/sub-token{i+1}_bold.nii.gz')
