import os
import pathlib  # Better path manipulation
import SimpleITK as sitk
import numpy as np
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

######################
#IMAGE MANIPULATION
######################


def get_sitk_header(file_path): # read a nifti header without loading the data
    reader = sitk.ImageFileReader()
    reader.SetFileName(file_path)
    reader.ReadImageInformation()
    return reader


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


class ResampleVolumesInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="Input 3D or 4D file to resample")
    ref_file = File(exists=True, mandatory=True,
                    desc="The reference 3D space to which the EPI will be warped.")
    transforms_3d_files = traits.List(desc="List of transform files to apply to every volume.")
    inverses_3d = traits.List(
        desc="Define whether some transforms must be inverse, with a boolean list where true defines inverse e.g.[0,1,0]")
    apply_motcorr = traits.Bool(
        default=True, desc="Whether to apply motion realignment, only for 4D file.")
    motcorr_params = File(
        exists=True, desc="xforms from head motion estimation .csv file")
    resampling_dim = traits.Str(
        desc="Specify the image dimension of post-resampling.")
    interpolation = traits.Int(
        desc="An SITK interpolator for resampling.")
    name_source = File(exists=True, mandatory=True,
                         desc='a Nifti file from which the header should be copied')
    clip_negative = traits.Bool(
        desc="Whether to clip out negative values after resampling.")
    rabies_data_type = traits.Int(mandatory=True,
                                  desc="Integer specifying SimpleITK data type.")


class ResampleVolumesOutputSpec(TraitedSpec):
    resampled_file = File(exists=True, desc='output 4D resampled file')


class ResampleVolumes(BaseInterface):
    """
    This interface will apply a set of transforms to an input 4D EPI as well as motion realignment if specified.
    Susceptibility distortion correction can be applied through the provided transforms. A list of the corrected
    single volumes will be provided as outputs, and these volumes require to be merged to recover timeseries.
    """

    input_spec = ResampleVolumesInputSpec
    output_spec = ResampleVolumesOutputSpec

    def _run_interface(self, runtime):
        # set default threader to platform to avoid freezing with MultiProc https://github.com/SimpleITK/SimpleITK/issues/1239
        sitk.ProcessObject_SetGlobalDefaultThreader('Platform')        
        # avoid crashes if receiving an empty file
        if pathlib.Path(self.inputs.in_file).name=='empty.nii.gz':
            setattr(self, 'resampled_file', self.inputs.in_file)
            return runtime

        # prepare output name
        filename_split = pathlib.Path(
            self.inputs.name_source).name.rsplit(".nii")
        resampled_file = os.path.abspath(
            f"{filename_split[0]}_resampled.nii.gz")

        ref_img = sitk.ReadImage(self.inputs.ref_file, self.inputs.rabies_data_type)

        # preparing the resampling dimensions of the reference space
        if not self.inputs.resampling_dim == 'ref_file': # with 'ref_file' the reference file is unaltered, and thus directly defines the commonspace resolution
            if self.inputs.resampling_dim == 'inputs_defined':
                spacing = get_sitk_header(self.inputs.in_file).GetSpacing()[:3]
            else:
                shape = self.inputs.resampling_dim.split('x')
                spacing = (float(shape[0]), float(shape[1]), float(shape[2]))
            ref_img = resample_image_spacing(ref_img, spacing)

        # replace undefined inputs for an empty list to avoid crashes
        from nipype.interfaces.base import isdefined
        transforms_3d_files = self.inputs.transforms_3d_files if isdefined(self.inputs.transforms_3d_files) else []
        inverses_3d = self.inputs.inverses_3d if isdefined(self.inputs.inverses_3d) else []
        n_procs = int(os.environ['RABIES_ITK_NUM_THREADS']) if "RABIES_ITK_NUM_THREADS" in os.environ else os.cpu_count() # default to number of CPUs

        motcorr_params_file = self.inputs.motcorr_params if self.inputs.apply_motcorr else None

        resampled_img = resample_volumes(self.inputs.in_file, ref_img, transforms_3d_files, inverses_3d, 
                                         motcorr_params_file = motcorr_params_file, interpolation=self.inputs.interpolation, 
                                         rabies_data_type=self.inputs.rabies_data_type, clip_negative=self.inputs.clip_negative, 
                                         n_procs=n_procs)
        sitk.WriteImage(resampled_img, resampled_file)

        setattr(self, 'resampled_file', resampled_file)
        return runtime

    def _list_outputs(self):
        return {'resampled_file': getattr(self, 'resampled_file')}


def resample_volumes(in_img, in_ref, transforms_3d_files = [], inverses_3d = [], motcorr_params_file = None, interpolation=sitk.sitkLinear, rabies_data_type=8, clip_negative=False, n_procs=os.cpu_count()):
    import SimpleITK as sitk
    from simpleitk_timeseries_motion_correction.apply_transforms import read_transforms_from_csv, resample_volume

    # the input can be either a nifti file or an SITK image
    if isinstance(in_img, sitk.Image):
        orig_img = in_img
    elif os.path.isfile(in_img):
        orig_img = sitk.ReadImage(in_img, rabies_data_type)
    if isinstance(in_ref, sitk.Image):
        ref_img = in_ref
    elif os.path.isfile(in_ref):
        ref_img = sitk.ReadImage(in_ref, rabies_data_type)

    transforms_3d_l = load_sitk_transforms(transforms_3d_files, inverses_3d)

    if orig_img.GetDimension()==3:
        composite = sitk.CompositeTransform(3)
        for transform in transforms_3d_l:
            composite.AddTransform(transform)
        resampled_img = resample_volume(orig_img, ref_img, composite, interpolation=interpolation, clip_negative=clip_negative, extrapolator=False)

    elif orig_img.GetDimension()==4:
        import concurrent.futures
        from tqdm import tqdm
        if motcorr_params_file is not None:
            hmc_transforms_l = read_transforms_from_csv(csv_file=motcorr_params_file)

        # create a resampling function that can create the composite transform on a per iteration basis to avoid memory overloads
        def resample_volume_composite(volume, reference, composite_transform_list, interpolation, clip_negative=False, extrapolator=False):
            composite = sitk.CompositeTransform(3)
            for transform in composite_transform_list:
                composite.AddTransform(transform)        
            return resample_volume(volume, reference, composite, interpolation, clip_negative, extrapolator)

        num_volumes = orig_img.GetSize()[3]
        # Convert to a list of 3D volumes from input to process in parallel
        input_volumes = []
        size_4d = orig_img.GetSize()
        for i in range(num_volumes):
            extractor = sitk.ExtractImageFilter()
            extractor.SetSize([size_4d[0], size_4d[1], size_4d[2], 0])
            extractor.SetIndex([0, 0, 0, i])
            input_volumes.append(extractor.Execute(orig_img))
        del orig_img, extractor # free memory

        resampled_volumes = [None] * num_volumes

        # Parallel processing
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_procs) as executor:    
            with tqdm(total=num_volumes) as pbar:
                futures = {}
                for i in range(num_volumes):
                    if motcorr_params_file is not None:
                        composite_transform_list = [hmc_transforms_l[i]] + transforms_3d_l
                    else:
                        composite_transform_list = transforms_3d_l
                        
                    future = executor.submit(
                        resample_volume_composite,
                        volume=input_volumes[i],
                        reference=ref_img,
                        composite_transform_list=composite_transform_list,
                        interpolation=interpolation,
                        clip_negative=clip_negative,
                        extrapolator=False,
                    )
                    futures[future] = i
                
                for future in concurrent.futures.as_completed(futures):
                    i = futures[future]
                    try:
                        resampled_volumes[i] = future.result()
                    except Exception as e:
                        print(f"Error processing volume {i}: {e}")
                        # Might want to abort or insert blank?
                        # For now re-raise to fail fast
                        raise e
                    pbar.update(1)
        resampled_img = sitk.JoinSeries(resampled_volumes)

    return resampled_img


def load_sitk_transforms(transforms_files, inverses):
    # input transforms are listed from last to first applied as per antsApplyTransforms convention; this need to be inverted since SITK applies from first to last
    transforms_files = transforms_files[::-1]
    inverses = inverses[::-1]

    # load as a list of SITK transforms
    transforms_l = []
    for transform_file,inverse in zip(transforms_files,inverses):
        if transform_file[-4:]=='.mat':
            transform = sitk.ReadTransform(transform_file)
            if inverse==1:
                transform = transform.GetInverse()
        elif transform_file[-7:]=='.nii.gz':
            # Read the displacement field image and convert to a transform
            transform = sitk.DisplacementFieldTransform(sitk.ReadImage(
                transform_file,
                sitk.sitkVectorFloat64
            ))
        elif transform_file=='NULL':
            continue # don't append if NULL, which is a possible RABIES input
        else:
            raise ValueError(f"Transform files must end with .mat or .nii.gz. Cannot read input {transform_file}.")
        transforms_l.append(transform)
    return transforms_l


class ResampleMaskInputSpec(BaseInterfaceInputSpec):
    mask_file = File(exists=True, mandatory=True,
                desc="Input mask.")
    ref_file = File(exists=True, mandatory=True,
                   desc="Target file defining resampling space.")
    transforms = traits.List(desc="List of transforms to resample the mask.")
    inverses = traits.List(
        desc="Define whether some transforms must be inverse, with a boolean list where true defines inverse e.g.[0,1,0]")
    name_suffix = traits.Str(desc="Suffix added at the output file.")
    name_source = File(exists=True, mandatory=True,
                       desc='Reference file for prefix of the output file.')


class ResampleMaskOutputSpec(TraitedSpec):
    resampled_file = traits.File(desc="The resampled mask file.")


class ResampleMask(BaseInterface):

    input_spec = ResampleMaskInputSpec
    output_spec = ResampleMaskOutputSpec

    def _run_interface(self, runtime):
        import os
        from rabies.utils import antsApplyTransforms
        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.name_source).name.rsplit(".nii")

        if self.inputs.name_suffix is None:
            new_mask_path = os.path.abspath(
                f'{filename_split[0]}_resampled.nii.gz')
        else:
            new_mask_path = os.path.abspath(f'{filename_split[0]}_{self.inputs.name_suffix}.nii.gz')

        # replace undefined inputs for an empty list to avoid crashes
        from nipype.interfaces.base import isdefined
        transforms = self.inputs.transforms if isdefined(self.inputs.transforms) else []
        inverses = self.inputs.inverses if isdefined(self.inputs.inverses) else []

        antsApplyTransforms(transforms = transforms, inverses = inverses, 
                        input_image = self.inputs.mask_file, ref_image = self.inputs.ref_file, output_filename = new_mask_path, interpolation='GenericLabel', rabies_data_type=sitk.sitkInt16, clip_negative=False)

        setattr(self, 'resampled_file', new_mask_path)
        return runtime

    def _list_outputs(self):
        return {'resampled_file': getattr(self, 'resampled_file')}


def antsApplyTransforms(transforms, inverses, input_image, ref_image, output_filename, interpolation, rabies_data_type=8, clip_negative=False):
    # tranforms is a list of transform files, set in order of call within antsApplyTransforms
    transform_string = ""
    for transform, inverse in zip(transforms, inverses):
        if transform=='NULL':
            continue
        elif bool(inverse):
            transform_string += f"-t [{transform},1] "
        else:
            transform_string += f"-t {transform} "

    command = f'antsApplyTransforms -i {input_image} {transform_string}-n {interpolation} -r {ref_image} -o {output_filename} -v'
    rc,c_out = run_command(command)

    if clip_negative or rabies_data_type is not None: # rabies_data_type can be set to None to avoid re-writing the image
        if rabies_data_type is not None:
            # need to reload/save the image to switch to rabies_data_type
            resampled_img = sitk.ReadImage(output_filename, rabies_data_type)
        else:
            resampled_img = sitk.ReadImage(output_filename)
        if clip_negative:
            # clip potential negative values
            array = sitk.GetArrayFromImage(resampled_img)
            array[(array < 0).astype(bool)] = 0
            pos_resampled_image = sitk.GetImageFromArray(array, isVector=False)
            pos_resampled_image.CopyInformation(resampled_img)
            resampled_img = pos_resampled_image
        sitk.WriteImage(resampled_img, output_filename)
    if not os.path.isfile(output_filename):
        raise ValueError(
            "Missing output image. Transform call failed: "+command)


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
        if os.environ.get('RABIES_ITK_THREADS_STATE') == 'ON':
            itk_threads = os.environ.get('RABIES_ITK_NUM_THREADS')
            if itk_threads is None:
                raise ValueError("RABIES_ITK_NUM_THREADS must be set when enabling ITK threading control.")
            os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = itk_threads        
        process = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=True,
            shell=True,
            env=os.environ,
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
            if match_targets[target] is None: # certain target may not have been generated, and needs to be attributed None
                split_dict[split][target] = None
            else:
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
    # create a new melodic with just 2 networks for low-dimensional dual regression
    melodic_networks = sitk.JoinSeries([melodic_img[:,:,:,5],melodic_img[:,:,:,19]]) 
    sitk.WriteImage(melodic_networks, tmppath+'/inputs/melodic_networks.nii.gz')
    
    dim = melodic_img.GetSize()[-1]
    for i in range(dim):
        IC_map = sitk.GetArrayFromImage(sitk.Resample(melodic_img[:,:,:,i], resampled_template))
        time = np.random.normal(0, array_4d.mean()/100, array_4d.shape[0]) # network timecourse; scale is 1% of image intensity
        array_4d+=(np.repeat(IC_map[np.newaxis, :, :, :], 15, axis=0).T*time).T

    for i in range(number_scans):
        # generate anatomical scan
        sitk.WriteImage(resampled_template, tmppath+f'/inputs/sub-token{i+1}_T1w.nii.gz')
        # generate functional scan
        array_4d_ = array_4d + np.random.normal(0, array_4d.mean()/ 100, array_4d.shape)  # add gaussian noise; scale is 1% of the mean intensity of the template
        sitk.WriteImage(sitk.GetImageFromArray(array_4d_, isVector=False),
                        tmppath+f'/inputs/sub-token{i+1}_bold.nii.gz')

        # necessary to read matrix orientation properly at the analysis stage
        sitk.WriteImage(copyInfo_4DImage(sitk.ReadImage(tmppath+f'/inputs/sub-token{i+1}_bold.nii.gz'), sitk.ReadImage(tmppath
                        + f'/inputs/sub-token{i+1}_T1w.nii.gz'), sitk.ReadImage(tmppath+f'/inputs/sub-token{i+1}_bold.nii.gz')), tmppath+f'/inputs/sub-token{i+1}_bold.nii.gz')
