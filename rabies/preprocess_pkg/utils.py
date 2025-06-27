from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

def prep_bids_iter(layout, bids_filter, bold_only=False, inclusion_list=['all'], exclusion_list=['none']):
    '''
    This function takes as input a BIDSLayout, and generates iteration lists
    for managing the workflow's iterables depending on whether --bold_only is
    selected or not.
    '''
    import pathlib

    scan_info = []
    split_name = []
    structural_scan_list = []
    run_iter = {}
    bold_scan_list = []

    subject_list = layout.get_subject()
    if len(subject_list) == 0:
        raise ValueError(
            "No subject information could be retrieved from the BIDS directory. The 'sub-' specification is mandatory.")

    if not 'subject' in list(bids_filter['func'].keys()):
        # enforce that only files with subject ID are read
        bids_filter['func']['subject']=subject_list

    # create the list for all functional images; this is applying all filters from bids_filter
    bold_bids = layout.get(extension=['nii', 'nii.gz'], **bids_filter['func'])
    if len(bold_bids) == 0:
        raise ValueError(
            f"No functional file were found respecting the functional BIDS spec: {bids_filter['func']}")
    
    # remove subject, session and run; these are used later to target single files
    bids_filter['func'].pop('subject', None)
    bids_filter['func'].pop('session', None)
    bids_filter['func'].pop('run', None)

    # filter inclusion/exclusion lists
    from rabies.utils import filter_scan_inclusion, filter_scan_exclusion
    boldname_list=[pathlib.Path(bold.filename).name.rsplit(".nii")[0] for bold in bold_bids]
    updated_split_name = filter_scan_inclusion(inclusion_list, boldname_list)
    updated_split_name = filter_scan_exclusion(exclusion_list, updated_split_name)
    
    filtered_bold_bids=[]
    for name in updated_split_name:
        for bold in bold_bids:
            if name in bold.filename:
                filtered_bold_bids.append(bold)
    bold_bids = filtered_bold_bids

    bold_dict = {}
    for bold in bold_bids:
        sub = bold.get_entities()['subject']
        try:
            ses = bold.get_entities()['session']
        except:
            ses = None

        try:
            run = bold.get_entities()['run']
        except:
            run = None

        if sub not in list(bold_dict.keys()):
            bold_dict[sub] = {}
        if ses not in list(bold_dict[sub].keys()):
            bold_dict[sub][ses] = {}

        bold_list = layout.get(subject=sub, session=ses, run=run, 
                               extension=['nii', 'nii.gz'], return_type='filename',**bids_filter['func'])
        bold_dict[sub][ses][run] = bold_list

    # if not bold_only, then the bold_list and run_iter will be a dictionary with keys being the anat filename
    # otherwise, it will be a list of bold scans themselves
    for sub in list(bold_dict.keys()):
        for ses in list(bold_dict[sub].keys()):
            if not bold_only:
                anat_list = layout.get(subject=sub, session=ses,
                                       extension=['nii', 'nii.gz'], return_type='filename',**bids_filter['anat'])
                if len(anat_list) == 0:
                    raise ValueError(
                        f'Missing an anatomical image for sub {sub} and ses- {ses}, and the following BIDS specs: {bids_filter["anat"]}')
                if len(anat_list) > 1:
                    raise ValueError(
                        f'Duplicate was found for the anatomical file sub- {sub}, ses- {ses}: {str(anat_list)}')
                file = anat_list[0]
                structural_scan_list.append(file)
                filename_template = pathlib.Path(file).name.rsplit(".nii")[0]
                split_name.append(filename_template)
                scan_info.append({'subject_id': sub, 'session': ses})
                run_iter[filename_template] = []

            for run in list(bold_dict[sub][ses].keys()):
                bold_list = bold_dict[sub][ses][run]
                if len(bold_list) > 1:
                    raise ValueError(
                        f'Duplicate was found for bold files sub- {sub}, ses- {ses} and run {run}: {str(bold_list)}')
                file = bold_list[0]
                bold_scan_list.append(file)
                if bold_only:
                    structural_scan_list.append(file)
                    filename_template = pathlib.Path(
                        file).name.rsplit(".nii")[0]
                    split_name.append(filename_template)
                    scan_info.append(
                        {'subject_id': sub, 'session': ses, 'run': run})
                else:
                    run_iter[filename_template].append(run)

    return split_name, scan_info, run_iter, structural_scan_list, bold_scan_list


class BIDSDataGraberInputSpec(BaseInterfaceInputSpec):
    bids_dir = traits.Str(exists=True, mandatory=True,
                          desc="BIDS data directory")
    bids_filter = traits.Dict(exists=True, mandatory=True,
                         desc="BIDS specs")
    scan_info = traits.Dict(exists=True, mandatory=True,
                            desc="Info required to find the scan")
    run = traits.Any(exists=True, desc="Run number")


class BIDSDataGraberOutputSpec(TraitedSpec):
    out_file = File(
        exists=True, desc="Selected file based on the provided parameters.")


class BIDSDataGraber(BaseInterface):
    """
    This interface will select a single scan from the BIDS directory based on the
    input specifications.
    """

    input_spec = BIDSDataGraberInputSpec
    output_spec = BIDSDataGraberOutputSpec

    def _run_interface(self, runtime):
        if 'run' in (self.inputs.scan_info.keys()):
            run = self.inputs.scan_info['run']
        else:
            run = self.inputs.run

        bids_filter = self.inputs.bids_filter.copy()
        bids_filter['subject'] = self.inputs.scan_info['subject_id']
        bids_filter['session'] = self.inputs.scan_info['session']
        if not run is None:
            bids_filter['run'] = run

        from bids.layout import BIDSLayout
        layout = BIDSLayout(self.inputs.bids_dir, validate=False)
        try:
            file_list = layout.get(extension=['nii', 'nii.gz'], return_type='filename', **bids_filter)
            if len(file_list) > 1:
                raise ValueError(f'Provided BIDS spec lead to duplicates: {bids_filter}')
            elif len(file_list)==0:
                raise ValueError(f'No file for found corresponding to the following BIDS spec: {bids_filter}')
        except:
            raise ValueError(f'Error with BIDS spec: {bids_filter}')

        setattr(self, 'out_file', file_list[0])

        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file')}


###WRAPPERS FOR AFNI'S FUNCTIONS; NECESSARY TO PREVENT ISSUES WHEN READING INPUTS/OUTPUTS FROM WORKFLOW GRAPH

def apply_despike(in_file):
    import pathlib
    import os
    from rabies.utils import run_command
    split = pathlib.Path(in_file).name.rsplit(".nii")[0]
    out_file = os.path.abspath(f"{split}_despike.nii.gz")
    command = f'3dDespike -prefix {out_file} {in_file}'
    rc,c_out = run_command(command)
    return out_file

def apply_autobox(in_file):
    import pathlib
    import os
    from rabies.utils import run_command
    split = pathlib.Path(in_file).name.rsplit(".nii")[0]
    out_file = os.path.abspath(f"{split}_autobox.nii.gz")
    command = f'3dAutobox -input {in_file} -prefix {out_file} -npad 1'
    rc,c_out = run_command(command)
    return out_file


def convert_to_RAS(img_file, out_dir=None):
    # convert the input image to the RAS orientation convention
    import os
    import nibabel as nb
    img = nb.load(img_file)
    if nb.aff2axcodes(img.affine) == ('R', 'A', 'S'):
        return img_file
    else:
        import pathlib  # Better path manipulation
        split = pathlib.Path(img_file).name.rsplit(".nii")
        if out_dir is None:
            out_file = os.path.abspath(split[0]+'_RAS.nii.gz')
        else:
            out_file = out_dir+'/'+split[0]+'_RAS.nii.gz'
            if not os.path.isdir(out_dir):
                os.makedirs(out_dir)
        nb.as_closest_canonical(img).to_filename(out_file)
        return out_file


def correct_oblique_affine(input):
    from nipype import logging
    log = logging.getLogger('nipype.workflow')

    import numpy as np
    import nibabel as nb
    img = nb.load(input)
    # check if the image is oblique
    if (nb.affines.obliquity(img.affine)!=0).sum()==0:
        # no axis is oblique, thus don't apply correction
        log.info(f"The file {input} is not oblique. No corrections applied.")
        return input
    
    log.info(f"The file {input} is considered oblique. The affine matrix is modified.")
    # nb.affines.obliquity assumes the largest dimension is the one that should correspond to voxel size
    aff = img.affine
    new_aff = np.zeros([4,4])
    vs = nb.affines.voxel_sizes(aff)
    vs_ratio = np.abs(aff[:-1, :-1] / vs)
    for ax in [0,1,2]:
        idx = np.where(vs_ratio[:,ax]==vs_ratio[:,ax].max())[0][0]
        new_aff[idx,ax]=vs[ax]*np.sign(aff[idx,ax])
        
    # copy 4th row/column
    new_aff[3,:] = aff[3,:]
    new_aff[:,3] = aff[:,3]

    import pathlib
    import os
    split = pathlib.Path(input).name.rsplit(".nii")[0]
    output = os.path.abspath(f"{split}_2card.nii.gz")
    nb.Nifti1Image(img.dataobj, new_aff, img.header).to_filename(output)
    return output


def convert_3dWarp(input):
    import pathlib
    import os
    from rabies.utils import run_command
    # apply AFNI's 3dWarp to convert from oblique to cartesian
    split = pathlib.Path(input).name.rsplit(".nii")[0]
    output = os.path.abspath(f"{split}_2card.nii.gz")
    command = f'3dWarp -oblique2card -prefix {output} {input}'
    rc = run_command(command)
    return output


def resample_template(opts, structural_scan_list, bold_scan_list):
    import os
    import SimpleITK as sitk
    import numpy as np
    from rabies.utils import resample_image_spacing, run_command
    from nipype import logging
    log = logging.getLogger('nipype.workflow')

    template_file = str(opts.anat_template)
    mask_file = str(opts.brain_mask)
    anatomical_resampling = opts.anatomical_resampling
    commonspace_resampling = opts.commonspace_resampling
    rabies_data_type = opts.data_type

    def sitk_get_spacing(f):
        reader = sitk.ImageFileReader()
        reader.SetFileName(f)
        reader.LoadPrivateTagsOn()
        reader.ReadImageInformation()
        return reader.GetSpacing()

    if anatomical_resampling == 'inputs_defined':
        spacing_list = [sitk_get_spacing(f)[:3] for f in structural_scan_list]
        
        # find the lowest dimension across all images and all axes
        low_dim = np.asarray([spacing for spacing in spacing_list]).flatten().min()
        registration_spacing = (low_dim, low_dim, low_dim)
    else:
        shape = anatomical_resampling.split('x')
        registration_spacing = (float(shape[0]), float(shape[1]), float(shape[2]))

    log.info(f"Resampling template to {registration_spacing[0]}x{registration_spacing[1]}x{registration_spacing[2]}mm dimensions for registration steps.")
    registration_template = os.path.abspath("registration_template.nii.gz")
    sitk.WriteImage(resample_image_spacing(sitk.ReadImage(
        template_file, rabies_data_type), registration_spacing), registration_template)

    # also resample the brain mask to ensure stable registrations further down
    registration_mask = os.path.abspath("registration_mask.nii.gz")
    command = f'antsApplyTransforms -d 3 -i {mask_file} -r {registration_template} -o {registration_mask} --verbose -n GenericLabel'
    rc,c_out = run_command(command)

    if commonspace_resampling == 'inputs_defined':
        # if bold_only, bold_scan_list and structural_scan_list are the same; avoid re-loading images if already computed above
        if not (anatomical_resampling == 'inputs_defined' and opts.bold_only):
            spacing_list = [sitk_get_spacing(f)[:3] for f in bold_scan_list]
        
        # create a list of all the types of spacing across listed files
        spacing_types = [spacing_list[0]]
        for spacing in spacing_list:
            if not spacing in spacing_types:
                spacing_types.append(spacing)

        if len(spacing_types)>1:
            raise ValueError(f"""The following list of image dimensions have been sampled across the input functional images: {spacing_types}. 
    Because of the conflicting image dimensions, a single resampling resolution cannot be inferred for commonspace. 
    Please manually define the desired commonspace resolution with the --commonspace_resampling parameter to pursue 
    preprocessing.""")
        else:
            commonspace_spacing = spacing_types[0]
        
    else:
        shape = commonspace_resampling.split('x')
        commonspace_spacing = (float(shape[0]), float(shape[1]), float(shape[2]))


    log.info(f"Images will be resampled to a commonspace resolution of {commonspace_spacing[0]}x{commonspace_spacing[1]}x{commonspace_spacing[2]}mm.")
    commonspace_template = os.path.abspath("commonspace_template.nii.gz")
    sitk.WriteImage(resample_image_spacing(sitk.ReadImage(
        template_file, rabies_data_type), commonspace_spacing), commonspace_template)

    return registration_template, registration_mask, commonspace_template

