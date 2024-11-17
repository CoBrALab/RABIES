from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from niworkflows.utils.connections import pop_file, listify

import pdb
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
    echo_iter = {}
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
        try:
            echo = bold.get_entities()['echo']
        except:
            echo = None

        if sub not in list(bold_dict.keys()):
            bold_dict[sub] = {}
        if ses not in list(bold_dict[sub].keys()):
            bold_dict[sub][ses] = {}
        if run not in list(bold_dict[sub][ses].keys()):
            bold_dict[sub][ses][run] = {}

        bold_list = layout.get(subject=sub, session=ses, run=run,echo=echo, 
                               extension=['nii', 'nii.gz'], return_type='filename',**bids_filter['func'])
        bold_dict[sub][ses][run][echo] = bold_list

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
                echo_iter[run] = []
                run_iter[filename_template].append(run)
                run_scans = []  # Initialize the run_scans list for each run
                for echo in list(bold_dict[sub][ses][run].keys()):
                    bold_list = bold_dict[sub][ses][run][echo]
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
                            {'subject_id': sub, 'session': ses, 'run': run, 'echo': echo})
                    else:
                        echo_iter[run].append(echo)


    number_functional_scans = len(bold_scan_list)
    return split_name, scan_info, run_iter,echo_iter, structural_scan_list, number_functional_scans,bold_scan_list


class BIDSDataGraberInputSpec(BaseInterfaceInputSpec):
    bids_dir = traits.Str(exists=True, mandatory=True,
                          desc="BIDS data directory")
    bids_filter = traits.Dict(exists=True, mandatory=True,
                         desc="BIDS specs")
    scan_info = traits.Dict(exists=True, mandatory=True,
                            desc="Info required to find the scan")
    run = traits.Any(exists=True, desc="Run number")
    echo = traits.Any(exists=True, desc="Echo number")

class BIDSDataGraberSingleEchoInputSpec(BaseInterfaceInputSpec):
    bids_dir = traits.Str(exists=True, mandatory=True,
                          desc="BIDS data directory")
    suffix = traits.List(exists=True, mandatory=True,
                         desc="Suffix to search for")
    scan_info = traits.Dict(exists=True, mandatory=True,
                            desc="Info required to find the scan")
    run = traits.Any(exists=True, desc="Run number")
    echo = traits.Any(exists=True, desc="Echo number")


class BIDSDataGraberOutputSpec(TraitedSpec):
    out_file = traits.Either(traits.List(File(exists=True)),
                             File(exists=True),
                             desc="Either a single file or a list of selected files based on the provided parameters.")


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
        if 'echo' in (self.inputs.scan_info.keys()):
            echo = self.inputs.scan_info['echo']
        else:
            echo = self.inputs.echo


        bids_filter = self.inputs.bids_filter.copy()
        bids_filter['subject'] = self.inputs.scan_info['subject_id']
        bids_filter['session'] = self.inputs.scan_info['session']
        if not run is None:
            bids_filter['run'] = run

        bids_filter = self.inputs.bids_filter.copy()
        bids_filter['subject'] = self.inputs.scan_info['subject_id']
        bids_filter['session'] = self.inputs.scan_info['session']
        if not run is None:
            bids_filter['run'] = run

        from bids.layout import BIDSLayout
        layout = BIDSLayout(self.inputs.bids_dir, validate=False)
        try:
            file_list = layout.get(extension=['nii', 'nii.gz'], return_type='filename', **bids_filter)
            '''if len(file_list) > 1:
                raise ValueError(f'Provided BIDS spec lead to duplicates: {bids_filter}')
            elif len(file_list)==0:
                raise ValueError(f'No file for found corresponding to the following BIDS spec: {bids_filter}')
                    '''
        except:
            raise ValueError(f'Error with BIDS spec: {bids_filter}')
        if len(file_list)<2:
            file_list = file_list[0]
        setattr(self, 'out_file', file_list)

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


def resample_template(template_file, mask_file, file_list, spacing='inputs_defined', rabies_data_type=8):
    import os
    import SimpleITK as sitk
    import numpy as np
    from rabies.utils import resample_image_spacing, run_command
    from nipype import logging
    log = logging.getLogger('nipype.workflow')

    if spacing == 'inputs_defined':
        file_list = list(np.asarray(file_list).flatten())
        img = sitk.ReadImage(file_list[0], rabies_data_type)
        low_dim = np.asarray(img.GetSpacing()[:3]).min()
        for file in file_list[1:]:
            img = sitk.ReadImage(file, rabies_data_type)
            new_low_dim = np.asarray(img.GetSpacing()[:3]).min()
            if new_low_dim < low_dim:
                low_dim = new_low_dim
        spacing = (low_dim, low_dim, low_dim)

        template_image = sitk.ReadImage(
            template_file, rabies_data_type)
        template_dim = template_image.GetSpacing()
        if np.asarray(template_dim[:3]).min() > low_dim:
            log.info("The template retains its original resolution.")
            return template_file, mask_file
    else:
        shape = spacing.split('x')
        spacing = (float(shape[0]), float(shape[1]), float(shape[2]))

    log.info(f"Resampling template to {spacing[0]}x{spacing[1]}x{spacing[2]}mm dimensions.")
    resampled_template = os.path.abspath("resampled_template.nii.gz")
    sitk.WriteImage(resample_image_spacing(sitk.ReadImage(
        template_file, rabies_data_type), spacing), resampled_template)

    # also resample the brain mask to ensure stable registrations further down
    resampled_mask = os.path.abspath("resampled_mask.nii.gz")
    command = f'antsApplyTransforms -d 3 -i {mask_file} -r {resampled_template} -o {resampled_mask} --verbose -n GenericLabel'
    rc,c_out = run_command(command)

    return resampled_template, resampled_mask

def extract_entities(file_list):
    """
    Return a dictionary of common entities given a list of files.

    Examples
    --------
    >>> extract_entities("sub-01/anat/sub-01_T1w.nii.gz")
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(["sub-01/anat/sub-01_T1w.nii.gz"] * 2)
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(["sub-01/anat/sub-01_run-1_T1w.nii.gz",
    ...                   "sub-01/anat/sub-01_run-2_T1w.nii.gz"])
    {'subject': '01', 'run': [1, 2], 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}

    """
    from collections import defaultdict
    from bids.layout import parse_file_entities

    # Flatten a list recursively
    def flatten(lis):
        for item in lis:
            if isinstance(item, list):
                for subitem in flatten(item):
                    yield subitem
            else:
                yield item

    file_list = list(flatten(file_list))  # Ensure the input is a flat list

    entities = defaultdict(list)
    for e, v in [
        ev_pair for f in file_list for ev_pair in parse_file_entities(f).items()
    ]:
        entities[e].append(v)

    def _unique(inlist):
        inlist = sorted(set(inlist))
        if len(inlist) == 1:
            return inlist[0]
        return inlist

    return {k: _unique(v) for k, v in entities.items()}


def get_img_orientation(imgf):
    """Return the image orientation as a string"""
    img = nb.load(imgf)
    return "".join(nb.aff2axcodes(img.affine))


def get_estimator(layout, fname):
    field_source = layout.get_metadata(fname).get("B0FieldSource")
    if isinstance(field_source, str):
        field_source = (field_source,)

    if field_source is None:
        import re
        from pathlib import Path

        from sdcflows.fieldmaps import get_identifier

        # Fallback to IntendedFor
        intended_rel = re.sub(r"^sub-[a-zA-Z0-9]*/", "", str(Path(fname).relative_to(layout.root)))
        field_source = get_identifier(intended_rel)

    return field_source

class BIDSDataGraberSingleEcho(BaseInterface):
    """
    This interface will select a single scan from the BIDS directory based on the
    input specifications.
    """

    input_spec = BIDSDataGraberSingleEchoInputSpec
    output_spec = BIDSDataGraberOutputSpec

    def _run_interface(self, runtime):
        subject_id = self.inputs.scan_info['subject_id']
        session = self.inputs.scan_info['session']
        if 'run' in (self.inputs.scan_info.keys()):
            run = self.inputs.scan_info['run']
        else:
            run = self.inputs.run
        if 'echo' in (self.inputs.scan_info.keys()):
            echo = self.inputs.scan_info['echo']
        else:
            echo = self.inputs.echo

        from bids.layout import BIDSLayout
        layout = BIDSLayout(self.inputs.bids_dir, validate=False)
        try:
            if run is None: # if there is no run spec to search, don't include it in the search
                file_list = layout.get(subject=subject_id, session=session, extension=[
                                  'nii', 'nii.gz'], suffix=self.inputs.suffix, return_type='filename')
            else:
                file_list = layout.get(subject=subject_id, session=session, run=run,echo=echo, extension=[
                                  'nii', 'nii.gz'], suffix=self.inputs.suffix, return_type='filename')
            '''if len(file_list) > 1:
                raise ValueError(f'Provided BIDS spec lead to duplicates: \
                    {str(self.inputs.suffix)} sub-{subject_id} ses-{session} run-{str(run)}')
            elif len(file_list)==0:
                raise ValueError(f'No file for found corresponding to the following BIDS spec: \
                    {str(self.inputs.suffix)} sub-{subject_id} ses-{session} run-{str(run)}')
                    '''
        except:
            raise ValueError(f'Error with BIDS spec: \
                    {str(self.inputs.suffix)} sub-{subject_id} ses-{session} run-{str(run)}')
        #if len(file_list)<2:
        file_list = file_list[0]
        setattr(self, 'out_file', file_list)

        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file')}

def remove_echo_string(filename):
    import re
    import shutil
    modified_filename = re.sub(r'echo-\d+_', '', filename)
    #shutil.copyfile(filename, modified_filename) Uncomment this to copy the final output to the initial func directory
    return modified_filename

#logging.getLogger('nipype.workflow').setLevel('DEBUG')

def extract_first_element(list):
    return list[0]

def extract_te_from_json(bold_file):
    import os
    import json
    """
    Extract the echo time (TE) from the JSON sidecar of a BOLD file.
    """
    # Replace .nii.gz or .nii extension with .json to get the JSON sidecar path
    json_file = os.path.splitext(bold_file)[0].rstrip('.nii.gz') + '.json'

    # Read the JSON file
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Extract the echo time
    te = data.get('EchoTime', None)

    return te
