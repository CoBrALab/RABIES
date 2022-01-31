import os
import pathlib
import SimpleITK as sitk
import numpy as np
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from rabies.utils import run_command

def prep_bids_iter(layout, bold_only=False):
    '''
    This function takes as input a BIDSLayout, and generates iteration lists
    for managing the workflow's iterables depending on whether --bold_only is
    selected or not.
    '''
    import pathlib

    scan_info = []
    split_name = []
    scan_list = []
    run_iter = {}
    bold_scan_list = []

    subject_list = layout.get_subject()
    if len(subject_list) == 0:
        raise ValueError(
            "No subject information could be retrieved from the BIDS directory. The 'sub-' specification is mandatory.")
    if not bold_only:
        anat_bids = layout.get(subject=subject_list, suffix=[
                               'T2w', 'T1w'], extension=['nii', 'nii.gz'])
        if len(anat_bids) == 0:
            raise ValueError(
                "No anatomical file with the suffix 'T2w' or 'T1w' were found among the BIDS directory.")
    bold_bids = layout.get(subject=subject_list, suffix=[
                           'bold', 'cbv'], extension=['nii', 'nii.gz'])
    if len(bold_bids) == 0:
        raise ValueError(
            "No functional file with the suffix 'bold' or 'cbv' were found among the BIDS directory.")

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

        bold_list = layout.get(subject=sub, session=ses, run=run, suffix=[
                               'bold', 'cbv'], extension=['nii', 'nii.gz'], return_type='filename')
        bold_dict[sub][ses][run] = bold_list

    # if not bold_only, then the bold_list and run_iter will be a dictionary with keys being the anat filename
    # otherwise, it will be a list of bold scans themselves
    for sub in list(bold_dict.keys()):
        for ses in list(bold_dict[sub].keys()):
            if not bold_only:
                anat_list = layout.get(subject=sub, session=ses, suffix=[
                                       'T2w', 'T1w'], extension=['nii', 'nii.gz'], return_type='filename')
                if len(anat_list) == 0:
                    raise ValueError(
                        f'Missing an anatomical image for sub {sub} and ses- {ses}')
                if len(anat_list) > 1:
                    raise ValueError(
                        f'Duplicate was found for the anatomical file sub- {sub}, ses- {ses}: {str(anat_list)}')
                file = anat_list[0]
                scan_list.append(file)
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
                    scan_list.append(file)
                    filename_template = pathlib.Path(
                        file).name.rsplit(".nii")[0]
                    split_name.append(filename_template)
                    scan_info.append(
                        {'subject_id': sub, 'session': ses, 'run': run})
                else:
                    run_iter[filename_template].append(run)

    return split_name, scan_info, run_iter, scan_list, bold_scan_list


class BIDSDataGraberInputSpec(BaseInterfaceInputSpec):
    bids_dir = traits.Str(exists=True, mandatory=True,
                          desc="BIDS data directory")
    suffix = traits.List(exists=True, mandatory=True,
                         desc="Suffix to search for")
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
        subject_id = self.inputs.scan_info['subject_id']
        session = self.inputs.scan_info['session']
        if 'run' in (self.inputs.scan_info.keys()):
            run = self.inputs.scan_info['run']
        else:
            run = self.inputs.run

        from bids.layout import BIDSLayout
        layout = BIDSLayout(self.inputs.bids_dir, validate=False)
        try:
            if run is None: # if there is no run spec to search, don't include it in the search
                file_list = layout.get(subject=subject_id, session=session, extension=[
                                  'nii', 'nii.gz'], suffix=self.inputs.suffix, return_type='filename')
            else:
                file_list = layout.get(subject=subject_id, session=session, run=run, extension=[
                                  'nii', 'nii.gz'], suffix=self.inputs.suffix, return_type='filename')
            if len(file_list) > 1:
                raise ValueError(f'Provided BIDS spec lead to duplicates: \
                    {str(self.inputs.suffix)} sub-{subject_id} ses-{session} run-{str(run)}')
            elif len(file_list)==0:
                raise ValueError(f'No file for found corresponding to the following BIDS spec: \
                    {str(self.inputs.suffix)} sub-{subject_id} ses-{session} run-{str(run)}')
        except:
            raise ValueError(f'Error with BIDS spec: \
                    {str(self.inputs.suffix)} sub-{subject_id} ses-{session} run-{str(run)}')

        setattr(self, 'out_file', file_list[0])

        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file')}


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
    rc = run_command(command)

    return resampled_template, resampled_mask

