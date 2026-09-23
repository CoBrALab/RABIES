#! /usr/bin/env python
# /// script
# requires-python = "==3.9.*"
# dependencies = [
#     "numpy==1.26.4",
#     "simpleitk==2.5.0",
# ]
# ///
"""Derive the RABIES rat template set from a SIGMA release.

RABIES requires binary white matter and CSF masks, which SIGMA distributes only as
probabilistic maps, and label descriptions as a CSV rather than in the ITK-SNAP
format. This script produces the files of the rat template set from an unpacked
SIGMA release, building the masks from structures of the anatomical atlas, and is
run once when building the distributed bundle rather than on a user's machine, so
that every installation gets identical files.

Usage: uv run scripts/gen_SIGMA_masks.py <unpacked SIGMA directory> <output directory>

The versions pinned above match rabies_environment.yml.
"""

import argparse
import csv
import os
import re
import shutil
import numpy as np
import SimpleITK as sitk

# the two SIGMA spaces used by RABIES: the anatomical template for the standard
# pipeline, and the functional template for --bold_only. The anatomical variant comes
# first, since the functional CSF mask is derived from the anatomical one.
VARIANTS = ['Anatomical', 'Functional']
# the anatomical template is named _template, the functional one _epi
TEMPLATE_SUFFIX = {'Anatomical': 'template', 'Functional': 'epi'}
# the functional atlas label description does not follow the image naming
LABEL_DESCRIPTION = {
    'Anatomical': 'SIGMA_InVivo_Anatomical_Brain_Atlas.txt',
    'Functional': 'SIGMA_Functional_Brain_Atlas_Labels.txt',
    }
MAPPING_NAME = {
    'Anatomical': 'SIGMA_InVivo_Anatomical_Brain_Atlas_mapping.csv',
    'Functional': 'SIGMA_Functional_Brain_Atlas_mapping.csv',
    }
# the white matter mask is built from atlas structures rather than by thresholding the
# probabilistic white matter map, which pulls in thalamus; regressing thalamic signal
# removes signal of interest. It is not eroded: rat tracts are a few voxels thick on
# the anatomical grid, and one iteration removes 70% of the mask. The functional atlas
# has no white matter structures, so that variant gets no white matter mask.
WM_STRUCTURES = {
    'Anatomical': [
        'corpus callosum and associated subcortical white matter',
        'corticofugal tract and corona radiata',
        'anterior commissure, anterior limb',
        'anterior commissure, posterior limb',
        'anterior commissure, intrabulbar part',
        'fimbria of the hippocampus',
        'alveus of the hippocampus',
        'fornix',
        'ventral hippocampal commissure',
        'optic tract and optic chiasm',
        'posterior commissure',
        'inferior cerebellar peduncle',
        'middle cerebellar peduncle',
        'superior cerebellar peduncle and prerubral field',
        ],
    }
# the CSF mask is the ventricular system of the atlas, as the mouse CSF mask is that of
# the DSURQE atlas. The probabilistic CSF map is not used: thresholded, it covers the
# surface of the brainstem and cerebellum rather than the ventricles, and overlaps the
# white matter mask. It is not eroded: rat ventricles are one or two voxels thick on
# the anatomical grid, and one iteration removes 81% of the mask, including nearly all
# of the lateral ventricles. The functional atlas has no ventricles, so that variant
# gets the anatomical mask resampled onto its grid.
CSF_STRUCTURES = {
    'Anatomical': [
        'Ventricular system, unspecified',
        '4th ventricle',
        'Central canal',
        ],
    }


ATTRIBUTION = """RABIES rat commonspace template set
===================================

Derived from the SIGMA Wistar rat brain templates and atlases, version 2.0,
distributed under the Creative Commons Attribution 4.0 International licence
(CC-BY-4.0).

  Source:  https://github.com/DavidBarriere/SIGMA-Rat-Brain-Templates-and-Atlases
  Archive: https://zenodo.org/records/10635831

Please cite:

  Barriere, D.A. et al. The SIGMA rat brain templates and atlases for multimodal
  MRI data analysis and visualization. Nature Communications 10, 5699 (2019).
  https://doi.org/10.1038/s41467-019-13575-7

The parcellation embedded in SIGMA is the Waxholm Space atlas of the rat brain,
normalized into SIGMA space by the SIGMA authors. Please also cite:

  Kleven, H. et al. Waxholm Space atlas of the rat brain: a 3D atlas supporting
  data analysis and integration. Nature Methods 20, 1822-1829 (2023).
  https://doi.org/10.1038/s41592-023-02034-3

Modifications made for RABIES
-----------------------------

Produced with scripts/gen_SIGMA_masks.py from the RABIES repository:

  * The anatomical white matter mask is the union of the following structures of the
    anatomical atlas. It is not eroded, since rat white matter tracts are only a few
    voxels thick on that grid:
{wm_structures}
    No white matter mask is provided for the functional template, whose atlas has no
    white matter structures.
  * The anatomical CSF mask is the union of the following structures of the anatomical
    atlas. It is not eroded, since rat ventricles are only one or two voxels thick on
    that grid:
{csf_structures}
    The functional CSF mask is the anatomical one resampled onto the functional
    template grid with linear interpolation and thresholded at 0.5, since the
    functional atlas has no ventricles.
  * The distributed brain masks were re-binarized at {mask_threshold}.
  * The ITK-SNAP label descriptions were converted to CSV.
  * The templates and atlas label images are unmodified copies.

Only the in-vivo anatomical and functional resources are included. The ex-vivo,
diffusion and CT/PET resources of the SIGMA release are not redistributed here.
"""


def find_file(sigma_dir, name):
    # SIGMA releases nest the files in modality folders whose layout has changed
    # between versions, so the files are located by name
    for root, dirs, files in os.walk(sigma_dir):
        if name in files:
            return os.path.join(root, name)
    raise ValueError(f"{name} was not found under {sigma_dir}.")


def binarize(in_file, out_file, threshold):
    img = sitk.ReadImage(in_file)
    array = sitk.GetArrayFromImage(img)
    # SIGMA probabilistic maps are scaled to [0,1], but a few are stored as percentages
    if array.max() > 1:
        array = array/array.max()
    mask = array >= threshold
    if mask.sum() == 0:
        raise ValueError(f"Thresholding {in_file} at {threshold} left an empty mask.")
    out_img = sitk.GetImageFromArray(mask.astype('int16'), isVector=False)
    out_img.CopyInformation(img)
    sitk.WriteImage(out_img, out_file)
    return out_file


def resample_mask(mask_file, reference_file, out_file):
    # the SIGMA templates share one world space, so a mask is carried onto the grid of
    # another template without registration; the interpolated mask is thresholded at
    # one half
    reference = sitk.ReadImage(reference_file)
    fraction = sitk.Resample(sitk.Cast(sitk.ReadImage(mask_file), sitk.sitkFloat32),
                             reference, sitk.Transform(), sitk.sitkLinear, 0.0,
                             sitk.sitkFloat32)
    mask = sitk.GetArrayFromImage(fraction) >= 0.5
    if mask.sum() == 0:
        raise ValueError(f"Resampling {mask_file} onto {reference_file} left an empty mask.")
    out_img = sitk.GetImageFromArray(mask.astype('int16'), isVector=False)
    out_img.CopyInformation(reference)
    sitk.WriteImage(out_img, out_file)
    return out_file


def read_label_description(in_file):
    # ITK-SNAP label descriptions list, per line:
    # IDX R G B A VIS MSH "LABEL"
    rows = []
    with open(in_file) as handle:
        for line in handle:
            line = line.strip()
            if len(line) == 0 or line.startswith('#'):
                continue
            match = re.match(r'^(\d+)\s+(?:\S+\s+){6}"(.*)"\s*$', line)
            if match is None:
                continue
            label, structure = int(match.group(1)), match.group(2)
            if label == 0: # the background label carries no structure
                continue
            rows.append([structure, label])
    if len(rows) == 0:
        raise ValueError(f"No labels could be read from {in_file}.")
    return rows


def roi_mask(atlas_file, label_description, structures, out_file):
    # names are matched exactly, so that a structure renamed or renumbered in a later
    # SIGMA release fails here instead of silently dropping out of the mask
    labels = {}
    for structure, label in read_label_description(label_description):
        labels.setdefault(structure, []).append(label)
    missing = [structure for structure in structures if structure not in labels]
    if len(missing) > 0:
        raise ValueError(f"{label_description} has no structure named {missing}.")
    img = sitk.ReadImage(atlas_file)
    array = sitk.GetArrayFromImage(img)
    mask = np.isin(array, [label for structure in structures for label in labels[structure]])
    if mask.sum() == 0:
        raise ValueError(f"The structures {structures} are empty in {atlas_file}.")
    out_img = sitk.GetImageFromArray(mask.astype('int16'), isVector=False)
    out_img.CopyInformation(img)
    sitk.WriteImage(out_img, out_file)
    return out_file


def convert_label_description(in_file, out_file):
    rows = read_label_description(in_file)
    with open(out_file, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['Structure', 'label'])
        writer.writerows(rows)
    return len(rows)


def check_overlap(template, file):
    # the same criterion RABIES applies to template files at runtime
    template_img = sitk.ReadImage(template)
    img = sitk.ReadImage(file)
    if not (template_img.GetOrigin() == img.GetOrigin()
            and template_img.GetDirection() == img.GetDirection()):
        raise ValueError(f"{file} does not overlap with {template}.")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('sigma_dir', help="path to an unpacked SIGMA release")
    parser.add_argument('out_dir', help="directory the derived files are written to")
    parser.add_argument('--mask_threshold', type=float, default=0.5,
                        help="probability above which a voxel belongs to the brain "
                             "(default: %(default)s)")
    opts = parser.parse_args()

    os.makedirs(opts.out_dir, exist_ok=True)

    # CC-BY requires attribution to travel with the redistributed files
    with open(os.path.join(opts.out_dir, 'ATTRIBUTION.txt'), 'w') as handle:
        handle.write(ATTRIBUTION.format(
            wm_structures='\n'.join(f'      - {structure}' for structure in WM_STRUCTURES['Anatomical']),
            csf_structures='\n'.join(f'      - {structure}' for structure in CSF_STRUCTURES['Anatomical']),
            mask_threshold=opts.mask_threshold))

    for variant in VARIANTS:
        prefix = f'SIGMA_InVivo_{variant}_Brain'

        template = f'{prefix}_{TEMPLATE_SUFFIX[variant]}.nii.gz'
        out_template = os.path.join(opts.out_dir, template)
        shutil.copyfile(find_file(opts.sigma_dir, template), out_template)

        # the distributed brain mask is binarized rather than trusted to be binary
        out_mask = binarize(find_file(opts.sigma_dir, f'{prefix}_mask.nii.gz'),
                            os.path.join(opts.out_dir, f'{prefix}_mask.nii.gz'),
                            opts.mask_threshold)

        outputs = [out_mask]

        atlas = f'{prefix}_Atlas.nii.gz'
        out_atlas = os.path.join(opts.out_dir, atlas)
        shutil.copyfile(find_file(opts.sigma_dir, atlas), out_atlas)
        outputs.append(out_atlas)

        if variant in WM_STRUCTURES:
            outputs.append(roi_mask(
                out_atlas, find_file(opts.sigma_dir, LABEL_DESCRIPTION[variant]),
                WM_STRUCTURES[variant], os.path.join(opts.out_dir, f'{prefix}_wm_mask.nii.gz')))

        out_csf = os.path.join(opts.out_dir, f'{prefix}_csf_mask.nii.gz')
        if variant in CSF_STRUCTURES:
            anatomical_csf = roi_mask(
                out_atlas, find_file(opts.sigma_dir, LABEL_DESCRIPTION[variant]),
                CSF_STRUCTURES[variant], out_csf)
            outputs.append(anatomical_csf)
        else:
            outputs.append(resample_mask(anatomical_csf, out_template, out_csf))

        n_labels = convert_label_description(
            find_file(opts.sigma_dir, LABEL_DESCRIPTION[variant]),
            os.path.join(opts.out_dir, MAPPING_NAME[variant]))

        for file in outputs:
            check_overlap(out_template, file)

        # the template and the label mapping are written alongside the checked outputs
        print(f"{variant}: wrote {len(outputs)+2} files and {n_labels} labels "
              f"to {opts.out_dir}")


if __name__ == '__main__':
    main()
