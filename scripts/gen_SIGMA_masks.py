#! /usr/bin/env python
# /// script
# requires-python = "==3.9.*"
# dependencies = [
#     "numpy==1.26.4",
#     "scipy==1.13.1",
#     "simpleitk==2.5.0",
# ]
# ///
"""Derive the RABIES rat template set from a SIGMA release.

SIGMA distributes probabilistic tissue maps, while RABIES requires binary masks,
and distributes its label descriptions in the ITK-SNAP format rather than as a
CSV. This script produces the files of the rat template set from an unpacked
SIGMA release, and is run once when building the distributed bundle rather than
on a user's machine, so that every installation gets identical files.

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
from scipy.ndimage import binary_erosion

# the two SIGMA spaces used by RABIES: the anatomical template for the standard
# pipeline, and the functional template for --bold_only
VARIANTS = ['Anatomical', 'Functional']
# the anatomical template is named _template, the functional one _epi
TEMPLATE_SUFFIX = {'Anatomical': 'template', 'Functional': 'epi'}
# the EPI grid is too coarse to erode: one iteration empties the CSF mask, so only
# the anatomical CSF mask is eroded, as in the mouse set
ERODED = {'Anatomical': True, 'Functional': False}
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
  * The probabilistic CSF maps were thresholded at {tissue_threshold} to give the binary
    masks RABIES requires. The anatomical CSF mask was then eroded by
    {erosion_iterations} iteration(s); the functional one was not, since that grid is too
    coarse to erode without emptying it.
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


def binarize(in_file, out_file, threshold, erosion_iterations=0):
    img = sitk.ReadImage(in_file)
    array = sitk.GetArrayFromImage(img)
    # SIGMA probabilistic maps are scaled to [0,1], but a few are stored as percentages
    if array.max() > 1:
        array = array/array.max()
    mask = array >= threshold
    if mask.sum() == 0:
        raise ValueError(f"Thresholding {in_file} at {threshold} left an empty mask.")
    if erosion_iterations > 0:
        # eroding keeps partial volume voxels at tissue boundaries out of the
        # nuisance timecourses
        mask = binary_erosion(mask, iterations=erosion_iterations)
        if mask.sum() == 0:
            raise ValueError(f"Eroding {in_file} by {erosion_iterations} iterations left "
                             "an empty mask; the grid is too coarse to erode.")
    out_img = sitk.GetImageFromArray(mask.astype('int16'), isVector=False)
    out_img.CopyInformation(img)
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
    parser.add_argument('--tissue_threshold', type=float, default=0.9,
                        help="probability above which a voxel belongs to CSF "
                             "(default: %(default)s)")
    parser.add_argument('--mask_threshold', type=float, default=0.5,
                        help="probability above which a voxel belongs to the brain "
                             "(default: %(default)s)")
    parser.add_argument('--erosion_iterations', type=int, default=1,
                        help="erosion applied to the CSF mask of the anatomical "
                             "template; the EPI grid is too coarse to erode "
                             "(default: %(default)s)")
    opts = parser.parse_args()

    os.makedirs(opts.out_dir, exist_ok=True)

    # CC-BY requires attribution to travel with the redistributed files
    with open(os.path.join(opts.out_dir, 'ATTRIBUTION.txt'), 'w') as handle:
        handle.write(ATTRIBUTION.format(
            wm_structures='\n'.join(f'      - {structure}' for structure in WM_STRUCTURES['Anatomical']),
            tissue_threshold=opts.tissue_threshold,
            mask_threshold=opts.mask_threshold,
            erosion_iterations=opts.erosion_iterations))

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
        erosion = opts.erosion_iterations if ERODED[variant] else 0
        name = 'eroded_csf_mask' if ERODED[variant] else 'csf_mask'
        outputs.append(binarize(
            find_file(opts.sigma_dir, f'{prefix}_csf.nii.gz'),
            os.path.join(opts.out_dir, f'{prefix}_{name}.nii.gz'),
            opts.tissue_threshold, erosion_iterations=erosion))

        atlas = f'{prefix}_Atlas.nii.gz'
        out_atlas = os.path.join(opts.out_dir, atlas)
        shutil.copyfile(find_file(opts.sigma_dir, atlas), out_atlas)
        outputs.append(out_atlas)

        if variant in WM_STRUCTURES:
            outputs.append(roi_mask(
                out_atlas, find_file(opts.sigma_dir, LABEL_DESCRIPTION[variant]),
                WM_STRUCTURES[variant], os.path.join(opts.out_dir, f'{prefix}_wm_mask.nii.gz')))

        n_labels = convert_label_description(
            find_file(opts.sigma_dir, LABEL_DESCRIPTION[variant]),
            os.path.join(opts.out_dir, MAPPING_NAME[variant]))

        for file in outputs:
            check_overlap(out_template, file)

        print(f"{variant}: wrote {len(outputs)+1} files and {n_labels} labels "
              f"to {opts.out_dir}")


if __name__ == '__main__':
    main()
