#! /usr/bin/env python
"""Derive the RABIES rat template set from a SIGMA release.

SIGMA distributes probabilistic tissue maps, while RABIES requires binary masks,
and distributes its label descriptions in the ITK-SNAP format rather than as a
CSV. This script produces the files of the rat template set from an unpacked
SIGMA release, and is run once when building the distributed bundle rather than
on a user's machine, so that every installation gets identical files.

Usage: gen_SIGMA_masks.py <unpacked SIGMA directory> <output directory>
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
# the anatomical masks are eroded, as in the mouse set
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


def convert_label_description(in_file, out_file):
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
                        help="probability above which a voxel belongs to a tissue "
                             "(default: %(default)s)")
    parser.add_argument('--mask_threshold', type=float, default=0.5,
                        help="probability above which a voxel belongs to the brain "
                             "(default: %(default)s)")
    parser.add_argument('--erosion_iterations', type=int, default=1,
                        help="erosion applied to the WM and CSF masks of the anatomical "
                             "template; the EPI grid is too coarse to erode "
                             "(default: %(default)s)")
    opts = parser.parse_args()

    os.makedirs(opts.out_dir, exist_ok=True)

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
        name = 'eroded_{tissue}_mask' if ERODED[variant] else '{tissue}_mask'
        for tissue in ['wm', 'csf']:
            outputs.append(binarize(
                find_file(opts.sigma_dir, f'{prefix}_{tissue}.nii.gz'),
                os.path.join(opts.out_dir,
                             f'{prefix}_' + name.format(tissue=tissue) + '.nii.gz'),
                opts.tissue_threshold, erosion_iterations=erosion))

        atlas = f'{prefix}_Atlas.nii.gz'
        out_atlas = os.path.join(opts.out_dir, atlas)
        shutil.copyfile(find_file(opts.sigma_dir, atlas), out_atlas)
        outputs.append(out_atlas)

        n_labels = convert_label_description(
            find_file(opts.sigma_dir, LABEL_DESCRIPTION[variant]),
            os.path.join(opts.out_dir, MAPPING_NAME[variant]))

        for file in outputs:
            check_overlap(out_template, file)

        print(f"{variant}: wrote {len(outputs)+1} files and {n_labels} labels "
              f"to {opts.out_dir}")


if __name__ == '__main__':
    main()
