# Input data requirements

The input dataset must be organised according to the
[BIDS data structure](https://bids.neuroimaging.io/) {cite}`Gorgolewski2016-zm`.

## How RABIES traverses a dataset

RABIES iterates through every subject found to contain a functional file, and
within each subject through the sessions and runs present.

When anatomical scans are used — that is, when `--bold_only` is not set — each
functional scan is matched to one anatomical scan **from the same subject and
session**.

## Scan identification

```{list-table}
:header-rows: 1
:widths: 30 30 40

* - Image type
  - Default BIDS suffixes
  - Parameter
* - Functional
  - `bold`, `cbv`
  - `--bids_filter`, key `func`
* - Structural
  - `T1w`, `T2w`
  - `--bids_filter`, key `anat`
```

Files matching neither filter are ignored. The default filter is equivalent to:

```{code-block} json
:caption: Default value of --bids_filter

{
    "func": {
        "suffix": ["bold", "cbv"]
    },
    "anat": {
        "suffix": ["T1w", "T2w"]
    }
}
```

```{seealso}
[How to select which scans get processed](../how_to/select_scans.md) for
customising the filter and for selecting individual scans.
```

## Image orientation

RABIES expects images in the NIfTI standard RAS+ orientation
(Right–Anterior–Superior). Incorrectly oriented images are a common source of
registration failures — see
[How to check image orientation](../how_to/check_orientation.md).

## Example dataset

The [RABIES example dataset](http://doi.org/10.5281/zenodo.8349029)
(`test_dataset.zip`) has the following structure:

```{code-block} text
:caption: Two subjects, one session each, with paired anatomical and functional scans

test_dataset/
├── sub-PHG001
│   └── ses-3
│       ├── anat
│       │   ├── sub-PHG001_ses-3_acq-RARE_T2w.json
│       │   └── sub-PHG001_ses-3_acq-RARE_T2w.nii.gz
│       └── func
│           ├── sub-PHG001_ses-3_task-rest_acq-EPI_run-1_bold.json
│           └── sub-PHG001_ses-3_task-rest_acq-EPI_run-1_bold.nii.gz
└── sub-PHG002
    └── ses-3
        ├── anat
        │   ├── sub-PHG002_ses-3_acq-RARE_T2w.json
        │   └── sub-PHG002_ses-3_acq-RARE_T2w.nii.gz
        └── func
            ├── sub-PHG002_ses-3_task-rest_acq-EPI_run-1_bold.json
            └── sub-PHG002_ses-3_task-rest_acq-EPI_run-1_bold.nii.gz

8 directories, 8 files
```

This is the dataset used in [the tutorial](../tutorials/first_run.md).

## Format conversion

Conversion from Bruker raw format to NIfTI can be handled with
[BrkRaw](https://brkraw.github.io/). The CoBrALab maintains
[notes on the conversion](https://github.com/CoBrALab/documentation/wiki/bruker2nifti-conversion).
