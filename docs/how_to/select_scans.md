# How to select which scans get processed

By default RABIES processes every functional scan it can find in the input BIDS
directory. There are two independent ways to narrow that down: change which
files count as functional or anatomical scans, or list the individual scans to
keep or drop.

## Change which files RABIES recognises

RABIES identifies functional scans by the `bold` or `cbv` suffix and structural
scans by the `T1w` or `T2w` suffix. Files matching neither are ignored.

Override this with `--bids_filter`, which takes a JSON file with the functional
filter under `func` and the structural filter under `anat`. This JSON reproduces
the defaults:

```{code-block} json
:caption: bids_filter.json

{
    "func": {
        "suffix": ["bold", "cbv"]
    },
    "anat": {
        "suffix": ["T1w", "T2w"]
    }
}
```

Pass it at the preprocessing stage:

```sh
rabies preprocess input_BIDS/ preprocess_outputs/ --bids_filter bids_filter.json
```

Add any BIDS entity to either filter to be more specific. To use only the RARE
anatomical acquisition and only the resting-state functional runs:

```{code-block} json
:caption: bids_filter.json

{
    "func": {
        "suffix": ["bold"],
        "task": ["rest"]
    },
    "anat": {
        "suffix": ["T2w"],
        "acquisition": ["RARE"]
    }
}
```

```{tip}
Reach for `--bids_filter` when your dataset contains scans that are not meant
for this pipeline at all — a second anatomical modality, a task run alongside
the resting-state runs. Reach for `--inclusion_ids` below when you want a
subset of otherwise-eligible scans.
```

## Include or exclude individual scans

`--inclusion_ids` and `--exclusion_ids` take the full paths of BOLD files.
They are execution options, so they go **before** the processing stage name,
and they can be given at any stage.

Process only two scans:

```sh
rabies --inclusion_ids input_BIDS/sub-001/ses-1/func/sub-001_ses-1_task-rest_bold.nii.gz \
                       input_BIDS/sub-002/ses-1/func/sub-002_ses-1_task-rest_bold.nii.gz \
  -p MultiProc preprocess input_BIDS/ preprocess_outputs/
```

Process everything except one scan:

```sh
rabies --exclusion_ids input_BIDS/sub-003/ses-1/func/sub-003_ses-1_task-rest_bold.nii.gz \
  -p MultiProc preprocess input_BIDS/ preprocess_outputs/
```

For longer lists, put one file path per row in a text file and pass the file:

```sh
rabies --inclusion_ids scans_to_process.txt -p MultiProc preprocess input_BIDS/ preprocess_outputs/
```

```{warning}
Do not put `--inclusion_ids` or `--exclusion_ids` immediately before the
processing stage name. The stage name gets swallowed into the list and
argument parsing fails. Put another option, such as `-p` or `--verbose`,
between them — as in every example above.

`--inclusion_ids` and `--exclusion_ids` cannot be used together.
```

## Drop scans after quality control

The two options above are also how you act on a quality control decision. To
carry out an analysis without scans you have judged unusable, re-run the
confound correction and analysis stages with `--exclusion_ids` listing them.
See [How to assess data quality](assess_data_quality.md) for how to arrive at
that decision, and for `--scan_QC_thresholds`, which excludes scans by
threshold rather than by name.

```{seealso}
[Input data requirements](../reference/bids_inputs.md) for how RABIES pairs
functional scans with anatomical scans, and the
[`preprocess` options](../reference/cli.md) for the full parameter list.
```
