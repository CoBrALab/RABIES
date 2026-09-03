# How to handle container syntax

A container has its own filesystem and cannot see your data unless you say so.
Running RABIES in a container is therefore the same as running it natively,
plus one rule: **every directory RABIES needs must be bound to a path inside
the container, and the RABIES arguments must use the container-side paths.**

Bind directories with `-B` for Apptainer and `-v` for Docker. Both take
`host_path:container_path`, with an optional `:ro` to make the bind read-only.

```{important}
Bind the same directories at the same container-side paths for **all three
stages**. Each stage reads the file paths recorded by the previous one, so a
path used during preprocessing must still resolve during confound correction
and analysis. Changing or dropping a bind between stages produces missing-file
errors.
```

## Apptainer

### Preprocessing

```sh
apptainer run -B $PWD/input_BIDS:/input_BIDS:ro \
  -B $PWD/preprocess_outputs:/preprocess_outputs/ \
  /path_to_apptainer_image/rabies.sif \
  -p MultiProc preprocess /input_BIDS/ /preprocess_outputs/ \
  --apply_STC --TR 1.2 \
  --commonspace_reg masking=true,brain_extraction=false,template_registration=SyN,fast_commonspace=false
```

`apptainer run /path_to_apptainer_image/rabies.sif` executes the image; every
argument after it is passed to RABIES and follows the ordinary
[command line syntax](../reference/cli.md). The two binds are what make the
data reachable:

`-B $PWD/input_BIDS:/input_BIDS:ro`
: Maps your BIDS folder to `/input_BIDS` inside the container, which is why the
  RABIES argument reads `/input_BIDS/`. `:ro` grants read-only access, so the
  container cannot modify your raw data.

`-B $PWD/preprocess_outputs:/preprocess_outputs/`
: Maps the desired output directory. There is no `:ro`, so the container can
  write here.

### Confound correction

```sh
apptainer run -B $PWD/input_BIDS:/input_BIDS:ro \
  -B $PWD/preprocess_outputs:/preprocess_outputs/ \
  -B $PWD/confound_correction_outputs:/confound_correction_outputs/ \
  /path_to_apptainer_image/rabies.sif \
  -p MultiProc confound_correction /preprocess_outputs/ /confound_correction_outputs/ \
  --nuisance_regressors WM_signal CSF_signal vascular_signal mot_6 \
  --smoothing_filter 0.3
```

`/input_BIDS` is still bound even though it does not appear in the RABIES
arguments — this is the rule stated above.

### Analysis

```sh
apptainer run -B $PWD/input_BIDS:/input_BIDS:ro \
  -B $PWD/preprocess_outputs:/preprocess_outputs/ \
  -B $PWD/confound_correction_outputs:/confound_correction_outputs/ \
  -B $PWD/analysis_outputs:/analysis_outputs/ \
  /path_to_apptainer_image/rabies.sif \
  -p MultiProc analysis /confound_correction_outputs /analysis_outputs/ \
  --group_ica apply=true,dim=30,random_seed=1
```

## Docker

The syntax mirrors Apptainer, with `-v` in place of `-B` and a few extra flags:

```sh
docker run -it --rm --user $(id -u) \
  -v $PWD/input_BIDS:/input_BIDS:ro \
  -v $PWD/preprocess_outputs:/preprocess_outputs/ \
  ghcr.io/cobralab/rabies:latest \
  -p MultiProc preprocess /input_BIDS/ /preprocess_outputs/ \
  --apply_STC --TR 1.2 \
  --commonspace_reg masking=true,brain_extraction=false,template_registration=SyN,fast_commonspace=false
```

`--user $(id -u)`
: Runs as your own user id, so output files are owned by you. Without it,
  Docker writes files as root and you may be unable to delete them.

`--rm`
: Removes the container when the run finishes.

Replace `latest` with a specific version tag for reproducible runs.

## Using a custom atlas or seed files

Template files, masks and seeds passed with `--anat_template`, `--brain_mask`,
`--WM_mask`, `--CSF_mask`, `--vascular_mask`, `--prior_maps` or `--seed_list`
live outside your input and output directories, so they need binds of their
own:

```sh
apptainer run -B $PWD/input_BIDS:/input_BIDS:ro \
  -B $PWD/preprocess_outputs:/preprocess_outputs/ \
  -B $PWD/my_atlas:/atlas:ro \
  /path_to_apptainer_image/rabies.sif \
  -p MultiProc preprocess /input_BIDS/ /preprocess_outputs/ \
  --anat_template /atlas/template.nii.gz \
  --brain_mask /atlas/brain_mask.nii.gz
```

```{seealso}
[CoBrALab recommendations](https://github.com/CoBrALab/documentation/wiki/Running-RABIES-on-niagara)
for running RABIES on Compute Canada clusters.
```
