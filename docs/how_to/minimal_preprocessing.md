# How to use RABIES on already-preprocessed data

If your functional data was preprocessed with your own pipeline and you only
want RABIES for confound correction and analysis, you cannot skip the
preprocessing stage — but you can run it with almost everything turned off.

Running `rabies confound_correction` directly on external data fails with
missing `.pkl` file errors. Those `.pkl` files are serialised states of the
internal preprocessing workflow and cannot be written by hand. Even with the
`.pkl` in place, the confound correction stage expects the full range of files
produced by preprocessing, so missing file errors would follow.

The practical solution is a *SHAM* preprocessing run, with the correction and
registration steps disabled. RABIES still computes the intermediary outputs that
the later stages require, while leaving the image data largely unchanged.

## Run a SHAM preprocessing

```sh
rabies preprocess bids_inputs/ preprocess_outputs/ \
  --anat_inho_cor method=disable \
  --bold_inho_cor method=disable \
  --commonspace_reg template_registration=no_reg,fast_commonspace=true \
  --bold2anat_coreg registration=no_reg \
  --no_HMC
```

```{list-table}
:header-rows: 1
:widths: 40 60

* - Parameter
  - What it turns off
* - `--anat_inho_cor method=disable`
  - inhomogeneity correction of the structural images
* - `--bold_inho_cor method=disable`
  - inhomogeneity correction of the EPI images
* - `--commonspace_reg template_registration=no_reg,fast_commonspace=true`
  - unbiased template generation, and the registration to the reference atlas (an identity transform is used instead)
* - `--bold2anat_coreg registration=no_reg`
  - EPI-to-structural coregistration, i.e. the susceptibility distortion correction
* - `--no_HMC`
  - the *application* of head motion correction to the resampled timeseries. Head motion parameters are still estimated and remain available to `--nuisance_regressors`, `--frame_censoring` and `--data_diagnosis`
```

If your dataset has no structural scans, add `--bold_only`, in which case
`--anat_inho_cor` and `--bold2anat_coreg` no longer apply.

## Your data must already be in commonspace

```{warning}
`template_registration=no_reg` does not skip the resampling to commonspace — it
replaces the estimated transform with an identity transform. The commonspace
outputs are therefore only meaningful if your input images **already overlap**
with the template given to `--anat_template`, which is the file that defines
the commonspace.

If they do not overlap, the commonspace timeseries and the atlas masks
(`--brain_mask`, `--WM_mask`, `--CSF_mask`, `--vascular_mask`) applied
downstream will not correspond to your data.
```

Confirm your image orientation before you start — see
[How to check image orientation](check_orientation.md).

## What still happens

A SHAM preprocessing is a minimal pass, not a strictly non-modifying one.

Operations which alter the data but are off by default stay off, and should not
be added: `--apply_STC`, `--apply_despiking`, `--detect_dummy`,
`--log_transform`, `--anat_autobox`, `--bold_autobox` and `--oblique2card`.

The timeseries are still resampled onto the output grid, using the identity
transforms described above. Use `--commonspace_resampling` and
`--anatomical_resampling` to control the output voxel dimensions, and
`--interpolation` to select the interpolator.

Two further operations are applied unconditionally and cannot be turned off:

**Negative values are clipped to zero**
: This happens when the preprocessed timeseries are written out. If your data
  legitimately contains negative values — because it was already demeaned or
  detrended by your own pipeline — those voxels will be set to zero. Bring in
  data on a positive scale, and leave centring to `--detrending_order` at the
  confound correction stage.

**The output is cast** to the type given by `--data_type`
: `float32` by default.

## Alternative: `--read_datasink`

`rabies confound_correction --read_datasink` reads the preprocessing outputs
from the datasink folders rather than from the saved workflow graph, which
removes the need for the `.pkl` file. This requires reproducing the RABIES
[output structure and file naming](../reference/outputs.md) exactly, and is
generally more work than running a SHAM preprocessing.

## Still not covered?

If your use case needs settings that are not exposed, open a
[discussion](https://github.com/CoBrALab/RABIES/discussions) describing
explicitly what you need. Providing example data lets us work out an
implementation supporting your use case.
