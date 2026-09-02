# How to troubleshoot registration

The RABIES defaults modify the images as little as possible and lean on the
quality of the images at acquisition. They are the right starting point, but
they do not produce a robust workflow for every dataset, and reaching good
outcomes on your own data will often mean tuning parameters.

This guide maps failures visible in the
[preprocessing QC report](../reference/qc_outputs.md) onto the parameters that
address them. Work through it in order: inhomogeneity correction feeds masking,
and masking feeds every registration that follows, so a problem there will
reappear downstream no matter what you change at the registration step.

```{important}
Before changing any parameter, confirm your images are correctly oriented —
see [How to check image orientation](check_orientation.md). Mis-oriented
images produce registration failures that no amount of parameter tuning will
fix.
```

## Inhomogeneity correction

Relevant parameters: `--anat_inho_cor`, `--bold_inho_cor`,
`--anat_robust_inho_cor`, `--bold_robust_inho_cor`. Inspect the
`anat_inho_cor/` and `bold_inho_cor/` folders of the QC report.

Only a subset of scans have failed masking, or the mask is partially misregistered
: Use `--anat_robust_inho_cor`/`--bold_robust_inho_cor`. These register all
  corrected images together to generate a temporary template representing the
  average of all scans, mask that template, and use it as the masking target in
  a second pass of inhomogeneity correction, which is a more robust target. The
  parameters for these options are the same as for `--commonspace_reg` below.

Inhomogeneity biases are not completely corrected, with signal drops remaining
: Apply `multiotsu=true`. This performs a staged correction, correcting low
  intensities first and iteratively including higher ones, which better handles
  images with strong inhomogeneity gradients and very low intensities.

Tissue outside the brain is causing registration failures
: If the initial correction enhanced the intensity of tissue outside the brain
  and masking then fails, use `--anat_autobox`/`--bold_autobox` to crop out the
  extra tissue automatically. You can also modify `otsu_thresh`, which sets the
  threshold for the automatic masking during the initial correction, to select
  a value more specific to brain tissue.

A large proportion of masking failures remain — mismatched brain sizes, non-linear warps, or the mask falling outside the brain
: Apply a less stringent registration `method`, stepping down through
  `SyN` → `Affine` → `Rigid` → `no_reg`. If you reach `no_reg`, you may also
  have to adjust `otsu_thresh` to obtain an automatically generated brain mask
  covering only brain tissue.

## Commonspace registration and susceptibility distortion correction

Relevant parameters: `--commonspace_reg`, `--bold2anat_coreg`. Inspect the
`commonspace_reg_wf.Anat2Unbiased/` and `commonspace_reg_wf.Unbiased2Atlas/`
folders of the QC report (or `commonspace_reg_wf.Anat2Atlas/` if you used
`fast_commonspace=true`), together with `EPI2Anat/`.

Many scans are misregistered, or brain edges are not well matched
: First inspect the quality of inhomogeneity correction for those scans, and
  follow the guidance above if the correction or brain masking was poor. If
  good quality masks were obtained during inhomogeneity correction, bring them
  into the registration with `masking=true`. If registration errors persist,
  particularly at the brain edges, `brain_extraction=true` and `keep_mask_after_extract=true` further constrains
  edge matching by removing tissue outside the brain.

  ```{warning}
  The quality of brain edge delineation depends on the masks derived during
  inhomogeneity correction, so `brain_extraction=true` is only as good as that
  earlier step. Fix masking first.
  ```

Scans have incomplete brain coverage, and surrounding tissue is stretched to fill the gap
: Non-linear registration assumes corresponding anatomy between the moving
  image and the target. When brain regions are missing — the cerebellum or
  olfactory bulbs are the usual cases — the surrounding tissue may be
  improperly stretched to fill the missing area. `brain_extraction=true,keep_mask_after_extract=true` can
  largely mitigate this.

## When registration cannot be salvaged

If a scan cannot be registered acceptably, exclude it rather than analysing it:
see [How to select which scans get processed](select_scans.md).

```{seealso}
- [Preprocessing QC outputs](../reference/qc_outputs.md) — what each QC folder shows
- [The preprocessing workflow](../explanation/preprocessing.md) — what each registration step is for
- [`preprocess` options](../reference/cli.md) — the complete parameter list, with all accepted values
```
