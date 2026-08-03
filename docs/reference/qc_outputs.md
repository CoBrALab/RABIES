# Preprocessing QC outputs

Several registration operations during preprocessing can fail to align images
accurately. RABIES automatically generates a set of PNG images allowing
efficient visual assessment of the key registration steps, so that failed
alignment is caught before it biases downstream analyses.

The images are written to `{output_folder}/preprocess_QC_report/`, one
subfolder per registration step or per supporting file group.

```{important}
Inspect this report on every run. A failed registration does not raise an
error — it produces plausible-looking output that is wrong.
```

## Which folders you get

Some folder names carry the name of the workflow that produced them as a
prefix, and which commonspace folders exist depends on `--commonspace_reg`.

```{list-table}
:header-rows: 1
:widths: 45 55

* - Folder
  - Present when
* - `anat_inho_cor/`
  - structural scans are used, i.e. not `--bold_only`
* - `bold_inho_cor/`
  - always
* - `commonspace_reg_wf.Anat2Unbiased/`
  - `fast_commonspace=false` (the default)
* - `commonspace_reg_wf.Unbiased2Atlas/`
  - `fast_commonspace=false` (the default)
* - `commonspace_reg_wf.unbiased_template_masking/`
  - `fast_commonspace=false` and `masking=true`
* - `commonspace_reg_wf.Anat2Atlas/`
  - `fast_commonspace=true`
* - `anat_robust_inho_cor_template.*/`
  - `--anat_robust_inho_cor apply=true`
* - `bold_robust_inho_cor_template.*/`
  - `--bold_robust_inho_cor apply=true`
* - `EPI2Anat/`
  - structural scans are used, i.e. not `--bold_only`
* - `template_files/`
  - always
* - `temporal_features/`
  - always
```

The `*_robust_inho_cor_template.*` folders contain the same registration
figures as the `commonspace_reg_wf.*` ones, for the template built during the
robust inhomogeneity correction pass rather than for the main commonspace
registration.

## `anat_inho_cor/`

Quality of the intensity inhomogeneity correction applied to the structural
image, which is performed before the important registration operations and is
crucial for their performance.

Each figure has 4 columns: **1** the raw image, **2** an initial correction of
the image, **3** an overlay of the anatomical mask used for the final
correction (by default obtained through a preliminary registration to the
commonspace template), and **4** the final corrected output.

```{figure} ../pics/sub-MFC067_ses-1_acq-FLASH_T1w_inho_cor.png
:alt: Four-column figure showing the stages of structural inhomogeneity correction

Structural inhomogeneity correction.
```

## `bold_inho_cor/`

The same as `anat_inho_cor/`, but for the 3D reference EPI image used to
estimate the alignment of the EPI.

```{figure} ../pics/sub-MFC068_ses-1_task-rest_acq-EPI_run-1_bold_inho_cor.png
:alt: Four-column figure showing the stages of EPI inhomogeneity correction

Functional inhomogeneity correction.
```

## `commonspace_reg_wf.Anat2Unbiased/`

Alignment between each anatomical image and the generated unbiased template.
This registration controls the overlap between different scanning sessions.

```{figure} ../pics/sub-MFC067_ses-1_acq-FLASH_T1w_inho_cor_registration.png
:alt: Overlap between a structural scan and the dataset-generated unbiased template

Structural scan (top) against the unbiased template (bottom).
```

## `commonspace_reg_wf.Unbiased2Atlas/`

Alignment of the generated unbiased template to the external anatomical
template in commonspace. This ensures proper alignment with the commonspace and
its associated brain parcellation.

```{figure} ../pics/atlas_registration.png
:alt: Overlap between the unbiased template and the reference atlas template

Unbiased template (top) against the reference atlas template (bottom).
```

## `commonspace_reg_wf.Anat2Atlas/`

Produced instead of the two folders above when `fast_commonspace=true`, which
skips the unbiased template and registers each scan directly to the reference
atlas. It shows that direct alignment, one figure per scan.

## `EPI2Anat/`

Alignment of the EPI image to the anatomical image from the same scanning
session. This step resamples the EPI into native space and corrects
susceptibility distortions through non-linear registration.

```{figure} ../pics/sub-MFC068_ses-1_task-rest_acq-EPI_run-1_bold_registration.png
:alt: Overlap between the volumetric EPI and the structural image

Volumetric EPI (top) against the structural image (bottom).
```

## `template_files/`

Overlap of the provided external anatomical template with its associated masks
and labels. Use it to confirm the correct template files were provided, and
share it alongside the RABIES report.

```{figure} ../pics/template_files.png
:alt: The anatomical template overlaid with its masks and labels

Template, masks and labels.
```

## `temporal_features/`

The timecourse of the head motion realignment parameters together with
framewise displacement, showing subject motion. Also includes a spatial map of
the signal variability at each voxel, and the temporal signal-to-noise ratio
(tSNR).

```{figure} ../pics/example_temporal_features.png
:alt: Motion parameter timecourses, framewise displacement, variability and tSNR maps

Temporal features.
```

```{seealso}
- [How to troubleshoot registration](../how_to/troubleshoot_registration.md) — what to change when these figures show a failure
- [Metric definitions](metrics.md) — how framewise displacement and the other quantities are computed
- [Data quality assessment](../explanation/data_quality.md) — the separate, analysis-stage quality reports
```
