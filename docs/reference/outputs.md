# Output files

Every processing stage writes its important outputs into `datasink/` folders,
created inside the output directory given at execution. This page describes
every file produced.

```{note}
The output directory also contains the Nipype working directory and the
serialised workflow state (`.pkl`). Those are internal, but they are what the
next stage reads to locate the previous stage's files — do not move or delete
them between stages.
```

## Preprocessing outputs

Preprocessing writes five datasinks: `anat_datasink/`, `bold_datasink/`,
`unbiased_template_datasink/`, `transforms_datasink/` and `motion_datasink/`.
The quality control images are written separately to `preprocess_QC_report/`,
described in [Preprocessing QC outputs](qc_outputs.md).

### `anat_datasink/`

Inhomogeneity-corrected anatomical scans.

- `anat_preproc/`: anatomical scans after inhomogeneity correction

### `bold_datasink/`

All outputs related to the functional scans. Files are resampled either onto
the native space or the commonspace of the EPI. Native space outputs are
resampled over the anatomical scan from the corresponding MRI session;
commonspace outputs are resampled over the reference atlas. The original EPI
voxel resolution is unchanged during resampling unless specified otherwise in
the RABIES command.

**Native space**

- `native_bold/`: preprocessed EPI timeseries resampled to native space
- `native_brain_mask/`: brain mask in native space
- `native_WM_mask/`: WM mask in native space
- `native_CSF_mask/`: CSF mask in native space
- `native_labels/`: atlas labels in native space
- `native_bold_ref/`: a volumetric 3D EPI average generated from the 4D `native_bold/`

**Commonspace**

- `commonspace_bold/`: preprocessed EPI timeseries resampled to commonspace
- `commonspace_mask/`: brain mask in commonspace
- `commonspace_WM_mask/`: WM mask in commonspace
- `commonspace_CSF_mask/`: CSF mask in commonspace
- `commonspace_vascular_mask/`: vascular mask in commonspace
- `commonspace_labels/`: atlas labels in commonspace
- `commonspace_resampled_template/`: the commonspace anatomical template, resampled to the EPI's dimensions

**Inputs and intermediates**

- `input_bold/`: the raw EPI scans provided as inputs in the BIDS data folder
- `initial_bold_ref/`: the initial volumetric 3D EPI average generated from the 4D `input_bold/`
- `raw_brain_mask/`: brain mask resampled onto the 4D `input_bold/`
- `inho_cor_bold/`: the volumetric 3D EPI (`initial_bold_ref/`) after inhomogeneity correction, later used for registration of the EPI
- `inho_cor_bold_warped2anat/`: `inho_cor_bold` after co-registration to the associated anatomical image (`anat_preproc/`)
- `std_map_preprocess/`: the temporal standard deviation at each voxel of `commonspace_bold/`
- `tSNR_map_preprocess/`: the temporal signal-to-noise ratio (tSNR) of `commonspace_bold/`

### `unbiased_template_datasink/`

Outputs from the generation of the unbiased template using
[optimized_antsMultivariateTemplateConstruction](https://github.com/CoBrALab/optimized_antsMultivariateTemplateConstruction).
The unbiased template is the average of all anatomical scans (or functional
scans with `--bold_only`) after their alignment.

- `unbiased_template/`: the unbiased template generated from the input dataset scans
- `warped_unbiased_template/`: the unbiased template, registered to the reference atlas in commonspace

### `transforms_datasink/`

All transform files for resampling between spaces.

The `bold_to_anat` registration transforms the raw EPI to overlap with the
anatomical image, correcting susceptibility distortions, which defines native
space. The `native_to_unbiased` registration overlaps every scan onto the
generated unbiased template. The `unbiased_to_atlas` registration aligns the
unbiased template with the reference atlas, which defines commonspace.

- `bold_to_anat_affine/`: affine transforms from the EPI co-registration to the anatomical image
- `bold_to_anat_warp/`: non-linear transforms from the EPI co-registration to the anatomical image
- `bold_to_anat_inverse_warp/`: inverse of `bold_to_anat_warp/`
- `native_to_unbiased_affine/`: affine transforms for the alignment between native space and the unbiased template
- `native_to_unbiased_warp/`: non-linear transforms for the same alignment
- `native_to_unbiased_inverse_warp/`: inverse of `native_to_unbiased_warp/`
- `unbiased_to_atlas_affine/`: affine transforms for the alignment between the unbiased template and the atlas in commonspace
- `unbiased_to_atlas_warp/`: non-linear transforms for the same alignment
- `unbiased_to_atlas_inverse_warp/`: inverse of `unbiased_to_atlas_warp/`

### `motion_datasink/`

Files derived from motion estimation.

- `motion_params_csv/`: the 24 motion parameters, usable as nuisance regressors at the confound correction stage
- `FD_csv/`: a CSV with timecourses for either the mean or maximal [framewise displacement](FD_target) estimations
- `FD_voxelwise/`: a NIfTI image containing framewise displacement evaluated at each voxel
- `pos_voxelwise/`: a NIfTI image tracking the displacement of each voxel across time, derived from the head motion realignment parameters

## Confound correction outputs

### `confound_correction_datasink/`

- `cleaned_timeseries/`: cleaned timeseries after the application of confound correction
- `frame_censoring_mask/`: CSV files recording, as a boolean vector, which timepoints were censored, if frame censoring was applied
- `aroma_out/`: outputs from running ICA-AROMA if `--ica_aroma` is applied, including the MELODIC ICA outputs and the component classification results
- `plot_CR_overfit/`: figures illustrating the variance explained by random regressors during confound correction, and the variance explained by the real regressors after subtracting the variance from random regressors

## Analysis outputs

### `commonspace_analysis_datasink/` and `nativespace_analysis_datasink/`

Which of these appears depends on the space the cleaned timeseries were
produced in:

```{list-table}
:header-rows: 1
:widths: 45 55

* - Datasink
  - Present when
* - `commonspace_analysis_datasink/`
  - confound correction produced commonspace timeseries, or `--resample_to_commonspace` was passed, or `--group_ica apply=true` was used
* - `nativespace_analysis_datasink/`
  - confound correction produced nativespace timeseries
```

Both contain the same set of files, listed below. `group_ICA_dir/` is always
written to `commonspace_analysis_datasink/`, since group ICA requires
commonspace alignment.

- `group_ICA_dir/`: complete output from MELODIC ICA, including the `melodic_IC.nii.gz` NIfTI giving all spatial components and a `report/` folder with an HTML visualisation
- `matrix_data_file/`: a `.pkl` file containing a 2D NumPy array representing the whole-brain correlation matrix. With `--ROI_type parcellated`, the row and column indices are matched in increasing order of the atlas ROI label number
- `matrix_fig/`: a `.png` displaying the correlation matrix
- `seed_correlation_maps/`: NIfTI files for [seed-based connectivity](SBC_target), one voxelwise correlation map per seed provided in `--seed_list`
- `dual_regression_nii/`: the spatial maps from [dual regression](DR_target), corresponding to the linear coefficients from the second regression. The 3D spatial maps are concatenated into a 4D NIfTI, with component order consistent with the priors provided in `--prior_maps`
- `dual_regression_timecourse_csv/`: a CSV storing the outputs from the first linear regression during dual regression — one timecourse per prior component from `--prior_maps`
- `NPR_prior_filename/`: spatial components fitted during NPR
- `NPR_prior_timecourse_csv/`: timecourses associated with each component in `NPR_prior_filename/`
- `NPR_extra_filename/`: the extra spatial components fitted during NPR which were not part of the priors
- `NPR_extra_timecourse_csv/`: timecourses associated with each component in `NPR_extra_filename/`

(diagnosis_datasink_target)=

### `data_diagnosis_datasink/`

Produced when `--data_diagnosis` is selected. See
[How to assess data quality](../how_to/assess_data_quality.md) for how to use
these.

- `figure_temporal_diagnosis/`: scan-level temporal features from the [spatiotemporal diagnosis](diagnosis_target)
- `figure_spatial_diagnosis/`: scan-level spatial features from the [spatiotemporal diagnosis](diagnosis_target)
- `temporal_info_csv/`: CSV containing the data plotted in `figure_temporal_diagnosis/`
- `spatial_VE_nii/`: NIfTI with the confound regression percentage variance explained ($R^2$) at each voxel
- `CR_prediction_std_nii/`: NIfTI with the confound regression variance explained at each voxel
- `random_CR_std_nii/`: NIfTI with the variance explained from random regressors at each voxel
- `corrected_CR_std_nii/`: NIfTI with the confound regression variance explained at each voxel after removing the variance explained by random regressors
- `temporal_std_nii/`: the standard deviation at each voxel after confound correction
- `GS_cov_nii/`: the covariance of each voxel with the global signal

`analysis_QC/` holds the group-level features of data quality:

- `sample_distributions/`: the [distribution plots](dist_plot_target)
    - `{analysis}_sample_distribution.png`: the distribution plot for a given network analysis
    - `{analysis}_outlier_detection.csv`: a CSV associating the measures displayed in the distribution plot with the corresponding scan IDs
- `parametric_stats/`: the [group statistical report](group_stats_target) for analysis quality control, using parametric measures
    - `DR{component #}_QC_maps.png`: statistical maps relevant to analysis quality control. `DR` refers to dual regression analysis, and `{component #}` relates the file to one of the BOLD components specified in `--prior_bold_idx`
    - `DR{component #}_QC_stats.csv`: a follow-up to `_QC_maps.png` allowing quantitative categorisation of data quality outcomes as in {cite}`Desrosiers-Gregoire2024-ou`
    - `seed_FC{seed #}_QC_maps.png`: the same statistical maps, for seed-based connectivity analysis
    - `seed_FC{seed #}_QC_stats.csv`: the same measures, for seed-based connectivity analysis
- `non_parametric_stats/`: as `parametric_stats/`, but using non-parametric measures

```{seealso}
[Metric definitions](metrics.md) gives the precise computation behind every
quantity named on this page.
```
