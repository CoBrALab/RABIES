# Understanding the Outputs

In this section, there is a description for all the output files provided at each processing stage. Important outputs from RABIES are stored into `datasink/` folders, which will be generated in the output folder specified at execution.

## Preprocessing Outputs

Multiple datasink folders are generated during preprocessing for different output types: `anat_datasink/`, `bold_datasink/`, `unbiased_template_datasink/`,  `transforms_datasink/` and `confounds_datasink/`. 

- `anat_datasink/`: Includes the inhomogeneity-correction anatomical scans.
    - `anat_preproc/`: anatomical scans after inhomogeneity correction

- `bold_datasink/`: Includes all outputs related to the functional scans, where files are either resampled onto the native or commonspace of the EPI. The native space outputs are resampled over the anatomical scan from each corresponding MRI session, whereas the commonspace outputs are resampled over the reference atlas (the original EPI voxel resolution is unchanged during resampling unless specified otherwise in the RABIES command).
    - `native_bold/`: preprocessed EPI timeseries resampled to nativespace
    - `native_brain_mask/`: brain mask in nativespace
    - `native_WM_mask/`: WM mask in nativespace
    - `native_CSF_mask/`: CSF mask in nativespace
    - `native_labels/`: atlas labels in nativespace
    - `native_bold_ref/`: a volumetric 3D EPI average generated from the 4D `native_bold/`
    - `commonspace_bold/`: preprocessed EPI timeseries resampled to commonspace
    - `commonspace_mask/`: brain mask in commonspace
    - `commonspace_WM_mask/`: WM mask in commonspace
    - `commonspace_CSF_mask/`: CSF mask in commonspace
    - `commonspace_vascular_mask/`: vascular mask in commonspace
    - `commonspace_labels/`: atlas labels in commonspace
    - `commonspace_resampled_template/`: the commonspace anatomical template, resampled to the EPI's dimensions
    - `input_bold/`: the raw EPI scans provided as inputs in the BIDS data folder
    - `initial_bold_ref/`: the initial volumetric 3D EPI average generated from the 4D `input_bold/`
    - `raw_brain_mask/`: brain mask resampled onto the 4D `input_bold/`
    - `inho_cor_bold/`: the volumetric 3D EPI (`initial_bold_ref/`) after inhomogeneity correction, which is later used for registration of the EPI
    - `inho_cor_bold_warped2anat/`: inho_cor_bold after co-registration to the associated anatomical image (`anat_preproc/`)
    - `std_map_preprocess/`: the temporal standard deviation at each voxel on the `commonspace_bold/`
    - `tSNR_map_preprocess/`: the temporal signal-to-noise ratio (tSNR) of the `commonspace_bold/`

- `unbiased_template_datasink/`: Outputs related to the generation of the unbiased template using https://github.com/CoBrALab/optimized_antsMultivariateTemplateConstruction. The unbiased template corresponds to the average of all anatomical (or functional with `--bold_only`) scans after their alignment.
    - `unbiased_template/`: the unbiased template generated from the input dataset scans
    - `warped_unbiased_template/`: the unbiased template, registered to the reference atlas in commonspace

- `transforms_datasink/`: datasink for all the relevant transform files resampling between the different spaces. The bold_to_anat registration transformed the raw EPI to overlap with the anatomical image, correcting for susceptibility distortions, which corresponds to the native space. The native_to_unbiased registration overlaps every scans to the generated unbiased template, and then the unbiased_to_atlas corresponds to the registration of the unbiased template with the reference atlas, which defines the commonspace.
    - `bold_to_anat_affine/`: affine transforms from the EPI co-registration to the anatomical image
    - `bold_to_anat_warp/`: non-linear transforms from the EPI co-registration to the anatomical image
    - `bold_to_anat_inverse_warp/`: inverse of the non-linear transforms from the EPI co-registration to the anatomical image
    - `native_to_unbiased_affine/`: affine transforms for the alignment between native space and the unbiased template
    - `native_to_unbiased_warp/`: non-linear transforms for the alignment between native space and the unbiased template
    - `native_to_unbiased_inverse_warp/`: inverse of the non-linear transforms  for the alignment between native space and the unbiased template
    - `unbiased_to_atlas_affine/`: affine transforms for the alignment between unbiased template and the atlas in commonspace
    - `unbiased_to_atlas_warp/`: non-linear transforms for the alignment between unbiased template and the atlas in commonspace
    - `unbiased_to_atlas_inverse_warp/`: inverse of the non-linear transforms for the alignment between unbiased template and the atlas in commonspace

- `motion_datasink/`: files derivated from motion estimation
    - `motion_params_csv/`: contains the 24 motion parameters which can be used as nuisance regressors at the confound correction pipeline stage.
    - `FD_csv/`: a CSV file with timescourses for either the mean or maximal framewise displacement (FD) estimations.
    - `FD_voxelwise/`: a Nifti image which contains framewise displacement evaluated at each voxel
    - `pos_voxelwise/`: a Nifti image which tracks the displacement (derived from the head motion realignment parameters) of each voxel across time


## Confound Correction Outputs
Important outputs from confound correction will be found in the `confound_correction_datasink/`:
- `confound_correction_datasink/`:
    - `cleaned_timeseries/`: cleaned timeseries after the application of confound correction
    - `frame_censoring_mask/`: contains CSV files each recording as a boolean vector which timepoints were censored if frame censoring was applied.
    - `aroma_out/`: if `--ica_aroma` is applied, this folder contains outputs from running ICA-AROMA, which includes the MELODIC ICA outputs and the component classification results
    - `plot_CR_overfit/`: will contain figures illustrating the variance explained by random regressors during confound correction, and the variance explained by the real regressors after substrating the variance from random regressors.
    - `background_masking_fig/`: will illustrate the image background masks automatically generated if using `--image_scaling background_noise`


## Analysis Outputs

Outputs from analyses will be found in the `analysis_datasink/`, whereas outputs relevant to the `--data_diagnosis` are found in `data_diagnosis_datasink/`:
- `analysis_datasink/`:
    - `group_ICA_dir/`: complete output from MELODIC ICA, which the melodic_IC.nii.gz Nifti which gives all spatial components, and `report/` folder which includes a HTML visualization. 
    - `matrix_data_file/`: .pkl file which contains a 2D numpy array representing the whole-brain correlation matrix. If `--ROI_type parcellated` is selected, the row/column indices of the array are matched in increasing order of the atlas ROI label number.
    - `matrix_fig/`: .png file which displays the correlation matrix
    - `seed_correlation_maps/`: nifti files for seed-based connectivity analysis, where each seed provided in `--seed_list` has an associated voxelwise correlation maps
    - `dual_regression_nii/`: the spatial maps from dual regression, which correspond to the linear coefficients from the second regression. The list of 3D spatial maps obtained are concatenated into a 4D Nifti file, where the order of component is consistent with the priors provided in `--prior_maps`.
    - `dual_regression_timecourse_csv/`: a CSV file which stores the outputs from the first linear regression during dual regression. This corresponds to a timecourse associated to each prior component from `--prior_maps`.
    - `NPR_prior_filename/`: spatial components fitted during NPR
    - `NPR_prior_timecourse_csv/`: timecourses associated to each components from NPR_prior_filename
    - `NPR_extra_filename/`: the extra spatial components fitted during NPR which were not part of priors
    - `NPR_extra_timecourse_csv/`: timecourses associated to each components from NPR_extra_filename
(diagnosis_datasink_target)=
- `data_diagnosis_datasink/`:
    - `figure_temporal_diagnosis/`: figure which displays scan-level temporal features from the [spatiotemporal diagnosis](diagnosis_target)
    - `figure_spatial_diagnosis/`: figure which displays scan-level spatial features from the [spatiotemporal diagnosis](diagnosis_target)
    - `analysis_QC/`: group-level features of data quality from `--data_diagnosis`
        - `sample_distributions/`: contains the [distribution plots](dist_plot_target)
            - `{analysis}_sample_distribution.png`: the distribution plot for a given network analysis
            - `{analysis}_outlier_detection.csv`: a CSV which associates the measures displayed in the distribution plot with corresponding scan IDs
        - `parametric_stats/`: [group statistical report](group_stats_target) for analysis quality control (using parametric measures)
            - `DR{component #}_QC_maps.png`: The _QC_maps.png files are displaying statistical maps relevant to analysis quality control. The DR refers to dual regression analysis, and the {component #} is relating the file to one of the BOLD components specified in `--prior_bold_idx`
            - `DR{component #}_QC_stats.csv`: a follow-up to _QC_maps.png which allows for the quantitative categorization of data quality outcomes as in {cite}`Desrosiers-Gregoire2024-ou`
            - `seed_FC{seed #}_QC_maps.png`: same statistical maps as with `DR{component #}_QC_maps.png`, but for seed-based connectivity analysis
            - `seed_FC{seed #}_QC_stats.csv`: same measures as with `DR{component #}_QC_maps.png`, but for seed-based connectivity analysis
        - `non_parametric_stats/`: same as `parametric_stats/`, but using non-parametric measures
    - `temporal_info_csv/`: CSV file containing the data plotted with `figure_temporal_diagnosis/`
    - `spatial_VE_nii/`: Nifti file with the confound regression percentage variance explained (R^2) at each voxel
    - `CR_prediction_std_nii/`: Nifti file with the confound regression variance explained at each voxel
    - `random_CR_std_nii/`: Nifti file with the variance explained from random regressors at each voxel
    - `corrected_CR_std_nii/`: Nifti file with the confound regression variance explained at each voxel after removing the variance explained by random regressors
    - `temporal_std_nii/`: the standard deviation at each voxel after confound correction
    - `GS_cov_nii/`: the covariance of each voxel with the global signal
