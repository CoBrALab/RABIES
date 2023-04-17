# Connectivity Analysis

Following the completion of the confound correction workflow, RABIES allows the estimation of resting-state connectivity using standard analyses: seed-based connectivity, whole-brain connectivity, group independent component analysis (ICA) and dual regression (DR). For each analysis (except for group-ICA), RABIES will compute individualized connectivity maps for each scan separately, which can then be exported for relevant statistical analyses (e.g. group comparison) conducted outside of RABIES.

## Correlation-based connectivity
Correlation-based analyses rely on computing a temporal correlation between different brain regions' BOLD fluctuations to estimate their functional coupling. The assessment of connectivity using correlation requires careful a priori cleaning of confounds (see [confound correction](confound_correction) and [data quality assessment](analysis_QC_target)), as various fMRI confounds introduce spurious correlations that won't be distinguished from neural activity.
(SBC_target)=
- **Seed-based connectivity** (`--seed_list`): Seed-based connectivity is the first technique developped for the mapping of connectivity during resting state {cite}`Biswal1995-vh`. The mean timecourse is first extracted from an anatomical seed of interest, and the correlation (Pearson’s r in RABIES) between this timecourse and every other voxel is computed to obtain a correlation map, representing the ‘connectivity strength’ between the seed and every other brain regions.
- **Whole-brain connectivity** (`--FC_matrix`/`--ROI_type`): This technique is an extension of the seed-based connectivity technique to encompass every brain region. That is, using the anatomical parcellation provided along with the atlas during preprocessing, the seed timecourse for every parcel is first extracted, and then the cross-correlation (Pearson’s r) is measured between every region pair. The correlation values are then re-organized into a whole-brain matrix representing the connectivity between every corresponding region pair.

## ICA-based connectivity
The second analysis approach available within RABIES rely on the spatial decomposition of BOLD timeseries using ICA, which models the data as a linear combination of independent sources. In contrast with correlation-based connectivity, which models a single linear relationship between regions, the ICA framework accounts for multiple potentially overlapping sources of BOLD fluctuations, which may further separate confound contributions from connectivity estimates. To obtain individualized connectivity estimates, this analysis framework consists of first deriving ICA components at the group level to define sources, and then recovering individual-specific versions of the sources with dual regression {cite}`Nickerson2017-gq`.
(ICA_target)=
- **Group ICA** (`--group_ica`): RABIES uses FSL’s MELODIC ICA algorithm {cite}`Beckmann2004-yw` to derive ICA components. For group-ICA, timeseries for all scans aligned in commonspace are concatenated to group all data, before computing the ICA decomposition, yielding
$$
Y_{concat} = A\hat{S}
$$
where $Y_{concat}$ are the concatenated timeseries, $\hat{S}$ are the set of spatial maps defining the independent sources, and $A$ is the mixing matrix storing timecourses associated to each component.
(DR_target)=
- **Dual regression (DR)** (`--prior_maps`/`--DR_ICA`) : DR builds on the group ICA decomposition to model scan-specific versions of the group-level components, thus allowing to estimate individualized connectivity for a given brain network first identified through group-ICA {cite}`Beckmann2009-cf,Nickerson2017-gq`. DR consists of two consecutive linear regression steps, where scan-specific timecourses are first derived for each ICA component, and then a scan-specific spatial map is obtained for each component timecourse. Using multivariate linear regression (with Ordinary least square (OLS)), component timecourses are obtained with 
$${\beta}_{TC} = OLS(\hat{S},Y)$$
describing $Y = \hat{S}{\beta}_{TC} + \epsilon$ where $Y$ are the scan timeseries, $\hat{S}$ are the ICA components and ${\beta}_{TC}$ corresponds to the estimated timecourses for each component. To accurately measure connectivity amplitude in the spatial maps derived from DR, the timecourses from the first regression step must be standardized prior to the second regression{cite}`Nickerson2017-gq`. In RABIES, timecourses are thus variance-normalized using an *L2-norm* 
$$
{\beta}^*_{TC} = \frac{{\beta}_{TC}}{||{\beta}_{TC}||_2}
$$
where $||x||_2 = \sqrt{\frac{1}{n}\sum_{i=1}^{n}x_i^2}$ is the *L2-norm* of $x$. The normalized timecourses ${\beta}^*_{TC}$ are then inputed into a second regression step to derive the spatial maps ${\beta}_{SM}$ with 
$${\beta}_{SM} = OLS({\beta}^*_{TC},Y^T)$$
where $Y = {\beta}^*_{TC}{\beta}_{SM} + \epsilon$, thus completing the linear model of the timeseries. The resulting scan-specific spatial maps ${\beta}_{SM}$ will comprise information about network amplitude and shape, which may be compared across subjects or groups with further statistical tests{cite}`Nickerson2017-gq`.
