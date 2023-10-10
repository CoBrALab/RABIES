# Metric definitions
(metrics_target)=
On this page, the L2-norm corresponds to $||x||_2 = \sqrt{\frac{1}{n}\sum_{i=1}^{n}x_i^2}$

(regressor_target)=
## Nuisance regressors for confound regression
* **mot_6**: Corresponds to the head motion translation and rotation parameters. Prior to the regression, the motion regressors are also subjected to the same frame censoring, detrending and frequency filtering which were applied to the BOLD timeseries to avoid the re-introduction of previously corrected confounds, as recommend in {cite}`Power2014-yf` and {cite}`Lindquist2019-lq`.
* **mot_24**: Corresponds to the 6 motion parameters together with their temporal derivatives, and 12 additional parameters are obtained by taking the squared terms (i.e. Friston 24 parameters {cite}`Friston1996-sa`)
$$
mot24_t = [mot6_t,(mot6_t-mot6_{t-1}),(mot6_t)^2,(mot6_t-mot6_{t-1})^2]
$$ 
with $mot24_t$ representing the list of 24 regressors for timepoint $t$. As with mot_6, the 24 regressors are additionally subjected to censoring, detrending and frequency filtering if applied on BOLD.
* **WM/CSF/vascular/global signal**: The mean signal is computed within the corresponding brain mask (WM,CSF,vascular or whole-brain mask) from the partially cleaned timeseries (i.e. after the confound correction steps 1-4 up to frequency filtering).
* **aCompCor_percent**: Principal component timecourses are derived from timeseries within the combined WM and CSF masks (aCompCor technique {cite}`Muschelli2014-vi`). From the timeseries within the WM/CSF masks $Y_{WM/CSF}$, a principal component analysis (PCA) decomposition is conducted to derive
$$
Y_{WM/CSF} = W_{aCompCor}C^T
$$
with $C$ corresponding to a set of spatial principal components, and $W$ to their associated loadings across time. The set of first components explaining 50% of the variance are kept, and their loadings $W_{aCompCor}$ provide the set of aCompCor nuisance regressors. PCA is conducted on the partially cleaned timeseries (i.e. after the confound correction steps 1-4 up to frequency filtering).
* **aCompCor_5**: Same as **aCompCor_percent**, but the first 5 components are kept instead of a set explaining 50% of the variance.


## Temporal scan diagnosis

(mot6_target)=
* **Head motion translation and rotation parameters**: Corresponds to 3 rotations (Euler angles in radians) and 3 translations (in mm) measured for head motion realignment at each timeframe.
(FD_target)=
* **Framewise displacement**: For each timepoint, corresponds to the displacement (mean across the brain voxels) between the current and the next frame. For each brain voxel within the referential space for head realignment (i.e. the [3D EPI](3D_EPI_target) which was provided as reference for realignment) and for each timepoint, the inverse transform of the head motion parameters (from the corresponding timepoint) is applied to obtain the voxel position pre-motion correction. Framewise displacement can then be computed for each voxel by computing the Euclidean distance between the positions pre-motion correction for the current and next timepoints. Thus, the mean framewise displacement $FD_t$ at timepoint $t$ is computed as 
$$
FD_t = \frac{1}{n}\sum_{i=1}^{n}\sqrt{(x_{i,t+1}-x_{i,t})^2+(y_{i,t+1}-y_{i,t})^2+(z_{i,t+1}-z_{i,t})^2}
$$
using the 3D $x$,$y$ and $z$ spatial coordinates in mm for timepoints $t$ and $t+1$ and for each voxel indices $i$. Framewise displacement for the last frame (which has no future timepoint) is set to 0.
(DVARS_target)=
* **DVARS**: represents the estimation of temporal shifts in global signal at each timepoint, measured as the root-mean-square of the timeseries’ temporal derivative
$$
DVARS_t = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(Y_{i,t}-Y_{i,t-1})^2}
$$
where $Y_{i,t}$ corresponds to the BOLD signal in brain voxel $i$ at timepoint $t$. The first timepoint is set to 0 (has no previous timepoint).
* **Edge/WM/CSF mask**: The L2-norm across voxels within a mask, at each timepoint.
* **$CR_{var}$**: The variance estimated by confound regression is computed for each timepoint. This is done by taking the L2-norm $CR_{var} = ||Y_{CR}||_2$ across voxels at each timepoints, where $Y_{CR}$ is the [predicted confound timeseries](CR_target).
* **CR $R^2$**: Represents the proportion of variance explained (and removed) by confound regression. This is obtained with $CR_{R^2}= 1-\frac{Var(\hat{Y})}{Var(Y)}$ at each timepoint, where $Y$ and $\hat{Y}$ are the timeseries pre- and post-regression respectively, and $Var(x) = \frac{1}{n}\sum_{i=1}^{n}(x_i - \mu_x)^2$ calculates the variance, with $\mu$ as the mean.
* **Mean amplitude**: A set of timecourse are averaged as $\frac{1}{n}\sum_{i=1}^{n}|X_i|$, where $X_i$ is the timecourse $i$. Timecourses can correspond to either of the following sets:
    * DR confounds: timecourses from the first stage of dual regression, using confound components provided to `--prior_confound_idx`.
    * DR networks: network timecourses from the first stage of dual regression as specified with  `--prior_bold_idx`.
    * SBC networks: network timecourses derived from the set of seeds provided in `--seed_list`. 

## Spatial scan diagnosis

* **BOLD<sub>SD</sub>**: The temporal standard deviation is computed for each voxel from the BOLD timeseries.
* **CR<sub>SD</sub>**: The temporal standard deviation computed on each voxel from the predicted confound timeseries during confound regression (i.e. [$Y_{CR}$](CR_target)).
* **CR R<sup>2</sup>**: The proportion of variance explained by confound regression at each voxel. This is obtained with $CR_{R^2}= 1-\frac{Var(\hat{Y})}{Var(Y)}$ at each voxel, where $Y$ and $\hat{Y}$ are the timeseries pre- and post-regression respectively, and $Var(x) = \frac{1}{n}\sum_{i=1}^{n}(x_i - \mu_x)^2$ is variance of $x$, with $\mu$ as the mean.
* **Global signal covariance (GS<sub>cov</sub>)**: The covariance between the global signal and the timeseries at each voxel is measured as $GS_{cov} = \frac{1}{n}\sum_{t=1}^{n}Y_t \times GS_t$, where $GS_t = \frac{1}{n}\sum_{i=1}^{n}Y_i$, i.e. the mean across all brain voxels for a given timepoint.
* **DR network X**: The linear coefficients resulting from the [second regression with dual regression](DR_target), corresponding to a network amplitude map (for the Xth network specified for analysis with `--prior_bold_idx`).
* **SBC network X**: The voxelwise correlation coefficients (pearson's r) estimated with seed-based connectivity (for the Xth seed provided for analysis with `--seed_list`).


## Distribution plot
(dist_plot_metrics)=

* **Network amplitude (not computed for seed-based connectivity)**: The overall network amplitude is summarized by computing the L2-norm across the network connectivity map outputed from dual regression (i.e. linear coefficients from the [second regression ${\beta}_{SM}$](DR_target))
* **Network specificity**: The network map (seed-based or dual regression) and the corresponding canonical network map are thresholded to include the top 4% of voxels with highest connectivity, and the overlap of the thresholded area is computed using Dice overlap. For dual regression, the 'canonical network' map will consist of the original ICA component corresponding to that network provided with `--prior_maps`, and for seed-based connectivity, the reference network maps are provided using the `--seed_prior_list` parameter.
* **Dual regression confound correlation**: The timecourse for a single network (from a seed or dual regression) is correlated with the timecourse from each confound component (provided using `--prior_confound_idx`) modelled through dual regression, then the absolute mean correlation is computed to obtain the average amplitude of confound correlations for this specific network analysis.
* **Total $CR_{SD}$**: The total standard deviation across the [predicted confound timeseries $Y_{CR}$](CR_target).
* **Mean framewise displacement**: The mean framewise displacement computed across time (only including frames after censoring applied for confound correction).
* **Temporal degrees of freedom**: The remaining degrees of freedom post-confound correction are calculated as `tDOF = Original number of timepoints - Number of censored timepoints - Number of AROMA components removed - Number of nuisance regressors`.

## Group statistical QC report
(group_QC_metrics)=
* **Specificity of network variability**: similarly to network specificity in the distribution plot, the network variability map  and the corresponding canonical network map are thresholded to include the top 4% of voxels, and then the overlap is estimated using Dice overlap.
* **Mean confound correlation**: for each confound correlation map (either $CR_{SD}$, mean FD or tDOF), the mean is computed across voxels included within the thresholded area of the canonical network map, to obtain a mean correlation within the network's core region.