(metrics_target)=

# Metric definitions

Precise definitions of every quantity RABIES computes and reports. Throughout
this page, the root-mean square (RMS) is
$||x||_2 = \sqrt{\frac{1}{n}\sum_{i=1}^{n}x_i^2}$.

```{seealso}
For where these values appear, see [Output files](outputs.md). For how to
interpret them, see [Data quality assessment](../explanation/data_quality.md).
```

(regressor_target)=

## Nuisance regressors for confound regression

Selected with `--nuisance_regressors` at the confound correction stage.

**mot_6**
: The head motion translation and rotation parameters. Prior to the regression,
  the motion regressors are subjected to the same frame censoring, detrending
  and frequency filtering applied to the BOLD timeseries, to avoid the
  re-introduction of previously corrected confounds, as recommended in
  {cite}`Power2014-yf` and {cite}`Lindquist2019-lq`.

**mot_24**
: The 6 motion parameters together with their temporal derivatives, plus 12
  additional parameters obtained by taking the squared terms — the Friston 24
  parameters {cite}`Friston1996-sa`:

  $$
  mot24_t = [mot6_t,(mot6_t-mot6_{t-1}),(mot6_t)^2,(mot6_t-mot6_{t-1})^2]
  $$

  with $mot24_t$ representing the list of 24 regressors for timepoint $t$. As
  with mot_6, the 24 regressors are additionally subjected to censoring,
  detrending and frequency filtering if applied on BOLD.

**WM/CSF/vascular/global signal**
: The mean signal computed within the corresponding brain mask (WM, CSF,
  vascular or whole-brain) from the partially cleaned timeseries, i.e. after
  confound correction steps 1–4 up to frequency filtering.

**aCompCor_percent**
: Principal component timecourses derived from timeseries within the combined
  WM and CSF masks — the aCompCor technique {cite}`Muschelli2014-vi`. From the
  timeseries within the WM/CSF masks $Y_{WM/CSF}$, a principal component
  analysis (PCA) decomposition is conducted to derive

  $$
  Y_{WM/CSF} = W_{aCompCor}C^T
  $$

  with $C$ a set of spatial principal components and $W$ their associated
  loadings across time. The first components explaining 50% of the variance are
  kept, and their loadings $W_{aCompCor}$ provide the aCompCor nuisance
  regressors. PCA is conducted on the partially cleaned timeseries, i.e. after
  confound correction steps 1–4 up to frequency filtering.

**aCompCor_5**
: As **aCompCor_percent**, but the first 5 components are kept instead of a set
  explaining 50% of the variance.

## Temporal scan diagnosis

(mot6_target)=

**Head motion translation and rotation parameters**
: 3 rotations (Euler angles in radians) and 3 translations (in mm) measured for
  head motion realignment at each timeframe.

(FD_target)=

**Framewise displacement**
: For each timepoint, the displacement — mean across the brain voxels — between
  the current and the next frame. For each brain voxel within the referential
  space for head realignment (the [3D EPI](3D_EPI_target) provided as reference
  for realignment) and for each timepoint, the inverse transform of the head
  motion parameters from the corresponding timepoint is applied to obtain the
  voxel position pre-motion correction. Framewise displacement is then computed
  for each voxel as the Euclidean distance between the pre-motion-correction
  positions for the current and next timepoints. The mean framewise
  displacement $FD_t$ at timepoint $t$ is therefore

  $$
  FD_t = \frac{1}{n}\sum_{i=1}^{n}\sqrt{(x_{i,t+1}-x_{i,t})^2+(y_{i,t+1}-y_{i,t})^2+(z_{i,t+1}-z_{i,t})^2}
  $$

  using the 3D $x$, $y$ and $z$ spatial coordinates in mm for timepoints $t$
  and $t+1$ and voxel indices $i$. Framewise displacement for the last frame,
  which has no future timepoint, is set to 0.

(DVARS_target)=

**DVARS**
: The estimation of temporal shifts in global signal at each timepoint,
  measured as the root-mean-square of the timeseries' temporal derivative

  $$
  DVARS_t = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(Y_{i,t}-Y_{i,t-1})^2}
  $$

  where $Y_{i,t}$ is the BOLD signal in brain voxel $i$ at timepoint $t$. The
  first timepoint is set to 0, having no previous timepoint.

**Framewise distance from mean**
: The mean square error (MSE) between each frame and the average EPI, with the
  average computed as the tri-mean across time voxelwise. The input is the
  resampled timeseries output from the `preprocess` stage.

**Whole-brain/Edge/WM/CSF mask**
: The mean signal across a given brain mask: whole-brain (the global signal),
  WM, CSF or brain edge.

**$CR_{var}$**
: The variance estimated by confound regression at each timepoint, computed as
  $CR_{var} = RMS(Y_{CR})$ across voxels, where $Y_{CR}$ is the
  [predicted confound timeseries](CR_target).

**CR $R^2$**
: The proportion of variance explained, and removed, by confound regression.
  Obtained with $CR_{R^2}= 1-\frac{Var(\hat{Y})}{Var(Y)}$ at each timepoint,
  where $Y$ and $\hat{Y}$ are the timeseries pre- and post-regression, and
  $Var(x) = \frac{1}{n}\sum_{i=1}^{n}(x_i - \mu_x)^2$ is the variance with
  $\mu$ the mean.

**Mean amplitude**
: A set of timecourses averaged as $\frac{1}{n}\sum_{i=1}^{n}|X_i|$, where
  $X_i$ is timecourse $i$. The timecourses correspond to one of:

  - *DR confounds*: timecourses from the first stage of dual regression, using
    the confound components provided to `--prior_confound_idx`
  - *DR networks*: network timecourses from the first stage of dual regression
    as specified with `--prior_bold_idx`
  - *SBC networks*: network timecourses derived from the seeds provided in
    `--seed_list`

## Spatial scan diagnosis

**BOLD<sub>SD</sub>**
: The temporal standard deviation computed for each voxel from the BOLD
  timeseries.

**CR<sub>SD</sub>**
: The temporal standard deviation computed for each voxel from the predicted
  confound timeseries during confound regression, i.e. [$Y_{CR}$](CR_target).

**CR R<sup>2</sup>**
: The proportion of variance explained by confound regression at each voxel.
  Obtained with $CR_{R^2}= 1-\frac{Var(\hat{Y})}{Var(Y)}$ at each voxel, where
  $Y$ and $\hat{Y}$ are the timeseries pre- and post-regression, and
  $Var(x) = \frac{1}{n}\sum_{i=1}^{n}(x_i - \mu_x)^2$ is the variance of $x$
  with $\mu$ the mean.

**Global signal covariance (GS<sub>cov</sub>)**
: The covariance between the global signal and the timeseries at each voxel,
  measured as $GS_{cov} = \frac{1}{n}\sum_{t=1}^{n}Y_t \times GS_t$, where
  $GS_t = \frac{1}{n}\sum_{i=1}^{n}Y_i$ is the mean across all brain voxels for
  a given timepoint.

**DR network X**
: The linear coefficients resulting from the
  [second regression with dual regression](DR_target), corresponding to a
  network amplitude map, for the Xth network specified with `--prior_bold_idx`.

**SBC network X**
: The voxelwise correlation coefficients (Pearson's r) estimated with
  seed-based connectivity, for the Xth seed provided in `--seed_list`.

(dist_plot_metrics)=

## Distribution plot

**Network amplitude**
: The overall network amplitude, summarised by computing the L2-norm across a
  network connectivity map from a subject-level analysis. Such a map can be
  derived from seed-based correlation, or correspond to the linear coefficients
  from the [second regression ${\beta}_{SM}$](DR_target) for dual regression.

**Network specificity**
: The network map (seed-based or dual regression) and the corresponding
  canonical network map are thresholded to include the top X% of voxels with
  highest connectivity, X% being defined by `--brainmap_percent_threshold`, and
  the overlap of the thresholded area is computed using Dice overlap. For dual
  regression, the canonical network map is the original ICA component
  corresponding to that network, provided with `--prior_maps`. For seed-based
  connectivity, the reference network maps are provided using
  `--seed_prior_list`.

**Dual regression confound correlation**
: The timecourse for a single network, from a seed or from dual regression, is
  correlated with the timecourse from each confound component (provided using
  `--prior_confound_idx`) modelled through dual regression. The absolute mean
  correlation is then computed to obtain the average amplitude of confound
  correlations for that network analysis.

**FD-DVARS corr.**
: For each scan, the correlation between the framewise displacement timecourse
  and the DVARS timecourse **computed post-confound correction** — this is not
  the DVARS plotted in the temporal diagnosis figure. Censored timeframes are
  excluded from both timecourses, and DVARS is recomputed after applying
  confound correction, so this metric represents *residual* associations
  between spontaneous motion and the cleaned global signal fluctuations.

**Total $CR_{SD}$**
: The total standard deviation across the
  [predicted confound timeseries $Y_{CR}$](CR_target).

**Mean framewise displacement**
: The mean framewise displacement computed across time, including only frames
  remaining after the censoring applied for confound correction.

**Temporal degrees of freedom**
: The degrees of freedom remaining after confound correction:

  ```text
  tDOF = Original number of timepoints
       - Number of censored timepoints
       - Number of AROMA components removed
       - Number of nuisance regressors
  ```

(group_QC_metrics)=

## Group statistical QC report

**Specificity of network variability**
: As with network specificity in the distribution plot, the network variability
  map and the corresponding canonical network map are thresholded to include
  the top X% of voxels (X% defined by `--brainmap_percent_threshold`), and the
  overlap is estimated using Dice overlap.

**Mean confound correlation**
: For each confound correlation map ($CR_{SD}$, mean FD or tDOF), the mean is
  computed across voxels within the thresholded area of the canonical network
  map, giving a mean correlation within the network's core region.
