# Connectivity analysis

Once confound correction is complete, RABIES estimates resting-state
connectivity using standard analyses: seed-based connectivity, whole-brain
connectivity, group independent component analysis (ICA) and dual regression
(DR).

For every analysis except group ICA, RABIES computes individualised
connectivity maps for each scan separately. These can be exported for
statistical analysis — group comparison and so on — conducted outside RABIES.

```{important}
Every analysis on this page assumes confounds have already been dealt with.
Various fMRI confounds introduce spurious correlations that are
indistinguishable from neural activity in the result, so the quality of a
connectivity estimate is bounded by the quality of the
[confound correction](confound_correction.md) that preceded it. See
[Data quality assessment](data_quality.md).
```

## Correlation-based connectivity

Correlation-based analyses compute a temporal correlation between different
brain regions' BOLD fluctuations to estimate their functional coupling.

(SBC_target)=

### Seed-based connectivity

`--seed_list`

Seed-based connectivity was the first technique developed for mapping
connectivity during rest {cite}`Biswal1995-vh`. The mean timecourse is
extracted from an anatomical seed of interest, and the correlation — Pearson's
r in RABIES — between that timecourse and every other voxel is computed,
producing a correlation map representing the connectivity strength between the
seed and every other brain region.

### Whole-brain connectivity

`--FC_matrix` / `--ROI_type`

An extension of seed-based connectivity to every brain region. Using the
anatomical parcellation provided with the atlas during preprocessing, the seed
timecourse for every parcel is extracted, then the cross-correlation (Pearson's
r) is measured between every region pair. The correlation values are
reorganised into a whole-brain matrix representing the connectivity between
every corresponding region pair.

## ICA-based connectivity

The second approach relies on the spatial decomposition of BOLD timeseries
using ICA, which models the data as a linear combination of independent
sources.

Where correlation-based connectivity models a single linear relationship
between regions, the ICA framework accounts for multiple, potentially
overlapping, sources of BOLD fluctuation. This can further separate confound
contributions from connectivity estimates.

To obtain individualised connectivity estimates, this framework first derives
ICA components at the group level to define the sources, then recovers
individual-specific versions of those sources with dual regression
{cite}`Nickerson2017-gq`.

(ICA_target)=

### Group ICA

`--group_ica`

RABIES uses FSL's MELODIC ICA algorithm {cite}`Beckmann2004-yw` to derive ICA
components. For group ICA, timeseries for all scans aligned in commonspace are
concatenated to group all data before computing the decomposition, yielding

$$
Y_{concat} = A\hat{S}
$$

where $Y_{concat}$ are the concatenated timeseries, $\hat{S}$ are the set of
spatial maps defining the independent sources, and $A$ is the mixing matrix
storing the timecourses associated with each component.

(DR_target)=

### Dual regression

`--prior_maps` / `--DR_ICA`

Dual regression builds on the group ICA decomposition to model scan-specific
versions of the group-level components, allowing individualised connectivity to
be estimated for a brain network first identified through group ICA
{cite}`Beckmann2009-cf,Nickerson2017-gq`.

It consists of two consecutive linear regression steps. First, scan-specific
timecourses are derived for each ICA component; second, a scan-specific spatial
map is obtained for each component timecourse.

Using multivariate OLS linear regression, component timecourses are obtained
with

$${\beta}_{TC} = OLS(\hat{S},Y)$$

describing $Y = \hat{S}{\beta}_{TC} + \epsilon$, where $Y$ are the scan
timeseries, $\hat{S}$ are the ICA components and ${\beta}_{TC}$ are the
estimated timecourses for each component.

To measure connectivity amplitude accurately in the spatial maps derived from
dual regression, the timecourses from the first regression step must be
standardised before the second regression {cite}`Nickerson2017-gq`. RABIES
variance-normalises them using root-mean square (RMS):

$$
{\beta}^*_{TC} = \frac{{\beta}_{TC}}{RMS({\beta}_{TC})}
$$

where $RMS(x) = \sqrt{\frac{1}{n}\sum_{i=1}^{n}x_i^2}$. The normalised
timecourses ${\beta}^*_{TC}$ are then fed into a second regression step to
derive the spatial maps ${\beta}_{SM}$:

$${\beta}_{SM} = OLS({\beta}^*_{TC},Y^T)$$

where $Y = {\beta}^*_{TC}{\beta}_{SM} + \epsilon$, completing the linear model
of the timeseries. The resulting scan-specific spatial maps ${\beta}_{SM}$
carry information about both network amplitude and network shape, which can be
compared across subjects or groups with further statistical tests
{cite}`Nickerson2017-gq`.

```{seealso}
- [Analysis outputs](../reference/outputs.md#analysis-outputs) — where each result is written
- [Metric definitions](../reference/metrics.md) — precise definitions of the derived quantities
- [How to assess data quality](../how_to/assess_data_quality.md) — checking these estimates are trustworthy
```
