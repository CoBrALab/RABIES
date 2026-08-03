(dist_plot_target)=

# The distribution plot

```{figure} ../pics/distribution_plot.png
:alt: Scatter plots of network connectivity measures against confound measures across scans

Each point is a scan. Measures of network connectivity — specificity and
amplitude — are contrasted with measures of confounds across the sample.
```

The distribution plot visualises the distribution of data quality measures
across the dataset, turning the per-scan judgements from the
[spatiotemporal diagnosis](diagnosis_target) into something you can threshold
and report.

Reading the plot:

- Points labelled in **grey** were removed using `--scan_QC_thresholds`. The
  grey dotted lines are the QC thresholds selected for network specificity
  (Dice overlap) and DR confound correlation.
- Among the remaining samples, and for each metric separately, scans presenting
  outlier values are labelled in **orange**. Outliers are detected with a
  modified Z-score threshold, set by `--outlier_threshold` and 3.5 by default.

The derivation of each quality metric is described in the
[metric definitions](dist_plot_metrics).

## What the report is for

The plot serves two distinct purposes.

**Setting scan inclusion criteria.** Inspect that network specificity is
sufficient and that the temporal correlation with confounds (DR confound corr.)
is minimal, then set thresholds for scan inclusion with `--scan_QC_thresholds`.
This is the top right subplot, discussed below.

**Complementing the group statistical report.** Visualise the association
between connectivity and the three confound measures included in the report —
$CR_{SD}$, mean FD and tDOF. This lets you determine, for instance, whether a
group-wise correlation in the [statistical report](group_stats_target) is
driven by a small number of outliers rather than by a dataset-wide effect.

## Scan-level thresholds

```{figure} ../pics/scan_QC_thresholds.png
:alt: Scan quality categories separated along network specificity and confound correlation axes

Reproduced from {cite}`Desrosiers-Gregoire2024-ou`: how the
[categories of scan quality outcome](quality_marker_target) separate along
these two measures.
```

The measures of network specificity (Dice overlap) and temporal correlation
with confounds — where confound timecourses are extracted using the confound
components specified with `--prior_confound_idx` and measured through dual
regression — were defined in {cite}`Desrosiers-Gregoire2024-ou` for conducting
scan-level QC.

They were selected as the measures best suited to quantifying network
detectability and spurious connectivity, and to applying inclusion thresholds
that select scans respecting the assumptions of network detectability and
minimal confound effects.

```{seealso}
[How to assess data quality](../how_to/assess_data_quality.md) gives the
`--scan_QC_thresholds` syntax and the procedure for choosing values.
```
