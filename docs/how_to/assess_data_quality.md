# How to assess data quality

This guide covers generating the RABIES data quality reports, using them to
decide which scans to keep, and reporting what you did in a publication. For
what the reports mean and why they exist, see
[Data quality assessment](../explanation/data_quality.md).

The guidance below is written for a standard resting-state fMRI design in which
you compare network connectivity between subjects or groups. The aim is to
identify features of spurious or absent connectivity, remove the scans where
those features dominate, and establish whether the remaining issues confound
your group statistics.

## Generate the reports

Pass `--data_diagnosis` at the analysis stage. It needs a set of ICA components
via `--prior_maps`, with the components corresponding to confounds identified
through `--prior_confound_idx`:

```sh
rabies -p MultiProc analysis confound_correction_outputs/ analysis_outputs/ \
  --data_diagnosis \
  --prior_maps melodic_IC.nii.gz \
  --prior_bold_idx 5 12 19 \
  --prior_confound_idx 0 1 2 6 7 8 \
  --DR_ICA
```

The reports appear in
[`data_diagnosis_datasink/`](../reference/outputs.md#analysis-outputs).

Connectivity can be evaluated for either analysis, or both:

For [dual regression](DR_target)
: Dual regression is always run using the full set of components from
  `--prior_maps`, because several report features are derived from the confound
  components named in `--prior_confound_idx`. Connectivity itself is evaluated for
  each network listed in `--prior_bold_idx`.

For [seed-based connectivity](SBC_target)
: A report is generated for each seed given to `--seed_list`. Each seed must be
  accompanied by a reference network map — a 3D NIfTI file per seed, passed
  through `--seed_prior_list` — representing the connectivity expected for the
  canonical network that seed belongs to.

### Classify your group ICA components

Ideally the components come from the dataset you are analysing, derived with
[group ICA](ICA_target). A
[pre-computed set](https://zenodo.org/records/19069284/files/melodic_IC.nii.gz)
for mice is used by default.

Newly generated components must be inspected visually to identify which
correspond to confound sources. Visualise `group_melodic.ica/melodic_IC.nii.gz`,
or use the FSL report generated automatically in `group_melodic.ica/report`.
Pass the confound components to `--prior_confound_idx` and the networks of interest
to `--prior_bold_idx`.

```{tip}
Classify conservatively. Not every component needs a label — include only
those with a clear feature delineating a network or a confound. The defaults
for `--prior_bold_idx` and `--prior_confound_idx` correspond to the classification
of the pre-computed set, which you can consult as a reference.

For guidance on classifying ICA components in rodents, see
{cite}`Zerbi2015-nl` and {cite}`Desrosiers-Gregoire2024-ou`.
```

## Work through the reports

```{figure} ../pics/QC_framework.png
:alt: The RABIES quality control framework, from scan-level diagnosis to group statistics

The quality control framework: scan-level diagnosis feeds scan inclusion
decisions, which in turn condition the validity of the group-level report.
```

### 1. Inspect each scan

Read the [spatiotemporal diagnosis](diagnosis_target) for every scan. Pay
particular attention to the four main quality markers that define the
[categories of scan quality](quality_marker_target), and judge whether features
of spurious or absent connectivity are prominent.

### 2. Remove scans with spurious or absent connectivity

If those features are prominent in a subset of scans, remove those scans to
mitigate false results. Set thresholds with `--scan_QC_thresholds` on the
scan-level measures of network specificity and confound correlation:

```sh
rabies -p MultiProc analysis confound_correction_outputs/ analysis_outputs/ \
  --data_diagnosis \
  --prior_maps melodic_IC.nii.gz \
  --prior_bold_idx 5 12 19 --prior_confound_idx 0 1 2 6 7 8 --DR_ICA \
  --scan_QC_thresholds '{DR:{Dice:[0.3,0.3,0.3],Conf:[0.25,0.25,0.25],Amp:false}}'
```

The value is a dictionary expression, quoted so the shell leaves it alone. Per
analysis (`DR`, `SBC` or `NPR`) you can set:

`Dice`
: Minimum network detectability, as Dice overlap with the prior. A list of
  values between 0 and 1, matched in order to `--prior_bold_idx` for DR and
  NPR, or to `--seed_list` for SBC. Either give an empty list, or give exactly
  as many thresholds as there are networks.

`Conf`
: Maximum temporal correlation with the dual regression confound timecourses.
  Same list rules as `Dice`.

`Amp`
: `true` to automatically remove scans with outlier network amplitude, which
  can indicate spurious connectivity {cite}`Nickerson2017-gq`.

To pick sensible threshold values, consult the
[distribution plots](dist_plot_target) and the accompanying CSV file, which
gives the measures per scan ID.

```{important}
Scans excluded by `--scan_QC_thresholds` are excluded from the group
statistical report, so the reports must be regenerated after you set the
thresholds.
```

### 3. Check the group level

Consult the [group statistical report](group_stats_target) to identify the main
driver of connectivity variability across scans, and whether it relates
primarily to network activity or to confounds.

### 4. Revisit confound correction if needed

If significant issues remain, redesign the confound correction stage —
see [How to optimise your confound correction strategy](optimise_confound_correction.md).

```{admonition} These guidelines are not prescriptive
:class: caution

They are meant to support identifying analysis pitfalls and improving research
transparency. The judgement of the experimenter is paramount: network
detectability is not always expected, for instance when studying the impact of
anaesthesia or inspecting a visual network in blind subjects. The conversation
about what should constitute proper standards for resting-state fMRI is still
evolving.
```

## Report your quality control in a publication

Every figure in the report is generated as PNG or SVG and can be shared
alongside a publication.

Share, at minimum:

- A spatiotemporal diagnosis for each scan used to derive connectivity results.
- A group statistical report and its affiliated distribution plot, for each
  group or dataset, if the analysis compares connectivity across subjects
  and/or groups.
- The set of ICA components classified as networks and confounds — for example
  the `melodic_IC.nii.gz` file, with its component classification.

If you excluded scans based on the guidelines above, describe the observations
that motivated your criteria and make the associated reports accessible: the
spatiotemporal diagnoses for the scans that displayed spurious or absent
connectivity and motivated a particular `--scan_QC_thresholds` value. If you
designed your confound correction using these tools, report that too.

```{seealso}
- [Data quality assessment](../explanation/data_quality.md) — what these reports are for
- [Metric definitions](../reference/metrics.md) — how every quantity is computed
- [Analysis outputs](../reference/outputs.md#analysis-outputs) — where each file lands
```
