# Command line interface

RABIES runs from a terminal. The software is divided into three processing
stages, and one of them must be selected on every invocation:

```text
rabies [execution options] {preprocess,confound_correction,analysis} [stage options] input output
```

```{list-table}
:header-rows: 1
:widths: 25 25 50

* - Stage
  - Takes as input
  - Produces
* - `preprocess`
  - A BIDS directory
  - Motion-corrected, distortion-corrected timeseries aligned to commonspace
* - `confound_correction`
  - A `preprocess` output directory
  - Cleaned timeseries
* - `analysis`
  - A `confound_correction` output directory
  - Connectivity estimates and data quality reports
```

Execution options — parallel processing, memory management, scan selection —
are given **before** the stage name. Stage options are given after it.

```{note}
The `--help` output below is generated at documentation build time from the
RABIES version these docs were built for. Run `rabies --help` locally to see
the options for the version you have installed.
```

## Argument syntax

RABIES uses diverse argument shapes:

Flags
: `--apply_STC`. Present or absent, no value.

Single values
: `--TR 1.2`. The value follows the parameter name after a space.

Lists
: `--nuisance_regressors WM_signal CSF_signal mot_6`. Values follow the
  parameter name, separated by spaces.

Key-value groups
: `--commonspace_reg masking=true,template_registration=SyN`. One parameter
  takes several settings as comma-separated `key=value` pairs, with no spaces.

## `rabies --help`

Execution options that apply to every stage: parallel processing plugin, thread
and memory limits, scan inclusion and exclusion, output data type and
interpolation.

:::{dropdown} rabies --help
:icon: terminal

```{program-output} rabies --help
```
:::

## `rabies preprocess --help`

Input selection, image corrections, registration options, resampling, slice
timing correction, and the reference atlas files.

:::{dropdown} rabies preprocess --help
:icon: terminal

```{program-output} rabies preprocess --help
```
:::

## `rabies confound_correction --help`

Frame censoring, detrending, ICA-AROMA, frequency filtering, confound
regression, intensity scaling and smoothing.

:::{dropdown} rabies confound_correction --help
:icon: terminal

```{program-output} rabies confound_correction --help
```
:::

## `rabies analysis --help`

Seed-based connectivity, whole-brain connectivity matrices, group ICA, dual
regression, and the `--data_diagnosis` quality reports.

:::{dropdown} rabies analysis --help
:icon: terminal

```{program-output} rabies analysis --help
```
:::
