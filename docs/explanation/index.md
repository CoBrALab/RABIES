# Explanation

Background and discussion of how RABIES works and why it was built this way.
These pages are for understanding, not for following along — read them away
from the keyboard.

## The processing stages

RABIES is structured into three stages that run in sequence, each consuming the
previous one's output.

```{toctree}
---
maxdepth: 1
---
preprocessing
confound_correction
analysis
```

## Assessing data quality

Confound correction cannot be designed blind. RABIES generates a set of reports
for characterising data quality and its impact on connectivity estimates.

```{toctree}
---
maxdepth: 1
---
data_quality
scan_diagnosis
distribution_plot
group_statistics
```

```{seealso}
The full methodological account is in the RABIES publication,
{cite}`Desrosiers-Gregoire2024-ou`.
```
