# How to optimise your confound correction strategy

(optim_CR)=

There is no single correct confound correction strategy. The right one is
dataset-specific, and the way to find it is to start minimal and add one
correction at a time, checking after each addition whether the data quality
reports improved.

This protocol comes from {cite}`Desrosiers-Gregoire2024-ou`. It assumes you can
already generate and read the data quality reports — if not, start with
[How to assess data quality](assess_data_quality.md).

```{important}
Start minimal and stay minimal for as long as the reports allow. Excessive
correction removes network activity along with the confounds, and
over-correction is harder to detect after the fact than under-correction.
```

## The protocol

1. **Start with a minimal correction** and generate the data quality reports at
   the analysis stage. A reasonable minimum is frame censoring on framewise
   displacement, regression of the 6 motion parameters, and spatial smoothing:

   ```sh
   rabies -p MultiProc confound_correction preprocess_outputs/ confound_correction_outputs/ \
     --frame_censoring FD_censoring=true,FD_threshold=0.05 \
     --nuisance_regressors mot_6 \
     --smoothing_filter 0.3
   ```

2. **Evaluate the reports**, following
   [How to assess data quality](assess_data_quality.md).

3. **Choose one additional correction**, using the table below to match what
   you observed to the correction that addresses it.

4. **Re-run confound correction with that one correction added**, regenerate
   the reports, and compare. Keep the addition only if it improved the quality
   outcomes. Adding one correction at a time is what makes its effect
   attributable.

5. **Repeat steps 3 and 4** until the quality outcomes are acceptable or you
   have run out of options.

## Matching observations to corrections

```{figure} ../pics/CR_optimization_table.svg
:alt: Table relating data quality observations to the corresponding confound correction options
:width: 100%

Guidance for prioritising additional corrections based on observations from the
data quality reports.
```

```{seealso}
[The confound correction workflow](confound_pipeline_target) describes every
correction step available and the order in which RABIES applies them.
```
