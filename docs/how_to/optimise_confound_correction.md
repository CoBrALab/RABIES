(optim_CR)=

# How to optimise your confound correction strategy

There is no agreed-upon single confound correction strategy for fMRI
functional connectivity analysis. Designing a confound correction 
pipeline is navigating a central trade-off: too little signal removal
risks leaving uncorrected artefacts that corrupt downstream analyses,
excessive signal removal risks removing signal of interest relating to 
network activity. The ideal strategy will depend on the extent and 
nature of artefacts present in a given dataset.

Here is described a protocol introduced by {cite}`Desrosiers-Gregoire2024-ou`
to navigate these decisions on a per-dataset basis. It assumes you can
already generate and read the data quality reports — if not, start with
[How to assess data quality](assess_data_quality.md).

## The protocol

1. **Start with a minimal correction** and generate the data quality reports at
   the analysis stage. A reasonable minimum can be frame censoring on framewise
   displacement, regression of the 6 motion parameters, and spatial smoothing:

   ```sh
   rabies -p MultiProc confound_correction preprocess_outputs/ confound_correction_outputs/ \
     --frame_censoring FD_censoring=true,FD_threshold=0.05 \
     --nuisance_regressors mot_6 \
     --smoothing_filter 0.3
   ```

```{Advice}
Starting minimal is ideal, as excessive
correction removes network activity along with the confounds, and
over-correction is harder to detect after the fact than under-correction.
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
