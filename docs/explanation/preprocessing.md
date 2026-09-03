# The preprocessing workflow

```{figure} ../pics/preprocessing.png
:alt: Diagram of the RABIES preprocessing workflow

The RABIES preprocessing workflow.
```

Preprocessing fMRI scans prior to analysis requires, at minimum, anatomical
alignment of scans to a common space, head realignment to correct for motion,
and correction of the susceptibility distortions arising from the echo-planar
imaging (EPI) acquisition of functional scans. The core preprocessing pipeline
in RABIES carries out each of these with state-of-the-art processing tools and
techniques.

## Alignment to the common space

Structural images, acquired alongside the EPI scans, are initially corrected
for inhomogeneities (**Structural inhomogeneity correction**) and then
registered together so that different MRI acquisitions can be aligned.

That registration works by generating an unbiased, data-driven template
(**Unbiased template generation**) through the iterative non-linear
registration of each image to the dataset consensus average, where the average
is updated at each iteration to provide an increasingly representative dataset
template ([optimized_antsMultivariateTemplateConstruction](https://github.com/CoBrALab/optimized_antsMultivariateTemplateConstruction);
{cite}`Avants2011-av`).

The finalised template, after the last iteration, provides a representative
alignment of each MRI session to a template sharing the acquisition properties
of the dataset — brain shape, field of view, anatomical contrast — which makes
it a stable registration target for cross-subject alignment. This
newly-generated unbiased template is then itself registered to an external
reference atlas (**Atlas registration**), which supplies both an anatomical
segmentation and a common space comparable across studies.

```{note}
This is why RABIES builds a study-specific template rather than registering
every scan directly to the atlas: the intermediate target resembles your data,
so each individual registration has less work to do and is less likely to
fail. The cost is computation time, which is why `fast_commonspace=true`
exists to skip it.
```

(3D_EPI_target)=

## EPI motion and distortion corrections

A volumetric EPI image is first derived using a trimmed mean across the EPI
frames, after an initial motion realignment step (**3D EPI generation**). Using
this volumetric EPI as a target, the head motion parameters are estimated by
realigning each EPI frame to the target with a rigid registration
(**Head motion estimation**).

To correct EPI susceptibility distortions, the volumetric EPI is first
subjected to an inhomogeneity correction step
(**Functional inhomogeneity correction**), then registered non-linearly to the
anatomical scan from the same MRI session, which yields the geometrical
transforms required to recover brain anatomy {cite}`Wang2017-ci`
(**Susceptibility distortion estimation**).

## Derivation of preprocessed EPI timeseries

The transforms providing head motion correction, susceptibility distortions 
and alignment to the common space are concatenated into a single resampling operation 
— avoiding multiple resampling — applied at each EPI frame {cite}`Esteban2019-rs` (**Frame-wise resampling**).
This generates the preprocessed EPI timeseries in common space, while alternatively 
the transforms to common space can be dropped to generate instead native space timeseries
using the `--resampling_space` parameter.

```{important}
Concatenating the transforms matters. Resampling an image is lossy, so
applying motion correction and then distortion correction as two separate
resampling steps blurs the data twice. RABIES composes the transforms first
and resamples once.
```

## Working without structural scans

Structural scans are recommended but not required. An alternative workflow,
selected with `--bold_only`, preprocesses an input dataset containing only EPI
functional images.

In this workflow the volumetric EPI corrected during
**Functional inhomogeneity correction** replaces the structural image for the
purpose of common space alignment, and is used to generate the unbiased
template, which is in turn registered to the reference atlas. This final
registration to the atlas accounts for the estimation of susceptibility
distortions, in place of the registration to a same-session structural image.

```{note}
When using the RABIES default mouse atlas, `--bold_only` also switches the
default template to an EPI reference template, which is a more robust target
for EPI registration than a structural reference.
```

```{seealso}
- [Preprocessing QC outputs](../reference/qc_outputs.md) — how to verify each step succeeded
- [How to troubleshoot registration](../how_to/troubleshoot_registration.md) — what to do when one did not
- [Workflow reference](../reference/workflows.md) — the source docstrings for every module above
```
