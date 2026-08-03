# Workflow reference

RABIES pipelines are built as [Nipype](https://nipype.readthedocs.io/en/latest/)
workflows: each node is a processing step, and the required inputs and outputs
define the links between nodes.

The docstrings below are included directly from the RABIES source, so they
describe the version of the code these docs were built from. For what each
workflow is for and why it exists, see
[The preprocessing workflow](../explanation/preprocessing.md) and
[The confound correction workflow](../explanation/confound_correction.md).

## Preprocessing

(wf_inho_correction)=

### `rabies.preprocess_pkg.inho_correction.init_inho_correction_wf`

[Source](https://github.com/CoBrALab/RABIES/blob/master/rabies/preprocess_pkg/inho_correction.py)

```{literalinclude} ../../rabies/preprocess_pkg/inho_correction.py
:start-after: inho_correction_head_start
:end-before: inho_correction_head_end
```

(wf_commonspace_reg)=

### `rabies.preprocess_pkg.commonspace_reg.init_commonspace_reg_wf`

[Source](https://github.com/CoBrALab/RABIES/blob/master/rabies/preprocess_pkg/commonspace_reg.py)

```{literalinclude} ../../rabies/preprocess_pkg/commonspace_reg.py
:start-after: commonspace_wf_head_start
:end-before: commonspace_wf_head_end
```

(wf_bold_ref)=

### `rabies.preprocess_pkg.bold_ref.init_bold_reference_wf`

[Source](https://github.com/CoBrALab/RABIES/blob/master/rabies/preprocess_pkg/bold_ref.py)

```{literalinclude} ../../rabies/preprocess_pkg/bold_ref.py
:start-after: gen_bold_ref_head_start
:end-before: gen_bold_ref_head_end
```

(wf_hmc)=

### `rabies.preprocess_pkg.hmc.init_bold_hmc_wf`

[Source](https://github.com/CoBrALab/RABIES/blob/master/rabies/preprocess_pkg/hmc.py)

```{literalinclude} ../../rabies/preprocess_pkg/hmc.py
:start-after: hmc_wf_head_start
:end-before: hmc_wf_head_end
```

(wf_motion_params)=

### `rabies.preprocess_pkg.hmc.EstimateMotionParams`

[Source](https://github.com/CoBrALab/RABIES/blob/master/rabies/preprocess_pkg/hmc.py)

```{literalinclude} ../../rabies/preprocess_pkg/hmc.py
:start-after: motion_param_head_start
:end-before: motion_param_head_end
```

(wf_cross_modal_reg)=

### `rabies.preprocess_pkg.registration.init_cross_modal_reg_wf`

[Source](https://github.com/CoBrALab/RABIES/blob/master/rabies/preprocess_pkg/registration.py)

```{literalinclude} ../../rabies/preprocess_pkg/registration.py
:start-after: cross_modal_reg_head_start
:end-before: cross_modal_reg_head_end
```

(wf_bold_resampling)=

### `rabies.preprocess_pkg.resampling.init_bold_preproc_trans_wf`

[Source](https://github.com/CoBrALab/RABIES/blob/master/rabies/preprocess_pkg/resampling.py)

```{literalinclude} ../../rabies/preprocess_pkg/resampling.py
:start-after: bold_resampling_head_start
:end-before: bold_resampling_head_end
```

(wf_mask_resampling)=

### `rabies.preprocess_pkg.resampling.init_mask_preproc_trans_wf`

[Source](https://github.com/CoBrALab/RABIES/blob/master/rabies/preprocess_pkg/resampling.py)

```{literalinclude} ../../rabies/preprocess_pkg/resampling.py
:start-after: mask_resampling_head_start
:end-before: mask_resampling_head_end
```

## Confound correction

(wf_confound_correction)=

### `rabies.confound_correction_pkg.confound_correction.init_confound_correction_wf`

[Source](https://github.com/CoBrALab/RABIES/blob/master/rabies/confound_correction_pkg/confound_correction.py)

```{literalinclude} ../../rabies/confound_correction_pkg/confound_correction.py
:start-after: confound_wf_head_start
:end-before: confound_wf_head_end
```

```{seealso}
[How to contribute to RABIES](../how_to/contribute.md) covers writing a new
Nipype interface and connecting it into one of these workflows.
```
