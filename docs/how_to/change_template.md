# How to override the default common space template

The files that will define the common space template and associated masks are 
controlled by the following parameters at the `preprocess` stage:
- `--anat_template`: the structural image that provides the target for common space
    alignment.
- `--brain_mask`: the brain mask file.
- `--WM_mask`: a mask for the white matter.
- `--CSF_mask`: a mask for the cerebrospinal fluid canals.
- `--vascular_mask`: a mask labelling major blood vessels.

RABIES will automatically input a file for each parameter by default (see [built-in commonspace template and atlas](../reference/template_files.md)),
but you may want to choose a separate template for your study (for instance if you are not inputting adult mouse data). To do so, you will need
to provide the full path to an adequate set of NIfTI-formatted files. At minimum, both `--anat_template` and `--brain_mask` must be provided,
as they are key components of the registration pipeline. The other brain parcellations are optional, but certain downstream options for 
confound correction or generating QC reports will be disabled.

```{important}
The template file `--anat_template` defines the common space coordinates for the entire pipeline. All other 
mask files are expected to overlap with this template file, as assessed from the metadata coordinates of the
NIfTI header. An error will be thrown if a mismatch is detected.
```

```{important}
At minimum, the image orientation of the template and input data must match. It is recommended to follow the
RAS+ orientation convention that RABIES inspects, see [How to check image orientation](./check_orientation.md).
```
