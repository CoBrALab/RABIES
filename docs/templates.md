# Commonspace template sets

RABIES registers every dataset to a reference atlas, which defines the common space
that outputs are reported in. The template, the masks derived from it, and the
parcellation used at the analysis stage all have to describe the same space, and the
same species. Selecting them one file at a time is error-prone, so they are grouped
into named **template sets**, selected with `--template_set` during preprocessing.

```sh
rabies preprocess bids_inputs/ preprocess_outputs/ --template_set rat
```

The default is `mouse`, which is what RABIES used before this option existed, so
existing commands keep behaving the same way.

## Available sets

| Set | Template | Source |
| --- | --- | --- |
| `mouse` | DSURQE, 40 micron | [Mouse Imaging Centre](https://wiki.mouseimaging.ca/display/MICePub/Mouse+Brain+Atlases) |
| `rat` | SIGMA in-vivo Wistar, 0.15 mm | [SIGMA v2.0](https://github.com/DavidBarriere/SIGMA-Rat-Brain-Templates-and-Atlases) {cite}`Barriere2019-sigma` |

With `--bold_only`, the pipeline registers functional images directly, and each set
switches to an EPI reference template, which is a more robust registration target
than a structural template: `EPICOMMON` for mouse, and the SIGMA in-vivo functional
template for rat.

## What each set provides

Not every set provides every file. A file a set does not provide is not substituted
from another set, since that would register data into the space of a different
species. Instead the operations depending on it are disabled, or refuse to run.

| File | `mouse` | `rat` |
| --- | --- | --- |
| Anatomical template and brain mask | yes | yes |
| White matter and CSF masks | yes | yes, derived by thresholding the SIGMA probabilistic tissue maps |
| Vascular mask | yes | **no** |
| Anatomical parcellation (`--ROI_labels_file`) | yes, 356 ROIs | yes, the Waxholm parcellation normalized to SIGMA, 222 ROIs {cite}`Kleven2023-whs` |
| Group ICA priors (`--prior_maps`) | yes | **no** |
| Pre-drawn seeds (`--seed_list`) | yes, 13 seeds | **no** |

For the rat set this means:

- Nuisance regression cannot use `vascular_signal`. The other regressors, including
  `WM_signal` and `CSF_signal`, are available.
- The analysis stage requires `--prior_maps` to be given explicitly. Dual regression,
  `--data_diagnosis` and the other prior-based analyses cannot run without it. If
  group ICA priors exist for your data, pass them; otherwise run group ICA on your
  own dataset first.
- `--seed_list` only accepts paths to your own seed images, not the pre-drawn names.

## Installing the files

Template files are downloaded into `$XDG_DATA_HOME/rabies`, or `~/.local/share/rabies`
when `XDG_DATA_HOME` is unset.

The two sets are not installed the same way:

- The **mouse** set is installed automatically the first time a run needs it, and is
  baked into the container image.
- The **rat** set is **not**. It has to be installed explicitly:

```sh
rabies install rat
```

Run that from a machine with network access before submitting jobs. Compute nodes on
a cluster frequently have no outbound network, and a run started there with a missing
set fails immediately, naming the command above, rather than partway through
processing.

`rabies install all` installs every set. A set that is already complete is left alone.

## Using your own template

Individual files still override the set, one at a time:

```sh
rabies preprocess bids_inputs/ preprocess_outputs/ --template_set rat \
  --vascular_mask my_vascular_mask.nii.gz
```

Providing `--anat_template` is different: the masks of a set are aligned with that
set's template and cannot be assumed to match a template from elsewhere. Passing your
own template therefore requires `--brain_mask` alongside it, and drops every mask you
do not provide yourself.

## The size check

Registering rat data to a mouse template does not fail. It runs to completion and
produces meaningless results, which is why the mismatch is caught up front: RABIES
compares the size of an input image against the selected template and refuses to run
when they differ by more than the difference between species.

If your data and template do match and the check rejects them anyway, which can happen
with an unusually large field of view, re-run with `--skip_scale_check`.

## Building the rat set from SIGMA

The distributed rat set is derived from the SIGMA release with
`scripts/gen_SIGMA_masks.py`, which binarizes and erodes the probabilistic white
matter and CSF maps and converts the ITK-SNAP label descriptions to CSV. It is run
once when the bundle is built, so that every installation gets identical files, and
is included in the repository so they can be reproduced.

The SIGMA resources are distributed under CC-BY-4.0. Work using the rat template set
should cite {cite}`Barriere2019-sigma` for SIGMA and {cite}`Kleven2023-whs` for the
Waxholm parcellation it embeds.
