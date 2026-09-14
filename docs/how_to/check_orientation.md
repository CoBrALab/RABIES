# How to check image orientation

RABIES assumes input data is oriented according to the NIfTI standard (RAS+,
i.e. Right–Anterior–Superior). Incorrectly oriented images are one of the most
common causes of registration failure and of results that look wrong for no
apparent reason, because the ANTs/ITK tools RABIES calls will happily register
a mis-oriented brain to the template and produce nonsense.

```{important}
Check orientation before you report a bug, and before you start
troubleshooting anything else. It costs five minutes and rules out the most
likely explanation.
```

## Check with ITK-SNAP

[ITK-SNAP](https://www.itksnap.org/pmwiki/pmwiki.php) is a free, open-source
medical image viewer that reads NIfTI orientation information correctly. Other
viewers do not always, which is why it is the recommended tool here.

1. **Open your image.** *File → Open Main Image…*, and load your NIfTI
   anatomical scan.

2. **Check the anatomy is where you expect it.** ITK-SNAP shows axial, coronal
   and sagittal views. The nose should be anterior, the top of the head
   superior.

3. **Compare against your reference template.** Open the atlas you intend to
   pass to `--anat_template` (SIGMA, DSURQE, Fischer rat, whichever applies) in
   a second ITK-SNAP window. Structures should appear in similar positions, and
   the orientation labels should match.

4. **Check the orientation labels.** Each view is labelled **R/L**, **A/P** and
   **S/I**. Confirm each corresponds to the real anatomical direction in your
   scan, then move the cursor: the crosshair should track consistently across
   all three views, so that dragging towards the label **R** moves towards
   anatomical right.

If any of these disagree, the orientation stored in your NIfTI header does not
describe your data, and you need to fix it before running RABIES. This is
usually introduced during conversion from the scanner format.

```{seealso}
For Bruker data, [BrkRaw](https://brkraw.github.io/) handles the raw-to-NIfTI
conversion; the CoBrALab maintains
[notes on the conversion](https://github.com/CoBrALab/documentation/wiki/bruker2nifti-conversion).
```

## If RABIES reports that a file does not overlap with the template

RABIES stops with `The file ... does not appear to overlap with provided template` when
a mask, labels file or seed does not share the origin and orientation of the template.
When mixing your own files with a template set, every file has to be aligned with the
template actually in use, not with the one it was originally drawn on.

The check compares only the image origin and the orientation of its axes. It does not
compare voxel size, image dimensions or field of view, so a file can pass and still not
line up with the template. Open the file over the template in ITK-SNAP to confirm it
sits where it should, as described in [Check with ITK-SNAP](#check-with-itk-snap).

Masks and labels files are converted to the RAS orientation before the check, so a file
whose axes are only flipped relative to the template passes. Seeds given to
`--seed_list` are not converted. The usual causes of this error are:

- a file drawn on a different template, which usually has a different origin;
- a file cropped or padded by another tool, which moves its origin;
- an oblique image, whose axes are rotated rather than flipped;
- a seed saved in a different orientation than the template.

Earlier versions of RABIES did not catch orientation mismatches in this check, so files
that used to be accepted may now stop a run.

## If orientation is correct and registration still fails

Move on to [How to troubleshoot registration](troubleshoot_registration.md),
which covers tuning inhomogeneity correction, masking and registration stages.
Read the [preprocessing QC report](../reference/qc_outputs.md) first, to
establish which registration step is the one that failed.

When you do report a problem, follow the
[issue template](https://github.com/CoBrALab/RABIES/blob/master/.github/ISSUE_TEMPLATE/standard-bug-report.md).
