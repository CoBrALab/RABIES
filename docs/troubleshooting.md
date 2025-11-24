# Troubleshooting

This page provides guidance on common issues encountered when using RABIES.

## Checking Image Orientation with ITK-SNAP

RABIES assumes that input data are properly oriented according to the NIfTI standard (RAS+ orientation).
Incorrectly oriented images are a common source of processing failures and unexpected results.
Before reporting bugs or troubleshooting other issues, it is critical to verify that your images have the correct orientation.

### Installing ITK-SNAP

[ITK-SNAP](https://www.itksnap.org/pmwiki/pmwiki.php) is a free, open-source medical image viewer that properly displays NIfTI orientation information.
RABIES uses ANTs/ITK tools under the hood.
Incorrect image orientation is one of the most common causes of registration failures in RABIES.
Since RABIES expects RAS orientation (Right–Anterior–Superior), you should always verify your anatomical scans in ITK-SNAP before running the pipeline.

### How to Check Orientation in ITK-SNAP

1. **Open your image in ITK-SNAP**
   - Go to **File → Open Main Image…**
   - Load your NIfTI anatomical scan.

2. **Verify anatomical orientation**
   - ITK-SNAP shows three orthogonal views (Axial, Coronal, Sagittal).
   - Ensure anatomical structures appear where you expect them (e.g., nose = anterior, top of head = superior).

3. **Compare your scan to a correctly oriented atlas/template**
   - Open your reference atlas (e.g., SIGMA, Fischer rat, etc.) in another ITK-SNAP window.
   - Verify that structures appear in similar positions and that the orientation labels match.

4. **Check the orientation labels**
   - Inspect axes labels around each view (**R/L**, **A/P**, **S/I**).
   - Confirm they correspond to the real anatomical directions in your scan.
   - Move the cursor: the crosshair should move consistently across all views (e.g., dragging right corresponds to anatomical right).

## See Also

- For registration-specific troubleshooting, see [Registration Troubleshooting](nested_docs/registration_troubleshoot.md)
- For QC-related guidance, see [Preprocessing QC](preproc_QC.md)
- When reporting bugs, refer to the [issue template](https://github.com/CoBrALab/RABIES/blob/master/.github/ISSUE_TEMPLATE/standard-bug-report.md)
