---
name: Standard bug report
about: Standard template to report bugs that are encountered
title: ""
labels: ""
assignees: ""
---

**Before submitting, please check our [Troubleshooting Guide](https://github.com/CoBrALab/RABIES/blob/master/docs/troubleshooting.md) for common issues and solutions.**

**Have you checked original image orientation?**
RABIES assumes properly oriented data.
Please ensure image orientation is correct by visualizing images with [ITK-snap](https://www.itksnap.org/pmwiki/pmwiki.php).
See the [troubleshooting guide](https://github.com/CoBrALab/RABIES/blob/master/docs/troubleshooting.md#checking-image-orientation-with-itk-snap) for detailed instructions on verifying image orientation.

**Describe the bug**
A clear and concise description of what the bug is.

**Describe your system**

- RABIES version (e.g., release number or container tag)
- Execution method (e.g., Docker, Singularity/Apptainer, or local install)
- Operating system (e.g., Ubuntu 22.04, macOS 13.5, Windows WSL2)
- (Optional) CPU threads and RAM available; container/runtime version

**Describe RABIES call**
Include a copy of the command call you executed from the terminal, and any additional information that could be relevant to your execution.

**Attach log file**
Attach to your issue the .log files present in the output folder. (e.g. rabies_out/rabies_preprocess.log)

**Attach QC_report**
Attach to your issue the QC_report folder present in the output folder. If QC_report is too large, mention that it couldn't be shared, and if possible, provide an alternative access to the files.

**Is the bug reproducible?**

**Additional context**
Add any other context about the problem here.
