# Third-party notices

The RABIES license (`LICENSE`, Academic Public License) applies to the RABIES
source code only. RABIES is distributed with third-party software and data that
remain under their own licenses. Those licenses are listed below.

RABIES starts the external tools below as separate programs. It does not link
them into its own process. Their presence alongside RABIES on a disk or in a
container image is aggregation, not combination.

## Corresponding source for GPL components

RABIES conveys the binaries below without changes. The corresponding source is
available at the locations named. This satisfies GNU GPL version 3 section 6(d).

For the AFNI binaries, RABIES elects version 3 of the GNU General Public
License, as permitted by the "or later" term on the GPL-licensed files inside
them. See the AFNI entry below.

| Component | Version | License | Corresponding source |
|---|---|---|---|
| AFNI (`3dDespike`, `3dTshift`, `3dWarp`, `3dAutobox`, `libmri.so`, `libf2c.so`) | see `AFNI_VERSION` below | Public domain and CC-BY-4.0, with GPL-2.0-or-later parts | <https://github.com/afni/afni> |
| MINC Toolkit v2 | 1.9.18-20200813 | GPL-3.0 | <https://github.com/BIC-MNI/minc-toolkit-v2> |
| GNU parallel | as packaged in Ubuntu 22.04 (jammy) | GPL-3.0-or-later | `apt-get source parallel` on Ubuntu 22.04 |
| rsync, patch, perl, imagemagick, curl, unzip, gdebi-core | as packaged in Ubuntu 22.04 (jammy) | various, including GPL | `apt-get source <package>` on Ubuntu 22.04 |

The AFNI build in this image is recorded in the file `AFNI_VERSION` at the image
root. Use that version to select the matching tag or commit in the AFNI
repository.

Most of AFNI is a United States Government work in the public domain, and the
portions copyrighted by the Medical College of Wisconsin are under CC-BY-4.0.
AFNI also bundles files under other licenses. One of them,
`src/niml/niml_md5.c` (Christophe Devine, GPL-2.0-or-later), is compiled into
`libmri`, so the AFNI binaries conveyed here are a combined work under
GPL-2.0-or-later. See <https://github.com/afni/afni/blob/master/LICENSE.txt>.

## FSL

FSL is installed from the FSL conda channel and is **not free for commercial
use**. The FSL license text is at `/opt/conda/LICENCE.FSL` in the image. The FSL
source is at `/opt/conda/src/` in the image, as the FSL conda packages ship it.
For commercial licensing, contact Oxford University Innovation at
<fsl@innovation.ox.ac.uk>, Reference Project 9564.

## ANTs and other Apache-2.0 components

ANTs 2.5.0, `nipype` and `SimpleITK` are under the Apache License 2.0. A copy of
that license is at <https://www.apache.org/licenses/LICENSE-2.0>.

## ICA-AROMA

`rabies/confound_correction_pkg/mod_ICA_AROMA/` is a modified copy of ICA-AROMA
(<https://github.com/maartenmennes/ICA-AROMA>), under the Apache License 2.0.
The license text is kept alongside the code at
`rabies/confound_correction_pkg/mod_ICA_AROMA/license.md`. RABIES changes to
those files are marked with `###RABIES modification` comments.

## Python dependencies

All Python dependencies are under permissive licenses: BSD-3-Clause (numpy,
scipy, pandas, scikit-learn, scikit-image, seaborn, pathos, traits, networkx),
MIT (nibabel, pybids, future), Apache-2.0 (nipype, SimpleITK), the Python
Software Foundation license (matplotlib), MPL-2.0 and MIT (tqdm), and the
Unlicense (qbatch).

## CoBrALab components

`minc-toolkit-extras` and `optimized_antsMultivariateTemplateConstruction` are
CoBrALab software under the same Academic Public License as RABIES. Their
license files are included with them.

## Atlases and data

The DSURQE 40-micron atlas comes from the Mouse Imaging Centre at The Centre for
Phenogenomics. See <https://www.mouseimaging.ca/>. Cite Dorr et al. 2008,
Steadman et al. 2014, Ullmann et al. 2013, Richards et al. 2011, Qiu et al.
2016, and Egan et al. 2015.

Other template and mask files come from Zenodo record 19069284 ("RABIES_extras",
Gabriel Desrosiers-Gregoire, CoBrALab), under CC-BY-4.0.
See <https://zenodo.org/records/19069284>.

## The N3 patch

`patch/nu_estimate_np_and_em.diff` changes a file from N3, copyright 1996 John
G. Sled, McConnell Brain Imaging Centre, Montreal Neurological Institute, McGill
University. That file carries a permissive notice which allows use, copying,
modification and distribution if the copyright notice stays. The patch keeps the
notice.
