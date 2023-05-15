# Metric definitions

(mot6_target)=
* **Head motion translation and rotation parameters**: Corresponds to 3 rotations (Euler angles in radians) and 3 translations (in mm) measured for head motion realignment at each timeframe.
(FD_target)=
* **Framewise displacement**: For each timepoint, corresponds to the displacement (mean across the brain voxels) between the current and the next frame. For each brain voxel within the referential space for head realignment (i.e. the [3D EPI](3D_EPI_target) which was provided as reference for realignment) and for each timepoint, the inverse transform of the head motion parameters (from the corresponding timepoint) is applied to obtain the voxel position pre-motion correction. Framewise displacement can then be computed for each voxel by computing the Euclidean distance between the positions pre-motion correction for the current and next timepoints. Thus, the mean framewise displacement $FD_t$ at timepoint $t$ is computed as 
$$
FD_t = \frac{1}{n}\sum_{i=1}^{n}\sqrt{(x_{i,t+1}-x_{i,t})^2+(y_{i,t+1}-y_{i,t})^2+(z_{i,t+1}-z_{i,t})^2}
$$
using the 3D $x$,$y$ and $z$ spatial coordinates in mm for timepoints $t$ and $t+1$ and for each voxel indices $i$. Framewise displacement for the last frame (which has no future timepoint) is set to 0.
(DVARS_target)=
* **DVARS**: represents the estimation of temporal shifts in global signal at each timepoint, measured as the root-mean-square of the timeseries’ temporal derivative
$$
DVARS_t = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(Y_{i,t}-Y_{i,t-1})^2}
$$
where $Y_{i,t}$ corresponds to the BOLD signal in brain voxel $i$ at timepoint $t$. The first timepoint is set to 0 (has no previous timepoint).
(tDOF_target)=
* **Temporal degrees of freedom (tDOF)**: The remaining degrees of freedom post-confound correction are calculated as *tDOF = Original number of timepoints - Number of censored timepoints - Number of AROMA components removed - Number of regressors*. 

(regressor_target)=
## Nuisance regressors for confound regression
* **mot_6**: Corresponds to the head motion translation and rotation parameters. Prior to the regression, the motion regressors are also subjected to the same frame censoring, detrending and frequency filtering which were applied to the BOLD timeseries to avoid the re-introduction of previously corrected confounds, as recommend in {cite}`Power2014-yf` and {cite}`Lindquist2019-lq`.
* **mot_24**: Corresponds to the 6 motion parameters together with their temporal derivatives, and 12 additional parameters are obtained by taking the squared terms (i.e. Friston 24 parameters {cite}`Friston1996-sa`)
$$
mot24_t = [mot6_t,(mot6_t-mot6_{t-1}),(mot6_t)^2,(mot6_t-mot6_{t-1})^2]
$$ 
with $mot24_t$ representing the list of 24 regressors for timepoint $t$. As with mot_6, the 24 regressors are additionally subjected to censoring, detrending and frequency filtering if applied on BOLD.
* **WM/CSF/vascular/global signal**: The mean signal is computed within the corresponding brain mask (WM,CSF,vascular or whole-brain mask) from the partially cleaned timeseries (i.e. after the confound correction steps 1-4 up to frequency filtering).
* **aCompCor_percent**: Principal component timecourses are derived from timeseries within the combined WM and CSF masks (aCompCor technique {cite}`Muschelli2014-vi`). From the timeseries within the WM/CSF masks $Y_{WM/CSF}$, a principal component analysis (PCA) decomposition is conducted to derive
$$
Y_{WM/CSF} = W_{aCompCor}C^T
$$
with $C$ corresponding to a set of spatial principal components, and $W$ to their associated loadings across time. The set of first components explaining 50% of the variance are kept, and their loadings $W_{aCompCor}$ provide the set of aCompCor nuisance regressors. PCA is conducted on the partially cleaned timeseries (i.e. after the confound correction steps 1-4 up to frequency filtering).
* **aCompCor_5**: Same as **aCompCor_percent**, but the first 5 components are kept instead of a set explaining 50% of the variance.

