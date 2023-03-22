# Metric description

Outputed from preprocessing outputs:
* 6 motion params
* FD
* DVARS
* tSNR

Regressors for confound regression:
(regressor_target)=
* mot6: If applied to BOLD timeseries for confound correction, frame censoring, detrending and frequency filtering are also appplied to the motion regressors prior to confound regression to avoid re-introduction of confounds, as recommend in {cite}`Power2014-yf` and {cite}`Lindquist2019-lq`.

to exclude motion spikes which may bias downstream corrections, in particular, detrending, frequency filtering and confound regression{cite}`Power2014-yf`. 

The nuisance regressors are also filtered to ensure orthogonality between the frequency filters and subsequent confound regression, which can otherwise re-introduce removed confounds{cite}`Lindquist2019-lq`. 

* mot24
* WM/CSF/vascular/global signal: 
* aCompCor_5
* aCompCor_percent


Spatiotemporal diagnosis:

