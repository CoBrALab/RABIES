# Scan diagnosis report

Description of the tool; summary of it can be used for; if possible share a link with examples

Unless specified otherwise, all metrics are computed from fMRI timeseries after the confound correction stage.

* Note that if frame censoring was applied, the time axis is discontinued (i.e. there are holes that are not shown.)
* Note that $CR_{var}$ and $CR_{R^2}$ are computing according to the specified list of regressors `--conf_list` during confound correction. If no regressor was specified, $CR_{var}$ and $CR_{R^2}$ are still estimated with a regression of the 6 motion parameters, but the regression isn't applied to remove signal from the timeseries.
* Below the L2-norm corresponds to $||x||_2 = \sqrt{\frac{1}{n}\sum_{i=1}^{n}x_i^2}$

## Temporal features
(image)




## Spatial features
(image)

## Key markers for scan quality categories

* Short description of what the observations are, what does the limited set of features consists of
* Image of 4 examples
* List of expectations for good quality outcome, or table of category X 4 quality features
* Comments on implications of quality outcomes (if good quality, minimal correction needed, if low network detectability not clear recovery post-processing, if confounded explore confound correction optimization)
