# Confound Correction pipeline

![Processing Schema](pics/confound_correction.png)

The workflow for confound correction is structured following best practices found in human litterature, largely based on recommendations from Power and colleagues (Power et al., 2014). 
1. `--FD_censoring`/`--FD_threshold`/`--DVARS_censoring` : First, frame censoring temporal masks are derived from FD and/or DVARS thresholds. The temporal masking is applied on both BOLD timeseries and nuisance regressors before any other correction step to exclude motion spikes which may bias downstream corrections, in particular, detrending, frequency filtering and confound regression (Power et al., 2014). 
2. Then, linear detrending is applied to both BOLD timeseries and nuisance regressors. 
3. `--run_aroma`/`--aroma_dim`/`--aroma_random_seed` : This is followed, if desired, with the cleaning of motion-related sources using ICA-AROMA (Pruim et al., 2015). The original code for the algorithm (https://github.com/maartenmennes/ICA-AROMA) was adapted to function without the hard-coded human priors for anatomical masking and linear coefficients for classification. ICA-AROMA is applied prior to frequency filtering to remove further effects of motion than can result in ringing after filtering (Carp, 2013; Pruim et al., 2015). 
4. `--TR`/`--highpass`/`--lowpass` : Next, frequency filtering requires particular considerations when applied after temporal censoring, which results in missing timepoints. Conventional filters cannot account for missing data, and to address this issue, Power and colleagues introduced a method for simulating missing timepoints (Power et al., 2014). This method relies on an adaptation of the Lomb-Scargle periodogram, which allows estimating the frequency composition of timeseries with missing data points, and from that estimation, missing timepoints can be simulated while preserving the frequency profile (Mathias et al., 2004). Missing data is thereby simulated with preserved frequency properties for both BOLD timeseries and nuisance regressors prior to frequency filtering. 
5. `--TR`/`--highpass`/`--lowpass`/`--edge_cutoff` : Following the simulation, frequency filtering (highpass, lowpass or bandpass) is applied using a 3rd-order Butterworth filter, and 30 seconds is removed at each edge of the timeseries to account for edge artefacts following filtering (Power et al., 2014). The removal of 30 seconds and the use of a 3rd order filter was selected based on the visualization of edge artefacts with simulated data. The nuisance regressors are also filtered to ensure orthogonality between the frequency filters and subsequent confound regression, which can otherwise re-introduce removed confounds (Lindquist et al., 2019). After frequency filtering, the temporal masks are re-applied to remove the simulated timepoints. 
6. `--conf_list` : After frequency filtering, remaining confound effects are modelled and removed using closed-form linear regression of the nuisance regressors onto the BOLD timeseries (i.e. confound regression). 
7. `--standardize` : Timeseries are temporally standardized (mean-substrated, variance-normalized)
8. `--smoothing_filter` : Timeseries are spatially smoothed using a Gaussian smoothing filter. 

Importantly, each confound correction step (with the exception of linear detrending) is optional when using RABIES. Currently, there are no generalizable confound correction workflow, and rather, the optimal confound correction strategy may well be largely study-specific. The issue of confound correction for resting-state fMRI remains largely unresolved among human literature, and is only beginning to be studied in rodents. Thus, a flexible workflow remains crucial, despite concerns regarding reproducibility.

A single Nipype workflow handles confound correction within RABIES: https://github.com/CoBrALab/RABIES/blob/master/rabies/confound_correction_pkg/confound_correction.py 
```python
"""
This workflow applies the RABIES confound correction pipeline to preprocessed EPI timeseries. The correction steps are 
orchestrated in line with recommendations from human litterature:   
#1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)
#2 - Linear detrending of fMRI timeseries and nuisance regressors
#3 - Apply ICA-AROMA.
#4 - If frequency filtering and frame censoring are applied, simulate data in censored timepoints using the Lomb-Scargle periodogram, 
        as suggested in Power et al. (2014, Neuroimage), for both the fMRI timeseries and nuisance regressors prior to filtering.
#5 - As recommended in Lindquist et al. (2019, Human brain mapping), make the nuisance regressors orthogonal
        to the temporal frequency filter.
#6 - Apply highpass and/or lowpass filtering on the fMRI timeseries (with simulated timepoints).
#7 - Re-apply the frame censoring mask onto filtered fMRI timeseries and nuisance regressors, taking out the
        simulated timepoints. Edge artefacts from frequency filtering can also be removed as recommended in Power et al. (2014, Neuroimage).
#8 - Apply confound regression using the selected nuisance regressors.
#9 - Standardize timeseries
#10 - Apply Gaussian spatial smoothing.

References:
    Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic 
        correlations in functional connectivity MRI networks arise from subject motion. Neuroimage, 59(3), 2142-2154.
    Power, J. D., Mitra, A., Laumann, T. O., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2014). Methods to detect, 
        characterize, and remove motion artifact in resting state fMRI. Neuroimage, 84, 320-341.
    Lindquist, M. A., Geuter, S., Wager, T. D., & Caffo, B. S. (2019). Modular preprocessing pipelines can reintroduce 
        artifacts into fMRI data. Human brain mapping, 40(8), 2358-2376.

Workflow:
    parameters
        cr_opts: command line interface parameters from confound_correction

    inputs
        bold_file: preprocessed EPI timeseries
        brain_mask: brain mask overlapping with EPI timeseries
        csf_mask: CSF mask overlapping with EPI timeseries
        confounds_file: CSV file with nuisance timecourses
        FD_file: CSV file with the framewise displacement

    outputs
        cleaned_path: the cleaned EPI timeseries
        aroma_out: folder with outputs from ICA-AROMA
        VE_file: variance explained (R^2) from confound regression at each voxel
        STD_file: standard deviation on the cleaned EPI timeseries
        CR_STD_file: standard deviation on the confound timeseries modelled during confound regression
        frame_mask_file: CSV file which records which frame were censored
        CR_data_dict: dictionary object storing extra data computed during confound correction
"""
```

### References
- Carp, J. (2013). Optimizing the order of operations for movement scrubbing: Comment on Power et al [Review of Optimizing the order of operations for movement scrubbing: Comment on Power et al]. NeuroImage, 76, 436–438.
- Lindquist, M. A., Geuter, S., Wager, T. D., & Caffo, B. S. (2019). Modular preprocessing pipelines can reintroduce artifacts into fMRI data. Human Brain Mapping, 40(8), 2358–2376.
- Mathias, A., Grond, F., Guardans, R., Seese, D., Canela, M., & Diebner, H. H. (2004). Algorithms for spectral analysis of irregularly sampled time series. Journal of Statistical Software, 11(1), 1–27.
- Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage, 59(3), 2142–2154.
- Power, J. D., Mitra, A., Laumann, T. O., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2014). Methods to detect, characterize, and remove motion artifact in resting state fMRI. NeuroImage, 84, 320–341.
- Pruim, R. H. R., Mennes, M., van Rooij, D., Llera, A., Buitelaar, J. K., & Beckmann, C. F. (2015). ICA-AROMA: A robust ICA-based strategy for removing motion artifacts from fMRI data. NeuroImage, 112, 267–277.

