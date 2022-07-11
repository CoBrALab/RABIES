# Confound Correction pipeline

![Processing Schema](pics/confound_correction.png)

The workflow for confound correction is structured following best practices found in human litterature, largely based on recommendations from Power and colleagues{cite}`Power2014-yf`. 
1. `--frame_censoring` : First, spike censoring temporal masks are derived from FD and/or DVARS thresholds. The temporal masking is applied on both BOLD timeseries and nuisance regressors before any other correction step to exclude motion spikes which may bias downstream corrections, in particular, detrending, frequency filtering and confound regression{cite}`Power2014-yf`. 
2. `--match_number_timepoints` : If spike censoring was applied in 1), all scans in the dataset can be set to retain the same final number of frames to account for downstream analysis biases due to unequal tDOF. To do so, a pre-set final number of frames is defined, and the number of extra frames is calculated after the application of spike censoring and taking into account potential edge removal in 5). A corresponding set of extra frames is then randomly selected and removed.
3. `--detrending_order` : Then, linear/quadratic detrending is applied to both BOLD timeseries and nuisance regressors. 
4. `--ica_aroma`: This is followed, if desired, with the cleaning of motion-related sources using ICA-AROMA{cite}`Pruim2015-nm`. The original code for the algorithm (https://github.com/maartenmennes/ICA-AROMA) was adapted to function without the hard-coded human priors for anatomical masking and linear coefficients for classification. ICA-AROMA is applied prior to frequency filtering to remove further effects of motion than can result in ringing after filtering{cite}`Carp2013-uf,Pruim2015-nm`. 
5. `--TR`/`--highpass`/`--lowpass` : Next, frequency filtering requires particular considerations when applied after temporal censoring, which results in missing timepoints. Conventional filters cannot account for missing data, and to address this issue, Power and colleagues introduced a method for simulating missing timepoints{cite}`Power2014-yf`. This method relies on an adaptation of the Lomb-Scargle periodogram, which allows estimating the frequency composition of timeseries with missing data points, and from that estimation, missing timepoints can be simulated while preserving the frequency profile {cite}`Mathias2004-rt`. Missing data is thereby simulated with preserved frequency properties for both BOLD timeseries and nuisance regressors prior to frequency filtering. 
6. `--TR`/`--highpass`/`--lowpass`/`--edge_cutoff` : Following the simulation, frequency filtering (highpass, lowpass or bandpass) is applied using a 3rd-order Butterworth filter, and 30 seconds is removed at each edge of the timeseries to account for edge artefacts following filtering{cite}`Power2014-yf`. The removal of 30 seconds and the use of a 3rd order filter was selected based on the visualization of edge artefacts with simulated data. The nuisance regressors are also filtered to ensure orthogonality between the frequency filters and subsequent confound regression, which can otherwise re-introduce removed confounds{cite}`Lindquist2019-lq`. After frequency filtering, the temporal masks are re-applied to remove the simulated timepoints. 
7. `--conf_list` : After frequency filtering, remaining confound effects are modelled and removed using ordinary least squares (OLS) linear regression of the nuisance regressors onto the BOLD timeseries (i.e. confound regression). 
8. `--image_scaling` : After previous confound correction steps, timeseries can be scaled according to one of three options: 1) based on background noise variance, 2) the global standard deviation calculated by including all voxels, or 3) by applying voxelwise standardization (each voxel is separately scaled to unit variance).
9. `--smoothing_filter` : Timeseries are spatially smoothed using a Gaussian smoothing filter. 

Importantly, each confound correction step (with the exception of linear detrending) is optional when using RABIES. Currently, there are no generalizable confound correction workflow, and rather, the optimal confound correction strategy may well be largely study-specific. The issue of confound correction for resting-state fMRI remains largely unresolved among human literature, and is only beginning to be studied in rodents. Thus, a flexible workflow remains crucial, despite concerns regarding reproducibility.

[[workflow source code](https://github.com/CoBrALab/RABIES/blob/master/rabies/confound_correction_pkg/confound_correction.py)]; **Ref.:** {cite}`Power2012-ji,Power2014-yf,Lindquist2019-lq`
```python
"""
This workflow applies the RABIES confound correction pipeline to preprocessed EPI timeseries. The correction steps are 
orchestrated in line with recommendations from human litterature:   
#1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)
#2 - If --match_number_timepoints is selected, each scan is matched to the defined minimum_timepoint number of frames.
#4 - Linear/Quadratic detrending of fMRI timeseries and nuisance regressors
#4 - Apply ICA-AROMA.
#5 - If frequency filtering and frame censoring are applied, simulate data in censored timepoints using the Lomb-Scargle periodogram, 
        as suggested in Power et al. (2014, Neuroimage), for both the fMRI timeseries and nuisance regressors prior to filtering.
#6 - As recommended in Lindquist et al. (2019, Human brain mapping), make the nuisance regressors orthogonal
        to the temporal frequency filter.
#7 - Apply highpass and/or lowpass filtering on the fMRI timeseries (with simulated timepoints).
#8 - Re-apply the frame censoring mask onto filtered fMRI timeseries and nuisance regressors, taking out the
        simulated timepoints. Edge artefacts from frequency filtering can also be removed as recommended in Power et al. (2014, Neuroimage).
#9 - Apply confound regression using the selected nuisance regressors.
#10 - Scaling of timeseries variance.
#11 - Apply Gaussian spatial smoothing.

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
        random_CR_STD_file_path: variance fitted by random regressors during confound regression
        corrected_CR_STD_file_path: CR_STD_file after substracting the variance fitted by random
            regressors.
        frame_mask_file: CSV file which records which frame were censored
        CR_data_dict: dictionary object storing extra data computed during confound correction
        background_mask_fig: a figure showing the automatically-generated mask of the image background 
            for --image_scaling background_noise.
"""
```
