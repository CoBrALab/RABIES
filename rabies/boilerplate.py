def define_registration(string):
    if string=='SyN':
        return 'nonlinear'
    elif string=='Affine':
        return 'affine'
    elif string=='Rigid':
        return 'rigid'


def preprocess_boilerplate(opts):
    methods=''
    references={}
    i=1

    # define references
    rabies_ref='Desrosiers-Gregoire et al., in preparation.'
    afni='Cox, R. W. (1996). AFNI: software for analysis and visualization of functional \
magnetic resonance neuroimages. Computers and Biomedical research, 29(3), 162-173.'
    avants_2011='Avants, B. B., Tustison, N. J., Song, G., Cook, P. A., Klein, A., & Gee, J. C. (2011). A reproducible evaluation of ANTs similarity metric performance in brain image registration. NeuroImage, 54(3), 2033–2044.'
    sdc_nonlinear='Wang, S. et al. Evaluation of field map and nonlinear registration methods for correction of susceptibility artifacts in diffusion MRI. Front.Neuroinform. 11, 17 (2017).'
    fmriprep='Esteban, O., Markiewicz, C. J., Blair, R. W., Moodie, C. A., Isik, A. I., Erramuzpe, A., Kent, J. D., Goncalves, M., DuPre, E., Snyder, M., Oya, H., Ghosh, S. S., Wright, J., Durnez, J., Poldrack, R. A., & Gorgolewski, K. J. (2019). fMRIPrep: a robust preprocessing pipeline for functional MRI. Nature Methods, 16(1), 111–116.'

    # citing RABIES
    references[rabies_ref]=i
    i+=1
    methods+=f"  The preprocessing of fMRI images was conducted using the open-source RABIES software (https://github.com/CoBrALab/RABIES)[{references[rabies_ref]}]. "

    # Applying autobox
    if opts.bold_autobox or opts.anat_autobox:
        references[afni]=i
        i+=1
        autobox=f"extra-space around the brain was automatically cropped (3dAutobox, AFNI)[{references[afni]}]. "
        if opts.bold_autobox and opts.anat_autobox:
            autobox="For both the anatomical and functional images, "+autobox
        elif opts.bold_autobox and not opts.anat_autobox:
            autobox="For the functional images, "+autobox
        elif not opts.bold_autobox and opts.anat_autobox:
            autobox="For the anatomical images, "+autobox
        methods+=autobox

    # Despiking
    if opts.apply_despiking:
        if not afni in references.keys():
            references[afni]=i
            i+=1

        methods+=f"Temporal spikes were corrected for at each voxel (3dDespike, AFNI)[{references[afni]}]. "

    # Detect dummy
    if opts.detect_dummy:
        methods+="Dummy scans were automatically detected and removed from each EPI. "
        detect_dummy="If dummy scans are detected, the median of these volumes provides a volumetric EPI image as reference, given their \
higher anatomical contrast. Otherwise, a"
    else:
        detect_dummy="A"

    # Generate 3D EPI ref
    methods+=f"{detect_dummy} volumetric EPI image was derived using a trimmed mean across the EPI frames, after an initial motion realignment step. "

    # HMC
    methods+="Using this volumetric EPI as a target, the head motion parameters are estimated by realigning each EPI frame to the target using a rigid registration. "

    # Slice MC
    if opts.apply_slice_mc:
        methods+="--apply_slice_mc NOT IMPLEMENTED YET "

    # common space alignment
    if not opts.bold_only:
        reg_image='structural images, which were acquired along the EPI scans,'
    else:
        reg_image='volumetric EPI scans'
    methods+=f"To conduct common space alignment, {reg_image} are initially corrected for \
inhomogeneities and "
    if not opts.fast_commonspace:
        references[avants_2011]=i
        i+=1
        methods+=f"then registered together to allow the alignment of different MRI acquisitions. This registration is conducted by generating \
an unbiased data-driven template through the iterative nonlinear registration of each image to the dataset consensus average, where the average \
gets updated at each iteration to provide an increasingly representative dataset template \
(https://github.com/CoBrALab/optimized_antsMultivariateTemplateConstruction) [{references[avants_2011]}]. The finalized template after the last \
iteration provides a representative alignment of each MRI session to a template that shares the acquisition properties of the dataset \
(e.g. brain shape, FOV, anatomical contrast, ...), which makes it a stable registration target for cross-subject alignment. After aligning the MRI \
sessions, this newly-generated unbiased template is then itself registered, using a {define_registration(opts.atlas_reg_script)} registration, to an external reference atlas to provide both an anatomical segmentation \
and a common space comparable across studies defined from the provided reference atlas."
    else:
        # Fast commonspace alignment
        methods+=f"then registered individually to an external reference atlas, with a {define_registration(opts.atlas_reg_script)} registration, to \
provide both an anatomical segmentation and a common space comparable across studies defined from the provided reference atlas. "

    # STC
    if opts.apply_STC:
        if not afni in references.keys():
            references[afni]=i
            i+=1
        stc_str=f"slice timing correction was applied to the timeseries (3dTshift)[{references[afni]}], then "
    else:
        stc_str=""

    # Resampling
    resampling_str=f"Finally, after calculating the transformations required to correct for head motion and susceptibility distortions, {stc_str}transforms were \
concatenated into a single resampling operation (avoiding multiple resampling) which is applied at each EPI frame, "

    if not opts.nativespace_resampling == 'inputs_defined':
        nativespace_resampling = f", resampled at a voxel resolution of {opts.nativespace_resampling}mm"
    else:
        nativespace_resampling = ""
    if not opts.commonspace_resampling == 'inputs_defined':
        commonspace_resampling = f", at a voxel resolution of {opts.commonspace_resampling}mm"
    else:
        commonspace_resampling = ""

    if not opts.bold_only: 
        if opts.coreg_script=='SyN':
            # Cross-modal alignment
            references[sdc_nonlinear]=i
            i+=1
            sdc_str=f"To correct for EPI susceptibility distortions, the volumetric EPI is also subjected to inhomogeneity correction, and then registered \
using a {define_registration(opts.coreg_script)} registration to the anatomical scan from the same MRI session, which allows to calculate the required \
geometrical transforms for recovering brain anatomy [{references[sdc_nonlinear]}]. "
        else:
            sdc_str="[NO SUSCEPTIBILITY DISTORTION CORRECTION CONDUCTED]"

        references[fmriprep]=i
        i+=1
        resampling_str+=f"generating the preprocessed EPI timeseries in native space [{references[fmriprep]}]{nativespace_resampling}. \
Preprocessed timeseries in common space are also generated by further concatenating the transforms allowing resampling to the reference atlas{commonspace_resampling}. "

    if opts.bold_only:
        if opts.atlas_reg_script=='SyN':
            references[sdc_nonlinear]=i
            i+=1
            sdc_str=f"The nonlinear registration to common space allows to calculate the required geometrical transforms for recovering brain anatomy, and \
is thus used to correct susceptibility distortions of the EPI[{references[sdc_nonlinear]}]. "
        else:
            sdc_str="[NO SUSCEPTIBILITY DISTORTION CORRECTION CONDUCTED]"

        references[fmriprep]=i
        i+=1
        resampling_str+=f"generating the preprocessed EPI timeseries [{references[fmriprep]}]. Since susceptibility distortion correction is completed with \
the registration to common space, the preprocessed timeseries are resampled to common space{commonspace_resampling}. "

    methods+=sdc_str
    methods+=resampling_str

    ref_string=''
    for key in references.keys():
        ref_string+=f"[{references[key]}] {key} \n"
    return methods,ref_string


def confound_correction_boilerplate(opts):
    methods=''
    references={}
    i=1

    # define references
    rabies_ref='Desrosiers-Gregoire et al., in preparation.'
    aroma='Pruim, R. H., Mennes, M., van Rooij, D., Llera, A., Buitelaar, J. K., & Beckmann, C. F. (2015). ICA-AROMA: A robust ICA-based strategy for removing motion artifacts from fMRI data. Neuroimage, 112, 267-277.'
    Lindquist='Lindquist, M. A., Geuter, S., Wager, T. D., & Caffo, B. S. (2019). Modular preprocessing pipelines can reintroduce artifacts into fMRI data. Human Brain Mapping, 40(8), 2358–2376.'
    mathias='Mathias, A., Grond, F., Guardans, R., Seese, D., Canela, M., & Diebner, H. H. (2004). Algorithms for spectral analysis of irregularly sampled time series. Journal of Statistical Software, 11(1), 1–27.'
    melodic='Beckmann, C. F., & Smith, S. M. (2004). Probabilistic independent component analysis for functional magnetic resonance imaging. IEEE transactions on medical imaging, 23(2), 137-152.'
    nilearn='Abraham, A., Pedregosa, F., Eickenberg, M., Gervais, P., Mueller, A., Kossaifi, J., ... & Varoquaux, G. (2014). Machine learning for neuroimaging with scikit-learn. Frontiers in neuroinformatics, 8, 14.'
    power2012='Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage, 59(3), 2142-2154.'
    power2014='Power, J. D., Mitra, A., Laumann, T. O., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2014). Methods to detect, characterize, and remove motion artifact in resting state fMRI. NeuroImage, 84, 320–341.'
    friston24='Friston, K. J., Williams, S., Howard, R., Frackowiak, R. S., & Turner, R. (1996). Movement‐related effects in fMRI time‐series. Magnetic resonance in medicine, 35(3), 346-355.'
    aCompCor='Muschelli, J., Nebel, M. B., Caffo, B. S., Barber, A. D., Pekar, J. J., & Mostofsky, S. H. (2014). Reduction of motion-related artifacts in resting state fMRI using aCompCor. Neuroimage, 96, 22-35.'

    # Commonspace VS native space
    if opts.nativespace_analysis:
        EPI_space=f"native space EPI timeseries"
    else:
        EPI_space=f"EPI timeseries resampled to commonspace"

    # citing RABIES
    references[rabies_ref]=i
    i+=1
    methods+=f"  Confound correction was executed using the RABIES software (https://github.com/CoBrALab/RABIES)[{references[rabies_ref]}] on {EPI_space}. "

    # Frame censoring
    if opts.FD_censoring or opts.DVARS_censoring:
        references[power2012]=i
        i+=1
        methods+=f"First, frames with prominent corruption were censored (i.e. scrubbing [{references[power2012]}]). "

        if opts.FD_censoring:
            methods+=f"Framewise displacement [{references[power2012]}] was measured across time and each frame surpassing {opts.FD_threshold}mm of motion, together with 1 backward and 2 forward frames, were removed. "

        if opts.DVARS_censoring:
            if not power2012 in references.keys():
                references[power2012]=i
                i+=1
            methods+=f"Using the DVARS measure[{references[power2012]}] of temporal fluctuations in global signal, timepoints presenting outlier DVARS values, characteristic of confounds, were removed. \
This was conducted by iteratively removing frames which present outlier DVARS values above or below 2.5 standard deviations until no more outliers are detected. "
        methods+="Next, "
    else:
        methods+="First, "

    # Detrending
    methods+=f"voxelwise linear detrending was first applied to remove first-order drifts and the average image. "

    # ICA-AROMA
    if opts.run_aroma:
        references[aroma]=i
        i+=1
        methods+=f"Then, motion sources were then automatically removed using a modified version of the ICA-AROMA classifier[{references[aroma]}], \
where classifier parameters and anatomical masks are instead adapted for rodent images. "
        if not opts.aroma_dim==0:
            references[melodic]=i
            i+=1
            methods+=f"Specifically {opts.aroma_dim} components were derived for each image independently using MELODIC-ICA algorithm[{melodic}] \
before classification. "

    # High/lowpass filtering
    highpass=opts.highpass is not None
    lowpass=opts.lowpass is not None
    if highpass or lowpass:
        references[nilearn]=i
        i+=1
        if highpass and lowpass:
            filter='bandpass'
            freq=f"{opts.highpass}-{opts.lowpass}"
        elif highpass:
            filter='highpass'
            freq=opts.highpass
        elif lowpass:
            filter='lowpass'
            freq=opts.lowpass
        if opts.FD_censoring or opts.DVARS_censoring:
            references[power2014]=i
            i+=1
            references[mathias]=i
            i+=1
            methods+=f"Next, frequency filtering requires particular considerations when frame censoring was \
previously applied, which results in missing timepoints. Conventional filters cannot account for missing data, \
and to address this issue, Power and colleagues introduced a method for simulating missing timepoints [{references[power2014]}]. \
This method relies on an adaptation of the Lomb-Scargle periodogram, which allows estimating the frequency composition of \
timeseries with missing data points, and from that estimation, missing timepoints can be simulated while preserving the \
frequency profile [{references[mathias]}]. Missing data is thereby simulated with preserved frequency properties for both \
BOLD timeseries and nuisance regressors prior to frequency filtering. "

            methods+=f"Following the simulation, "
        else:
            methods+=f"Next, "

        if not power2014 in references.keys():
            references[power2014]=i
            i+=1
        methods+=f"{filter} filtering({freq}Hz) was applied using a 3rd-order Butterworth filter, and {opts.edge_cutoff} seconds \
is removed at each edge of the timeseries to account for edge artefacts following filtering [{references[power2014]}]. The removal \
of 30 seconds and the use of a 3rd order filter was selected based on the visualization of edge artefacts with simulated data. "

        if len(opts.conf_list) > 0:
            references[Lindquist]=i
            i+=1
            methods+=f"The nuisance regressors are also filtered to ensure orthogonality between the frequency filters and subsequent confound regression, which can otherwise re-introduce removed confounds [{references[Lindquist]}]."

        if opts.FD_censoring or opts.DVARS_censoring:
            methods+=f"After frequency filtering, the temporal masks are re-applied to remove the simulated timepoints. "

    # Confound Regression
    if len(opts.conf_list)>0:
        str_list=[]
        if 'mot_6' in opts.conf_list:
            str_list.append('the 6 rigid motion parameters')
        if 'mot_24' in opts.conf_list:
            references[friston24]=i
            i+=1
            str_list.append(f'24 motion parameters (the 6 rigid motion parameters, their temporal derivative, together with all 12 parameters squared[{references[friston24]}])')
        if 'aCompCor' in opts.conf_list:
            references[aCompCor]=i
            i+=1
            str_list.append(f'the timecourses from principal components explaining 50% of the variance among the voxels from the combined WM and CSF masks (aCompCor)[{references[aCompCor]}]')
        if 'mean_FD' in opts.conf_list:
            str_list.append(f'the mean framewise displacement')
        if 'WM_signal' in opts.conf_list or 'CSF_signal' in opts.conf_list or 'vascular_signal' in opts.conf_list:
            mask_str='the mean signal from the '

            masks=[]
            if 'WM_signal' in opts.conf_list:
                masks.append('WM')
            if 'CSF_signal' in opts.conf_list:
                masks.append('CSF')
            if 'vascular_signal' in opts.conf_list:
                masks.append('vascular')

            mask_str+=masks[0]
            for str_ in masks[1:-1]:
                mask_str+=', '+str_
            if len(masks)>1:
                mask_str+=f' and {masks[-1]} masks'
            else:
                mask_str+=' mask'

            str_list.append(mask_str)

        if 'global_signal' in opts.conf_list:
            str_list.append(f'the global signal')

        conf_list_str=str_list[0]
        for str_ in str_list[1:-1]:
            conf_list_str+=', '+str_
        if len(str_list)>1:
            conf_list_str+=f' and {str_list[-1]}'


        methods+=f"Selected nuisance regressors were then used for confound regression. More specifically, using ordinary least square regression, \
{conf_list_str} were modelled at each voxel and regressed from the data. "


    if opts.standardize or opts.smoothing_filter is not None:
        methods+="Before analysis, "
    # Standardize
    if opts.standardize:
        methods+=f"timeseries were temporally standardized (mean-substrated, variance-normalized) and "

    # Spatial smoothing
    if opts.smoothing_filter is not None:
        if not nilearn in references.keys():
            references[nilearn]=i
            i+=1
        methods+=f"a spatial Gaussian smoothing filter(nilearn.image.smooth_img)[{references[nilearn]}] was applied at {opts.smoothing_filter}mm full-width at half maximum (FWHM). "

    methods+='\n'
    ref_string=''
    for key in references.keys():
        ref_string+=f"[{references[key]}] {key} \n"
    return methods,ref_string
