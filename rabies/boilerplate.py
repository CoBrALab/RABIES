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
    denoising='Jose V. Manjon, Pierrick Coupe, Luis Marti-Bonmati, D. Louis Collins, Montserrat Robles "Adaptive non-local means denoising of MR images with spatially varying noise levels" Journal of Magnetic Resonance Imaging Volume 31, Issue 1, pages 192–203, January 2010 DOI: 10.1002/jmri.2200'
    ants='Avants, B. B., Tustison, N., & Song, G. (2009). Advanced normalization tools (ANTS). Insight j, 2(365), 1-35.'
    inho_cor='J.G. Sled, A.P. Zijdenbos and A.C. Evans, "A non-parametric method for automatic correction of intensity non-uniformity in MRI data", in "IEEE Transactions on Medical Imaging", vol. 17, n. 1, pp. 87-97, 1998'
    reg_script='https://github.com/CoBrALab/minc-toolkit-extras/blob/master/ants_generate_iterations.py'
    sdc_nonlinear='Wang, S. et al. Evaluation of field map and nonlinear registration methods for correction of susceptibility artifacts in diffusion MRI. Front.Neuroinform. 11, 17 (2017).'


    # citing RABIES
    references[rabies_ref]=i
    i+=1
    methods+=f"  The preprocessing of fMRI images was conducted using the open-source RABIES software (https://github.com/CoBrALab/RABIES)[{references[rabies_ref]}]. "

    # RAS re-orientation
    methods+="Before preprocessing, each image was re-oriented to match the RAS convention. "

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
        detect_dummy="If dummy scans are detected, the median of these volumes is selected as reference, given their higher anatomical contrast. Otherwise, t"
    else:
        detect_dummy="T"

    # Generate 3D EPI ref
    references[denoising]=i
    i+=1
    methods+=f"A representative 3D volume (EPI ref.) was generated from the EPI for downstream motion estimation and \
registration steps. {detect_dummy}wo iterations of motion realignment to a median of the volumes \
was conducted, and a trimmed mean ignoring 5% extreme values was taken as final reference. The final image was then corrected using \
non-local means denoising[{references[denoising]}]. "


    # STC
    if opts.apply_STC:
        if not afni in references.keys():
            references[afni]=i
            i+=1
        methods+=f"Slice timing correction was applied to the timeseries (3dTshift)[{references[afni]}], and r"
    else:
        methods+="R"

    # HMC
    references[ants]=i
    i+=1
    methods+=f"igid head motion realignment parameters were estimated through rigid registration to this EPI ref. (antsMotionCorr)[{references[ants]}]. "

    # Slice MC
    if opts.apply_slice_mc:
        methods+="SLICE-MC NOT IMPLEMENTED YET "


    # EPI inhomogeneity correction
    references[inho_cor]=i
    i+=1
    methods+=f"\n   Prior to registration, the EPI ref. image was then denoised with non-local mean denoising[{references[denoising]}] followed by iterative correction for intensity inhomogeneities [{references[inho_cor]}]. \
Initial masking was achieved via intensity thresholding giving an initial correction of the image, and {define_registration(opts.bold_inho_cor_method)} registration was then \
conducted to align a brain mask derived from the reference atlas for a final round of correction. "

    if not opts.bold_only:
        # Anat inhomogeneity correction
        methods+=f"The anatomical image was subjected to similar inhomogeneity correction and denoising, with a {define_registration(opts.anat_inho_cor_method)} registration. "
        reg_image='anatomical image'
    else:
        reg_image='EPI ref. image'

    if not opts.fast_commonspace:
        # Unbiased template generation
        methods+=f"Following inhomogeneity correction, a dataset-specific unbiased template was generated from the input {reg_image}s. \
To do so, through a set of iterations, images were registered to a consensus average generated from the overlap of the previous iteration. Registrations were \
increasingly stringent, executing 2 iterations each for a rigid, then affine and finally nonlinear template generation. The alignment to the \
finalized template provides the individual transforms to align each scan. "

        # Commonspace alignment
        methods+=f"The unbiased template was then registered to the provided reference atlas with a {define_registration(opts.atlas_reg_script)} registration, providing the \
transforms to resample brain masks and atlases to each scan. "
    else:
        # Fast commonspace alignment
        methods+=f"Each {reg_image} was then registered individually to the provided reference atlas with a {define_registration(opts.atlas_reg_script)} registration, providing the \
transforms to resample brain masks and atlases to each scan. "

    if not opts.bold_only:
        # Cross-modal alignment
        methods+=f"Finally, the EPI ref. and anatomical images were co-registered within subject using a {define_registration(opts.coreg_script)} registration. "

    # Scale-independent registration method
    references[reg_script]=i
    i+=1
    methods+=f"For every registration operation, RABIES uses an adaptive scale-independent registration framework built on top of ANTs (antsRegistration)[{references[reg_script]}]. \n"

    if not opts.bold_only:
        transform_resampling="cross-modal alignment"
        if opts.coreg_script=='SyN':
            references[sdc_nonlinear]=i
            i+=1
            transform_resampling+=f", correcting for EPI distortions from the nonlinear registration [{references[sdc_nonlinear]}],"
    else:
        references[sdc_nonlinear]=i
        i+=1
        transform_resampling=f"alignment to common reference space, correcting for EPI distortions from the nonlinear registration [{references[sdc_nonlinear]}],"


    # EPI resampling for HMC and SDC
    methods+=f"  Each volume from the EPI was finally resampled using the transforms from the {transform_resampling} and using \
the head motion realignment parameters, hence providing minimally preprocessed EPI timeseries. Transforms were concatenated before application to conduct \
resampling with a single operation to mitigate interpolation effects from repeated resampling. \n"

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
    nilearn='Abraham, A., Pedregosa, F., Eickenberg, M., Gervais, P., Mueller, A., Kossaifi, J., ... & Varoquaux, G. (2014). Machine learning for neuroimaging with scikit-learn. Frontiers in neuroinformatics, 8, 14.'
    power2012='Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage, 59(3), 2142-2154.'
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

    # Detrending
    methods+=f"Voxelwise linear detrending was first applied to remove first-order drifts and the average image. "

    # ICA-AROMA
    if opts.run_aroma:
        references[aroma]=i
        i+=1
        methods+=f"Motion sources were then automatically removed using a modified version of the ICA-AROMA classifier[{references[aroma]}], \
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
        methods+=f"Next, {filter} filtering({freq}Hz) was applied(nilearn.signal.clean)[{references[nilearn]}]. "

    # Frame censoring
    if opts.FD_censoring:
        references[power2012]=i
        i+=1
        methods+=f"Timepoints with prominent motion corruption were censored from the data (scrubbing [{references[power2012]}]). \
Framewise displacement was measured across time and each frame surpassing {opts.FD_threshold}mm of motion, together with 1 backward and 4 forward frames, were removed. "

    if opts.DVARS_censoring:
        if not power2012 in references.keys():
            references[power2012]=i
            i+=1
        methods+=f"Frame censoring was conducted using the DVARS measure[{references[power2012]}] of temporal fluctuations in global signal. This built-in strategy \
removes timepoints presenting outlier values in DVARS, which is characteristic of confounds that are not necessarily well accounted for by framewise displacement estimates[{references[rabies_ref]}]. \
This is conducted by iteratively removing frames which present outlier DVARS values above or below 2.5 standard deviations until no more outliers are detected. "

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


        methods+=f"Estimated nuisance timecourses during preprocessing were then used for confound regression. More specifically, using ordinary least square regression, \
{conf_list_str} were modelled at each voxel and regressed from the data. "


    if opts.standardize or opts.smoothing_filter is not None:
        methods+="Before analysis, "
    # Standardize
    if opts.standardize:
        methods+=f"timeseries were normalized to unit variance and "

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
