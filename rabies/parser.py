import os
import argparse
from pathlib import Path
import pathos.multiprocessing as multiprocessing  # Better multiprocessing

if 'XDG_DATA_HOME' in os.environ.keys():
    rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
else:
    rabies_path = os.environ['HOME']+'/.local/share/rabies'


def get_parser():
    """Build parser object"""
    parser = argparse.ArgumentParser(
        description=
            "RABIES performs multiple stages of rodent fMRI image processing, including preprocessing, \n"
            "confound correction, simple analyses and data quality assessment.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparsers = parser.add_subparsers(
        title='Processing options',
        description=
            "The RABIES workflow is seperated into three main processing stages: preprocessing, \n"
            "confound correction and analysis. Outputs from the preprocessing provide the inputs for\n"
            "the subsequent confound correction, and finally analysis.",
        help='Description',
        dest='rabies_stage',
        metavar='Processing stage')

    preprocess = subparsers.add_parser("preprocess",
        help=
            "\n"
            "Conducts preprocessing on an input dataset in BIDS format. Preprocessing includes \n"
            "motion realignment, susceptibility distortions correction through non-linear \n"
            "registration, alignment to commonspace, anatomical parcellation and evaluation of \n"
            "nuisance timecourses.\n"
            "\n",
        formatter_class=argparse.RawTextHelpFormatter)
    confound_correction = subparsers.add_parser("confound_correction",
        help=
            "\n"
            "Flexible options for confound correction are applied directly on preprocessing outputs\n"
            "from RABIES to derive cleaned timeseries. Various correction strategies, if selected, are\n"
            "applied in the following order, following best practices from human litterature:\n"
            "   #1 - Compute and apply frame censoring mask (from FD and/or DVARS thresholds)\n"
            "   #2 - If --match_number_timepoints is selected, each scan is matched to the \n"
            "       defined minimum_timepoint number of frames.\n"
            "   #3 - Linear/Quadratic detrending of fMRI timeseries and nuisance regressors\n"
            "   #4 - Apply ICA-AROMA.\n"
            "   #5 - If frequency filtering and frame censoring are applied, simulate data in censored\n" 
            "       timepoints using the Lomb-Scargle periodogram, as suggested in Power et al. (2014, \n"
            "       Neuroimage), for both the fMRI timeseries and nuisance regressors prior to filtering.\n"
            "   #6 - As recommended in Lindquist et al. (2019, Human brain mapping), make the nuisance \n"
            "       regressors orthogonal to the temporal frequency filter.\n"
            "   #7 - Apply highpass and/or lowpass filtering on the fMRI timeseries (with simulated \n"
            "       timepoints).\n"
            "   #8 - Re-apply the frame censoring mask onto filtered fMRI timeseries and nuisance \n"
            "       regressors, taking out the simulated timepoints. Edge artefacts from frequency \n"
            "       filtering can also be removed as recommended in Power et al. (2014, Neuroimage).\n"
            "   #9 - Apply confound regression using the selected nuisance regressors (see --conf_list\n" 
            "       options).\n"
            "   #10 - Scaling of timeseries variance\n"
            "   #11 - Apply Gaussian spatial smoothing.\n"
            "\n",
        formatter_class=argparse.RawTextHelpFormatter)
    analysis = subparsers.add_parser("analysis",
        help=
            "\n"
            "Conduct simple resting-state functional connectivity (FC) analysis, or data quality\n" 
            "diagnosis, on cleaned timeseries after confound correction. Analysis options include\n"
            "seed-based FC, whole-brain FC matrix, group-ICA and dual regression. --data_diagnosis\n" 
            "computes features of data quality at the individual scan and group levels, as in \n"
            "Desrosiers-Gregoire et al. (in prep)\n" 
            "\n",
        formatter_class=argparse.RawTextHelpFormatter)

    ####Execution
    g_execution = parser.add_argument_group(
        title='Execution Options', 
        description=
            "Options for parallel execution and memory management."
        )
    g_execution.add_argument(
        "-p", "--plugin", default='Linear',
        choices=['Linear', 'MultiProc', 'SGE', 'SGEGraph',
                'PBS', 'LSF', 'SLURM', 'SLURMGraph'],
        help=
            "Specify the nipype plugin for workflow execution.\n"
            "Consult https://nipype.readthedocs.io/en/0.11.0/users/plugins.html for details.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_execution.add_argument(
        '--local_threads', type=int, default=multiprocessing.cpu_count(),
        help=
            "For --plugin MultiProc, set the maximum number of processors run in parallel.\n"
            "Defaults to number of CPUs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_execution.add_argument(
        "--scale_min_memory", type=float, default=1.0,
        help=
            "For --plugin MultiProc, set the memory scaling factor attributed to nodes during\n"
            "execution. Increase the scaling if memory crashes are reported.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_execution.add_argument(
        "--min_proc", type=int, default=1,
        help=
            "For --plugin SGE/SGEGraph, scale the number of nodes attributed to jobs to\n"
            "avoid memory crashes.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_execution.add_argument(
        "--verbose", type=int, default=1,
        help=
            "Set the verbose level. 0=WARNING, 1=INFO, 2 or above=DEBUG.\n"
            "(default: %(default)s)\n"
            "\n"
        )


    ####Preprocessing
    preprocess.add_argument(
        'bids_dir', action='store', type=Path,
        help=
            "The root folder of the BIDS-formated input data directory.\n"
            "\n"
        )
    preprocess.add_argument(
        'output_dir', action='store', type=Path,
        help=
            "the output path to drop outputs from major preprocessing steps.\n"
            "\n"
        )
    preprocess.add_argument(
        "--bold_only", dest='bold_only', action='store_true',
        help=
            "Apply preprocessing with only EPI scans. Commonspace registration is executed directly using\n"
            "the corrected EPI 3D reference images. The commonspace registration simultaneously applies\n"
            "distortion correction, this option will produce only commonspace outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--anat_autobox', dest='anat_autobox', action='store_true',
        help=
            "Crops out extra space around the brain on the structural image using AFNI's 3dAutobox\n"
            "https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutobox.html.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--bold_autobox', dest='bold_autobox', action='store_true',
        help=
            "Crops out extra space around the brain on the EPI image using AFNI's 3dAutobox\n"
            "https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutobox.html.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--apply_despiking', dest='apply_despiking', action='store_true',
        help=
            "Applies AFNI's 3dDespike https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        "--HMC_option", type=str, default='intraSubjectBOLD',
        choices=['intraSubjectBOLD', '0', '1', '2', '3'],
        help=
            "Select an option for head motion realignment among the pre-built options from\n"
            "https://github.com/ANTsX/ANTsR/blob/master/R/ants_motion_estimation.R.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--apply_slice_mc', dest='apply_slice_mc', action='store_true',
        help=
            "Whether to apply a slice-specific motion correction after initial volumetric HMC. This can \n"
            "correct for interslice misalignment resulting from within-TR motion. With this option, \n"
            "motion corrections and the subsequent resampling from registration are applied sequentially\n"
            "since the 2D slice registrations cannot be concatenate with 3D transforms. \n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        '--detect_dummy', dest='detect_dummy', action='store_true',
        help=
            "Detect and remove initial dummy volumes from the EPI, and generate a reference EPI based on\n"
            "these volumes if detected. Dummy volumes will be removed from the output preprocessed EPI.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    preprocess.add_argument(
        "--data_type", type=str, default='float32',
        choices=['int16', 'int32', 'float32', 'float64'],
        help=
            "Specify data format outputs to control for file size.\n"
            "(default: %(default)s)\n"
            "\n"
        )

    g_registration = preprocess.add_argument_group(
        title='Registration Options', 
        description=
            "Customize registration operations and troubleshoot registration failures.\n"
        )
    g_registration.add_argument(
        "--anat_inho_cor", type=str, 
        default='method=SyN,otsu_thresh=2,multiotsu=false',
        help=
            "Select options for the inhomogeneity correction of the structural image.\n"
            "* method: specify which registration strategy is employed for providing a brain mask.\n"
            "*** Rigid: conducts only rigid registration.\n"
            "*** Affine: conducts Rigid then Affine registration.\n"
            "*** SyN: conducts Rigid, Affine then non-linear registration.\n"
            "*** no_reg: skip registration.\n"
            "*** N4_reg: previous correction script prior to version 0.3.1.\n"
            "*** disable: disables the inhomogeneity correction.\n"
            "* otsu_thresh: The inhomogeneity correction script necessitates an initial correction with a \n"
            " Otsu masking strategy (prior to registration of an anatomical mask). This option sets the \n"
            " Otsu threshold level to capture the right intensity distribution. \n"
            "*** Specify an integer among [0,1,2,3,4]. \n"
            "* multiotsu: Select this option to perform a staged inhomogeneity correction, where only \n"
            " lower intensities are initially corrected, then higher intensities are iteratively \n"
            " included to eventually correct the whole image. This technique may help with images with \n"
            " particularly strong inhomogeneity gradients and very low intensities.\n"
            "*** Specify 'true' or 'false'. \n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        '--anat_robust_inho_cor', type=str,
        default='apply=false,masking=false,brain_extraction=false,template_registration=SyN',
        help=
            "When selecting this option, inhomogeneity correction is executed twice to optimize \n"
            "outcomes. After completing an initial inhomogeneity correction step, the corrected outputs \n"
            "are co-registered to generate an unbiased template, using the same method as the commonspace \n"
            "registration. This template is then masked, and is used as a new target for masking during a \n"
            "second iteration of inhomogeneity correction. Using this dataset-specific template should \n"
            "improve the robustness of masking for inhomogeneity correction.\n"
            "* apply: select 'true' to apply this option. \n"
            " *** Specify 'true' or 'false'. \n"
            "* masking: Combine masks derived from the inhomogeneity correction step to support \n"
            " registration during the generation of the unbiased template, and then during template \n"            
            " registration.\n"
            " *** Specify 'true' or 'false'. \n"
            "* brain_extraction: conducts brain extraction prior to template registration based on the \n"
            " combined masks from inhomogeneity correction. This will enhance brain edge-matching, but \n"
            " requires good quality masks. This should be selected along the 'masking' option.\n"
            " *** Specify 'true' or 'false'. \n"
            "* template_registration: Specify a registration script for the alignment of the \n"
            " dataset-generated unbiased template to a reference template for masking.\n"
            "*** Rigid: conducts only rigid registration.\n"
            "*** Affine: conducts Rigid then Affine registration.\n"
            "*** SyN: conducts Rigid, Affine then non-linear registration.\n"
            "*** no_reg: skip registration.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        "--bold_inho_cor", type=str, 
        default='method=Rigid,otsu_thresh=2,multiotsu=false',
        help=
            "Same as --anat_inho_cor, but for the EPI images.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        '--bold_robust_inho_cor', type=str,
        default='apply=false,masking=false,brain_extraction=false,template_registration=SyN',
        help=
            "Same as --anat_robust_inho_cor, but for the EPI images.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        '--commonspace_reg', type=str,
        default='masking=false,brain_extraction=false,template_registration=SyN,fast_commonspace=false',
        help=
            "Specify registration options for the commonspace registration.\n"
            "* masking: Combine masks derived from the inhomogeneity correction step to support \n"
            " registration during the generation of the unbiased template, and then during template \n"            
            " registration.\n"
            "*** Specify 'true' or 'false'. \n"
            "* brain_extraction: conducts brain extraction prior to template registration based on the \n"
            " combined masks from inhomogeneity correction. This will enhance brain edge-matching, but \n"
            " requires good quality masks. This should be selected along the 'masking' option.\n"
            "*** Specify 'true' or 'false'. \n"
            "* template_registration: Specify a registration script for the alignment of the \n"
            " dataset-generated unbiased template to the commonspace atlas.\n"
            "*** Rigid: conducts only rigid registration.\n"
            "*** Affine: conducts Rigid then Affine registration.\n"
            "*** SyN: conducts Rigid, Affine then non-linear registration.\n"
            "*** no_reg: skip registration.\n"
            "* fast_commonspace: Skip the generation of a dataset-generated unbiased template, and \n"
            " instead, register each scan independently directly onto the commonspace atlas, using the \n"
            " template_registration. This option can be faster, but may decrease the quality of \n"
            " alignment between subjects. \n"
            "*** Specify 'true' or 'false'. \n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_registration.add_argument(
        "--bold2anat_coreg", type=str, 
        default='masking=false,brain_extraction=false,registration=SyN',
        help=
            "Specify the registration script for cross-modal alignment between the EPI and structural\n"
            "images. This operation is responsible for correcting EPI susceptibility distortions.\n"
            "* masking: With this option, the brain masks obtained from the EPI inhomogeneity correction \n"
            " step are used to support registration.\n"
            "*** Specify 'true' or 'false'. \n"
            "* brain_extraction: conducts brain extraction prior to registration using the EPI masks from \n"
            " inhomogeneity correction. This will enhance brain edge-matching, but requires good quality \n"
            " masks. This should be selected along the 'masking' option.\n"
            "*** Specify 'true' or 'false'. \n"
            "* registration: Specify a registration script.\n"
            "*** Rigid: conducts only rigid registration.\n"
            "*** Affine: conducts Rigid then Affine registration.\n"
            "*** SyN: conducts Rigid, Affine then non-linear registration.\n"
            "*** no_reg: skip registration.\n"
            "(default: %(default)s)\n"
            "\n"

        )

    g_resampling = preprocess.add_argument_group(
        title='Resampling Options', 
        description=
            "The following options allow to resample the voxel dimensions for the preprocessed EPIs\n"
            "or for the anatomical images during registration.\n"
            "The resampling syntax must be 'dim1xdim2xdim3' (in mm), follwing the RAS axis convention\n"
            "(dim1=Right-Left, dim2=Anterior-Posterior, dim3=Superior-Inferior). If 'inputs_defined'\n"
            "is provided instead of axis dimensions, the original dimensions are preserved.\n"
        )
    g_resampling.add_argument(
        '--nativespace_resampling', type=str, default='inputs_defined',
        help=
            "Can specify a resampling dimension for the nativespace fMRI outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_resampling.add_argument(
        '--commonspace_resampling', type=str, default='inputs_defined',
        help=
            "Can specify a resampling dimension for the commonspace fMRI outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_resampling.add_argument(
        '--anatomical_resampling', type=str, default='inputs_defined',
        help=
            "\n"
            "This specifies resampling dimensions for the anatomical registration targets. By \n"
            "default, images are resampled to isotropic resolution based on the smallest dimension\n"
            "among the provided anatomical images (EPI images instead if --bold_only is True). \n"
            "Increasing voxel resampling size will increase registration speed at the cost of \n"
            "accuracy.\n"
            "(default: %(default)s)\n"
            "\n"
        )

    g_stc = preprocess.add_argument_group(
        title='STC Options', 
        description=
            "Specify Slice Timing Correction (STC) info that is fed to AFNI's 3dTshift\n"
            "(https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html). The STC is applied\n"
            "in the anterior-posterior orientation, and thus RABIES assumes slices were acquired in\n"
            "this direction.\n"
        )
    g_stc.add_argument(
        '--apply_STC', dest='apply_STC', action='store_true',
        help=
            "Select this option to apply the STC step.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_stc.add_argument(
        '--TR', type=str, default='auto',
        help=
            "Specify repetition time (TR) in seconds. (e.g. --TR 1.2)\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_stc.add_argument(
        '--tpattern', type=str, default='alt-z',
        choices=['alt-z', 'seq-z', 'alt+z', 'seq+z'],
        help=
            "Specify if interleaved ('alt') or sequential ('seq') acquisition, and specify in which \n"
            "direction (- or +) to apply the correction. If slices were acquired from front to back, \n"
            "the correction should be in the negative (-) direction. Refer to this discussion on the \n"
            "topic for more information https://github.com/CoBrALab/RABIES/discussions/217.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_stc.add_argument(
        '--stc_axis', type=str, default='Y',
        choices=['X', 'Y', 'Z'],
        help=
            "Can specify over which axis of the image the STC must be applied. Generally, the correction \n"
            "should be over the Y axis, which corresponds to the anteroposterior axis in RAS convention. \n"
            "(default: %(default)s)\n"
            "\n"
        )

    g_atlas = preprocess.add_argument_group(
        title='Template Files', 
        description=
            "Specify commonspace template and associated mask/label files. By default, RABIES\n"
            "provides the mouse DSURQE atlas\n"
            "https://wiki.mouseimaging.ca/display/MICePub/Mouse+Brain+Atlases.\n"
        )
    g_atlas.add_argument(
        '--anat_template', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_average.nii.gz",
        help=
            "Anatomical file for the commonspace atlas.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--brain_mask', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_mask.nii.gz",
        help=
            "Brain mask aligned with the template.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--WM_mask', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_eroded_WM_mask.nii.gz",
        help=
            "White matter mask aligned with the template.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--CSF_mask', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_eroded_CSF_mask.nii.gz",
        help=
            "CSF mask aligned with the template.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--vascular_mask', action='store', type=Path,
        default=f"{rabies_path}/vascular_mask.nii.gz",
        help=
            "Can provide a mask of major blood vessels to compute associated nuisance timeseries.\n"
            "The default mask was generated by applying MELODIC ICA and selecting the resulting \n"
            "component mapping onto major brain vessels.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    g_atlas.add_argument(
        '--labels', action='store', type=Path,
        default=f"{rabies_path}/DSURQE_40micron_labels.nii.gz",
        help=
            "Labels file providing the atlas anatomical annotations.\n"
            "(default: %(default)s)\n"
            "\n"
        )



    ####Confound correction
    confound_correction.add_argument(
        'preprocess_out', action='store', type=Path,
        help=
            "path to RABIES preprocessing output directory.\n"
            "\n"
        )
    confound_correction.add_argument(
        'output_dir', action='store', type=Path,
        help=
            "path for confound correction output directory.\n"
            "\n"
        )
    confound_correction.add_argument(
        '--nativespace_analysis', dest='nativespace_analysis', action='store_true',
        help=
            "Conduct confound correction and analysis in native space.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--image_scaling', type=str,
        default="None",
        choices=["None", "background_noise", "global_variance", "voxelwise_standardization"],
        help=
            "Select an option for scaling the image variance to match the intensity profile of \n"
            "different scans and avoid biases in data variance and amplitude estimation during analysis.\n"
            "The variance explained from confound regression is also scaled accordingly for later use with \n"
            "--data_diagnosis. \n"
            "*** None: No scaling is applied, only detrending.\n"
            "*** background_noise: a mask is derived to map background noise, and scale the image \n"
            "   intensity relative to the noise standard deviation. \n"
            "*** global_variance: After applying confound correction, the cleaned timeseries are scaled \n"
            "   according to the total standard deviation of all voxels, to scale total variance to 1. \n"
            "*** voxelwise_standardization: After applying confound correction, each voxel is separately \n"
            "   scaled to unit variance (z-scoring). \n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--detrending_order', type=str,
        default="linear",
        choices=["linear", "quadratic"],
        help=
            "Select between linear or quadratic (second-order) detrending of voxel timeseries.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--conf_list', type=str,
        nargs="*",  # 0 or more values expected => creates a list
        default=[],
        choices=["WM_signal", "CSF_signal", "vascular_signal",
                "global_signal", "aCompCor", "mot_6", "mot_24", "mean_FD"],
        help=
            "Select list of nuisance regressors that will be applied on voxel timeseries, i.e., confound\n"
            "regression.\n"
            "*** WM/CSF/vascular/global_signal: correspond to mean signal from WM/CSF/vascular/brain \n"
            "   masks.\n"
            "*** mot_6: 6 rigid head motion correction parameters.\n"
            "*** mot_24: mot_6 + their temporal derivative, then all 12 parameters squared, as in \n"
            "   Friston et al. (1996, Magnetic Resonance in Medicine).\n"
            "*** aCompCor: method from Muschelli et al. (2014, Neuroimage), where component timeseries\n"
            "   are obtained using PCA, conducted on the combined WM and CSF masks voxel timeseries. \n"
            "   Components adding up to 50 percent of the variance are included.\n"
            "*** mean_FD: the mean framewise displacement timecourse.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--frame_censoring', type=str, default='FD_censoring=false,FD_threshold=0.05,DVARS_censoring=false,minimum_timepoint=3',
        help=
            "Censor frames that are highly corrupted (i.e. 'scrubbing'). \n"
            "* FD_censoring: Apply frame censoring based on a framewise displacement threshold. The frames \n"
            " that exceed the given threshold, together with 1 back and 2 forward frames will be masked \n"
            " out, as in Power et al. (2012, Neuroimage).\n"
            "*** Specify 'true' or 'false'. \n"
            "* FD_threshold: the FD threshold in mm. \n"
            "* DVARS_censoring: Will remove timepoints that present outlier values on the DVARS metric \n"
            " (temporal derivative of global signal). This method will censor timepoints until the \n" 
            " distribution of DVARS values across time does not contain outliers values above or below 2.5 \n" 
            " standard deviations.\n"
            "*** Specify 'true' or 'false'. \n"
            "* minimum_timepoint: Can set a minimum number of timepoints remaining after frame censoring. \n" 
            " If the threshold is not met, an empty file is generated and the scan is not considered in \n" 
            " further steps. \n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--TR', type=str, default='auto',
        help=
            "Specify repetition time (TR) in seconds. (e.g. --TR 1.2)\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--highpass', type=float, default=None,
        help=
            "Specify highpass filter frequency.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--lowpass', type=float, default=None,
        help=
            "Specify lowpass filter frequency.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--edge_cutoff', type=float, default=0,
        help=
            "Specify the number of seconds to cut at beginning and end of acquisition if applying a\n"
            "frequency filter. Frequency filters generate edge effects at begining and end of the\n" 
            "timeseries. We recommend to cut those timepoints (around 30sec at both end for 0.01Hz \n" 
            "highpass.).\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--smoothing_filter', type=float, default=None,
        help=
            "Specify filter size in mm for spatial smoothing. Will apply nilearn's function \n"
            "https://nilearn.github.io/modules/generated/nilearn.image.smooth_img.html\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--match_number_timepoints', dest='match_number_timepoints', action='store_true', default=False,
        help=
            "With this option, only a subset of the timepoints are kept post-censoring to match the \n"
            "--minimum_timepoint number for all scans. This can be conducted to avoid inconsistent \n" 
            "temporal degrees of freedom (tDOF) between scans during downstream analysis. We recommend \n" 
            "selecting this option if a significant confounding effect of tDOF is detected during --data_diagnosis.\n" 
            "The extra timepoints removed are randomly selected among the set available post-censoring.\n" 
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--ica_aroma', type=str, default='apply=false,dim=0,random_seed=1',
        help=
            "Apply ICA-AROMA denoising (Pruim et al. 2015). The original classifier was modified to incorporate \n"
            "rodent-adapted masks and classification hyperparameters.\n"
            "* apply: apply the denoising.\n"
            "*** Specify 'true' or 'false'. \n"
            "* dim: Specify a pre-determined number of MELODIC components to derive. '0' will use an automatic \n"
            " estimator. \n"
            "* random_seed: For reproducibility, this option sets a fixed random seed for MELODIC. \n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--read_datasink', dest='read_datasink', action='store_true', default=False,
        help=
            "\n"
            "Choose this option to read preprocessing outputs from datasinks instead of the saved \n"
            "preprocessing workflow graph. This allows to run confound correction without having \n"
            "available RABIES preprocessing folders, but the targetted datasink folders must follow the\n"
            "structure of RABIES preprocessing.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    confound_correction.add_argument(
        '--timeseries_interval', type=str, default='all',
        help=
            "Before confound correction, can crop the timeseries within a specific interval.\n"
            "e.g. '0,80' for timepoint 0 to 80.\n"
            "(default: %(default)s)\n"
            "\n"
        )


    ####Analysis
    analysis.add_argument(
        'confound_correction_out', action='store', type=Path,
        help=
            "path to RABIES confound correction output directory.\n"
            "\n"
        )
    analysis.add_argument(
        'output_dir', action='store', type=Path,
        help=
            "path for analysis outputs.\n"
            "\n"
        )
    analysis.add_argument(
        '--scan_list', type=str,
        nargs="*",  # 0 or more values expected => creates a list
        default=['all'],
        help=
            "This option offers to run the analysis on a subset of the scans. The scans are selected by\n"
            "providing the full path to the corresponding EPI file in the input BIDS folder. The list \n"
            "of scan can be specified manually as a list of file name '--scan_list scan1.nii.gz \n"
            "scan2.nii.gz ...' or the files can be imbedded into a .txt file with one filename per row.\n"
            "By default, 'all' will use all the scans previously processed.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--prior_maps', action='store', type=Path,
        default=f"{rabies_path}/melodic_IC.nii.gz",
        help=
            "Provide a 4D nifti image with a series of spatial priors representing common sources of\n"
            "signal (e.g. ICA components from a group-ICA run). This 4D prior map file will be used for \n"
            "Dual regression, Dual ICA and --data_diagnosis. The RABIES default corresponds to a MELODIC \n"
            "run on a combined group of anesthetized-ventilated and awake mice. Confound correction \n"
            "consisted of highpass at 0.01 Hz, FD censoring at 0.03mm, DVARS censoring, and \n"
            "mot_6,WM_signal,CSF_signal as regressors.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--prior_bold_idx', type=int,
        nargs="*",  # 0 or more values expected => creates a list
        default=[5, 12, 19],
        help=
            "Specify the indices for the priors corresponding to BOLD sources from --prior_maps. These will\n"
            "be fitted during Dual ICA and provide the BOLD components during --data_diagnosis.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--prior_confound_idx', type=int,
        nargs="*",  # 0 or more values expected => creates a list
        default=[0, 1, 2, 6, 7, 8, 9, 10, 11,
                    13, 14, 21, 22, 24, 26, 28, 29],
        help=
            "Specify the indices for the confound components from --prior_maps. This is pertinent for the\n" 
            "--data_diagnosis outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        "--data_diagnosis", dest='data_diagnosis', action='store_true',
        help=
            "This option carries out the spatiotemporal diagnosis as described in Desrosiers-Gregoire et al. \n"
            "The diagnosis generates key temporal and spatial features both at the scan level and the group\n"
            "level, allowing the identification of sources of confounds and data quality issues. We recommend \n"
            "using this data diagnosis workflow, more detailed in the publication, to improve the control for \n"
            "data quality issues and prevent the corruptions of analysis outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )

    analysis.add_argument(
        '--seed_list', type=str,
        nargs="*",  # 0 or more values expected => creates a list
        default=[],
        help=
            "Can provide a list of Nifti files providing a mask for an anatomical seed, which will be used\n"
            "to evaluate seed-based connectivity maps using on Pearson's r. Each seed must consist of \n"
            "a binary mask representing the ROI in commonspace.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--seed_prior_list', type=str,
        nargs="*",  # 0 or more values expected => creates a list
        default=[],
        help=
            "For analysis QC of seed-based FC during --data_diagnosis, prior network maps are required for \n"
            "each seed provided in --seed_list. Provide the list of prior files in matching order of the \n"
            "--seed_list arguments to match corresponding seed maps.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument("--FC_matrix", dest='FC_matrix', action='store_true',
        help=
            "Compute whole-brain connectivity matrices using Pearson's r between ROI timeseries.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        "--ROI_type", type=str, default='parcellated',
        choices=['parcellated', 'voxelwise'],
        help=
            "Define ROIs for --FC_matrix between 'parcellated' from the provided atlas during preprocessing,\n"
            "or 'voxelwise' to derive the correlations between every voxel."
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        "--ROI_csv", action='store', type=Path, 
        default=f"{rabies_path}/DSURQE_40micron_labels.nii.gz",
        help=
            "A CSV file with the ROI names matching the ROI index numbers in the atlas labels Nifti file. \n"
            "A copy of this file is provided along the FC matrix generated for each subject.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--group_ica', type=str, default='apply=false,dim=0,random_seed=1',
        help=
            "Perform group-ICA using FSL's MELODIC on the whole dataset's cleaned timeseries.\n"
            "Note that confound correction must have been conducted on commonspace outputs.\n"
            "* apply: compute group-ICA.\n"
            "*** Specify 'true' or 'false'. \n"
            "* dim: Specify a pre-determined number of MELODIC components to derive. '0' will use an automatic \n"
            " estimator. \n"
            "* random_seed: For reproducibility, this option sets a fixed random seed for MELODIC. \n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        "--DR_ICA", dest='DR_ICA', action='store_true',
        help=
            "Conduct dual regression on each subject timeseries, using the priors from --prior_maps. The\n"
            "linear coefficients from both the first and second regressions will be provided as outputs.\n"
            "Requires that confound correction was conducted on commonspace outputs.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--NPR_temporal_comp', type=int, default=-1,
        help=
            "Option for performing Neural Prior Recovery (NPR). Specify with this option how many extra \n"
            "subject-specific sources will be computed to account for non-prior confounds. This options \n"
            "specifies the number of temporal components to compute. After computing \n"
            "these sources, NPR will provide a fit for each prior in --prior_maps indexed by --prior_bold_idx.\n"
            "Specify at least 0 extra sources to run NPR.\n"
            "(default: %(default)s)\n"
            "\n"
        )
    analysis.add_argument(
        '--NPR_spatial_comp', type=int, default=-1,
        help=
            "Same as --NPR_temporal_comp, but specify how many spatial components to compute (which are \n"
            "additioned to the temporal components).\n"
            "(default: %(default)s)\n"
            "\n"
        )

    return parser


def read_parser(parser):
    opts = parser.parse_args()

    if opts.rabies_stage == 'preprocess':
        opts.anat_inho_cor = parse_argument(opt=opts.anat_inho_cor, 
            key_value_pairs = {'method':['Rigid', 'Affine', 'SyN', 'no_reg', 'N4_reg', 'disable'], 
                'otsu_thresh':['0','1','2','3','4'], 'multiotsu':['true', 'false']},
            name='anat_inho_cor')

        opts.bold_inho_cor = parse_argument(opt=opts.bold_inho_cor, 
            key_value_pairs = {'method':['Rigid', 'Affine', 'SyN', 'no_reg', 'N4_reg', 'disable'], 
                'otsu_thresh':['0','1','2','3','4'], 'multiotsu':['true', 'false']},
            name='bold_inho_cor')

        opts.commonspace_reg = parse_argument(opt=opts.commonspace_reg, 
            key_value_pairs = {'masking':['true', 'false'], 'brain_extraction':['true', 'false'], 
                'template_registration':['Rigid', 'Affine', 'SyN', 'no_reg'], 'fast_commonspace':['true', 'false']},
            name='commonspace_reg')

        opts.bold2anat_coreg = parse_argument(opt=opts.bold2anat_coreg, 
            key_value_pairs = {'masking':['true', 'false'], 'brain_extraction':['true', 'false'], 
                'registration':['Rigid', 'Affine', 'SyN', 'no_reg']},
            name='bold2anat_coreg')

        opts.anat_robust_inho_cor = parse_argument(opt=opts.anat_robust_inho_cor, 
            key_value_pairs = {'apply':['true', 'false'], 'masking':['true', 'false'], 'brain_extraction':['true', 'false'], 
                'template_registration':['Rigid', 'Affine', 'SyN', 'no_reg']},
            name='anat_robust_inho_cor')

        opts.bold_robust_inho_cor = parse_argument(opt=opts.bold_robust_inho_cor, 
            key_value_pairs = {'apply':['true', 'false'], 'masking':['true', 'false'], 'brain_extraction':['true', 'false'], 
                'template_registration':['Rigid', 'Affine', 'SyN', 'no_reg']},
            name='bold_robust_inho_cor')

    elif opts.rabies_stage == 'confound_correction':
        opts.frame_censoring = parse_argument(opt=opts.frame_censoring, 
            key_value_pairs = {'FD_censoring':['true', 'false'], 'FD_threshold':float, 'DVARS_censoring':['true', 'false'],
                'minimum_timepoint':int},
            name='frame_censoring')

        opts.ica_aroma = parse_argument(opt=opts.ica_aroma, 
            key_value_pairs = {'apply':['true', 'false'], 'dim':int, 'random_seed':int},
            name='ica_aroma')

    elif opts.rabies_stage == 'analysis':
        opts.group_ica = parse_argument(opt=opts.group_ica, 
            key_value_pairs = {'apply':['true', 'false'], 'dim':int, 'random_seed':int},
            name='group_ica')

    return opts

def parse_argument(opt, key_value_pairs, name):
    key_list = list(key_value_pairs.keys())
    l = opt.split(',')
    opt_dict = {}
    for e in l:
        if not '=' in e:
            raise ValueError(f"Provided option must follow the 'key=value' syntax, {e} was found instead.")
        s = e.split('=')
        if not len(s)==2:
            raise ValueError(f"Provided option must follow the 'key=value' syntax, {e} was found instead.")
        [key,value] = s
        if not key in key_list:
            raise ValueError(f"The provided key {key} is not part of the available options {key_list}.")
        if key_value_pairs[key] in [int,float]:
            value = key_value_pairs[key](value)
        else:
            if not value in key_value_pairs[key]:
                raise ValueError(f"The provided value {value} is not part of the available options {key_value_pairs[key]} for the key {key}.")
            if value=='true':
                value=True
            elif value=='false':
                value=False
        opt_dict[key]=value

    for key in key_list:
        if not key in list(opt_dict.keys()):
            raise ValueError(f"The key {key} is missing from the necessary attributes for --{name}.")
    return opt_dict
