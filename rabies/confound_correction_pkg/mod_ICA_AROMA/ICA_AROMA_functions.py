#!/usr/bin/env python

# Functions for ICA-AROMA v0.3 beta

from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np

###RABIES modification
def run_ICA_AROMA(outDir,inFile,mc,TR,mask="",mask_csf="",denType="nonaggr",melDir="",dim=0,overwrite=False, random_seed=1):
    ###additional function for the execution of ICA_AROMA within RABIES
    import os
    #import subprocess
    from rabies.utils import run_command
    import shutil
    import rabies.confound_correction_pkg.mod_ICA_AROMA.classification_plots as classification_plots
    import rabies.confound_correction_pkg.mod_ICA_AROMA.ICA_AROMA_functions as aromafunc

    # Change to script directory
    cwd = os.path.realpath(os.path.curdir)
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(scriptDir)

    print('\n------------------------------- RUNNING ICA-AROMA ------------------------------- ')
    print('--------------- \'ICA-based Automatic Removal Of Motion Artifacts\' --------------- \n')

    # Define variables based on the type of input (i.e. Feat directory or specific input arguments), and check whether the specified files exist.
    cancel = False

    # Check whether the files exist
    if not inFile:
        print('No input file specified.')
    else:
        if not os.path.isfile(inFile):
            print('The specified input file does not exist.')
            cancel = True
    if not mc:
        print('No mc file specified.')
    else:
        if not os.path.isfile(mc):
            print('The specified mc file does does not exist.')
            cancel = True

    # Check if the mask exists, when specified.
    if mask:
        if not os.path.isfile(mask):
            print('The specified mask does not exist.')
            cancel = True

    # Check if the type of denoising is correctly specified, when specified
    if not (denType == 'nonaggr') and not (denType == 'aggr') and not (denType == 'both') and not (denType == 'no'):
        print('Type of denoising was not correctly specified. Non-aggressive denoising will be run.')
        denType = 'nonaggr'

    # If the criteria for file/directory specifications have not been met. Cancel ICA-AROMA.
    if cancel:
        print('\n----------------------------- ICA-AROMA IS CANCELED -----------------------------\n')
        exit()

    #------------------------------------------- PREPARE -------------------------------------------#

    # Define the FSL-bin directory
    fslDir = os.path.join(os.environ["FSLDIR"], 'bin', '')

    # Create output directory if needed
    if os.path.isdir(outDir) and overwrite is False:
        print('Output directory', outDir, """already exists.
              AROMA will not continue.
              Rerun with the -overwrite option to explicitly overwrite existing output.""")
        exit()
    elif os.path.isdir(outDir) and overwrite is True:
        print('Warning! Output directory', outDir, 'exists and will be overwritten.\n')
        shutil.rmtree(outDir)
        os.makedirs(outDir)
    else:
        os.makedirs(outDir)

    # Get TR of the fMRI data, if not specified
    if TR:
        TR = TR
    else:
        cmd = ' '.join([os.path.join(fslDir, 'fslinfo'),
                        inFile,
                        '| grep pixdim4 | awk \'{print $2}\''])
        rc,c_out = run_command(cmd)
        TR = float(c_out)

    # Check TR
    if TR == 0:
        print('TR is zero. ICA-AROMA requires a valid TR and will therefore exit. Please check the header, or define the TR as an additional argument.\n----------------------------- ICA-AROMA IS CANCELED -----------------------------\n')
        exit()

    # Define mask.
    mask_cp = os.path.join(outDir, 'mask.nii.gz')
    shutil.copyfile(mask, mask_cp)
    mask=mask_cp


    #---------------------------------------- Run ICA-AROMA ----------------------------------------#

    print('Step 1) MELODIC')
    try:
        aromafunc.runICA(fslDir, inFile, outDir, melDir, mask, dim, TR, random_seed=random_seed)
        melIC = os.path.join(outDir, 'melodic_IC_thr.nii.gz')
    except Exception as e:
        raise ValueError(f"MELODIC FAILED DURING ICA-AROMA: {e}")

    print('Step 2) Automatic classification of the components')

    print('  - *modified version skips commonspace registration')

    print('  - computing edge and out masks')
    mask_edge = os.path.join(outDir, 'mask_edge.nii.gz')
    mask_out = os.path.join(outDir, 'mask_out.nii.gz')
    aromafunc.compute_edge_mask(mask,mask_edge, num_edge_voxels=1)
    aromafunc.compute_out_mask(mask,mask_out)

    print('  - extracting the CSF & Edge fraction features')
    #modified inputs for the spatial features, by providing the required masks manually
    edgeFract, csfFract = aromafunc.mod_feature_spatial(fslDir, outDir, melIC, mask_csf, mask_edge, mask_out)

    print('  - extracting the Maximum RP correlation feature')
    melmix = os.path.join(outDir, 'melodic.ica', 'melodic_mix')
    maxRPcorr = aromafunc.feature_time_series(melmix, mc)

    print('  - extracting the High-frequency content feature')
    melFTmix = os.path.join(outDir, 'melodic.ica', 'melodic_FTmix')
    HFC = aromafunc.feature_frequency(melFTmix, TR)

    print('  - classification')
    motionICs = aromafunc.classification(outDir, maxRPcorr, edgeFract, HFC, csfFract)
    classification_plots.classification_plot(os.path.join(outDir, 'classification_overview.txt'),
                                             outDir)


    if (denType != 'no'):
        print('Step 3) Data denoising')
        denoising_applied = aromafunc.denoising(fslDir, inFile, outDir, melmix, denType, motionICs)
    else:
        denoising_applied=False

    # Revert to old directory
    os.chdir(cwd)

    print('\n----------------------------------- Finished -----------------------------------\n')
    return True, denoising_applied
###end of RABIES modification


def runICA(fslDir, inFile, outDir, melDirIn, mask, dim, TR, random_seed):
    """ This function runs MELODIC and merges the mixture modeled thresholded ICs into a single 4D nifti file

    Parameters
    ---------------------------------------------------------------------------------
    fslDir:     Full path of the bin-directory of FSL
    inFile:     Full path to the fMRI data file (nii.gz) on which MELODIC should be run
    outDir:     Full path of the output directory
    melDirIn:   Full path of the MELODIC directory in case it has been run before, otherwise define empty string
    mask:       Full path of the mask to be applied during MELODIC
    dim:        Dimensionality of ICA
    TR:     TR (in seconds) of the fMRI data

    Output (within the requested output directory)
    ---------------------------------------------------------------------------------
    melodic.ica     MELODIC directory
    melodic_IC_thr.nii.gz   merged file containing the mixture modeling thresholded Z-statistical maps located in melodic.ica/stats/ """

    # Import needed modules
    import os
    #import subprocess
    from rabies.utils import run_command

    # Define the 'new' MELODIC directory and predefine some associated files
    melDir = os.path.join(outDir, 'melodic.ica')
    melIC = os.path.join(melDir, 'melodic_IC.nii.gz')
    melICmix = os.path.join(melDir, 'melodic_mix')
    melICthr = os.path.join(outDir, 'melodic_IC_thr.nii.gz')

    # When a MELODIC directory is specified,
    # check whether all needed files are present.
    # Otherwise... run MELODIC again
    if len(melDir) != 0 and os.path.isfile(os.path.join(melDirIn, 'melodic_IC.nii.gz')) and os.path.isfile(os.path.join(melDirIn, 'melodic_FTmix')) and os.path.isfile(os.path.join(melDirIn, 'melodic_mix')):

        print('  - The existing/specified MELODIC directory will be used.')

        # If a 'stats' directory is present (contains thresholded spatial maps)
        # create a symbolic link to the MELODIC directory.
        # Otherwise create specific links and
        # run mixture modeling to obtain thresholded maps.
        if os.path.isdir(os.path.join(melDirIn, 'stats')):
            os.symlink(melDirIn, melDir)
        else:
            print('  - The MELODIC directory does not contain the required \'stats\' folder. Mixture modeling on the Z-statistical maps will be run.')

            # Create symbolic links to the items in the specified melodic directory
            os.makedirs(melDir)
            for item in os.listdir(melDirIn):
                os.symlink(os.path.join(melDirIn, item),
                           os.path.join(melDir, item))

            # Run mixture modeling
            run_command(' '.join([os.path.join(fslDir, 'melodic'),
                                '--in=' + melIC,
                                '--ICs=' + melIC,
                                '--mix=' + melICmix,
                                '--outdir=' + melDir,
                                '--Ostats --mmthresh=0.5 --seed='+str(random_seed)]))

    else:
        # If a melodic directory was specified, display that it did not contain all files needed for ICA-AROMA (or that the directory does not exist at all)
        if len(melDirIn) != 0:
            if not os.path.isdir(melDirIn):
                print('  - The specified MELODIC directory does not exist. MELODIC will be run seperately.')
            else:
                print('  - The specified MELODIC directory does not contain the required files to run ICA-AROMA. MELODIC will be run seperately.')

        # Run MELODIC
        run_command(' '.join([os.path.join(fslDir, 'melodic'),
                            '--in=' + inFile,
                            '--outdir=' + melDir,
                            '--mask=' + mask,
                            '--dim=' + str(dim),
                            '--Ostats --nobet --mmthresh=0.5 --report --seed='+str(random_seed),
                            '--tr=' + str(TR)]))

    # Get number of components
    cmd = ' '.join([os.path.join(fslDir, 'fslinfo'),
                    melIC,
                    '| grep dim4 | head -n1 | awk \'{print $2}\''])
    rc,c_out = run_command(cmd)
    nrICs = int(float(c_out))

    # Merge mixture modeled thresholded spatial maps. Note! In case that mixture modeling did not converge, the file will contain two spatial maps. The latter being the results from a simple null hypothesis test. In that case, this map will have to be used (first one will be empty).
    for i in range(1, nrICs + 1):
        # Define thresholded zstat-map file
        zTemp = os.path.join(melDir, 'stats', 'thresh_zstat' + str(i) + '.nii.gz')
        cmd = ' '.join([os.path.join(fslDir, 'fslinfo'),
                        zTemp,
                        '| grep dim4 | head -n1 | awk \'{print $2}\''])
        rc,c_out = run_command(cmd)
        lenIC = int(float(c_out))

        # Define zeropad for this IC-number and new zstat file
        cmd = ' '.join([os.path.join(fslDir, 'zeropad'),
                        str(i),
                        '4'])
        rc,c_out = run_command(cmd)
        ICnum = c_out
        zstat = os.path.join(outDir, 'thr_zstat' + ICnum)

        # Extract last spatial map within the thresh_zstat file
        run_command(' '.join([os.path.join(fslDir, 'fslroi'),
                            zTemp,      # input
                            zstat,      # output
                            str(lenIC - 1),   # first frame
                            '1']))      # number of frames

    # Merge and subsequently remove all mixture modeled Z-maps within the output directory
    run_command(' '.join([os.path.join(fslDir, 'fslmerge'),
                        '-t',                       # concatenate in time
                        melICthr,                   # output
                        os.path.join(outDir, 'thr_zstat????.nii.gz')]))  # inputs

    run_command('rm ' + os.path.join(outDir, 'thr_zstat????.nii.gz'))

    # Apply the mask to the merged file (in case a melodic-directory was predefined and run with a different mask)
    run_command(' '.join([os.path.join(fslDir, 'fslmaths'),
                        melICthr,
                        '-mas ' + mask,
                        melICthr]))


def register2MNI(fslDir, inFile, outFile, affmat, warp):
    """ This function registers an image (or time-series of images) to MNI152 T1 2mm. If no affmat is defined, it only warps (i.e. it assumes that the data has been registerd to the structural scan associated with the warp-file already). If no warp is defined either, it only resamples the data to 2mm isotropic if needed (i.e. it assumes that the data has been registered to a MNI152 template). In case only an affmat file is defined, it assumes that the data has to be linearly registered to MNI152 (i.e. the user has a reason not to use non-linear registration on the data).

    Parameters
    ---------------------------------------------------------------------------------
    fslDir:     Full path of the bin-directory of FSL
    inFile:     Full path to the data file (nii.gz) which has to be registerd to MNI152 T1 2mm
    outFile:    Full path of the output file
    affmat:     Full path of the mat file describing the linear registration (if data is still in native space)
    warp:       Full path of the warp file describing the non-linear registration (if data has not been registered to MNI152 space yet)

    Output (within the requested output directory)
    ---------------------------------------------------------------------------------
    melodic_IC_mm_MNI2mm.nii.gz merged file containing the mixture modeling thresholded Z-statistical maps registered to MNI152 2mm """


    # Import needed modules
    import os
    import subprocess

    # Define the MNI152 T1 2mm template
    fslnobin = fslDir.rsplit('/', 2)[0]
    ref = os.path.join(fslnobin, 'data', 'standard', 'MNI152_T1_2mm_brain.nii.gz')

    # If the no affmat- or warp-file has been specified, assume that the data is already in MNI152 space. In that case only check if resampling to 2mm is needed
    if (len(affmat) == 0) and (len(warp) == 0):
        # Get 3D voxel size
        pixdim1 = float(subprocess.getoutput('%sfslinfo %s | grep pixdim1 | awk \'{print $2}\'' % (fslDir, inFile)))
        pixdim2 = float(subprocess.getoutput('%sfslinfo %s | grep pixdim2 | awk \'{print $2}\'' % (fslDir, inFile)))
        pixdim3 = float(subprocess.getoutput('%sfslinfo %s | grep pixdim3 | awk \'{print $2}\'' % (fslDir, inFile)))

        # If voxel size is not 2mm isotropic, resample the data, otherwise copy the file
        if (pixdim1 != 2) or (pixdim2 != 2) or (pixdim3 != 2):
            os.system(' '.join([os.path.join(fslDir, 'flirt'),
                                ' -ref ' + ref,
                                ' -in ' + inFile,
                                ' -out ' + outFile,
                                ' -applyisoxfm 2 -interp trilinear']))
        else:
            os.system('cp ' + inFile + ' ' + outFile)

    # If only a warp-file has been specified, assume that the data has already been registered to the structural scan. In that case apply the warping without a affmat
    elif (len(affmat) == 0) and (len(warp) != 0):
        # Apply warp
        os.system(' '.join([os.path.join(fslDir, 'applywarp'),
                            '--ref=' + ref,
                            '--in=' + inFile,
                            '--out=' + outFile,
                            '--warp=' + warp,
                            '--interp=trilinear']))

    # If only a affmat-file has been specified perform affine registration to MNI
    elif (len(affmat) != 0) and (len(warp) == 0):
        os.system(' '.join([os.path.join(fslDir, 'flirt'),
                            '-ref ' + ref,
                            '-in ' + inFile,
                            '-out ' + outFile,
                            '-applyxfm -init ' + affmat,
                            '-interp trilinear']))

    # If both a affmat- and warp-file have been defined, apply the warping accordingly
    else:
        os.system(' '.join([os.path.join(fslDir, 'applywarp'),
                            '--ref=' + ref,
                            '--in=' + inFile,
                            '--out=' + outFile,
                            '--warp=' + warp,
                            '--premat=' + affmat,
                            '--interp=trilinear']))

def cross_correlation(a, b):
    """Cross Correlations between columns of two matrices"""
    assert a.ndim == b.ndim == 2
    _, ncols_a = a.shape
    # nb variables in columns rather than rows hence transpose
    # extract just the cross terms between cols in a and cols in b
    return np.corrcoef(a.T, b.T)[:ncols_a, ncols_a:]


def feature_time_series(melmix, mc):
    """ This function extracts the maximum RP correlation feature scores.
    It determines the maximum robust correlation of each component time-series
    with a model of 72 realignment parameters.

    Parameters
    ---------------------------------------------------------------------------------
    melmix:     Full path of the melodic_mix text file
    mc:     Full path of the text file containing the realignment parameters

    Returns
    ---------------------------------------------------------------------------------
    maxRPcorr:  Array of the maximum RP correlation feature scores for the components
    of the melodic_mix file"""

    # Import required modules
    import numpy as np
    import random

    # Read melodic mix file (IC time-series), subsequently define a set of squared time-series
    mix = np.loadtxt(melmix)
    if len(mix.shape)==1:
        mix = mix[:,np.newaxis]

    # Read motion parameter file
    rp6 = np.loadtxt(mc)
    _, nparams = rp6.shape

    # Determine the derivatives of the RPs (add zeros at time-point zero)
    rp6_der = np.vstack((np.zeros(nparams),
                         np.diff(rp6, axis=0)
                         ))

    # Create an RP-model including the RPs and its derivatives
    rp12 = np.hstack((rp6, rp6_der))

    # Add the squared RP-terms to the model
    # add the fw and bw shifted versions
    rp12_1fw = np.vstack((
        np.zeros(2 * nparams),
        rp12[:-1]
    ))
    rp12_1bw = np.vstack((
        rp12[1:],
        np.zeros(2 * nparams)
    ))
    rp_model = np.hstack((rp12, rp12_1fw, rp12_1bw))

    # Determine the maximum correlation between RPs and IC time-series
    nsplits = 1000
    nmixrows, nmixcols = mix.shape
    nrows_to_choose = int(round(0.9 * nmixrows))

    # Max correlations for multiple splits of the dataset (for a robust estimate)
    max_correls = np.empty((nsplits, nmixcols))
    for i in range(nsplits):
        # Select a random subset of 90% of the dataset rows (*without* replacement)
        chosen_rows = random.sample(population=range(nmixrows),
                                    k=nrows_to_choose)

        # Combined correlations between RP and IC time-series, squared and non squared
        correl_nonsquared = cross_correlation(mix[chosen_rows],
                                              rp_model[chosen_rows])
        correl_squared = cross_correlation(mix[chosen_rows]**2,
                                           rp_model[chosen_rows]**2)
        correl_both = np.hstack((correl_squared, correl_nonsquared))

        # Maximum absolute temporal correlation for every IC
        max_correls[i] = np.abs(correl_both).max(axis=1)

    # Feature score is the mean of the maximum correlation over all the random splits
    # Avoid propagating occasional nans that arise in artificial test cases
    return np.nanmean(max_correls, axis=0)


def feature_frequency(melFTmix, TR):
    """ This function extracts the high-frequency content feature scores.
    It determines the frequency, as fraction of the Nyquist frequency,
    at which the higher and lower frequencies explain half
    of the total power between 0.01Hz and Nyquist.

    Parameters
    ---------------------------------------------------------------------------------
    melFTmix:   Full path of the melodic_FTmix text file
    TR:     TR (in seconds) of the fMRI data (float)

    Returns
    ---------------------------------------------------------------------------------
    HFC:        Array of the HFC ('High-frequency content') feature scores
    for the components of the melodic_FTmix file"""

    # Import required modules
    import numpy as np

    # Determine sample frequency
    Fs = old_div(1, TR)

    # Determine Nyquist-frequency
    Ny = old_div(Fs, 2)

    # Load melodic_FTmix file
    FT = np.loadtxt(melFTmix)
    if len(FT.shape)==1:
        FT = FT[:,np.newaxis]

    # Determine which frequencies are associated with every row in the melodic_FTmix file  (assuming the rows range from 0Hz to Nyquist)
    f = Ny * (np.array(list(range(1, FT.shape[0] + 1)))) / (FT.shape[0])

    # Only include frequencies higher than 0.01Hz
    fincl = np.squeeze(np.array(np.where(f > 0.01)))
    FT = FT[fincl, :]
    f = f[fincl]

    # Set frequency range to [0-1]
    f_norm = old_div((f - 0.01), (Ny - 0.01))

    # For every IC; get the cumulative sum as a fraction of the total sum
    fcumsum_fract = old_div(np.cumsum(FT, axis=0), np.sum(FT, axis=0))

    # Determine the index of the frequency with the fractional cumulative sum closest to 0.5
    idx_cutoff = np.argmin(np.abs(fcumsum_fract - 0.5), axis=0)

    # Now get the fractions associated with those indices index, these are the final feature scores
    HFC = f_norm[idx_cutoff]

    # Return feature score
    return HFC

###RABIES modification
def compute_edge_mask(in_mask,out_file, num_edge_voxels):
    #custom function for computing edge mask from an input brain mask
    import numpy as np
    import nibabel as nb
    img=nb.load(in_mask)
    mask_array=np.asarray(img.dataobj)
    shape=mask_array.shape

    #iterate through all voxels from the three dimensions and look if it contains surrounding voxels
    edge_mask=np.zeros(shape, dtype=bool)
    num_voxel=0
    while num_voxel<num_edge_voxels:
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    #only look if the voxel is part of the mask
                    if mask_array[x,y,z]:
                        if (mask_array[x-1:x+2,y-1:y+2,z-1:z+2]==0).sum()>0:
                            edge_mask[x,y,z]=1
        mask_array=mask_array-edge_mask
        num_voxel+=1

    nb.Nifti1Image(edge_mask, img.affine, img.header).to_filename(out_file)

def compute_out_mask(in_mask,out_file):
    #custom function for computing a mask for the outside of the brain
    import numpy as np
    import nibabel as nb
    img=nb.load(in_mask)
    mask_array=np.asarray(img.dataobj)

    #iterate through all voxels from the three dimensions and look if it contains surrounding voxels
    out_mask=np.asarray(mask_array==0,dtype=bool)
    nb.Nifti1Image(out_mask, img.affine, img.header).to_filename(out_file)

def mod_feature_spatial(fslDir, tempDir, melIC, mask_csf, mask_edge, mask_out):
    #This is a modified version of the orginial ICA-AROMA function where the CSF and edge masks are provided manually.
    """
    This function extracts the spatial feature scores. For each IC it determines the fraction of the mixture modeled thresholded Z-maps respecitvely located within the CSF or at the brain edges, using predefined standardized masks.

    Parameters
    ---------------------------------------------------------------------------------
    fslDir:     Full path of the bin-directory of FSL
    tempDir:    Full path of a directory where temporary files can be stored (called 'temp_IC.nii.gz')
    melIC:      Full path of the nii.gz file containing mixture-modeled threholded (p>0.5) Z-maps, which overlays with the provided masks.

    Returns
    ---------------------------------------------------------------------------------
    edgeFract:  Array of the edge fraction feature scores for the components of the melIC file
    csfFract:   Array of the CSF fraction feature scores for the components of the melIC file"""

    # Import required modules
    import numpy as np
    import os
    #import subprocess
    from rabies.utils import run_command

    # Get the number of ICs
    cmd = '%sfslinfo %s | grep dim4 | head -n1 | awk \'{print $2}\'' % (fslDir, melIC)
    rc,c_out = run_command(cmd)
    numICs = int(c_out)

    # Loop over ICs
    edgeFract = np.zeros(numICs)
    csfFract = np.zeros(numICs)
    for i in range(0, numICs):
        # Define temporary IC-file
        tempIC = os.path.join(tempDir, 'temp_IC.nii.gz')

        # Extract IC from the merged melodic_IC_thr2MNI2mm file
        run_command(' '.join([os.path.join(fslDir, 'fslroi'),
                  melIC,
                  tempIC,
                  str(i),
                  '1']))

        # Change to absolute Z-values
        run_command(' '.join([os.path.join(fslDir, 'fslmaths'),
                  tempIC,
                  '-abs',
                  tempIC]))

        # Get sum of Z-values within the total Z-map (calculate via the mean and number of non-zero voxels)
        cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                    tempIC,
                                                    '-V | awk \'{print $1}\''])
        rc,c_out = run_command(cmd)
        totVox = int(c_out)

        if not (totVox == 0):
            cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                           tempIC,
                                                           '-M'])
            rc,c_out = run_command(cmd)
            totMean = float(c_out)
        else:
            print('     - The spatial map of component ' + str(i + 1) + ' is empty. Please check!')
            totMean = 0

        totSum = totMean * totVox

        # Get sum of Z-values of the voxels located within the CSF (calculate via the mean and number of non-zero voxels)
        cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                    tempIC,
                                                    '-k ', mask_csf,
                                                    '-V | awk \'{print $1}\''])
        rc,c_out = run_command(cmd)
        csfVox = int(c_out)

        if not (csfVox == 0):
            cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                           tempIC,
                                                           '-k ', mask_csf,
                                                           '-M'])
            rc,c_out = run_command(cmd)
            csfMean = float(c_out)
        else:
            csfMean = 0

        csfSum = csfMean * csfVox

        # Get sum of Z-values of the voxels located within the Edge (calculate via the mean and number of non-zero voxels)
        cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                     tempIC,
                                                     '-k ', mask_edge,
                                                     '-V | awk \'{print $1}\''])
        rc,c_out = run_command(cmd)
        edgeVox = int(c_out)
        if not (edgeVox == 0):
            cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                            tempIC,
                                                            '-k ', mask_edge,
                                                            '-M'])
            rc,c_out = run_command(cmd)
            edgeMean = float(c_out)
        else:
            edgeMean = 0

        edgeSum = edgeMean * edgeVox

        # Get sum of Z-values of the voxels located outside the brain (calculate via the mean and number of non-zero voxels)
        cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                    tempIC,
                                                    '-k ', mask_out,
                                                    '-V | awk \'{print $1}\''])
        rc,c_out = run_command(cmd)
        outVox = int(c_out)
        if not (outVox == 0):
            cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                           tempIC,
                                                           '-k ', mask_out,
                                                           '-M'])
            rc,c_out = run_command(cmd)
            outMean = float(c_out)
        else:
            outMean = 0

        outSum = outMean * outVox

        # Determine edge and CSF fraction
        if totSum==csfSum:
            edgeFract[i] = 1
            csfFract[i] = 1
        elif not (totSum == 0):
            edgeFract[i] = old_div((outSum + edgeSum), (totSum - csfSum))
            csfFract[i] = old_div(csfSum, totSum)
        else:
            edgeFract[i] = 0
            csfFract[i] = 0

    # Remove the temporary IC-file
    os.remove(tempIC)

    # Return feature scores
    return edgeFract, csfFract
###end of RABIES modification


def feature_spatial(fslDir, tempDir, aromaDir, melIC):
    """ This function extracts the spatial feature scores. For each IC it determines the fraction of the mixture modeled thresholded Z-maps respecitvely located within the CSF or at the brain edges, using predefined standardized masks.

    Parameters
    ---------------------------------------------------------------------------------
    fslDir:     Full path of the bin-directory of FSL
    tempDir:    Full path of a directory where temporary files can be stored (called 'temp_IC.nii.gz')
    aromaDir:   Full path of the ICA-AROMA directory, containing the mask-files (mask_edge.nii.gz, mask_csf.nii.gz & mask_out.nii.gz)
    melIC:      Full path of the nii.gz file containing mixture-modeled threholded (p>0.5) Z-maps, registered to the MNI152 2mm template

    Returns
    ---------------------------------------------------------------------------------
    edgeFract:  Array of the edge fraction feature scores for the components of the melIC file
    csfFract:   Array of the CSF fraction feature scores for the components of the melIC file"""

    # Import required modules
    import numpy as np
    import os
    #import subprocess
    from rabies.utils import run_command

    # Get the number of ICs
    cmd = '%sfslinfo %s | grep dim4 | head -n1 | awk \'{print $2}\'' % (fslDir, melIC)
    rc,c_out = run_command(cmd)
    numICs = int(c_out)

    # Loop over ICs
    edgeFract = np.zeros(numICs)
    csfFract = np.zeros(numICs)
    for i in range(0, numICs):
        # Define temporary IC-file
        tempIC = os.path.join(tempDir, 'temp_IC.nii.gz')

        # Extract IC from the merged melodic_IC_thr2MNI2mm file
        run_command(' '.join([os.path.join(fslDir, 'fslroi'),
                  melIC,
                  tempIC,
                  str(i),
                  '1']))

        # Change to absolute Z-values
        run_command(' '.join([os.path.join(fslDir, 'fslmaths'),
                  tempIC,
                  '-abs',
                  tempIC]))

        # Get sum of Z-values within the total Z-map (calculate via the mean and number of non-zero voxels)
        cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                    tempIC,
                                                    '-V | awk \'{print $1}\''])
        rc,c_out = run_command(cmd)
        totVox = int(c_out)

        if not (totVox == 0):
            cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                           tempIC,
                                                           '-M'])
            rc,c_out = run_command(cmd)
            totMean = float(c_out)
        else:
            print('     - The spatial map of component ' + str(i + 1) + ' is empty. Please check!')
            totMean = 0

        totSum = totMean * totVox

        # Get sum of Z-values of the voxels located within the CSF (calculate via the mean and number of non-zero voxels)
        cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                           tempIC,
                                                           '-k mask_csf.nii.gz',
                                                           '-M'])
        rc,c_out = run_command(cmd)
        csfVox = int(c_out)

        if not (csfVox == 0):
            cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                           tempIC,
                                                           '-k mask_csf.nii.gz',
                                                           '-M'])
            rc,c_out = run_command(cmd)
            csfMean = float(c_out)
        else:
            csfMean = 0

        csfSum = csfMean * csfVox

        # Get sum of Z-values of the voxels located within the Edge (calculate via the mean and number of non-zero voxels)
        cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                     tempIC,
                                                     '-k mask_edge.nii.gz',
                                                     '-V | awk \'{print $1}\''])
        rc,c_out = run_command(cmd)
        edgeVox = int(c_out)
        if not (edgeVox == 0):
            cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                            tempIC,
                                                            '-k mask_edge.nii.gz',
                                                            '-M'])
            rc,c_out = run_command(cmd)
            edgeMean = float(c_out)
        else:
            edgeMean = 0

        edgeSum = edgeMean * edgeVox

        # Get sum of Z-values of the voxels located outside the brain (calculate via the mean and number of non-zero voxels)
        cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                    tempIC,
                                                    '-k mask_out.nii.gz',
                                                    '-V | awk \'{print $1}\''])
        rc,c_out = run_command(cmd)
        outVox = int(c_out)
        if not (outVox == 0):
            cmd = ' '.join([os.path.join(fslDir, 'fslstats'),
                                                           tempIC,
                                                           '-k mask_out.nii.gz',
                                                           '-M'])
            rc,c_out = run_command(cmd)
            outMean = float(c_out)
        else:
            outMean = 0

        outSum = outMean * outVox

        # Determine edge and CSF fraction
        if not (totSum == 0):
            edgeFract[i] = old_div((outSum + edgeSum), (totSum - csfSum))
            csfFract[i] = old_div(csfSum, totSum)
        else:
            edgeFract[i] = 0
            csfFract[i] = 0

    # Remove the temporary IC-file
    os.remove(tempIC)

    # Return feature scores
    return edgeFract, csfFract


def classification(outDir, maxRPcorr, edgeFract, HFC, csfFract):
    """ This function classifies a set of components into motion and
    non-motion components based on four features;
    maximum RP correlation, high-frequency content, edge-fraction and CSF-fraction

    Parameters
    ---------------------------------------------------------------------------------
    outDir:     Full path of the output directory
    maxRPcorr:  Array of the 'maximum RP correlation' feature scores of the components
    edgeFract:  Array of the 'edge fraction' feature scores of the components
    HFC:        Array of the 'high-frequency content' feature scores of the components
    csfFract:   Array of the 'CSF fraction' feature scores of the components

    Return
    ---------------------------------------------------------------------------------
    motionICs   Array containing the indices of the components identified as motion components

    Output (within the requested output directory)
    ---------------------------------------------------------------------------------
    classified_motion_ICs.txt   A text file containing the indices of the components identified as motion components """

    # Import required modules
    import numpy as np
    import os

    # Classify the ICs as motion or non-motion

    # Define criteria needed for classification (thresholds and hyperplane-parameters)
    thr_csf = 0.10
    thr_HFC = 0.60 #modified threshold for RABIES
    hyp = np.array([-14, 9.95127547670627, 24.8333160239175]) #modified hyperparameters for RABIES

    # Project edge & maxRPcorr feature scores to new 1D space
    x = np.array([maxRPcorr, edgeFract])
    proj = hyp[0] + np.dot(x.T, hyp[1:])

    # Classify the ICs
    motionICs = np.squeeze(np.array(np.where((proj > 0) + (csfFract > thr_csf) + (HFC > thr_HFC))))

    # Put the feature scores in a text file
    np.savetxt(os.path.join(outDir, 'feature_scores.txt'),
               np.vstack((maxRPcorr, edgeFract, HFC, csfFract)).T)

    # Put the indices of motion-classified ICs in a text file
    txt = open(os.path.join(outDir, 'classified_motion_ICs.txt'), 'w')
    if motionICs.size > 1:  # and len(motionICs) != 0: if motionICs is not None and
        txt.write(','.join(['{:.0f}'.format(num) for num in (motionICs + 1)]))
    elif motionICs.size == 1:
        txt.write('{:.0f}'.format(motionICs + 1))
    txt.close()

    # Create a summary overview of the classification
    txt = open(os.path.join(outDir, 'classification_overview.txt'), 'w')
    txt.write('\t'.join(['IC',
                         'Motion/noise',
                         'maximum RP correlation',
                         'Edge-fraction',
                         'High-frequency content',
                         'CSF-fraction']))
    txt.write('\n')
    for i in range(0, len(csfFract)):
        if (proj[i] > 0) or (csfFract[i] > thr_csf) or (HFC[i] > thr_HFC):
            classif = "True"
        else:
            classif = "False"
        txt.write('\t'.join(['{:d}'.format(i + 1),
                             classif,
                             '{:.2f}'.format(maxRPcorr[i]),
                             '{:.2f}'.format(edgeFract[i]),
                             '{:.2f}'.format(HFC[i]),
                             '{:.2f}'.format(csfFract[i])]))
        txt.write('\n')
    txt.close()

    return motionICs


def denoising(fslDir, inFile, outDir, melmix, denType, denIdx):
    """ This function classifies the ICs based on the four features;
    maximum RP correlation, high-frequency content, edge-fraction and CSF-fraction

    Parameters
    ---------------------------------------------------------------------------------
    fslDir:     Full path of the bin-directory of FSL
    inFile:     Full path to the data file (nii.gz) which has to be denoised
    outDir:     Full path of the output directory
    melmix:     Full path of the melodic_mix text file
    denType:    Type of requested denoising ('aggr': aggressive, 'nonaggr': non-aggressive, 'both': both aggressive and non-aggressive
    denIdx:     Indices of the components that should be regressed out

    Output (within the requested output directory)
    ---------------------------------------------------------------------------------
    denoised_func_data_<denType>.nii.gz:        A nii.gz file of the denoised fMRI data"""

    # Import required modules
    import os
    import numpy as np
    from rabies.utils import run_command

    # Check if denoising is needed (i.e. are there components classified as motion)
    check = denIdx.size > 0

    if check == 1:
        # Put IC indices into a char array
        if denIdx.size == 1:
            denIdxStrJoin = "%d"%(denIdx + 1)
        else:
            denIdxStr = np.char.mod('%i', (denIdx + 1))
            denIdxStrJoin = ','.join(denIdxStr)

        # Non-aggressive denoising of the data using fsl_regfilt (partial regression), if requested
        if (denType == 'nonaggr') or (denType == 'both'):
            run_command(' '.join([os.path.join(fslDir, 'fsl_regfilt'),
                                '--in=' + inFile,
                                '--design=' + melmix,
                                '--filter="' + denIdxStrJoin + '"',
                                '--out=' + os.path.join(outDir, 'denoised_func_data_nonaggr.nii.gz')]))

        # Aggressive denoising of the data using fsl_regfilt (full regression)
        if (denType == 'aggr') or (denType == 'both'):
            run_command(' '.join([os.path.join(fslDir, 'fsl_regfilt'),
                                '--in=' + inFile,
                                '--design=' + melmix,
                                '--filter="' + denIdxStrJoin + '"',
                                '--out=' + os.path.join(outDir, 'denoised_func_data_aggr.nii.gz'),
                                '-a']))
        return True # we return whether denoising was applied or not
    else:
        print("  - None of the components were classified as motion, so no denoising is applied.")
        return False # we return whether denoising was applied or not
