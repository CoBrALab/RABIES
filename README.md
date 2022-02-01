# RABIES: Rodent Automated Bold Improvement of EPI Sequences.

RABIES is an open source image processing pipeline for rodent fMRI. It conducts state-of-the-art preprocessing and confound correction, and supplies standard resting-state functional connectivity analyses. Visit our documentation at <https://rabies.readthedocs.io/en/stable/>.

![RABIES Schema](https://raw.githubusercontent.com/CoBrALab/RABIES/master/docs/pics/RABIES_schema.png)

## What you can do with RABIES

The primary purpose of RABIES consists of equipping independent rodent fMRI research with a standard, yet flexible and reliable, image processing platform, 
and complement these tools with guidelines supporting best practices for reproducibility and quality control. The RABIES software is structured into three main processing stages: **preprocessing**, **confound correction** and **analysis**. 

### Preprocessing
The preprocessing workflow regroups essential fMRI preprocessing steps prior to analysis. It includes a robust registration workflow with automatically-adapting parameters allowing to succesfully process diverse acquisition types (i.e. rodent species, scan field strenght, coil type, ...), and can conduct the following preprocessing steps:
- head motion correction
- susceptibility distortion correction
- resampling to native or common space
- brain parcellation
- slice timing correction (optional)
- despiking (optional)
- visual assessment of registration for quality control

### Confound correction
Following preprocessing, a range of strategies to correct fMRI confounds (e.g. motion) can then be conducted within RABIES:
- confound regression (with several options for nuisance regressors)
- frequency filtering (highpass, lowpass, bandpass)
- frame censoring (or scrubbing)
- ICA-AROMA
- spatial smoothing

### Analysis
Simple resting-state connectivity analyses are made available after preprocessing and confound correction. RABIES also provides a 'data diagnosis' workflow, which generates several indices of data quality and potential confounds, and conversaly, aims to improve the correction of confounds and transparency with regards to data quality:
- seed-based functional connectivity
- whole-brain connectivity matrix
- group-ICA
- dual regression
- data diagnosis


## Notes on software design

**Nipype workflows**: The image processing pipelines are structured using the [Nipype library](https://nipype.readthedocs.io/en/latest/), which allows to build dynamic workflows in the form of a computational graph. Each node in the graph consists of a processing step, and the required inputs/outputs define the links between nodes.

**Reproducible and transparent research**: RABIES aims to follow best practices for reproducible and transparent research, including the following:
- open source code <https://github.com/CoBrALab/RABIES>
- standardized input data format with [BIDS](https://bids.neuroimaging.io/)
- easily shared, automatically-generated visual outputs for quality control
- containerized distribution of the software hosted on [Docker Hub](https://hub.docker.com/r/gabdesgreg/rabies) which can be downloaded via Docker and Singularity platforms

## Citation

**Citing RABIES**: to cite the software, we currently ask to reference its [Github page](https://github.com/CoBrALab/RABIES).

**Boilerplate**: a boilerplate summarizing the preprocessing and confound correction operations is automatically generated in the output folder. You can use the boilerplate as reference to described your methods, or reproduce it verbatim, in a paper.

## License
The [RABIES license](https://github.com/CoBrALab/RABIES/blob/master/LICENSE) allows for uses in academic and educational environments only. Commercial use requires a commercial license from CoBrALab <contact@cobralab.ca>, <http://cobralab.ca>

## Acknowledgements
This software was developped by the [CoBrALab](https://cobralab.ca/), located at the Cerebral Imaging Center of the Douglas Mental Health University Institute, Montreal, Canada, in affiliation with McGill University, Montreal, Canada. This work was supported by funding from Healthy Brains, Healthy Lives (HBHL), the Fonds de recherche du Québec - Santé (FRQS) and - Nature et technologies (FRQNT), and the Natural Sciences and Engineering Research Council (NSERC) of Canada. [fMRIPrep](https://fmriprep.org/en/stable/) was an important inspirational source for this project, in particular with regards to best practices for software reproducibility and code design using Nipype. We also thank the organizers of [BrainHack School Montreal](https://school.brainhackmtl.org/), which guided the initial steps of this project in 2018.


## Ask for help
If you need support in using the software or experience issues that are not documented, we'll provide support on the [Github discussion](https://github.com/CoBrALab/RABIES/discussions).

## Contributing to RABIES

RABIES is under continuous development to keep up with evolving demand and ongoing research. We welcome suggestions for improvements using the [Github issue system](https://github.com/CoBrALab/RABIES/issues). If you're interested in contributing code, you can reach out on the [Github discussion](https://github.com/CoBrALab/RABIES/discussions) and we welcome contributions as pull requests.
