## Details on the installation

* install.sh: 1. will install the RABIES package in $HOME directory, creates an executable script, and setup permanent paths in .bashrc 2. Downloads the DSURQE mouse template as default template, and generates white matter and CSF masks from gen_masks.py 3. Also installs https://github.com/CobraLab/twolevel_ants_dbm
    * gen_masks.py requires numpy, pandas and nibabel to run
* conda environment: can use the rabies_environment.yml to directly create the environment by running "conda env create -f rabies_environment.yml". Should append channels conda-forge, bioconda, simpleitk ("conda config --append channels conda-forge")
* Software Dependencies: can refer to the Dockerfile for installing dependencies.
