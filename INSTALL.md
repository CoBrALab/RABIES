## Details on the installation

* Minimal Dependencies: before installing RABIES, the minimal dependencies listed in dependencies.txt must be functional. You might want to refer to the Dockerfile installation build to get started.
* python environment: we highly recommend installing a anaconda environment for the python dependencies (also listed in dependencies.txt). You can do so by installing Miniconda3 from https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html for linux.
    * rabies_environment.yml: A new "rabies" conda environment can be created directly from the rabies_environment.yml file by running "conda env create -f rabies_environment.yml".
* install.sh: After the dependencies are met, you can execute install.sh to install RABIES on the local computer. 1. will install the RABIES package in $HOME directory, creates an executable script, and setup permanent paths in .bashrc 2. Downloads the DSURQE mouse template as default template, and generates white matter and CSF masks from gen_masks.py 3. Also installs https://github.com/CobraLab/twolevel_ants_dbm
    * gen_masks.py requires numpy, pandas and nibabel to run
* Docker container: as an alternative to a local installation, a Docker version is also available to run a containerized version of RABIES, which creates a complete virtual environment that meets all dependencies.