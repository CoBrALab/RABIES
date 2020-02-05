## Containerized installation (Recommended)
Containers are virtual computing environments where all dependencies required to execute a software can be controlled. Thus, containers allow to run RABIES without having to install and control dependencies manually. However, containers can only be executed in parallel through the 'MultiProc' option.
* Singularity container: RABIES can be installed as a Singularity container to run locally or on remote computing clusters. After installing Singularity (https://singularity.lbl.gov), the container can be generated as a singularity image with "singularity build rabies.sif docker://gabdesgreg/rabies:tagname". See README for execution instructions.
* Docker container: similarly to Singularity, a Docker version is also available (https://www.docker.com). You can directly install the container from docker hub with "docker pull gabdesgreg/rabies:tagname", to create a new image tagged "gabdesgreg/rabies:tagname" among your docker images. To test the installation, run "docker run gabdesgreg/rabies:tagname -h" and the help message from the command line interface should be printed. The usage of Docker requires root privileges, and thus cannot be used on remote computing clusters. See README for execution instructions.

## Local installation

* Minimal Dependencies: before installing RABIES, the minimal dependencies listed in dependencies.txt must be functional. You might want to refer to the Dockerfile installation build to get started.
* python environment: we highly recommend installing a anaconda environment for the python dependencies (also listed in dependencies.txt). You can do so by installing Miniconda3 from https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html for linux. A conda environment with all required python dependencies can be created from the rabies_environment.yml file by running "conda env create -f rabies_environment.yml". The environment must be activated through "source activate rabies" before running RABIES.
* install.sh: After the dependencies are met, you can execute install.sh to install RABIES on the local computer.
<br/>
1. will install the RABIES package in $HOME directory, creates an executable script, and setup permanent paths in .bashrc
<br/>
2. Downloads the DSURQE mouse template as default template, and generates white matter and CSF masks from gen_masks.py
<br/>
3. Also installs https://github.com/CobraLab/twolevel_ants_dbm
