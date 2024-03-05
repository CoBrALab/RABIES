# Installation

## Container (Singularity/Docker) \*\*RECOMMENDED\*\*
For most uses, we recommend instead using a containerized installation with [Singularity](https://singularity.lbl.gov) or [Docker](https://www.docker.com). Containers allow to build entire computing environments, grouping all dependencies required to run the software. This in turn reduces the burden of installing dependencies manually and ensures reproducible behavior of the software. Singularity is generally preferred over Docker since it requires less permission, and can thus be imported from most computing environment (e.g. high performance computing clusters such as Compute Canada.)

A [containerized version](https://github.com/CoBrALab/RABIES/pkgs/container/rabies) of RABIES is available from Github. After installing Singularity or Docker, the following command will pull and build the container:
* Install Singularity .sif file: 
```
singularity build rabies.sif docker://ghcr.io/cobralab/rabies:latest
```
* Install Docker image: 
```
docker pull ghcr.io/cobralab/rabies:latest
```
A specific tag version can be selected (instead of `latest`) from the [list online](https://github.com/CoBrALab/RABIES/pkgs/container/rabies). Versions prior to 0.5.0 are found on [Docker Hub](https://hub.docker.com/r/gabdesgreg/rabies).

## PyPi
The software is available on [PyPi](https://pypi.org/project/rabies/), which makes the rabies python package widely accessible with
```
pip install rabies
```
However, this does not account for non-python dependencies found in `dependencies.txt`.

## Neurodesk
RABIES is also made available on the [Neurodesk platform](https://neurodesk.github.io/), as part of the [built-in tools](https://neurodesk.github.io/applications/) for neuroimaging. The Neurodesk platform allows for an entirely browser-based neuroimaging computing environment, with pre-built neuroimaging tools from the community, and aims at reducing needs for manual development of computing environments and at improving reproducible neuroimaging. More details on Neurodesk here <https://neurodesk.github.io/docs/>. 