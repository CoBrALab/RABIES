# How to install RABIES

Pick the installation route that matches where you will run RABIES.

::::{grid} 1 1 3 3
:gutter: 2

:::{grid-item-card} Container
**Recommended.** Apptainer on Linux and HPC clusters, Docker elsewhere. All
dependencies included.
:::

:::{grid-item-card} PyPI
Python package only. You install the non-Python dependencies yourself.
:::

:::{grid-item-card} Neurodesk
Browser-based neuroimaging environment. Nothing to install locally.
:::

::::

## Install as a container

Containers package the entire computing environment, so you do not install
dependencies by hand and the software behaves identically wherever you run it.
[Apptainer](https://apptainer.org/) is generally preferred over
[Docker](https://www.docker.com) because it does not require root permissions,
which makes it usable on high performance computing clusters.

Install Apptainer or Docker first (Apptainer publishes
[quick start guidelines](https://apptainer.org/docs/user/main/quick_start.html)),
then pull the [RABIES image](https://github.com/CoBrALab/RABIES/pkgs/container/rabies):

::::{tab-set}

:::{tab-item} Apptainer
:sync: apptainer

```sh
apptainer build rabies-latest.sif docker://ghcr.io/cobralab/rabies:latest
```

This produces a single `.sif` file containing the whole environment.
:::

:::{tab-item} Docker
:sync: docker

```sh
docker pull ghcr.io/cobralab/rabies:latest
```
:::

::::

To pin a version, replace `latest` with a tag from the
[list of published images](https://github.com/CoBrALab/RABIES/pkgs/container/rabies).

```{note}
Versions prior to 0.5.0 are not on the GitHub container registry. They remain
available on [Docker Hub](https://hub.docker.com/r/gabdesgreg/rabies).
```

For the execution syntax once the image is built, see
[How to handle container syntax](run_with_containers.md).

## Install from PyPI

RABIES is published on [PyPI](https://pypi.org/project/rabies/):

```sh
pip install rabies
```

```{warning}
`pip install` gives you the Python package only. RABIES also calls out to
external neuroimaging tools, which are listed in
[`dependencies.txt`](https://github.com/CoBrALab/RABIES/blob/master/dependencies.txt)
and which you must install yourself. If you are not prepared to manage those,
use a container.
```

## Use RABIES on Neurodesk

RABIES is one of the [built-in tools](https://neurodesk.github.io/applications/)
on the [Neurodesk platform](https://neurodesk.github.io/), a browser-based
neuroimaging computing environment with community-maintained prebuilt tools.
Nothing is installed on your own machine. See the
[Neurodesk documentation](https://neurodesk.github.io/docs/) to get started.
