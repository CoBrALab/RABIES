### initial setup from https://github.com/BIC-MNI/build_packages/blob/master/build_ubuntu_18.04_x64/Dockerfile
FROM ubuntu:bionic

# install basic system packages
RUN apt-get -y update && \
    apt-get -y dist-upgrade && \
    apt-get install -y --no-install-recommends \
         sudo \
         build-essential g++ gfortran bc \
         bison flex \
         libx11-dev x11proto-core-dev \
         libxi6 libxi-dev \
         libxmu6 libxmu-dev libxmu-headers \
         libgl1-mesa-dev libglu1-mesa-dev \
         libjpeg-dev \
         libssl-dev ccache libapt-inst2.0 git lsb-release \
         curl ca-certificates && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/*

# add user to build all tools
RUN useradd -ms /bin/bash rabies && \
    echo "rabies ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/rabies && \
    chmod 0440 /etc/sudoers.d/rabies

ENV PATH=/usr/lib/ccache:$PATH

# add new cmake
RUN mkdir src && \
    cd src && \
    curl -L --output cmake-3.14.5.tar.gz https://github.com/Kitware/CMake/releases/download/v3.14.5/cmake-3.14.5.tar.gz  && \
    tar zxf cmake-3.14.5.tar.gz && \
    cd cmake-3.14.5 && \
    ./configure --prefix=/usr --no-qt-gui && \
    make && \
    make install && \
    cd ../../ && \
    rm -rf src


WORKDIR /home/rabies
ENV HOME="/home/rabies"

### install ANTs/AFNI softwares
RUN apt-get update -qq \
    && sudo apt-get install -y -q --no-install-recommends \
           gcc \
           g++ \
           graphviz \
           tree \
           git \
           vim \
           emacs-nox \
           nano \
           less \
           ncdu \
           tig \
           git-annex-remote-rclone \
           netbase \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Installing ANTs 2.3.1 based on neurodocker
ENV ANTSPATH="$HOME/ants-v2.3.1/bin" \
    PATH="$HOME/ants-v2.3.1/bin:$PATH" \
    LD_LIBRARY_PATH="$HOME/ants-v2.3.1/lib:$LD_LIBRARY_PATH"
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           cmake \
           g++ \
           gcc \
           git \
           make \
           zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /tmp/ants/build \
    && git clone https://github.com/ANTsX/ANTs.git /tmp/ants/source \
    && cd /tmp/ants/source \
    && git fetch --tags \
    && git checkout v2.3.1 \
    && cd /tmp/ants/build \
    && cmake -DBUILD_SHARED_LIBS=ON /tmp/ants/source \
    && make -j1 \
    && mkdir -p $HOME/ants-v2.3.1 \
    && mv bin lib $HOME/ants-v2.3.1/ \
    && mv /tmp/ants/source/Scripts/* $HOME/ants-v2.3.1/bin \
    && rm -rf /tmp/ants


# install AFNI latest
RUN wget afni.nimh.nih.gov/pub/dist/tgz/linux_ubuntu_16_64.tgz \
    && mkdir $HOME/afni \
    && tar -xzf linux_ubuntu_16_64.tgz -C $HOME/afni \
    && rm -rf linux_ubuntu_16_64.tgz

ENV AFNIPATH="$HOME/afni/linux_ubuntu_16_64/" \
    PATH="$HOME/afni/linux_ubuntu_16_64:$PATH"

RUN apt-get update && \
  apt-get install -y --no-install-recommends htop nano wget imagemagick parallel zram-config debconf

#Build tools and dependencies
RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections && \
  apt install -y --no-install-recommends build-essential gdebi-core \
    git imagemagick libssl-dev cmake autotools-dev automake \
    ed zlib1g-dev libxml2-dev libxslt-dev openjdk-8-jre \
    zenity libcurl4-openssl-dev bc gawk libxkbcommon-x11-0 \
    ttf-mscorefonts-installer bc

#Install python environment

ENV CONDA_DIR="$HOME/miniconda-latest" \
    PATH="$HOME/miniconda-latest/bin:$PATH" \
    ND_ENTRYPOINT="$HOME/startup.sh"

RUN export PATH="$HOME/miniconda-latest/bin:$PATH" \
    && echo "Downloading Miniconda installer ..." \
    && conda_installer="/tmp/miniconda.sh" \
    && curl -fsSL --retry 5 -o "$conda_installer" https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash "$conda_installer" -b -p $HOME/miniconda-latest \
    && rm -f "$conda_installer" \
    && conda update -yq -nbase conda


#### install RABIES
ENV export RABIES_VERSION=0.1.1 \
    export RABIES=$HOME/RABIES-${RABIES_VERSION} \
    export PYTHONPATH="${PYTHONPATH}:$RABIES"

RUN export RABIES_VERSION=0.1.1 && \
  export RABIES=$HOME/RABIES-${RABIES_VERSION} && \
  mkdir -p temp && \
  curl -L --retry 5 -o temp/RABIES.tar.gz https://github.com/CoBrALab/RABIES/archive/${RABIES_VERSION}.tar.gz && \
  cd temp && \
  tar zxf RABIES.tar.gz && \
  cd .. && \
  conda env create -f temp/RABIES-${RABIES_VERSION}/rabies_environment.yml && \
  bash temp/RABIES-${RABIES_VERSION}/install.sh && \
  rm -r temp && \
  DSURQE_100micron_labels=${RABIES}/template_files/DSURQE_100micron_labels.nii.gz && \
  csv_labels=${RABIES}/template_files/DSURQE_40micron_R_mapping.csv && \
  /home/rabies/miniconda-latest/envs/rabies/bin/python ${RABIES}/gen_masks.py $DSURQE_100micron_labels $csv_labels ${RABIES}/template_files/DSURQE_100micron && \
  echo "#! /home/rabies/miniconda-latest/envs/rabies/bin/python" > ${RABIES}/exec.py && \
  echo "import os" >> ${RABIES}/exec.py && \
  echo "import sys" >> ${RABIES}/exec.py && \
  echo "os.environ['RABIES'] = '${RABIES}'" >> ${RABIES}/exec.py && \
  echo "sys.path.insert(0,os.environ['RABIES'])" >> ${RABIES}/exec.py && \
  echo "from rabies.run_main import execute_workflow" >> ${RABIES}/exec.py && \
  echo "execute_workflow()" >> ${RABIES}/exec.py && \
  chmod +x ${RABIES}/exec.py

ENV QBATCH_SYSTEM local

WORKDIR /tmp/
ENTRYPOINT ["/home/rabies/RABIES-0.1.1/exec.py"]
