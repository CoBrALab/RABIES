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


# install FSL
RUN sudo apt-get update && \
  sudo apt-get install -y --no-install-recommends gnupg gnupg2 gnupg1
RUN git clone https://github.com/poldracklab/mriqc.git && \
  mv mriqc/docker/files/neurodebian.gpg $HOME && \
  rm -rf mriqc

RUN curl -sSL "http://neuro.debian.net/lists/$( lsb_release -c | cut -f2 ).us-ca.full" >> /etc/apt/sources.list.d/neurodebian.sources.list && \
  sudo apt-key add $HOME/neurodebian.gpg && \
  (apt-key adv --refresh-keys --keyserver hkp://ha.pool.sks-keyservers.net 0xA5D32F012649A5A9 || true)

RUN sudo ln -fs /usr/share/zoneinfo/America/Montreal /etc/localtime && \
  sudo apt-get install -y --no-install-recommends tzdata && \
  sudo dpkg-reconfigure -f noninteractive tzdata && \
  sudo apt-get update && \
  sudo apt-get install -y --no-install-recommends fsl-core && \
  sudo apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Configure FSL environment
ENV FSLDIR="/usr/share/fsl/5.0/"
ENV FSL_DIR="${FSLDIR}" \
  FSLOUTPUTTYPE=NIFTI_GZ \
  PATH="/usr/share/fsl/5.0/bin:$PATH" \
  LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH

# Install minc-toolkit
RUN curl -L --output $HOME/minc-toolkit-1.9.18.deb http://packages.bic.mni.mcgill.ca/minc-toolkit/Debian/minc-toolkit-1.9.18-20200813-Ubuntu_18.04-x86_64.deb && \
  sudo dpkg -i $HOME/minc-toolkit-1.9.18.deb

# minc-toolkit configuration parameters for 1.9.18-20200813
ENV MINC_TOOLKIT=/opt/minc/1.9.18 \
  MINC_TOOLKIT_VERSION="1.9.18-20200813"
ENV PATH=${MINC_TOOLKIT}/bin:${MINC_TOOLKIT}/pipeline:${PATH} \
  PERL5LIB=${MINC_TOOLKIT}/perl:${MINC_TOOLKIT}/pipeline${PERL5LIB:+:$PERL5LIB} \
  LD_LIBRARY_PATH=${MINC_TOOLKIT}/lib:${MINC_TOOLKIT}/lib/InsightToolkit${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} \
  MNI_DATAPATH=${MINC_TOOLKIT}/../share:${MINC_TOOLKIT}/share \
  MINC_FORCE_V2=1 \
  MINC_COMPRESS=4 \
  VOLUME_CACHE_THRESHOLD=-1 \
  MANPATH=${MINC_TOOLKIT}/man${MANPATH:+:$MANPATH} \
  # integrated ANTs tools
  ANTSPATH=${MINC_TOOLKIT}/bin

USER rabies

# install miniconda
ENV MINICONDA_VERSION=4.8.2
ENV CONDA_DIR=${HOME}/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_${MINICONDA_VERSION}-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH

# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile

# make conda activate command available from /bin/bash --interative shells
RUN conda init bash

#### install RABIES
ENV RABIES_VERSION=0.2.1-dev
ENV RABIES=$HOME/RABIES-${RABIES_VERSION}
ENV PYTHONPATH="${PYTHONPATH}:$RABIES" \
  PATH=${PATH}:${RABIES}/minc-toolkit-extras:${RABIES}/twolevel_ants_dbm:${RABIES}/rabies/shell_scripts

# download code and create conda environment
RUN mkdir -p temp && \
  #curl -L --retry 5 -o temp/RABIES.tar.gz https://github.com/CoBrALab/RABIES/archive/${RABIES_VERSION}.tar.gz && \
  #cd temp && \
  #tar zxf RABIES.tar.gz && \
  #cd .. && \
  #conda env create -f temp/RABIES-${RABIES_VERSION}/rabies_environment.yml && \
  #bash temp/RABIES-${RABIES_VERSION}/install.sh && \
  rm -r temp && \
  git clone https://github.com/CoBrALab/RABIES && \
  mv RABIES $RABIES && \
  conda env create -f $RABIES/rabies_environment.yml && \
  conda run -n rabies /bin/bash $RABIES/install.sh

# overide the path with the rabies environment so that it becomes the default
ENV PATH=$CONDA_DIR/envs/rabies/bin:$PATH

WORKDIR /tmp/

ENTRYPOINT ["/home/rabies/RABIES-0.2.1-dev/bin/rabies"]
