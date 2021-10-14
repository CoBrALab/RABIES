FROM ubuntu:18.04 as base
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    perl \
    imagemagick \
    parallel \
    locales \
    gdebi-core curl unzip \
    tcsh \
    xfonts-base \
    gsl-bin\
    netpbm \
    libjpeg62 \
    xvfb \
    libglw1-mesa \
    libxm4 \
    libgfortran4 \
    sudo \
    ca-certificates \
    rsync \
    gnupg software-properties-common \
  && rm -rf /var/lib/apt/lists/* \
  && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8

# add user to build all tools
RUN useradd -ms /bin/bash rabies && \
    echo "rabies ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/rabies && \
    chmod 0440 /etc/sudoers.d/rabies

#Fix silly AFNI not being properly built for modern ubuntu
RUN ln -sf /usr/lib/x86_64-linux-gnu/libgsl.so.23 /usr/lib/x86_64-linux-gnu/libgsl.so.19

####################################################################################################
FROM base as builder
RUN apt-get update && apt-get install -y gnupg software-properties-common --no-install-recommends \
    && curl -sSL https://apt.kitware.com/keys/kitware-archive-latest.asc | apt-key add - \
    && apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main' \
    && apt-get update && apt-get install -y --no-install-recommends \
      git cmake \
      build-essential automake libtool bison \
      libz-dev libjpeg-dev libpng-dev libtiff-dev \
      liblcms2-dev flex libx11-dev freeglut3-dev libxmu-dev \
      libxi-dev libqt4-dev libxml2-dev ninja-build \
    && rm -rf /var/lib/apt/lists/*


#Download and build ANTs
RUN mkdir -p /opt/ANTs/build && git clone https://github.com/ANTsX/ANTs.git /opt/ANTs/src \
    && cd /opt/ANTs/src \
    && git checkout 1759e5e23772e114a490cfa33a5764b400307b9d \
    && cd /opt/ANTs/build \
    && cmake -GNinja -DITK_BUILD_MINC_SUPPORT=ON -DBUILD_TESTING=OFF ../src \
    && cmake --build . \
    && cd ANTS-build \
    && cmake --install .


####################################################################################################
FROM base
#We only copy the ANTs commands we use, otherwise the container is huge
COPY --from=builder /opt/ANTs/bin/ /opt/quarantine/ANTs/bin/

#Install afni
RUN curl -L -O https://afni.nimh.nih.gov/pub/dist/bin/misc/@update.afni.binaries && \
    tcsh @update.afni.binaries -package linux_ubuntu_16_64 -apsearch yes -bindir /opt/quarantine/afni && \
    rm -f @update.afni.binaries

ENV PATH=/opt/quarantine/afni${PATH:+:$PATH}

#Install FSL
RUN curl -sSL https://raw.githubusercontent.com/nipy/nipype/master/docker/files/neurodebian.gpg | apt-key add - && \
    curl -sSL http://neuro.debian.net/lists/bionic.us-nh.full > /etc/apt/sources.list.d/neurodebian.sources.list && \
    apt-get update && apt-get install -y --no-install-recommends fsl-core && \
    rm -rf /var/lib/apt/lists/*

# Configure FSL environment
ENV FSLDIR="/usr/share/fsl/5.0/"
ENV FSL_DIR="${FSLDIR}" \
  FSLOUTPUTTYPE=NIFTI_GZ \
  PATH="/usr/share/fsl/5.0/bin:$PATH" \
  LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH
ENV PATH=/usr/lib/fsl/5.0:$PATH

#Install minc-toolkit
RUN curl -L --output /tmp/minc-toolkit-1.9.18.deb https://packages.bic.mni.mcgill.ca/minc-toolkit/min/minc-toolkit-1.9.18-20200813-Ubuntu_18.04-x86_64.deb && \
  gdebi -n /tmp/minc-toolkit-1.9.18.deb && \
  rm -f /tmp/minc-toolkit-1.9.18.deb

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
  MANPATH=${MINC_TOOLKIT}/man${MANPATH:+:$MANPATH}

# add a patch to nu_estimate_np_and_em
COPY patch/nu_estimate_np_and_em.diff nu_estimate_np_and_em.diff
RUN apt-get update -y && apt-get install -y --no-install-recommends patch
RUN (cd / && sudo patch -p0) < nu_estimate_np_and_em.diff && rm nu_estimate_np_and_em.diff

#Enable ANTs
ENV PATH=/opt/quarantine/ANTs/bin${PATH:+:$PATH} \
  ANTSPATH=/opt/quarantine/ANTs/bin
ENV PATH=/opt/ANTs/bin${PATH:+:$PATH}

USER rabies
WORKDIR /home/rabies
ENV HOME="/home/rabies"

#install conda
RUN curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash Miniforge3-Linux-x86_64.sh -b -p ${HOME}/miniforge && \
    rm -f Miniforge3-Linux-x86_64.sh
ENV CONDA_DIR=${HOME}/miniforge

#enable conda
RUN echo '. ${HOME}/miniforge/etc/profile.d/conda.sh' >> ${HOME}/.bashrc
RUN echo 'conda activate' >> ${HOME}/.bashrc

ENV CONDA_EXE='${CONDA_DIR}/bin/conda' \
  CONDA_PYTHON_EXE='${CONDA_DIR}/bin/python' \
  # override the path with the conda environment so that it becomes the default
  PATH=${CONDA_DIR}/bin:$PATH

# install RABIES
ENV RABIES=${HOME}/RABIES
RUN mkdir $RABIES

COPY rabies_environment.yml setup.py MANIFEST.in README.md LICENSE dependencies.txt $RABIES/

RUN . ${HOME}/miniforge/etc/profile.d/conda.sh && \
  conda activate && \
  conda env update -f $RABIES/rabies_environment.yml

COPY rabies $RABIES/rabies
COPY minc-toolkit-extras $RABIES/minc-toolkit-extras
COPY optimized_antsMultivariateTemplateConstruction $RABIES/optimized_antsMultivariateTemplateConstruction
COPY scripts $RABIES/scripts

RUN . ${HOME}/miniforge/etc/profile.d/conda.sh && \
  conda activate && \
  pip install -e $RABIES && \
  conda clean --all -y

RUN . ${HOME}/miniforge/etc/profile.d/conda.sh && conda activate && conda config --set auto_activate_base true

# adding 'agg' as default backend to avoid matplotlib errors
ENV MPLBACKEND agg

# pre-install the template defaults
ENV XDG_DATA_HOME=$HOME/.local/share

RUN install_DSURQE.sh $XDG_DATA_HOME/rabies

SHELL ["/bin/bash", "--login", "-c"]
ENTRYPOINT ["rabies"]
