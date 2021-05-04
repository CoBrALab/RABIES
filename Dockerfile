
FROM ubuntu:18.04 as base
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    perl \
    imagemagick \
    parallel \
    locales \
    python \
    gdebi-core curl unzip \
    tcsh \
    xfonts-base \
    gsl-bin\
    netpbm \
    libjpeg62 \
    xvfb \
    libglu1-mesa-dev \
    libglw1-mesa \
    libxm4 \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libgfortran4 \
    sudo \
    ca-certificates \
    rsync \
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
    && apt-get update && apt-get install -y \
      git cmake \
      build-essential automake libtool bison \
      libz-dev libjpeg-dev libpng-dev libtiff-dev \
      liblcms2-dev flex libx11-dev freeglut3-dev libxmu-dev \
      libxi-dev libqt4-dev libxml2-dev ninja-build  \
    && rm -rf /var/lib/apt/lists/*


#Download and build ANTs
RUN mkdir -p /opt/ANTs/build && git clone https://github.com/ANTsX/ANTs.git /opt/ANTs/src \
    && cd /opt/ANTs/src \
    && git checkout 5f1fab66da8ccfb30d242109aefb8df636698e9d \
    && cd /opt/ANTs/build \
    && cmake -GNinja -DITK_BUILD_MINC_SUPPORT=ON ../src \
    && cmake --build . \
    && cd ANTS-build \
    && cmake --install .


####################################################################################################
FROM base
#We only copy the ANTs commands we use, otherwise the container is huge
COPY --from=builder /opt/ANTs/bin/ /opt/ANTs/bin/

#Install afni
RUN curl -L -O https://afni.nimh.nih.gov/pub/dist/bin/misc/@update.afni.binaries && \
    tcsh @update.afni.binaries -package linux_ubuntu_16_64 -apsearch yes -bindir /opt/quarantine/afni && \
    rm -f @update.afni.binaries

RUN echo 'export PATH=/opt/quarantine/afni${PATH:+:$PATH}' > /etc/profile.d/99afni.sh

#Install FSL
RUN curl -L -O https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py && \
    python fslinstaller.py -d /opt/quarantine/fsl -p -q && \
    rm -f fslinstaller.py

RUN echo 'export FSLDIR=/opt/quarantine/fsl && export PATH=/opt/quarantine/fsl/bin${PATH:+:$PATH} && . ${FSLDIR}/etc/fslconf/fsl.sh' > /etc/profile.d/99fsl.sh

#Install minc-toolkit
RUN curl -L --output /tmp/minc-toolkit-1.9.18.deb http://packages.bic.mni.mcgill.ca/minc-toolkit/Debian/minc-toolkit-1.9.18-20200813-Ubuntu_18.04-x86_64.deb && \
  gdebi -n /tmp/minc-toolkit-1.9.18.deb && \
  rm -f /tmp/minc-toolkit-1.9.18.deb

COPY minc-toolkit-extras /opt/quarantine
COPY twolevel_ants_dbm /opt/quarantine

#Enable minc-toolkit
RUN echo '. /opt/minc/1.9.18/minc-toolkit-config.sh' > /etc/profile.d/98minc.sh
RUN echo 'export PATH=/opt/quarantine/minc-toolkit-extras${PATH:+:$PATH}' >> /etc/profile.d/98minc.sh
RUN echo 'export PATH=/opt/quarantine/twolevel_ants_dbm${PATH:+:$PATH}' >> /etc/profile.d/98minc.sh

#Enable ANTs
RUN echo 'export PATH=/opt/ANTs/bin${PATH:+:$PATH}' > /etc/profile.d/99ANTS.sh
RUN echo 'export ANTSPATH=/opt/ANTs/bin' >> /etc/profile.d/99ANTS.sh

#install conda
RUN curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash Miniforge3-Linux-x86_64.sh -b -p /opt/quarantine/miniforge && \
    rm -f Miniforge3-Linux-x86_64.sh

#enable conda
RUN cp /opt/quarantine/miniforge/etc/profile.d/conda.sh /etc/profile.d/99conda.sh
RUN echo 'conda activate' >> /etc/profile.d/99conda.sh
RUN echo 'export PATH=/opt/ANTs/bin${PATH:+:$PATH}' >> /etc/profile.d/99ANTs.sh

RUN . /etc/profile.d/99conda.sh && conda config --append channels simpleitk && \
  conda install -y 'networkx>=2.4' 'matplotlib>=3.1.1' 'nibabel>=2.3.1' 'nilearn>=0.4.2' 'nipype>=1.1.4' 'numpy>=1.16.2' 'pandas' 'scikit-learn>=0.20.0' 'scipy' 'simpleitk>=1.2.2' 'tqdm' 'pathos' && \
  conda activate && \
  pip install rabies==0.2.4

USER rabies
WORKDIR /home/rabies
ENV HOME="/home/rabies"

RUN . /etc/profile.d/99conda.sh && conda config --set auto_activate_base true

# pre-install the template defaults
ENV XDG_DATA_HOME=$HOME/.local/share

RUN . /etc/profile.d/99conda.sh && . /etc/profile.d/98minc.sh && install_DSURQE.sh $XDG_DATA_HOME/rabies

SHELL ["/bin/bash", "--login", "-c"]
ENTRYPOINT ["/bin/bash" , "--login", "-c", "rabies \"$@\"", "bash"]
