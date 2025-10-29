FROM mambaorg/micromamba:jammy

USER root

# System-level dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    perl \
    imagemagick \
    parallel \
    gdebi-core \
    curl  \
    ca-certificates \
    patch \
    rsync \
    unzip \
  && rm -rf /var/lib/apt/lists/*

# Install a few minimal AFNI components
RUN curl -L -O https://afni.nimh.nih.gov/pub/dist/tgz/linux_ubuntu_16_64.tgz && \
    mkdir -p /opt/afni && \
    tar xzvf linux_ubuntu_16_64.tgz -C /opt/afni linux_ubuntu_16_64/{libmri.so,libf2c.so,3dDespike,3dTshift,3dWarp,3dAutobox} --strip-components=1 && \
    rm -f linux_ubuntu_16_64.tgz
ENV PATH=/opt/afni${PATH:+:$PATH}

# Install minc-toolkit
RUN curl -L --output /tmp/minc-toolkit-1.9.18.deb \
  https://packages.bic.mni.mcgill.ca/minc-toolkit/min/minc-toolkit-1.9.18-20200813-Ubuntu_18.04-x86_64.deb && \
  gdebi -n /tmp/minc-toolkit-1.9.18.deb && \
  rm -f /tmp/minc-toolkit-1.9.18.deb && \
  rm -f /opt/minc/1.9.18/bin/{ants*,ANTS*,ImageMath,AverageImages,ThresholdImage,ExtractRegionFromImageByMask,ConvertImage,AverageAffine*,ResampleImage}

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
RUN (cd / && patch -p0) < nu_estimate_np_and_em.diff && rm nu_estimate_np_and_em.diff

# ANTs install
RUN curl -L --output /tmp/ants.zip https://github.com/ANTsX/ANTs/releases/download/v2.5.0/ants-2.5.0-ubuntu-22.04-X64-gcc.zip && \
    unzip -d /opt /tmp/ants.zip && \
    rm -rf /opt/ants-2.5.0/lib && \
    rm -f /tmp/ants.zip

ENV PATH=/opt/ants-2.5.0/bin:${PATH}

USER $MAMBA_USER
ENV HOME=/home/$MAMBA_USER

# install RABIES
ENV RABIES=${HOME}/RABIES
RUN mkdir $RABIES

COPY rabies_environment.yml setup.py MANIFEST.in README.md LICENSE dependencies.txt $RABIES/

COPY rabies $RABIES/rabies
COPY minc-toolkit-extras $RABIES/minc-toolkit-extras
COPY optimized_antsMultivariateTemplateConstruction $RABIES/optimized_antsMultivariateTemplateConstruction
COPY scripts $RABIES/scripts

RUN micromamba install -y -n base -f $RABIES/rabies_environment.yml && \
    micromamba run -n base pip install -e $RABIES && \
    micromamba clean --all --yes

# FSL conda packages don't properly setup FSL, do it manually
ENV FSLDIR=/opt/conda
ENV FSLWISH=/opt/conda/bin/fslwish
ENV FSLTCLSH=/opt/conda/bin/fsltclsh
ENV FSLMULTIFILEQUIT=TRUE
ENV FSL_LOAD_NIFTI_EXTENSIONS=0
ENV FSLGECUDAQ=
ENV FSL_SKIP_GLOBAL=0
ENV FSLOUTPUTTYPE=NIFTI_GZ

# adding 'agg' as default backend to avoid matplotlib errors
ENV MPLBACKEND agg

# pre-install the template defaults
ENV XDG_DATA_HOME=${HOME}/.local/share

RUN micromamba run -n base install_DSURQE.sh $XDG_DATA_HOME/rabies

# Run a basic test
RUN micromamba run -n base error_check_rabies.py --complete

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "rabies"]
