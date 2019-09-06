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
RUN useradd -ms /bin/bash nistmni && \
    echo "nistmni ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/nistmni && \
    chmod 0440 /etc/sudoers.d/nistmni

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

ENV AFNIPATH="$HOME/ants-v2.3.1/bin" \
    PATH="$HOME/afni/linux_ubuntu_16_64:$PATH"


### Local CIC VM from https://github.com/CobraLab/MINC-VM/blob/master/provision.sh

RUN apt-get update && \
  apt-get install -y --no-install-recommends htop nano wget imagemagick parallel zram-config debconf

#Build tools and dependencies
RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections && \
  apt install -y --no-install-recommends build-essential gdebi-core \
    git imagemagick libssl-dev cmake autotools-dev automake \
    ed zlib1g-dev libxml2-dev libxslt-dev openjdk-8-jre \
    zenity libcurl4-openssl-dev bc gawk libxkbcommon-x11-0 \
    ttf-mscorefonts-installer bc

USER nistmni
ENV HOME /home/nistmni
WORKDIR /home/nistmni

ENV CONDA_DIR="$HOME/miniconda-latest" \
    PATH="$HOME/miniconda-latest/bin:$PATH" \
    ND_ENTRYPOINT="$HOME/startup.sh"

RUN export PATH="$HOME/miniconda-latest/bin:$PATH" \
    && echo "Downloading Miniconda installer ..." \
    && conda_installer="/tmp/miniconda.sh" \
    && curl -fsSL --retry 5 -o "$conda_installer" https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash "$conda_installer" -b -p $HOME/miniconda-latest \
    && rm -f "$conda_installer" \
    && conda update -yq -nbase conda \
    && conda config --append channels conda-forge \
    && conda config --append channels bioconda \
    && conda config --append channels simpleitk \
    && sync && conda clean --all && sync \
    && conda install -y -q python=3.6.8 numpy==1.16.2 pandas scipy python-graphviz pip cython setuptools simpleitk scikit-image scikit-learn==0.20.0 nibabel==2.3.1 nilearn==0.4.2 nipype==1.1.4 \
    && sync && conda clean --all && sync \
    && rm -rf ~/.cache/pip/* \
    && sync

USER root

# set paths
ENV minc_toolkit_v2=https://packages.bic.mni.mcgill.ca/minc-toolkit/Debian/minc-toolkit-1.9.17-20190313-Ubuntu_18.04-x86_64.deb \
  minc_toolkit_v1=https://packages.bic.mni.mcgill.ca/minc-toolkit/Debian/minc-toolkit-1.0.09-20170529-Ubuntu_18.04-x86_64.deb \
  bic_mni_models=http://packages.bic.mni.mcgill.ca/minc-toolkit/Debian/bic-mni-models-0.1.1-20120421.deb \
  beast_library=http://packages.bic.mni.mcgill.ca/minc-toolkit/Debian/beast-library-1.1.0-20121212.deb \
  pyminc=https://github.com/Mouse-Imaging-Centre/pyminc/archive/v0.52.tar.gz \
  minc_stuffs=https://github.com/Mouse-Imaging-Centre/minc-stuffs/archive/v0.1.24.tar.gz \
  pyezminc=https://github.com/BIC-MNI/pyezminc/archive/release-1.2.01.tar.gz \
  quarter=https://github.com/Alexpux/Quarter/archive/master.tar.gz \
  bicinventor=https://github.com/BIC-MNI/bicInventor/archive/master.tar.gz \
  brain_view2=https://github.com/Mouse-Imaging-Centre/brain-view2/archive/master.tar.gz \
  itksnap_minc=http://www.bic.mni.mcgill.ca/~vfonov/temp/itksnap-3.4.0-20151130-Linux-x86_64-qt4.tar.gz \
  mni_cortical_statistics=https://github.com/BIC-MNI/mni.cortical.statistics/archive/ver-0_9_5.tar.gz \
  generate_deformation_fields=https://github.com/Mouse-Imaging-Centre/generate_deformation_fields/archive/1.0.1.tar.gz \
  pydpiper=https://github.com/Mouse-Imaging-Centre/pydpiper/archive/v2.0.13.tar.gz \
  bpipe=https://github.com/ssadedin/bpipe/releases/download/0.9.9.6/bpipe-0.9.9.6.tar.gz


#Download and install external debs
RUN wget --progress=dot:mega $minc_toolkit_v2 && \
  wget --progress=dot:mega $minc_toolkit_v1 && \
  wget --progress=dot:mega $bic_mni_models && \
  wget --progress=dot:mega $beast_library && \
  for file in *.deb; do gdebi --n $file; done && \
  rm -f *.deb


RUN echo '. /opt/minc/1.9.17/minc-toolkit-config.sh' >> $HOME/.bashrc && \
  echo 'export LD_LIBRARY_PATH=/opt/minc/1.9.17/lib' >> $HOME/.bashrc && \
  echo 'export PATH=/opt/minc-toolkit-extras/:$PATH' >> $HOME/.bashrc

#Enable minc-toolkit in this script, need to escape error checking
RUN set +u && \
  . /opt/minc/1.9.17/minc-toolkit-config.sh && \
  set -u

#Download other packages
RUN wget --progress=dot:mega $pyminc -O pyminc.tar.gz && \
  wget --progress=dot:mega $minc_stuffs -O minc-stuffs.tar.gz && \
  wget --progress=dot:mega $pyezminc -O pyezminc.tar.gz && \
  wget --progress=dot:mega $generate_deformation_fields -O generate_deformation_fields.tar.gz && \
  wget --progress=dot:mega $pydpiper -O pydpiper.tar.gz && \
  wget --progress=dot:mega $bpipe -O bpipe.tar.gz && \
  wget https://raw.githubusercontent.com/andrewjanke/volgenmodel/master/volgenmodel -O /usr/local/bin/volgenmodel && \
  git clone https://github.com/CobraLab/minc-toolkit-extras.git /opt/minc-toolkit-extras


#Do this so that we don't need to keep track of version numbers for build
RUN mkdir -p pyminc && tar xzvf pyminc.tar.gz -C pyminc --strip-components 1 && \
  mkdir -p minc-stuffs && tar xzvf minc-stuffs.tar.gz -C minc-stuffs --strip-components 1 && \
  mkdir -p generate_deformation_fields && tar xzvf generate_deformation_fields.tar.gz -C generate_deformation_fields --strip-components 1 && \
  mkdir -p pyezminc && tar xzvf pyezminc.tar.gz -C pyezminc --strip-components 1 && \
  mkdir -p pydpiper && tar xzvf pydpiper.tar.gz -C pydpiper --strip-components 1 && \
  mkdir -p /opt/bpipe && tar xzvf bpipe.tar.gz -C /opt/bpipe --strip-components 1 && ln -s /opt/bpipe/bin/bpipe /usr/local/bin/bpipe


#Build and install packages
RUN ( cd pyezminc && python setup.py install --mincdir /opt/minc/1.9.17 ) && \
  ( cd pyminc && python setup.py install ) && \
  ( cd minc-stuffs && ./autogen.sh && ./configure --with-build-path=/opt/minc/1.9.17 && make && make install && python setup.py install ) && \
  ( cd generate_deformation_fields && ./autogen.sh && ./configure --with-minc2 --with-build-path=/opt/minc/1.9.17 && make && make install) && \
  ( cd generate_deformation_fields/scripts && python setup.py build_ext --inplace && python setup.py install) && \
  ( cd pydpiper && python setup.py install)

RUN pip install -U qbatch==2.1.5 && \
  rm -rf pyezminc* pyminc* minc-stuffs* generate_deformation_fields* pydpiper* bpipe*

#Installing brain-view2
RUN apt install -y --no-install-recommends libcoin80-dev libpcre++-dev qt4-default libqt4-opengl-dev libtool && \
  wget $quarter -O quarter.tar.gz && \
  wget $bicinventor -O bicinventor.tar.gz && \
  wget $brain_view2 -O brain-view2.tar.gz && \
  mkdir quarter && tar xzvf quarter.tar.gz -C quarter --strip-components 1 && \
  mkdir bicinventor && tar xzvf bicinventor.tar.gz -C bicinventor --strip-components 1 && \
  mkdir brain-view2 && tar xzvf brain-view2.tar.gz -C brain-view2 --strip-components 1 && \
  ( cd quarter && cmake . && make && make install ) && \
  ( cd bicinventor && ./autogen.sh && ./configure --with-build-path=/opt/minc/1.9.17 --prefix=/opt/minc/1.9.17 --with-minc2 && make && make install ) && \
  ( cd brain-view2 && /usr/bin/qmake-qt4 MINCDIR=/opt/minc/1.9.17 HDF5DIR=/opt/minc/1.9.17 INVENTORDIR=/opt/minc/1.9.17 && make && cp brain-view2 /opt/minc/1.9.17/bin ) && \
  rm -rf quarter* bicinventor* brain-view2*

#Install itksnap-MINC
RUN wget $itksnap_minc -O itksnap_minc.tar.gz && \
  tar xzvf itksnap_minc.tar.gz -C /usr/local --strip-components 1 && \
  rm -f itksnap_minc.tar.gz

ENV export MINC_PATH=/opt/minc/1.9.17 \
  export PATH=${OLDPATH}

#Purge unneeded packages
RUN apt-get purge $(dpkg -l | tr -s ' ' | cut -d" " -f2 | sed 's/:amd64//g' | grep -e -E '(-dev|-doc)$')

#Remove a hunk of useless packages which seem to be safe to remove
RUN apt-get -y purge printer-driver.* xserver-xorg-video.* xscreensaver.* wpasupplicant wireless-tools .*vdpau.* \
  bluez-cups cups-browsed cups-bsd cups-client cups-common cups-core-drivers cups-daemon cups-filters \
  cups-filters-core-drivers cups-ppdc cups-server-common linux-headers.* snapd bluez linux-firmware .*sane.* .*ppds.* && \
  apt-get -y clean && \
  apt-get -y --purge autoremove

#Cleanup to ensure extra files aren't packed into VM
RUN cd ~ && \
  rm -rf /tmp/provision && \
  rm -f /var/cache/apt/archives/*.deb && \
  rm -rf /var/lib/apt/lists/*


#### install Mouse_fmriPype

FROM test:latest

ENV export LD_LIBRARY_PATH=/opt/minc/1.9.17/lib \
  export PATH=/opt/minc-toolkit-extras/:$PATH

RUN . /opt/minc/1.9.17/minc-toolkit-config.sh && \
  cd $HOME && \
  git clone https://github.com/Gab-D-G/Mouse_fmriPype && \
  bash $HOME/Mouse_fmriPype/setup.sh && \
  git clone https://github.com/CobraLab/twolevel_ants_dbm && \
  echo 'export PATH=$HOME/twolevel_ants_dbm/twolevel_dbm.py:$PATH' >> $HOME/.bashrc

#write container execution script
RUN echo "#! /usr/bin/env python" > /home/nistmni/Mouse_fmriPype/exec.py && \
  echo "import os" >> /home/nistmni/Mouse_fmriPype/exec.py && \
  echo "import sys" >> /home/nistmni/Mouse_fmriPype/exec.py && \
  echo "os.environ['MFP'] = '/home/nistmni/Mouse_fmriPype'" >> /home/nistmni/Mouse_fmriPype/exec.py && \
  echo "sys.path.insert(0,os.environ['MFP'])" >> /home/nistmni/Mouse_fmriPype/exec.py && \
  echo "from mfp.run_main import execute_workflow" >> /home/nistmni/Mouse_fmriPype/exec.py && \
  echo "execute_workflow()" >> /home/nistmni/Mouse_fmriPype/exec.py && \
  chmod +x /home/nistmni/Mouse_fmriPype/exec.py

ENV QBATCH_SYSTEM local

USER nistmni
ENV HOME /home/nistmni
WORKDIR /home/nistmni

ENTRYPOINT ["/home/nistmni/Mouse_fmriPype/exec.py"]
