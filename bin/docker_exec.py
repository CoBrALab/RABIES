#! /home/rabies/miniconda-latest/envs/rabies/bin/python
import os
import sys

# minc-toolkit configuration parameters for 1.9.18-20200813
os.environ['PATH'] = '/opt/minc/1.9.18/bin:${PATH}'
os.environ['LD_LIBRARY_PATH'] = '/opt/minc/1.9.18/lib:${LD_LIBRARY_PATH}'

os.environ['MINC_TOOLKIT']='/opt/minc/1.9.18'
os.environ['MINC_TOOLKIT_VERSION']="1.9.18-20200813"
os.environ['PATH']="${MINC_TOOLKIT}/bin:${MINC_TOOLKIT}/pipeline:${PATH}"
os.environ['PERL5LIB']="${MINC_TOOLKIT}/perl:${MINC_TOOLKIT}/pipeline${PERL5LIB:+:$PERL5LIB}"
os.environ['LD_LIBRARY_PATH']="${MINC_TOOLKIT}/lib:${MINC_TOOLKIT}/lib/InsightToolkit${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
os.environ['MNI_DATAPATH']="${MINC_TOOLKIT}/../share:${MINC_TOOLKIT}/share"
os.environ['MINC_FORCE_V2']="1"
os.environ['MINC_COMPRESS']="4"
os.environ['VOLUME_CACHE_THRESHOLD']="-1"
os.environ['MANPATH']="${MINC_TOOLKIT}/man${MANPATH:+:$MANPATH}"
# integrated ANTs tools
os.environ['ANTSPATH']="${MINC_TOOLKIT}/bin"


# FSL
os.environ['FSLDIR'] = '/usr/share/fsl/5.0/'
os.environ['FSL_DIR'] = '${FSLDIR}'
os.environ['FSLOUTPUTTYPE'] = 'NIFTI_GZ'
os.environ['PATH'] = '/usr/share/fsl/5.0/bin:${PATH}'
os.environ['LD_LIBRARY_PATH'] = '/usr/lib/fsl/5.0:${LD_LIBRARY_PATH}'
# RABIES
os.environ['RABIES_VERSION'] = '0.2.1-dev'
os.environ['RABIES'] = '${HOME}/RABIES-${RABIES_VERSION}'
os.environ['PYTHONPATH'] = '${PYTHONPATH}:${RABIES}'
sys.path.insert(0,os.environ['RABIES'])
os.environ['PATH'] = '${RABIES}/rabies/shell_scripts:${RABIES}/twolevel_ants_dbm:${RABIES}/minc-toolkit-extras:/home/rabies/miniconda-latest/envs/rabies/bin:${PATH}'
from rabies.run_main import execute_workflow
execute_workflow()
