#! /home/rabies/miniconda-latest/envs/rabies/bin/python
import os
import sys

# minc-toolkit configuration parameters for 1.9.18-20200813
os.environ['PATH'] = '/opt/minc/1.9.18/bin:%s' % (os.environ['PATH'])
os.environ['LD_LIBRARY_PATH'] = '/opt/minc/1.9.18/lib:%s' % (os.environ['LD_LIBRARY_PATH'])

os.environ['MINC_TOOLKIT']='/opt/minc/1.9.18'
os.environ['MINC_TOOLKIT_VERSION']="1.9.18-20200813"
os.environ['PATH']="%s/bin:%s/pipeline:%s" % (os.environ['MINC_TOOLKIT'],os.environ['MINC_TOOLKIT'],os.environ['PATH'])
os.environ['PERL5LIB']="%s/perl:%s/pipeline${PERL5LIB:+:$PERL5LIB}" % (os.environ['MINC_TOOLKIT'],os.environ['MINC_TOOLKIT'],)
os.environ['LD_LIBRARY_PATH']="%s/lib:%s/lib/InsightToolkit${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" % (os.environ['MINC_TOOLKIT'],os.environ['MINC_TOOLKIT'],)
os.environ['MNI_DATAPATH']="%s/../share:%s/share" % (os.environ['MINC_TOOLKIT'],os.environ['MINC_TOOLKIT'],)
os.environ['MINC_FORCE_V2']="1"
os.environ['MINC_COMPRESS']="4"
os.environ['VOLUME_CACHE_THRESHOLD']="-1"
os.environ['MANPATH']="%s/man${MANPATH:+:$MANPATH}" % (os.environ['MINC_TOOLKIT'],)
# integrated ANTs tools
os.environ['ANTSPATH']="%s/bin" % (os.environ['MINC_TOOLKIT'],)


# FSL
os.environ['FSLDIR'] = '/usr/share/fsl/5.0/'
os.environ['FSL_DIR'] = os.environ['FSLDIR']
os.environ['FSLOUTPUTTYPE'] = 'NIFTI_GZ'
os.environ['PATH'] = '/usr/share/fsl/5.0/bin:%s' % (os.environ['PATH'])
os.environ['LD_LIBRARY_PATH'] = '/usr/lib/fsl/5.0:%s' % (os.environ['LD_LIBRARY_PATH'])
# RABIES
os.environ['RABIES_VERSION'] = '0.2.1-dev'
os.environ['RABIES'] = "%s/RABIES-%s" % (os.environ['HOME'],os.environ['RABIES_VERSION'])
if 'PYTHONPATH' in os.environ.keys():
    os.environ['PYTHONPATH'] = '%s:%s' % (os.environ['PYTHONPATH'],os.environ['RABIES'])
else:
    os.environ['PYTHONPATH'] = os.environ['RABIES']
sys.path.insert(0,os.environ['RABIES'])
os.environ['PATH'] = '%s/rabies/shell_scripts:%s/twolevel_ants_dbm:%s/minc-toolkit-extras:/home/rabies/miniconda-latest/envs/rabies/bin:%s' % (os.environ['RABIES'],os.environ['RABIES'],os.environ['RABIES'],os.environ['PATH'])
from rabies.run_main import execute_workflow
execute_workflow()
