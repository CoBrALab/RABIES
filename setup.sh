module load anaconda/5.1.0-python3 minc-toolkit/1.9.17 minc-stuffs/0.1.24^minc-toolkit-1.9.17 qbatch/git pydpiper/2.0.13
module load AFNI ANTs/20190211
source activate /home/cic/desgab/miniconda3/envs/custom_nipype
export PYTHONPATH="${PYTHONPATH}:/data/chamal/projects/Gabriel_DG/software/Mouse_fmriPype/"

mkdir -p $HOME/bin
echo -e '#! /usr/bin/env python \nfrom mfp.run_main import execute_workflow \nexecute_workflow()' > $HOME/bin/mfp
chmod +x $HOME/bin/mfp
export PATH=$PATH":$HOME/bin"
