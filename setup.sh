module load anaconda/5.0.1-python3 minc-toolkit/1.9.15 minc-stuffs/0.1.21^minc-toolkit-1.9.15 qbatch/git pydpiper/2.0.10
module load AFNI ANTs
source activate /home/cic/desgab/miniconda3/envs/custom_nipype
export PYTHONPATH="${PYTHONPATH}:$PWD"
