mkdir -p anat_preproc
mkdir -p flip

for file in ../minc/*anat.mnc; do var=$(basename $file); echo volflip $file flip/${var:0:8}flipped.mnc; done > parallel_flip.sh
bash parallel_flip.sh
for file in flip/*flipped.mnc; do var=$(basename $file); echo /data/chamal/projects/gabriel/src/minc-toolkit-extras/mouse-preprocessing-v3.sh $file anat_preproc/${var:0:8}preproc_anat.mnc; done > parallel_preproc.sh
bash parallel_preproc.sh

#qbatch --ppj 4 parallel_flip.sh
