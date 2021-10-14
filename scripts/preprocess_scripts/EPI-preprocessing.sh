#!/bin/bash
#Preprocessing script for mouse brains, using Dorr-steadman template and MINC files converted from DICOM from Paravision 5.x

#If using conversion from PV6, might need to remove the volflip -y command


set -euo pipefail
set -x

calc(){ awk "BEGIN { print "$*" }"; }

tmpdir=$(mktemp -d)

input=$1
output=$2
modelfile=$3
modelmask=$4
reg_option=$5
otsu_thresh=$6

# convert inputs to mnc
nii2mnc ${input} ${tmpdir}/input.mnc
input=${tmpdir}/input.mnc
nii2mnc ${modelfile} ${tmpdir}/modelfile.mnc
modelfile=${tmpdir}/modelfile.mnc
nii2mnc ${modelmask} ${tmpdir}/modelmask.mnc
modelmask=${tmpdir}/modelmask.mnc

#Calculate Maximum FOV using the foreground/background of the fixed image
ThresholdImage 3 ${modelfile} ${tmpdir}/bgmask.h5 1e-12 Inf 1 0
ThresholdImage 3 ${modelfile} ${tmpdir}/otsu.h5 Otsu 4 ${tmpdir}/bgmask.h5
ThresholdImage 3 ${tmpdir}/otsu.h5 ${tmpdir}/otsu.h5 2 Inf 1 0
LabelGeometryMeasures 3 ${tmpdir}/otsu.h5 none ${tmpdir}/geometry.csv
fixed_maximum_resolution=$(python -c "print(max([ a*b for a,b in zip( [ a-b for a,b in zip( [float(x) for x in \"$(tail -1 ${tmpdir}/geometry.csv | cut -d, -f 14,16,18)\".split(\",\") ],[float(x) for x in \"$(tail -1 ${tmpdir}/geometry.csv | cut -d, -f 13,15,17)\".split(\",\") ])],[abs(x) for x in [float(x) for x in \"$(PrintHeader ${modelfile} 1)\".split(\"x\")]])]))")

origdistance=$fixed_maximum_resolution
distance=${origdistance}
levels=4
cycles=3
iters=100
lambda=2e-6
shrink=1.0
fwhm=0.1
stop=1e-5

do_correct() {
  distance=${origdistance}
  j=0
  while (( j < levels )); do
    i=0
    while (( i < cycles )); do
      nu_correct -clobber -normalize_field \
        -stop ${stop} -distance ${distance} -iterations ${iters} -fwhm ${fwhm} -shrink ${shrink} -lambda ${lambda} \
        -mask ${tmpdir}/weight.mnc ${n3input} ${tmpdir}/corrected_${distance}_${i}.mnc
      n3input=${tmpdir}/corrected_${distance}_${i}.mnc
      ((++i))
    done
    distance=$(calc  ${distance}/2)
    ((++j))
  done

}


minc_modify_header $input -sinsert :history=''

#Forceably convert to MINC2, and clamp range to avoid negative numbers, rescale to 0-65535
mincconvert -2 ${input} ${tmpdir}/originput.mnc
#Rescale initial data into entirely positive range (fix for completely negative data)
ImageMath 3 ${tmpdir}/originput.mnc RescaleImage ${tmpdir}/originput.mnc 0 65535


#Construct Otsu Mask of entire image
input=${tmpdir}/originput.mnc

ImageMath 3 ${tmpdir}/originput.mnc PadImage ${tmpdir}/originput.mnc 20

ThresholdImage 3 ${input} ${tmpdir}/bgmask.mnc 1e-12 Inf 1 0
ThresholdImage 3 ${input} ${tmpdir}/weight.mnc Otsu 4 ${tmpdir}/bgmask.mnc
ThresholdImage 3 ${tmpdir}/weight.mnc ${tmpdir}/weight.mnc $otsu_thresh Inf 1 0
iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 1 1 ball 1
ExtractRegionFromImageByMask 3 ${input} ${tmpdir}/repad.mnc ${tmpdir}/weight.mnc 1 10
ExtractRegionFromImageByMask 3 ${tmpdir}/originput.mnc ${tmpdir}/repad_orig.mnc ${tmpdir}/weight.mnc 1 10
mv -f ${tmpdir}/repad_orig.mnc ${tmpdir}/originput.mnc
mv -f ${tmpdir}/repad.mnc ${input}
ThresholdImage 3 ${input} ${tmpdir}/bgmask.mnc 1e-12 Inf 1 0
mincresample -like ${input} -keep -near -labels ${tmpdir}/weight.mnc ${tmpdir}/weight2.mnc
mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc

n3input=${tmpdir}/originput.mnc
do_correct

for file in ${tmpdir}/*imp; do
  echo nu_evaluate -clobber -mapping ${file} -mask ${tmpdir}/weight.mnc -field ${tmpdir}/$(basename $file .imp)_field.mnc ${input} ${tmpdir}/$(basename $file .imp).mnc
done | parallel

mincmath -clobber -mult ${tmpdir}/*field.mnc ${tmpdir}/precorrect_field_combined.mnc
mincmath -clobber -copy_header -zero -div ${tmpdir}/originput.mnc ${tmpdir}/precorrect_field_combined.mnc ${tmpdir}/precorrect.mnc

if [ $reg_option == 'Rigid' ]; then
  reg_type="--linear-type rigid --skip-nonlinear"
  antsRegistration_affine_SyN.sh $reg_type --convergence 1e-7 \
    ${tmpdir}/precorrect.mnc ${modelfile} ${tmpdir}/tomodel

  antsApplyTransforms -d 3 -i ${modelmask} -t [ ${tmpdir}/tomodel0_GenericAffine.xfm,1 ] \
    -o ${tmpdir}/newmask.mnc -n GenericLabel -r ${n3input} --verbose

elif [ $reg_option == 'Affine' ]; then
  reg_type="--linear-type affine --skip-nonlinear"
  antsRegistration_affine_SyN.sh $reg_type --convergence 1e-7 \
    ${tmpdir}/precorrect.mnc ${modelfile} ${tmpdir}/tomodel

  antsApplyTransforms -d 3 -i ${modelmask} -t [ ${tmpdir}/tomodel0_GenericAffine.xfm,1 ] \
    -o ${tmpdir}/newmask.mnc -n GenericLabel -r ${n3input} --verbose

elif [ $reg_option == 'SyN' ]; then
  reg_type="--linear-type affine"
  antsRegistration_affine_SyN.sh $reg_type --convergence 1e-7 \
    ${tmpdir}/precorrect.mnc ${modelfile} ${tmpdir}/tomodel

  antsApplyTransforms -d 3 -i ${modelmask} -t [ ${tmpdir}/tomodel0_GenericAffine.xfm,1 ] -t ${tmpdir}/tomodel1_inverse_NL.xfm \
    -o ${tmpdir}/newmask.mnc -n GenericLabel -r ${n3input} --verbose

fi

if [ $reg_option == 'no_reg' ]; then
  cp ${tmpdir}/precorrect.mnc ${tmpdir}/correct.mnc
else

  iMath 3 ${tmpdir}/newmask.mnc MD ${tmpdir}/newmask.mnc 1 1 ball 1
  ThresholdImage 3 ${n3input} ${tmpdir}/weight.mnc Otsu 4 ${tmpdir}/newmask.mnc
  ThresholdImage 3 ${tmpdir}/weight.mnc ${tmpdir}/weight.mnc 2 Inf 1 0
  ImageMath 3 ${tmpdir}/weight.mnc m ${tmpdir}/newmask.mnc ${tmpdir}/weight.mnc
  iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 1 1 ball 1
  ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
  iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 1 1 ball 1
  mincresample -like ${input} -keep -near -labels ${tmpdir}/weight.mnc ${tmpdir}/weight2.mnc
  mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc

  n3input=${tmpdir}/precorrect.mnc
  do_correct

  for file in ${tmpdir}/*imp; do
    echo nu_evaluate -clobber -mapping ${file} -mask ${tmpdir}/weight.mnc -field ${tmpdir}/$(basename $file .imp)_field.mnc ${input} ${tmpdir}/$(basename $file .imp).mnc
  done | parallel

  mincmath -clobber -mult ${tmpdir}/*field.mnc ${tmpdir}/precorrect_field_combined.mnc ${tmpdir}/field_final.mnc
  mincmath -clobber -copy_header -zero -div ${tmpdir}/originput.mnc ${tmpdir}/field_final.mnc ${tmpdir}/correct.mnc
fi

ExtractRegionFromImageByMask 3 ${tmpdir}/correct.mnc ${tmpdir}/recrop.mnc ${tmpdir}/weight.mnc 1 10
mv -f ${tmpdir}/recrop.mnc ${tmpdir}/correct.mnc
mincresample -keep -near -like ${tmpdir}/correct.mnc ${tmpdir}/weight.mnc ${tmpdir}/weight2.mnc
mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc

# reconvert to .nii
mnc2nii ${tmpdir}/correct.mnc ${tmpdir}/correct.nii
mnc2nii ${tmpdir}/weight.mnc ${tmpdir}/weight.nii
mnc2nii ${tmpdir}/precorrect.mnc ${tmpdir}/init_denoise.nii
gzip ${tmpdir}/*.nii
cp ${tmpdir}/correct.nii.gz $(dirname ${output})/$(basename ${output} .nii.gz).nii.gz
cp ${tmpdir}/weight.nii.gz $(dirname ${output})/$(basename ${output} .nii.gz)_mask.nii.gz
cp ${tmpdir}/init_denoise.nii.gz $(dirname ${output})/$(basename ${output} .nii.gz)_init_denoise.nii.gz

rm -rf ${tmpdir}
