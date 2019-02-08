mask_path=$1
output=$2

mkdir -p tmp
nii2mnc $mask_path tmp/mask.mnc

mincmorph -successive E tmp/mask.mnc tmp/eroded_mask.mnc

cp tmp/eroded_mask.mnc $output

rm -r tmp
