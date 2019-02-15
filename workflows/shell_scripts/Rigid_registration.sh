EPI=$1
anat_file=$2
mask=$3
filename_template=$4

antsRegistration -d 3 \
--verbose -o [${filename_template}_output_,${filename_template}_output_warped_image.nii.gz] \
-t Rigid[0.1] -m Mattes[$EPI,$anat_file,1,64,None] \
-c 1000x500x250x100x50x25 -s 8x4x2x1x0.5x0 -f 6x5x4x3x2x1 --masks [NULL,NULL] \
--interpolation BSpline[5] -z 1 -u 0 -a 1
