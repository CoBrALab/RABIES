mkdir -p nii_files

for file in realigned_files/*.mnc; do mnc2nii $file nii_files/$(basename $file ${file:(-4)}).nii; done

gzip nii_files/*.nii
