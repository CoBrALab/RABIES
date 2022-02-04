#! /usr/bin/env python
import sys
import SimpleITK as sitk
from rabies.utils import run_command

def otsu_bias_cor(target, otsu_ref, out_dir, b_value, mask=None, n_iter=200):
    command = f'ImageMath 3 {out_dir}/null_mask.mnc ThresholdAtMean {otsu_ref} 0'
    rc = run_command(command)
    command = f'ThresholdImage 3 {otsu_ref} {out_dir}/otsu_weight.mnc Otsu 4'
    rc = run_command(command)

    otsu_img = sitk.ReadImage(
        f'{out_dir}/otsu_weight.mnc', sitk.sitkUInt8)
    otsu_array = sitk.GetArrayFromImage(otsu_img)

    if mask is not None:
        resampled_mask_img = sitk.ReadImage(
            mask, sitk.sitkUInt8)
        resampled_mask_array = sitk.GetArrayFromImage(resampled_mask_img)

        otsu_array = otsu_array*resampled_mask_array

    combined_mask=(otsu_array==1.0)+(otsu_array==2.0)
    mask_img=sitk.GetImageFromArray(combined_mask.astype('uint8'), isVector=False)
    mask_img.CopyInformation(otsu_img)
    sitk.WriteImage(mask_img, f'{out_dir}/mask12.mnc')

    combined_mask=(otsu_array==1.0)+(otsu_array==2.0)+(otsu_array==3.0)
    mask_img=sitk.GetImageFromArray(combined_mask.astype('uint8'), isVector=False)
    mask_img.CopyInformation(otsu_img)
    sitk.WriteImage(mask_img, f'{out_dir}/mask123.mnc')

    combined_mask=(otsu_array==1.0)+(otsu_array==2.0)+(otsu_array==3.0)+(otsu_array==4.0)
    mask_img=sitk.GetImageFromArray(combined_mask.astype('uint8'), isVector=False)
    mask_img.CopyInformation(otsu_img)
    sitk.WriteImage(mask_img, f'{out_dir}/mask1234.mnc')

    extend_mask12(target,out_dir, b_value=b_value, mask=mask, n_iter=200)

    command = f'N4BiasFieldCorrection -v -d 3 -i {target} -b {str(b_value)} -s 4 -c [{str(n_iter)}x{str(n_iter)}x{str(n_iter)},1e-4] -w {out_dir}/mask12.mnc -x {out_dir}/null_mask.mnc -o {out_dir}/corrected1_.mnc'
    rc = run_command(command)

    #command = 'N4BiasFieldCorrection -v -d 3 -i corrected1_.mnc -b %s -s 4 -c [%sx%sx%s,1e-4] -w mask123.mnc -x null_mask.mnc -o corrected2_.mnc' % (str(b_value), str(n_iter),str(n_iter),str(n_iter),)
    #rc = run_command(command)

    command = f'N4BiasFieldCorrection -v -d 3 -i {out_dir}/corrected1_.mnc -b {str(b_value)} -s 4 -c [{str(n_iter)}x{str(n_iter)}x{str(n_iter)},1e-4] -w {out_dir}/mask1234.mnc -x {out_dir}/null_mask.mnc -o {out_dir}/multistage_corrected.mnc'
    rc = run_command(command)

def extend_mask12(target, out_dir, b_value, mask=None, n_iter=200):

    command = f'N4BiasFieldCorrection -d 3 -i {target} -b {str(b_value)} -s 4 -c [{str(n_iter)}x{str(n_iter)}x{str(n_iter)},1e-4] -w {out_dir}/mask12.mnc -x {out_dir}/null_mask.mnc -o {out_dir}/corrected_lower.mnc'
    rc = run_command(command)

    command = f'ThresholdImage 3 {out_dir}/corrected_lower.mnc {out_dir}/otsu_weight.mnc Otsu 4'
    rc = run_command(command)

    otsu_img = sitk.ReadImage(
        f'{out_dir}/otsu_weight.mnc', sitk.sitkUInt8)
    otsu_array = sitk.GetArrayFromImage(otsu_img)

    if mask is not None:
        resampled_mask_img = sitk.ReadImage(
            mask, sitk.sitkUInt8)
        resampled_mask_array = sitk.GetArrayFromImage(resampled_mask_img)

        otsu_array = otsu_array*resampled_mask_array

    lower_mask=(otsu_array==1.0)+(otsu_array==2.0)

    def extend_mask(mask, lower_mask):
        otsu_img = sitk.ReadImage(
            mask, sitk.sitkUInt8)
        otsu_array = sitk.GetArrayFromImage(otsu_img)
        combined_mask = (otsu_array+lower_mask)>0

        mask_img=sitk.GetImageFromArray(combined_mask.astype('uint8'), isVector=False)
        mask_img.CopyInformation(otsu_img)
        sitk.WriteImage(mask_img, mask)

    extend_mask(f'{out_dir}/mask12.mnc', lower_mask)
    extend_mask(f'{out_dir}/mask123.mnc', lower_mask)
    extend_mask(f'{out_dir}/mask1234.mnc', lower_mask)


target=sys.argv[1]
out_dir=sys.argv[2]
b_value=sys.argv[3]
otsu_ref=target
otsu_bias_cor(target=target, otsu_ref=otsu_ref, out_dir=out_dir, b_value=b_value, n_iter=100)
