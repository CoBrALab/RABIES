#!/home/murosevic/miniconda3/envs/rabies/bin/python
from argparse import ArgumentParser
import os
import numpy as np
import SimpleITK as sitk


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-o", "--output", type=str,
                        help="""
                        Name of output average file.
                        """)
    parser.add_argument('--file_list', type=str,
                        nargs="*",  # 0 or more values expected => creates a list
                        help="""
                        Specify a list of input files, space-seperated (i.e. file1 file2 ...).
                        """)
    parser.add_argument("--image_type", default='image',
                        choices=['image', 'affine', 'warp'],
                        help="""
                        Specify whether the type of image is a nifti structural image,
                        or a set of affine or non-linear (warp) transforms.
                        """)
    parser.add_argument("--method", default='trimmed_mean',
                        choices=['mean', 'median', 'trimmed_mean', 'huber'],
                        help="""
                        Specify of average method to create from the image list.
                        """)
    parser.add_argument("--trim_percent", type=float, default=15,
                        help="""
                        Specify % to trim off if using trimmed_mean.
                        """)
    parser.add_argument("--winsorize_lower_bound", type=float, default=0.005,
                        help="""
                        Specify % of lowest outlier intensities values to not consider during registration.
                        """)
    parser.add_argument("--winsorize_upper_bound", type=float, default=0.995,
                        help="""
                        Specify % of highest outlier intensities values to not consider during registration.
                        """)
    parser.add_argument("--normalize", dest='normalize', action='store_true',
                        help="""
                        Whether to divide each image by its mean before computing average.
                        """)
    opts = parser.parse_args()
    output = os.path.abspath(opts.output)

    if len(opts.file_list)==1:
        print("ONLY ONE INPUT PROVIDED TO --file_list. THE OUTPUT IS THE INPUT.")
        sitk.WriteImage(sitk.ReadImage(opts.file_list[0]), output)
        quit()

    # takes a average out of the array values from a list of Niftis
    array_list = []
    for file in opts.file_list:
        if os.path.isfile(file) and (opts.image_type=='warp' or opts.image_type=='image'):
            array = sitk.GetArrayFromImage(sitk.ReadImage(file))
        elif os.path.isfile(file) and opts.image_type=='affine':
            transform = sitk.ReadTransform(file)
            array = np.array(transform.GetParameters())
        else:
            continue
        shape = array.shape # we assume all inputs have the same shape
        array = array.flatten()
        if opts.normalize: # divide the image values by its mean
            array /= array.mean()
        array_list.append(array)

    concat_array = np.array(array_list)
    if opts.method == 'mean':
        average = np.mean(concat_array, axis=0)
    elif opts.method == 'median':
        average = np.median(concat_array,axis=0)
    elif opts.method == 'trimmed_mean':
        from scipy import stats
        average = stats.trim_mean(concat_array, opts.trim_percent/100, axis=0)
    elif opts.method == 'huber':
        import statsmodels.api as sm
        average = sm.robust.scale.huber(concat_array)

    average = average.reshape(shape)

    if opts.image_type=='image':
        average_img = sitk.GetImageFromArray(average, isVector=False)
        average_img.CopyInformation(sitk.ReadImage(opts.file_list[0]))
        sitk.WriteImage(average_img, output)
    elif opts.image_type=='warp':
        average_img = sitk.GetImageFromArray(average, isVector=True)
        average_img.CopyInformation(sitk.ReadImage(opts.file_list[0]))
        sitk.WriteImage(average_img, output)
    elif opts.image_type=='affine':
        transform.SetParameters(average)
        sitk.WriteTransform(transform, output)
