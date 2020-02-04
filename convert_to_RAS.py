import sys
import os

in_file=os.path.abspath(sys.argv[1])

def convert_to_RAS(img_file, out_dir):
    #convert the input image to the RAS orientation convention
    import nibabel as nb
    img=nb.load(img_file)
    if not nb.aff2axcodes(img.affine)==('R','A','S'):
        nb.as_closest_canonical(img).to_filename(img_file)

convert_to_RAS(in_file, os.path.dirname(in_file))
