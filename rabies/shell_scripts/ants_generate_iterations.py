#!/usr/bin/env python

# This file generates steps of registration between two images and attempts to compensate
# For ANTs' dependency on the resolution of the file

# We do this by defining two scales to step over
# blur_scale, which is the real-space steps in blurring we will do
# shrink_scale, which is the subsampling scale that is 1/2 the fwhm blur scale, adjusted for file minimum resolution and max size

from __future__ import division, print_function

import argparse
import numpy as np
import sys


def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


parser.add_argument(
    '--min', help='minimum resolution of fixed file (mm)', type=float, required=True)
parser.add_argument(
    '--max', help='max size of fixed file (mm)', type=float, required=True)
parser.add_argument(
    '--start-scale', help='set starting scale (mm), default calculated from max size', type=float)
parser.add_argument(
    '--final-iterations', help='total number of iterations at lowest scale', type=int, default=25)
parser.add_argument(
    '--output', help='type of output to generate', default='generic', choices=['generic', 'affine', 'modelbuild', 'twolevel_dbm', 'multilevel-halving'])
parser.add_argument('--step-size', help='step mode for generation', default=1)
parser.add_argument(
    '--convergence', help='set convergence for generated stages', default='1e-6')

args = parser.parse_args()

# Setup inital inputs
min_resolution = args.min
max_size = args.max

if RepresentsInt(args.step_size):
    step_size = int(args.step_size)
elif args.step_size == "power2":
    step_size = args.step_size
else:
    sys.exit("Unrecognized step size")

# Make empty arrays
shrinks = []
blurs = []
iterations = []
bins = []

if args.output == "affine" or args.output == "multilevel-halving" and args.final_iterations == 25:
  min_iterations = 50
else:
  min_iterations = args.final_iterations

# Converter
fwhm_to_sigma = 2 * np.sqrt(2 * np.log(2))

# Inital resolution scaling
if args.start_scale:
  start_shrink = args.start_scale / min_resolution
else:
  start_shrink = max_size / 28 / min_resolution

if isinstance(step_size, int):
    for shrink_scale in range(int(np.around(start_shrink)), 0, -1 * step_size):
        shrinks.append(
            str(int(min(max_size / 32 / min_resolution, max(1.0, np.around(shrink_scale))))))
        blurs.append(str(shrink_scale * 2 * min_resolution / fwhm_to_sigma))
        iterations.append(str(min(2025, int(min_iterations * 3**(max(0,shrink_scale - 1))))))
        bins.append(
            str(int(np.around((max(32, 256 / max(1, int(shrinks[-1]))))))))
else:
    blur_scale = start_shrink * 2 * min_resolution
    shrink_scale = start_shrink
    while (blur_scale > 0.5 * min_resolution):
        shrinks.append(
            str(int(min(max_size / 32 / min_resolution, max(1.0, np.around(shrink_scale))))))
        blurs.append(str(blur_scale / fwhm_to_sigma))
        iterations.append(str(min(2025, int(min_iterations * 3**(max(0,shrink_scale-1))))))
        bins.append(
            str(int(np.around((max(32, 256 / max(1,  int(shrinks[-1]))))))))
        blur_scale = blur_scale / 2
        shrink_scale = shrink_scale / 2

shrinks.append("1")
blurs.append("0")
iterations.append(str(min_iterations))
bins.append(str(int(np.around((max(32, 256 / max(1, 1 * min_resolution)))))))

if args.output == 'affine':
    transforms = ["--transform Translation[ 0.1 ]",
                  "--transform Rigid[ 0.1 ]",
                  "--transform Similarity[ 0.1 ]",
                  "--transform Affine[ 0.1 ]"]
    masks = ["--masks [ NOMASK,NOMASK ]",
             "--masks [ NOMASK,NOMASK ]",
             "--masks [ NOMASK,NOMASK ]",
             "--masks [ ${fixedmask},${movingmask} ]" ]
    repeatmask = [ False,
                   False,
                   "--masks [ ${fixedmask},${movingmask} ]",
                   False ] 

    slicestart = 0
    slicesize = max(3,int(np.around(len(shrinks) / (len(transforms)) * 1.5)))
    
    for i, transform in enumerate(transforms):
        if i == len(transforms)-1:
          print(transform, end=' \\\n')
          print("\t--metric Mattes[ ${{fixedfile}},${{movingfile}},1,{},None ]".format(bins[slicestart]), end=' \\\n')
          print("\t--convergence [ {},{},10 ]".format("x".join(iterations[slicestart:]), args.convergence), end=' \\\n')
          print("\t--shrink-factors {}".format("x".join(shrinks[slicestart:])), end=' \\\n')
          print("\t--smoothing-sigmas {}mm".format("x".join(blurs[slicestart:])), end=' \\\n')
          print("\t" + masks[i], end= ' ')
          slicestart += int(np.floor(slicesize / 2))
        else:
          print(transform, end=' \\\n')
          print("\t--metric Mattes[ ${{fixedfile}},${{movingfile}},1,{},None ]".format(bins[slicestart]), end=' \\\n')
          print("\t--convergence [ {},{},10 ]".format("x".join(iterations[slicestart:slicestart + slicesize]), args.convergence), end=' \\\n')
          print("\t--shrink-factors {}".format("x".join(shrinks[slicestart:slicestart + slicesize])), end=' \\\n')
          print("\t--smoothing-sigmas {}mm".format("x".join(blurs[slicestart:slicestart + slicesize])), end=' \\\n')
          print("\t" + masks[i], end=' \\\n')
          if repeatmask[i]:
            print(transform, end=' \\\n')
            print("\t--metric Mattes[ ${{fixedfile}},${{movingfile}},1,{},None ]".format(bins[slicestart]), end=' \\\n')
            print("\t--convergence [ {},{},10 ]".format("x".join(iterations[slicestart:slicestart + slicesize]), args.convergence), end=' \\\n')
            print("\t--shrink-factors {}".format("x".join(shrinks[slicestart:slicestart + slicesize])), end=' \\\n')
            print("\t--smoothing-sigmas {}mm".format("x".join(blurs[slicestart:slicestart + slicesize])), end=' \\\n')
            print("\t" + repeatmask[i], end=' \\\n')
          slicestart += int(np.floor(slicesize / 2))
          
elif args.output == 'twolevel_dbm':
    print("--reg-iterations {}".format("x".join(iterations)), end=' \\\n')
    print("--reg-shrinks {}".format("x".join(shrinks)), end=' \\\n')
    print("--reg-smoothing {}mm".format("x".join(blurs)), end=' ')

elif args.output == 'modelbuild':
    print("-q {}".format("x".join(iterations)), end=' \\\n')
    print("-f {}".format("x".join(shrinks)), end=' \\\n')
    print("-s {}mm".format("x".join(blurs)), end=' ')

elif args.output == 'generic':
    print("--convergence [ {},{},10 ]".format("x".join(iterations), args.convergence), end=' \\\n')
    print("--shrink-factors {}".format("x".join(shrinks)), end=' \\\n')
    print("--smoothing-sigmas {}mm".format("x".join(blurs)), end=' ')

elif args.output == "multilevel-halving":
    transforms = ["--transform Translation[ 0.1 ]",
                  "--transform Rigid[ 0.1 ]",
                  "--transform Similarity[ 0.1 ]",
                  "--transform Affine[ 0.1 ]"]
    masks = ["--masks [ NOMASK,NOMASK ]",
             "--masks [ NOMASK,NOMASK ]",
             "--masks [ NOMASK,NOMASK ]",
             "--masks [ ${fixedmask},${movingmask} ]" ]
    repeatmask = [ False,
                   False,
                   "--masks [ ${fixedmask},${movingmask} ]",
                   False ] 

    slicestart=0
    for i, transform in enumerate(transforms):
      if i == len(transforms) - 1:
        print(transform, end=' \\\n')
        print("\t--metric Mattes[ ${{fixedfile}},${{movingfile}},1,{},None ]".format(bins[slicestart]), end=' \\\n')
        print("\t--convergence [ {},{},10 ]".format("x".join(iterations[slicestart:]), args.convergence), end=' \\\n')
        print("\t--shrink-factors {}".format("x".join(shrinks[slicestart:])), end=' \\\n')
        print("\t--smoothing-sigmas {}mm".format("x".join(blurs[slicestart:])), end=' \\\n')
        print("\t" + masks[i], end=' ')
      else:
        print(transform, end=' \\\n')
        print("\t--metric Mattes[ ${{fixedfile}},${{movingfile}},1,{},None ]".format(bins[slicestart]), end=' \\\n')
        print("\t--convergence [ {},{},10 ]".format("x".join(iterations[slicestart:len(blurs)-(4-i)]), args.convergence), end=' \\\n')
        print("\t--shrink-factors {}".format("x".join(shrinks[slicestart:len(blurs)-(4-i)])), end=' \\\n')
        print("\t--smoothing-sigmas {}mm".format("x".join(blurs[slicestart:len(blurs)-(4-i)])), end=' \\\n')
        print("\t" + masks[i], end=' \\\n')
        if repeatmask[i]:
          print(transform, end=' \\\n')
          print("\t--metric Mattes[ ${{fixedfile}},${{movingfile}},1,{},None ]".format(bins[slicestart]), end=' \\\n')
          print("\t--convergence [ {},{},10 ]".format("x".join(iterations[slicestart:len(blurs)-(4-i)]), args.convergence), end=' \\\n')
          print("\t--shrink-factors {}".format("x".join(shrinks[slicestart:len(blurs)-(4-i)])), end=' \\\n')
          print("\t--smoothing-sigmas {}mm".format("x".join(blurs[slicestart:len(blurs)-(4-i)])), end=' \\\n')
          print("\t" + repeatmask[i], end=' \\\n')
      slicestart = slicestart + int(len(blurs)/2**(i+1) - 1)
