#!/home/murosevic/miniconda3/envs/rabies/bin/python

# This file generates steps of registration between two images and attempts to compensate
# For ANTs' dependency on the resolution of the file

# We do this by defining two scales to step over
# blur_scale, which is the real-space steps in blurring we will do
# shrink_scale, which is the subsampling scale that is 1/2 the fwhm blur scale, adjusted for file minimum resolution and max size

from __future__ import division, print_function

import argparse
import sys
import re


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

def check_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

class SplitArgsComma(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.rstrip(',').split(','))

parser.add_argument(
    '--min', help='minimum resolution of fixed image (mm)', type=float, required=True)
parser.add_argument(
    '--max', help='maximum dimension of fixed image features (mm)', type=float, required=True)
parser.add_argument(
    '--start-scale', help='set starting sigma blur scale (mm), default calculated from max size', type=float)
parser.add_argument(
    '--final-iterations', help='total number of iterations at lowest scale', type=int, default=20)
parser.add_argument(
    '--output', help='type of output to generate', default='generic',
      choices=['generic', 'affine', 'modelbuild', 'twolevel_dbm', 'multilevel-halving', 'exhaustive-affine',
                'lsq6', 'lsq9', 'lsq12', 'rigid', 'similarity','volgenmodel','affine-plain'])
parser.add_argument('--volgen-iteration', help='for volgenmodel mode, control which iteration is output', default=0, type=int)
parser.add_argument('--step-size', help='step size for fwhm scale space, default 1/2 of voxel size', type=float)
parser.add_argument(
    '--convergence', help='set convergence for generated stages', default='1e-6')
parser.add_argument(
    '--convergence-window', help='set convergence window for generated stages', default='10')
parser.add_argument(
    '--close', help='images are already close, skip large scales of pyramid and translation and rigid steps for affine', action='store_true')
parser.add_argument(
    '--rough', help='skip full-resolution iterations for a "rough" alignment', action='store_true')
parser.add_argument(
    '--affine-metric', help='which metric to use for affine stages, use comma separated list for multiple image pairs (MI, Mattes, GC, MeanSquares, Demons)', default='Mattes')
parser.add_argument('--reg-pairs', help='number of image pairs for affine output', default=1, type=check_positive)
parser.add_argument('--reg-pairs-weights', help='either a single number for all weights, or a comma separated list of weights equal to reg_pairs', default=['1'], action=SplitArgsComma)
parser.add_argument('--no-masks', help='for linear registration outputs skip repeat stages with masks', action='store_true')

parser.add_argument('--override-shrink-factors', help='override calculation of optimal image pyramid with specified settings')
parser.add_argument('--override-smoothing-sigmas', help='override calculation of optimal image pyramid with specified settings')
parser.add_argument('--override-convergence', help='override calculation of optimal image pyramid with specified settings')


args = parser.parse_args()

# Setup inital inputs
min_resolution = args.min
max_size = args.max

affinemetric = args.affine_metric.rstrip(',').split(",")
affineweights= args.reg_pairs_weights

if len(affinemetric) != args.reg_pairs:
  if len(affinemetric) != 1:
    sys.exit("Number of metrics provided not equal to the number of image pairs, or 1")
  else:
    affinemetric = affinemetric * args.reg_pairs

if len(args.reg_pairs_weights) != args.reg_pairs:
  if len(args.reg_pairs_weights) != 1:
    sys.exit("Number of affine weights provided not equal to the number of image pairs, or 1")
  else:
    affineweights = affineweights * args.reg_pairs


if (args.output in ['affine', 'multilevel-halving', 'exhaustive-affine','lsq6', 'lsq9', 'lsq12', 'rigid', 'similarity']):
  if args.step_size:
     step_size = args.step_size
  else:
    step_size = min_resolution / 2
  if args.start_scale:
     starting_sigma = args.start_scale
  else:
    starting_sigma = max_size / 32
else:
  if args.step_size:
     step_size = args.step_size
  else:
    step_size = min_resolution / 2

  if args.start_scale:
     starting_sigma = args.start_scale
  else:
    starting_sigma = max_size / 12.3125 / 4

if args.close:
  starting_sigma = starting_sigma / 2

# Step down sigmas by step_size from starting_sigma to zero
sigmas = [ x * step_size for x in list(range(int(round(starting_sigma/step_size)),-1,-1)) ]

# Shrinks are 2x sigma, with a capped maximum shrink
shrinks = [ round(max(min(max_size/min_resolution/32, 2*x/min_resolution),1)) for x in sigmas ]

iterations = [ min(500,int(args.final_iterations * 3**(max(0,x - 1)))) for x in shrinks ]

if args.rough:
  to_remove = [i for i in range(len(shrinks)) if shrinks[i]==1]
  sigmas = [i for j, i in enumerate(sigmas) if j not in to_remove]
  shrinks = [i for j, i in enumerate(shrinks) if j not in to_remove]
  iterations = [i for j, i in enumerate(iterations) if j not in to_remove]
  sigmas.append('0')
  shrinks.append('1')
  iterations.append('0')

# Convert to strings
sigmas = [ str(x) for x in sigmas ]
shrinks = [ str(int(x)) for x in shrinks ]
iterations = [ str(x) for x in iterations ]

suffix = "mm"

if args.override_shrink_factors and args.override_smoothing_sigmas and args.override_convergence:
   if re.search("mm", args.override_smoothing_sigmas):
      suffix="mm"
      sigmas = args.override_smoothing_sigmas.strip("mm").split("x")
   elif re.search("vox", args.override_smoothing_sigmas):
      suffix="vox"
      sigmas = args.override_smoothing_sigmas.strip("vox").split("x")
   else:
      sigmas = args.override_smoothing_sigmas.split("x")
      suffix = ""
   shrinks = args.override_shrink_factors.split("x")
   iterations = args.override_convergence.split("x")

# Setup transforms
if args.output in ["multilevel-halving", "affine", "lsq12","exhaustive-affine"]:
  transforms = ["--transform Translation[ ",
                "--transform Rigid[ ",
                "--transform Similarity[ ",
                "--transform Affine[ "]
elif args.output in ["lsq9","similarity"]:
  transforms = ["--transform Translation[ ",
                "--transform Rigid[ ",
                "--transform Similarity[ ",
                "--transform Similarity[ "]
elif args.output in ["lsq6","rigid"]:
  transforms = ["--transform Translation[ ",
                "--transform Rigid[ ",
                "--transform Rigid[ ",
                "--transform Rigid[ "]
elif args.output in ['affine-plain']:
  transforms = ["--transform Rigid[ ",
                "--transform Similarity[ ",
                "--transform Affine[ "]
if not args.close:
  gradient_steps = [ 0.5, 0.33437015, 0.2236068, 0.1 ]
  gradient_steps_repeat = [ 0.5, 0.33437015, 0.14953488, 0.1 ]
else:
  gradient_steps = [ 0.1, 0.1, 0.1, 0.1 ]
  gradient_steps_repeat = [ 0.1, 0.1, 0.1, 0.1 ]

masks = ["--masks [ NOMASK,NOMASK ]",
         "--masks [ NOMASK,NOMASK ]",
         "--masks [ NOMASK,NOMASK ]",
         "--masks [ ${fixedmask},${movingmask} ]" ]

# If defined, repeat a stage previously done with a mask
repeatmask = [ False,
               False,
               "--masks [ ${fixedmask},${movingmask} ]",
               False ]

if args.output == 'exhaustive-affine' or args.output == 'affine-plain':
    for i, transform in enumerate(transforms):
            if i == len(transforms) - 1:
              print(transform + str(gradient_steps[i]) + " ]", end=' \\\n')
              for j in range(1, args.reg_pairs+1):
                  print("\t--metric {affinemetric}[ ${{fixedfile{j}}},${{movingfile{j}}},{affineweights},32,None,1,1 ]".format(j=j, affinemetric=affinemetric[j-1], affineweights=affineweights[j-1]), end=' \\\n')
              print("\t--convergence [ {},{},{} ]".format("x".join(iterations), args.convergence, args.convergence_window), end=' \\\n')
              print("\t--shrink-factors {}".format("x".join(shrinks)), end=' \\\n')
              if args.no_masks:
                print("\t--smoothing-sigmas {}{}".format("x".join(sigmas),suffix), end=' ')
              else:
                print("\t--smoothing-sigmas {}{}".format("x".join(sigmas),suffix), end=' \\\n')
                print("\t" + masks[i], end=' ')
            else:
              print(transform + str(gradient_steps[i]) + " ]", end=' \\\n')
              for j in range(1, args.reg_pairs+1):
                  print("\t--metric {affinemetric}[ ${{fixedfile{j}}},${{movingfile{j}}},{affineweights},32,None,1,1 ]".format(j=j, affinemetric=affinemetric[j-1], affineweights=affineweights[j-1]), end=' \\\n')
              print("\t--convergence [ {},{},{} ]".format("x".join(iterations), args.convergence, args.convergence_window), end=' \\\n')
              print("\t--shrink-factors {}".format("x".join(shrinks)), end=' \\\n')
              print("\t--smoothing-sigmas {}{}".format("x".join(sigmas),suffix), end=' \\\n')
              if not args.no_masks:
                print("\t" + masks[i], end=' \\\n')
                if repeatmask[i]:
                  print(transform + str(gradient_steps[i]) + " ]", end=' \\\n')
                  for j in range(1, args.reg_pairs+1):
                    print("\t--metric {affinemetric}[ ${{fixedfile{j}}},${{movingfile{j}}},{affineweights},32,None,1,1 ]".format(j=j, affinemetric=affinemetric[j-1], affineweights=affineweights[j-1]), end=' \\\n')
                  print("\t--convergence [ {},{},{} ]".format("x".join(iterations), args.convergence, args.convergence_window), end=' \\\n')
                  print("\t--shrink-factors {}".format("x".join(shrinks)), end=' \\\n')
                  print("\t--smoothing-sigmas {}{}".format("x".join(sigmas),suffix), end=' \\\n')
                  print("\t" + repeatmask[i], end=' \\\n')

elif args.output == 'twolevel_dbm':
    print("--reg-iterations {}".format("x".join(iterations)), end=' \\\n')
    print("--reg-shrinks {}".format("x".join(shrinks)), end=' \\\n')
    print("--reg-smoothing {}{}".format("x".join(sigmas),suffix), end=' ')

elif args.output == 'modelbuild':
    print("-q {}".format("x".join(iterations)), end=' \\\n')
    print("-f {}".format("x".join(shrinks)), end=' \\\n')
    print("-s {}{}".format("x".join(sigmas),suffix), end=' ')

elif args.output == 'generic':
      print("--convergence [ {},{},{} ]".format("x".join(iterations), args.convergence, args.convergence_window), end=' \\\n')
      print("--shrink-factors {}".format("x".join(shrinks)), end=' \\\n')
      print("--smoothing-sigmas {}{}".format("x".join(sigmas),suffix), end=' ')

elif args.output == 'volgenmodel':
    print("--convergence [ {}x0,{},{} ]".format("x".join(iterations[0:min(args.volgen_iteration+1,len(iterations))]), args.convergence, args.convergence_window), end=' \\\n')
    print("--shrink-factors {}x1".format("x".join(shrinks[0:min(args.volgen_iteration+1,len(shrinks))])), end=' \\\n')
    print("--smoothing-sigmas {}x0{}".format("x".join(sigmas[0:min(args.volgen_iteration+1,len(sigmas))]),suffix), end=' ')

else:

    slicestart = [ 0,
                   int(round(0.25*len(sigmas))),
                   int(round(0.50*len(sigmas))),
                   int(round(0.75*len(sigmas)))]

    sliceend = [   int(round(0.5*len(sigmas))),
                   int(round(0.75*len(sigmas))),
                   -2,
                   -1]

    if args.close:
      transforms = transforms[2:]
      slicestart = slicestart[0:3]
      sliceend = sliceend[1:]
      masks = masks[2:]
      repeatmask = repeatmask[2:]


    for i, transform in enumerate(transforms):
      if i == len(transforms) - 1:
        print(transform + str(gradient_steps[i]) + " ]", end=' \\\n')
        for j in range(1, args.reg_pairs+1):
          print("\t--metric {affinemetric}[ ${{fixedfile{j}}},${{movingfile{j}}},{affineweights},64,None,1,1 ]".format(j=j, affinemetric=affinemetric[j-1], affineweights=affineweights[j-1]), end=' \\\n')
        print("\t--convergence [ {},{},{} ]".format("x".join(iterations[slicestart[i]:]), args.convergence, args.convergence_window), end=' \\\n')
        print("\t--shrink-factors {}".format("x".join(shrinks[slicestart[i]:])), end=' \\\n')
        if args.no_masks:
          print("\t--smoothing-sigmas {}{}".format("x".join(sigmas[slicestart[i]:]),suffix), end=' ')
        else:
          print("\t--smoothing-sigmas {}{}".format("x".join(sigmas[slicestart[i]:]),suffix), end=' \\\n')
          print("\t" + masks[i], end=' ')
      else:
        print(transform + str(gradient_steps[i]) + " ]", end=' \\\n')
        for j in range(1, args.reg_pairs+1):
          print("\t--metric {affinemetric}[ ${{fixedfile{j}}},${{movingfile{j}}},{affineweights},32,None,1,1 ]".format(j=j, affinemetric=affinemetric[j-1], affineweights=affineweights[j-1]), end=' \\\n')
        print("\t--convergence [ {},{},{} ]".format("x".join(iterations[slicestart[i]:max(-slicestart[i]+1,sliceend[i])]), args.convergence, args.convergence_window), end=' \\\n')
        print("\t--shrink-factors {}".format("x".join(shrinks[slicestart[i]:max(-slicestart[i]+1,sliceend[i])])), end=' \\\n')
        print("\t--smoothing-sigmas {}{}".format("x".join(sigmas[slicestart[i]:max(-slicestart[i]+1,sliceend[i])]),suffix), end=' \\\n')
        if not args.no_masks:
          print("\t" + masks[i], end=' \\\n')
          if repeatmask[i]:
            print(transform + str(gradient_steps_repeat[i]) + " ]", end=' \\\n')
            for j in range(1, args.reg_pairs+1):
                print("\t--metric {affinemetric}[ ${{fixedfile{j}}},${{movingfile{j}}},{affineweights},32,None,1,1 ]".format(j=j, affinemetric=affinemetric[j-1], affineweights=affineweights[j-1]), end=' \\\n')
            print("\t--convergence [ {},{},{} ]".format("x".join(iterations[slicestart[i]:max(-slicestart[i]+1,sliceend[i])]), args.convergence, args.convergence_window), end=' \\\n')
            print("\t--shrink-factors {}".format("x".join(shrinks[slicestart[i]:max(-slicestart[i]+1,sliceend[i])])), end=' \\\n')
            print("\t--smoothing-sigmas {}{}".format("x".join(sigmas[slicestart[i]:max(-slicestart[i]+1,sliceend[i])]),suffix), end=' \\\n')
            print("\t" + repeatmask[i], end=' \\\n')
