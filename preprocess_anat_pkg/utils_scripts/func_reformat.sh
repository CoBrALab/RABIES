#!/bin/bash
#Script adapted from /opt/quarantine/minc-toolkit-extras/1.0/build/mouse-preprocessing-v2.sh only for flipping and centering
#usage:
#func_reformat.sh input.mnc OUTPUT_DIR

#Operations
# swaps zy to re-orient mouse
# flips x to fix left-right mixup
# centers brain in space

set -euo pipefail
set -v

tmpdir=$(mktemp -d)

cp $1 $tmpdir/input.mnc

input=$tmpdir/input.mnc
OUTPUT_DIR=$2

minc_modify_header $input -sinsert :history=‘’

volmash -swap zy $input $tmpdir/mash.mnc
volflip $tmpdir/mash.mnc $tmpdir/flip.mnc

clean_and_center_minc.pl $tmpdir/flip.mnc $tmpdir/centered.mnc

cp $tmpdir/centered.mnc $OUTPUT_DIR/reformat_$(basename $1)

rm -rf $tmpdir
