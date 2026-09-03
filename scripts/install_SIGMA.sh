#!/bin/bash
set -euo pipefail

# Installs the rat template set, derived from the SIGMA Wistar rat brain templates
# and atlases (Barriere et al. 2019, CC-BY-4.0) with scripts/gen_SIGMA_masks.py.
# Unlike the mouse set, this one is not installed automatically and is not baked
# into the container; `rabies install rat` runs this script.

out_dir=$1
mkdir -p "$out_dir"

bundle_url="https://github.com/CoBrALab/RABIES/releases/download/SIGMA-1.0/SIGMA.zip"

curl -L --retry 5 --fail --silent --show-error "${bundle_url}" -o "${out_dir}/SIGMA.zip"

unzip -o "${out_dir}/SIGMA.zip" -d "${out_dir}"
rm "${out_dir}/SIGMA.zip"
