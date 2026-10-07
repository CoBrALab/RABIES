#!/bin/bash
set -euo pipefail

# Installs the rat template set, derived from the SIGMA Wistar rat brain templates
# and atlases (Barriere et al. 2019, CC-BY-4.0) with scripts/gen_SIGMA_masks.py.
# Unlike the mouse set, this one is not installed automatically when a run needs it;
# `rabies install rat` runs this script. The container image runs it at build time.
#
# The bundle is downloaded rather than the SIGMA archive itself
# (https://zenodo.org/records/10635831) because the WM and CSF masks and the label CSV
# are derived from it by gen_SIGMA_masks.py, so every installation gets identical files
# without running the derivation.

out_dir=$1
mkdir -p "$out_dir"

bundle_url="https://github.com/CoBrALab/RABIES/releases/download/SIGMA-0.1/SIGMA.zip"
bundle_sha256="b280a61ad9f676a5959a63e693f6ff2707ea8f8f9ff4941badc978d344c70c4e"

curl -L --retry 5 --fail --silent --show-error "${bundle_url}" -o "${out_dir}/SIGMA.zip"

# macOS has shasum rather than sha256sum
if command -v sha256sum >/dev/null; then
    checksum=$(sha256sum "${out_dir}/SIGMA.zip" | cut -d' ' -f1)
else
    checksum=$(shasum -a 256 "${out_dir}/SIGMA.zip" | cut -d' ' -f1)
fi
if [ "${checksum}" != "${bundle_sha256}" ]; then
    rm "${out_dir}/SIGMA.zip"
    echo "The downloaded rat template set does not match the expected checksum." >&2
    exit 1
fi

unzip -o "${out_dir}/SIGMA.zip" -d "${out_dir}"
rm "${out_dir}/SIGMA.zip"
