#!/bin/bash
set -euo pipefail

out_dir=$1
mkdir -p "$out_dir"


files=(
    "DSURQE_40micron_average.nii"
    "DSURQE_40micron_labels.nii"
    "DSURQE_40micron_mask.nii"
    "DSURQE_40micron_R_mapping.csv"
)
fallback_files=(
    "DSURQE_40micron_average.nii.gz"
    "DSURQE_40micron_labels.nii.gz"
    "DSURQE_40micron_mask.nii.gz"
    "DSURQE_40micron_R_mapping.csv"
)

# Primary URLs (mouseimaging.ca)
primary_urls=(
    "https://www.mouseimaging.ca/repo/DSURQE_40micron/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/nifti/DSURQE_40micron_average.nii"
    "https://www.mouseimaging.ca/repo/DSURQE_40micron/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/nifti/DSURQE_40micron_labels.nii"
    "https://www.mouseimaging.ca/repo/DSURQE_40micron/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/nifti/DSURQE_40micron_mask.nii"
    "https://www.mouseimaging.ca/repo/DSURQE_40micron/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/mappings/DSURQE_40micron_R_mapping.csv"
)

# Fallback URLs (github)
fallback_urls=(
    "https://github.com/CoBrALab/RABIES/releases/download/0.5.1/DSURQE_40micron.nii.gz"
    "https://github.com/CoBrALab/RABIES/releases/download/0.5.1/DSURQE_40micron_labels.nii.gz"
    "https://github.com/CoBrALab/RABIES/releases/download/0.5.1/DSURQE_40micron_mask.nii.gz"
    "https://github.com/CoBrALab/RABIES/releases/download/0.5.1/DSURQE_40micron_R_mapping.csv"
)

# # Try to download from mouseimaging first, if that doesn't work download from github
# primary_success=true
# for i in {0..3}; do
#     if ! curl -L --retry 5 --fail --silent --show-error -o "${out_dir}/${files[$i]}" "${primary_urls[$i]}"; then
#         primary_success=false
#         break
#     fi
# done
#
# if [ "$primary_success" = true ]; then
#   for f in "${out_dir}"/DSURQE_40micron_*.nii; do
#     gzip -f "$f"
#   done
# fi


# if these fail too exit 1
#
# temporary hotfix to force download from fallback_urls
# remove line 57 when primary_urls can be used again and uncomment lines 38-50
primary_success=false
if [ "$primary_success" = false ]; then
    for i in {0..3}; do
        if ! curl -L --retry 5 --fail --silent --show-error -o "${out_dir}/${fallback_files[$i]}" "${fallback_urls[$i]}"; then
            exit 1
        fi
    done
fi

# create regional masks
gen_DSURQE_masks.py ${out_dir}/DSURQE_40micron_labels.nii.gz ${out_dir}/DSURQE_40micron_R_mapping.csv ${out_dir}/DSURQE_40micron

# RABIES_extras: try Zenodo first, then a GitHub release mirror.
# Zenodo returns 403 to some clients (e.g. GitHub Actions runners), so the
# same files are mirrored as assets of a GitHub release.
extras_zenodo="https://zenodo.org/records/19069284/files"
extras_github="https://github.com/CoBrALab/RABIES/releases/download/RABIES_extras_19069284"

fetch_extra() {
    local name=$1
    for base in "$extras_zenodo" "$extras_github"; do
        echo "Downloading ${base}/${name}"
        if curl -L --retry 5 --retry-all-errors --fail --silent --show-error \
            -o "${out_dir}/${name}" "${base}/${name}"; then
            return 0
        fi
        echo "  failed: ${base}/${name}" >&2
    done
    echo "ERROR: could not download ${name} from any source" >&2
    exit 1
}

for f in melodic_IC.nii.gz vascular_mask.nii.gz EPICOMMON.zip DSURQE_seeds.zip; do
    fetch_extra "$f"
done

unzip "${out_dir}/EPICOMMON.zip" -d "${out_dir}"
unzip "${out_dir}/DSURQE_seeds.zip" -d "${out_dir}"
rm "${out_dir}/EPICOMMON.zip"
rm "${out_dir}/DSURQE_seeds.zip"
rm -r ${out_dir}/EPICOMMON/EPICOMMON2DSURQE/ # remove unnecessary transforms