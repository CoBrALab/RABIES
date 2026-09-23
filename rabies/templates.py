"""Registry of the commonspace template sets distributed with RABIES.

Each template set groups the anatomical template and its associated masks and
atlas files into a single named unit selected with --template_set. A set holds
one 'anat' variant, used for the standard pipeline, and one 'epi' variant, used
when --bold_only bypasses the structural images.

Roles that a set does not provide are declared as None. Downstream operations
depending on a missing role are disabled rather than falling back to another
set, since mixing files across sets silently registers data to the wrong space.
"""

import os

# setting all default template files
if 'XDG_DATA_HOME' in os.environ.keys():
    rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
else:
    rabies_path = os.environ['HOME']+'/.local/share/rabies'

# roles resolved during preprocessing; brain_mask and the masks after it must be binary
PREPROCESS_ROLES = ['anat_template', 'brain_mask', 'WM_mask', 'CSF_mask', 'vascular_mask']
# roles resolved during analysis
ANALYSIS_ROLES = ['labels', 'mapping', 'prior_maps']

MOUSE_SEED_NAMES = [
        'ACA_seed', 'ECT_AI_seed', 'HY_seed', 'ORB_limbic_seed', 'SS_frontal_seed',
        'AMYG_seed', 'RSP_seed', 'THAL_seed', 'basal_ganglia_seed', 'HIP_seed',
        'MO_seed', 'SS_dorsal_seed', 'VIS_seed']

TEMPLATE_SETS = {
    'mouse': {
        'description': "DSURQE mouse atlas (https://www.mouseimaging.ca/repo/DSURQE_40micron/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/)",
        # the mouse set is installed on demand at the start of every run
        'auto_install': True,
        'install_script': 'install_DSURQE.sh',
        'seed_names': MOUSE_SEED_NAMES,
        'anat': {
            'anat_template': f"{rabies_path}/DSURQE_40micron_average.nii.gz",
            'brain_mask': f"{rabies_path}/DSURQE_40micron_mask.nii.gz",
            'WM_mask': f"{rabies_path}/DSURQE_40micron_eroded_WM_mask.nii.gz",
            'CSF_mask': f"{rabies_path}/DSURQE_40micron_eroded_CSF_mask.nii.gz",
            'vascular_mask': f"{rabies_path}/vascular_mask.nii.gz",
            'labels': f"{rabies_path}/DSURQE_40micron_labels.nii.gz",
            'mapping': f"{rabies_path}/DSURQE_40micron_R_mapping.csv",
            'prior_maps': f"{rabies_path}/melodic_IC.nii.gz",
            'seed_dir': f"{rabies_path}/DSURQE_seeds",
            'seed_suffix': '_left.nii.gz',
            },
        'epi': {
            'anat_template': f"{rabies_path}/EPICOMMON/EPICOMMON_atlas/EPICOMMON_template.nii.gz",
            'brain_mask': f"{rabies_path}/EPICOMMON/EPICOMMON_atlas/EPICOMMON_brain_mask.nii.gz",
            'WM_mask': f"{rabies_path}/EPICOMMON/EPICOMMON_atlas/EPICOMMON_WM_mask.nii.gz",
            'CSF_mask': f"{rabies_path}/EPICOMMON/EPICOMMON_atlas/EPICOMMON_CSF_mask.nii.gz",
            'vascular_mask': f"{rabies_path}/EPICOMMON/EPICOMMON_atlas/EPICOMMON_vascular_mask.nii.gz",
            'labels': f"{rabies_path}/EPICOMMON/EPICOMMON_atlas/EPICOMMON_labels.nii.gz",
            'mapping': f"{rabies_path}/DSURQE_40micron_R_mapping.csv",
            'prior_maps': f"{rabies_path}/EPICOMMON/EPICOMMON_atlas/melodic_IC_resampled.nii.gz",
            'seed_dir': f"{rabies_path}/DSURQE_seeds/EPICOMMON_resampled",
            'seed_suffix': '_left_EPICOMMON_resampled.nii.gz',
            },
        },
    'rat': {
        'description': "SIGMA Wistar rat atlas (https://github.com/DavidBarriere/SIGMA-Rat-Brain-Templates-and-Atlases)",
        # the rat set is not installed automatically; `rabies install rat` fetches it
        'auto_install': False,
        'install_script': 'install_SIGMA.sh',
        'seed_names': [],
        'anat': {
            'anat_template': f"{rabies_path}/SIGMA/SIGMA_InVivo_Anatomical_Brain_template.nii.gz",
            'brain_mask': f"{rabies_path}/SIGMA/SIGMA_InVivo_Anatomical_Brain_mask.nii.gz",
            # union of the white matter structures of the SIGMA anatomical atlas
            'WM_mask': f"{rabies_path}/SIGMA/SIGMA_InVivo_Anatomical_Brain_wm_mask.nii.gz",
            # union of the ventricles of the SIGMA anatomical atlas, not eroded
            'CSF_mask': f"{rabies_path}/SIGMA/SIGMA_InVivo_Anatomical_Brain_csf_mask.nii.gz",
            # SIGMA provides no vessel segmentation
            'vascular_mask': None,
            'labels': f"{rabies_path}/SIGMA/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz",
            # derived from the ITK-SNAP label description shipped with the atlas
            'mapping': f"{rabies_path}/SIGMA/SIGMA_InVivo_Anatomical_Brain_Atlas_mapping.csv",
            # no group ICA priors have been published in SIGMA space
            'prior_maps': None,
            'seed_dir': None,
            'seed_suffix': None,
            },
        'epi': {
            'anat_template': f"{rabies_path}/SIGMA/SIGMA_InVivo_Functional_Brain_epi.nii.gz",
            'brain_mask': f"{rabies_path}/SIGMA/SIGMA_InVivo_Functional_Brain_mask.nii.gz",
            # the functional atlas has no white matter structures to build a mask from, and
            # thresholding the probabilistic white matter map pulls in thalamus
            'WM_mask': None,
            # the anatomical CSF mask resampled onto this grid, since the functional atlas
            # has no ventricles
            'CSF_mask': f"{rabies_path}/SIGMA/SIGMA_InVivo_Functional_Brain_csf_mask.nii.gz",
            'vascular_mask': None,
            # the 59-ROI SIGMA functional parcellation, the only atlas in this space
            'labels': f"{rabies_path}/SIGMA/SIGMA_InVivo_Functional_Brain_Atlas.nii.gz",
            'mapping': f"{rabies_path}/SIGMA/SIGMA_Functional_Brain_Atlas_mapping.csv",
            'prior_maps': None,
            'seed_dir': None,
            'seed_suffix': None,
            },
        },
    }

TEMPLATE_SET_NAMES = list(TEMPLATE_SETS.keys())

# the set used when --template_set is not given; RABIES does not detect the species of
# the data, so falling back to it is logged as a warning
DEFAULT_TEMPLATE_SET = 'mouse'


def describe_sets():
    """Return a help-text listing of the available template sets."""
    return ''.join([f"* {name}: {TEMPLATE_SETS[name]['description']}\n"
                    for name in TEMPLATE_SET_NAMES])


def get_template_set(preprocess_opts):
    """Return the template set a preprocessing run was registered to."""
    # runs preprocessed before --template_set existed carry no such attribute, and
    # were necessarily run with the mouse files that were the only ones available
    return getattr(preprocess_opts, 'template_set', 'mouse')


def get_variant(template_set, bold_only=False):
    """Return the dictionary of files for a template set, given the pipeline variant."""
    if template_set not in TEMPLATE_SETS.keys():
        raise ValueError(
            f"Unknown --template_set {template_set}. Available sets: {TEMPLATE_SET_NAMES}.")
    return TEMPLATE_SETS[template_set]['epi' if bold_only else 'anat']


def resolve(template_set, role, bold_only=False):
    """Return the file a template set provides for a role, or None if it provides none."""
    variant = get_variant(template_set, bold_only=bold_only)
    if role not in variant.keys():
        raise ValueError(f"Unknown template role {role}.")
    return variant[role]


def registered_to_set_template(preprocess_opts):
    """Return whether a preprocessing run registered its data to its set's own template.

    The set's atlas files are only aligned with that template. The template file of
    the run is compared rather than recorded when the options are resolved, since
    --inherit_unbiased_template replaces it afterwards with the template of an earlier
    run, and runs preprocessed before template sets existed record nothing.
    """
    set_template = resolve(get_template_set(preprocess_opts), 'anat_template',
                           bold_only=preprocess_opts.bold_only)
    return str(preprocess_opts.anat_template) == set_template


def analysis_default(preprocess_opts, role):
    """Return the default file for an analysis role, given the preprocessing run.

    The file comes from the set the run was registered to. None is returned when the
    set does not provide the role, or when the run registered its data to a template
    from elsewhere, which the set's files are not aligned with.
    """
    if not registered_to_set_template(preprocess_opts):
        return None
    return resolve(get_template_set(preprocess_opts), role,
                   bold_only=preprocess_opts.bold_only)


def seed_names(template_set):
    """Return the pre-built seed names available for a template set."""
    return TEMPLATE_SETS[template_set]['seed_names']


def install_script(template_set):
    """Return the name of the script that downloads the files of a template set."""
    return TEMPLATE_SETS[template_set]['install_script']


def auto_install(template_set):
    """Return whether a template set is installed on demand when a run needs it."""
    return TEMPLATE_SETS[template_set]['auto_install']


def seed_file(template_set, seed, bold_only=False):
    """Return the path of a pre-built seed, or None if the set ships no seeds."""
    variant = get_variant(template_set, bold_only=bold_only)
    if variant['seed_dir'] is None:
        return None
    return f"{variant['seed_dir']}/{seed}{variant['seed_suffix']}"


def required_files(template_set):
    """Return every file a complete installation of a template set must contain."""
    files = []
    for variant_name in ['anat', 'epi']:
        variant = TEMPLATE_SETS[template_set][variant_name]
        for role in PREPROCESS_ROLES+ANALYSIS_ROLES:
            f = variant.get(role)
            if f is not None and f not in files:
                files.append(f)
        for seed in seed_names(template_set):
            files.append(seed_file(template_set, seed,
                                   bold_only=(variant_name == 'epi')))
    return files


def missing_files(template_set):
    """Return the files of a template set that are not installed."""
    return [f for f in required_files(template_set) if not os.path.isfile(f)]


def resolve_options(opts, log):
    # Fill in the template files that the user did not provide from the selected
    # --template_set. A file given explicitly overrides the one from the set; a role
    # that neither provides is left as None, which blocks the downstream operations
    # depending on it while still allowing preprocessing to run.
    # The roles the user provided are recorded first, since the loop below overwrites
    # them with the resolved files.
    provided = [role for role in PREPROCESS_ROLES if getattr(opts, role) is not None]

    if 'anat_template' in provided and 'brain_mask' not in provided:
        raise ValueError("--anat_template was provided, but not --brain_mask "
                         "- it is necessary to provide a brain mask matching the template.")

    for role in PREPROCESS_ROLES:
        if role in provided:
            # make sure we have absolute paths
            setattr(opts, role, os.path.abspath(getattr(opts, role)))
            log.info(f'--{role} overrides the file from --template_set {opts.template_set}.')
            continue
        if 'anat_template' in provided:
            # the set's masks do not match a user-provided template, so they are not used
            setattr(opts, role, None)
            continue
        setattr(opts, role, resolve(opts.template_set, role, bold_only=opts.bold_only))

    log.info(f"COMMONSPACE TEMPLATE SET: {opts.template_set}"
             + (" (--bold_only EPI variant)" if opts.bold_only else ""))
    for role in PREPROCESS_ROLES:
        opt_file = getattr(opts, role)
        log.info(f"    --{role}: {opt_file if opt_file is not None else 'not available'}")

    # the default only matters when the set's own files are used: an inherited run fixes
    # the set, and a user-provided template replaces its files
    if not opts.explicit_template_set and opts.inherit_unbiased_template == 'none' \
            and 'anat_template' not in provided:
        log.warning(
            "\n############################################# WARNING\n"
            f"No --template_set was selected, so the {opts.template_set} template set is used. "
            "RABIES does not detect the species of the data: data registered to the template "
            "of another species runs to completion, but produces meaningless outputs. Select "
            f"the set matching the species with --template_set ({', '.join(TEMPLATE_SET_NAMES)}), "
            "which also silences this warning."
            "\n############################################# WARNING\n")


def prior_maps_required(opts):
    """Return whether the analyses selected in `opts` require --prior_maps.

    Dual regression and neural prior recovery fit the prior maps to the data.
    --data_diagnosis uses them when available but runs without them, and seed-based
    connectivity, FC matrices and group ICA do not use them, so a set that provides no
    priors does not block those.
    """
    return bool(opts.DR_ICA
                or opts.NPR_temporal_comp > -1 or opts.NPR_spatial_comp > -1
                or opts.optimize_NPR['apply'])
