"""Tests for the commonspace template set registry and option resolution.

These cover the logic that decides which template files a run uses. Getting it
wrong does not crash the pipeline, it registers data to the wrong space, so the
resolution rules are checked directly rather than through a pipeline run.

Run with: python -m unittest discover tests
"""

import os
import unittest
from argparse import Namespace

from rabies import templates


class StubLog:
    def __init__(self):
        self.messages = []
        self.warnings = []

    def info(self, message):
        self.messages.append(message)

    def warning(self, message):
        self.messages.append(message)
        self.warnings.append(message)


def preprocess_options(**overrides):
    opts = Namespace(template_set='mouse', bold_only=False, explicit_template_set=True,
                     inherit_unbiased_template='none')
    for role in templates.PREPROCESS_ROLES:
        setattr(opts, role, None)
    for key, value in overrides.items():
        setattr(opts, key, value)
    return opts


def analysis_options(**overrides):
    opts = Namespace(DR_ICA=False, data_diagnosis=False, NPR_temporal_comp=-1,
                     NPR_spatial_comp=-1, optimize_NPR={'apply': False})
    for key, value in overrides.items():
        setattr(opts, key, value)
    return opts


class TestRegistry(unittest.TestCase):

    def test_every_set_defines_every_role(self):
        for name in templates.TEMPLATE_SET_NAMES:
            for variant in ['anat', 'epi']:
                for role in templates.PREPROCESS_ROLES+templates.ANALYSIS_ROLES:
                    self.assertIn(role, templates.TEMPLATE_SETS[name][variant],
                                  f"{name}/{variant} is missing {role}")

    def test_a_set_always_provides_a_template_and_a_brain_mask(self):
        for name in templates.TEMPLATE_SET_NAMES:
            for bold_only in [False, True]:
                self.assertIsNotNone(templates.resolve(name, 'anat_template', bold_only=bold_only))
                self.assertIsNotNone(templates.resolve(name, 'brain_mask', bold_only=bold_only))

    def test_the_epi_variant_differs_from_the_anat_variant(self):
        # --bold_only must not silently reuse the structural template
        for name in templates.TEMPLATE_SET_NAMES:
            self.assertNotEqual(templates.resolve(name, 'anat_template'),
                                templates.resolve(name, 'anat_template', bold_only=True))

    def test_sets_share_no_template_file(self):
        # a file resolving to another set's would register data to the wrong species
        mouse = {f for f in templates.required_files('mouse')}
        rat = {f for f in templates.required_files('rat')}
        self.assertEqual(mouse & rat, set())

    def test_rat_provides_no_vascular_mask_priors_or_seeds(self):
        self.assertIsNone(templates.resolve('rat', 'vascular_mask'))
        self.assertIsNone(templates.resolve('rat', 'prior_maps'))
        self.assertEqual(templates.seed_names('rat'), [])
        self.assertIsNone(templates.seed_file('rat', 'HIP_seed'))

    def test_rat_provides_a_white_matter_mask_only_for_the_anat_variant(self):
        # the functional atlas has no white matter structures to build the mask from
        self.assertIsNotNone(templates.resolve('rat', 'WM_mask'))
        self.assertIsNone(templates.resolve('rat', 'WM_mask', bold_only=True))

    def test_mouse_seeds_differ_between_variants(self):
        anat_seed = templates.seed_file('mouse', 'HIP_seed')
        epi_seed = templates.seed_file('mouse', 'HIP_seed', bold_only=True)
        self.assertNotEqual(anat_seed, epi_seed)

    def test_required_files_excludes_roles_a_set_does_not_provide(self):
        self.assertNotIn(None, templates.required_files('rat'))

    def test_the_set_of_a_preprocessing_run_is_read_back(self):
        self.assertEqual(templates.get_template_set(Namespace(template_set='rat')), 'rat')

    def test_a_run_from_before_template_sets_used_the_mouse_files(self):
        self.assertEqual(templates.get_template_set(Namespace()), 'mouse')

    def test_required_files_include_the_label_mapping(self):
        for name in templates.TEMPLATE_SET_NAMES:
            self.assertIn(templates.resolve(name, 'mapping'), templates.required_files(name))

    def test_only_the_mouse_set_installs_itself(self):
        # compute nodes often have no network, so the rat set is installed explicitly
        self.assertTrue(templates.auto_install('mouse'))
        self.assertFalse(templates.auto_install('rat'))

    def test_each_set_names_its_install_script(self):
        self.assertEqual(templates.install_script('mouse'), 'install_DSURQE.sh')
        self.assertEqual(templates.install_script('rat'), 'install_SIGMA.sh')

    def test_unknown_set_and_role_are_rejected(self):
        with self.assertRaises(ValueError):
            templates.resolve('hamster', 'anat_template')
        with self.assertRaises(ValueError):
            templates.resolve('mouse', 'not_a_role')

    def test_describe_sets_names_every_set(self):
        description = templates.describe_sets()
        for name in templates.TEMPLATE_SET_NAMES:
            self.assertIn(name, description)


class TestResolveOptions(unittest.TestCase):

    def test_defaults_come_from_the_selected_set(self):
        for name in templates.TEMPLATE_SET_NAMES:
            opts = preprocess_options(template_set=name)
            templates.resolve_options(opts, StubLog())
            for role in templates.PREPROCESS_ROLES:
                self.assertEqual(getattr(opts, role),
                                 templates.resolve(name, role))

    def test_bold_only_selects_the_epi_variant(self):
        opts = preprocess_options(template_set='rat', bold_only=True)
        templates.resolve_options(opts, StubLog())
        self.assertEqual(opts.anat_template,
                         templates.resolve('rat', 'anat_template', bold_only=True))

    def test_an_individual_file_overrides_only_itself(self):
        # the case the sentinel comparisons could not express: keep the set's
        # template while replacing one of its masks
        opts = preprocess_options(template_set='rat', brain_mask='my_mask.nii.gz')
        templates.resolve_options(opts, StubLog())
        self.assertEqual(opts.anat_template, templates.resolve('rat', 'anat_template'))
        self.assertEqual(opts.brain_mask, os.path.abspath('my_mask.nii.gz'))
        self.assertEqual(opts.WM_mask, templates.resolve('rat', 'WM_mask'))

    def test_a_user_template_drops_the_sets_masks(self):
        # the set's masks are not aligned with a template from elsewhere
        opts = preprocess_options(anat_template='t.nii.gz', brain_mask='m.nii.gz')
        templates.resolve_options(opts, StubLog())
        self.assertEqual(opts.anat_template, os.path.abspath('t.nii.gz'))
        self.assertEqual(opts.brain_mask, os.path.abspath('m.nii.gz'))
        for role in ['WM_mask', 'CSF_mask', 'vascular_mask']:
            self.assertIsNone(getattr(opts, role))

    def test_a_user_template_keeps_the_masks_given_with_it(self):
        opts = preprocess_options(anat_template='t.nii.gz', brain_mask='m.nii.gz',
                                  WM_mask='wm.nii.gz')
        templates.resolve_options(opts, StubLog())
        self.assertEqual(opts.WM_mask, os.path.abspath('wm.nii.gz'))
        self.assertIsNone(opts.CSF_mask)

    def test_a_user_template_without_a_brain_mask_is_rejected(self):
        opts = preprocess_options(anat_template='t.nii.gz')
        with self.assertRaises(ValueError):
            templates.resolve_options(opts, StubLog())

    def test_the_resolved_set_is_logged(self):
        log = StubLog()
        opts = preprocess_options(template_set='rat')
        templates.resolve_options(opts, log)
        self.assertTrue(any('rat' in message for message in log.messages))

    def test_falling_back_to_the_default_set_is_warned_about(self):
        # rat data run without --template_set is registered to the mouse template
        log = StubLog()
        opts = preprocess_options(template_set=templates.DEFAULT_TEMPLATE_SET,
                                  explicit_template_set=False)
        templates.resolve_options(opts, log)
        self.assertEqual(len(log.warnings), 1)
        self.assertIn('--template_set', log.warnings[0])

    def test_a_selected_set_is_not_warned_about(self):
        log = StubLog()
        templates.resolve_options(preprocess_options(template_set='mouse'), log)
        self.assertEqual(log.warnings, [])

    def test_the_default_is_not_warned_about_when_its_files_are_not_used(self):
        # an inherited run fixes the set, and a user template replaces the set's files
        for overrides in [{'inherit_unbiased_template': '/previous/run'},
                          {'anat_template': 't.nii.gz', 'brain_mask': 'm.nii.gz'}]:
            log = StubLog()
            opts = preprocess_options(explicit_template_set=False, **overrides)
            templates.resolve_options(opts, log)
            self.assertEqual(log.warnings, [], f"{overrides} should not warn")


class TestRegisteredToSetTemplate(unittest.TestCase):
    # the analysis stage applies a set's atlas only to data registered to that set's
    # template, since the atlas is not aligned with any other

    def test_a_run_on_the_set_template_is_registered_to_it(self):
        for name in templates.TEMPLATE_SET_NAMES:
            for bold_only in [False, True]:
                opts = preprocess_options(template_set=name, bold_only=bold_only)
                templates.resolve_options(opts, StubLog())
                self.assertTrue(templates.registered_to_set_template(opts),
                                f"{name} bold_only={bold_only}")

    def test_a_run_on_a_user_template_is_not(self):
        opts = preprocess_options(anat_template='t.nii.gz', brain_mask='m.nii.gz')
        templates.resolve_options(opts, StubLog())
        self.assertFalse(templates.registered_to_set_template(opts))

    def test_a_run_inheriting_a_user_template_is_not(self):
        # --inherit_unbiased_template replaces the template after the options are resolved
        opts = preprocess_options()
        templates.resolve_options(opts, StubLog())
        opts.anat_template = os.path.abspath('inherited_template.nii.gz')
        self.assertFalse(templates.registered_to_set_template(opts))

    def test_a_run_from_before_template_sets_is_checked_against_the_mouse_set(self):
        opts = Namespace(bold_only=False,
                         anat_template=templates.resolve('mouse', 'anat_template'))
        self.assertTrue(templates.registered_to_set_template(opts))
        opts.anat_template = os.path.abspath('t.nii.gz')
        self.assertFalse(templates.registered_to_set_template(opts))


class TestPriorMapsRequired(unittest.TestCase):

    def test_analyses_that_do_not_fit_priors_do_not_require_them(self):
        # seed-based connectivity, FC matrices and --data_diagnosis must run on a set that
        # ships no priors
        self.assertFalse(templates.prior_maps_required(analysis_options()))
        self.assertFalse(templates.prior_maps_required(analysis_options(data_diagnosis=True)))

    def test_each_prior_based_analysis_requires_them(self):
        for overrides in [{'DR_ICA': True},
                          {'NPR_temporal_comp': 0}, {'NPR_spatial_comp': 0},
                          {'optimize_NPR': {'apply': True}}]:
            self.assertTrue(templates.prior_maps_required(analysis_options(**overrides)),
                            f"{overrides} should require prior maps")


if __name__ == '__main__':
    unittest.main()
