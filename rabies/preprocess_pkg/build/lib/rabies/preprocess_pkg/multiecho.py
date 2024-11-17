# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright 2023 The NiPreps Developers <nipreps@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# We support and encourage derived works from this project, please read
# about our expectations at
#
#     https://www.nipreps.org/community/licensing/
#
"""
Multi-echo EPI
~~~~~~~~~~~~~~

For using multi-echo EPI data.

"""
import os
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.utility import IdentityInterface

from nipype import logging
from nipype.interfaces.base import (
    CommandLine,
    CommandLineInputSpec,
    File,
    TraitedSpec,
    traits,
)

LOGGER = logging.getLogger('nipype.interface')


class T2SMapInputSpec(CommandLineInputSpec):
    in_files = traits.List(
        File(exists=True),
        argstr='-d %s',
        position=1,
        mandatory=True,
        minlen=3,
        desc='multi-echo BOLD EPIs',
    )
    echo_times = traits.List(
        traits.Float, argstr='-e %s', position=2, mandatory=True, minlen=3, desc='echo times'
    )
    mask_file = File(argstr='--mask %s', position=3, desc='mask file', exists=True)
    fittype = traits.Enum(
        'curvefit',
        'loglin',
        argstr='--fittype %s',
        position=4,
        usedefault=True,
        desc=(
            'Desired fitting method: '
            '"loglin" means that a linear model is fit to the log of the data. '
            '"curvefit" means that a more computationally demanding '
            'monoexponential model is fit to the raw data.'
        ),
    )


class T2SMapOutputSpec(TraitedSpec):
    t2star_map = File(exists=True, desc='limited T2* map')
    s0_map = File(exists=True, desc='limited S0 map')
    optimal_comb = File(exists=True, desc='optimally combined ME-EPI time series')


class T2SMap(CommandLine):
    """
    Runs the tedana T2* workflow to generate an adaptive T2* map and create
    an optimally combined ME-EPI time series.

    Example
    =======

    >>> from fmriprep.interfaces import multiecho
    >>> t2smap = multiecho.T2SMap()
    >>> t2smap.inputs.in_files = ['sub-01_run-01_echo-1_bold.nii.gz',
    ...                           'sub-01_run-01_echo-2_bold.nii.gz',
    ...                           'sub-01_run-01_echo-3_bold.nii.gz']
    >>> t2smap.inputs.echo_times = [0.013, 0.027, 0.043]
    >>> t2smap.cmdline  # doctest: +ELLIPSIS
    't2smap -d sub-01_run-01_echo-1_bold.nii.gz sub-01_run-01_echo-2_bold.nii.gz \
sub-01_run-01_echo-3_bold.nii.gz -e 13.0 27.0 43.0 --fittype curvefit'

    """

    _cmd = 't2smap'
    input_spec = T2SMapInputSpec
    output_spec = T2SMapOutputSpec

    def _format_arg(self, name, trait_spec, value):
        if name == 'echo_times':
            value = [te * 1000 for te in value]
            return trait_spec.argstr % ' '.join([str(te) for te in value])
        elif name == 'in_files':
            return trait_spec.argstr % ' '.join(value)
        else:
            return super()._format_arg(name, trait_spec, value)

    def _list_outputs(self):
        outputs = self._outputs().get()
        out_dir = os.getcwd()
        outputs['t2star_map'] = os.path.join(out_dir, 'T2starmap.nii.gz')
        outputs['s0_map'] = os.path.join(out_dir, 'S0map.nii.gz')
        outputs['optimal_comb'] = os.path.join(out_dir, 'desc-optcom_bold.nii.gz')
        return outputs


class TedanaInputSpec(CommandLineInputSpec):
    in_files = traits.List(
        File(exists=True),
        argstr='-d %s',
        mandatory=True,
        desc='Multi-echo datasets',
    )
    echo_times = traits.List(
        traits.Float,
        argstr='-e %s',
        mandatory=True,
        desc='Echo times in seconds',
    )
    mask_file = File(argstr='--mask %s', exists=True, desc='Binary mask')
    fittype = traits.Enum(
        'loglin', 'curvefit', argstr='--fittype %s', desc='Desired T2*/S0 fitting method'
    )
    combmode = traits.Enum('t2s', argstr='--combmode %s', desc='Combination scheme for TEs')
    tedpca = traits.Either(
        traits.Enum('mdl', 'kic', 'aic', 'kundu', 'kundu-stabilize'),
        traits.Float(),
        traits.Int(),
        argstr='--tedpca %s',
        desc='Method with which to select components in TEDPCA',
    )
    out_dir = traits.Str(argstr='--out-dir %s', desc='Output directory')
    prefix = traits.Str(argstr='--prefix %s', desc='Prefix for filenames generated')
    # Include other options as needed


class TedanaOutputSpec(TraitedSpec):
    t2star_map = File(exists=True, desc='T2* map')
    s0_map = File(exists=True, desc='S0 map')
    optimal_comb = File(exists=True, desc='Optimally combined ME-EPI time series')
    denoised_data = File(exists=True, desc='Denoised ME-EPI time series')
    mixing_matrix = File(exists=True, desc='Mixing matrix')
    component_maps = File(exists=True, desc='Component maps')
    ica_metrics = File(exists=True, desc='ICA metrics')


class Tedana(CommandLine):
    """
    Runs the tedana workflow to perform multi-echo ICA on ME-EPI data.

    Example
    =======

    >>> from fmriprep.interfaces import multiecho
    >>> tedana = multiecho.Tedana()
    >>> tedana.inputs.in_files = ['sub-01_run-01_echo-1_bold.nii.gz',
    ...                           'sub-01_run-01_echo-2_bold.nii.gz',
    ...                           'sub-01_run-01_echo-3_bold.nii.gz']
    >>> tedana.inputs.echo_times = [0.013, 0.027, 0.043]
    >>> tedana.cmdline  # doctest: +ELLIPSIS
    'tedana -d sub-01_run-01_echo-1_bold.nii.gz sub-01_run-01_echo-2_bold.nii.gz \
sub-01_run-01_echo-3_bold.nii.gz -e 13.0 27.0 43.0'

    """

    _cmd = 'tedana'
    input_spec = TedanaInputSpec
    output_spec = TedanaOutputSpec

    def _format_arg(self, name, trait_spec, value):
        if name == 'echo_times':
            # tedana expects echo times in ms, so multiply by 1000
            value = [te * 1000 for te in value]
            return trait_spec.argstr % ' '.join([str(te) for te in value])
        elif name == 'in_files':
            return trait_spec.argstr % ' '.join(value)
        else:
            return super()._format_arg(name, trait_spec, value)

    def _list_outputs(self):
        outputs = self._outputs().get()
        out_dir = self.inputs.out_dir or os.getcwd()
        prefix = self.inputs.prefix or ''
        outputs['t2star_map'] = os.path.join(out_dir, f'{prefix}desc-t2s_map.nii.gz')
        outputs['s0_map'] = os.path.join(out_dir, f'{prefix}desc-s0_map.nii.gz')
        outputs['optimal_comb'] = os.path.join(out_dir, f'{prefix}desc-optcom_bold.nii.gz')
        outputs['denoised_data'] = os.path.join(out_dir, f'{prefix}desc-denoised_bold.nii.gz')
        outputs['mixing_matrix'] = os.path.join(out_dir, f'{prefix}desc-mixing.tsv')
        outputs['component_maps'] = os.path.join(out_dir, f'{prefix}desc-component_maps.nii.gz')
        outputs['ica_metrics'] = os.path.join(out_dir, f'{prefix}ica_metrics.tsv')
        # Include other outputs as needed
        return outputs


def create_multiecho_wf(opts):
    # Instantiate a workflow
    multiecho_wf = Workflow(name='multiecho_wf')

    # Create input nodes
    inputnode = Node(IdentityInterface(fields=['bold_list', 'te_list', 'mask']), name='inputnode')

    # Define output fields based on use_meica
    output_fields = ['output_bold', 't2star_map', 's0_map', 'optimal_comb']
    if getattr(opts, 'use_meica', False):
        output_fields.extend(['denoised_data', 'mixing_matrix', 'component_maps', 'ica_metrics', 'csf_mask', 'wm_mask'])
    outputnode = Node(IdentityInterface(fields=output_fields), name='outputnode')

    if getattr(opts, 'use_meica', False):
        # Create the Tedana node
        tedana_node = Node(Tedana(), name='tedana_node')
        # Set options if any
        if hasattr(opts, 'fittype') and opts.fittype:
            tedana_node.inputs.fittype = opts.fittype
        if hasattr(opts, 'combmode') and opts.combmode:
            tedana_node.inputs.combmode = opts.combmode
        if hasattr(opts, 'tedpca') and opts.tedpca:
            tedana_node.inputs.tedpca = opts.tedpca
        if hasattr(opts, 'out_dir') and opts.out_dir:
            tedana_node.inputs.out_dir = opts.out_dir
        if hasattr(opts, 'prefix') and opts.prefix:
            tedana_node.inputs.prefix = opts.prefix
        # Connect inputs
        multiecho_wf.connect([
            (inputnode, tedana_node, [('bold_list', 'in_files'),
                                      ('te_list', 'echo_times'),
                                      ('mask', 'mask_file')])
        ])
        # Connect the outputs to the outputnode
        multiecho_wf.connect([
            (tedana_node, outputnode, [
                ('denoised_data', 'output_bold'),
                ('t2star_map', 't2star_map'),
                ('s0_map', 's0_map'),
                ('optimal_comb', 'optimal_comb'),
                ('denoised_data', 'denoised_data'),
                ('mixing_matrix', 'mixing_matrix'),
                ('component_maps', 'component_maps'),
                ('ica_metrics', 'ica_metrics'),
                # Include other outputs as needed
            ])
        ])
    else:
        # Create the T2SMap node to get the optimally combined output
        optimal_combination = Node(T2SMap(fittype='curvefit'), name='optimal_combination')
        # Connect input nodes to the T2SMap node
        multiecho_wf.connect([
            (inputnode, optimal_combination, [('bold_list', 'in_files'),
                                              ('te_list', 'echo_times'),
                                              ('mask', 'mask_file')])
        ])
        # Connect the outputs to the outputnode
        multiecho_wf.connect([
            (optimal_combination, outputnode, [
                ('optimal_comb', 'output_bold'),
                ('t2star_map', 't2star_map'),
                ('s0_map', 's0_map'),
                ('optimal_comb', 'optimal_comb'),
            ])
        ])

    return multiecho_wf
