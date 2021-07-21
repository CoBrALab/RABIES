from nipype.interfaces import utility as niu
import nipype.interfaces.ants as ants
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.utility import Function
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from .registration import run_antsRegistration
from .preprocess_visual_QC import PlotOverlap


def init_commonspace_reg_wf(opts, output_folder, transforms_datasink, num_scan, name='commonspace_reg_wf'):
    inputnode_iterable = pe.Node(niu.IdentityInterface(fields=['iter_name']),
                                        name="inputnode_iterable")
    inputnode = pe.Node(niu.IdentityInterface(fields=['moving_image', 'moving_mask', 'atlas_anat', 'atlas_mask']),
                                        name="inputnode")
    outputnode = pe.Node(niu.IdentityInterface(fields=['unbiased_template', 'native_mask', 'to_atlas_affine', 'to_atlas_warp', 'to_atlas_inverse_warp',
                                                       'native_to_unbiased_affine', 'native_to_unbiased_warp', 'native_to_unbiased_inverse_warp',
                                                       'native_to_commonspace_transform_list','native_to_commonspace_inverse_list',
                                                       'commonspace_to_native_transform_list','commonspace_to_native_inverse_list']),
                                        name="outputnode")
    workflow = pe.Workflow(name=name)

    atlas_reg = pe.Node(Function(input_names=['reg_method', 'moving_image', 'moving_mask', 'fixed_image', 'fixed_mask', 'rabies_data_type'],
                                    output_names=['affine', 'warp',
                                                  'inverse_warp', 'warped_image'],
                                    function=run_antsRegistration),
                           name='atlas_reg', mem_gb=2*opts.scale_min_memory)
    atlas_reg.inputs.rabies_data_type = opts.data_type
    atlas_reg.plugin_args = {
        'qsub_args': f'-pe smp {str(3*opts.min_proc)}', 'overwrite': True}
    atlas_reg.inputs.reg_method = str(opts.atlas_reg_script)

    prep_commonspace_transform_node = pe.Node(Function(input_names=['native_ref', 'atlas_mask', 'native_to_unbiased_affine',
                                                               'native_to_unbiased_warp','native_to_unbiased_inverse_warp',
                                                               'to_atlas_affine','to_atlas_warp','to_atlas_inverse_warp',
                                                               'fast_commonspace'],
                                               output_names=[
                                                   'native_mask', 'native_to_commonspace_transform_list','native_to_commonspace_inverse_list','commonspace_to_native_transform_list','commonspace_to_native_inverse_list'],
                                               function=prep_commonspace_transform),
                                      name='prep_commonspace_transform')
    prep_commonspace_transform_node.inputs.atlas_mask = str(opts.brain_mask)
    prep_commonspace_transform_node.inputs.fast_commonspace = opts.fast_commonspace

    workflow.connect([
        (inputnode, atlas_reg, [
            ("atlas_anat", "fixed_image"),
            ("atlas_mask", "fixed_mask"),
            ]),
        (atlas_reg, prep_commonspace_transform_node, [
            ("affine", "to_atlas_affine"),
            ("warp", "to_atlas_warp"),
            ("inverse_warp", "to_atlas_inverse_warp"),
            ]),
        (atlas_reg, outputnode, [
            ("affine", "to_atlas_affine"),
            ("warp", "to_atlas_warp"),
            ("inverse_warp", "to_atlas_inverse_warp"),
            ]),
        (prep_commonspace_transform_node, outputnode, [
            ('native_to_commonspace_transform_list', 'native_to_commonspace_transform_list'),
            ('native_to_commonspace_inverse_list', 'native_to_commonspace_inverse_list'),
            ('commonspace_to_native_transform_list', 'commonspace_to_native_transform_list'),
            ('commonspace_to_native_inverse_list', 'commonspace_to_native_inverse_list'),
            ('native_mask', 'native_mask'),
            ]),
        ])

    if opts.fast_commonspace:
        prep_commonspace_transform_node.inputs.native_to_unbiased_affine = None
        prep_commonspace_transform_node.inputs.native_to_unbiased_warp = None
        prep_commonspace_transform_node.inputs.native_to_unbiased_inverse_warp = None

        PlotOverlap_Native2Atlas_node = pe.Node(
            PlotOverlap(), name='PlotOverlap_Native2Atlas')
        PlotOverlap_Native2Atlas_node.inputs.out_dir = output_folder+f'/preprocess_QC_report/{name}.Native2Atlas/'

        workflow.connect([
            (inputnode, atlas_reg, [
                ("moving_image", "moving_image"),
                ]),
            (inputnode, prep_commonspace_transform_node, [
                ("moving_image", "native_ref"),
                ]),
            (inputnode_iterable, PlotOverlap_Native2Atlas_node, [
                ("iter_name", "name_source"),
                ]),
            (inputnode, PlotOverlap_Native2Atlas_node,[
                ("atlas_anat", "fixed"),
                ]),
            (atlas_reg, PlotOverlap_Native2Atlas_node, [
                ("warped_image", "moving"),
                ]),
            ])
        if opts.commonspace_masking:
            workflow.connect([
                (inputnode, atlas_reg, [
                    ("moving_mask", "moving_mask"),
                    ]),
                ])

    else:
        # setup a node to select the proper files associated with a given input scan for commonspace registration
        commonspace_selectfiles = pe.Node(Function(input_names=['filename', 'native_list', 'affine_list', 'warp_list', 'inverse_warp_list', 'warped_native_list'],
                                                   output_names=[
                                                       'native_ref', 'warped_native', 'native_to_unbiased_affine','native_to_unbiased_warp','native_to_unbiased_inverse_warp'],
                                                   function=select_commonspace_outputs),
                                          name='commonspace_selectfiles')

        generate_template = pe.Node(GenerateTemplate(masking=opts.commonspace_masking, output_folder=output_folder+f'/{name}.commonspace_datasink/', cluster_type=opts.plugin,
                                              ),
                                      name='generate_template', n_procs=opts.local_threads, mem_gb=1*num_scan*opts.scale_min_memory)

        PlotOverlap_Native2Unbiased_node = pe.Node(
            PlotOverlap(), name='PlotOverlap_Native2Unbiased')
        PlotOverlap_Native2Unbiased_node.inputs.out_dir = output_folder+f'/preprocess_QC_report/{name}.Native2Unbiased/'
        PlotOverlap_Unbiased2Atlas_node = pe.Node(
            PlotOverlap(), name='PlotOverlap_Unbiased2Atlas')
        PlotOverlap_Unbiased2Atlas_node.inputs.out_dir = output_folder+f'/preprocess_QC_report/{name}.Unbiased2Atlas'
        PlotOverlap_Unbiased2Atlas_node.inputs.name_source = ''

        workflow.connect([
            (inputnode, generate_template, [
                ("moving_image", "moving_image_list"),
                ("atlas_anat", "template_anat"),
                ]),
            (inputnode, commonspace_selectfiles, [
                ("moving_image", "native_list"),
                ]),
            (inputnode, PlotOverlap_Unbiased2Atlas_node,[
                ("atlas_anat", "fixed"),
                ]),
            (inputnode_iterable, commonspace_selectfiles, [
                ("iter_name", "filename"),
                ]),
            (inputnode_iterable, PlotOverlap_Native2Unbiased_node, [
                ("iter_name", "name_source"),
                ]),
            (generate_template, atlas_reg, [
                ("unbiased_template", "moving_image"),
                ]),
            (generate_template, commonspace_selectfiles, [
                ("affine_list", "affine_list"),
                ("warp_list", "warp_list"),
                ("inverse_warp_list", "inverse_warp_list"),
                ("warped_image_list", "warped_native_list"),
                ]),
            (generate_template, outputnode, [
                ("unbiased_template", "unbiased_template"),
                ]),
            (commonspace_selectfiles, prep_commonspace_transform_node, [
                ("native_ref", "native_ref"),
                ("native_to_unbiased_affine", "native_to_unbiased_affine"),
                ("native_to_unbiased_warp", "native_to_unbiased_warp"),
                ("native_to_unbiased_inverse_warp", "native_to_unbiased_inverse_warp"),
                ]),
            (commonspace_selectfiles, outputnode, [
                ("native_to_unbiased_affine", "native_to_unbiased_affine"),
                ("native_to_unbiased_warp", "native_to_unbiased_warp"),
                ("native_to_unbiased_inverse_warp", "native_to_unbiased_inverse_warp"),
                ]),
            (commonspace_selectfiles, PlotOverlap_Native2Unbiased_node, [
                ("warped_native", "moving"),
                ]),
            (generate_template, PlotOverlap_Native2Unbiased_node, [
                ("unbiased_template", "fixed"),
                ]),
            (atlas_reg, PlotOverlap_Unbiased2Atlas_node, [
                ("warped_image", "moving"),
                ]),
            ])
        if opts.commonspace_masking:
            workflow.connect([
                (inputnode, generate_template, [
                    ("moving_mask", "moving_mask_list"),
                    ]),
                (generate_template, atlas_reg, [
                    ("unbiased_mask", "moving_mask"),
                    ]),
                ])

    if opts.rabies_step == 'preprocess':
        if opts.fast_commonspace:
            workflow.connect([
                (atlas_reg, transforms_datasink, [
                    ("affine", "native_to_atlas_affine"),
                    ("warp", "native_to_atlas_warp"),
                    ("inverse_warp", "native_to_atlas_inverse_warp"),
                    ]),
                ])
        else:
            workflow.connect([
                (atlas_reg, transforms_datasink, [
                    ("affine", "unbiased_to_atlas_affine"),
                    ("warp", "unbiased_to_atlas_warp"),
                    ("inverse_warp", "unbiased_to_atlas_inverse_warp"),
                    ]),
                (commonspace_selectfiles, transforms_datasink, [
                    ("native_to_unbiased_affine", "native_to_unbiased_affine"),
                    ("native_to_unbiased_warp", "native_to_unbiased_warp"),
                    ("native_to_unbiased_inverse_warp", "native_to_unbiased_inverse_warp"),
                    ]),
                ])

    return workflow


def select_commonspace_outputs(filename, native_list, affine_list, warp_list, inverse_warp_list, warped_native_list):
    from rabies.preprocess_pkg.utils import select_from_list
    native_to_unbiased_affine = select_from_list(filename, affine_list)
    native_to_unbiased_warp = select_from_list(filename, warp_list)
    native_to_unbiased_inverse_warp = select_from_list(
        filename, inverse_warp_list)
    warped_native = select_from_list(filename, warped_native_list)
    native_ref = select_from_list(filename, native_list)
    return native_ref, warped_native, native_to_unbiased_affine, native_to_unbiased_warp, native_to_unbiased_inverse_warp


def prep_commonspace_transform(native_ref, atlas_mask, native_to_unbiased_affine,
                               native_to_unbiased_warp,native_to_unbiased_inverse_warp,
                               to_atlas_affine,to_atlas_warp,to_atlas_inverse_warp, fast_commonspace=False):

    if fast_commonspace:
        native_to_commonspace_transform_list=[to_atlas_warp,to_atlas_affine]
        native_to_commonspace_inverse_list=[0,0]
        commonspace_to_native_transform_list=[to_atlas_affine,to_atlas_inverse_warp]
        commonspace_to_native_inverse_list=[1,0]
    else:
        native_to_commonspace_transform_list=[to_atlas_warp,to_atlas_affine,native_to_unbiased_warp,native_to_unbiased_affine]
        native_to_commonspace_inverse_list=[0,0,0,0]
        commonspace_to_native_transform_list=[native_to_unbiased_affine,native_to_unbiased_inverse_warp,to_atlas_affine,to_atlas_inverse_warp]
        commonspace_to_native_inverse_list=[1,0,1,0]

    import os
    import pathlib  # Better path manipulation
    filename_split = pathlib.Path(
        native_ref).name.rsplit(".nii")
    native_mask = os.path.abspath(filename_split[0]+'_mask.nii.gz')
    from rabies.preprocess_pkg.utils import exec_applyTransforms
    exec_applyTransforms(commonspace_to_native_transform_list, commonspace_to_native_inverse_list, atlas_mask, native_ref, native_mask, mask=True)

    return native_mask, native_to_commonspace_transform_list,native_to_commonspace_inverse_list,commonspace_to_native_transform_list,commonspace_to_native_inverse_list


class GenerateTemplateInputSpec(BaseInterfaceInputSpec):
    moving_image_list = traits.List(exists=True, mandatory=True,
                            desc="List of anatomical images used for commonspace registration.")
    moving_mask_list = traits.List(exists=True,
                            desc="List of masks accompanying each image. (optional)")
    masking = traits.Bool(
        desc="Whether to use the masking option.")
    output_folder = traits.Str(
        exists=True, mandatory=True, desc="Path to output folder.")
    template_anat = File(exists=True, mandatory=True,
                         desc="Reference anatomical template to define the target space.")
    cluster_type = traits.Str(
        exists=True, mandatory=True, desc="Choose the type of cluster system to submit jobs to. Choices are local, sge, pbs, slurm.")

class GenerateTemplateOutputSpec(TraitedSpec):
    unbiased_template = File(
        exists=True, mandatory=True, desc="Output template generated from commonspace registration.")
    unbiased_mask = File(desc="Output mask generated along with the template. (optional)")
    affine_list = traits.List(exists=True, mandatory=True,
                              desc="List of affine transforms from anat to template space.")
    warp_list = traits.List(exists=True, mandatory=True,
                            desc="List of non-linear transforms from anat to template space.")
    inverse_warp_list = traits.List(exists=True, mandatory=True,
                                    desc="List of the inverse non-linear transforms from anat to template space.")
    warped_image_list = traits.List(exists=True, mandatory=True,
                                   desc="List of anatomical images warped to template space..")


class GenerateTemplate(BaseInterface):
    """

    """

    input_spec = GenerateTemplateInputSpec
    output_spec = GenerateTemplateOutputSpec

    def _run_interface(self, runtime):
        import os
        import pandas as pd
        import pathlib
        import multiprocessing
        from rabies.preprocess_pkg.utils import run_command

        cwd = os.getcwd()
        template_folder = self.inputs.output_folder

        command = f'mkdir -p {template_folder}'
        rc = run_command(command)

        from rabies.preprocess_pkg.utils import flatten_list
        merged = flatten_list(list(self.inputs.moving_image_list))
        # create a csv file of the input image list
        csv_path = cwd+'/commonspace_input_files.csv'
        df = pd.DataFrame(data=merged)
        df.to_csv(csv_path, header=False, sep=',', index=False)

        if self.inputs.masking:
            merged_masks = flatten_list(list(self.inputs.moving_mask_list))
            mask_csv_path = cwd+'/commonspace_input_masks.csv'
            df = pd.DataFrame(data=merged_masks)
            df.to_csv(mask_csv_path, header=False, sep=',', index=False)
            masks = f'--masks {mask_csv_path}'
        else:
            merged_masks = ['NULL']
            masks=''
            unbiased_mask='NULL'

        if len(merged) == 1:
            import logging
            log = logging.getLogger('root')
            log.info("Only a single scan was provided as input for commonspace registration. Commonspace registration "
                  "won't be run, and the output template will be the input scan.")

            # create an identity transform as a surrogate for the commonspace transforms
            import SimpleITK as sitk
            dimension = 3
            identity = sitk.Transform(dimension, sitk.sitkIdentity)

            file = merged[0]
            mask = merged_masks[0]
            filename_template = pathlib.Path(file).name.rsplit(".nii")[0]
            transform_file = template_folder+filename_template+'_identity.mat'
            sitk.WriteTransform(identity, transform_file)

            setattr(self, 'unbiased_template', file)
            setattr(self, 'unbiased_mask', mask)
            setattr(self, 'affine_list', [transform_file])
            setattr(self, 'warp_list', [transform_file])
            setattr(self, 'inverse_warp_list', [transform_file])
            setattr(self, 'warped_image_list', [file])

            return runtime


        # convert nipype plugin spec to match QBATCH
        plugin = self.inputs.cluster_type
        if plugin=='MultiProc' or plugin=='Linear':
            cluster_type='local'
            num_threads = multiprocessing.cpu_count()
        elif plugin=='SGE' or plugin=='SGEGraph':
            cluster_type='sge'
            num_threads = 1
        elif plugin=='PBS':
            cluster_type='PBS'
            num_threads = 1
        elif plugin=='SLURM' or plugin=='SLURMGraph':
            cluster_type='slurm'
            num_threads = 1
        else:
            raise ValueError("Plugin option must correspond to one of 'local', 'sge', 'pbs' or 'slurm'")

        command = f'QBATCH_SYSTEM={cluster_type} QBATCH_CORES={num_threads} modelbuild.sh \
            --float --average-type median --gradient-step 0.25 --iterations 2 --starting-target {self.inputs.template_anat} --stages rigid,affine,nlin \
            --output-dir {template_folder} --sharpen-type unsharp --block --debug {masks} {csv_path}'
        rc = run_command(command)


        unbiased_template = template_folder + \
            '/nlin/1/average/template_sharpen_shapeupdate.nii.gz'
        # verify that all outputs are present
        if not os.path.isfile(unbiased_template):
            raise ValueError(unbiased_template+" doesn't exists.")

        if self.inputs.masking:
            unbiased_mask = template_folder + \
                '/nlin/1/average/mask_shapeupdate.nii.gz'
            # verify that all outputs are present
            if not os.path.isfile(unbiased_mask):
                raise ValueError(unbiased_mask+" doesn't exists.")

        affine_list = []
        warp_list = []
        inverse_warp_list = []
        warped_image_list = []

        i = 0
        for file in merged:
            file = str(file)
            filename_template = pathlib.Path(file).name.rsplit(".nii")[0]
            native_to_unbiased_inverse_warp = f'{template_folder}/nlin/1/transforms/{filename_template}_1InverseWarp.nii.gz'
            if not os.path.isfile(native_to_unbiased_inverse_warp):
                raise ValueError(
                    native_to_unbiased_inverse_warp+" file doesn't exists.")
            native_to_unbiased_warp = f'{template_folder}/nlin/1/transforms/{filename_template}_1Warp.nii.gz'
            if not os.path.isfile(native_to_unbiased_warp):
                raise ValueError(native_to_unbiased_warp+" file doesn't exists.")
            native_to_unbiased_affine = f'{template_folder}/nlin/1/transforms/{filename_template}_0GenericAffine.mat'
            if not os.path.isfile(native_to_unbiased_affine):
                raise ValueError(native_to_unbiased_affine
                                 + " file doesn't exists.")
            warped_image = f'{template_folder}/nlin/1/resample/{filename_template}.nii.gz'
            if not os.path.isfile(warped_image):
                raise ValueError(warped_image
                                 + " file doesn't exists.")
            inverse_warp_list.append(native_to_unbiased_inverse_warp)
            warp_list.append(native_to_unbiased_warp)
            affine_list.append(native_to_unbiased_affine)
            warped_image_list.append(warped_image)
            i += 1

        setattr(self, 'unbiased_template', unbiased_template)
        setattr(self, 'unbiased_mask', unbiased_mask)
        setattr(self, 'affine_list', affine_list)
        setattr(self, 'warp_list', warp_list)
        setattr(self, 'inverse_warp_list', inverse_warp_list)
        setattr(self, 'warped_image_list', warped_image_list)

        return runtime

    def _list_outputs(self):
        return {'unbiased_template': getattr(self, 'unbiased_template'),
                'unbiased_mask': getattr(self, 'unbiased_mask'),
                'affine_list': getattr(self, 'affine_list'),
                'warp_list': getattr(self, 'warp_list'),
                'inverse_warp_list': getattr(self, 'inverse_warp_list'),
                'warped_image_list': getattr(self, 'warped_image_list'), }
