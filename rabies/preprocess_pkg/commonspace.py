from nipype.interfaces import utility as niu
import nipype.interfaces.ants as ants
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.utility import Function
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


class GenerateTemplateInputSpec(BaseInterfaceInputSpec):
    moving_image = traits.List(exists=True, mandatory=True,
                            desc="List of anatomical images used for commonspace registration.")
    output_folder = traits.Str(
        exists=True, mandatory=True, desc="Path to output folder.")
    template_anat = File(exists=True, mandatory=True,
                         desc="Reference anatomical template to define the target space.")
    cluster_type = traits.Str(
        exists=True, mandatory=True, desc="Choose the type of cluster system to submit jobs to. Choices are local, sge, pbs, slurm.")
    #walltime = traits.Str(
    #    exists=True, mandatory=True, desc="Option for job submission specifying requested time per pairwise registration.")

class GenerateTemplateOutputSpec(TraitedSpec):
    warped_image = File(
        exists=True, desc="Output template generated from commonspace registration.")
    affine_list = traits.List(exists=True, mandatory=True,
                              desc="List of affine transforms from anat to template space.")
    warp_list = traits.List(exists=True, mandatory=True,
                            desc="List of non-linear transforms from anat to template space.")
    inverse_warp_list = traits.List(exists=True, mandatory=True,
                                    desc="List of the inverse non-linear transforms from anat to template space.")
    warped_anat_list = traits.List(exists=True, mandatory=True,
                                   desc="List of anatomical images warped to template space..")


class GenerateTemplate(BaseInterface):
    """
    Runs commonspace registration using ants_dbm.
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

        command = 'mkdir -p %s' % (template_folder,)
        rc = run_command(command)

        # create a csv file of the input image list
        csv_path = cwd+'/commonspace_input_files.csv'
        from rabies.preprocess_pkg.utils import flatten_list
        merged = flatten_list(list(self.inputs.moving_image))

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
            filename_template = pathlib.Path(file).name.rsplit(".nii")[0]
            transform_file = template_folder+filename_template+'_identity.mat'
            sitk.WriteTransform(identity, transform_file)

            setattr(self, 'warped_image', file)
            setattr(self, 'affine_list', [transform_file])
            setattr(self, 'warp_list', [transform_file])
            setattr(self, 'inverse_warp_list', [transform_file])
            setattr(self, 'warped_anat_list', [file])

            return runtime

        df = pd.DataFrame(data=merged)
        df.to_csv(csv_path, header=False, sep=',', index=False)

        global OPTIONS
        OPTIONS = [os.environ.get("QBATCH_OPTIONS")] if os.environ.get(
            "QBATCH_OPTIONS") else []

        # convert nipype plugin spec to match QBATCH
        plugin = self.inputs.cluster_type
        if plugin=='MultiProc' or plugin=='Linear':
            cluster_type='local'
            num_threads = multiprocessing.cpu_count()
        elif plugin=='SGE' or plugin=='SGEGraph':
            cluster_type='sge'
            num_threads = 1
            OPTIONS.append("'-sync y'")
        elif plugin=='PBS':
            cluster_type='PBS'
            num_threads = 1
            OPTIONS.append("'-Wblock=true'")
        elif plugin=='SLURM' or plugin=='SLURMGraph':
            cluster_type='slurm'
            num_threads = 1
            OPTIONS.append("'--wait'")
        else:
            raise ValueError("Plugin option must correspond to one of 'local', 'sge', 'pbs' or 'slurm'")

        qbatch_extras='%s QBATCH_SYSTEM=%s QBATCH_CORES=%s' % ('QBATCH_OPTIONS='+','.join(OPTIONS), cluster_type, num_threads,)
        command = '%s modelbuild.sh \
            --float --average-type mean --gradient-step 0.25 --iterations 3 --starting-target %s --stages nlin \
            --output-dir %s --sharpen-type none --debug %s' % (
            qbatch_extras, self.inputs.template_anat, template_folder, csv_path)
        rc = run_command(command)


        # PBS and SLURM were not tested, so the execution remains in a while loop waiting for the final output file
        ants_dbm_template = template_folder + \
            '/nlin/2/average/template_sharpen_shapeupdate.nii.gz'
        if not cluster_type in ['local', 'sge']:
            import time
            filename = ants_dbm_template
            while True:
                try:
                    with open(filename, 'rb') as _:
                        break
                except IOError:
                    time.sleep(3)

        # verify that all outputs are present
        if not os.path.isfile(ants_dbm_template):
            raise ValueError(ants_dbm_template+" doesn't exists.")

        affine_list = []
        warp_list = []
        inverse_warp_list = []
        warped_anat_list = []

        i = 0
        for file in merged:
            file = str(file)
            filename_template = pathlib.Path(file).name.rsplit(".nii")[0]
            anat_to_template_inverse_warp = '%s/nlin/2/transforms/%s_1InverseWarp.nii.gz' % (
                template_folder, filename_template, )
            if not os.path.isfile(anat_to_template_inverse_warp):
                raise ValueError(
                    anat_to_template_inverse_warp+" file doesn't exists.")
            anat_to_template_warp = '%s/nlin/2/transforms/%s_1Warp.nii.gz' % (
                template_folder, filename_template, )
            if not os.path.isfile(anat_to_template_warp):
                raise ValueError(anat_to_template_warp+" file doesn't exists.")
            anat_to_template_affine = '%s/nlin/2/transforms/%s_0GenericAffine.mat' % (
                template_folder, filename_template, )
            if not os.path.isfile(anat_to_template_affine):
                raise ValueError(anat_to_template_affine
                                 + " file doesn't exists.")
            warped_anat = '%s/nlin/2/resample/%s.nii.gz' % (
                template_folder, filename_template, )
            if not os.path.isfile(warped_anat):
                raise ValueError(warped_anat
                                 + " file doesn't exists.")
            inverse_warp_list.append(anat_to_template_inverse_warp)
            warp_list.append(anat_to_template_warp)
            affine_list.append(anat_to_template_affine)
            warped_anat_list.append(warped_anat)
            i += 1

        setattr(self, 'warped_image', ants_dbm_template)
        setattr(self, 'affine_list', affine_list)
        setattr(self, 'warp_list', warp_list)
        setattr(self, 'inverse_warp_list', inverse_warp_list)
        setattr(self, 'warped_anat_list', warped_anat_list)

        return runtime

    def _list_outputs(self):
        return {'warped_image': getattr(self, 'warped_image'),
                'affine_list': getattr(self, 'affine_list'),
                'warp_list': getattr(self, 'warp_list'),
                'inverse_warp_list': getattr(self, 'inverse_warp_list'),
                'warped_anat_list': getattr(self, 'warped_anat_list'), }
