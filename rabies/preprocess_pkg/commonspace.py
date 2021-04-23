from nipype.interfaces import utility as niu
import nipype.interfaces.ants as ants
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.utility import Function
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


class ANTsDBMInputSpec(BaseInterfaceInputSpec):
    moving_image = traits.List(exists=True, mandatory=True,
                            desc="List of anatomical images used for commonspace registration.")
    output_folder = traits.Str(
        exists=True, mandatory=True, desc="Path to output folder.")
    template_anat = File(exists=True, mandatory=True,
                         desc="Reference anatomical template to define the target space.")
    cluster_type = traits.Str(
        exists=True, mandatory=True, desc="Choose the type of cluster system to submit jobs to. Choices are local, sge, pbs, slurm.")
    walltime = traits.Str(
        exists=True, mandatory=True, desc="Option for job submission specifying requested time per pairwise registration.")
    memory_request = traits.Str(
        exists=True, mandatory=True, desc="Option for job submission specifying requested memory per pairwise registration.")
    local_threads = traits.Int(
        exists=True, mandatory=True, desc="Number of threads to run in parallel if cluster_type is local.")


class ANTsDBMOutputSpec(TraitedSpec):
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


class ANTsDBM(BaseInterface):
    """
    Runs commonspace registration using ants_dbm.
    """

    input_spec = ANTsDBMInputSpec
    output_spec = ANTsDBMOutputSpec

    def _run_interface(self, runtime):
        import os
        import pandas as pd
        import pathlib
        from rabies.preprocess_pkg.utils import run_command

        cwd = os.getcwd()
        template_folder = self.inputs.output_folder+'/ants_dbm_outputs/'

        if os.path.isdir(template_folder):
            # remove previous run
            command = 'rm -r %s' % (template_folder,)
            rc = run_command(command)
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

        command = 'cd %s ; ants_dbm.sh %s %s %s %s %s %s' % (
            template_folder, csv_path, self.inputs.template_anat, self.inputs.cluster_type, self.inputs.walltime, self.inputs.memory_request, self.inputs.local_threads)
        rc = run_command(command)

        # verify that all outputs are present
        ants_dbm_template = template_folder + \
            '/output/secondlevel/secondlevel_template0.nii.gz'
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
            anat_to_template_inverse_warp = '%s/output/secondlevel/secondlevel_%s%s1InverseWarp.nii.gz' % (
                template_folder, filename_template, str(i),)
            if not os.path.isfile(anat_to_template_inverse_warp):
                raise ValueError(
                    anat_to_template_inverse_warp+" file doesn't exists.")
            anat_to_template_warp = '%s/output/secondlevel/secondlevel_%s%s1Warp.nii.gz' % (
                template_folder, filename_template, str(i),)
            if not os.path.isfile(anat_to_template_warp):
                raise ValueError(anat_to_template_warp+" file doesn't exists.")
            anat_to_template_affine = '%s/output/secondlevel/secondlevel_%s%s0GenericAffine.mat' % (
                template_folder, filename_template, str(i),)
            if not os.path.isfile(anat_to_template_affine):
                raise ValueError(anat_to_template_affine
                                 + " file doesn't exists.")
            warped_anat = '%s/output/secondlevel/secondlevel_template0%s%sWarpedToTemplate.nii.gz' % (
                template_folder, filename_template, str(i),)
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


# workflow inspired from https://nipype.readthedocs.io/en/latest/users/examples/smri_antsregistration_build_template.html
def init_commonspace_wf(name="antsRegistrationTemplateBuilder"):
    # from nipype.workflows.smri.ants import antsRegistrationTemplateBuildSingleIterationWF

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['file_list']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=[
                         'PrimaryTemplate', 'PassiveTemplate', 'Transforms', 'PreRegisterAverage']), name='outputnode')

    datasource = pe.Node(Function(input_names=['InitialTemplateInputs'],
                                  output_names=['InitialTemplateInputs', 'ListOfImagesDictionaries',
                                                'registrationImageTypes', 'interpolationMapping'],
                                  function=prep_data),
                         name='datasource')

    # creates an average from the input images as initial target template
    initAvg = pe.Node(interface=ants.AverageImages(), name='initAvg')
    initAvg.inputs.dimension = 3
    initAvg.inputs.normalize = True

    # Define the iterations for template building
    buildTemplateIteration1 = antsRegistrationTemplateBuildSingleIterationWF(
        'iteration01')
    buildTemplateIteration2 = antsRegistrationTemplateBuildSingleIterationWF(
        'iteration02')
    buildTemplateIteration3 = antsRegistrationTemplateBuildSingleIterationWF(
        'iteration03')

    workflow.connect(inputnode, "file_list", datasource,
                     "InitialTemplateInputs")
    workflow.connect(datasource, "InitialTemplateInputs", initAvg, "images")

    workflow.connect(initAvg, 'output_average_image', buildTemplateIteration1,
                     'inputspec.fixed_image')
    workflow.connect(datasource, 'ListOfImagesDictionaries',
                     buildTemplateIteration1, 'inputspec.ListOfImagesDictionaries')
    workflow.connect(datasource, 'registrationImageTypes', buildTemplateIteration1,
                     'inputspec.registrationImageTypes')
    workflow.connect(datasource, 'interpolationMapping', buildTemplateIteration1,
                     'inputspec.interpolationMapping')

    '''
    #the template created from the previous iteration becomes the new target template
    workflow.connect(buildTemplateIteration1, 'outputspec.template',
                     buildTemplateIteration2, 'inputspec.fixed_image')
    workflow.connect(datasource, 'ListOfImagesDictionaries',
                     buildTemplateIteration2, 'inputspec.ListOfImagesDictionaries')
    workflow.connect(datasource, 'registrationImageTypes', buildTemplateIteration2,
                     'inputspec.registrationImageTypes')
    workflow.connect(datasource, 'interpolationMapping', buildTemplateIteration2,
                     'inputspec.interpolationMapping')
    #the template created from the previous iteration becomes the new target template
    workflow.connect(buildTemplateIteration2, 'outputspec.template',
                     buildTemplateIteration3, 'inputspec.fixed_image')
    workflow.connect(datasource, 'ListOfImagesDictionaries',
                     buildTemplateIteration3, 'inputspec.ListOfImagesDictionaries')
    workflow.connect(datasource, 'registrationImageTypes', buildTemplateIteration3,
                     'inputspec.registrationImageTypes')
    workflow.connect(datasource, 'interpolationMapping', buildTemplateIteration3,
                     'inputspec.interpolationMapping')
    '''

    workflow.connect(buildTemplateIteration1, 'outputspec.template', outputnode,
                     'PrimaryTemplate')
    workflow.connect(buildTemplateIteration1,
                     'outputspec.passive_deformed_templates', outputnode,
                     'PassiveTemplate')
    workflow.connect(buildTemplateIteration1,
                     'outputspec.transforms_list', outputnode,
                     'Transforms')
    workflow.connect(initAvg, 'output_average_image', outputnode,
                     'PreRegisterAverage')

    return workflow


def prep_data(InitialTemplateInputs):
    interpolationMapping = {
        'anat': 'Linear'
    }

    registrationImageTypes = ['anat']

    # create a list of dictionaries of the input files
    ListOfImagesDictionaries = []
    for file in InitialTemplateInputs:
        ListOfImagesDictionaries.append({'anat': file})

    return InitialTemplateInputs, ListOfImagesDictionaries, registrationImageTypes, interpolationMapping
