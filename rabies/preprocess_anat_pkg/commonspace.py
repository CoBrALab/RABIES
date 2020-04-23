#workflow inspired from https://nipype.readthedocs.io/en/latest/users/examples/smri_antsregistration_build_template.html
import os
from nipype.interfaces import utility as niu
import nipype.interfaces.ants as ants
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.workflows.smri.ants import antsRegistrationTemplateBuildSingleIterationWF
from nipype.interfaces.utility import Function

def init_commonspace_wf(name="antsRegistrationTemplateBuilder"):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['file_list']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['PrimaryTemplate', 'PassiveTemplate', 'Transforms', 'PreRegisterAverage']), name='outputnode')

    datasource = pe.Node(Function(input_names=['InitialTemplateInputs'],
                              output_names=['InitialTemplateInputs', 'ListOfImagesDictionaries','registrationImageTypes', 'interpolationMapping'],
                              function=prep_data),
                     name='datasource')

    #creates an average from the input images as initial target template
    initAvg = pe.Node(interface=ants.AverageImages(), name='initAvg')
    initAvg.inputs.dimension = 3
    initAvg.inputs.normalize = True

    #Define the iterations for template building
    buildTemplateIteration1 = antsRegistrationTemplateBuildSingleIterationWF(
        'iteration01')
    buildTemplateIteration2 = antsRegistrationTemplateBuildSingleIterationWF(
        'iteration02')
    buildTemplateIteration3 = antsRegistrationTemplateBuildSingleIterationWF(
        'iteration03')

    workflow.connect(inputnode, "file_list", datasource, "InitialTemplateInputs")
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

    #create a list of dictionaries of the input files
    ListOfImagesDictionaries = []
    for file in InitialTemplateInputs:
        ListOfImagesDictionaries.append({'anat':file})

    return InitialTemplateInputs, ListOfImagesDictionaries,registrationImageTypes, interpolationMapping
