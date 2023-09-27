# Contributing to RABIES

RABIES aims to provide an accessible tool responding to growing needs across the preclinical fMRI community. This effort should be community-driven, and community involvement will be commensurate in achieving this goal in several respects:
- adapting and maintaining accessibility for uses across the broader community
- reproducibility and transparency, as well as scientific scrutiny and rigor
- maintaining best practices across the different aspects of image processing and analysis, as well as quality control
- bringing in appropriate expertise for the integration of new tools

Suggestions for improvements can be shared using the Github [issues system](https://github.com/CoBrALab/RABIES/issues) and [discussion board](https://github.com/CoBrALab/RABIES/discussions). Additionally, contributions from developpers will be important in incorporating novel features as well as addressing revisions. This page provides preliminary guidelines for getting started as a RABIES developper, and covers: setting a developper environment, submitting a pull request, testing and debugging, and basic instructions for adding a new module to the pipeline. We recommend discussing your proposed updates on the Github discussion board or issues prior to creating a pull request. Thank you for your support!

## Dev environment

For development, it is recommend to install RABIES locally, as this will make the testing and debugging process smoother. This requires installing the dependencies listed in dependencies.txt, and then installing RABIES in an appropriate python environment (e.g. using anaconda) from the Github repository. This can be done by cloning the repository, and then running ```python setup.py install```.

### ...using a container

It is possible to do development using an interactive container to avoid installing dependencies manually. An interactive shell can be opened through Docker with ```docker run -it --rm --entrypoint bash ...``` or through Singularity with ```singularity shell ...```. From inside the container, the RABIES package can be upgraded by running ```python setup.py install``` on your local version of the package. It is then possible to run scripts, or test the whole pipeline with ```error_check_rabies.py```. However, it won't be possible to use an interface for debugging (e.g. Spyder).

## Instructions to create a pull request

1. On github, fork the RABIES repository to have your own copy. 
2. Clone your repository to carry out local modifications and testing. Use the `--recursive` option to download the submodules together with the main RABIES package.
3. Make your modifications to the code.
4. Testing and debugging: install your updated version of the package with ```python setup.py install```, using a proper dev environment (see above). Your can test the workflow with specific parameters by editing the ```debug_workflow.py``` script, and executing in debug mode with Spyder (see below). Before commiting changes, make sure that running ```error_check_rabies.py``` completes with no error.
5. Commit and push your modifications to Github, and create a pull request from your forked repo to the original.

### Interactive debugging with Spyder and debug_workflow.py

Here are some recommendations for debugging using Spyder:
1. open the debug_workflow.py file in Spyder
2. find the scrips with your local installation to add breakpoints for debugging. Using `import rabies; os.path.abspath(rabies.__file__)` will provide the path to the __init__.py file of your installed package, and from there you can find file of interest and add a breakpoint where desired.
3. execute debug_workflow.py in debug mode, and run until it finds the breakpoint, and debug from there.


## Creation of a new module and integration within a Nipype workflow

RABIES' workflow is structured using Nipype (for more info on Nipype, see online [documentation](https://nipype.readthedocs.io/en/latest/) and [tutorial](https://miykael.github.io/nipype_tutorial/)). Preferably, a new function should be created as a Nipype interface, which has the following syntax:

```python
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
class NewInterfaceInputSpec(BaseInterfaceInputSpec):
    # you must select an appropriate input type with traits.type (can be Dict, File, Int, ...)
    input_str = traits.Str(exists=True, mandatory=True,
                          desc="An input string.")

class NewInterfaceOutputSpec(TraitedSpec):
    out_file = File(
        exists=True, desc="An output file.")


class NewInterface(BaseInterface):
    """
    Describe your module.
    """

    input_spec = NewInterfaceInputSpec
    output_spec = NewInterfaceOutputSpec

    def _run_interface(self, runtime):
        input_str = self.inputs.input_str

        '''
        YOUR CODE
        '''

        setattr(self, 'out_file', out_file)

        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file')}


```

You can then create a Nipype node for your interface: 
```python
from .other_script import NewInterface # import your interface if from a different script
from nipype.pipeline import engine as pe

new_interface_node = pe.Node(NewInterface(),
                            name='new_interface')
```

Instead of an interface, it is also possible to create a Nipype node from any python function:
```python
from nipype.pipeline import engine as pe
from nipype.interfaces.utility import Function

new_function_node = pe.Node(Function(input_names=['input_1', 'input_2', ...],
                                                output_names=['output_1', 'output_2', ...],
                                                function=NewFunction),
                                        name='new_function')
```

After creating a node which can carry the desired operation, it must be integrated within a workflow by linking up the inputs and outputs with other nodes. Below is an example of a simple workflow which conducts slice-timing correction:

```python
from nipype.pipeline import engine as pe
from nipype.interfaces.utility import Function
from nipype.interfaces import utility as niu

# this function creates and return a Nipype workflow which conducts slice timing correction
def init_bold_stc_wf(name='bold_stc_wf'):

    workflow = pe.Workflow(name=name) # creating a new Nipype workflow
    # creating an intermediate node for storing inputs to the workflow
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file']), name='inputnode')
    # creating an intermediate node for storing outputs to the workflow
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['stc_file']), name='outputnode')

    # preparing the node conducting STC
    slice_timing_correction_node = pe.Node(Function(input_names=['in_file', 'tr', 'tpattern', 'stc_axis', 
                                                        'interp_method', 'rabies_data_type'],
                                                    output_names=[
                                                        'out_file'],
                                                    function=slice_timing_correction),
                                            name='slice_timing_correction', mem_gb=1.5*opts.scale_min_memory)

    # linking up the inputnode to provide inputs to the STC node, and outputs from STC to the outputnode of the workflow
    workflow.connect([
        (inputnode, slice_timing_correction_node, [('bold_file', 'in_file')]),
        (slice_timing_correction_node,
            outputnode, [('out_file', 'stc_file')]),
    ])
    return workflow

```

This example demonstrates the basic syntax of a Nipype workflow. Most likely, a new interface will be integrated as part of a pre-existing workflow (instead of creating a new one), in which case the right nodes must be linked up with the new interface.
