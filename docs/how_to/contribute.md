# How to contribute to RABIES

RABIES aims to provide an accessible tool responding to growing needs across
the preclinical fMRI community. This effort should be community-driven, and
community involvement is paramount in several respects:

- adapting and maintaining **accessibility** for users across the broader community
- **reproducibility and transparency**, as well as scientific scrutiny and rigour
- defining and incorporating **best practices** across image processing, analysis and quality control
- leveraging appropriate expertise for the **integration of new tools**

Suggestions for improvements can be shared through the GitHub
[issues system](https://github.com/CoBrALab/RABIES/issues) and
[discussion board](https://github.com/CoBrALab/RABIES/discussions).

## Set up a development environment

Install RABIES locally rather than working through a container — testing and
debugging are much smoother that way. Install the dependencies listed in
[`dependencies.txt`](https://github.com/CoBrALab/RABIES/blob/master/dependencies.txt),
then install RABIES from a clone of the repository into a Python environment of
your choice:

```sh
git clone --recursive https://github.com/CoBrALab/RABIES.git
cd RABIES
python setup.py install
```

```{note}
Use `--recursive`. RABIES pulls in submodules, and a clone without them will
not run.
```

### Working inside a container instead

You can run operations in a container to avoid installing dependencies by hand,
at the cost of losing interactive debugging (Spyder and similar will not be
available).

Open a persistent container instance with `-d --entrypoint sh --name mycontainer`,
binding the paths you need with `-v`:

```sh
docker run -it -v $PWD:/work_dir -v /path_to_local_RABIES_package:/RABIES:ro \
  --rm --entrypoint sh -d --name mycontainer rabies:local_testing
```

Then execute commands inside it:

```sh
docker exec mycontainer micromamba run $COMMAND
```

To check for errors, `$COMMAND` can be `error_check_rabies.py --complete`.

To test your changes you must reinstall RABIES inside the container:

```{code-block} sh
:caption: Reinstalling the package inside a running container

mkdir -p /tmp/RABIES
# copy the upgraded package into the container
# (/RABIES must be bound to your local package with -v)
rsync -avz /RABIES/* /tmp/RABIES/.
cd /tmp/RABIES
python setup.py install
```

Compile these into a `.sh` script and run it in place of `$COMMAND` above.

## Submit a pull request

1. **Fork** the RABIES repository on GitHub.
2. **Clone** your fork, with `--recursive` to pull the submodules.
3. **Branch**: `git checkout -b my_new_branch`, with a name describing the
   change. You are ready to modify the code.
4. **Test and debug.** Install your updated package with
   `python setup.py install` in a proper development environment. Test the
   workflow with specific parameters by editing `debug_workflow.py` and running
   it in debug mode (see below).

   ```{important}
   Before committing, confirm that `error_check_rabies.py --complete` finishes
   with no errors.
   ```
5. **Commit, push and open a pull request** from your fork to the original
   repository.

## Debug interactively with Spyder

1. Open `debug_workflow.py` in Spyder.
2. Find the scripts in your local installation and add breakpoints. Running
   `import rabies; os.path.abspath(rabies.__file__)` gives you the path to the
   installed package's `__init__.py`, and from there you can locate the file of
   interest.
3. Run `debug_workflow.py` in debug mode until it reaches the breakpoint.

## Add a new module to the pipeline

RABIES workflows are structured with [Nipype](https://nipype.readthedocs.io/en/latest/)
(see also the [Nipype tutorial](https://miykael.github.io/nipype_tutorial/)).

### Write the function as a Nipype interface

```python
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


class NewInterfaceInputSpec(BaseInterfaceInputSpec):
    # select an appropriate input type with traits.type (Dict, File, Int, ...)
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

### Wrap it in a node

```python
from .other_script import NewInterface  # if the interface is in a different script
from nipype.pipeline import engine as pe

new_interface_node = pe.Node(NewInterface(),
                             name='new_interface')
```

A node can also be built from any plain Python function, without writing an
interface:

```python
from nipype.pipeline import engine as pe
from nipype.interfaces.utility import Function

new_function_node = pe.Node(Function(input_names=['input_1', 'input_2'],
                                     output_names=['output_1', 'output_2'],
                                     function=NewFunction),
                            name='new_function')
```

### Connect the node into a workflow

Once the node carries out the operation you want, integrate it by linking its
inputs and outputs to other nodes. Here is a complete minimal workflow, the one
that performs slice timing correction:

```{code-block} python
:caption: A Nipype workflow conducting slice timing correction
:linenos:

from nipype.pipeline import engine as pe
from nipype.interfaces.utility import Function
from nipype.interfaces import utility as niu


def init_bold_stc_wf(name='bold_stc_wf'):

    workflow = pe.Workflow(name=name)
    # intermediate node storing the workflow inputs
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_file']), name='inputnode')
    # intermediate node storing the workflow outputs
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['stc_file']), name='outputnode')

    slice_timing_correction_node = pe.Node(
        Function(input_names=['in_file', 'tr', 'tpattern', 'stc_axis',
                              'interp_method', 'rabies_data_type'],
                 output_names=['out_file'],
                 function=slice_timing_correction),
        name='slice_timing_correction', mem_gb=1.5 * opts.scale_min_memory)

    # feed the inputnode into the STC node, and STC outputs into the outputnode
    workflow.connect([
        (inputnode, slice_timing_correction_node, [('bold_file', 'in_file')]),
        (slice_timing_correction_node,
            outputnode, [('out_file', 'stc_file')]),
    ])
    return workflow
```

Most contributions integrate a new interface into a pre-existing workflow
rather than creating a new one, in which case the work is connecting the right
nodes to your interface.

```{seealso}
[Workflow reference](../reference/workflows.md) documents the existing
workflows and links to their source.
```
