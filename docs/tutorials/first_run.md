# Your first RABIES run

In this tutorial you will take two raw mouse fMRI scans and carry them all the
way through RABIES: preprocessing, confound correction, and a whole-brain
connectivity analysis. Along the way you will look at the quality control
images RABIES produces, so that by the end you will have seen every part of the
software in action.

You do not need to understand every option you type. The point is to build a
feel for the shape of a RABIES run. Everything you use here is explained
elsewhere, and this page links to the relevant pages as it goes.

## What you will end up with

Three output folders that mirror the three stages of the software, and a
connectivity matrix computed from the two scans:

```text
rabies_tutorial/
├── bids_inputs/                  the data you download in step 2
├── preprocess_outputs/           step 3
├── confound_correction_outputs/  step 5
└── analysis_outputs/             step 6
```

```{admonition} How long this takes
:class: important

Steps 1 and 2 take a few minutes. Step 3 is real image processing on real
brains, and is the long one — expect it to run for an hour or more, depending
on how many cores your machine has. Steps 5 and 6 take a couple of minutes
each. Start step 3 and come back to it.
```

## Step 1: Install the container

RABIES has a lot of non-Python dependencies, so you will use the prebuilt
container. Install [Apptainer](https://apptainer.org/docs/user/main/quick_start.html),
then build the RABIES image:

```sh
apptainer build rabies.sif docker://ghcr.io/cobralab/rabies:latest
```

This downloads roughly a gigabyte and takes a few minutes. When it finishes you
will have a single file, `rabies.sif`, which contains the entire software
environment.

Check that it runs:

```sh
apptainer run rabies.sif --help
```

You should see the RABIES usage message, listing the three processing stages:
`preprocess`, `confound_correction` and `analysis`. That is the structure of
the whole tutorial.

```{tip}
If you are on macOS or Windows, or Apptainer is not available to you, use
Docker instead — see [How to install RABIES](../how_to/install.md). Every
`apptainer run ... -B src:dst` below becomes `docker run ... -v src:dst`.
```

## Step 2: Get the example data

Make a working directory and download the two-subject example dataset:

```sh
mkdir -p rabies_tutorial && cd rabies_tutorial
curl -L -o test_dataset.zip https://zenodo.org/records/8349029/files/test_dataset.zip
unzip test_dataset.zip -d bids_inputs
```

Look at what you got:

```sh
find bids_inputs -name '*.nii.gz'
```

Two subjects, each with one anatomical scan and one resting-state functional
scan:

```text
bids_inputs/
├── sub-PHG001/ses-3/
│   ├── anat/sub-PHG001_ses-3_acq-RARE_T2w.nii.gz
│   └── func/sub-PHG001_ses-3_task-rest_acq-EPI_run-1_bold.nii.gz
└── sub-PHG002/ses-3/
    ├── anat/sub-PHG002_ses-3_acq-RARE_T2w.nii.gz
    └── func/sub-PHG002_ses-3_task-rest_acq-EPI_run-1_bold.nii.gz
```

Notice the file names. The `sub-`, `ses-`, `task-` and `T2w`/`bold` parts are
not decoration — this is the [BIDS](https://bids.neuroimaging.io/) naming
standard, and it is how RABIES finds the scans and pairs each functional image
with the anatomical image from the same session. Your own data will need to
look like this too. See [Input data requirements](../reference/bids_inputs.md).

## Step 3: Preprocess

Create the output folder, then run the preprocessing stage:

```sh
mkdir -p preprocess_outputs

apptainer run \
  -B $PWD/bids_inputs:/inputs:ro \
  -B $PWD/preprocess_outputs:/outputs \
  rabies.sif -p MultiProc \
  preprocess /inputs /outputs \
  --commonspace_reg masking=true,fast_commonspace=true
```

Three things are happening in that command:

- `-B $PWD/bids_inputs:/inputs:ro` makes your input folder visible inside the
  container under the name `/inputs`, read-only. The container cannot see any
  path you do not bind, which is why the RABIES arguments say `/inputs` and
  `/outputs` rather than the paths on your machine.
- `-p MultiProc` runs the pipeline in parallel across the cores of your
  machine. Leave it out and everything runs one step at a time.
- `--commonspace_reg masking=true,fast_commonspace=true` registers each scan
  straight to the reference atlas. RABIES would otherwise first build a
  study-specific template out of your scans, which is the better choice for a
  real study but is not worth it for two subjects.

You did not specify an atlas. The container comes with the DSURQE mouse atlas
already installed, and RABIES uses it by default — which is what makes this
command so short.

Let it run. You will see a stream of Nipype log lines as each node in the
workflow completes.

## Step 4: Look at the quality control report

This is the step people skip, and it is the one that catches problems. Open the
report folder:

```sh
ls preprocess_outputs/preprocess_QC_report/
```

Each subfolder holds PNG images for one registration step. Open
`preprocess_QC_report/commonspace_reg_wf.Anat2Atlas/` first — this is where the
alignment to the reference atlas landed, because you passed
`fast_commonspace=true`:

```{figure} ../pics/atlas_registration.png
:alt: Overlap of the study template on the reference atlas

The alignment of your data (top) with the reference atlas (bottom). The
outlines should follow the same anatomy.
```

Then open `preprocess_QC_report/EPI2Anat/`, which shows each functional scan
aligned onto its own anatomical scan:

```{figure} ../pics/sub-MFC068_ses-1_task-rest_acq-EPI_run-1_bold_registration.png
:alt: Overlap of a functional scan on its anatomical scan

The functional image (top) matched to the anatomical image (bottom) from the
same session. This registration is also what corrects susceptibility
distortion.
```

In both figures the brain edges should line up. On this example dataset they
will. When they do not on your own data, that is a registration failure, and
[How to troubleshoot registration](../how_to/troubleshoot_registration.md)
tells you which parameters to change.

Every folder in the report is described in
[Preprocessing QC outputs](../reference/qc_outputs.md).

## Step 5: Correct confounds

Preprocessed data is not yet ready to analyse: head motion and physiological
noise produce correlations that look exactly like connectivity. The confound
correction stage removes them.

```sh
mkdir -p confound_correction_outputs

apptainer run \
  -B $PWD/bids_inputs:/inputs:ro \
  -B $PWD/preprocess_outputs:/outputs \
  -B $PWD/confound_correction_outputs:/cc_outputs \
  rabies.sif -p MultiProc \
  confound_correction /outputs /cc_outputs \
  --nuisance_regressors mot_6 WM_signal CSF_signal \
  --smoothing_filter 0.3
```

Note that you bind `bids_inputs` again even though it is not in the RABIES
arguments. Each stage reads the file paths recorded by the previous one, so
every path used in step 3 has to stay reachable, at the same location, for the
rest of the pipeline.

The two options you passed:

- `--nuisance_regressors mot_6 WM_signal CSF_signal` models the signal using
  the six head motion parameters plus the mean white-matter and CSF signals,
  and subtracts what it modelled.
- `--smoothing_filter 0.3` applies 0.3 mm Gaussian spatial smoothing.

That is a deliberately modest correction. Choosing a correction strategy for a
real dataset is its own task, covered in
[How to optimise your confound correction strategy](../how_to/optimise_confound_correction.md).

The cleaned timeseries land in
`confound_correction_outputs/confound_correction_datasink/cleaned_timeseries/`.

## Step 6: Analyse

Now compute a whole-brain connectivity matrix. RABIES will extract a timecourse
from every region of the atlas parcellation and correlate every pair:

```sh
mkdir -p analysis_outputs

apptainer run \
  -B $PWD/bids_inputs:/inputs:ro \
  -B $PWD/preprocess_outputs:/outputs \
  -B $PWD/confound_correction_outputs:/cc_outputs \
  -B $PWD/analysis_outputs:/analysis_outputs \
  rabies.sif -p MultiProc \
  analysis /cc_outputs /analysis_outputs \
  --FC_matrix --ROI_type parcellated
```

Open the result:

```sh
ls analysis_outputs/commonspace_analysis_datasink/matrix_fig/
```

There is one PNG per scan, showing the correlation between every pair of atlas
regions. The bright block structure along the diagonal is what functional
organisation looks like in this representation: regions near each other, and
regions belonging to the same network, fluctuate together.

The same values are in `commonspace_analysis_datasink/matrix_data_file/` as a `.pkl` file
holding a 2D NumPy array, with rows and columns ordered by atlas label number,
ready to take into your own statistics.

## What you have done

You have run all three RABIES stages and produced a connectivity estimate from
raw scans. Specifically, you:

- built the container and confirmed it runs
- preprocessed two subjects, registering them into a common atlas space
- inspected the registration quality control images
- removed motion and physiological confounds from the timeseries
- computed a whole-brain connectivity matrix

Run through it a second time on the same data. The commands will make more
sense the second time, and the whole sequence will take you a few minutes of
typing.

## Where to go next

For your own data, the decisions the tutorial made for you become yours to
make:

- Your data has to be in BIDS format — see
  [Input data requirements](../reference/bids_inputs.md).
- Your registrations will not always succeed on the first try — see
  [How to troubleshoot registration](../how_to/troubleshoot_registration.md).
- Confound correction should be tuned to your data, not copied from a tutorial
  — see [How to optimise your confound correction strategy](../how_to/optimise_confound_correction.md)
  and [How to assess data quality](../how_to/assess_data_quality.md).
- To understand what the preprocessing actually did, read
  [The preprocessing workflow](../explanation/preprocessing.md).
```
